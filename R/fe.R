## Linear models with absorbed fixed effects: design_lm(y ~ x | fe1 + fe2).
##
## Frisch-Waugh-Lovell on aggregates instead of demeaning. Model:
## y = X b + F a + e, weights w, F = [F_1 ... F_K] the dummies of K fixed
## effects (FE). With S = F'WX, s = F'Wy and K = F'WF (diagonal blocks
## diag(n_k), off-diagonal blocks the weighted crosstabs N_kl),
##
##   X'MX = X'WX - S' K^- S,   X'My = X'Wy - S' K^- s,
##
## for any generalized inverse K^- (the columns of S and s lie in the range
## of K). Every input is a sum over rows: X'WX and X'Wy are products with
## the design, S and s grouped sums, K counts and crosstabs; the rows are
## never demeaned and the dense design is never formed. The same estimator
## (dm_fe_core()) runs on data in memory (C kernels, one pass for all
## grouped sums: dm_fe_struct_gsum() in src/fe.c) and on remote data
## (the remote methods of the same primitives): one code path, see
## dm_fe_core().
##
## The solve with K eliminates the largest FE (diagonal block D_1)
## analytically and works on the Schur complement over the other FEs,
##
##   R = blockdiag(n_2, ..., n_K) + offdiag(N_kl) - N_1' D_1^-1 N_1,
##
## (G_2 + ... + G_K square, sparse), with right-hand sides
## B_r - N_1' D_1^-1 B_1, B = [S, s]. R is singular: one level per connected
## component of the (FE_1, FE_k) graph is redundant and set to zero (a valid
## normalization, since shifting FE_k and FE_1 in opposite directions within
## such a component leaves the fit unchanged). Each diagonal block of R is a
## weighted graph Laplacian (its rows sum to zero before the reference
## levels are removed), so the reduced system is a symmetric M-matrix. With
## three or more FE, R can be singular even so (FE collinear beyond the
## pairwise graphs, such as firms assigned round robin over years); the
## solver then drops a further level (dm_fe_solve()).
##
## R is solved densely when small (at most 1000 levels: beyond, the dense
## factorization costs more than conjugate gradients), by a pivoted
## Cholesky factorization that sets the levels of pivots at rounding level
## to zero (and gives the rank of the FE), otherwise by
## conjugate gradients, on R formed or applied without forming it:
##
## * Formed (dm_fe_schur() in src/fe.c, the lower triangle by Gustavson's
##   product): when R is sparse (low mobility between the levels of the
##   other FE, as on worker x firm panels where workers rarely change
##   firms), preconditioned by the incomplete Cholesky factor IC(0) of R
##   (two FE) or of its block diagonal (three or more: block Jacobi over
##   the FEs, each block a Laplacian). On such panels Jacobi needs hundreds
##   of iterations where IC(0) needs about ten.
## * Matrix-free: when forming R would cost more than applying it (high
##   mobility: R near dense), R v = A v - N_1' D_1^-1 (N_1 v), two passes
##   over the sparse crosstab N_1 per product, preconditioned by the
##   diagonal of R (computed without forming R); there it converges in
##   about 15 iterations. Forming gives up as soon as R has more nonzeros
##   than N_1 or the cost of IC(0) (the sum of the squared column counts
##   of R) exceeds form_max (option design.matrix.fe_form_max, default 25,
##   calibrated on worker x firm panels: the two paths take the same time
##   near 20 to 40) times nnz(N_1); a probe of 64 columns estimates both
##   first, so that giving up is cheap.
##
## The stopping rule is a normwise backward error below 1e-12 (see
## dm_fe_pcg() in src/fe.c), with the norm of the right-hand side floored
## at that of the first solve of the system (a refinement step's right-hand
## side is at rounding level), and a warm start (Z0) is accepted, for
## sequences of nearby systems (the steps of a GLM). Sparse Cholesky
## (CHOLMOD) is used only on request (option design.matrix.fe_solver =
## "chol"): on worker x firm graphs its fill is close to dense.
##
## A = K^- B gives the FE coefficients without another solve: a = A_y - A_x b.
## The coefficients are refined from the residuals (computed row by row):
## the normal equations X'M e = 0 correct the rounding of the one-pass
## aggregates. Numeric variables far from zero are centered first when that
## is an exact linear transform of the design (as for remote fits, see
## dm_center() in qr.R), which keeps X'MX, hence the covariance
## matrix, accurate.
##
## With FE, the score of cluster c is
## sum_{i in c} w_i e_i (MX)_i = colSums(w e X, group = cl)_c
## - sum_k E_ck Gamma_k, Gamma = K^- S and E_ck the crosstab of w e over
## (cluster, FE_k); the correction vanishes for FEs nested in the clusters
## (their normal equations make each level's sum of w e zero).
##
## Degrees of freedom and small sample corrections follow fixest's defaults
## (ssc(K.adj = TRUE, K.fixef = "nonnested", K.exact = FALSE, G.adj = TRUE),
## t(G - 1) for clustered tests), and singletons are dropped iteratively as
## feols() does by default.

## ---------------------------------------------------------------------------
## Formula and arguments

## Does the formula have a fixed-effect part (y ~ x | fe)?
dm_fe_has_bar <- function(formula) {
    if (!inherits(formula, "formula")) return(FALSE)
    rhs <- formula[[length(formula)]]
    is.call(rhs) && identical(rhs[[1L]], as.name("|"))
}

## Split `y ~ x | a + b^c` into the main formula and the FE terms (each a
## character vector of variable names, several for an interaction a^b).
dm_fe_parse <- function(formula) {
    rhs <- formula[[length(formula)]]
    if (!dm_fe_has_bar(formula))
        stop("The formula must have a fixed-effect part: y ~ x | fe.", call. = FALSE)
    main <- formula
    main[[length(formula)]] <- rhs[[2L]]
    split_plus <- function(e) {
        if (is.call(e) && identical(e[[1L]], as.name("+"))) {
            c(split_plus(e[[2L]]), split_plus(e[[3L]]))
        } else {
            list(e)
        }
    }
    split_hat <- function(e) {
        if (is.name(e)) return(as.character(e))
        if (is.call(e) && identical(e[[1L]], as.name("^")))
            return(c(split_hat(e[[2L]]), split_hat(e[[3L]])))
        stop("Fixed effects must be variable names or interactions a^b.", call. = FALSE)
    }
    fe <- lapply(split_plus(rhs[[3L]]), split_hat)
    names(fe) <- vapply(fe, paste, "", collapse = "^")
    fe <- fe[!duplicated(names(fe))]
    list(main = main, fe = fe)
}

## The cluster variable of `vcov`: NULL for "iid", else a column name.
dm_fe_vcov_arg <- function(vcov) {
    msg <- "`vcov` must be \"iid\" or a one-sided formula naming a cluster variable, such as ~firm."
    if (is.character(vcov) && length(vcov) == 1L) {
        if (identical(vcov, "iid")) return(NULL)
        stop(msg, call. = FALSE)
    }
    if (inherits(vcov, "formula") && length(vcov) == 2L) {
        if (is.name(vcov[[2L]])) return(as.character(vcov[[2L]]))
        stop("Only one-way clustering is supported: `vcov = ~cluster`, one variable.", call. = FALSE)
    }
    stop(msg, call. = FALSE)
}

dm_fe_rm_arg <- function(fixef.rm) {
    if (!checkmate::test_choice(fixef.rm, c("singleton", "none")))
        stop("`fixef.rm` must be \"singleton\" or \"none\".", call. = FALSE)
    fixef.rm
}

## The formula with variables added as terms (FE and cluster variables, for
## the model frame and the rows used).
dm_fe_augment <- function(main, vars) {
    rhs <- main[[length(main)]]
    for (v in vars) rhs <- call("+", rhs, as.name(v))
    main[[length(main)]] <- rhs
    main
}

## Options of the solver: design.matrix.fe_solver ("auto", "dense", "chol"
## or "pcg"), design.matrix.fe_precond ("auto", default: IC(0) when the
## reduced system is formed, Jacobi when it is applied matrix-free;
## "jacobi"; or "ichol", which forms it), design.matrix.fe_form_max (when
## to form the reduced system, see dm_fe_system()) and
## design.matrix.fe_refine (steps of iterative refinement, default 1).
dm_fe_opts <- function() {
    solver <- getOption("design.matrix.fe_solver", "auto")
    if (!checkmate::test_choice(solver, c("auto", "dense", "chol", "pcg")))
        stop("Option `design.matrix.fe_solver` must be \"auto\", \"dense\", \"chol\" or \"pcg\".", call. = FALSE)
    precond <- getOption("design.matrix.fe_precond", "auto")
    if (!checkmate::test_choice(precond, c("auto", "ichol", "jacobi")))
        stop("Option `design.matrix.fe_precond` must be \"auto\", \"ichol\" or \"jacobi\".", call. = FALSE)
    form_max <- getOption("design.matrix.fe_form_max", 25)
    if (!checkmate::test_number(form_max, lower = 0))
        stop("Option `design.matrix.fe_form_max` must be a non-negative number.", call. = FALSE)
    refine <- getOption("design.matrix.fe_refine", 1L)
    if (!checkmate::test_count(refine))
        stop("Option `design.matrix.fe_refine` must be a count.", call. = FALSE)
    list(solver = solver, precond = precond, pcg_tol = 1e-12, pcg_maxit = 10000L,
         refine = as.integer(refine), form_max = form_max, dense_max = 1000L)
}

## ---------------------------------------------------------------------------
## Codes of fixed effects (in memory)

## Integer codes 1..G of a vector without missing values, and its levels
## (the values that occur, sorted, or in level order for a factor).
dm_fe_code <- function(v) {
    if (is.factor(v)) {
        cd <- dm_fe_compact(as.integer(v), nlevels(v))
        return(list(codes = cd$codes, G = cd$G, levels = levels(v)[cd$used]))
    }
    if (is.logical(v)) v <- as.integer(v)
    if (is.numeric(v)) {
        r <- .Call(C_dm_fe_code_int, v, max(4 * length(v), 1e6))
        if (!is.null(r)) return(list(codes = r[[1L]], G = r[[2L]], levels = r[[3L]]))
    }
    u <- sort(unique(v), method = "radix")
    list(codes = match(v, u), G = length(u), levels = u)
}

## Renumber codes so that the levels that occur are 1..G.
dm_fe_compact <- function(cd, G) {
    used <- tabulate(cd, nbins = G) > 0L
    if (all(used)) return(list(codes = cd, G = G, used = used))
    map <- cumsum(used)
    list(codes = map[cd], G = sum(used), used = used)
}

## Codes of an FE (one or several variables, columns of `d`), with its
## levels as a data frame with one column per variable.
dm_fe_code_term <- function(d, vars) {
    if (length(vars) == 1L) {
        cd <- dm_fe_code(d[[vars]])
        lv <- list(cd$levels)
        names(lv) <- vars
        return(list(codes = cd$codes, G = cd$G, levels = as.data.frame(lv, stringsAsFactors = FALSE)))
    }
    parts <- lapply(vars, function(v) dm_fe_code(d[[v]]))
    cd <- parts[[1L]]$codes
    for (k in seq_along(parts)[-1L]) cd <- (cd - 1) * parts[[k]]$G + parts[[k]]$codes
    comb <- dm_fe_code(cd)
    ## the components of each combined level
    code <- comb$levels - 1
    lv <- vector("list", length(vars))
    for (k in rev(seq_along(parts))) {
        Gk <- parts[[k]]$G
        lv[[k]] <- parts[[k]]$levels[code %% Gk + 1]
        code <- code %/% Gk
    }
    names(lv) <- vars
    list(codes = comb$codes, G = comb$G, levels = as.data.frame(lv, stringsAsFactors = FALSE))
}

## Labels of FE levels, as fixest names them (a_b for interactions).
dm_fe_labels <- function(lv) {
    if (ncol(lv) == 1L) return(as.character(lv[[1L]]))
    do.call(paste, c(lapply(lv, as.character), sep = "_"))
}

## ---------------------------------------------------------------------------
## Small helpers on sums

## Component of each level of b in the bipartite graph of t(N_ab).
dm_fe_components <- function(Nt) {
    i <- rep.int(seq_len(ncol(Nt)) - 1L, diff(Nt@p))
    .Call(C_dm_fe_components, i, Nt@i, ncol(Nt), nrow(Nt))
}

## ---------------------------------------------------------------------------
## The reduced system

## Conjugate gradients for R Z = B, preconditioned by IC(0) of the block
## diagonal of R (blocks: the FEs; the whole of R with two FE), or by its
## diagonal when IC(0) breaks down or with precond = "jacobi". R: a sparse
## symmetric matrix (formed), or a "dm_fe_op" (R = A - N diag(d) N',
## applied without being formed, the rows of N mapped to those of R by
## `map`: Jacobi only). Z0: a starting value (a warm
## start; dropped for columns where it is worse than zero) or NULL. bmin: a
## floor for the norm of the right-hand sides in the stopping rule. The
## preconditioner is built once per `cache`.
dm_fe_pcg <- function(R, B, tol, maxit, block, precond = "jacobi", cache = new.env(), Z0 = NULL,
                      bmin = 0) {
    if (is.null(cache$A)) {
        if (inherits(R, "dm_fe_op")) {
            cache$A <- R$A
            cache$N <- R$N
            cache$map <- R$map
            cache$d <- R$d
            precond <- "jacobi"
        } else {
            cache$A <- dm_fe_lower(R)
        }
        cache$precond <- if (precond == "ichol") "ichol" else "jacobi"
        if (precond == "ichol") {
            A <- cache$A
            P <- list(p = A@p, i = A@i, x = A@x)
            if (length(unique(block)) > 1L) {
                ## block diagonal: the entries within an FE
                j <- rep.int(seq_len(ncol(A)), diff(A@p))
                same <- block[A@i + 1L] == block[j]
                P <- list(p = c(0L, cumsum(tabulate(j[same], ncol(A)))), i = A@i[same], x = A@x[same])
            }
            ## plain IC(0): the modified factor (MIC, omega > 0) took 2.5 to 5
            ## times as many iterations on worker x firm graphs
            Lx <- .Call(C_dm_fe_ichol, P$p, P$i, P$x, ncol(A), 0)
            if (!is.null(Lx)) {
                cache$L <- list(p = P$p, i = P$i, x = Lx)
            } else {
                cache$precond <- "jacobi"
            }
        }
    }
    A <- cache$A
    N <- cache$N
    L <- cache$L
    storage.mode(B) <- "double"
    if (!is.null(Z0)) {
        Z0 <- as.matrix(Z0)
        storage.mode(Z0) <- "double"
        if (!identical(dim(Z0), dim(B))) Z0 <- NULL
    }
    r <- .Call(C_dm_fe_pcg, A@p, A@i, A@x, if (!is.null(N)) N@p, if (!is.null(N)) N@i,
               if (!is.null(N)) N@x, cache$map, cache$d, L$p, L$i, L$x, B, Z0,
               as.double(bmin), as.double(tol), as.integer(maxit))
    ## (nullv: on a breakdown, a direction close to a null vector of R)
    list(x = r[[1L]], iter = r[[2L]], resid = r[[3L]], converged = r[[3L]] <= tol,
         nullv = r[[4L]], precond = cache$precond)
}

## The lower triangle of a symmetric sparse matrix, as a dgCMatrix.
dm_fe_lower <- function(R) {
    if (methods::is(R, "symmetricMatrix") && R@uplo == "U") R <- Matrix::t(R)
    R <- methods::as(methods::as(R, "generalMatrix"), "CsparseMatrix")
    if (!methods::is(R, "dsparseMatrix")) R <- methods::as(R, "dMatrix")
    Matrix::tril(R)
}

## Solve the reduced system R Z = B (R symmetric positive definite after the
## reference levels are removed, or semidefinite with a consistent
## right-hand side). R: a sparse matrix, or a "dm_fe_op" (not formed: by
## conjugate gradients with the Jacobi preconditioner). Z0: a starting
## value for conjugate gradients (a warm start: the solution of a nearby
## system), or NULL. bmin: a floor for the norm of the right-hand sides in
## the stopping rule of conjugate gradients; by default the largest norm of
## the right-hand sides solved before with the same `cache` (a refinement
## step's right-hand side is at rounding level: its solution is needed to
## the accuracy of the first solve, not to tol relative to itself).
## Returns Z and solver details.
dm_fe_solve <- function(R, B, opts, block, cache = new.env(), dense_max = opts$dense_max %||% 1000L,
                        Z0 = NULL, bmin = NULL) {
    m <- if (inherits(R, "dm_fe_op")) R$m else nrow(R)
    if (m == 0L) return(list(z = B, solver = "none"))
    solver <- opts$solver
    if (solver == "auto") solver <- if (m <= dense_max) "dense" else "pcg"
    if (inherits(R, "dm_fe_op")) solver <- "pcg"
    if (solver == "dense") {
        ## pivoted Cholesky: R can be singular with three or more FE (FE
        ## collinear beyond the references removed); the levels of pivots
        ## below rounding are set to zero too (a valid normalization)
        M <- as.matrix(R)
        ch <- tryCatch(suppressWarnings(chol(M, pivot = TRUE, tol = 1e-11 * max(diag(M)))),
                       error = function(e) NULL)
        if (!is.null(ch)) {
            rk <- attr(ch, "rank")
            pv <- attr(ch, "pivot")[seq_len(rk)]
            Rk <- ch[seq_len(rk), seq_len(rk), drop = FALSE]
            z <- matrix(0, m, ncol(B))
            z[pv, ] <- backsolve(Rk, backsolve(Rk, B[pv, , drop = FALSE], transpose = TRUE))
            if (all(is.finite(z))) return(list(z = z, solver = "dense", rank = rk))
        }
        solver <- "pcg"
    }
    if (solver == "chol") {
        L <- tryCatch(Matrix::Cholesky(R, perm = TRUE, LDL = FALSE), error = function(e) NULL)
        if (!is.null(L)) {
            z <- as.matrix(Matrix::solve(L, B, system = "A"))
            ok <- all(is.finite(z)) &&
                max(abs(as.matrix(R %*% z) - B)) <= 1e-8 * max(1, max(abs(B)))
            if (ok) return(list(z = z, solver = "chol", nnz_L = length(L@x)))
        }
        solver <- "pcg"
    }
    precond <- opts$precond
    if (precond == "auto") precond <- if (inherits(R, "dm_fe_op")) "jacobi" else "ichol"
    ## R singular (FE collinear beyond the references removed: with three
    ## or more FE, the components of the pairwise graphs do not find every
    ## redundant level): conjugate gradients break down once rounding
    ## excites the null direction. The level with the largest entry of the
    ## direction of the breakdown (close to a null vector) is then set to
    ## zero too, a valid normalization, and the system solved again; the
    ## levels dropped are kept in the cache for later solves.
    if (is.null(bmin)) bmin <- cache$bnorm %||% 0
    cache$bnorm <- max(cache$bnorm, sqrt(colSums(B^2)))
    drop <- cache$drop
    repeat {
        if (is.null(cache$pc)) cache$pc <- new.env(parent = emptyenv())
        use <- if (length(drop)) seq_len(m)[-drop] else seq_len(m)
        sol <- dm_fe_pcg(dm_fe_sub(R, use), B[use, , drop = FALSE], opts$pcg_tol, opts$pcg_maxit,
                         block[use], precond, cache$pc,
                         Z0 = if (!is.null(Z0) && nrow(Z0) == m) Z0[use, , drop = FALSE], bmin = bmin)
        if (sol$converged || is.null(sol$nullv) || length(drop) >= 20L) break
        drop <- sort(c(drop, use[which.max(abs(sol$nullv))]))
        cache$drop <- drop
        cache$pc <- NULL
    }
    if (!sol$converged)
        warning(sprintf("Conjugate gradients did not converge in %d iterations (relative residual %.2g).",
                        sol$iter, sol$resid), call. = FALSE)
    z <- matrix(0, m, ncol(B))
    z[use, ] <- sol$x
    list(z = z, solver = "pcg", precond = sol$precond, iter = sol$iter, resid = sol$resid,
         ndrop = length(drop))
}

## The rows and columns `use` of the reduced system (a sparse matrix or a
## "dm_fe_op").
dm_fe_sub <- function(R, use) {
    m <- if (inherits(R, "dm_fe_op")) R$m else nrow(R)
    if (length(use) == m) return(R)
    if (!inherits(R, "dm_fe_op")) return(R[use, use, drop = FALSE])
    R$A <- R$A[use, use, drop = FALSE]
    mp <- match(R$map + 1L, use) - 1L
    mp[is.na(mp)] <- -1L
    R$map <- mp
    R$m <- length(use)
    R
}

## The FE system K = F'WF, reduced by eliminating the largest FE.
## G: levels per FE; nk: weighted counts per level; crosstab(a, b): t(N_ab).
##
## R = A - N_1' D_1^-1 N_1 with A = blockdiag(n_k) + offdiag(N_kl) is formed
## (dm_fe_schur(), lower triangle) unless forming it would cost more than
## applying it: forming gives up as soon as R has more nonzeros than N_1
## (min(form_max, 1) times), or the cost of its IC(0) factor (the sum of
## its squared column counts) exceeds form_max times nnz(N_1). R is then
## applied matrix-free (a "dm_fe_op"), two passes over N_1 per product.
## R is always formed for the dense solver (m <= dense_max), for "chol"
## and for precond = "ichol".
dm_fe_system <- function(G, nk, crosstab, opts) {
    nfe <- length(G)
    o <- order(G, decreasing = TRUE)
    sys <- list(nfe = nfe, f1 = o[1L], D1 = nk[[o[1L]]], rest = o[-1L], opts = opts)
    if (nfe == 1L) return(sys)
    rest <- sys$rest
    sys$Gr <- G[rest]
    sys$off <- off <- c(0L, cumsum(sys$Gr))
    m <- off[length(off)]
    ## t(N_1k) for each other FE, and references from the (FE_1, FE_k) graphs
    Nt <- lapply(rest, function(k) crosstab(sys$f1, k))
    ncomp <- integer(length(rest))
    ref <- integer(0)
    for (r in seq_along(rest)) {
        comp <- dm_fe_components(Nt[[r]])
        ncomp[r] <- max(comp)
        ref <- c(ref, off[r] + which(!duplicated(comp)))
    }
    sys$Ntr <- if (length(Nt) == 1L) Nt[[1L]] else do.call(rbind, Nt)
    keep <- setdiff(seq_len(m), ref)
    mk <- length(keep)
    ## A = blockdiag(n_k) + offdiag(N_kl), lower triangle, on the levels kept
    ii <- jj <- seq_len(m)
    xx <- unlist(nk[rest], use.names = FALSE)
    if (length(rest) > 1L) {
        for (r in seq_along(rest)[-length(rest)]) for (s in (r + 1L):length(rest)) {
            Ts <- Matrix::summary(crosstab(rest[r], rest[s])) # G_s x G_r
            ii <- c(ii, Ts$i + off[s])
            jj <- c(jj, Ts$j + off[r])
            xx <- c(xx, Ts$x)
        }
    }
    A <- Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(m, m), repr = "C")[keep, keep, drop = FALSE]
    ## N_1' (rows: the levels of the other FE), its rows mapped to those of R
    N <- methods::as(methods::as(sys$Ntr, "generalMatrix"), "CsparseMatrix")
    map <- rep(-1L, m)
    map[keep] <- seq_len(mk) - 1L
    nnzN <- max(length(N@x), mk)
    form <- opts$solver %in% c("dense", "chol") || opts$precond == "ichol" ||
        (opts$solver == "auto" && mk <= (opts$dense_max %||% 1000L))
    fm <- if (form) Inf else opts$form_max %||% 25
    S <- .Call(C_dm_fe_schur, A@p, A@i, A@x, N@p, N@i, N@x, map, 1 / sys$D1,
               as.double(if (is.finite(fm)) min(fm, 1) * nnzN else Inf), as.double(fm * nnzN))
    if (!is.null(S)) {
        sys$R <- methods::new("dsCMatrix", Dim = c(mk, mk), uplo = "L", p = S[[1L]], i = S[[2L]], x = S[[3L]])
    } else {
        sys$R <- structure(list(A = A, N = N, map = map, d = 1 / sys$D1, m = mk), class = "dm_fe_op")
    }
    sys$keep <- keep
    sys$cache <- new.env(parent = emptyenv()) # the preconditioner, built once
    sys$block <- rep(seq_along(rest), sys$Gr)[keep]
    sys$m <- m
    sys$info <- list(m = mk, nnz_R = if (is.null(S)) NA_integer_ else length(S[[3L]]),
                     nnz_N = length(N@x), formed = !is.null(S), ncomp = ncomp)
    sys
}

## A = K^- B (split by FE) and Q = B' K^- B, for B a list of G_k x q
## matrices (grouped sums by FE). `Z0`: a starting value of the solution of
## the reduced system (`Z`, returned, of an earlier call: a warm start), or
## NULL.
dm_fe_apply <- function(sys, B, Z0 = NULL) {
    f1 <- sys$f1
    A1p <- B[[f1]] / sys$D1
    Q <- crossprod(B[[f1]], A1p)
    A <- vector("list", sys$nfe)
    if (sys$nfe == 1L) {
        A[[f1]] <- A1p
        return(list(Q = Q, A = A, info = list(solver = "direct", m = 0L)))
    }
    rest <- sys$rest
    Brr <- do.call(rbind, B[rest]) - as.matrix(sys$Ntr %*% A1p)
    sol <- dm_fe_solve(sys$R, Brr[sys$keep, , drop = FALSE], sys$opts, sys$block, sys$cache, Z0 = Z0)
    Z <- matrix(0, sys$m, ncol(Brr))
    Z[sys$keep, ] <- sol$z
    Zk <- sol$z
    sol$z <- NULL
    Q <- Q + crossprod(Brr, Z)
    A[[f1]] <- (B[[f1]] - as.matrix(Matrix::crossprod(sys$Ntr, Z))) / sys$D1
    for (r in seq_along(rest)) A[[rest[r]]] <- Z[sys$off[r] + seq_len(sys$Gr[r]), , drop = FALSE]
    info <- c(sol, sys$info)
    ## (factored densely: m, the rank of the reduced system)
    if (!is.null(sol$rank)) info$m <- sol$rank
    info$rank <- NULL
    list(Q = Q, A = A, info = info, Z = Zk)
}

## Sequential pivoted Cholesky of a p x p matrix M in column order: a column
## is aliased when its residual variance after the earlier kept columns is
## below tol times `scale` (its variance before the FE are absorbed).
## Returns the kept columns and the upper factor of M[keep, keep].
dm_fe_pivchol <- function(M, scale, tol) {
    p <- ncol(M)
    keep <- integer(0)
    Rf <- matrix(0, 0, 0)
    for (j in seq_len(p)) {
        if (length(keep)) {
            r <- backsolve(Rf, M[keep, j], transpose = TRUE)
            d <- M[j, j] - sum(r^2)
        } else {
            r <- numeric(0)
            d <- M[j, j]
        }
        if (is.finite(d) && scale[j] > 0 && d > tol * scale[j]) {
            k <- length(keep)
            Rn <- matrix(0, k + 1L, k + 1L)
            if (k) {
                Rn[seq_len(k), seq_len(k)] <- Rf
                Rn[seq_len(k), k + 1L] <- r
            }
            Rn[k + 1L, k + 1L] <- sqrt(d)
            Rf <- Rn
            keep <- c(keep, j)
        }
    }
    list(keep = keep, R = Rf)
}

## ---------------------------------------------------------------------------
## The estimator, on the primitives
##
## One implementation for designs in memory and remote designs: every sum
## is a primitive of the fitting core (primitives.R, primitives_lookup.R):
## X'WX and X'Wz (dmp_crossprod), F_k'WX, F_k'W1 and F_k'Wz for all FE
## (dmp_colsums_groups, one pass in memory), crosstabs N_kl (dmp_crosstab),
## the residual column (dmp_mult and lookups of the FE values), F_k'We
## (dmp_groupsum), sums over rows (dm_colsum).
##
## Arguments of dm_fe_core(): the design x, and `cols`, the columns of the
## problem for that design: z (response), w (weights or NULL), codes (per
## FE, a dmc_codes() grouping) and, for clustered standard errors, cl (a
## dmc_codes() grouping of the clusters).

dm_fe_core <- function(x, cols, tol, opts, cluster = FALSE) {
    z <- cols$z
    w <- cols$w
    codes <- cols$codes
    G <- vapply(codes, function(g) g$G, 1L)
    p <- ncol(x)
    icol <- which(colnames(x) == "(Intercept)")
    xcols <- setdiff(seq_len(p), icol)
    px <- length(xcols)
    nfe <- length(G)
    n <- dm_nrows(x)
    ## moments
    wz <- if (is.null(w)) z else w * z
    xtx <- unname(if (is.null(w)) dmp_crossprod(x) else dmp_crossprod(x, dmp_scale(x, w)))
    xty <- unname(dmp_crossprod(x, wz))[, 1L]
    yty <- dm_colsum(wz * z)
    sw <- if (is.null(w)) n else dm_colsum(w)
    ybar <- dm_colsum(wz) / sw
    tss <- dm_colsum(if (is.null(w)) (z - ybar)^2 else w * (z - ybar)^2)
    if (!all(is.finite(xtx)) || !all(is.finite(xty))) stop("NA/NaN/Inf in 'x'", call. = FALSE)
    ## grouped sums and the FE system
    gs <- dmp_colsums_groups(if (is.null(w)) x else dmp_scale(x, w), codes, extra = z)
    nk <- lapply(gs, `[[`, "n")
    if (any(vapply(nk, function(v) any(!(v > 0)), NA)))
        stop("internal error: a fixed-effect level without weight", call. = FALSE)
    B <- lapply(gs, function(g) cbind(g$S[, xcols, drop = FALSE], g$extra[, 1L]))
    ## t(N_ab), G_b x G_a
    crosstab <- function(a, b) dmp_crosstab(codes[[b]], codes[[a]], w)
    sys <- dm_fe_system(G, nk, crosstab, opts)
    ab <- dm_fe_apply(sys, B)
    XtX <- xtx[xcols, xcols, drop = FALSE]
    Xty <- xty[xcols]
    Mfull <- rbind(cbind(XtX, Xty), c(Xty, yty)) - ab$Q
    Mxx <- Mfull[seq_len(px), seq_len(px), drop = FALSE]
    Mxy <- Mfull[seq_len(px), px + 1L]
    sx <- colSums(B[[1L]])[seq_len(px)]
    scale <- diag(XtX) - sx^2 / sum(nk[[1L]])
    pc <- dm_fe_pivchol(Mxx, scale, tol)
    kp <- pc$keep
    b <- rep(NA_real_, px)
    if (length(kp)) b[kp] <- backsolve(pc$R, backsolve(pc$R, Mxy[kp], transpose = TRUE))
    b0 <- b
    b0[is.na(b0)] <- 0
    fe <- lapply(ab$A, function(Ak) drop(Ak[, px + 1L] - Ak[, seq_len(px), drop = FALSE] %*% b0))
    ## residuals e = z - X b - sum_k fe_k[code_k], and their sums
    resid <- function(b0, fe) {
        beta <- numeric(p)
        beta[xcols] <- b0
        f <- dm_fe_fitted(x, beta, fe, codes)
        e <- z - f
        we <- if (is.null(w)) e else w * e
        list(xte = unname(dmp_crossprod(x, we))[, 1L],
             Be = lapply(seq_len(nfe), function(k) dmp_groupsum(we, codes[[k]])[, 1L]),
             rss = dm_colsum(we * e), we = we)
    }
    r <- resid(b0, fe)
    ## iterative refinement: the normal equations X'M e = 0 from residuals
    ## computed row by row correct the rounding of the one-pass aggregates
    for (it in seq_len(opts$refine)) {
        Ae <- dm_fe_apply(sys, lapply(r$Be, function(v) matrix(v, ncol = 1L)))$A
        delta <- numeric(0)
        if (length(kp)) {
            g <- r$xte[xcols]
            for (k in seq_len(nfe)) g <- g - drop(crossprod(B[[k]][, seq_len(px), drop = FALSE], Ae[[k]]))
            delta <- backsolve(pc$R, backsolve(pc$R, g[kp], transpose = TRUE))
            b0[kp] <- b0[kp] + delta
            b[kp] <- b0[kp]
        }
        for (k in seq_len(nfe))
            fe[[k]] <- fe[[k]] + Ae[[k]][, 1L] - drop(ab$A[[k]][, kp, drop = FALSE] %*% delta)
        r <- resid(b0, fe)
    }
    ## degrees of freedom (fixest: K.fixef = "nonnested", K.exact = FALSE)
    fe_dof <- sum(G) - (nfe - 1L)
    ## the exact rank of K (fixest's K.exact = TRUE): with two FE, one level
    ## per connected component is redundant; with more, the count is exact
    ## when the reduced system was factored (a pivoted dense Cholesky, which
    ## gives its rank, or a sparse Cholesky that succeeded); otherwise NA
    fe_dof_exact <- if (nfe == 1L) {
        G
    } else if (nfe == 2L || ab$info$solver %in% c("dense", "chol")) {
        max(G) + ab$info$m
    } else {
        NA_integer_
    }
    rank <- length(kp)
    bread <- matrix(NA_real_, px, px)
    if (rank) bread[kp, kp] <- chol2inv(pc$R)
    out <- list(b = b, kp = kp, xcols = xcols, bread = bread, rss = r$rss, fe = fe,
                rank = rank, G = G, fe_dof = fe_dof, fe_dof_exact = fe_dof_exact,
                df.residual = n - rank - fe_dof, nobs = n, tss = tss,
                tss_within = max(Mfull[px + 1L, px + 1L], 0), info = ab$info)
    if (cluster) {
        cl <- cols$cl
        Gc <- cl$G
        U <- dmp_colsums_groups(dmp_scale(x, r$we), list(cl))[[1L]]$S[, xcols, drop = FALSE]
        ## an FE is nested in the clusters when each of its levels is in one
        ## cluster; otherwise its scores need a correction
        nested <- vapply(seq_len(nfe), function(k) dmp_nested(cl, codes[[k]]), NA)
        for (k in which(!nested)) {
            E <- dmp_crosstab(codes[[k]], cl, r$we) # G_k x Gc
            U <- U - as.matrix(Matrix::crossprod(E, ab$A[[k]][, seq_len(px), drop = FALSE]))
        }
        U <- U[, kp, drop = FALSE]
        nn <- !nested
        K_cl <- rank + if (any(nn)) sum(G[nn]) - (sum(nn) - 1L) else 1L
        adj <- Gc / (Gc - 1) * (n - 1) / (n - K_cl)
        Vc <- matrix(NA_real_, px, px)
        if (rank) {
            bk <- bread[kp, kp, drop = FALSE]
            Vc[kp, kp] <- adj * bk %*% crossprod(U) %*% bk
        }
        out$vcov_cluster <- Vc
        out$G_cluster <- Gc
    }
    out
}

## X beta plus the FE values of the rows' levels, as a column (codes: per
## FE, integer codes, or, with `levels`, key columns matched to levels).
dm_fe_fitted <- function(x, beta, fe, codes, levels = NULL) {
    f <- dmp_mult(x, beta)
    for (k in seq_along(fe))
        f <- f + dmc_lookup(fe[[k]], codes[[k]], levels[[k]], along = list(f))
    f
}

## Map a fit on the centered design X~ = X T back to X: b = T b~ (but the
## intercept column, absorbed by the FE, which take the constant t0'b~),
## bread T bread~ T'. NULL when the aliased columns do not allow it (a kept
## column of X~ uses an aliased column of X).
dm_fe_uncenter <- function(est, Tm, icol) {
    xc <- est$xcols
    Tx <- Tm[xc, xc, drop = FALSE]
    kp <- est$kp
    al <- setdiff(seq_along(xc), kp)
    if (length(al) && length(kp) && any(Tx[al, kp] != 0)) return(NULL)
    bt <- est$b
    bt[is.na(bt)] <- 0
    b <- drop(Tx %*% bt)
    b[al] <- NA_real_
    est$b <- b
    if (length(kp)) {
        Tk <- Tx[kp, kp, drop = FALSE]
        est$bread[kp, kp] <- Tk %*% est$bread[kp, kp, drop = FALSE] %*% t(Tk)
        if (!is.null(est$vcov_cluster))
            est$vcov_cluster[kp, kp] <- Tk %*% est$vcov_cluster[kp, kp, drop = FALSE] %*% t(Tk)
    }
    ## X~ b~ = X b + (t0'b~) 1: the constant goes to the first FE
    if (length(icol)) est$fe[[1L]] <- est$fe[[1L]] + sum(Tm[icol, xc] * bt)
    est
}


## Centering of numeric variables (in memory; see dmp_center()): list(shift, T), with the
## design of the shifted frame X~ = X T and T unit upper triangular, or
## NULL. As dm_center() does for remote designs: only variables whose mean
## is large compared with their spread, and only when such a T exists
## (checked on prototype rows).
dm_fe_center <- function(mt, mf, contrasts) {
    if (!isTRUE(getOption("design.matrix.fe_center", TRUE))) return(NULL)
    if (attr(mt, "intercept") == 0L) return(NULL)
    ## the decision and the shift from a systematic sample of rows (any
    ## shift is exact; a shift near the mean keeps the sums accurate)
    tv <- dm_terms_vars(mt)$names
    nr <- nrow(mf)
    idx <- if (nr > 1e5) unique(round(seq(1, nr, length.out = 1e5))) else seq_len(nr)
    st <- mf[idx, intersect(tv, names(mf)), drop = FALSE]
    attr(st, "terms") <- mt
    m0 <- mf[seq_len(min(nrow(mf), 2000L)), , drop = FALSE]
    attr(m0, "terms") <- NULL
    p <- tryCatch(ncol(design_matrix(mt, m0, contrasts.arg = contrasts)), error = function(e) NA)
    if (is.na(p) || p > 2000L) return(NULL)
    dm_center_shift(st, m0, function(fr) design_matrix(mt, fr, contrasts.arg = contrasts), p)
}

## ---------------------------------------------------------------------------
## design_lm() with fixed effects
##
## Only the first step differs between data in memory and remote data, as
## it does without FE: building the design of the rows with complete model
## variables, with the response and weights as columns (a model frame in
## memory; a remote design and dmp_resolve() remotely). The rows used, the
## FE codes, the removal of singletons, the estimator, the methods and the
## predictions are one code path on the primitives.

#' @section Fixed effects:
#' A formula `y ~ x | fe1 + fe2` absorbs fixed effects, as
#' [fixest::feols()]: each FE is a variable (a factor, character, integer or
#' numeric column) or an interaction `a^b`. The fit never builds the FE
#' dummies or demeans the rows: it solves the normal equations from sums
#' over rows (the design's products, sums within FE levels and crosstabs of
#' the FE), eliminating the largest FE exactly and solving for the others
#' densely (up to 1000 levels) or by conjugate gradients: on the reduced
#' system formed, with an incomplete Cholesky preconditioner, when it is
#' sparse (low mobility, such as workers who rarely change firms), and
#' without forming it, with a Jacobi preconditioner, when it would be near
#' dense (high mobility). The same code fits data in memory and remote
#' data (DuckDB tables or queries, Parquet), with identical results up to
#' rounding.
#'
#' * Rows: rows with a missing value in the model, FE or cluster variables
#'   are dropped, as are rows with zero weight and, with `fixef.rm =
#'   "singleton"` (default), singletons (rows alone in a level of some FE),
#'   iteratively, as `feols()` does. `fit$nsingletons` counts them.
#' * Coefficients: the intercept is absorbed; a column collinear with the
#'   FE or earlier columns is `NA` (fixest's rule: its residual variance is
#'   below `tol` times its variance; `tol` defaults to `1e-10` with FE, as
#'   `collin.tol`). One step of iterative refinement from the residuals is
#'   taken, and numeric variables far from zero are centered when that is
#'   an exact transform of the design, so that coefficients and standard
#'   errors stay accurate.
#' * Standard errors: `vcov = "iid"` (default) or `vcov = ~cl` (clustered by
#'   one variable), with fixest's default small sample corrections
#'   (`ssc(K.adj = TRUE, K.fixef = "nonnested", G.adj = TRUE)`); clustered
#'   tests in `summary()` use `t(G - 1)`. The cluster scores are corrected
#'   for FE not nested in the clusters. `fit$fe_dof` is the number of FE
#'   parameters counted (`sum(levels) - (number of FE - 1)`); `fit$fe_dof_exact`
#'   the exact rank of the FE (fixest's `K.exact = TRUE`), exact for one or
#'   two FE and for three or more when the reduced system was factored
#'   (otherwise `NA`).
#' * Methods: those of a `design_lm` fit (class `c("design_lm_fe",
#'   "design_lm")`); `summary()` reports R-squared and within R-squared.
#'   `fixef()` (the generic of fixest or nlme, when loaded) returns the FE
#'   values, also in `fit$fixef`: one named vector per FE, normalized with
#'   one level per connected set of levels at zero in every FE but the
#'   largest (differences within a connected set are those of fixest).
#'   `predict(fit, newdata)` is `X b` plus the FE values of the rows' levels
#'   (`NA` for levels not in the fit). With remote data, `fitted()` and
#'   `predict()` return tables and `residuals()` a remote design, as for
#'   remote fits without FE.
#' * Not supported with FE: `subset`, offsets, `compress = TRUE`, `se.fit`
#'   and intervals in `predict()`, clustering by more than one variable.
#'   For remote data, the response and `weights` are column names,
#'   `DBI::SQL()` or R expressions, as without FE (expressions with no SQL
#'   translation are evaluated in R, in batches).
#' * Options: `design.matrix.fe_solver` (`"auto"`, `"dense"`, `"pcg"`, or
#'   `"chol"` for a sparse Cholesky factor), `design.matrix.fe_precond`
#'   (`"auto"`, default: incomplete Cholesky when the reduced system is
#'   formed, Jacobi otherwise; `"jacobi"`; or `"ichol"`, which always forms
#'   it), `design.matrix.fe_form_max` (when to form the reduced system:
#'   default 25, `0` never, `Inf` always), `design.matrix.fe_refine` (refinement steps,
#'   default 1) and `design.matrix.fe_center` (default `TRUE`).
#' @param vcov With fixed effects: `"iid"` or a one-sided formula naming a
#'   cluster variable (`~firm`). Without fixed effects only `"iid"`.
#' @param fixef.rm With fixed effects: `"singleton"` (drop singletons
#'   iteratively, as [fixest::feols()]) or `"none"`.
#' @examples
#' ## fixed effects (examples below: without)
#' set.seed(2)
#' panel <- data.frame(id = rep(1:100, each = 10), t = rep(1:10, 100), x = rnorm(1000))
#' panel$y <- panel$x + panel$id / 50 + rnorm(1000)
#' fe <- design_lm(y ~ x | id + t, panel, vcov = ~id)
#' summary(fe)
#' head(fe$fixef$id)
#'
#' @rdname design_lm
#' @name design_lm
NULL

## Entry point from design_lm(): `call` is its matched call, `env` its
## parent frame.
dm_fe_lm <- function(call, env, formula, contrasts, compress, tol, tol_missing, model,
                     vcov, fixef.rm, dots) {
    if (isTRUE(compress)) stop("`compress = TRUE` is not supported with fixed effects.", call. = FALSE)
    if (!is.null(call$offset)) stop("Offsets are not supported with fixed effects.", call. = FALSE)
    if (!is.null(call$subset)) stop("`subset` is not supported with fixed effects.", call. = FALSE)
    checkmate::assert_number(tol, lower = 0, finite = TRUE)
    checkmate::assert_flag(model)
    ## aliasing: fixest's criterion (residual variance below tol times the
    ## variance before the FE), with fixest's default collin.tol = 1e-10
    if (tol_missing) tol <- 1e-10
    pf <- dm_fe_parse(formula)
    clv <- dm_fe_vcov_arg(vcov)
    fixef.rm <- dm_fe_rm_arg(fixef.rm)
    data <- if (is.null(call$data)) NULL else eval(call$data, env)
    ## the design of the rows with complete model variables, with the
    ## response and weights as columns (as for fits without FE, the
    ## construction differs: a model frame, or a remote design)
    base <- if (!is.null(data) && dm_is_remote_source(data)) {
        dm_fe_base_remote(pf, data, call, env, contrasts, dots, clv)
    } else {
        if (length(dots)) stop("Unknown arguments: ", paste(names(dots), collapse = ", "), call. = FALSE)
        dm_fe_base_memory(pf, call, env, contrasts, clv)
    }
    prob <- dm_fe_prepare(base, pf, clv, fixef.rm)
    est <- dm_fe_estimate(prob$x, prob$cols, tol, !is.null(clv))
    fit <- dm_fe_finish(est, colnames(prob$x), pf, prob$levels, clv, prob$nsingle, base$nall)
    w <- prob$cols$w
    fit$weighted <- !is.null(w)
    fit$sumlogw <- if (is.null(w)) 0 else dm_colsum(log(w))
    fit$x <- prob$x
    fit$cols <- prob$cols
    fit$terms <- base$terms
    fit$contrasts <- base$contrasts
    fit$xlevels <- base$xlevels
    fit$spec <- base$spec
    ## (remotely, na.exclude pads the results with every row not used)
    fit$na.action <- if (isTRUE(base$exclude)) dm_remote_exclude(prob$x)
                     else dm_fe_na_action(base$na_act, prob$rows, base$nall)
    fit$method <- switch(dm_engine(prob$x), sql = "SQL aggregates", chunked = "aggregates computed in batches", NULL)
    ## (the data of a remote fit stay in the database: always kept)
    if (!is.null(fit$method)) {
        ## nothing to drop
    } else if (model) {
        fit$model <- prob$x$vars
    } else {
        fit$x <- fit$cols <- NULL
    }
    fit$call <- call
    fit
}

## In memory: the model frame of the main terms plus the FE and cluster
## variables (so that rows with a missing value in any of them are
## dropped), the design of the main terms (with the frame's predvars; the
## design keeps the frame, FE variables included), response and weights.
dm_fe_base_memory <- function(pf, call, env, contrasts, clv) {
    fe_vars <- unique(unlist(pf$fe, use.names = FALSE))
    aug <- dm_fe_augment(pf$main, unique(c(fe_vars, clv)))
    environment(aug) <- environment(pf$main)
    mf <- dmfit_model_frame(call, c("weights", "na.action"), env, formula = aug)
    mt_aug <- attr(mf, "terms")
    labs <- attr(mt_aug, "term.labels")
    drop <- which(!labs %in% attr(stats::terms(pf$main), "term.labels"))
    mt <- if (!length(drop)) {
        mt_aug
    } else if (length(drop) < length(labs)) {
        stats::drop.terms(mt_aug, drop, keep.response = TRUE)
    } else {
        ## no main term left: a formula with the response and the intercept
        tt <- stats::terms(pf$main)
        attr(tt, "predvars") <- attr(mt_aug, "predvars")[1:2]
        attr(tt, "dataClasses") <- attr(mt_aug, "dataClasses")[1L]
        tt
    }
    environment(mt) <- environment(pf$main)
    y <- dmfit_response(mf)
    if (is.null(y)) stop("The formula has no response.", call. = FALSE)
    if (is.matrix(y)) stop("Multivariate responses are not supported.", call. = FALSE)
    w <- as.vector(stats::model.weights(mf))
    if (!is.null(w) && !is.numeric(w)) stop("'weights' must be a numeric vector", call. = FALSE)
    if (!is.null(w)) w <- as.double(w)
    D <- design_matrix(mt, mf, contrasts.arg = contrasts)
    list(x = D, z = y, w = w, keys = dm_fe_keys(D, pf, clv, mf),
         rows = seq_len(nrow(mf)), nall = nrow(mf), na_act = attr(mf, "na.action"),
         terms = mt, contrasts = D$contrasts, xlevels = stats::.getXlevels(mt, mf), spec = dmfit_spec(D))
}

## Remotely: the remote design of the right-hand side, and the response and
## weights resolved as columns (a column name, DBI::SQL(), or an R
## expression: in SQL when translatable, else evaluated in R), on the rows
## where the model, FE and cluster variables, the response and the weights
## are present, as the model frame in memory (dm_remote_fit()).
dm_fe_base_remote <- function(pf, source, call, env, contrasts, dots, clv) {
    na <- if (is.null(call$na.action)) "omit" else dm_na_kind(eval(call$na.action, env), c("omit", "exclude"), "For remote data, ")
    nms <- names(dots)
    if (is.null(nms)) nms <- rep("", length(dots))
    keys <- c("table", "sql", "engine", "chunk_size")
    if (any(!nms %in% keys))
        stop("Unknown arguments: ", paste(nms[!nms %in% keys], collapse = ", "), call. = FALSE)
    if (length(attr(stats::terms(pf$main), "offset")))
        stop("Offsets are not supported with fixed effects.", call. = FALSE)
    dm_need_duckdb()
    src <- dm_remote_source(source, dots$table, dots$sql)
    r <- dm_remote_fit(pf$main, src, list(y = pf$main[[2L]], weights = call$weights), env, dots, contrasts,
                       present = unique(c(unlist(pf$fe, use.names = FALSE), clv)))
    D <- r$x
    list(x = D, z = r$cols$y, w = r$cols$weights, keys = dm_fe_keys(D, pf, clv), rows = NULL, nall = D$n,
         na_act = NULL, exclude = na == "exclude", terms = stats::terms(pf$main), contrasts = D$plan$contrasts, xlevels = D$plan$xlevels,
         spec = D$plan)
}

## The keys of the FE (and of the clusters, when they are not an FE): the
## values of their variables over the rows of x (from `data` in memory).
dm_fe_keys <- function(x, pf, clv, data = NULL) {
    list(fe = lapply(pf$fe, function(v) dmc_key(x, v, data)),
         cl = if (!is.null(clv) && is.na(match(clv, names(pf$fe)))) dmc_key(x, clv, data))
}

## The rows used and the FE codes, on the primitives (one code path): rows
## with a missing or non-finite response or weight, a zero weight or a
## missing FE or cluster value are dropped; then, with fixef.rm =
## "singleton", rows alone in a level of some FE, iteratively (grouped
## counts, a lookup of each row's count, dmp_subset()), until none is left.
## Every column is computed on the design of the rows kept. Returns
## list(x, cols = list(z, w, codes, cl), levels, nsingle, rows).
dm_fe_prepare <- function(base, pf, clv, fixef.rm) {
    x <- base$x
    z <- base$z
    w <- base$w
    rows <- base$rows
    keys <- base$keys$fe
    clkey <- base$keys$cl
    if (!is.null(w) && dm_colsum(dmc_is_finite(w) & w < 0) > 0)
        stop("missing or negative weights not allowed", call. = FALSE)
    restrict <- function(keep) {
        ## (keep is never missing: the logical columns are NA-safe)
        nd <- dm_colsum(!keep)
        if (nd == 0) return(0)
        x2 <- dmp_subset(x, keep)
        z <<- dmc_subset(z, keep, x2)
        w <<- dmc_subset(w, keep, x2)
        rows <<- dmc_subset(rows, keep, x2)
        keys <<- lapply(keys, dmc_subset, keep, x2)
        clkey <<- dmc_subset(clkey, keep, x2)
        x <<- x2
        if (dm_nrows(x) == 0) stop("0 (non-NA) cases", call. = FALSE)
        nd
    }
    keep <- dmc_is_finite(z)
    if (!is.null(w)) keep <- keep & dmc_is_finite(w) & w != 0
    for (k in c(keys, list(clkey))) if (!is.null(k)) keep <- keep & dmc_not_na(k)
    restrict(keep)
    nsingle <- 0
    repeat {
        last <- list()
        cds <- vector("list", length(keys))
        for (k in seq_along(keys)) {
            cds[[k]] <- dmc_dense_codes(keys[[k]], along = last)
            last <- list(cds[[k]]$codes$codes)
        }
        if (fixef.rm == "none") break
        one <- dmc_const(z, 1)
        cnt <- lapply(cds, function(cd) dmp_groupsum(one, cd$codes)[, 1L])
        ## no level with one row: no singleton
        if (!any(vapply(cnt, function(v) any(v == 1), NA))) break
        single <- NULL
        for (k in seq_along(keys)) {
            s <- dmc_lookup(cnt[[k]], cds[[k]]$codes, along = last) == 1
            single <- if (is.null(single)) s else single | s
            last <- list(single)
        }
        nd <- restrict(!single)
        if (nd == 0) break
        nsingle <- nsingle + nd
    }
    cl <- NULL
    if (!is.null(clv)) {
        m <- match(clv, names(pf$fe))
        if (is.na(m)) {
            cc <- dmc_dense_codes(clkey, along = last)
            last <- list(cc$codes$codes)
            cl <- cc$codes
        } else {
            cl <- cds[[m]]$codes
        }
        if (cl$G < 2L) stop("Clustered standard errors need at least two clusters.", call. = FALSE)
    }
    ## the design over the frame of the last columns (which keeps the others)
    x <- dm_rebase(x, if (length(last)) dm_rcols_frame(last[[1L]]))
    list(x = x, cols = list(z = z, w = w, codes = lapply(cds, `[[`, "codes"), cl = cl),
         levels = lapply(cds, `[[`, "levels"), nsingle = as.integer(nsingle), rows = rows)
}

## The fit: on the design with numeric variables centered when that is
## exact (dmp_center()), else on the design.
dm_fe_estimate <- function(x, cols, tol, cluster) {
    opts <- dm_fe_opts()
    ctr <- if (isTRUE(getOption("design.matrix.fe_center", TRUE))) dmp_center(x) else NULL
    if (!is.null(ctr)) {
        xc <- ctr$x
        on <- function(g) if (is.null(g)) NULL else dmc_codes(dmc_on(g$codes, xc), g$G)
        cc <- list(z = dmc_on(cols$z, xc), w = dmc_on(cols$w, xc), codes = lapply(cols$codes, on),
                   cl = on(cols$cl))
        est <- dm_fe_core(xc, cc, tol, opts, cluster)
        est <- dm_fe_uncenter(est, ctr$T, which(colnames(xc) == "(Intercept)"))
        if (!is.null(est)) return(est)
    }
    dm_fe_core(x, cols, tol, opts, cluster)
}

## The na.action of the fit: rows dropped for missing values (those of the
## model frame's na.action) and rows of the frame not used (zero weights,
## singletons), as indices in the data, with the class of the na.action
## ("omit" by default).
dm_fe_na_action <- function(na_act, rows, nall) {
    if (is.null(rows) || length(rows) == nall) return(na_act)
    orig <- if (is.null(na_act)) seq_len(nall) else seq_len(nall + length(na_act))[-na_act]
    dropped <- sort(c(as.integer(na_act), orig[-rows]))
    names(dropped) <- as.character(dropped)
    class(dropped) <- if (is.null(na_act)) "omit" else class(na_act)
    dropped
}

## The design_lm object of an FE fit.
dm_fe_finish <- function(est, cn, pf, lv, clv, nsingle, nall) {
    nm <- cn[est$xcols]
    b <- est$b
    names(b) <- nm
    sigma2 <- est$rss / est$df.residual
    V <- est$bread * sigma2
    dimnames(V) <- list(nm, nm)
    cu <- est$bread
    dimnames(cu) <- list(nm, nm)
    fe <- est$fe
    names(fe) <- names(pf$fe)
    names(lv) <- names(pf$fe)
    for (k in seq_along(fe)) names(fe[[k]]) <- dm_fe_labels(lv[[k]])
    G <- est$G
    names(G) <- names(pf$fe)
    fit <- list(coefficients = b, vcov = V, vcov_type = "iid",
                sigma = sqrt(sigma2), rss = est$rss, df.residual = est$df.residual,
                nobs = est$nobs, rank = est$rank, fixef = fe, fe_levels = lv,
                fe_vars = pf$fe, fe_sizes = G, fe_dof = est$fe_dof,
                fe_dof_exact = est$fe_dof_exact, nsingletons = nsingle,
                cov.unscaled = cu, tss = est$tss, tss_within = est$tss_within,
                compressed = FALSE, patterns = est$nobs, rows = nall, fe_info = est$info,
                intercept = TRUE)
    if (!is.null(clv)) {
        Vc <- est$vcov_cluster
        dimnames(Vc) <- list(nm, nm)
        fit$vcov <- Vc
        fit$vcov_type <- "cluster"
        fit$cluster <- clv
        fit$G_cluster <- est$G_cluster
    }
    class(fit) <- c("design_lm_fe", "design_lm")
    fit
}

## ---------------------------------------------------------------------------
## Methods

#' @export
vcov.design_lm_fe <- function(object, ...) object$vcov

#' @export
logLik.design_lm_fe <- function(object, ...) {
    N <- object$nobs
    val <- 0.5 * (object$sumlogw - N * (log(2 * pi) + 1 - log(N) + log(object$rss)))
    structure(val, nall = N, nobs = N, df = object$rank + object$fe_dof + 1L, class = "logLik")
}

#' @export
print.design_lm_fe <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    if (!is.null(x$call)) cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
    cat("\nFixed effects: ", paste(sprintf("%s (%s)", names(x$fe_sizes),
                                            dm_fmt_count(x$fe_sizes)), collapse = ", "), "\n", sep = "")
    dmfit_print_rows(x, prefix = "")
    dmfit_print_coef(x, digits)
    cat("\n")
    invisible(x)
}

#' @export
summary.design_lm_fe <- function(object, ...) {
    est <- object$coefficients
    aliased <- is.na(est)
    keep <- !aliased
    se <- sqrt(diag(object$vcov)[keep])
    tval <- est[keep] / se
    ## clustered: t(G - 1), as fixest
    tdf <- if (identical(object$vcov_type, "cluster")) object$G_cluster - 1L else object$df.residual
    coefs <- cbind(Estimate = est[keep], "Std. Error" = se, "t value" = tval,
                   "Pr(>|t|)" = 2 * stats::pt(abs(tval), tdf, lower.tail = FALSE))
    ans <- list(call = object$call, terms = object$terms, coefficients = coefs,
                aliased = aliased, sigma = object$sigma,
                df = c(object$rank, object$df.residual, length(est)),
                r.squared = 1 - object$rss / object$tss,
                within.r.squared = if (object$tss_within > 0) 1 - object$rss / object$tss_within else NA_real_,
                vcov_type = object$vcov_type, cluster = object$cluster, G_cluster = object$G_cluster,
                fe_sizes = object$fe_sizes, nobs = object$nobs, nsingletons = object$nsingletons,
                compressed = FALSE, method = object$method,
                rows = object$rows)
    n <- object$nobs
    ans$adj.r.squared <- 1 - (1 - ans$r.squared) * (n - 1) / object$df.residual
    class(ans) <- c("summary.design_lm_fe", "summary.design_lm")
    ans
}

#' @export
print.summary.design_lm_fe <- function(x, digits = max(3L, getOption("digits") - 3L),
                                       signif.stars = getOption("show.signif.stars"), ...) {
    if (!is.null(x$call)) cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
    cat("\nObservations: ", dm_fmt_count(x$nobs),
        if (x$nsingletons) sprintf(" (%s singletons removed)", dm_fmt_count(x$nsingletons)), "\n", sep = "")
    cat("Fixed effects: ", paste(sprintf("%s (%s)", names(x$fe_sizes), dm_fmt_count(x$fe_sizes)),
                                 collapse = ", "), "\n", sep = "")
    cat("Standard errors: ", if (identical(x$vcov_type, "cluster"))
        sprintf("clustered by %s (%s clusters)", x$cluster, dm_fmt_count(x$G_cluster)) else "iid",
        "\n", sep = "")
    dmfit_print_coefmat(x, digits, signif.stars, ...)
    cat("\nResidual standard error:", format(signif(x$sigma, digits)),
        "on", x$df[2L], "degrees of freedom\n")
    cat("Multiple R-squared: ", formatC(x$r.squared, digits = digits),
        ",\tWithin R-squared: ", formatC(x$within.r.squared, digits = digits), "\n\n", sep = "")
    invisible(x)
}

## The fitted values and residuals of the rows used, as columns.
dm_fe_columns <- function(object) {
    x <- object$x
    if (is.null(x))
        stop("This needs the fit's data: refit with `model = TRUE`, or use `predict(fit, newdata)`.",
             call. = FALSE)
    cols <- object$cols
    f <- dm_fe_fitted(x, dm_fe_beta(object, x), object$fixef, cols$codes)
    list(fitted = f, resid = cols$z - f)
}

## Coefficients for the columns of design x (0 for the intercept, absorbed,
## and for aliased columns).
dm_fe_beta <- function(object, x) {
    b <- object$coefficients
    b[is.na(b)] <- 0
    beta <- numeric(ncol(x))
    beta[match(names(b), colnames(x))] <- b
    beta
}

#' @export
residuals.design_lm_fe <- function(object, ...) {
    ## remote: the remote design with the columns .fitted and .resid
    r <- dm_fe_columns(object)
    dmp_output(object$x, list(.fitted = r$fitted, .resid = r$resid), design = TRUE,
               memory = function(d) stats::naresid(object$na.action, d$.resid))
}

#' @export
fitted.design_lm_fe <- function(object, ...) {
    ## remote: a table with .fitted and .resid
    r <- dm_fe_columns(object)
    dmp_output(object$x, list(.fitted = r$fitted, .resid = r$resid),
               memory = function(d) stats::napredict(object$na.action, d$.fitted))
}

#' @export
predict.design_lm_fe <- function(object, newdata, se.fit = FALSE,
                                 interval = c("none", "confidence", "prediction"),
                                 level = 0.95, na.action = stats::na.pass, ...) {
    interval <- match.arg(interval)
    if (isTRUE(se.fit) || interval != "none")
        stop("With fixed effects, predict() gives the fit only (no `se.fit` or intervals).", call. = FALSE)
    if (missing(newdata) || is.null(newdata)) return(fitted(object))
    if (!identical(na.action, stats::na.pass))
        stop("With fixed effects, predict() keeps every row of `newdata` (na.action = na.pass).", call. = FALSE)
    ## the design of the new rows (as predict() without FE)
    Dn <- dmfit_newdata(object, newdata, stats::na.pass, list(...))$D
    keys <- lapply(object$fe_vars, function(vars) {
        k <- tryCatch(dmc_key(Dn, vars, data = if (is.data.frame(newdata)) newdata),
                      error = function(e) NULL)
        if (is.null(k))
            stop(sprintf("`newdata` has no fixed-effect variable among `%s`.", paste(vars, collapse = "`, `")),
                 call. = FALSE)
        k
    })
    f <- dm_fe_fitted(Dn, dm_fe_beta(object, Dn), object$fixef, keys, object$fe_levels)
    ## a table remotely, a named vector in memory
    dmp_output(Dn, list(.fitted = f), memory = function(d) stats::setNames(d$.fitted, dimnames(Dn)[[1L]]))
}

## fixef(): the FE values, a list with one named vector per FE (fixest's
## normalization differs: compare sums over FE, or differences within a
## connected set of levels). The method is registered for fixest's and
## nlme's fixef() generics when those packages are loaded, without a
## dependency on them; the values are also `fit$fixef`.
#' @rawNamespace S3method(fixest::fixef, design_lm_fe)
#' @rawNamespace S3method(nlme::fixef, design_lm_fe)
fixef.design_lm_fe <- function(object, ...) object$fixef


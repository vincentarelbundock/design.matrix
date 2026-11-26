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
## levels are removed), so the reduced system is a symmetric M-matrix. It is
## solved densely when small (at most 3000 levels), otherwise by conjugate
## gradients preconditioned by the incomplete Cholesky factor IC(0) of R
## (two FE) or of its block diagonal (three or more: block Jacobi over the
## FEs, each block a Laplacian). Sparse Cholesky (CHOLMOD) is used only on
## request (option design.matrix.fe_solver = "chol"): on worker x firm
## graphs its fill is close to dense.
##
## A = K^- B gives the FE coefficients without another solve: a = A_y - A_x b.
## The coefficients are refined from the residuals (computed row by row):
## the normal equations X'M e = 0 correct the rounding of the one-pass
## aggregates. Numeric variables far from zero are centered first when that
## is an exact linear transform of the design (as for remote fits, see
## dm_center() in remote_lm.R), which keeps X'MX, hence the covariance
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
## or "pcg"), design.matrix.fe_precond ("ichol" or "jacobi") and
## design.matrix.fe_refine (steps of iterative refinement, default 1).
dm_fe_opts <- function() {
    solver <- getOption("design.matrix.fe_solver", "auto")
    if (!checkmate::test_choice(solver, c("auto", "dense", "chol", "pcg")))
        stop("Option `design.matrix.fe_solver` must be \"auto\", \"dense\", \"chol\" or \"pcg\".", call. = FALSE)
    precond <- getOption("design.matrix.fe_precond", "ichol")
    if (!checkmate::test_choice(precond, c("ichol", "jacobi")))
        stop("Option `design.matrix.fe_precond` must be \"ichol\" or \"jacobi\".", call. = FALSE)
    refine <- getOption("design.matrix.fe_refine", 1L)
    if (!checkmate::test_count(refine))
        stop("Option `design.matrix.fe_refine` must be a count.", call. = FALSE)
    list(solver = solver, precond = precond, pcg_tol = 1e-12, pcg_maxit = 10000L,
         refine = as.integer(refine))
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

## Grouped sums of the columns of a numeric matrix.
dm_fe_gsum <- function(codes, G, W) {
    if (!is.matrix(W)) W <- matrix(W, ncol = 1L)
    storage.mode(W) <- "double"
    .Call(C_dm_fe_gsum, codes, as.integer(G), W)
}

## The weighted crosstab of codes a and b as t(N_ab), a G_b x G_a dgCMatrix.
dm_fe_crosstab_t <- function(a, b, Ga, Gb, w = NULL) {
    ct <- .Call(C_dm_fe_crosstab, a, b, if (is.null(w)) NULL else as.double(w),
                as.integer(Ga), as.integer(Gb))
    methods::new("dgCMatrix", i = ct[[2L]], p = ct[[1L]], x = ct[[3L]],
                 Dim = c(as.integer(Gb), as.integer(Ga)))
}

## Component of each level of b in the bipartite graph of t(N_ab).
dm_fe_components <- function(Nt) {
    i <- rep.int(seq_len(ncol(Nt)) - 1L, diff(Nt@p))
    .Call(C_dm_fe_components, i, Nt@i, ncol(Nt), nrow(Nt))
}

## ---------------------------------------------------------------------------
## The reduced system

## Conjugate gradients for R Z = B, preconditioned by IC(0) of the block
## diagonal of R (blocks: the FEs; the whole of R with two FE), or by its
## diagonal when IC(0) breaks down or with precond = "jacobi".
dm_fe_pcg <- function(R, B, tol, maxit, block, precond = "ichol", cache = new.env()) {
    if (is.null(cache$A)) {
        cache$A <- methods::as(methods::as(R, "generalMatrix"), "CsparseMatrix")
        cache$precond <- precond
        if (precond == "ichol") {
            P <- Matrix::tril(cache$A)
            if (length(unique(block)) > 1L) {
                Ts <- Matrix::summary(P)
                same <- block[Ts$i] == block[Ts$j]
                P <- Matrix::sparseMatrix(i = Ts$i[same], j = Ts$j[same], x = Ts$x[same],
                                          dims = dim(P), repr = "C")
            }
            P <- methods::as(methods::as(P, "generalMatrix"), "CsparseMatrix")
            ## plain IC(0): the modified factor (MIC, omega > 0) took 2.5 to 5
            ## times as many iterations on worker x firm graphs
            Lx <- .Call(C_dm_fe_ichol, P@p, P@i, P@x, nrow(P), 0)
            if (!is.null(Lx)) {
                cache$L <- list(p = P@p, i = P@i, x = Lx)
            } else {
                cache$precond <- "jacobi"
            }
        }
    }
    A <- cache$A
    L <- cache$L
    storage.mode(B) <- "double"
    r <- .Call(C_dm_fe_pcg, A@p, A@i, A@x, L$p, L$i, L$x, B, as.double(tol), as.integer(maxit))
    list(x = r[[1L]], iter = r[[2L]], resid = r[[3L]], converged = r[[3L]] <= tol,
         precond = cache$precond)
}

## Solve the reduced system R Z = B (R symmetric positive definite after the
## reference levels are removed, or semidefinite with a consistent
## right-hand side). Returns Z and solver details.
dm_fe_solve <- function(R, B, opts, block, cache = new.env(), dense_max = 3000L) {
    m <- nrow(R)
    if (m == 0L) return(list(z = B, solver = "none"))
    solver <- opts$solver
    if (solver == "auto") solver <- if (m <= dense_max) "dense" else "pcg"
    if (solver == "dense") {
        ch <- tryCatch(chol(as.matrix(R)), error = function(e) NULL)
        if (!is.null(ch)) {
            z <- backsolve(ch, backsolve(ch, B, transpose = TRUE))
            return(list(z = z, solver = "dense"))
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
    sol <- dm_fe_pcg(R, B, opts$pcg_tol, opts$pcg_maxit, block, opts$precond, cache)
    if (!sol$converged)
        warning(sprintf("Conjugate gradients did not converge in %d iterations (relative residual %.2g).",
                        sol$iter, sol$resid), call. = FALSE)
    list(z = sol$x, solver = "pcg", precond = sol$precond, iter = sol$iter, resid = sol$resid)
}

## The FE system K = F'WF, reduced by eliminating the largest FE.
## G: levels per FE; nk: weighted counts per level; crosstab(a, b): t(N_ab).
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
    ## R = blockdiag(n_k) + offdiag(N_kl) - N_1' D_1^-1 N_1
    R <- Matrix::Diagonal(x = unlist(nk[rest], use.names = FALSE))
    if (length(rest) > 1L) {
        blocks <- list()
        for (r in seq_along(rest)[-length(rest)]) for (s in (r + 1L):length(rest)) {
            Nrs <- crosstab(rest[r], rest[s]) # G_s x G_r
            Ts <- Matrix::summary(Nrs)
            blocks[[length(blocks) + 1L]] <- data.frame(i = Ts$i + off[s], j = Ts$j + off[r], x = Ts$x)
        }
        bl <- do.call(rbind, blocks)
        Off <- Matrix::sparseMatrix(i = bl$i, j = bl$j, x = bl$x, dims = c(m, m))
        R <- R + Off + Matrix::t(Off)
    }
    R <- R - Matrix::tcrossprod(sys$Ntr %*% Matrix::Diagonal(x = 1 / sqrt(sys$D1)))
    R <- Matrix::forceSymmetric(methods::as(R, "CsparseMatrix"), uplo = "U")
    sys$keep <- setdiff(seq_len(m), ref)
    sys$cache <- new.env(parent = emptyenv()) # the preconditioner, built once
    sys$R <- R[sys$keep, sys$keep]
    sys$block <- rep(seq_along(rest), sys$Gr)[sys$keep]
    sys$m <- m
    sys$info <- list(m = length(sys$keep), nnz_R = length(R@x), ncomp = ncomp)
    sys
}

## A = K^- B (split by FE) and Q = B' K^- B, for B a list of G_k x q
## matrices (grouped sums by FE).
dm_fe_apply <- function(sys, B) {
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
    sol <- dm_fe_solve(sys$R, Brr[sys$keep, , drop = FALSE], sys$opts, sys$block, sys$cache)
    Z <- matrix(0, sys$m, ncol(Brr))
    Z[sys$keep, ] <- sol$z
    sol$z <- NULL
    Q <- Q + crossprod(Brr, Z)
    A[[f1]] <- (B[[f1]] - as.matrix(Matrix::crossprod(sys$Ntr, Z))) / sys$D1
    for (r in seq_along(rest)) A[[rest[r]]] <- Z[sys$off[r] + seq_len(sys$Gr[r]), , drop = FALSE]
    list(Q = Q, A = A, info = c(sol, sys$info))
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
## is a primitive of the fitting core (primitives.R, primitives_lookup.R)
## or one of the three sums below, whose in-memory methods use the C
## kernels of src/fe.c and whose remote methods are primitives:
##
##   dm_fe_groupsums(x, w, z, codes, G)  per FE: weighted counts, F_k'WX
##                                       and F_k'Wz (in memory: one pass
##                                       for every FE; remotely
##                                       dmp_colsums(group =) and
##                                       dmp_groupsum())
##   dm_fe_groupsum(col, code, G)        F_k'col (dmp_groupsum())
##   dm_fe_crosstab(a, b, w, Ga, Gb)     t(N_ab), a G_b x G_a dgCMatrix
##                                       (dmp_crosstab())
##
## Arguments of dm_fe_core(): the design x, and `cols`, the columns of the
## problem for that design: z (response), w (weights or NULL), codes (per
## FE, integer codes 1..G_k), G, and for clustered standard errors cl (the
## cluster codes 1..Gc) and Gc.

dm_fe_core <- function(x, cols, tol, opts, cluster = FALSE) {
    z <- cols$z
    w <- cols$w
    codes <- cols$codes
    G <- cols$G
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
    gs <- dm_fe_groupsums(x, w, z, codes, G)
    nk <- lapply(gs, `[[`, "nk")
    if (any(vapply(nk, function(v) any(!(v > 0)), NA)))
        stop("internal error: a fixed-effect level without weight", call. = FALSE)
    B <- lapply(gs, function(g) cbind(g$S[, xcols, drop = FALSE], g$s))
    crosstab <- function(a, b) dm_fe_crosstab(codes[[a]], codes[[b]], w, G[a], G[b])
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
             Be = lapply(seq_len(nfe), function(k) dm_fe_groupsum(we, codes[[k]], G[k])),
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
    ## when the reduced system was positive definite once references were
    ## removed (a dense or sparse Cholesky succeeded); otherwise NA
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
        Gc <- cols$Gc
        U <- dm_fe_groupsums(x, r$we, NULL, list(cl), Gc)[[1L]]$S[, xcols, drop = FALSE]
        ## an FE is nested in the clusters when each of its levels is in one
        ## cluster; otherwise its scores need a correction
        nested <- vapply(seq_len(nfe), function(k) dm_fe_nested(cl, codes[[k]], Gc, G[k]), NA)
        for (k in which(!nested)) {
            E <- dm_fe_crosstab(cl, codes[[k]], r$we, Gc, G[k])
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
    if (!dm_is_rcols(f)) f <- drop(f)
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

## The fit of a problem (see dm_fe_lm()): on the centered design when that
## is exact, else on the design.
dm_fe_estimate <- function(prob, tol, cluster) {
    opts <- dm_fe_opts()
    if (!is.null(prob$ctr)) {
        xc <- prob$x_centered
        est <- dm_fe_core(xc, prob$cols(xc), tol, opts, cluster)
        est <- dm_fe_uncenter(est, prob$ctr$T, which(colnames(xc) == "(Intercept)"))
        if (!is.null(est)) return(est)
    }
    dm_fe_core(prob$x, prob$cols(prob$x), tol, opts, cluster)
}

## ---------------------------------------------------------------------------
## The three sums with in-memory kernels (see above)

dm_fe_groupsums <- function(x, w, z, codes, G) {
    if (!inherits(x, "dm_remote")) return(dm_fe_gsums(x, w, z, codes, G))
    one <- if (is.null(w)) dmc_const(codes[[1L]], 1) else w
    xs <- if (is.null(w)) x else dmp_scale(x, w)
    lapply(seq_along(codes), function(k) {
        S <- dm_fe_place(dmp_colsums(xs, group = codes[[k]]), G[k])
        ns <- dm_fe_place(dmp_groupsum(if (is.null(z)) one else dmc_cbind(one, one * z), codes[[k]]), G[k])
        list(nk = ns[, 1L], S = S, s = if (is.null(z)) NULL else ns[, 2L])
    })
}

dm_fe_groupsum <- function(col, code, G) {
    if (!dm_is_rcols(col)) {
        return(.Call(C_dm_fe_struct_gsum, NULL, NULL, as.double(col), list(code), as.integer(G),
                     NULL, FALSE)[[1L]][1L, ])
    }
    dm_fe_place(dmp_groupsum(col, code), G)[, 1L]
}

dm_fe_crosstab <- function(a, b, w, Ga, Gb) {
    if (!dm_is_rcols(a) && !dm_is_rcols(b) && !dm_is_rcols(w))
        return(dm_fe_crosstab_t(a, b, Ga, Gb, w))
    ct <- dmp_crosstab(a, b, w)
    st <- Matrix::summary(ct)
    ia <- as.integer(rownames(ct))[st$i]
    ib <- as.integer(colnames(ct))[st$j]
    Matrix::sparseMatrix(i = ib, j = ia, x = st$x, dims = c(as.integer(Gb), as.integer(Ga)), repr = "C")
}

## Is code b nested in code a (each level of b in one level of a)?
dm_fe_nested <- function(a, b, Ga, Gb) {
    if (identical(a, b)) return(TRUE)
    if (!dm_is_rcols(a) && !dm_is_rcols(b)) return(.Call(C_dm_fe_nested, a, b, as.integer(Gb)))
    ct <- dm_fe_crosstab(a, b, NULL, Ga, Gb) # G_b x G_a
    all(tabulate(ct@i + 1L, nbins = Gb) == 1L)
}

## Rows of grouped sums labelled by codes 1..G, in code order.
dm_fe_place <- function(S, G) {
    out <- matrix(0, G, ncol(S))
    out[as.integer(rownames(S)), ] <- unname(S)
    out
}

## For each FE k: list(nk = sums of w by level, S = grouped sums of w X
## (G_k x p), s = grouped sums of w z (NULL when z is NULL)), in memory:
## one pass over the rows with the structured kernel when the design
## allows it, the generic row kernel otherwise.
dm_fe_gsums <- function(D, w, y, codes, G) {
    nfe <- length(codes)
    wd <- if (is.null(w)) NULL else as.double(w)
    p <- ncol(D)
    Y <- if (is.null(y)) NULL else matrix(as.double(y), ncol = 1L)
    if (dm_fe_struct_ok(D, G)) {
        M <- .Call(C_dm_fe_struct_gsum, D$setup, D$vars, wd, codes, as.integer(G), Y, TRUE)
        info <- dm_struct_info(D)
        icpt <- isTRUE(D$setup$intercept > 0L)
        return(lapply(seq_len(nfe), function(k) {
            m <- M[[k]]
            S <- matrix(0, G[k], p)
            q <- 1L
            for (t in seq_along(info)) {
                Tt <- info[[t]]
                if (t == 1L && icpt) {
                    S[, Tt$cols[1L, 1L]] <- m[1L, ]
                    next
                }
                L <- nrow(Tt$C)
                ns <- nrow(Tt$cols)
                for (u in seq_len(ns)) {
                    rows <- q + (u - 1L) * L + seq_len(L)
                    S[, Tt$cols[u, ]] <- crossprod(m[rows, , drop = FALSE], Tt$C)
                }
                q <- q + L * ns
            }
            list(nk = m[1L, ], S = S, s = if (is.null(Y)) NULL else m[q + 1L, ])
        }))
    }
    one <- if (is.null(wd)) rep(1, nrow(D)) else wd
    lapply(seq_len(nfe), function(k) {
        list(nk = dm_fe_gsum(codes[[k]], G[k], one)[, 1L],
             S = dm_fe_dsum(D, wd, codes[[k]], G[k]),
             s = if (is.null(Y)) NULL else dm_fe_gsum(codes[[k]], G[k], one * Y)[, 1L])
    })
}

## Can the structured kernel compute the grouped sums? (A plain design that
## the structured kernels accept, and accumulators of at most 2^25 doubles
## per FE.)
dm_fe_struct_ok <- function(D, G) {
    if (!isTRUE(getOption("design.matrix.structured", TRUE))) return(FALSE)
    if (!dm_struct_ok(D)) return(FALSE)
    info <- dm_struct_info(D)
    Q <- 2 + sum(vapply(info, function(Tt) nrow(Tt$C) * nrow(Tt$cols), 0))
    all(Q * as.double(G) <= 2^25)
}

## colSums(w * D, group = codes) for codes 1..G: a G x p matrix (row kernel).
dm_fe_dsum <- function(D, w, codes, G) {
    J <- as_jac(D)
    W <- matrix(if (is.null(w)) 1 else as.double(w), nrow = nrow(D), ncol = 1L)
    if (jac_kernel_ok(J)) return(t(dm_crossprod_kernel(J, W, codes, G)))
    X <- as.matrix(D)
    storage.mode(X) <- "double"
    dm_fe_gsum(codes, G, W[, 1L] * X)
}

## Centering of numeric variables (in memory): list(shift, T), with the
## design of the shifted frame X~ = X T and T unit upper triangular, or
## NULL. As dm_center() does for remote designs: only variables whose mean
## is large compared with their spread, and only when such a T exists
## (checked on prototype rows).
dm_fe_center <- function(mt, mf, contrasts, w = NULL) {
    if (!isTRUE(getOption("design.matrix.fe_center", TRUE))) return(NULL)
    if (attr(mt, "intercept") == 0L) return(NULL)
    tv <- dm_terms_vars(mt)$names
    resp <- if (attr(mt, "response") > 0L) tv[attr(mt, "response")] else character()
    cand <- setdiff(intersect(tv, names(mf)), resp)
    num <- cand[vapply(cand, function(nm) identical(stats::.MFclass(mf[[nm]]), "numeric"), NA)]
    if (!length(num)) return(NULL)
    ## the decision and the shift from a systematic sample of rows (any
    ## shift is exact; a shift near the mean keeps the sums accurate)
    nr <- nrow(mf)
    idx <- if (nr > 1e5) unique(round(seq(1, nr, length.out = 1e5))) else seq_len(nr)
    mu <- vapply(num, function(nm) {
        x <- mf[[nm]][idx]
        if (is.null(w)) mean(x) else sum(w[idx] * x) / sum(w[idx])
    }, 0)
    sdv <- vapply(num, function(nm) stats::sd(mf[[nm]][idx]), 0)
    sdv[!is.finite(sdv)] <- 0
    big <- is.finite(mu) & abs(mu) > sdv
    if (!any(big)) return(NULL)
    shift <- mu[big]
    m0 <- mf[seq_len(min(nrow(mf), 2000L)), , drop = FALSE]
    attr(m0, "terms") <- NULL
    p <- tryCatch(ncol(design_matrix(mt, m0, contrasts.arg = contrasts)), error = function(e) NA)
    if (is.na(p) || p > 2000L) return(NULL)
    m <- max(4L * p, 200L)
    proto <- m0[sample.int(nrow(m0), m, replace = TRUE), , drop = FALSE]
    for (nm in names(proto)) {
        x <- proto[[nm]]
        if (is.factor(x)) {
            proto[[nm]] <- factor(sample(levels(x), m, TRUE), levels = levels(x))
            attr(proto[[nm]], "contrasts") <- attr(x, "contrasts")
        } else if (is.logical(x)) {
            proto[[nm]] <- sample(c(FALSE, TRUE), m, TRUE)
        } else if (nm %in% num) {
            proto[[nm]] <- stats::rnorm(m)
        }
    }
    shifted <- proto
    for (nm in names(shift)) shifted[[nm]] <- proto[[nm]] - shift[[nm]]
    X0 <- tryCatch(as.matrix(design_matrix(mt, proto, contrasts.arg = contrasts)), error = function(e) NULL)
    X1 <- tryCatch(as.matrix(design_matrix(mt, shifted, contrasts.arg = contrasts)), error = function(e) NULL)
    if (is.null(X0) || is.null(X1) || ncol(X0) != p || ncol(X1) != p) return(NULL)
    dm_center_T(X0, X1, shift)
}

## T with X1 = X0 T, unit upper triangular, or NULL (see dm_center()).
dm_center_T <- function(X0, X1, shift) {
    p <- ncol(X0)
    q <- qr(X0)
    if (q$rank < p) return(NULL)
    Tm <- unname(qr.coef(q, X1))
    if (!all(is.finite(Tm)) || max(abs(X0 %*% Tm - X1)) > 1e-8 * max(1, abs(X1))) return(NULL)
    ts <- max(1, abs(Tm))
    if (max(abs(diag(Tm) - 1)) > 1e-8 * ts || max(abs(Tm[lower.tri(Tm)])) > 1e-8 * ts) return(NULL)
    Tm[lower.tri(Tm)] <- 0
    diag(Tm) <- 1
    Tm[abs(Tm) < 1e-10 * ts] <- 0
    list(shift = shift, T = Tm)
}

## ---------------------------------------------------------------------------
## design_lm() with fixed effects
##
## Preparing the problem differs between data in memory and remote data,
## as it does without FE: the rows used (model frame, missing values,
## weights), the FE codes and the iterative removal of singletons (in
## memory: integer codes and a C loop before the design is built;
## remotely: a temporary table of the rows used with dense-rank codes, and
## DELETE ... HAVING COUNT(*) = 1 until no row goes). Both give a problem:
## list(x, x_centered, ctr, cols(x), levels, nsingle, nall, ...), and the
## estimator, the methods and predictions are one code path.

#' @section Fixed effects:
#' A formula `y ~ x | fe1 + fe2` absorbs fixed effects, as
#' [fixest::feols()]: each FE is a variable (a factor, character, integer or
#' numeric column) or an interaction `a^b`. The fit never builds the FE
#' dummies or demeans the rows: it solves the normal equations from sums
#' over rows (the design's products, sums within FE levels and crosstabs of
#' the FE), eliminating the largest FE exactly and solving for the others
#' densely (up to 3000 levels) or by conjugate gradients with an incomplete
#' Cholesky preconditioner. The same code fits data in memory and remote
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
#'   For remote data, `weights` names a column, and the response must have
#'   a SQL translation.
#' * Options: `design.matrix.fe_solver` (`"auto"`, `"dense"`, `"pcg"`, or
#'   `"chol"` for a sparse Cholesky factor), `design.matrix.fe_precond`
#'   (`"ichol"` or `"jacobi"`), `design.matrix.fe_refine` (refinement steps,
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
    prob <- if (!is.null(data) && dm_is_remote_source(data)) {
        if (!is.null(call$na.action) && dm_remote_na(eval(call$na.action, env)) != "omit")
            stop("For remote data, `na.action` must be na.omit.", call. = FALSE)
        w <- call$weights
        if (is.name(w)) w <- as.character(w)
        if (!is.null(w) && !checkmate::test_string(w))
            stop("For remote data, `weights` must name a column.", call. = FALSE)
        dm_fe_remote_problem(pf, data, w, contrasts, clv, fixef.rm, dots)
    } else {
        if (length(dots)) stop("Unknown arguments: ", paste(names(dots), collapse = ", "), call. = FALSE)
        dm_fe_memory_problem(pf, call, env, contrasts, clv, fixef.rm)
    }
    est <- dm_fe_estimate(prob, tol, !is.null(clv))
    fit <- dm_fe_finish(est, colnames(prob$x), pf, prob$levels, clv, prob$nsingle, prob$nall)
    fit$weighted <- prob$weighted
    fit$sumlogw <- prob$sumlogw
    fit$x <- prob$x
    fit$cols <- prob$cols(prob$x)
    fit$terms <- prob$terms
    fit$contrasts <- prob$contrasts
    fit$xlevels <- prob$xlevels
    fit$spec <- prob$spec
    fit$na.action <- prob$na.action
    if (inherits(prob$x, "dm_remote")) {
        fit$remote <- TRUE
        fit$remote_design <- prob$x
        fit$method <- if (prob$x$engine == "sql") "SQL aggregates" else "aggregates computed in batches"
    } else if (model) {
        fit$model <- prob$mf
    } else {
        fit$x <- fit$cols <- NULL
    }
    fit$call <- call
    fit
}

## The model frame of an FE fit: the main terms plus the FE and cluster
## variables (so that rows with a missing value in any of them are
## dropped), and the terms of the main formula (with the frame's predvars).
dm_fe_frame <- function(pf, call, env, clv) {
    fe_vars <- unique(unlist(pf$fe, use.names = FALSE))
    aug <- dm_fe_augment(pf$main, unique(c(fe_vars, clv)))
    environment(aug) <- environment(pf$main)
    mf <- call[c(1L, match(c("data", "weights", "na.action"), names(call), 0L))]
    mf$formula <- aug
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(design.matrix::design_frame)
    mf <- eval(mf, env)
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
    list(mf = mf, mt = mt, fe_vars = fe_vars)
}

dm_fe_memory_problem <- function(pf, call, env, contrasts, clv, fixef.rm) {
    fr <- dm_fe_frame(pf, call, env, clv)
    mf <- fr$mf
    mt <- fr$mt
    fe_vars <- fr$fe_vars
    na_act <- attr(mf, "na.action")
    y <- dmfit_response(mf)
    if (is.null(y)) stop("The formula has no response.", call. = FALSE)
    if (is.matrix(y)) stop("Multivariate responses are not supported.", call. = FALSE)
    w <- as.vector(stats::model.weights(mf))
    if (!is.null(w) && !is.numeric(w)) stop("'weights' must be a numeric vector", call. = FALSE)
    nall <- nrow(mf)
    dmfit_check(nall, y, w, NULL)
    if (!is.null(w)) w <- as.double(w)
    ## rows used: positive weights, no singletons
    rows <- if (!is.null(w) && any(w == 0)) which(w != 0) else NULL
    fe_codes <- function(r) {
        d <- if (is.null(r)) mf else mf[r, fe_vars, drop = FALSE]
        lapply(pf$fe, function(v) dm_fe_code_term(d, v))
    }
    cd <- fe_codes(rows)
    nsingle <- 0L
    if (fixef.rm == "singleton") {
        keep <- .Call(C_dm_fe_singletons, lapply(cd, `[[`, "codes"), vapply(cd, `[[`, 0L, "G"))
        if (!is.null(keep)) {
            nsingle <- sum(!keep)
            rows <- if (is.null(rows)) which(keep) else rows[keep]
            if (!length(rows)) stop("No observation left after removing singletons.", call. = FALSE)
            cd <- fe_codes(rows)
        }
    }
    if (!is.null(rows)) {
        tt <- attr(mf, "terms")
        mf <- mf[rows, , drop = FALSE]
        attr(mf, "terms") <- tt
        ## levels of the design's factors absent from the rows used
        for (v in intersect(dm_terms_vars(mt)$names, names(mf))) if (is.factor(mf[[v]])) {
            ct <- attr(mf[[v]], "contrasts")
            mf[[v]] <- droplevels(mf[[v]])
            if (!is.null(ct)) attr(mf[[v]], "contrasts") <- ct
        }
        y <- y[rows]
        if (!is.null(w)) w <- w[rows]
    }
    codes <- lapply(cd, `[[`, "codes")
    G <- vapply(cd, `[[`, 0L, "G")
    cl <- Gc <- NULL
    if (!is.null(clv)) {
        same <- match(clv, names(pf$fe))
        ccd <- if (is.na(same)) dm_fe_code(mf[[clv]]) else cd[[same]]
        cl <- ccd$codes
        Gc <- ccd$G
        if (Gc < 2L) stop("Clustered standard errors need at least two clusters.", call. = FALSE)
    }
    D <- design_matrix(mt, mf, contrasts.arg = contrasts)
    ctr <- dm_fe_center(mt, mf, contrasts, w)
    D2 <- NULL
    if (!is.null(ctr)) {
        mf2 <- mf
        for (nm in names(ctr$shift)) mf2[[nm]] <- mf[[nm]] - ctr$shift[[nm]]
        D2 <- design_matrix(mt, mf2, contrasts.arg = contrasts)
    }
    list(x = D, x_centered = D2, ctr = ctr,
         cols = function(x) list(z = y, w = w, codes = codes, G = G, cl = cl, Gc = Gc),
         levels = lapply(cd, `[[`, "levels"), nsingle = nsingle, nall = nall,
         weighted = !is.null(w), sumlogw = if (is.null(w)) 0 else sum(log(w)),
         terms = mt, contrasts = D$contrasts, xlevels = stats::.getXlevels(mt, mf),
         spec = dmfit_spec(D), na.action = dm_fe_na_action(na_act, rows, nall), mf = mf)
}

dm_fe_remote_problem <- function(pf, source, w, contrasts, clv, fixef.rm, dots) {
    dm_need_duckdb()
    nms <- names(dots)
    if (is.null(nms)) nms <- rep("", length(dots))
    keys <- c("table", "sql", "engine", "chunk_size")
    if (any(!nms %in% keys))
        stop("Unknown arguments: ", paste(nms[!nms %in% keys], collapse = ", "), call. = FALSE)
    if (length(attr(stats::terms(pf$main), "offset")))
        stop("Offsets are not supported with fixed effects.", call. = FALSE)
    ## the rows used by the main formula (and the weights)
    D0 <- suppressMessages(dm_remote(pf$main, source, contrasts, NULL, dots, weights = w, drop_unused = TRUE))
    if (is.null(D0$response)) stop("The formula has no response.", call. = FALSE)
    if (D0$vars[[D0$response]]$kind != "numeric")
        stop("The response of a remote fit must be one numeric variable with a SQL translation.", call. = FALSE)
    fe_vars <- unique(unlist(pf$fe, use.names = FALSE))
    miss <- setdiff(c(fe_vars, clv), names(D0$src$types))
    if (length(miss)) stop(sprintf("No column `%s` in the remote data.", miss[1L]), call. = FALSE)
    tab <- dm_fe_remote_table(D0, pf, fe_vars, w, clv, fixef.rm)
    rest <- dots[!nms %in% c("table", "sql")]
    D <- dm_remote(pf$main, dm_con(D0), contrasts, NULL, c(list(table = tab$name), rest),
                   weights = w, drop_unused = TRUE)
    ## the temporary table lives on D0's connection (a private one for Parquet)
    D$src$holder$parent <- D0$src$holder
    q <- function(x) dm_qi(D, x)
    nfe <- length(pf$fe)
    fecols <- paste0("__dm_fe", seq_len(nfe))
    clcol <- NULL
    if (!is.null(clv)) {
        m <- match(clv, names(pf$fe))
        clcol <- if (is.na(m)) "__dm_cl" else fecols[m]
    }
    gmax <- function(col) as.integer(dm_query(D, sprintf("SELECT MAX(%s) AS g FROM %s", q(col), q(tab$name)))$g)
    G <- vapply(fecols, gmax, 0L, USE.NAMES = FALSE)
    Gc <- NULL
    if (!is.null(clcol)) {
        Gc <- gmax(clcol)
        if (Gc < 2L) stop("Clustered standard errors need at least two clusters.", call. = FALSE)
    }
    levels <- lapply(seq_len(nfe), function(k) {
        vars <- pf$fe[[k]]
        r <- dm_query(D, sprintf("SELECT DISTINCT %s, %s FROM %s ORDER BY 1", q(fecols[k]),
                                 paste(q(vars), collapse = ", "), q(tab$name)))
        out <- r[vars]
        for (v in vars) if (is.factor(out[[v]])) out[[v]] <- as.character(out[[v]])
        out
    })
    ## centering (as dm_remote_lm_moments())
    ctr <- if (isTRUE(getOption("design.matrix.fe_center", TRUE)) &&
               isTRUE(getOption("design.matrix.remote_center", TRUE))) dm_center(D) else NULL
    D2 <- NULL
    if (!is.null(ctr)) {
        D2 <- D
        D2$cache <- new.env(parent = emptyenv())
        D2$id <- dm_tmp_name("dm")
        if (D$engine == "sql") {
            for (nm in names(ctr$shift))
                D2$vars[[nm]]$sql <- sprintf("(%s - %s)", D$vars[[nm]]$sql, dm_sql_double(ctr$shift[[nm]]))
        } else {
            D2$shift <- ctr$shift
        }
    }
    ysql <- D$vars[[D$response]]$sql
    wsql <- if (is.null(w)) NULL else sprintf("CAST(s.%s AS DOUBLE)", q(w))
    cols <- function(x) list(z = dmp_column(x, ysql), w = dmp_column(x, wsql),
                             codes = lapply(fecols, function(col) dmp_column(x, sprintf("s.%s", q(col)))),
                             G = G, cl = if (!is.null(clcol)) dmp_column(x, sprintf("s.%s", q(clcol))), Gc = Gc)
    sumlogw <- if (is.null(w)) 0 else as.double(dm_query(D, sprintf("SELECT SUM(ln(%s)) AS s FROM %s AS s",
                                                                       wsql, q(tab$name)))$s)
    list(x = D, x_centered = D2, ctr = ctr, cols = cols, levels = levels, nsingle = tab$nsingle,
         nall = D0$n, weighted = !is.null(w), sumlogw = sumlogw, terms = D$terms_full,
         contrasts = D$plan$contrasts, xlevels = D$plan$xlevels, spec = D$plan, na.action = NULL)
}

## The temporary table of the rows used: the columns the model reads, the
## FE codes __dm_fe<k> (dense ranks of the FE variables) and, for a cluster
## variable that is not an FE, its codes __dm_cl; singletons deleted
## iteratively, then the codes ranked again.
dm_fe_remote_table <- function(D0, pf, fe_vars, w, clv, fixef.rm) {
    q <- function(x) dm_qi(D0, x)
    cols <- unique(c(dm_remote_inputs(D0), fe_vars, clv, w))
    types <- D0$src$types
    present <- vapply(unique(c(fe_vars, clv)), function(col) {
        if (types[[col]] == "double") sprintf("s.%s IS NOT NULL AND NOT isnan(s.%s)", q(col), q(col))
        else sprintf("s.%s IS NOT NULL", q(col))
    }, "")
    where <- c(sprintf("(%s)", D0$filter), present,
               if (!is.null(w)) sprintf("CAST(s.%s AS DOUBLE) <> 0", q(w)))
    nfe <- length(pf$fe)
    fecols <- paste0("__dm_fe", seq_len(nfe))
    ranks <- vapply(seq_len(nfe), function(k) {
        sprintf("DENSE_RANK() OVER (ORDER BY %s) AS %s", paste0("s.", q(pf$fe[[k]]), collapse = ", "), q(fecols[k]))
    }, "")
    clrank <- if (!is.null(clv) && is.na(match(clv, names(pf$fe))))
        sprintf("DENSE_RANK() OVER (ORDER BY s.%s) AS %s", q(clv), q("__dm_cl"))
    name <- dm_tmp_name("dm_fe")
    dm_execute(D0, sprintf("CREATE TEMPORARY TABLE %s AS SELECT %s FROM (%s) AS s WHERE %s",
                           q(name), paste(c(paste0("s.", q(cols)), ranks, clrank), collapse = ", "),
                           D0$src$select, paste(where, collapse = " AND ")))
    nsingle <- 0
    if (fixef.rm == "singleton") {
        conds <- vapply(fecols, function(col) {
            sprintf("%s IN (SELECT %s FROM %s GROUP BY %s HAVING COUNT(*) = 1)", q(col), q(col), q(name), q(col))
        }, "")
        repeat {
            k <- dm_execute(D0, sprintf("DELETE FROM %s WHERE %s", q(name), paste(conds, collapse = " OR ")))
            if (!k) break
            nsingle <- nsingle + k
        }
        if (nsingle) {
            rr <- c(sprintf("DENSE_RANK() OVER (ORDER BY %s) AS %s", q(fecols), q(fecols)),
                    if (!is.null(clrank)) sprintf("DENSE_RANK() OVER (ORDER BY %s) AS %s", q("__dm_cl"), q("__dm_cl")))
            name2 <- dm_tmp_name("dm_fe")
            dm_execute(D0, sprintf("CREATE TEMPORARY TABLE %s AS SELECT * REPLACE (%s) FROM %s",
                                   q(name2), paste(rr, collapse = ", "), q(name)))
            dm_execute(D0, sprintf("DROP TABLE %s", q(name)))
            name <- name2
        }
    }
    n <- as.numeric(dm_query(D0, sprintf("SELECT COUNT(*) AS n FROM %s", q(name)))$n)
    if (n == 0) stop("No observation left after removing singletons.", call. = FALSE)
    list(name = name, nsingle = as.integer(nsingle))
}

## The na.action of the fit: rows dropped for missing values (those of the
## model frame's na.action) and rows of the frame not used (zero weights,
## singletons), as indices in the data, with the class of the na.action
## ("omit" by default).
dm_fe_na_action <- function(na_act, rows, nall) {
    if (is.null(rows)) return(na_act)
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
    if (isTRUE(x$remote)) cat(sprintf("Fitted from %s of %s rows.\n", x$method, dm_fmt_count(x$rows)))
    if (length(x$coefficients)) {
        cat("\nCoefficients:\n")
        print.default(format(x$coefficients, digits = digits), print.gap = 2L, quote = FALSE)
    } else {
        cat("\nNo coefficients\n")
    }
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
                compressed = FALSE, remote = object$remote, method = object$method,
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
    nsing <- sum(x$aliased)
    if (length(x$aliased) == 0L) {
        cat("\nNo Coefficients\n")
    } else {
        cat(if (nsing) sprintf("\nCoefficients: (%d not defined because of singularities)\n", nsing)
            else "\nCoefficients:\n")
        coefs <- matrix(NA_real_, length(x$aliased), 4L,
                        dimnames = list(names(x$aliased), colnames(x$coefficients)))
        coefs[!x$aliased, ] <- x$coefficients
        stats::printCoefmat(coefs, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
    }
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
    r <- dm_fe_columns(object)
    ## remote: the remote design with the columns .fitted and .resid
    if (isTRUE(object$remote)) return(dm_with_columns(object$x, dmc_cbind(r$fitted, r$resid), c(".fitted", ".resid")))
    stats::naresid(object$na.action, as.vector(r$resid))
}

#' @export
fitted.design_lm_fe <- function(object, ...) {
    r <- dm_fe_columns(object)
    ## remote: a table with .fitted and .resid
    if (isTRUE(object$remote)) return(dmp_output(object$x, list(.fitted = r$fitted, .resid = r$resid)))
    stats::napredict(object$na.action, as.vector(r$fitted))
}

#' @export
predict.design_lm_fe <- function(object, newdata, se.fit = FALSE,
                                 interval = c("none", "confidence", "prediction"),
                                 level = 0.95, na.action = stats::na.pass, ...) {
    interval <- match.arg(interval)
    if (isTRUE(se.fit) || interval != "none")
        stop("With fixed effects, predict() gives the fit only (no `se.fit` or intervals).", call. = FALSE)
    if (missing(newdata) || is.null(newdata)) return(fitted(object))
    if (is.null(object$spec)) stop("This design cannot be applied to new data.", call. = FALSE)
    if (!identical(na.action, stats::na.pass))
        stop("With fixed effects, predict() keeps every row of `newdata` (na.action = na.pass).", call. = FALSE)
    remote <- dm_is_remote_source(newdata)
    Dn <- if (remote) {
        do.call(design_matrix, c(list(object$spec, newdata, na.action = stats::na.pass), list(...)))
    } else {
        design_matrix(object$spec, newdata, na.action = stats::na.pass)
    }
    have <- if (remote) names(Dn$src$types) else names(newdata)
    keys <- lapply(object$fe_vars, function(vars) {
        miss <- setdiff(vars, have)
        if (length(miss))
            stop(sprintf("`newdata` has no fixed-effect variable `%s`.", miss[1L]), call. = FALSE)
        if (remote) dmp_column(Dn, sprintf("s.%s", dm_qi(Dn, vars))) else newdata[vars]
    })
    f <- dm_fe_fitted(Dn, dm_fe_beta(object, Dn), object$fixef, keys, object$fe_levels)
    if (remote) {
        if (".fitted" %in% names(Dn$src$types))
            stop("The remote data already has a column `.fitted`.", call. = FALSE)
        return(dmp_output(Dn, list(.fitted = f)))
    }
    f <- as.vector(f)
    names(f) <- dimnames(Dn)[[1L]]
    f
}

## fixef(): the FE values, a list with one named vector per FE (fixest's
## normalization differs: compare sums over FE, or differences within a
## connected set of levels). The method is registered for fixest's and
## nlme's fixef() generics when those packages are loaded, without a
## dependency on them; the values are also `fit$fixef`.
#' @rawNamespace S3method(fixest::fixef, design_lm_fe)
#' @rawNamespace S3method(nlme::fixef, design_lm_fe)
fixef.design_lm_fe <- function(object, ...) object$fixef


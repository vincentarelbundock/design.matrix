## Absorbed fixed effects: the estimator of design_lm(y ~ x | fe1 + fe2) and of
## the weighted least squares steps of design_glm() with fixed effects.
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
## (xtx_fe_core()) runs on data in memory (C kernels: one pass for all
## grouped sums and the small crosstabs, xtx_fe_struct_gsum() in src/fe.c,
## and one for the crosstabs of the largest FE with the others,
## xtx_fe_crosstabs(); on several threads with design_threads()) and on
## remote data (the remote methods of the same primitives): one code path,
## see xtx_fe_core().
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
## solver then drops a further level (xtx_fe_solve()).
##
## R is solved densely when small (at most 1000 levels: beyond, the dense
## factorization costs more than conjugate gradients), by a pivoted
## Cholesky factorization that sets the levels of pivots at rounding level
## to zero (and gives the rank of the FE), otherwise by
## conjugate gradients, on R formed or applied without forming it:
##
## * Formed (xtx_fe_schur() in src/fe.c, the lower triangle by Gustavson's
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
##   of R) exceeds form_max (option xtx.fe_form_max, default 25,
##   calibrated on worker x firm panels: the two paths take the same time
##   near 20 to 40) times nnz(N_1); a probe of 64 columns estimates both
##   first, so that giving up is cheap.
##
## The stopping rule is a normwise backward error below 1e-12 (see
## xtx_fe_pcg() in src/fe.c), with the norm of the right-hand side floored
## at that of the first solve of the system (a refinement step's right-hand
## side is at rounding level), and a warm start (Z0) is accepted, for
## sequences of nearby systems (the steps of a GLM). Sparse Cholesky
## (CHOLMOD) is used only on request (option xtx.fe_solver =
## "chol"): on worker x firm graphs its fill is close to dense.
##
## A = K^- B gives the FE coefficients without another solve: a = A_y - A_x b.
## The coefficients are refined from the residuals (computed row by row):
## the normal equations X'M e = 0 correct the rounding of the one-pass
## aggregates. Numeric variables far from zero are centered first when that
## is an exact linear transform of the design (as for remote fits, see
## xtx_center() in qr.R), which keeps X'MX, hence the covariance
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
## The estimator, on the primitives
##
## One implementation for designs in memory and remote designs: every sum
## is a primitive of the fitting core (primitives.R, primitives_lookup.R,
## primitives_fused.R): X'WX and X'Wz (xtx_p_crossprod), the moments of z
## (xtx_c_moments), F_k'WX, F_k'W1 and F_k'Wz for all FE and the crosstabs
## N_kl (xtx_p_fe_sums, two passes in memory), the residuals, F_k'We and
## the rss (xtx_p_mult and xtx_c_fe_resid, one pass in memory).
##
## Arguments of xtx_fe_core(): the design x, and `cols`, the columns of the
## problem for that design: z (response), w (weights or NULL), codes (per
## FE, a xtx_c_codes() grouping) and, for clustered standard errors, cl (a
## xtx_c_codes() grouping of the clusters). `warm`: NULL, or an environment
## that keeps the solution of the reduced system from one call to the next
## (the weighted least squares steps of a GLM, whose solutions change
## little: a warm start for the solver of the reduced system) and the
## numbers of iterations. With resid_final = FALSE, the residuals after the
## last refinement are not computed (no rss, `we`; for the steps of a GLM,
## which need them only for the last step: see xtx_fe_glm_estimate()).
## The fit keeps its absorption operator (`absorb`), for the cluster scores
## (xtx_fe_cluster()).
##
## Steps: project (the moments, the FE absorbed) -> pivoted Cholesky ->
## solve -> refine (from the residuals) -> degrees of freedom.

#' Options of the fixed-effects solver
#'
#' `absorb_control()` collects the options of the solver of fixed effects
#' (FE) that [design_absorb()] and [design_absorb_wls()] use: by default the
#' options `xtx.fe_solver` (`"auto"`, `"dense"`, `"chol"` or `"pcg"`),
#' `xtx.fe_precond` (`"auto"`, default: incomplete Cholesky when the reduced
#' system is formed, Jacobi when it is applied without being formed;
#' `"jacobi"`; or `"ichol"`, which always forms it), `xtx.fe_form_max` (when
#' to form the reduced system: default 25, `0` never, `Inf` always) and
#' `xtx.fe_refine` (refinement steps: `"auto"`, default, one step whose final
#' pass over the residuals is skipped when its correction is negligible; or a
#' count of steps, each followed by that pass).
#'
#' Names may change before the first release.
#'
#' @param solver,precond,form_max,refine See Details.
#' @return A list of options for [design_absorb()] and [design_absorb_wls()].
#' @examples
#' absorb_control()$solver
#' @export
absorb_control <- function(solver = getOption("xtx.fe_solver", "auto"),
                           precond = getOption("xtx.fe_precond", "auto"),
                           form_max = getOption("xtx.fe_form_max", 25),
                           refine = getOption("xtx.fe_refine", "auto")) {
    if (!checkmate::test_choice(solver, c("auto", "dense", "chol", "pcg")))
        stop("Option `xtx.fe_solver` must be \"auto\", \"dense\", \"chol\" or \"pcg\".", call. = FALSE)
    if (!checkmate::test_choice(precond, c("auto", "ichol", "jacobi")))
        stop("Option `xtx.fe_precond` must be \"auto\", \"ichol\" or \"jacobi\".", call. = FALSE)
    if (!checkmate::test_number(form_max, lower = 0))
        stop("Option `xtx.fe_form_max` must be a non-negative number.", call. = FALSE)
    auto <- identical(refine, "auto")
    if (!auto && !checkmate::test_count(refine))
        stop("Option `xtx.fe_refine` must be \"auto\" or a count.", call. = FALSE)
    list(solver = solver, precond = precond, pcg_tol = 1e-12, pcg_maxit = 10000L,
         refine = if (auto) 1L else as.integer(refine), refine_auto = auto, form_max = form_max,
         dense_max = 1000L)
}

## Component of each level of b in the bipartite graph of t(N_ab).
xtx_fe_components <- function(Nt) {
    i <- rep.int(seq_len(ncol(Nt)) - 1L, diff(Nt@p))
    .Call(C_xtx_fe_components, i, Nt@i, ncol(Nt), nrow(Nt))
}

## Conjugate gradients for R Z = B, preconditioned by IC(0) of the block
## diagonal of R (blocks: the FEs; the whole of R with two FE), or by its
## diagonal when IC(0) breaks down or with precond = "jacobi". R: a sparse
## symmetric matrix (formed), or a "xtx_fe_op" (R = A - N diag(d) N',
## applied without being formed, the rows of N mapped to those of R by
## `map`: Jacobi only). Z0: a starting value (a warm
## start; dropped for columns where it is worse than zero) or NULL. bmin: a
## floor for the norm of the right-hand sides in the stopping rule. The
## preconditioner is built once per `cache`. nthreads: threads for the
## products and sums of the recursion (identical results for a given
## number; one thread is the single-threaded code; small systems run on one
## thread; a negative count -T forces T threads, for the tests).
xtx_fe_pcg <- function(R, B, tol, maxit, block, precond = "jacobi", cache = new.env(), Z0 = NULL,
                      bmin = 0, nthreads = design_threads()) {
    if (is.null(cache$A)) {
        if (inherits(R, "xtx_fe_op")) {
            cache$A <- R$A
            cache$N <- R$N
            cache$map <- R$map
            cache$d <- R$d
            precond <- "jacobi"
        } else {
            cache$A <- xtx_fe_lower(R)
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
            Lx <- .Call(C_xtx_fe_ichol, P$p, P$i, P$x, ncol(A), 0)
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
    r <- .Call(C_xtx_fe_pcg, A@p, A@i, A@x, if (!is.null(N)) N@p, if (!is.null(N)) N@i,
               if (!is.null(N)) N@x, cache$map, cache$d, L$p, L$i, L$x, B, Z0,
               as.double(bmin), as.double(tol), as.integer(maxit), as.integer(nthreads))
    ## (nullv: on a breakdown, a direction close to a null vector of R)
    list(x = r[[1L]], iter = r[[2L]], resid = r[[3L]], converged = r[[3L]] <= tol,
         nullv = r[[4L]], precond = cache$precond)
}

## The lower triangle of a symmetric sparse matrix, as a dgCMatrix.
xtx_fe_lower <- function(R) {
    if (methods::is(R, "symmetricMatrix") && R@uplo == "U") R <- Matrix::t(R)
    R <- methods::as(methods::as(R, "generalMatrix"), "CsparseMatrix")
    if (!methods::is(R, "dsparseMatrix")) R <- methods::as(R, "dMatrix")
    Matrix::tril(R)
}

## Solve the reduced system R Z = B (R symmetric positive definite after the
## reference levels are removed, or semidefinite with a consistent
## right-hand side). R: a sparse matrix, or a "xtx_fe_op" (not formed: by
## conjugate gradients with the Jacobi preconditioner). Z0: a starting
## value for conjugate gradients (a warm start: the solution of a nearby
## system), or NULL. bmin: a floor for the norm of the right-hand sides in
## the stopping rule of conjugate gradients; by default the largest norm of
## the right-hand sides solved before with the same `cache` (a refinement
## step's right-hand side is at rounding level: its solution is needed to
## the accuracy of the first solve, not to tol relative to itself).
## Returns Z and solver details.
xtx_fe_solve <- function(R, B, opts, block, cache = new.env(), dense_max = opts$dense_max %||% 1000L,
                        Z0 = NULL, bmin = NULL) {
    m <- if (inherits(R, "xtx_fe_op")) R$m else nrow(R)
    if (m == 0L) return(list(z = B, solver = "none"))
    solver <- opts$solver
    if (solver == "auto") solver <- if (m <= dense_max) "dense" else "pcg"
    if (inherits(R, "xtx_fe_op")) solver <- "pcg"
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
    if (precond == "auto") precond <- if (inherits(R, "xtx_fe_op")) "jacobi" else "ichol"
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
        sol <- xtx_fe_pcg(xtx_fe_sub(R, use), B[use, , drop = FALSE], opts$pcg_tol, opts$pcg_maxit,
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
## "xtx_fe_op").
xtx_fe_sub <- function(R, use) {
    m <- if (inherits(R, "xtx_fe_op")) R$m else nrow(R)
    if (length(use) == m) return(R)
    if (!inherits(R, "xtx_fe_op")) return(R[use, use, drop = FALSE])
    R$A <- R$A[use, use, drop = FALSE]
    mp <- match(R$map + 1L, use) - 1L
    mp[is.na(mp)] <- -1L
    R$map <- mp
    R$m <- length(use)
    R
}

## The crosstabs that xtx_fe_system() takes, as pairs c(a, b) (crosstab(a,
## b) = t(N_ab)): the largest FE with each other, then the pairs of the
## others.
xtx_fe_pairs <- function(G) {
    o <- order(G, decreasing = TRUE)
    rest <- o[-1L]
    pairs <- lapply(rest, function(k) c(o[1L], k))
    if (length(rest) > 1L)
        for (r in seq_along(rest)[-length(rest)]) for (s in (r + 1L):length(rest))
            pairs <- c(pairs, list(c(rest[r], rest[s])))
    pairs
}

## The FE system K = F'WF, reduced by eliminating the largest FE.
## G: levels per FE; nk: weighted counts per level; crosstab(a, b): t(N_ab).
##
## R = A - N_1' D_1^-1 N_1 with A = blockdiag(n_k) + offdiag(N_kl) is formed
## (xtx_fe_schur(), lower triangle) unless forming it would cost more than
## applying it: forming gives up as soon as R has more nonzeros than N_1
## (min(form_max, 1) times), or the cost of its IC(0) factor (the sum of
## its squared column counts) exceeds form_max times nnz(N_1). R is then
## applied matrix-free (a "xtx_fe_op"), two passes over N_1 per product.
## R is always formed for the dense solver (m <= dense_max), for "chol"
## and for precond = "ichol".
xtx_fe_system <- function(G, nk, crosstab, opts) {
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
        comp <- xtx_fe_components(Nt[[r]])
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
    S <- .Call(C_xtx_fe_schur, A@p, A@i, A@x, N@p, N@i, N@x, map, 1 / sys$D1,
               as.double(if (is.finite(fm)) min(fm, 1) * nnzN else Inf), as.double(fm * nnzN))
    if (!is.null(S)) {
        sys$R <- methods::new("dsCMatrix", Dim = c(mk, mk), uplo = "L", p = S[[1L]], i = S[[2L]], x = S[[3L]])
    } else {
        sys$R <- structure(list(A = A, N = N, map = map, d = 1 / sys$D1, m = mk), class = "xtx_fe_op")
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
xtx_fe_apply <- function(sys, B, Z0 = NULL) {
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
    sol <- xtx_fe_solve(sys$R, Brr[sys$keep, , drop = FALSE], sys$opts, sys$block, sys$cache, Z0 = Z0)
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

## The kept columns and the factor of the FWL moments X'MX (fixest's rule:
## a column is aliased when its residual variance after the FE and the
## earlier kept columns is below tol times `scale`, its variance before
## the FE are absorbed), by the pivoted Cholesky kernel of qr()
## (xtx_pivchol(), ties drop the column). `band` and `recheck`: suspect
## columns, see xtx_pivchol(). Returns the kept columns and the upper factor
## of M[keep, keep].
xtx_fe_pivchol <- function(M, scale, tol, band = NULL, recheck = NULL) {
    ch <- xtx_pivchol(M, scale, tol, "drop", band = band, recheck = recheck)
    k <- ch$rank
    list(keep = ch$pivot[seq_len(k)], R = ch$R[seq_len(k), seq_len(k), drop = FALSE])
}

## ---------------------------------------------------------------------------
## The absorption operator
##
## The FE of a fit, absorbed from sums over rows. xtx_absorb() builds it from
## the one fused pass of xtx_p_fe_sums() (F_k'W[X z] and F_k'W1 for every FE
## and the crosstabs N_kl: the pass must see x and z, so the operator is
## built for them), and reduces and sets up the FE system K = F'WF for
## solves once (xtx_fe_system()). It is a list of functions:
##
##   solve(B, warm)    K^- B (split by FE) and Q = B'K^-B, for B a list of
##                     G_k x q matrices (xtx_fe_apply(); warm: a starting
##                     value of the reduced solve, the `Z` of an earlier
##                     call, or NULL)
##   project(M, warm)  the FWL moments M - Q of [X z] (M = [X z]'W[X z] on
##                     the columns xcols and z, B = F'W[X z]): list(M, info,
##                     Z); keeps Gamma = K^-B for the functions below
##   effects(b, j)     the FE values of column j of [X z] (default: z) less
##                     those of X b: Gamma_j - Gamma_X b, one vector per FE
##   update(fe, Ae, kp, delta)  FE values after a refinement step: fe +
##                     Ae - Gamma_X[, kp] delta
##   cross(g, Ae)      g - sum_k S_k' Ae_k (S_k = F_k'WX): X'M v from X'W v
##                     and Ae = K^- F'W v
##   variance(dxx)     the variances of the columns of X before absorption
##                     (sums of squares around the weighted mean), from
##                     their sums of squares dxx
##   scores(we, cl)    the cluster scores sum_{i in c} we_i (MX)_i, corrected
##                     for FE not nested in the clusters (cl: a xtx_c_codes()
##                     grouping; see the header): list(U, nested)
##
## and G, the levels per FE.
xtx_absorb <- function(x, z, w, codes, xcols, opts) {
    G <- vapply(codes, function(g) g$G, 1L)
    px <- length(xcols)
    xi <- seq_len(px)
    pairs <- xtx_fe_pairs(G)
    ## (with the crosstabs of the FE system: one or two passes in memory)
    fs <- xtx_p_fe_sums(x, w, codes, extra = z, pairs = pairs)
    nk <- lapply(fs$gs, `[[`, "n")
    if (any(vapply(nk, function(v) any(!(v > 0)), NA)))
        stop("internal error: a fixed-effect level without weight", call. = FALSE)
    B <- lapply(fs$gs, function(g) cbind(g$S[, xcols, drop = FALSE], g$extra[, 1L]))
    ## t(N_ab), G_b x G_a
    pkey <- vapply(pairs, paste, "", collapse = ",")
    sys <- xtx_fe_system(G, nk, function(a, b) fs$ct[[match(paste(a, b, sep = ","), pkey)]], opts)
    fs <- NULL
    Gam <- NULL
    solve <- function(B, warm = NULL) xtx_fe_apply(sys, B, Z0 = warm)
    project <- function(M, warm = NULL) {
        ab <- solve(B, warm)
        Gam <<- ab$A
        list(M = M - ab$Q, info = ab$info, Z = ab$Z)
    }
    effects <- function(b, j = px + 1L) lapply(Gam, function(Ak) drop(Ak[, j] - Ak[, xi, drop = FALSE] %*% b))
    update <- function(fe, Ae, kp, delta) {
        for (k in seq_along(fe)) fe[[k]] <- fe[[k]] + Ae[[k]][, 1L] - drop(Gam[[k]][, kp, drop = FALSE] %*% delta)
        fe
    }
    cross <- function(g, Ae) {
        for (k in seq_along(B)) g <- g - drop(crossprod(B[[k]][, xi, drop = FALSE], Ae[[k]]))
        g
    }
    variance <- function(dxx) dxx - colSums(B[[1L]])[xi]^2 / sum(nk[[1L]])
    scores <- function(we, cl) {
        ## (a grouping of design_groups() as codes)
        cl <- xtx_absorb_codes(cl)[[1L]]
        U <- xtx_p_colsums_groups(xtx_p_scale(x, we), list(cl))[[1L]]$S[, xcols, drop = FALSE]
        ## an FE is nested in the clusters when each of its levels is in
        ## one cluster; otherwise its scores need a correction
        nested <- vapply(seq_along(codes), function(k) xtx_p_nested(cl, codes[[k]]), NA)
        for (k in which(!nested)) {
            E <- xtx_p_crosstab(codes[[k]], cl, we) # G_k x Gc
            U <- U - as.matrix(Matrix::crossprod(E, Gam[[k]][, xi, drop = FALSE]))
        }
        list(U = U, nested = nested)
    }
    structure(list(G = G, solve = solve, project = project, effects = effects, update = update, cross = cross,
                   variance = variance, scores = scores), class = "xtx_absorb")
}

xtx_fe_core <- function(x, cols, tol, opts, warm = NULL, resid_final = TRUE) {
    xcols <- setdiff(seq_len(ncol(x)), which(colnames(x) == "(Intercept)"))
    px <- length(xcols)
    xi <- seq_len(px)
    n <- xtx_nrows(x)
    mo <- xtx_fe_moments(x, cols$z, cols$w, xcols, n)
    A <- xtx_absorb(x, cols$z, cols$w, cols$codes, xcols, opts)
    pj <- A$project(mo$M, warm = if (!is.null(warm)) warm$Z)
    if (!is.null(warm)) {
        warm$Z <- pj$Z
        warm$iter <- c(warm$iter, pj$info$iter %||% 0L)
    }
    M <- pj$M
    dxx <- diag(mo$M)[xi]
    scale <- A$variance(dxx)
    su <- xtx_fe_suspects(x, cols, A, xcols, M, scale, dxx, tol)
    pc <- xtx_fe_pivchol(M[xi, xi, drop = FALSE], scale, tol, band = su$band, recheck = su$recheck)
    kp <- pc$keep
    b0 <- numeric(px)
    if (length(kp)) b0[kp] <- backsolve(pc$R, backsolve(pc$R, M[kp, px + 1L], transpose = TRUE))
    rf <- xtx_fe_refine(x, cols, A, xcols, pc, b0, A$effects(b0), opts, resid_final)
    b <- rep(NA_real_, px)
    b[kp] <- rf$b[kp]
    bread <- matrix(NA_real_, px, px)
    if (length(kp)) bread[kp, kp] <- chol2inv(pc$R)
    dof <- xtx_fe_dof(A$G, pj$info)
    out <- list(b = b, kp = kp, xcols = xcols, bread = bread, rss = rf$r$rss, fe = rf$fe,
                rank = length(kp), G = A$G, fe_dof = dof$fe_dof, fe_dof_exact = dof$fe_dof_exact,
                df.residual = n - length(kp) - dof$fe_dof, nobs = n, tss = mo$tss,
                tss_within = max(M[px + 1L, px + 1L], 0), info = c(pj$info, list(resid_pass = rf$info)),
                we = rf$r$we, absorb = A)
    class(out) <- "xtx_absorb_fit"
    out
}

## The suspect band of the pivoted Cholesky of a fit with FE (xtx_pivchol()):
## the residual variance d of a column after the FE and the earlier kept
## columns comes from the aggregates, X'WX - S'K^-S, accurate only to about
## c eps diag(X'WX). A column with |d - tol * scale| <= c eps diag(X'WX)
## cannot be judged from them: its d is computed again from the data, the
## sum of w e^2 of its residual e = x_j - X_K c - F f (c from the factor of
## the columns K kept before it, f = Gamma_j - Gamma_K c), one pass over
## the rows per suspect (xtx_c_fe_resid()). That residual is not refined, but
## the error of c and f from the aggregates changes its sum of squares only
## to second order, far below the band.
##
## c = 1e4 (c eps ~ 2.2e-12): the aggregates are one-pass sums (rounding of
## order eps each), their difference cancels (absolute error of order eps
## diag(X'WX)), and with conjugate gradients K^-S is accurate only to the
## stopping rule, a normwise backward error of 1e-12 (~4500 eps), which
## enters S'K^-S at that order relative to diag(X'WX); c covers that with a
## margin of 2. With centered columns and tol = 1e-10, a column is a
## suspect only when its residual variance is within about 2% of the
## threshold: rare, so the extra passes cost nothing in general. Without
## a positive scale a column is never taken, nor a suspect.
## Returns list(band (lo, hi), recheck) for xtx_fe_pivchol().
xtx_fe_suspects <- function(x, cols, A, xcols, M, scale, dxx, tol) {
    width <- 1e4 * .Machine$double.eps * abs(dxx)
    ok <- is.finite(scale) & scale > 0 & is.finite(width)
    thr <- tol * scale
    band <- list(lo = ifelse(ok, thr - width, Inf), hi = ifelse(ok, thr + width, -Inf))
    px <- length(xcols)
    recheck <- function(ch, sus) {
        ## (a zero column, made only when there are suspects)
        z0 <- xtx_c_const(cols$z, 0)
        acc <- ch$pivot[seq_len(ch$rank)]
        d <- vapply(sus, function(j) {
            K <- acc[acc < j]
            a <- length(K)
            Ra <- ch$R[seq_len(a), seq_len(a), drop = FALSE]
            bk <- numeric(px)
            if (a) bk[K] <- backsolve(Ra, backsolve(Ra, M[K, j], transpose = TRUE))
            beta <- numeric(ncol(x))
            beta[xcols] <- bk
            beta[xcols[j]] <- -1
            ## e = 0 - (X_K c - x_j) - F f
            r <- xtx_c_fe_resid(z0, xtx_p_mult(x, beta), cols$w, A$effects(bk, j), cols$codes)
            r[[3L]]
        }, 0)
        list(dres = d, taken = is.finite(d) & d > 0 & d > thr[sus])
    }
    list(band = band, recheck = recheck)
}

## The moments of a fit with FE: M = [X z]'W[X z] on the columns xcols of
## x and the response z (weights w, NULL for none; n rows), and the total
## sum of squares of z around its weighted mean.
xtx_fe_moments <- function(x, z, w, xcols, n) {
    mo <- lm_moments(x, z, w, n)
    xtx <- unname(mo$xtx)
    xty <- unname(mo$xty)[, 1L]
    if (!all(is.finite(xtx)) || !all(is.finite(xty))) stop("NA/NaN/Inf in 'x'", call. = FALSE)
    XtX <- xtx[xcols, xcols, drop = FALSE]
    Xty <- xty[xcols]
    list(M = rbind(cbind(XtX, Xty), c(Xty, mo$yty)), tss = mo$tss)
}

## The residuals e = z - X b - sum_k fe_k[code_k] of a fit with FE (b: the
## coefficients of the columns xcols of x) and their sums: X'We, F_k'We,
## rss = sum(w e^2) and we = w e (one pass in memory).
xtx_fe_residuals <- function(x, cols, xcols, b, fe) {
    beta <- numeric(ncol(x))
    beta[xcols] <- b
    ## (remotely the product and the lookups of the FE bind on the frame of
    ## z and w, at the sums of xtx_c_fe_resid())
    r <- xtx_c_fe_resid(cols$z, xtx_p_mult(x, beta), cols$w, fe, cols$codes)
    xte <- xtx_p_crossprod(x, r[[1L]], compress = FALSE)
    list(xte = unname(xte)[, 1L], Be = lapply(r[[2L]], function(m) m[, 1L]),
         rss = r[[3L]], we = drop(r[[1L]]))
}

## Iterative refinement of a fit with FE (xtx_refine()): the normal
## equations X'M e = 0 from residuals computed row by row correct the
## rounding of the one-pass aggregates. b, fe: the solution from the
## aggregates (b zero for aliased columns), pc its pivoted factor, A the
## absorption operator. Stopping rule: opts$refine steps; with fe_refine =
## "auto", the pass over the rows after the last step is skipped when the
## step's correction is negligible (xtx_fe_negligible()), and the residuals
## of the last pass stand for those after the step. Returns list(b, fe, r
## (the residuals; with resid_final = FALSE, after at least one step, rss
## NA and no `we`), info).
xtx_fe_refine <- function(x, cols, A, xcols, pc, b, fe, opts, resid_final) {
    kp <- pc$keep
    resid <- function(s) xtx_fe_residuals(x, cols, xcols, s$b, s$fe)
    correct <- function(s, r) {
        ## (its right-hand side is at rounding level: the solver floors its
        ## norm at that of the first solve with the system, see
        ## xtx_fe_solve())
        Ae <- A$solve(lapply(r$Be, function(v) matrix(v, ncol = 1L)))$A
        b <- s$b
        delta <- numeric(0)
        if (length(kp)) {
            g <- A$cross(r$xte[xcols], Ae)
            delta <- backsolve(pc$R, backsolve(pc$R, g[kp], transpose = TRUE))
            b[kp] <- b[kp] + delta
        }
        list(sol = list(b = b, fe = A$update(s$fe, Ae, kp, delta)), Ae = Ae, delta = delta)
    }
    stop <- function(step, r, it) {
        it == opts$refine && isTRUE(opts$refine_auto) && xtx_fe_negligible(r, step$Ae, step$delta, pc$R)
    }
    rf <- xtx_refine(list(b = b, fe = fe), resid, correct, opts$refine, stop)
    r <- if (rf$stopped) rf$r
         else if (resid_final || opts$refine == 0L) resid(rf$sol)
         else list(rss = NA_real_, we = NULL)
    list(b = rf$sol$b, fe = rf$sol$fe, r = r, info = if (rf$stopped) "skipped" else "done")
}

## Degrees of freedom of the FE (fixest: K.fixef = "nonnested", K.exact =
## FALSE), and the exact rank of K (fixest's K.exact = TRUE): with two FE,
## one level per connected component is redundant; with more, the count is
## exact when the reduced system was factored (a pivoted dense Cholesky,
## which gives its rank, or a sparse Cholesky that succeeded); otherwise NA.
## info: the solver details of the absorption.
xtx_fe_dof <- function(G, info) {
    nfe <- length(G)
    exact <- if (nfe == 1L) {
        G
    } else if (nfe == 2L || info$solver %in% c("dense", "chol")) {
        max(G) + info$m
    } else {
        NA_integer_
    }
    list(fe_dof = sum(G) - (nfe - 1L), fe_dof_exact = exact)
}

## Is the correction of a refinement step negligible for the residuals,
## so that the pass over the rows that recomputes them can be skipped?
## The step moves the residuals e by d = F Ae + X delta (Ae = K^- F'We,
## delta the change of b, `Rc` the Cholesky factor of X'MX on the kept
## columns), which lowers sum(w e^2) by Be'Ae + delta' X'MX delta to first
## order (Be = F'We: e is orthogonal to the columns of [X F] at the
## solution) and changes e by about the square root of that, relative to
## the residuals. Negligible: that decrease is below 1e-20 times the
## residual sum of squares (residuals, hence sigma and the clustered
## scores, accurate to 1e-10, the rss to 1e-20); otherwise (ill
## conditioned or near collinear systems, uncentered variables, which the
## step corrects) the residuals are recomputed. The coefficients and FE
## values are those of the step either way.
xtx_fe_negligible <- function(r, Ae, delta, Rc) {
    if (!is.finite(r$rss) || !(r$rss > 0)) return(FALSE)
    q <- sum(vapply(seq_along(Ae), function(k) sum(r$Be[[k]] * Ae[[k]][, 1L]), 0))
    if (length(delta)) q <- q + sum(drop(Rc %*% delta)^2)
    is.finite(q) && abs(q) <= 1e-20 * r$rss
}

## X beta plus the FE values of the rows' levels, as a column (codes: per
## FE, integer codes, or, with `levels`, key columns matched to levels).
## With beta a p x k matrix and each fe[[k]] a G_k x k matrix, k columns.
xtx_fe_fitted <- function(x, beta, fe, codes, levels = NULL) {
    xtx_c_fe_add(xtx_p_mult(x, beta), fe, codes, levels)
}

## Map a fit on the centered design X~ = X T back to X: b = T b~ (but the
## intercept column, absorbed by the FE, which take the constant t0'b~),
## bread T bread~ T'. NULL when the aliased columns do not allow it (a kept
## column of X~ uses an aliased column of X).
xtx_fe_uncenter <- function(est, Tm, icol) {
    xc <- est$xcols
    Tx <- Tm[xc, xc, drop = FALSE]
    kp <- est$kp
    al <- setdiff(seq_along(xc), kp)
    if (!xtx_uncenter_ok(Tx, kp, al)) return(NULL)
    bt <- est$b
    bt[is.na(bt)] <- 0
    b <- drop(Tx %*% bt)
    b[al] <- NA_real_
    est$b <- b
    if (length(kp)) {
        Tk <- Tx[kp, kp, drop = FALSE]
        est$bread[kp, kp] <- xtx_uncenter_cov(Tk, est$bread[kp, kp, drop = FALSE])
        if (!is.null(est$vcov_cluster))
            est$vcov_cluster[kp, kp] <- xtx_uncenter_cov(Tk, est$vcov_cluster[kp, kp, drop = FALSE])
    }
    ## X~ b~ = X b + (t0'b~) 1: the constant goes to the first FE
    if (length(icol)) est$fe[[1L]] <- est$fe[[1L]] + sum(Tm[icol, xc] * bt)
    est
}

#' The rows of a fit with fixed effects
#'
#' `design_absorb_rows()` selects the rows a fit with fixed effects (FE) uses
#' and codes the FE on them, with one code path for data in memory and remote
#' data: rows with a missing or non-finite response or weight, a zero weight,
#' or a missing FE or cluster value are dropped; then, iteratively until none
#' is left, rows in a level of some FE that the FE fit perfectly: levels with
#' one row (singletons, `singletons = TRUE`), and levels whose response is all
#' 0 (`zero = TRUE`) or all 1 (`one = TRUE`), whose FE would be infinite in a
#' GLM (fixest's `fixef.rm`).
#'
#' Names may change before the first release.
#'
#' @param problem The input stage of the fit, a [design_problem()] made with
#'   `fe` (and `cluster`): its design `x`, response `y` (the response the FE
#'   are fitted to, such as the response minus an offset), weights `w`, rows
#'   `rows`, FE and cluster keys `keys`, and `more`, a list of other columns
#'   of the rows to keep (offset, GLM columns), all of which the fit may set
#'   before the call.
#' @param singletons,zero,one Which levels to drop, see Details.
#' @return A list: `x` (the design of the rows kept; `x_rows` the same),
#'   `cols` (`z`, the response; `w`, the weights or `NULL`; `codes`, the FE
#'   groupings, one per FE; `cl`, the cluster grouping or `NULL`), `more` (the
#'   columns of `problem$more` on those rows), `levels` (per FE, a data frame
#'   of the levels in code order), `nsingle` and `nperfect` (the rows dropped
#'   as singletons and as perfectly fit), `rows` (in memory, the rows of
#'   the model frame kept) and `frame` (in memory, the model frame of those
#'   rows; `NULL` remotely). The groupings are for [design_absorb_wls()],
#'   [design_absorb()] and [design_absorb_fitted()]; `g$G` is the number of
#'   levels of grouping `g`.
#' @export
design_absorb_rows <- function(problem, singletons = TRUE, zero = FALSE, one = FALSE) {
    checkmate::assert_flag(singletons)
    checkmate::assert_flag(zero)
    checkmate::assert_flag(one)
    base <- problem
    clv <- problem$cluster
    fe_names <- names(problem$keys$fe)
    x <- base$x
    z <- base$y
    w <- base$w
    rows <- base$rows
    more <- base$more
    keys <- base$keys$fe
    clkey <- base$keys$cl
    restrict <- function(keep) {
        ## (keep is never missing: the logical columns are NA-safe)
        nd <- xtx_colsum(!keep)
        if (nd == 0) return(0)
        x2 <- xtx_p_subset(x, keep)
        z <<- xtx_c_subset(z, keep, x2)
        w <<- xtx_c_subset(w, keep, x2)
        rows <<- xtx_c_subset(rows, keep, x2)
        more <<- lapply(more, xtx_c_subset, keep, x2)
        keys <<- lapply(keys, xtx_c_subset, keep, x2)
        clkey <<- xtx_c_subset(clkey, keep, x2)
        x <<- x2
        if (xtx_nrows(x) == 0) stop("0 (non-NA) cases", call. = FALSE)
        nd
    }
    if (xtx_nrows(x) == 0) stop("0 (non-NA) cases", call. = FALSE)
    keep <- xtx_c_is_finite(z)
    if (!is.null(w)) keep <- keep & xtx_c_is_finite(w) & w != 0
    for (k in c(keys, list(clkey))) if (!is.null(k)) keep <- keep & xtx_c_complete(k)
    restrict(keep)
    rm_single <- singletons
    rm0 <- zero
    rm1 <- one
    nsingle <- nperfect <- 0
    repeat {
        cds <- unname(lapply(keys, xtx_c_dense_codes, raw = TRUE))
        if (!rm_single && !rm0 && !rm1) break
        ## per level: rows, rows with a response not 0, not 1 (without
        ## bounds, the rows per level: no columns to form; in memory they
        ## come with the codes, xtx_fe_code())
        if (!rm0 && !rm1) {
            cnt <- lapply(cds, function(cd) {
                m <- xtx_p_groupcount(cd$codes, cd$n)
                matrix(m, length(m), 3L)
            })
        } else {
            one <- xtx_c_const(z, 1)
            cols <- xtx_c_cbind(one, if (rm0) xtx_c_count(z != 0) else one, if (rm1) xtx_c_count(z != 1) else one)
            cnt <- lapply(cds, function(cd) xtx_p_groupsum(cols, cd$codes))
        }
        ## per level, once: perfectly fit (counted as such when it is also
        ## a singleton), or a singleton
        masks <- lapply(cnt, function(m) {
            pf_k <- (rm0 & m[, 2L] == 0) | (rm1 & m[, 3L] == 0)
            cbind(pf_k, rm_single & m[, 1L] == 1 & !pf_k)
        })
        ## no level to drop (no lookup of the rows' levels)
        if (!any(vapply(masks, any, NA))) break
        single <- perfect <- NULL
        for (k in seq_along(keys)) {
            ## (the flags of the rows' levels)
            lk <- xtx_c_lookup_flags(masks[[k]], cds[[k]]$codes)
            p1 <- lk[[1L]]
            s1 <- lk[[2L]]
            perfect <- if (is.null(perfect)) p1 else perfect | p1
            single <- if (is.null(single)) s1 else single | s1
        }
        drop <- perfect | single
        s <- xtx_c_sums(xtx_c_count(perfect), xtx_c_count(drop))
        nd <- restrict(!drop)
        if (nd == 0) break
        nperfect <- nperfect + s[1L]
        nsingle <- nsingle + (s[2L] - s[1L])
    }
    cl <- NULL
    if (!is.null(clv)) {
        m <- match(clv, fe_names)
        if (is.na(m)) {
            cl <- xtx_c_dense_codes(clkey, raw = TRUE)$codes
        } else {
            cl <- cds[[m]]$codes
        }
        if (cl$G < 2L) stop("Clustered standard errors need at least two clusters.", call. = FALSE)
    }
    ## Keep the filtered source for grouped summaries and on-demand row outputs.
    x <- xtx_p_rebase(x, list(z, w, more, lapply(cds, `[[`, "codes"), cl))
    on <- function(col) xtx_c_on(col, x)
    oncodes <- function(g) xtx_cols_map(g, on)
    list(x = x, x_rows = x,
         cols = list(z = on(z), w = on(w), codes = lapply(cds, function(cd) oncodes(cd$codes)), cl = oncodes(cl)),
         more = lapply(more, on), levels = lapply(cds, `[[`, "levels"), nsingle = as.integer(nsingle),
         nperfect = as.integer(nperfect), rows = rows, frame = if (!inherits(x, "xtx_remote")) x$vars)
}


## ---------------------------------------------------------------------------
## Exports: the operator, the estimator, fitted values, centering

## FE groupings as xtx_c_codes(): those of design_absorb_rows(), or
## groupings of design_groups() (whose codes, with `levels`, can be NA).
xtx_absorb_codes <- function(groups) {
    if (xtx_is_codes(groups) || inherits(groups, "xtx_groups")) groups <- list(groups)
    if (!is.list(groups) || !length(groups))
        stop("`groups` must be a list of groupings (design_absorb_rows(), design_groups()).", call. = FALSE)
    lapply(groups, function(g) {
        if (xtx_is_codes(g)) return(g)
        if (inherits(g, "xtx_groups")) {
            xtx_groups_check(g, "groups", na_ok = TRUE)
            return(xtx_groups_codes(g))
        }
        stop("`groups` must be a list of groupings (design_absorb_rows(), design_groups()).", call. = FALSE)
    })
}

#' Absorbed fixed effects
#'
#' The estimator of linear models with fixed effects (FE) of [design_lm()],
#' which [design_glm()] also uses for the weighted least squares steps of its
#' iterations, from sums over rows (Frisch-Waugh-Lovell on aggregates: the
#' dummies of the FE are never formed and the rows never demeaned), with one
#' code path for designs in memory and remote designs.
#'
#' * `design_absorb()` builds the absorption operator of the FE for a
#'   response `z` and weights `w` (one fused pass over the rows for every
#'   grouped sum and crosstab), with the FE system reduced and set up for
#'   solves: an object of class `xtx_absorb`, a list of functions (see
#'   *The absorption operator*).
#' * `design_absorb_wls()` fits `z` on the columns of `x` (but the
#'   intercept, absorbed) and the FE by weighted least squares: the moments
#'   of `[X z]` with the FE absorbed, columns judged collinear by fixest's
#'   rule (residual variance below `tol` times the variance before
#'   absorption), the solution refined from residuals computed row by row,
#'   and the FE values.
#' * `design_absorb_fitted()` is `X beta` plus the FE values of the rows'
#'   levels, as a column.
#'
#' Names may change before the first release.
#'
#' @section The absorption operator:
#' A list of class `xtx_absorb` with `G` (the levels per FE) and the
#' functions
#' \describe{
#'   \item{`solve(B, warm)`}{`K^- B` (split by FE) and `Q = B'K^-B`, `K = F'WF`,
#'     for `B` a list of `G_k x q` matrices (`warm`: a starting value of the
#'     reduced solve, the `Z` of an earlier call, or `NULL`).}
#'   \item{`project(M, warm)`}{The moments `M - Q` of `[X z]` with the FE
#'     absorbed, `M = [X z]'W[X z]` on the columns of `x` but the intercept
#'     and `z`: `list(M, info, Z)`; keeps `Gamma = K^- F'W[X z]` for the
#'     functions below.}
#'   \item{`effects(b, j)`}{The FE values of column `j` of `[X z]` (default:
#'     `z`) less those of `X b`, one vector per FE.}
#'   \item{`update(fe, Ae, kp, delta)`}{FE values after a refinement step.}
#'   \item{`cross(g, Ae)`}{`X'M v` from `X'W v` and `Ae = K^- F'W v`.}
#'   \item{`variance(dxx)`}{The variances of the columns of `X` before
#'     absorption, from their sums of squares.}
#'   \item{`scores(we, cl)`}{The cluster scores `sum_{i in c} we_i (MX)_i`
#'     (the raw scores, with no small sample correction), corrected for FE
#'     not nested in the clusters (`cl` a grouping): `list(U, nested)`.}
#' }
#'
#' @param x A design (`design_matrix` or remote design), unscaled. Its
#'   intercept column, if any, is absorbed by the FE.
#' @param z The response: a numeric vector with one value per row of `x` in
#'   memory, a [design_column()] of `x` remotely (as [design_absorb_rows()]
#'   returns it).
#' @param w Weights (positive), in the same form as `z`, or `NULL`.
#' @param groups The FE: a list of groupings, from [design_absorb_rows()]
#'   (`cols$codes`) or [design_groups()]. For `design_absorb_fitted()`,
#'   groupings of new rows made with `design_groups(x, key, levels =)` give
#'   `NA` for levels the fit does not have.
#' @param tol The tolerance for collinear columns (fixest's `collin.tol`).
#' @param control Options of the solver, [absorb_control()].
#' @param warm `NULL`, or an environment that keeps the solution of the
#'   reduced system from one call to the next (`warm$Z`, a warm start for
#'   sequences of nearby problems, such as the steps of a GLM) and the
#'   numbers of iterations of conjugate gradients (`warm$iter`).
#' @param resid Compute the residuals after the last refinement step (the
#'   residual sum of squares `rss` and the weighted residuals `we`, for
#'   clustered standard errors). `FALSE` skips that pass when the refinement
#'   step did not compute them.
#' @param beta The coefficients of the columns of `x` (a vector, or a matrix
#'   with one column per vector).
#' @param fe The FE values: a list with one vector per FE (in code order),
#'   or one matrix per FE with as many columns as `beta`.
#' @return `design_absorb()`: a `xtx_absorb` (see *The absorption
#'   operator*). `design_absorb_wls()`: a list of class `xtx_absorb_fit`:
#'   `b` (the coefficients of the columns of `x` but the intercept, `NA` for
#'   collinear ones), `kp` (the columns kept), `xcols` (the columns of `x`
#'   fitted), `bread` (the inverse of `X'MX` on the kept columns, `NA`
#'   elsewhere), `rss`, `fe` (the FE values, one vector per FE), `rank`, `G`,
#'   `fe_dof` (the FE parameters, `sum(G) - (number of FE - 1)`),
#'   `fe_dof_exact` (the exact rank of the FE, or `NA`), `df.residual`,
#'   `nobs`, `tss` (total sum of squares of `z` around its weighted mean),
#'   `tss_within` (after the FE), `info` (solver details), `we` (the weighted
#'   residuals) and `absorb` (the operator, for `absorb$scores()`).
#'   `design_absorb_fitted()`: a column (a numeric vector or matrix in
#'   memory, a [design_column()] remotely).
#' @examples
#' set.seed(1)
#' d <- data.frame(y = rnorm(200), x = rnorm(200), id = sample(1:20, 200, TRUE))
#' D <- design_matrix(y ~ x, d)
#' g <- design_groups(D, d["id"])
#' fit <- design_absorb_wls(D, d$y, NULL, list(g))
#' all.equal(fit$b, unname(coef(lm(y ~ x + factor(id), d))["x"]))
#' @export
design_absorb <- function(x, z, w = NULL, groups, control = absorb_control()) {
    checkmate::assert_multi_class(x, c("xtx_matrix", "xtx_remote"))
    codes <- xtx_absorb_codes(groups)
    xcols <- setdiff(seq_len(ncol(x)), which(colnames(x) == "(Intercept)"))
    xtx_absorb(x, z, w, codes, xcols, control)
}

#' @rdname design_absorb
#' @export
design_absorb_wls <- function(x, z, w = NULL, groups, tol = 1e-10, control = absorb_control(), warm = NULL,
                              resid = TRUE) {
    checkmate::assert_multi_class(x, c("xtx_matrix", "xtx_remote"))
    if (!is.null(warm) && !is.environment(warm)) stop("`warm` must be NULL or an environment.", call. = FALSE)
    xtx_fe_core(x, list(z = z, w = w, codes = xtx_absorb_codes(groups)), tol, control, warm = warm,
                resid_final = isTRUE(resid))
}

#' @rdname design_absorb
#' @export
design_absorb_fitted <- function(x, beta, fe, groups) {
    xtx_fe_fitted(x, beta, fe, xtx_absorb_codes(groups))
}

#' @rdname design_absorb
#' @param ... Not used.
#' @export
print.xtx_absorb <- function(x, ...) {
    cat(sprintf("<design_absorb: %d fixed effect%s (%s levels)>\n", length(x$G), if (length(x$G) == 1L) "" else "s",
                paste(x$G, collapse = ", ")))
    invisible(x)
}

#' Centered designs
#'
#' `design_center()` shifts the numeric variables of a design by their mean
#' (in a sample of rows) when that is an exact linear transform of the
#' design, `X~ = X T` with `T` unit upper triangular: a fit on the centered
#' design keeps the accuracy of sums of squares and cross products when a
#' variable is far from zero. In memory the design is rebuilt from its model
#' frame; remotely the variables' SQL changes. `design_uncenter()` maps a fit
#' of [design_absorb_wls()] on the centered design back to the design: `b = T
#' b~` (the intercept column, absorbed, moves to the first fixed effect),
#' `bread` (and `vcov_cluster`, when the fit has one) `T V T'`.
#'
#' Names may change before the first release.
#'
#' @param x A design: a `design_matrix` made from a formula and a data frame
#'   (one that keeps its model frame), or a remote design.
#' @param cols Columns of `x` to carry over to the centered design (a list,
#'   possibly nested, of columns and groupings, as [design_absorb_rows()]
#'   returns them): remotely re-made over the new design, in memory unchanged.
#' @param fit A fit of [design_absorb_wls()] on `center$x`.
#' @param center The result of `design_center()`.
#' @return `design_center()`: `NULL` when no exact centering applies, else a
#'   list with `x` (the centered design), `T` (the transform) and `cols`.
#'   `design_uncenter()`: the fit mapped back, or `NULL` when the collinear
#'   columns do not allow it (a kept column of the centered design uses a
#'   collinear column of the design): then fit the design itself.
#' @examples
#' d <- data.frame(y = rnorm(100), x = rnorm(100) + 1e6, id = rep(1:10, 10))
#' D <- design_matrix(y ~ x, d)
#' ctr <- design_center(D)
#' colMeans(as.matrix(ctr$x))
#' @export
design_center <- function(x, cols = list()) {
    checkmate::assert_multi_class(x, c("xtx_matrix", "xtx_remote"))
    ctr <- xtx_p_center(x)
    if (is.null(ctr)) return(NULL)
    xc <- ctr$x
    list(x = xc, T = ctr$T, cols = xtx_cols_map(cols, function(col) xtx_c_on(col, xc)))
}

#' @rdname design_center
#' @export
design_uncenter <- function(fit, center) {
    xtx_fe_uncenter(fit, center$T, which(colnames(center$x) == "(Intercept)"))
}

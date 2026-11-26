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
## (GROUP BY queries and the products of remote designs, fe_remote.R); the
## two differ only in how they compute these sums (a "backend").
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

## Keys of FE levels, for matching new data: one string per row of `lv`
## (a data frame of level values).
dm_fe_keys <- function(lv) {
    v <- lapply(lv, function(x) if (is.double(x)) format(x, digits = 17, trim = TRUE) else as.character(x))
    if (length(v) == 1L) return(v[[1L]])
    do.call(paste, c(v, sep = "\r"))
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
dm_fe_pcg <- function(R, B, tol, maxit, block, precond = "ichol") {
    A <- methods::as(methods::as(R, "generalMatrix"), "CsparseMatrix")
    Lp <- Li <- Lx <- NULL
    if (precond == "ichol") {
        P <- Matrix::tril(A)
        if (length(unique(block)) > 1L) {
            Ts <- Matrix::summary(P)
            same <- block[Ts$i] == block[Ts$j]
            P <- Matrix::sparseMatrix(i = Ts$i[same], j = Ts$j[same], x = Ts$x[same],
                                      dims = dim(P), repr = "C")
        }
        P <- methods::as(methods::as(P, "generalMatrix"), "CsparseMatrix")
        Lx <- .Call(C_dm_fe_ichol, P@p, P@i, P@x, nrow(P))
        if (!is.null(Lx)) {
            Lp <- P@p
            Li <- P@i
        } else {
            precond <- "jacobi"
        }
    }
    storage.mode(B) <- "double"
    r <- .Call(C_dm_fe_pcg, A@p, A@i, A@x, Lp, Li, Lx, B, as.double(tol), as.integer(maxit))
    list(x = r[[1L]], iter = r[[2L]], resid = r[[3L]], converged = r[[3L]] <= tol,
         precond = precond)
}

## Solve the reduced system R Z = B (R symmetric positive definite after the
## reference levels are removed, or semidefinite with a consistent
## right-hand side). Returns Z and solver details.
dm_fe_solve <- function(R, B, opts, block, dense_max = 3000L) {
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
    sol <- dm_fe_pcg(R, B, opts$pcg_tol, opts$pcg_maxit, block, opts$precond)
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
    sol <- dm_fe_solve(sys$R, Brr[sys$keep, , drop = FALSE], sys$opts, sys$block)
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
## The estimator, on a backend
##
## A backend is a list:
##   n (rows), p (design columns), icol (the intercept column, or empty),
##   G (levels per FE);
##   moments: list(xtx, xty, yty, tss) for the whole design (tss: weighted
##       sum of squares of y around its weighted mean);
##   gs: per FE, list(nk = weighted counts, S = G_k x p grouped sums of w X,
##       s = grouped sums of w y);
##   crosstab(a, b): t(N_ab), weighted, a G_b x G_a dgCMatrix;
##   resid(beta, fe): residuals e = y - X beta - sum_k fe_k[level], kept
##       for cluster(); returns list(xte = X'We, Be = per FE, grouped sums
##       of w e, rss = sum w e^2);
##   cluster(): after resid(), list(G = number of clusters,
##       U = colSums(w e X, group = cl), nested = per FE, whether it is
##       nested in the clusters, E(k) = t(crosstab(cl, FE_k, w e)),
##       a G_k x G dgCMatrix).

dm_fe_core <- function(be, tol, opts, cluster = FALSE) {
    p <- be$p
    xcols <- setdiff(seq_len(p), be$icol)
    px <- length(xcols)
    G <- be$G
    nfe <- length(G)
    mom <- be$moments
    gs <- be$gs
    nk <- lapply(gs, `[[`, "nk")
    B <- lapply(gs, function(g) cbind(g$S[, xcols, drop = FALSE], g$s))
    if (any(vapply(nk, function(v) any(!(v > 0)), NA)))
        stop("internal error: a fixed-effect level without weight", call. = FALSE)
    sys <- dm_fe_system(G, nk, be$crosstab, opts)
    ab <- dm_fe_apply(sys, B)
    XtX <- mom$xtx[xcols, xcols, drop = FALSE]
    Xty <- mom$xty[xcols]
    Mfull <- rbind(cbind(XtX, Xty), c(Xty, mom$yty)) - ab$Q
    Mxx <- Mfull[seq_len(px), seq_len(px), drop = FALSE]
    Mxy <- Mfull[seq_len(px), px + 1L]
    sw <- sum(nk[[1L]])
    sx <- colSums(B[[1L]])[seq_len(px)]
    scale <- diag(XtX) - sx^2 / sw
    pc <- dm_fe_pivchol(Mxx, scale, tol)
    kp <- pc$keep
    b <- rep(NA_real_, px)
    if (length(kp)) b[kp] <- backsolve(pc$R, backsolve(pc$R, Mxy[kp], transpose = TRUE))
    b0 <- b
    b0[is.na(b0)] <- 0
    fe <- lapply(ab$A, function(Ak) drop(Ak[, px + 1L] - Ak[, seq_len(px), drop = FALSE] %*% b0))
    full <- function(b0) {
        bf <- numeric(p)
        bf[xcols] <- b0
        bf
    }
    r <- be$resid(full(b0), fe)
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
        r <- be$resid(full(b0), fe)
    }
    n <- be$n
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
    K <- rank + fe_dof
    bread <- matrix(NA_real_, px, px)
    if (rank) bread[kp, kp] <- chol2inv(pc$R)
    out <- list(b = b, kp = kp, xcols = xcols, bread = bread, rss = r$rss, fe = fe,
                rank = rank, G = G, fe_dof = fe_dof, fe_dof_exact = fe_dof_exact,
                df.residual = n - K, nobs = n, tss = mom$tss,
                tss_within = max(Mfull[px + 1L, px + 1L], 0), info = ab$info)
    if (cluster) {
        cl <- be$cluster()
        U <- cl$U[, xcols, drop = FALSE]
        for (k in which(!cl$nested))
            U <- U - as.matrix(Matrix::crossprod(cl$E(k), ab$A[[k]][, seq_len(px), drop = FALSE]))
        U <- U[, kp, drop = FALSE]
        nn <- !cl$nested
        K_cl <- rank + if (any(nn)) sum(G[nn]) - (sum(nn) - 1L) else 1L
        Gc <- cl$G
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

## ---------------------------------------------------------------------------
## In-memory backend

dm_fe_mem_backend <- function(D, y, w, codes, G, clcodes = NULL, Gc = NULL) {
    n <- nrow(D)
    p <- ncol(D)
    cn <- colnames(D)
    wy <- if (is.null(w)) y else w * y
    xtx <- if (is.null(w)) crossprod(D) else crossprod(D, w * D)
    xty <- drop(crossprod(D, wy))
    sw <- if (is.null(w)) n else sum(w)
    ybar <- sum(wy) / sw
    tss <- if (is.null(w)) sum((y - ybar)^2) else sum(w * (y - ybar)^2)
    mom <- list(xtx = unname(xtx), xty = unname(xty), yty = sum(wy * y), tss = tss)
    gs <- dm_fe_gsums(D, w, y, codes, G)
    state <- new.env(parent = emptyenv())
    be <- list(n = n, p = p, icol = which(cn == "(Intercept)"), G = G, moments = mom, gs = gs)
    be$crosstab <- function(a, b) dm_fe_crosstab_t(codes[[a]], codes[[b]], G[a], G[b], w)
    be$resid <- function(beta, fe) {
        e <- y - drop(D %*% beta)
        for (k in seq_along(codes)) e <- e - fe[[k]][codes[[k]]]
        we <- if (is.null(w)) e else w * e
        state$we <- we
        Be <- dm_fe_gsums(NULL, we, NULL, codes, G)
        list(xte = drop(crossprod(D, we)), Be = lapply(Be, `[[`, "nk"), rss = sum(we * e))
    }
    be$cluster <- function() {
        we <- state$we
        U <- dm_fe_gsums(D, we, NULL, list(clcodes), Gc)[[1L]]$S
        nested <- vapply(seq_along(codes), function(k) {
            cd <- codes[[k]]
            first <- clcodes[match(seq_len(G[k]), cd)]
            all(first[cd] == clcodes)
        }, NA)
        list(G = Gc, U = U, nested = nested,
             E = function(k) dm_fe_crosstab_t(clcodes, codes[[k]], Gc, G[k], we))
    }
    be
}

## For each FE k: list(nk = sums of w by level, S = grouped sums of w X
## (G_k x p; NULL when D is NULL), s = grouped sums of w y (NULL when y is
## NULL)). One pass over the rows with the structured kernel when the
## design allows it, the generic row kernel otherwise.
dm_fe_gsums <- function(D, w, y, codes, G) {
    nfe <- length(codes)
    wd <- if (is.null(w)) NULL else as.double(w)
    if (is.null(D)) {
        M <- .Call(C_dm_fe_struct_gsum, NULL, NULL, wd, codes, as.integer(G), NULL, FALSE)
        return(lapply(M, function(m) list(nk = m[1L, ])))
    }
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
    mu <- vapply(num, function(nm) {
        x <- mf[[nm]]
        if (is.null(w)) mean(x) else sum(w * x) / sum(w)
    }, 0)
    sdv <- vapply(num, function(nm) stats::sd(mf[[nm]]), 0)
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

## Fit on a backend built by make_be(design), centered when exact.
dm_fe_estimate <- function(make_be, D, ctr, D_centered, tol, cluster) {
    opts <- dm_fe_opts()
    if (!is.null(ctr)) {
        be <- make_be(D_centered)
        est <- dm_fe_uncenter(dm_fe_core(be, tol, opts, cluster), ctr$T, be$icol)
        if (!is.null(est)) return(est)
    }
    dm_fe_core(make_be(D), tol, opts, cluster)
}

## ---------------------------------------------------------------------------
## design_lm() with fixed effects

## Entry point from design_lm(): `call` is its matched call, `env` its
## parent frame.
dm_fe_lm <- function(call, env, formula, contrasts, compress, tol, tol_missing, model,
                     vcov, fixef.rm, dots) {
    if (isTRUE(compress)) stop("`compress = TRUE` is not supported with fixed effects.", call. = FALSE)
    if (!is.null(call$offset)) stop("Offsets are not supported with fixed effects.", call. = FALSE)
    checkmate::assert_number(tol, lower = 0, finite = TRUE)
    checkmate::assert_flag(model)
    ## aliasing: fixest's criterion (residual variance below tol times the
    ## variance before the FE), with fixest's default collin.tol = 1e-10
    if (tol_missing) tol <- 1e-10
    pf <- dm_fe_parse(formula)
    clv <- dm_fe_vcov_arg(vcov)
    fixef.rm <- dm_fe_rm_arg(fixef.rm)
    data <- if (is.null(call$data)) NULL else eval(call$data, env)
    if (!is.null(data) && dm_is_remote_source(data)) {
        if (!is.null(call$subset)) stop("`subset` is not supported for remote data.", call. = FALSE)
        if (!is.null(call$na.action) && dm_remote_na(eval(call$na.action, env)) != "omit")
            stop("For remote data, `na.action` must be na.omit.", call. = FALSE)
        w <- call$weights
        if (is.name(w)) w <- as.character(w)
        if (!is.null(w) && !checkmate::test_string(w))
            stop("For remote data, `weights` must name a column.", call. = FALSE)
        fit <- dm_fe_remote(pf, data, w, contrasts, tol, clv, fixef.rm, dots)
    } else {
        if (length(dots)) stop("Unknown arguments: ", paste(names(dots), collapse = ", "), call. = FALSE)
        fit <- dm_fe_memory(pf, call, env, contrasts, tol, clv, fixef.rm, model)
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
    mf <- call[c(1L, match(c("data", "subset", "weights", "na.action"), names(call), 0L))]
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

dm_fe_memory <- function(pf, call, env, contrasts, tol, clv, fixef.rm, model) {
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
    clcodes <- Gc <- NULL
    if (!is.null(clv)) {
        same <- match(clv, names(pf$fe))
        ccd <- if (is.na(same)) dm_fe_code(mf[[clv]]) else cd[[same]]
        clcodes <- ccd$codes
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
    est <- dm_fe_estimate(function(X) dm_fe_mem_backend(X, y, w, codes, G, clcodes, Gc),
                          D, ctr, D2, tol, !is.null(clv))
    fit <- dm_fe_finish(est, colnames(D), pf, lapply(cd, `[[`, "levels"), clv, nsingle, w, nall)
    fit$terms <- mt
    fit$contrasts <- D$contrasts
    fit$xlevels <- stats::.getXlevels(mt, mf)
    fit$spec <- dmfit_spec(D)
    fit$na.action <- dm_fe_na_action(na_act, rows, nall)
    if (model) {
        fit$model <- mf
        fit$data <- list(design = D, y = y, weights = w, codes = codes)
    }
    fit
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
dm_fe_finish <- function(est, cn, pf, lv, clv, nsingle, w, nall) {
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
                weighted = !is.null(w), compressed = FALSE, patterns = est$nobs,
                rows = nall, fe_info = est$info, intercept = TRUE,
                sumlogw = if (is.null(w)) 0 else sum(log(w)))
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

## Sum of the FE values of the rows: codes as a list (one integer vector per
## FE), NA for unseen levels.
dm_fe_sum <- function(fe, codes) {
    out <- 0
    for (k in seq_along(fe)) out <- out + fe[[k]][codes[[k]]]
    out
}

#' @export
residuals.design_lm_fe <- function(object, ...) {
    if (isTRUE(object$remote)) return(dm_fe_remote_residuals(object))
    dat <- dmfit_data(object, "residuals()")
    r <- dat$y - dmfit_linpred_fe(object, dat$design) - dm_fe_sum(object$fixef, dat$codes)
    stats::naresid(object$na.action, as.vector(r))
}

#' @export
fitted.design_lm_fe <- function(object, ...) {
    if (isTRUE(object$remote)) stop("For a remote fit, use residuals(fit): its `.fitted` column.", call. = FALSE)
    dat <- dmfit_data(object, "fitted()")
    f <- dmfit_linpred_fe(object, dat$design) + dm_fe_sum(object$fixef, dat$codes)
    stats::napredict(object$na.action, as.vector(f))
}

## X b for the design columns but the intercept (absorbed).
dmfit_linpred_fe <- function(object, D) {
    b <- object$coefficients
    b[is.na(b)] <- 0
    beta <- numeric(ncol(D))
    beta[match(names(b), colnames(D))] <- b
    if (!ncol(D)) return(rep(0, nrow(D)))
    drop(D %*% beta)
}

## Codes of new data's FE levels in the fit's levels (NA when unseen).
dm_fe_match <- function(object, newdata) {
    lapply(seq_along(object$fe_vars), function(k) {
        vars <- object$fe_vars[[k]]
        miss <- setdiff(vars, names(newdata))
        if (length(miss))
            stop(sprintf("`newdata` has no fixed-effect variable `%s`.", miss[1L]), call. = FALSE)
        lv <- object$fe_levels[[k]]
        nd <- as.data.frame(lapply(vars, function(v) {
            x <- newdata[[v]]
            if (is.factor(x)) as.character(x) else x
        }), stringsAsFactors = FALSE)
        names(nd) <- vars
        match(dm_fe_keys(nd), dm_fe_keys(lv))
    })
}

#' @export
predict.design_lm_fe <- function(object, newdata, se.fit = FALSE,
                                 interval = c("none", "confidence", "prediction"),
                                 level = 0.95, na.action = stats::na.pass, ...) {
    interval <- match.arg(interval)
    if (isTRUE(se.fit) || interval != "none")
        stop("With fixed effects, predict() gives the fit only (no `se.fit` or intervals).", call. = FALSE)
    if (missing(newdata) || is.null(newdata)) {
        if (isTRUE(object$remote)) stop("For a remote fit, give `newdata`, or use residuals(fit).", call. = FALSE)
        return(fitted(object))
    }
    if (dm_is_remote_source(newdata)) return(dm_fe_remote_predict(object, newdata, list(...)))
    if (is.null(object$spec)) stop("This design cannot be applied to new data.", call. = FALSE)
    if (!identical(na.action, stats::na.pass))
        stop("With fixed effects, predict() keeps every row of `newdata` (na.action = na.pass).", call. = FALSE)
    D <- design_matrix(object$spec, newdata, na.action = stats::na.pass)
    f <- dmfit_linpred_fe(object, D) + dm_fe_sum(object$fixef, dm_fe_match(object, newdata))
    names(f) <- dimnames(D)[[1L]]
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

## (remote fits: fe_remote.R)

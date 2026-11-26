## Prototype (internal, not exported): linear models with absorbed fixed
## effects, by Frisch-Waugh-Lovell on aggregates instead of demeaning.
##
## Model: y = X b + F a + e, weights w, F = [F_1 ... F_K] the dummies of
## K fixed effects (FE). With S = F'WX, s = F'Wy and K = F'WF (diagonal
## blocks diag(n_k), off-diagonal blocks the weighted crosstabs N_kl),
##
##   X'MX = X'WX - S' K^- S,   X'My = X'Wy - S' K^- s,
##
## for any generalized inverse K^- (the columns of S and s lie in the range
## of K). Every input is a sum over rows: X'WX and X'Wy are products with
## the design, S and s grouped sums, K counts and crosstabs; the rows are
## never demeaned and the dense design is never formed.
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
## such a component leaves the fit unchanged). The reduced system is solved
## densely when small, by sparse Cholesky (CHOLMOD) otherwise, or by
## Jacobi-preconditioned conjugate gradients with p + 1 right-hand sides.
##
## A = K^- B gives the FE coefficients without another solve: a = A_y - A_x b.
## Residuals per row (one product with the design) give the residual sum of
## squares and the cluster scores. With FE, the score of cluster c is
## sum_{i in c} w_i e_i (MX)_i = colSums(w e X, group = cl)_c
## - sum_k E_ck Gamma_k, Gamma = K^- S and E_ck the crosstab of w e over
## (cluster, FE_k); the correction vanishes for FEs nested in the clusters
## (their normal equations make each level's sum of w e zero).
##
## Degrees of freedom and small sample corrections follow fixest's defaults
## (ssc(K.adj = TRUE, K.fixef = "nonnested", K.exact = FALSE, G.adj = TRUE)),
## and singletons are dropped iteratively as feols() does by default.

## Split `y ~ x | a + b^c` into the main formula and the FE terms (each a
## character vector of variable names, several for an interaction a^b).
dm_fe_parse <- function(formula) {
    rhs <- formula[[length(formula)]]
    if (!is.call(rhs) || !identical(rhs[[1L]], as.name("|")))
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
    list(main = main, fe = fe)
}

## Integer codes 1..G of a vector (no missing values), in order of first
## appearance (factors: levels that occur, in level order).
dm_fe_code <- function(v) {
    if (is.factor(v)) {
        cd <- as.integer(v)
        return(dm_fe_compact(cd, nlevels(v)))
    }
    if (is.numeric(v) && length(v)) {
        ## integer-valued with a small range: offset and compact, no hashing
        r <- range(v)
        if (r[2L] - r[1L] < 4 * length(v) + 1e6 && all(v == round(v)))
            return(dm_fe_compact(as.integer(v - r[1L] + 1), as.integer(r[2L] - r[1L] + 1)))
    }
    cd <- match(v, unique(v))
    list(codes = cd, G = max(0L, cd))
}

## Renumber codes so that the levels that occur are 1..G.
dm_fe_compact <- function(cd, G) {
    used <- tabulate(cd, nbins = G) > 0L
    if (all(used)) return(list(codes = cd, G = G))
    map <- cumsum(used)
    list(codes = map[cd], G = sum(used))
}

## Codes of an interaction of several code vectors.
dm_fe_interact <- function(cl) {
    cd <- cl[[1L]]$codes
    for (k in seq_along(cl)[-1L]) cd <- (cd - 1) * cl[[k]]$G + cl[[k]]$codes
    dm_fe_code(cd)
}

## Rows kept after dropping singletons iteratively (rows alone in a level of
## some FE), as feols() does by default.
dm_fe_singletons <- function(codes) {
    n <- length(codes[[1L]]$codes)
    keep <- rep(TRUE, n)
    repeat {
        drop <- rep(FALSE, n)
        for (cd in codes) {
            cnt <- tabulate(cd$codes[keep], nbins = cd$G)
            drop <- drop | (keep & cnt[cd$codes] == 1L)
        }
        if (!any(drop)) break
        keep <- keep & !drop
    }
    keep
}

## Grouped sums of the columns of a numeric matrix.
dm_fe_gsum <- function(codes, G, W) {
    if (!is.matrix(W)) W <- matrix(W, ncol = 1L)
    storage.mode(W) <- "double"
    .Call(C_dm_fe_gsum, codes, as.integer(G), W)
}

## colSums(w * D, group = codes) for codes 1..G: a G x p matrix.
dm_fe_dsum <- function(D, w, codes, G) {
    J <- as_jac(D)
    W <- matrix(if (is.null(w)) 1 else as.double(w), nrow = nrow(D), ncol = 1L)
    if (jac_kernel_ok(J)) return(t(dm_crossprod_kernel(J, W, codes, G)))
    X <- as.matrix(D)
    storage.mode(X) <- "double"
    dm_fe_gsum(codes, G, W[, 1L] * X)
}

## The weighted crosstab of codes a and b as t(N_ab), a G_b x G_a dgCMatrix.
dm_fe_crosstab_t <- function(a, b, w = NULL) {
    ct <- .Call(C_dm_fe_crosstab, a$codes, b$codes,
                if (is.null(w)) NULL else as.double(w), as.integer(a$G), as.integer(b$G))
    methods::new("dgCMatrix", i = ct[[2L]], p = ct[[1L]], x = ct[[3L]],
                 Dim = c(as.integer(b$G), as.integer(a$G)))
}

## Component of each level of b in the bipartite graph of t(N_ab).
dm_fe_components <- function(Nt) {
    i <- rep.int(seq_len(ncol(Nt)) - 1L, diff(Nt@p))
    .Call(C_dm_fe_components, i, Nt@i, ncol(Nt), nrow(Nt))
}

## Jacobi-preconditioned conjugate gradients for R Z = B, one independent
## recursion per column of B, until every relative residual is below tol.
dm_fe_pcg <- function(R, B, tol = 1e-12, maxit = 10000L) {
    dinv <- 1 / Matrix::diag(R)
    m <- ncol(B)
    X <- matrix(0, nrow(B), m)
    Rr <- B
    bn <- sqrt(colSums(B^2))
    bn[bn == 0] <- 1
    Z <- dinv * Rr
    P <- Z
    rz <- colSums(Rr * Z)
    it <- 0L
    res <- sqrt(colSums(Rr^2)) / bn
    while (max(res) > tol && it < maxit) {
        it <- it + 1L
        AP <- as.matrix(R %*% P)
        pap <- colSums(P * AP)
        alpha <- ifelse(pap > 0, rz / pap, 0)
        X <- X + P * rep(alpha, each = nrow(P))
        Rr <- Rr - AP * rep(alpha, each = nrow(P))
        Z <- dinv * Rr
        rz_new <- colSums(Rr * Z)
        beta <- ifelse(rz > 0, rz_new / rz, 0)
        P <- Z + P * rep(beta, each = nrow(P))
        rz <- rz_new
        res <- sqrt(colSums(Rr^2)) / bn
    }
    list(x = X, iter = it, resid = max(res), converged = max(res) <= tol)
}

## Solve the reduced system R Z = B (R symmetric positive definite after the
## reference levels are removed). Returns Z and solver details.
dm_fe_solve <- function(R, B, solver, pcg_tol, pcg_maxit, dense_max = 3000L) {
    m <- nrow(R)
    if (m == 0L) return(list(z = B, solver = "none"))
    ## sparse Cholesky is never chosen automatically: on worker x firm
    ## graphs its fill is close to dense (see the b3fm experiment)
    if (solver == "auto") solver <- if (m <= dense_max) "dense" else "pcg"
    if (solver == "dense") {
        Rd <- as.matrix(R)
        ch <- tryCatch(chol(Rd), error = function(e) NULL)
        if (!is.null(ch)) {
            z <- backsolve(ch, backsolve(ch, B, transpose = TRUE))
            return(list(z = z, solver = "dense", nnz_L = m * (m + 1) / 2))
        }
        solver <- "pcg"
    }
    if (solver == "chol") {
        L <- tryCatch(Matrix::Cholesky(R, perm = TRUE, LDL = FALSE),
                      error = function(e) NULL)
        if (!is.null(L)) {
            z <- as.matrix(Matrix::solve(L, B, system = "A"))
            ok <- all(is.finite(z)) &&
                max(abs(as.matrix(R %*% z) - B)) <= 1e-8 * max(1, max(abs(B)))
            if (ok) return(list(z = z, solver = "chol", nnz_L = length(L@x)))
        }
        solver <- "pcg"
    }
    sol <- dm_fe_pcg(R, B, pcg_tol, pcg_maxit)
    if (!sol$converged)
        warning(sprintf("PCG did not converge in %d iterations (relative residual %.2g).",
                        sol$iter, sol$resid), call. = FALSE)
    list(z = sol$x, solver = "pcg", iter = sol$iter, resid = sol$resid)
}

## The FE system K = F'WF, reduced by eliminating the largest FE.
## codes: list of list(codes, G); nk: list of weighted counts per level.
dm_fe_system <- function(codes, nk, w, solver = "auto", pcg_tol = 1e-12,
                         pcg_maxit = 10000L) {
    nfe <- length(codes)
    G <- vapply(codes, `[[`, 0L, "G")
    o <- order(G, decreasing = TRUE)
    sys <- list(nfe = nfe, f1 = o[1L], D1 = nk[[o[1L]]], rest = o[-1L],
                solver = solver, pcg_tol = pcg_tol, pcg_maxit = pcg_maxit)
    if (nfe == 1L) return(sys)
    rest <- sys$rest
    sys$Gr <- G[rest]
    sys$off <- off <- c(0L, cumsum(sys$Gr))
    m <- off[length(off)]
    ## t(N_1k) for each other FE, and references from the (FE_1, FE_k) graphs
    Nt <- lapply(rest, function(k) dm_fe_crosstab_t(codes[[sys$f1]], codes[[k]], w))
    ncomp <- integer(length(rest))
    ref <- integer(0)
    for (r in seq_along(rest)) {
        comp <- dm_fe_components(Nt[[r]])
        ncomp[r] <- max(comp)
        ref <- c(ref, off[r] + which(!duplicated(comp)))
    }
    sys$Ntr <- do.call(rbind, Nt)
    ## R = blockdiag(n_k) + offdiag(N_kl) - N_1' D_1^-1 N_1
    R <- Matrix::Diagonal(x = unlist(nk[rest], use.names = FALSE))
    if (length(rest) > 1L) {
        blocks <- list()
        for (r in seq_along(rest)[-length(rest)]) for (s in (r + 1L):length(rest)) {
            Nrs <- dm_fe_crosstab_t(codes[[rest[r]]], codes[[rest[s]]], w) # G_s x G_r
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
    t0 <- proc.time()[[3L]]
    sol <- dm_fe_solve(sys$R, Brr[sys$keep, , drop = FALSE], sys$solver,
                       sys$pcg_tol, sys$pcg_maxit)
    Z <- matrix(0, sys$m, ncol(Brr))
    Z[sys$keep, ] <- sol$z
    sol$z <- NULL
    sol$time_solve <- proc.time()[[3L]] - t0
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

## Evaluate a weights or cluster argument: NULL, a column name, a one-sided
## formula, or a vector with one value per row of data.
dm_fe_arg <- function(x, data, env) {
    if (is.null(x)) return(NULL)
    if (inherits(x, "formula")) return(eval(x[[2L]], data, env))
    if (is.character(x) && length(x) == 1L && nrow(data) != 1L) return(data[[x]])
    x
}

dm_fe_fit <- function(formula, data, cluster = NULL, weights = NULL,
                      singletons = "drop", solver = c("auto", "dense", "chol", "pcg"),
                      tol = 1e-10, pcg_tol = 1e-12, pcg_maxit = 10000L,
                      refine = 0L) {
    solver <- match.arg(solver)
    singletons <- match.arg(singletons, c("drop", "keep"))
    times <- numeric(0)
    tick <- function(name) {
        now <- proc.time()[[3L]]
        times[name] <<- now - t_last
        t_last <<- now
    }
    t_last <- proc.time()[[3L]]
    env <- environment(formula)
    pf <- dm_fe_parse(formula)
    fe_vars <- unique(unlist(pf$fe, use.names = FALSE))
    vars <- unique(c(all.vars(pf$main), fe_vars))
    df <- as.data.frame(data)[, intersect(vars, names(data)), drop = FALSE]
    w <- dm_fe_arg(weights, data, env)
    cl <- dm_fe_arg(cluster, data, env)
    ok <- stats::complete.cases(df)
    if (!is.null(w)) ok <- ok & !is.na(w) & w != 0
    if (!is.null(cl)) ok <- ok & !is.na(cl)
    if (!is.null(w) && any(w[ok] < 0)) stop("Negative weights.", call. = FALSE)
    ## FE codes, singletons
    fe_codes <- function(d) lapply(pf$fe, function(v) {
        if (length(v) == 1L) dm_fe_code(d[[v]]) else dm_fe_interact(lapply(v, function(u) dm_fe_code(d[[u]])))
    })
    rows <- which(ok)
    codes <- fe_codes(if (all(ok)) df else df[rows, fe_vars, drop = FALSE])
    nsingle <- 0L
    if (singletons == "drop") {
        keep <- dm_fe_singletons(codes)
        nsingle <- sum(!keep)
        if (nsingle) {
            rows <- rows[keep]
            codes <- lapply(codes, function(cd) dm_fe_compact(cd$codes[keep], cd$G))
        }
    }
    if (length(rows) < nrow(df)) df <- df[rows, , drop = FALSE]
    if (!is.null(w)) w <- as.double(w[rows])
    if (!is.null(cl)) cl <- cl[rows]
    for (v in names(df)) if (is.factor(df[[v]])) df[[v]] <- droplevels(df[[v]])
    n <- length(rows)
    tick("prepare")
    ## design and response
    y <- as.double(eval(pf$main[[2L]], df, env))
    D <- design_matrix(pf$main[-2L], df)
    cn <- colnames(D)
    ic <- which(cn == "(Intercept)")
    xcols <- setdiff(seq_along(cn), ic)
    p <- length(xcols)
    tick("design")
    ## aggregates: X'WX, X'Wy, y'Wy, grouped sums and counts
    XtX <- if (is.null(w)) crossprod(D) else crossprod(sqrt(w) * D)
    wy <- if (is.null(w)) y else w * y
    Xty <- crossprod(D, wy)[, 1L]
    yty <- sum(wy * y)
    tick("crossprod")
    nk <- lapply(codes, function(cd) {
        if (is.null(w)) as.double(tabulate(cd$codes, nbins = cd$G)) else dm_fe_gsum(cd$codes, cd$G, w)[, 1L]
    })
    B <- lapply(codes, function(cd) {
        S <- dm_fe_dsum(D, w, cd$codes, cd$G)[, xcols, drop = FALSE]
        cbind(S, dm_fe_gsum(cd$codes, cd$G, wy))
    })
    tick("group_sums")
    sys <- dm_fe_system(codes, nk, w, solver, pcg_tol, pcg_maxit)
    ab <- dm_fe_apply(sys, B)
    tick("absorb")
    ## FWL normal equations
    Mfull <- rbind(cbind(XtX[xcols, xcols, drop = FALSE], Xty[xcols]), c(Xty[xcols], yty)) - ab$Q
    Mxx <- Mfull[seq_len(p), seq_len(p), drop = FALSE]
    Mxy <- Mfull[seq_len(p), p + 1L]
    sw <- sum(nk[[1L]])
    sx <- colSums(B[[1L]])[seq_len(p)]
    scale <- diag(XtX)[xcols] - sx^2 / sw
    pc <- dm_fe_pivchol(Mxx, scale, tol)
    kp <- pc$keep
    b <- rep(NA_real_, p)
    names(b) <- cn[xcols]
    if (length(kp)) b[kp] <- backsolve(pc$R, backsolve(pc$R, Mxy[kp], transpose = TRUE))
    b0 <- b
    b0[is.na(b0)] <- 0
    ## FE coefficients and residuals
    fe <- lapply(ab$A, function(Ak) drop(Ak[, p + 1L] - Ak[, seq_len(p), drop = FALSE] %*% b0))
    names(fe) <- names(pf$fe)
    bfull <- numeric(length(cn))
    bfull[xcols] <- b0
    resid <- function() {
        bfull[xcols] <- b0
        e <- y - drop(D %*% bfull)
        for (k in seq_along(codes)) e <- e - fe[[k]][codes[[k]]$codes]
        e
    }
    e <- resid()
    we <- if (is.null(w)) e else w * e
    tick("solve_resid")
    ## iterative refinement: the normal equations X'M e = 0 from the
    ## residuals, which are computed row by row (two-pass accuracy), correct
    ## the rounding of the one-pass aggregates in X'My
    for (it in seq_len(refine)) {
        if (!length(kp)) break
        g <- crossprod(D, we)[xcols, 1L]
        Be <- lapply(codes, function(cd) dm_fe_gsum(cd$codes, cd$G, we))
        Ae <- dm_fe_apply(sys, Be)$A
        r <- g
        for (k in seq_along(codes)) r <- r - drop(crossprod(B[[k]][, seq_len(p), drop = FALSE], Ae[[k]]))
        delta <- backsolve(pc$R, backsolve(pc$R, r[kp], transpose = TRUE))
        b0[kp] <- b0[kp] + delta
        b[kp] <- b0[kp]
        for (k in seq_along(codes))
            fe[[k]] <- fe[[k]] + Ae[[k]][, 1L] - drop(ab$A[[k]][, kp, drop = FALSE] %*% delta)
        e <- resid()
        we <- if (is.null(w)) e else w * e
    }
    rss <- sum(we * e)
    tick("refine")
    ## degrees of freedom (fixest: K.exact = FALSE)
    G <- vapply(codes, `[[`, 0L, "G")
    fe_dof <- sum(G) - (length(G) - 1L)
    ## the exact rank of K (fixest's K.exact = TRUE): with two FE, one level
    ## per connected component is redundant; with more, the count is exact
    ## when the reduced system was positive definite once references were
    ## removed (a dense or sparse Cholesky succeeded)
    fe_dof_exact <- if (length(G) == 1L) {
        G
    } else if (length(G) == 2L || ab$info$solver %in% c("dense", "chol")) {
        max(G) + ab$info$m
    } else {
        NA_integer_
    }
    rank <- length(kp)
    K <- rank + fe_dof
    bread <- if (rank) chol2inv(pc$R) else matrix(0, 0, 0)
    V <- matrix(NA_real_, p, p, dimnames = list(cn[xcols], cn[xcols]))
    V_iid <- V
    V_iid[kp, kp] <- bread * rss / (n - K)
    V_cl <- NULL
    if (!is.null(cl)) {
        ## a cluster variable that is one of the FE: reuse its codes
        same <- if (inherits(cluster, "formula") && is.name(cluster[[2L]])) {
            match(as.character(cluster[[2L]]), names(pf$fe))
        } else if (is.character(cluster) && length(cluster) == 1L) {
            match(cluster, names(pf$fe))
        } else {
            NA
        }
        ccode <- if (is.na(same)) dm_fe_code(cl) else codes[[same]]
        U <- dm_fe_dsum(D, we, ccode$codes, ccode$G)[, xcols, drop = FALSE]
        nested <- vapply(seq_along(codes), function(k) {
            if (identical(k, same)) return(TRUE)
            cd <- codes[[k]]$codes
            first <- ccode$codes[match(seq_len(codes[[k]]$G), cd)]
            all(first[cd] == ccode$codes)
        }, NA)
        for (k in which(!nested)) {
            Et <- dm_fe_crosstab_t(ccode, codes[[k]], we) # G_k x G_cl
            U <- U - as.matrix(Matrix::crossprod(Et, ab$A[[k]][, seq_len(p), drop = FALSE]))
        }
        U <- U[, kp, drop = FALSE]
        Gc <- ccode$G
        K_cl <- rank + if (any(!nested)) sum(G[!nested]) - (sum(!nested) - 1L) else 1L
        adj <- Gc / (Gc - 1) * (n - 1) / (n - K_cl)
        V_cl <- V
        V_cl[kp, kp] <- adj * bread %*% crossprod(U) %*% bread
        tick("cluster")
    }
    list(coefficients = b, vcov_iid = V_iid, vcov_cluster = V_cl, fixef = fe,
         residuals = e, rss = rss, nobs = n, nsingletons = nsingle, rank = rank,
         fe_sizes = G, fe_dof = fe_dof, fe_dof_exact = fe_dof_exact, df.residual = n - K,
         solver = ab$info, times = times)
}

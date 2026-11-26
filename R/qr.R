#' QR decomposition of a design matrix without forming the dense matrix
#'
#' `qr(x)` decomposes a `design_matrix` as `X P = Q R` without forming the
#' dense matrix `X` or the `n` by `p` matrix `Q`. It returns an object of
#' class `"design_qr"`, which the usual base R functions accept: `qr.coef()`,
#' `qr.fitted()`, `qr.resid()`, `qr.qty()`, `qr.qy()`, `qr.Q()`, `qr.R()` and
#' `solve()`. This gives least squares coefficients, fitted values and
#' residuals at the cost of a few passes over the rows.
#'
#' **Method.** `R` is the Cholesky factor of `X'X`, computed by
#' `crossprod(x)` without forming `X`. Columns are pivoted as in
#' `qr(as.matrix(x))` (LINPACK's limited pivoting, used by [lm()]): they
#' are taken in order, and a column whose residual norm is below `tol`
#' times its norm is moved to the end and counts as aliased. From `X'X`, a
#' residual norm is accurate only to about `sqrt(.Machine$double.eps)`
#' times the column norm, which cannot separate an aliased column from one
#' just above `tol`; the residual norms of such columns are computed from
#' the data (one product with `x` per block of columns). `Q` is
#' `X[, pivot] %*% solve(R)`, never formed. Coefficients solve the normal
#' equations and are then refined from the residuals `y - X b`, each step
#' costing one `x %*% b` and one `crossprod(x, r)`. Refinement recovers
#' much of the accuracy that the normal equations lose on badly conditioned
#' designs.
#'
#' **Agreement with `qr()`.** Results agree with those of
#' `qr(as.matrix(x))` up to floating point rounding, not bit for bit. On
#' badly conditioned designs the coefficients agree to fewer digits than
#' with a Householder QR. A column whose residual norm is within rounding
#' of `tol` can be aliased by one decomposition and not the other. Rows of `R` can
#' differ in sign (the diagonal here is positive), and so can the columns
#' of `Q` and the values of `qr.qty()`. For an aliased design, the block of
#' `R` for the aliased columns is zero.
#'
#' **What is not formed.** The complete `Q` (`n` by `n`) is never formed,
#' so:
#' * `qr.Q(q)` returns the `n` by `rank` matrix of the first `rank` columns
#'   of `Q`, which spans the columns of `X`; `complete = TRUE` is not
#'   supported.
#' * `qr.qty(q, y)` returns the first `rank` rows of `t(Q) %*% y`.
#' * `qr.qy(q, y)` takes `y` with `rank` rows and returns `Q[, 1:rank] %*% y`.
#' * `qr.fitted()` only supports `k = q$rank`.
#'
#' Non-finite entries in `x` are an error, as in [qr()]. `y` must be finite.
#' When `X'X` cannot represent the columns' scale (entries whose squares
#' overflow or underflow, around `1e154` and `1e-154`), `qr()` falls back
#' to the dense QR of `qr(as.matrix(x))`, kept in the object, and the
#' methods give its results.
#'
#' **Calling the methods.** `qr.coef()`, `qr.fitted()`, `qr.resid()`,
#' `qr.qty()`, `qr.qy()`, `qr.Q()` and `qr.R()` are S4 generics created by
#' this package (as the Matrix package does for its decompositions): they
#' dispatch when design.matrix is attached, or when called as
#' `design.matrix::qr.coef()`. Calls to `base::qr.coef()` (for example from
#' code in another package that imports base's function) do not dispatch,
#' and stop because the object is not a base QR decomposition. `qr()` and
#' `solve()` are S3 generics and always dispatch.
#'
#' @param x A `design_matrix`.
#' @param tol The tolerance for detecting aliased columns, as in [qr()].
#' @param ... Not used.
#' @return An object of class `"design_qr"`: a list with `design` (`x`),
#'   `R` (the `p` by `p` upper triangular factor, columns in pivoted
#'   order), `rank`, `pivot` and `tol`.
#' @examples
#' set.seed(1)
#' d <- data.frame(y = rnorm(100), x = rnorm(100), f = gl(4, 25))
#' D <- design_matrix(y ~ x * f, d)
#' q <- qr(D)
#' all.equal(qr.coef(q, d$y), coef(lm(y ~ x * f, d)))
#' all.equal(qr.fitted(q, d$y), unname(fitted(lm(y ~ x * f, d))))
#' @name design_qr
NULL

#' @rdname design_qr
#' @method qr design_matrix
#' @export
qr.design_matrix <- function(x, tol = 1e-07, ...) {
    checkmate::assert_number(tol, lower = 0, finite = TRUE)
    dots <- list(...)
    if (isTRUE(dots$LAPACK))
        stop("`LAPACK = TRUE` is not supported for design matrices.", call. = FALSE)
    G <- crossprod(x)
    storage.mode(G) <- "double"
    if (!all(is.finite(G))) {
        if (!all(finite_rows(x))) stop("NA/NaN/Inf in 'x'", call. = FALSE)
        ## finite entries whose squares overflow
        return(dmqr_dense(x, tol))
    }
    dg <- diag(G)
    ## squares near the limits of double precision lose the column's scale
    if (any(dg > 1e280 | (dg > 0 & dg < 1e-280))) return(dmqr_dense(x, tol))
    ch <- dmqr_chol(x, G, tol)
    if (is.null(ch)) return(dmqr_dense(x, tol))
    R <- ch$R
    cn <- colnames(x)
    if (!is.null(cn)) colnames(R) <- cn[ch$pivot]
    structure(list(design = x, R = R, rank = ch$rank, pivot = ch$pivot, tol = tol),
              class = "design_qr")
}

## The same factorization for a remote design, from the primitives (see
## primitives.R): X'X in the database, and the residual norms of suspect
## columns from products with the data. There is no dense fallback.
#' @method qr dm_remote
#' @export
qr.dm_remote <- function(x, tol = 1e-07, ...) dmqr_remote(x, tol)

## (`Y`, remote columns: X'Y comes from the scan that computes X'X, and is
## kept for dmqr_coef(), see dm_remote_gram_xty())
##
## The columns kept and their order come from X'X (lm()'s rule, as in
## memory), but the values of the factor of the kept columns come from the
## X'X of the design with numeric variables shifted by their mean in a
## sample, when that design is X~ = X T with T unit upper triangular
## (dm_center()): R = R~ T^-1. A numeric variable far from zero (x + 1e7)
## makes X'X nearly singular, while the shifted one is well conditioned; the
## database sums X'X with another rounding than the in-memory kernels, so
## without the shift the factor would lose the accuracy that in-memory
## designs keep. The sums are then those of the shifted design only (one
## scan), and X'X = T^-T (X~'X~) T^-1 follows from them.
dmqr_remote <- function(x, tol, Y = NULL) {
    checkmate::assert_number(tol, lower = 0, finite = TRUE)
    ctr <- if (isTRUE(getOption("design.matrix.remote_center", TRUE))) dm_center(x)
    mc <- xs <- NULL
    if (!is.null(ctr)) {
        xs <- dm_shifted(x, ctr$shift)
        ## (not compressed: the design of the patterns would not be shifted)
        op <- options(design.matrix.remote_compress = FALSE)
        mc <- tryCatch(dm_remote_gram_xty(xs, Y), finally = options(op))
        if (!all(is.finite(mc$xtx))) mc <- NULL
    }
    if (is.null(mc)) {
        mo <- dm_remote_gram_xty(x, Y)
    } else {
        Ti <- backsolve(ctr$T, diag(ncol(ctr$T)))
        G <- crossprod(Ti, mc$xtx %*% Ti)
        mo <- list(xtx = (G + t(G)) / 2, xty = if (!is.null(mc$xty)) crossprod(Ti, mc$xty))
        dimnames(mo$xtx) <- dimnames(mc$xtx)
    }
    G <- mo$xtx
    storage.mode(G) <- "double"
    if (!all(is.finite(G))) stop("NA/NaN/Inf in 'x'", call. = FALSE)
    ch <- dmqr_chol(x, G, tol)
    if (is.null(ch)) stop("The squares of some entries of the design underflow.", call. = FALSE)
    if (!is.null(mc)) ch <- dmqr_centered(ch, G, ctr$T, xs, mc, Y)
    R <- ch$R
    cn <- colnames(x)
    if (!is.null(cn)) colnames(R) <- cn[ch$pivot]
    structure(list(design = x, R = R, rank = ch$rank, pivot = ch$pivot, tol = tol, G = unname(G),
                   centered = ch$centered, xty = mo$xty, xty_of = if (!is.null(mo$xty)) unclass(Y)$sql),
              class = "design_qr")
}

## The factorization of a design for a least squares fit of the columns Y:
## qr(x), which remotely also computes X'Y in the same scan.
dmqr_for <- function(x, tol, Y) {
    if (inherits(x, "dm_remote") && dm_is_rcols(Y)) dmqr_remote(x, tol, Y) else qr(x, tol = tol)
}

## X'Y for the columns Y, from the scan of the factorization when it was
## made for them (dmqr_for()), else computed.
dmqr_xty <- function(q, D, Y) {
    if (!is.null(q$xty) && dm_is_rcols(Y) && identical(q$xty_of, unclass(Y)$sql)) return(q$xty)
    dmp_crossprod(D, Y)
}

## The factor of a remote design with the values of the shifted design
## (see dmqr_remote()): `xs` the shifted design, `mc` its X'X and X'Y, `Tm`
## with X~ = X T. Unchanged when the kept columns do not allow it.
dmqr_centered <- function(ch, G, Tm, xs, mc, Y = NULL) {
    if (ch$rank == 0L) return(ch)
    k <- ch$rank
    p <- length(ch$pivot)
    keep <- ch$pivot[seq_len(k)]
    al <- ch$pivot[seq_len(p) > k]
    if (!identical(keep, sort(keep)) || (length(al) && any(Tm[al, keep] != 0))) return(ch)
    Rt <- tryCatch(chol(unname(mc$xtx[keep, keep, drop = FALSE])), error = function(e) NULL)
    if (is.null(Rt)) return(ch)
    R11 <- Rt %*% backsolve(Tm[keep, keep, drop = FALSE], diag(k))
    R <- ch$R
    R[seq_len(k), seq_len(k)] <- R11
    if (length(al)) {
        R[seq_len(k), k + seq_along(al)] <- backsolve(R11, unname(G[keep, al, drop = FALSE]), transpose = TRUE)
    }
    ch$R <- R
    ch$centered <- list(design = xs, T = Tm[keep, keep, drop = FALSE], Rt = Rt, xty = mc$xty,
                        xty_of = if (!is.null(mc$xty)) unclass(Y)$sql)
    ch
}

## The remote design with numeric variables shifted (their SQL, or the
## batches of the chunked engine).
dm_shifted <- function(x, shift) {
    if (x$engine == "sql") {
        for (nm in names(shift))
            x$vars[[nm]]$sql <- sprintf("(%s - %s)", x$vars[[nm]]$sql, dm_sql_double(shift[[nm]]))
    } else {
        x$shift <- shift
    }
    x
}

## Centering of the numeric variables of a remote design: list(shift (named
## by variable), T) with X(shifted) = X T, T unit upper triangular, or NULL
## (nothing to center, or no such T); from a sample of its rows.
dm_center <- function(D) {
    ## (once for a design: its scaled copies share the cache)
    if (is.null(D$cache$center)) D$cache$center <- list(dm_center_sample(D))
    D$cache$center[[1L]]
}

dm_center_sample <- function(D) {
    plan <- D$plan
    sample <- dm_remote_sample(D)
    if (!nrow(sample)) return(NULL)
    mf <- tryCatch(design_frame(plan$terms, sample, xlev = plan$xlevels, na.action = stats::na.omit),
                   error = function(e) NULL)
    if (is.null(mf) || !nrow(mf)) return(NULL)
    if (D$p > 2000L) return(NULL)
    dm_center_shift(mf, mf, function(fr) dm_frame_design(D, fr), D$p)
}

## The shift of the numeric variables of a model frame whose mean (in the
## rows of `stats_frame`) is large compared with their spread, and T with
## X(shifted) = X T unit upper triangular, checked on prototype rows drawn
## from `proto_frame` (levels and logical values at random, numeric values
## near zero, which keeps the prototype design well conditioned), the
## designs built by build(frame) with p columns: list(shift, T), or NULL.
## (For remote designs, dm_center(); in memory, dm_fe_center().)
dm_center_shift <- function(stats_frame, proto_frame, build, p) {
    tt <- attr(stats_frame, "terms")
    resp <- if (!is.null(tt) && attr(tt, "response") > 0L) dm_terms_vars(tt)$names[attr(tt, "response")]
    cand <- setdiff(names(stats_frame), resp)
    num <- cand[vapply(cand, function(nm) identical(stats::.MFclass(stats_frame[[nm]]), "numeric"), NA)]
    if (!length(num)) return(NULL)
    mu <- vapply(num, function(nm) mean(stats_frame[[nm]]), 0)
    sdv <- vapply(num, function(nm) stats::sd(stats_frame[[nm]]), 0)
    sdv[!is.finite(sdv)] <- 0
    ## only variables whose mean is large compared with their spread
    big <- is.finite(mu) & abs(mu) > sdv
    if (!any(big)) return(NULL)
    shift <- mu[big]
    m <- max(4L * p, 200L)
    proto <- proto_frame[sample.int(nrow(proto_frame), m, replace = TRUE), , drop = FALSE]
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
    X0 <- tryCatch(as.matrix(build(proto)), error = function(e) NULL)
    X1 <- tryCatch(as.matrix(build(shifted)), error = function(e) NULL)
    if (is.null(X0) || is.null(X1) || ncol(X0) != p || ncol(X1) != p) return(NULL)
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

## The in-memory design of a model frame of the plan's terms.
dm_frame_design <- function(D, mf) {
    plan <- D$plan
    d <- design_matrix(plan$terms, mf, contrasts.arg = plan$contrasts)
    if (isTRUE(plan$drop_intercept)) d <- dm_drop_intercept(d)
    d
}


## A design whose X'X loses its scale: the dense QR of qr(), kept in the
## object, so that every method gives the dense results.
dmqr_dense <- function(x, tol) {
    q <- base::qr(as.matrix(x), tol = tol)
    p <- ncol(x)
    R <- matrix(0, p, p)
    Rd <- base::qr.R(q)
    R[seq_len(nrow(Rd)), ] <- Rd
    cn <- colnames(x)
    if (!is.null(cn)) colnames(R) <- cn[q$pivot]
    structure(list(design = x, R = R, rank = q$rank, pivot = q$pivot, tol = tol, dense = q),
              class = "design_qr")
}

## Pivoted Cholesky factor of G = X'X (see dm_chol_pivot() in kernels.c).
## Columns whose residual norm is too small to judge from G, or whose norm
## in G is zero, are suspects: their residual norms are computed from the
## data, and the factorization is repeated until every suspect was judged
## against the columns actually taken before it. Returns NULL when a column
## with a zero norm in G is not zero in the data (its square underflowed).
dmqr_chol <- function(x, G, tol) {
    p <- ncol(G)
    ## columns taken from G alone are clearly above tol
    susp <- max(1e-8, 100 * tol^2)
    dg <- diag(G)
    dres <- rep(NA_real_, p)
    ## A zero norm in G is a zero column unless the column's non-zero
    ## entries can be so small that their squares underflow: when the
    ## design's smallest possible non-zero entry rules that out, zero
    ## columns need no check against the data.
    if (any(dg == 0) && dmp_min_nonzero(x)^2 > 1e-300) dres[dg == 0] <- 0
    repeat {
        ch <- .Call(C_dm_chol_pivot, G, tol, susp, as.double(nrow(x)), dres)
        names(ch) <- c("R", "rank", "pivot", "suspects")
        new <- ch$suspects[is.na(dres[ch$suspects])]
        if (length(new) == 0L) return(ch)
        v <- dmqr_resid2(x, G, ch, new)
        if (any(v$nonzero & !(dg[new] > 0))) return(NULL)
        ## the C code takes a column when dres >= tol^2 * G[j, j], that is,
        ## when its residual norm is at least tol times its norm
        dres[new] <- ifelse(v$nonzero, v$ratio^2 * dg[new], 0)
        ## as dqrdc2, a tie with tol keeps the column
        taken <- new[v$nonzero & v$ratio > 0 & v$ratio >= tol]
        if (length(taken) == 0L) return(ch)
        ## a suspect turned out to be independent: every later column,
        ## suspects included, was judged without it
        dres[seq_len(p) > min(taken)] <- NA_real_
    }
}

## For columns `cols` of x: the norm of the residual after projection on
## the columns taken before each of them, relative to the column's norm,
## both computed from the data (one refinement step) and scaled so that
## neither overflows nor underflows; and whether the column has a non-zero
## entry. In blocks of columns.
dmqr_resid2 <- function(x, G, ch, cols) {
    acc <- ch$pivot[seq_len(ch$rank)]
    R <- ch$R
    p <- ncol(G)
    ratio <- numeric(length(cols))
    nonzero <- logical(length(cols))
    for (blk in split(seq_along(cols), (seq_along(cols) - 1L) %/% 32L)) {
        k <- length(blk)
        M <- matrix(0, p, k)
        E <- matrix(0, p, k)
        A <- vector("list", k)
        for (u in seq_len(k)) {
            j <- cols[blk[u]]
            a <- sum(acc < j)
            A[[u]] <- acc[seq_len(a)]
            Ra <- R[seq_len(a), seq_len(a), drop = FALSE]
            cf <- dmqr_bs(Ra, dmqr_bs(Ra, G[A[[u]], j, drop = FALSE], transpose = TRUE))
            M[j, u] <- 1
            M[A[[u]], u] <- -cf
            E[j, u] <- 1
        }
        res <- dmp_mult(x, M)
        h <- crossprod(x, res)
        for (u in seq_len(k)) {
            a <- length(A[[u]])
            Ra <- R[seq_len(a), seq_len(a), drop = FALSE]
            M[A[[u]], u] <- M[A[[u]], u] -
                dmqr_bs(Ra, dmqr_bs(Ra, h[A[[u]], u, drop = FALSE], transpose = TRUE))
        }
        ## a column's norm is sqrt(G[j, j]), accurate to rounding, unless
        ## G[j, j] is zero (possibly underflow): then it comes from the data
        dz <- diag(G)[cols[blk]] == 0
        nrm <- dmp_colnorms(dmp_mult(x, cbind(M, E[, dz, drop = FALSE])))
        res_n <- nrm[seq_len(k)]
        col_n <- sqrt(diag(G)[cols[blk]])
        col_n[dz] <- nrm[k + seq_len(sum(dz))]
        nonzero[blk] <- col_n > 0
        ratio[blk] <- ifelse(col_n > 0, res_n / col_n, 0)
    }
    list(ratio = ratio, nonzero = nonzero)
}

#' @export
print.design_qr <- function(x, ...) {
    cat(sprintf("<design_qr: %d x %d design, rank %d>\n",
                nrow(x$design), ncol(x$design), x$rank))
    invisible(x)
}

## y as a finite double matrix with `rows` rows; remembers whether it was a
## vector.
dmqr_y <- function(y, rows, what = "y") {
    vec <- !is.matrix(y)
    Y <- if (vec) matrix(y, ncol = 1L) else y
    if (!is.numeric(Y) && !is.logical(Y))
        stop(sprintf("`%s` must be numeric.", what), call. = FALSE)
    if (nrow(Y) != rows)
        stop(sprintf("`%s` must have %d rows.", what, rows), call. = FALSE)
    if (!all(is.finite(Y)))
        stop(sprintf("`%s` must be finite.", what), call. = FALSE)
    storage.mode(Y) <- "double"
    list(Y = unname(Y), vec = vec,
         rows = if (vec) names(y) else rownames(y),
         cols = if (vec) NULL else colnames(y))
}

## R11 (the factor of the first `rank` pivoted columns) and triangular
## solves with it, which also work at rank 0.
dmqr_R11 <- function(qr) {
    k <- qr$rank
    unname(qr$R[seq_len(k), seq_len(k), drop = FALSE])
}

dmqr_bs <- function(R11, rhs, transpose = FALSE) {
    if (nrow(R11) == 0L) return(matrix(0, 0L, ncol(rhs)))
    backsolve(R11, rhs, transpose = transpose)
}

## Solve R11' R11 b = rhs.
dmqr_solve <- function(qr, rhs) {
    R11 <- dmqr_R11(qr)
    dmqr_bs(R11, dmqr_bs(R11, rhs, transpose = TRUE))
}

## Coefficients of the full design (zero for aliased columns, original
## column order) from those of the first `rank` pivoted columns.
dmqr_full <- function(qr, b, fill = 0) {
    p <- ncol(qr$design)
    out <- matrix(fill, p, ncol(b))
    out[qr$pivot[seq_len(qr$rank)], ] <- b
    out
}

## Least squares coefficients of the first `rank` pivoted columns, with
## iterative refinement on the residuals. Each step costs one product and
## one crossproduct with the design. Refinement stops once the correction
## is at the level of rounding, or stops shrinking (the attainable accuracy
## is reached), and after at most `steps` steps.
dmqr_coef <- function(qr, Y, steps = 3L) {
    k <- qr$rank
    D <- qr$design
    keep <- qr$pivot[seq_len(k)]
    if (!is.null(qr$dense))
        return(base::qr.coef(qr$dense, as.matrix(Y))[keep, , drop = FALSE])
    if (!is.null(qr$centered)) {
        ## a remote design with numeric variables shifted (dmqr_centered()):
        ## X~'Y on the shifted design, solved with its well-conditioned
        ## factor, and mapped back (beta = T beta~)
        cc <- qr$centered
        op <- options(design.matrix.remote_compress = FALSE)
        on.exit(options(op))
        ht <- dmqr_xty(cc, cc$design, Y)[keep, , drop = FALSE]
        bt <- backsolve(cc$Rt, backsolve(cc$Rt, ht, transpose = TRUE))
        b <- cc$T %*% bt
        ## the moments of the solution, for sums of squares (dmfit_ls())
        attr(b, "moments") <- list(b = bt, h = ht, G = crossprod(cc$Rt))
        return(b)
    }
    h <- dmqr_xty(qr, D, Y)
    b <- dmqr_solve(qr, h[keep, , drop = FALSE])
    if (k == 0L) return(b)
    ## The normal equations lose about kappa^2 * eps of accuracy, where
    ## kappa is the condition number of X with columns scaled to unit norm
    ## (Cholesky is insensitive to column scaling). When that is already at
    ## the level of rounding, refinement cannot help.
    R11 <- dmqr_R11(qr)
    kap <- 1 / rcond(R11 / rep(sqrt(colSums(R11^2)), each = k), triangular = TRUE)
    if (kap^2 <= 64) {
        if (!is.null(qr$G)) attr(b, "moments") <- list(b = b, h = h[keep, , drop = FALSE], G = qr$G[keep, keep, drop = FALSE])
        return(b)
    }
    last <- Inf
    for (s in seq_len(steps)) {
        ## X'(Y - X b): from the data in memory; for a remote design, from
        ## the moments already computed (X'Y - X'X b), which refines the
        ## solution of the normal equations without another pass
        g <- if (!is.null(qr$G)) h - qr$G %*% dmqr_full(qr, b)
             else dmp_crossprod(D, Y - dmp_mult(D, dmqr_full(qr, b), along = list(Y)))
        db <- dmqr_solve(qr, g[keep, , drop = FALSE])
        b <- b + db
        size <- max(abs(db)) / max(abs(b), .Machine$double.xmin)
        if (!(size > 64 * .Machine$double.eps) || size > 0.25 * last) break
        last <- size
    }
    if (!is.null(qr$G)) attr(b, "moments") <- list(b = b, h = h[keep, , drop = FALSE], G = qr$G[keep, keep, drop = FALSE])
    b
}

## The result as a vector or a matrix like `y`, with row names `rows` and
## (for a matrix) the column names of `y`.
dmqr_out <- function(out, yy, rows = yy$rows) {
    if (yy$vec) {
        out <- out[, 1L]
        names(out) <- rows
    } else {
        dimnames(out) <- if (is.null(rows) && is.null(yy$cols)) NULL else list(rows, yy$cols)
    }
    out
}

## A row-wise result of the methods (the columns `cols`, one or several):
## a matrix or vector like `y` in memory (dmqr_out()), a table with columns
## `prefix` (or prefix1, prefix2, ...) remotely (dmp_output()).
dmqr_result <- function(qr, cols, prefix, yy, rows = yy$rows) {
    k <- dmc_ncol(cols)
    nm <- if (k == 1L) prefix else paste0(prefix, seq_len(k))
    dmp_output(qr$design, stats::setNames(lapply(seq_len(k), function(j) dmc_col(cols, j)), nm),
               memory = function(d) dmqr_out(unname(as.matrix(d)), yy, rows))
}

## The fitted values of y (a column, one or several).
dmqr_fitted <- function(qr, Y) {
    if (!is.null(qr$dense)) return(unname(base::qr.fitted(qr$dense, Y)))
    dmp_mult(qr$design, dmqr_full(qr, dmqr_coef(qr, Y)), along = list(Y))
}

#' @importFrom methods setOldClass
NULL
setOldClass("design_qr")

#' @rdname design_qr
#' @param qr A `"design_qr"` object, from `qr()` on a `design_matrix`.
#' @param y A numeric vector or matrix: the response for `qr.coef()`,
#'   `qr.fitted()`, `qr.resid()` and `qr.qty()` (`nrow(x)` rows), or the
#'   values to multiply by `Q` for `qr.qy()` (`rank` rows).
#' @param k The rank for `qr.fitted()`; only `qr$rank` is supported.
#' @param complete,Dvec Not supported; only the defaults are accepted.
#' @param a,b For `solve()`: a `"design_qr"` object and the right-hand side,
#'   as in [solve.qr()].
#' @exportMethod qr.coef
setMethod("qr.coef", "design_qr", function(qr, y) {
    yy <- dmp_y(qr$design, y)
    b <- dmqr_coef(qr, yy$Y)
    out <- dmqr_full(qr, b, fill = NA_real_)
    dmqr_out(out, yy, colnames(qr$design))
})

#' @rdname design_qr
#' @exportMethod qr.fitted
setMethod("qr.fitted", "design_qr", function(qr, y, k = qr$rank) {
    if (!identical(as.integer(k), as.integer(qr$rank)))
        stop("Only `k = qr$rank` is supported for design matrices.", call. = FALSE)
    yy <- dmp_y(qr$design, y)
    dmqr_result(qr, dmqr_fitted(qr, yy$Y), ".fitted", yy)
})

#' @rdname design_qr
#' @exportMethod qr.resid
setMethod("qr.resid", "design_qr", function(qr, y) {
    yy <- dmp_y(qr$design, y)
    dmqr_result(qr, yy$Y - dmqr_fitted(qr, yy$Y), ".resid", yy)
})

#' @rdname design_qr
#' @exportMethod qr.qty
setMethod("qr.qty", "design_qr", function(qr, y) {
    k <- qr$rank
    keep <- qr$pivot[seq_len(k)]
    yy <- dmp_y(qr$design, y)
    if (!is.null(qr$dense))
        return(dmqr_out(base::qr.qty(qr$dense, yy$Y)[seq_len(k), , drop = FALSE], yy,
                        yy$rows[seq_len(k)]))
    out <- dmqr_bs(dmqr_R11(qr), dmp_crossprod(qr$design, yy$Y)[keep, , drop = FALSE],
                   transpose = TRUE)
    dmqr_out(out, yy, yy$rows[seq_len(k)])
})

#' @rdname design_qr
#' @exportMethod qr.qy
setMethod("qr.qy", "design_qr", function(qr, y) {
    k <- qr$rank
    yy <- dmqr_y(y, k)
    if (!is.null(qr$dense)) {
        Z <- matrix(0, nrow(qr$design), ncol(yy$Y))
        Z[seq_len(k), ] <- yy$Y
        return(dmqr_out(unname(base::qr.qy(qr$dense, Z)), yy, NULL))
    }
    b <- dmqr_bs(dmqr_R11(qr), yy$Y)
    dmqr_result(qr, dmp_mult(qr$design, dmqr_full(qr, b)), ".qy", yy, NULL)
})

#' @rdname design_qr
#' @exportMethod qr.Q
setMethod("qr.Q", "design_qr", function(qr, complete = FALSE, Dvec) {
    if (!isFALSE(complete))
        stop("`complete = TRUE` is not supported for design matrices.", call. = FALSE)
    if (!missing(Dvec))
        stop("`Dvec` is not supported for design matrices.", call. = FALSE)
    k <- qr$rank
    if (!is.null(qr$dense)) return(unname(base::qr.Q(qr$dense)[, seq_len(k), drop = FALSE]))
    dmqr_result(qr, dmp_mult(qr$design, dmqr_full(qr, dmqr_bs(dmqr_R11(qr), diag(1, k)))), ".Q",
                list(vec = FALSE, cols = NULL), NULL)
})

#' @rdname design_qr
#' @exportMethod qr.R
setMethod("qr.R", "design_qr", function(qr, complete = FALSE) {
    if (!isFALSE(complete))
        stop("`complete = TRUE` is not supported for design matrices.", call. = FALSE)
    qr$R
})

#' @rdname design_qr
#' @method solve design_qr
#' @export
solve.design_qr <- function(a, b, ...) {
    n <- nrow(a$design)
    p <- ncol(a$design)
    if (a$rank != min(n, p)) stop("singular matrix 'a' in solve", call. = FALSE)
    if (missing(b)) {
        if (n != p) stop("only square matrices can be inverted", call. = FALSE)
        b <- diag(1, n)
    }
    ## as solve.qr(): coefficients outside the basis are zero, so that the
    ## result solves the system
    res <- qr.coef(a, b)
    res[is.na(res)] <- 0
    res
}

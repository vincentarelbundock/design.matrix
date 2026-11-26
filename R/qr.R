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
    if (!all(is.finite(G))) stop("NA/NaN/Inf in 'x'", call. = FALSE)
    storage.mode(G) <- "double"
    ch <- dmqr_chol(x, G, tol)
    R <- ch$R
    cn <- colnames(x)
    if (!is.null(cn)) colnames(R) <- cn[ch$pivot]
    structure(list(design = x, R = R, rank = ch$rank, pivot = ch$pivot, tol = tol),
              class = "design_qr")
}

## Pivoted Cholesky factor of G = X'X (see dm_chol_pivot() in kernels.c).
## Columns whose residual norm is too small to judge from G are suspects;
## their residual norms are computed from the data, and the factorization
## is repeated until every suspect was judged against the columns actually
## taken before it.
dmqr_chol <- function(x, G, tol, susp = 1e-8) {
    p <- ncol(G)
    dres <- rep(NA_real_, p)
    repeat {
        ch <- .Call(C_dm_chol_pivot, G, tol, susp, as.double(nrow(x)), dres)
        names(ch) <- c("R", "rank", "pivot", "suspects")
        new <- ch$suspects[is.na(dres[ch$suspects])]
        if (length(new) == 0L) return(ch)
        dres[new] <- dmqr_resid2(x, G, ch, new)
        taken <- new[dres[new] > tol^2 * diag(G)[new]]
        if (length(taken) == 0L) return(ch)
        ## a suspect turned out to be independent: later columns were
        ## judged without it
        later <- seq_len(p) > min(taken)
        later[taken] <- FALSE
        dres[later] <- NA_real_
    }
}

## Squared residual norms of columns `cols` of x, after projection on the
## columns taken before each of them, computed from the data (one
## refinement step), in blocks of columns.
dmqr_resid2 <- function(x, G, ch, cols) {
    acc <- ch$pivot[seq_len(ch$rank)]
    R <- ch$R
    out <- numeric(length(cols))
    for (blk in split(seq_along(cols), (seq_along(cols) - 1L) %/% 32L)) {
        M <- matrix(0, ncol(G), length(blk))
        A <- vector("list", length(blk))
        for (u in seq_along(blk)) {
            j <- cols[blk[u]]
            a <- sum(acc < j)
            A[[u]] <- acc[seq_len(a)]
            Ra <- R[seq_len(a), seq_len(a), drop = FALSE]
            cf <- dmqr_bs(Ra, dmqr_bs(Ra, G[A[[u]], j, drop = FALSE], transpose = TRUE))
            M[j, u] <- 1
            M[A[[u]], u] <- -cf
        }
        res <- x %*% M
        h <- crossprod(x, res)
        for (u in seq_along(blk)) {
            a <- length(A[[u]])
            Ra <- R[seq_len(a), seq_len(a), drop = FALSE]
            M[A[[u]], u] <- M[A[[u]], u] -
                dmqr_bs(Ra, dmqr_bs(Ra, h[A[[u]], u, drop = FALSE], transpose = TRUE))
        }
        res <- x %*% M
        out[blk] <- colSums(res^2)
    }
    out
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
    b <- dmqr_solve(qr, crossprod(D, Y)[keep, , drop = FALSE])
    if (k == 0L) return(b)
    ## The normal equations lose about kappa^2 * eps of accuracy, where
    ## kappa is the condition number of X with columns scaled to unit norm
    ## (Cholesky is insensitive to column scaling). When that is already at
    ## the level of rounding, refinement cannot help.
    R11 <- dmqr_R11(qr)
    kap <- 1 / rcond(R11 / rep(sqrt(colSums(R11^2)), each = k), triangular = TRUE)
    if (kap^2 <= 64) return(b)
    last <- Inf
    for (s in seq_len(steps)) {
        r <- Y - D %*% dmqr_full(qr, b)
        db <- dmqr_solve(qr, crossprod(D, r)[keep, , drop = FALSE])
        b <- b + db
        size <- max(abs(db)) / max(abs(b), .Machine$double.xmin)
        if (!(size > 64 * .Machine$double.eps) || size > 0.25 * last) break
        last <- size
    }
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

dmqr_fitted <- function(qr, y) {
    yy <- dmqr_y(y, nrow(qr$design))
    b <- dmqr_coef(qr, yy$Y)
    list(fit = unname(qr$design %*% dmqr_full(qr, b)), yy = yy)
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
    yy <- dmqr_y(y, nrow(qr$design))
    b <- dmqr_coef(qr, yy$Y)
    out <- dmqr_full(qr, b, fill = NA_real_)
    dmqr_out(out, yy, colnames(qr$design))
})

#' @rdname design_qr
#' @exportMethod qr.fitted
setMethod("qr.fitted", "design_qr", function(qr, y, k = qr$rank) {
    if (!identical(as.integer(k), as.integer(qr$rank)))
        stop("Only `k = qr$rank` is supported for design matrices.", call. = FALSE)
    f <- dmqr_fitted(qr, y)
    dmqr_out(f$fit, f$yy)
})

#' @rdname design_qr
#' @exportMethod qr.resid
setMethod("qr.resid", "design_qr", function(qr, y) {
    f <- dmqr_fitted(qr, y)
    dmqr_out(f$yy$Y - f$fit, f$yy)
})

#' @rdname design_qr
#' @exportMethod qr.qty
setMethod("qr.qty", "design_qr", function(qr, y) {
    yy <- dmqr_y(y, nrow(qr$design))
    k <- qr$rank
    keep <- qr$pivot[seq_len(k)]
    out <- dmqr_bs(dmqr_R11(qr), crossprod(qr$design, yy$Y)[keep, , drop = FALSE],
                   transpose = TRUE)
    dmqr_out(out, yy, yy$rows[seq_len(k)])
})

#' @rdname design_qr
#' @exportMethod qr.qy
setMethod("qr.qy", "design_qr", function(qr, y) {
    k <- qr$rank
    yy <- dmqr_y(y, k)
    b <- dmqr_bs(dmqr_R11(qr), yy$Y)
    out <- unname(qr$design %*% dmqr_full(qr, b))
    dmqr_out(out, yy, NULL)
})

#' @rdname design_qr
#' @exportMethod qr.Q
setMethod("qr.Q", "design_qr", function(qr, complete = FALSE, Dvec) {
    if (!isFALSE(complete))
        stop("`complete = TRUE` is not supported for design matrices.", call. = FALSE)
    if (!missing(Dvec))
        stop("`Dvec` is not supported for design matrices.", call. = FALSE)
    k <- qr$rank
    unname(qr$design %*% dmqr_full(qr, dmqr_bs(dmqr_R11(qr), diag(1, k))))
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
    qr.coef(a, b)
}

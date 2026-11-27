## Fused primitives: each computes, from columns, several sums in one call.
## In memory, on plain double vectors (no attributes), one pass in C
## (src/fused.c, whose arithmetic is that of the R code below, to the last
## bit); otherwise (other vectors in memory, remote columns) the operations
## of primitives.R, which build SQL remotely. The fused primitives of the
## fits themselves (residual sums, the working weights of IRLS, deviances)
## are in xtx.fit (its R/kernels.R).

## a plain double vector: the kernels read it directly, without coercion
xtx_plain <- function(x) is.double(x) && is.null(attributes(x))

#' Moments of a design and a response
#'
#' `design_moments()` computes the weighted moments of a least squares
#' problem, `X'WX`, `X'Wz` and `z'Wz`, and the total sum of squares of `z`
#' around its weighted mean (computed in a separate pass around that mean,
#' which keeps its accuracy when the mean is large next to the variation).
#' In memory these are products with the design and one pass over `z`; for
#' a remote design, one scan for the moments and one for the total sum of
#' squares.
#'
#' @param x A design: a `design_matrix` or a remote design, unscaled.
#' @param z The response: a numeric vector with one value per row of `x` in
#'   memory, a [design_column()] of `x` remotely.
#' @param w Weights, in the same form as `z`, or `NULL`.
#' @param n The number of rows (default: those of `x`).
#' @param centered Compute the total sum of squares around the mean
#'   (`FALSE`: the moments only, for a remote design).
#' @return A list with `xtx` (`X'WX`), `xty` (`X'Wz`, one column), `yty`
#'   (`z'Wz`) and `tss`.
#' @examples
#' d <- data.frame(y = rnorm(50), x = rnorm(50), f = gl(5, 10))
#' D <- design_matrix(y ~ x + f, d)
#' m <- design_moments(D, d$y)
#' all.equal(unname(solve(m$xtx, m$xty)[, 1]), unname(coef(lm(y ~ x + f, d))))
#' @keywords internal
#' @export
design_moments <- function(x, z, w = NULL, n = NULL, centered = TRUE) {
    if (is.null(n)) n <- xtx_nrows(x)
    if (inherits(x, "xtx_remote") && is.null(x$scale) && xtx_is_column(z) &&
        (is.null(w) || xtx_is_column(w))) {
        x <- xtx_p_rebase(x, list(z, w))
        y <- xtx_remote_cols(x, z, "response")
        weight <- if (!is.null(w)) xtx_remote_cols(y$D, w, "weights")$sql
        mo <- xtx_sql_moments(y$D, y$sql, weight, sum_y = centered)
        if (!centered) return(mo)
        deviation2 <- (z - mo$sumy / mo$sumw)^2
        mo$tss <- xtx_colsum(if (is.null(w)) deviation2 else w * deviation2)
        return(mo)
    }
    wz <- if (is.null(w)) z else w * z
    sm <- xtx_c_moments(z, w, n)
    list(xtx = if (is.null(w)) xtx_p_crossprod(x) else xtx_p_crossprod(x, xtx_p_scale(x, w)), xty = xtx_p_crossprod(x, wz),
         yty = sm[1L], tss = sm[4L])
}

## The moments of the response z of a fit with fixed effects, weights w
## (NULL for none), n rows: c(sum(w z z), sum(w), sum(w z), sum(w (z -
## zbar)^2)), zbar = sum(w z) / sum(w), without weights w = 1 (and then
## sum(w) is given as 0, zbar = sum(z) / n).
xtx_c_moments <- function(z, w, n) {
    if (xtx_plain(z) && (is.null(w) || xtx_plain(w))) {
        ## (the kernel checks the lengths: one weight per row)
        return(.Call(C_xtx_fe_moments, z, w))
    }
    wz <- if (is.null(w)) z else w * z
    sm <- xtx_c_sums(wz * z, if (is.null(w)) 0 else w, wz)
    zbar <- sm[3L] / (if (is.null(w)) n else sm[2L])
    unname(c(sm, xtx_colsum(if (is.null(w)) (z - zbar)^2 else w * (z - zbar)^2)))
}

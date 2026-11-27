## Weighted least squares on a design (in memory or remote) with lm()'s
## pivoting, and the row patterns that compress it: the least squares core
## of the fits (design_wls(), design_patterns()).

## Groups of rows that share a design pattern, or NULL when the fit should
## not be compressed. `compress = "auto"` follows the option
## `xtx.dedup` (its "auto" search gives up beyond n / 2
## patterns); TRUE searches in any case. Either way, a fit on more than
## half as many patterns as rows is not worth compressing.
xtx_fit_groups <- function(x, compress) {
    if (isFALSE(compress) || nrow(x) < 2L) return(NULL)
    J <- as_jac(x)
    if (!jac_kernel_ok(J)) return(NULL)
    g <- if (isTRUE(compress)) xtx_dedup_groups(J, force = TRUE) else xtx_dedup_groups(J, cost = Inf)
    if (is.null(g)) return(NULL)
    if (length(g$first) == nrow(x)) return(NULL)
    if (!isTRUE(compress) && length(g$first) > 0.5 * nrow(x)) return(NULL)
    g
}

## Per-pattern sufficient statistics of z with weights w (NULL for 1):
## the sum of the weights, the weighted mean, and the weighted sum of
## squares within patterns (computed around the means, which avoids the
## cancellation of sum(w z^2) - sum(w) mean^2).
xtx_fit_suffstats <- function(z, w, map, G) {
    if (is.null(w)) {
        sw <- tabulate(map, nbins = G)
        swz <- drop(rowsum(z, map, reorder = TRUE))
    } else {
        ## doubles: integer weights could overflow in the sums
        w <- as.double(w)
        sw <- drop(rowsum(w, map, reorder = TRUE))
        swz <- drop(rowsum(w * z, map, reorder = TRUE))
    }
    m <- ifelse(sw > 0, swz / sw, 0)
    within <- if (is.null(w)) {
        sum((z - m[map])^2)
    } else {
        ## rows with zero weight do not enter the fit (their squared
        ## deviation can overflow, and 0 * Inf is NaN)
        pos <- w > 0
        sum(w[pos] * (z[pos] - m[map[pos]])^2)
    }
    list(w = as.double(sw), mean = m, within = within)
}

#' Weighted least squares on a design
#'
#' `design_wls()` fits `z` on the columns of a design by weighted least
#' squares, with the pivoting and the rule for aliased columns of
#' [stats::lm.fit()] (same tolerance, same `NA` pattern in the
#' coefficients), without forming the dense matrix: the factor comes from
#' [qr()][design_qr] of the scaled design (`sqrt(w) * x`), and the
#' coefficients are refined from the residuals. The same call works on a
#' design in memory (columns are numeric vectors) and on a remote design
#' (columns are [design_column()]s of it). This is the least squares core of
#' `xtx.fit::lm_xtx()` and of the iterations of `xtx.fit::glm_xtx()`.
#'
#' Names may change before the first release.
#'
#' @param x A design: a `design_matrix`, or a remote design (see
#'   [design_remote]), unscaled.
#' @param z The response: a numeric vector with one value per row of `x`
#'   in memory, a [design_column()] of `x` for a remote design. Finite.
#' @param w Optional weights, in the same form as `z`: finite and
#'   non-negative; rows with zero weight do not enter the fit. `NULL` for
#'   none.
#' @param tol The tolerance for aliased columns, as in [stats::lm.fit()].
#' @param intercept Whether the model has an intercept: the explained sum
#'   of squares `mss` is then around the weighted mean (as in
#'   [summary.lm()]), otherwise around zero.
#' @param sums Also compute the residual and explained sums of squares
#'   (`FALSE`: the solution only, as the steps of a GLM need it).
#' @return A list: `coefficients` (named by the columns, `NA` for aliased
#'   columns), `rank`, `pivot`, `R` (the `p x p` upper factor of `X'WX`,
#'   columns in pivoted order), `tol` and, with `sums = TRUE`, `rss` (the
#'   weighted residual sum of squares), `mss` and `intercept`.
#' @examples
#' d <- data.frame(y = rnorm(50), x = rnorm(50), w = runif(50), f = gl(5, 10))
#' D <- design_matrix(y ~ x + f, d)
#' fit <- design_wls(D, d$y, d$w)
#' all.equal(fit$coefficients, coef(lm(y ~ x + f, d, weights = w)))
#' @export
design_wls <- function(x, z, w = NULL, tol = 1e-07, intercept = TRUE, sums = TRUE) {
    p <- ncol(x)
    cn <- colnames(x)
    if (p == 0L) {
        rss <- xtx_colsum(if (is.null(w)) z^2 else xtx_c_ifelse(w > 0, w * z^2, 0))
        return(list(coefficients = numeric(0),
                    rank = 0L, pivot = integer(0), R = matrix(0, 0, 0),
                    tol = tol, rss = rss, mss = 0, intercept = intercept))
    }
    ## (the primitives of primitives.R: columns are numeric vectors in
    ## memory, remote columns for a remote design)
    sw <- if (is.null(w)) NULL else sqrt(w)
    Y <- if (is.null(sw)) z else sw * z
    q <- xtx_qr_for(xtx_p_scale(x, sw), tol, Y)
    b <- xtx_qr_coef(q, Y)
    beta <- xtx_qr_full(q, b)[, 1L]
    coef <- beta
    coef[setdiff(seq_len(p), q$pivot[seq_len(q$rank)])] <- NA_real_
    names(coef) <- cn
    ## (the weighted least squares steps of a GLM need only the solution)
    if (!sums) return(list(coefficients = coef, rank = q$rank, pivot = q$pivot, R = q$R, tol = tol))
    ## (remotely pending: its view is made only if xtx_c_ls_sums() reads it,
    ## on the frame of z and w)
    f <- xtx_p_mult(x, beta)
    ww <- if (is.null(w)) 1 else w
    ## w v^2, 0 for rows with zero weight (whose squares can overflow)
    wsq <- function(v) if (is.null(w)) v^2 else xtx_c_ifelse(w > 0, w * v^2, 0)
    mo <- attr(b, "moments")
    rss <- NULL
    if (!is.null(mo)) {
        ## from the moments of the solution, when the factorization has them
        ## (a remote design: no pass over its products): rss = z'Wz - b.h,
        ## mss = b'Gb - (sum w f)^2 / sum w, with the intercept column of G
        ## giving sum w f; xtx_c_ls_sums() below when cancellation would
        ## cost accuracy
        zz <- if (is.null(w) && length(q$response_ss) == 1L) {
            c(q$response_ss, q$weight_sum)
        } else {
            xtx_c_sums(wsq(z), ww + 0 * z)
        }
        bh <- sum(mo$b * mo$h)
        rss_m <- zz[1L] - bh
        ic <- match(which(colnames(x) == "(Intercept)"), q$pivot[seq_len(q$rank)])
        if (rss_m > 1e-6 * zz[1L] && (!intercept || !is.na(ic))) {
            rss <- rss_m
            mss <- if (intercept) bh - sum(mo$b * mo$G[, ic])^2 / zz[2L] else bh
            mss <- max(mss, 0)
        }
    }
    if (is.null(rss)) {
        ## (primitives_fused.R: in memory one pass in C)
        s <- xtx_c_ls_sums(z, f, w, intercept)
        rss <- s[1L]
        mss <- s[2L]
    }
    R <- q$R
    list(coefficients = coef, rank = q$rank, pivot = q$pivot, R = R, tol = tol,
         rss = rss, mss = mss, intercept = intercept)
}


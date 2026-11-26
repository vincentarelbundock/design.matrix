## Fused primitives of the fitting code (Kata design.matrix#b89r, FE
## design.matrix#xyk5): each
## computes, from columns, what the fit needs in one call. In memory, on
## plain double vectors (no attributes), one pass in C (src/fused.c, whose
## arithmetic is that of the R code below, to the last bit); otherwise
## (other vectors in memory, remote columns) the operations of
## primitives.R, which build SQL remotely. The fits call these and contain
## no branch on the backend (website/architecture.typ, "One fitting code,
## two backends").
##
## Tests: test_alloc.R (in memory and remote, against the generic code).

## a plain double vector: the kernels read it directly, without coercion
xtx_plain <- function(x) is.double(x) && is.null(attributes(x))

## a xtx_c_codes() grouping in memory whose codes the kernels read directly
## (an integer vector without attributes; the kernels that index by code
## check that the codes are in 1..G)
xtx_plain_codes <- function(g) xtx_is_codes(g) && is.integer(g$codes) && is.null(attributes(g$codes))

## The residual and explained sums of squares of a least squares fit, from
## the response z, the fitted values f and the weights w (NULL for none):
## rss = sum(w r^2) with r = z - f (0 for rows with zero weight, whose
## squares can overflow), mss = sum(w (f - fbar)^2) with fbar the weighted
## mean of f with an intercept, sum(w f^2) without one. c(rss, mss).
xtx_c_ls_sums <- function(z, f, w, intercept) {
    if (xtx_plain(z) && xtx_plain(f) && (is.null(w) || xtx_plain(w))) {
        ## (the kernel checks the lengths: one value per row)
        return(.Call(C_xtx_ls_sums, z, f, w, isTRUE(intercept)))
    }
    ww <- if (is.null(w)) 1 else w
    wsq <- function(v) if (is.null(w)) v^2 else xtx_c_ifelse(w > 0, w * v^2, 0)
    r <- z - f
    rss <- xtx_colsum(wsq(r))
    mss <- if (intercept) {
        fbar <- if (is.null(w)) xtx_c_mean(f) else xtx_colsum(w * f) / xtx_colsum(w)
        xtx_colsum(ww * (f - fbar)^2)
    } else {
        xtx_colsum(ww * f^2)
    }
    c(rss, mss)
}

#' The working weights and responses of IRLS
#'
#' `col_irls_working()` computes, from columns, the working weights and
#' responses of an iteration of iteratively reweighted least squares, with
#' the checks of [stats::glm.fit()], in one call: in memory one pass in C
#' over plain double vectors (the arithmetic of `glm.fit()`, to the last
#' bit), remotely SQL expressions and one query for the checks. With
#' `good = W > 0 & me != 0`, the working weights are `W me^2 / v` and the
#' working responses `(eta - off) + (y - mu) / me` where `good`, 0
#' elsewhere.
#'
#' Names may change before the first release.
#'
#' @param W The prior weights.
#' @param eta The linear predictor (with the offset).
#' @param mu The fitted means.
#' @param y The response.
#' @param off The offset (a column, or a number such as 0).
#' @param v The variance at `mu`.
#' @param me `d mu / d eta` at `eta`.
#' @param table Remotely, store the two columns in a temporary table (for
#'   columns that many queries read, such as those of a fit with fixed
#'   effects); otherwise they are computed once in a projection (see
#'   [col_compute()]).
#' @return A list: the working weights, the working responses (columns),
#'   and the counts of rows where `W > 0` and `v` is missing, `W > 0` and
#'   `v == 0`, `W > 0` and `me` is missing, `good`, and where a working weight
#'   or response is not finite.
#' @examples
#' cf <- column_family(poisson())
#' y <- c(0, 2, 5); eta <- log(c(1, 2, 4)); mu <- cf$linkinv(eta)
#' col_irls_working(1, eta, mu, y, 0, cf$variance(mu), cf$mu_eta(eta, mu))
#' @export
col_irls_working <- function(W, eta, mu, y, off, v, me, table = FALSE) {
    if (xtx_plain(W) && xtx_plain(eta) && xtx_plain(mu) && xtx_plain(y) && xtx_plain(off) &&
        xtx_plain(v) && xtx_plain(me)) {
        ## (the kernel checks the lengths: 1 or length(eta))
        return(.Call(C_xtx_glm_working, W, eta, mu, y, off, v, me))
    }
    good <- (W > 0) & (me != 0)
    wz <- xtx_c_freeze(xtx_c_cbind(xtx_c_ifelse(good, W * me * me / v, 0),
                               xtx_c_ifelse(good, (eta - off) + (y - mu) / me, 0)),
                     table = isTRUE(table))
    w2 <- xtx_c_col(wz, 1L)
    z <- xtx_c_col(wz, 2L)
    s <- xtx_c_sums(xtx_c_count((W > 0) & !xtx_c_present(v)), xtx_c_count((W > 0) & v == 0),
                  xtx_c_count((W > 0) & !xtx_c_present(me)), xtx_c_count(good),
                  xtx_c_count(xtx_c_is_finite(w2) & xtx_c_is_finite(z), negate = TRUE))
    list(w2, z, s)
}

## The deviance residuals of poisson(): 2 w (y log(y / m) - (y - m)) where
## y > 0, 2 m w elsewhere (`positive`: m > 0 by construction, so that the
## logarithm, reached only for y > 0, takes the fast ln of xtx_fast_log()).
xtx_c_poisson_dev <- function(y, m, w, positive = FALSE) {
    if (xtx_plain(y) && xtx_plain(m) && xtx_plain(w)) {
        ## (the kernel checks the lengths: 1 or the longest)
        return(.Call(C_xtx_poisson_dev, y, m, w))
    }
    lg <- if (positive) xtx_fast_log else log
    xtx_c_ifelse(y > 0, 2 * w * (y * lg(y / m) - (y - m)), 2 * m * w)
}

## The weighted binomial log densities of the AIC (as binomial()$aic): m
## trials, y the proportion of successes, w the prior weights, mu the
## means.
xtx_c_binomial_ll <- function(m, y, w, mu) {
    if (xtx_plain(m) && xtx_plain(y) && xtx_plain(w) && xtx_plain(mu)) {
        ## (the kernel checks the lengths: 1 or the longest)
        return(.Call(C_xtx_binomial_ll, m, y, w, mu))
    }
    kk <- xtx_c_round(m * y)
    mm <- xtx_c_round(m)
    xtx_c_ifelse(m > 0, w / m, 0) *
        (xtx_c_lgamma(mm + 1) - xtx_c_lgamma(kk + 1) - xtx_c_lgamma(mm - kk + 1) +
         xtx_c_ifelse(kk > 0, kk * log(mu), 0) + xtx_c_ifelse(mm - kk > 0, (mm - kk) * log(1 - mu), 0))
}

## ---- fixed effects (fe.R) ----

## Weighted regression moments. Center the response in a separate pass to
## retain accuracy when its mean is large relative to its variation.
lm_moments <- function(x, z, w, n = NULL, centered = TRUE) {
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

## The column f plus the values of the fixed effects at the rows' levels:
## f + fe_1[codes_1] + ... + fe_K[codes_K] (codes: per FE, a xtx_c_codes()
## grouping, or, with `levels`, key columns matched to levels, as
## xtx_c_lookup() takes them; NA for a level out of range). With f an n x k
## column and each fe[[k]] a G_k x k matrix, k columns.
xtx_c_fe_add <- function(f, fe, codes, levels = NULL) {
    if (is.null(levels) && length(fe) && xtx_plain(f) && all(vapply(fe, xtx_plain, NA)) &&
        all(vapply(codes, xtx_plain_codes, NA))) {
        ## (the kernel checks the lengths: one code per row)
        return(.Call(C_xtx_fe_fitted_add, f, fe, lapply(codes, `[[`, "codes")))
    }
    for (k in seq_along(fe))
        f <- f + xtx_c_lookup(fe[[k]], codes[[k]], levels[[k]])
    f
}

## The residuals of a fit with fixed effects and their sums, from the
## response z, the fitted values of the design xb, the weights w (NULL for
## none), the FE values fe and codes (per FE, a xtx_c_codes() grouping):
## with e = z - xtx_c_fe_add(xb, fe, codes) and we = w e (e without
## weights), list(we, the sums of we within the levels of each FE (one
## G_k x 1 matrix per FE, as xtx_p_groupsum()), sum(we e)). In memory one
## pass (we comes as an n x 1 matrix, which products take without a copy);
## remotely one query for the sums (xtx_c_groupsums()).
xtx_c_fe_resid <- function(z, xb, w, fe, codes) {
    if (xtx_plain(z) && xtx_plain(xb) && (is.null(w) || xtx_plain(w)) && all(vapply(fe, xtx_plain, NA)) &&
        all(vapply(codes, xtx_plain_codes, NA))) {
        ## (the kernel checks the lengths and the codes)
        return(.Call(C_xtx_fe_resid, z, xb, w, fe, lapply(codes, `[[`, "codes"),
                     vapply(codes, function(g) g$G, 1L)))
    }
    e <- z - xtx_c_fe_add(xb, fe, codes)
    we <- if (is.null(w)) e else w * e
    gs <- xtx_c_groupsums(we, codes, total = we * e)
    list(we, gs$groups, gs$total)
}

## The rows per level of a xtx_c_codes() grouping, as doubles (`n`: those
## counts when known, as xtx_c_dense_codes() gives them in memory).
xtx_p_groupcount <- function(group, n = NULL) {
    if (!is.null(n)) return(as.double(n))
    if (xtx_plain_codes(group)) return(as.double(tabulate(group$codes, nbins = group$G)))
    xtx_p_groupsum(xtx_c_const(group$codes, 1), group)[, 1L]
}

## Columns computed once, together (xtx_c_freeze() of their xtx_c_cbind(): one
## frame remotely, so that they combine with each other): a list of the
## columns. In memory, the columns themselves (no copy).
xtx_c_freeze_all <- function(..., table = FALSE) {
    a <- list(...)
    if (is.null(do.call(xtx_column_frame, a))) return(a)
    m <- xtx_c_freeze(do.call(xtx_c_cbind, a), table = table)
    lapply(seq_along(a), function(j) xtx_c_col(m, j))
}

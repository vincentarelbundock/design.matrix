## Fused primitives of the fitting code (Kata design.matrix#b89r): each
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
dm_plain <- function(x) is.double(x) && is.null(attributes(x))

## a dmc_codes() grouping in memory whose codes the kernels read directly
## (an integer vector without attributes; the kernels that index by code
## check that the codes are in 1..G)
dm_plain_codes <- function(g) dm_is_codes(g) && is.integer(g$codes) && is.null(attributes(g$codes))

## The residual and explained sums of squares of a least squares fit, from
## the response z, the fitted values f and the weights w (NULL for none):
## rss = sum(w r^2) with r = z - f (0 for rows with zero weight, whose
## squares can overflow), mss = sum(w (f - fbar)^2) with fbar the weighted
## mean of f with an intercept, sum(w f^2) without one. c(rss, mss).
dmc_ls_sums <- function(z, f, w, intercept) {
    if (dm_plain(z) && dm_plain(f) && (is.null(w) || dm_plain(w))) {
        ## (the kernel checks the lengths: one value per row)
        return(.Call(C_dm_ls_sums, z, f, w, isTRUE(intercept)))
    }
    ww <- if (is.null(w)) 1 else w
    wsq <- function(v) if (is.null(w)) v^2 else dmc_ifelse(w > 0, w * v^2, 0)
    r <- z - f
    rss <- dm_colsum(wsq(r))
    mss <- if (intercept) {
        fbar <- if (is.null(w)) dmc_mean(f) else dm_colsum(w * f) / dm_colsum(w)
        dm_colsum(ww * (f - fbar)^2)
    } else {
        dm_colsum(ww * f^2)
    }
    c(rss, mss)
}

## The working weights and responses of IRLS, from the prior weights W, the
## linear predictor eta, the means mu, the response y, the offset off, and
## the family's variance v and d mu / d eta me at mu: with good = W > 0 &
## me != 0, w2 = W me^2 / v and z = (eta - off) + (y - mu) / me where good,
## 0 elsewhere. `table`: remotely, the two columns as a temporary table
## (dmc_freeze()). list(w2, z, counts): the rows where W > 0 and v is
## missing, W > 0 and v == 0, W > 0 and me is missing, good, and w2 or z is
## not finite.
dmc_irls_working <- function(W, eta, mu, y, off, v, me, table = FALSE) {
    if (dm_plain(W) && dm_plain(eta) && dm_plain(mu) && dm_plain(y) && dm_plain(off) &&
        dm_plain(v) && dm_plain(me)) {
        ## (the kernel checks the lengths: 1 or length(eta))
        return(.Call(C_dm_glm_working, W, eta, mu, y, off, v, me))
    }
    good <- (W > 0) & (me != 0)
    wz <- dmc_freeze(dmc_cbind(dmc_ifelse(good, W * me * me / v, 0),
                               dmc_ifelse(good, (eta - off) + (y - mu) / me, 0)),
                     table = isTRUE(table))
    w2 <- dmc_col(wz, 1L)
    z <- dmc_col(wz, 2L)
    s <- dmc_sums(dmc_count((W > 0) & !dmc_notnan(v)), dmc_count((W > 0) & v == 0),
                  dmc_count((W > 0) & !dmc_notnan(me)), dmc_count(good),
                  dmc_count(dmc_is_finite(w2) & dmc_is_finite(z), negate = TRUE))
    list(w2, z, s)
}

## The deviance residuals of poisson(): 2 w (y log(y / m) - (y - m)) where
## y > 0, 2 m w elsewhere.
dmc_poisson_dev <- function(y, m, w) {
    if (dm_plain(y) && dm_plain(m) && dm_plain(w)) {
        ## (the kernel checks the lengths: 1 or the longest)
        return(.Call(C_dm_poisson_dev, y, m, w))
    }
    dmc_ifelse(y > 0, 2 * w * (y * log(y / m) - (y - m)), 2 * m * w)
}

## The weighted binomial log densities of the AIC (as binomial()$aic): m
## trials, y the proportion of successes, w the prior weights, mu the
## means.
dmc_binomial_ll <- function(m, y, w, mu) {
    if (dm_plain(m) && dm_plain(y) && dm_plain(w) && dm_plain(mu)) {
        ## (the kernel checks the lengths: 1 or the longest)
        return(.Call(C_dm_binomial_ll, m, y, w, mu))
    }
    kk <- dmc_round(m * y)
    mm <- dmc_round(m)
    dmc_ifelse(m > 0, w / m, 0) *
        (dmc_lgamma(mm + 1) - dmc_lgamma(kk + 1) - dmc_lgamma(mm - kk + 1) +
         dmc_ifelse(kk > 0, kk * log(mu), 0) + dmc_ifelse(mm - kk > 0, (mm - kk) * log(1 - mu), 0))
}

## ---- fixed effects (fe.R) ----

## The moments of the response z of a fit with fixed effects, weights w
## (NULL for none), n rows: c(sum(w z z), sum(w), sum(w z), sum(w (z -
## zbar)^2)), zbar = sum(w z) / sum(w), without weights w = 1 (and then
## sum(w) is given as 0, zbar = sum(z) / n).
dmc_moments <- function(z, w, n) {
    if (dm_plain(z) && (is.null(w) || dm_plain(w))) {
        ## (the kernel checks the lengths: one weight per row)
        return(.Call(C_dm_fe_moments, z, w))
    }
    wz <- if (is.null(w)) z else w * z
    sm <- dmc_sums(wz * z, if (is.null(w)) 0 else w, wz)
    zbar <- sm[3L] / (if (is.null(w)) n else sm[2L])
    unname(c(sm, dm_colsum(if (is.null(w)) (z - zbar)^2 else w * (z - zbar)^2)))
}

## The column f plus the values of the fixed effects at the rows' levels:
## f + fe_1[codes_1] + ... + fe_K[codes_K] (codes: per FE, a dmc_codes()
## grouping, or, with `levels`, key columns matched to levels, as
## dmc_lookup() takes them; NA for a level out of range). With f an n x k
## column and each fe[[k]] a G_k x k matrix, k columns.
dmc_fe_add <- function(f, fe, codes, levels = NULL) {
    if (is.null(levels) && length(fe) && dm_plain(f) && all(vapply(fe, dm_plain, NA)) &&
        all(vapply(codes, dm_plain_codes, NA))) {
        ## (the kernel checks the lengths: one code per row)
        return(.Call(C_dm_fe_fitted_add, f, fe, lapply(codes, `[[`, "codes")))
    }
    for (k in seq_along(fe))
        f <- f + dmc_lookup(fe[[k]], codes[[k]], levels[[k]], along = list(f))
    f
}

## The residuals of a fit with fixed effects and their sums, from the
## response z, the fitted values of the design xb, the weights w (NULL for
## none), the FE values fe and codes (per FE, a dmc_codes() grouping):
## with e = z - dmc_fe_add(xb, fe, codes) and we = w e (e without
## weights), list(we, the sums of we within the levels of each FE (one
## G_k x 1 matrix per FE, as dmp_groupsum()), sum(we e)). In memory one
## pass (we comes as an n x 1 matrix, which products take without a copy);
## remotely one query for the sums (dmc_groupsums()).
dmc_fe_resid <- function(z, xb, w, fe, codes) {
    if (dm_plain(z) && dm_plain(xb) && (is.null(w) || dm_plain(w)) && all(vapply(fe, dm_plain, NA)) &&
        all(vapply(codes, dm_plain_codes, NA))) {
        ## (the kernel checks the lengths and the codes)
        return(.Call(C_dm_fe_resid, z, xb, w, fe, lapply(codes, `[[`, "codes"),
                     vapply(codes, function(g) g$G, 1L)))
    }
    e <- z - dmc_fe_add(xb, fe, codes)
    we <- if (is.null(w)) e else w * e
    gs <- dmc_groupsums(we, codes, total = we * e)
    list(we, gs$groups, gs$total)
}

## The rows per level of a dmc_codes() grouping, as doubles (`n`: those
## counts when known, as dmc_dense_codes() gives them in memory).
dmp_groupcount <- function(group, n = NULL) {
    if (!is.null(n)) return(as.double(n))
    if (dm_plain_codes(group)) return(as.double(tabulate(group$codes, nbins = group$G)))
    dmp_groupsum(dmc_const(group$codes, 1), group)[, 1L]
}

## Columns computed once, together (dmc_freeze() of their dmc_cbind(): one
## frame remotely, so that they combine with each other): a list of the
## columns. In memory, the columns themselves (no copy).
dmc_freeze_all <- function(..., table = FALSE) {
    a <- list(...)
    if (is.null(do.call(dm_rcols_frame, a))) return(a)
    m <- dmc_freeze(do.call(dmc_cbind, a), table = table)
    lapply(seq_along(a), function(j) dmc_col(m, j))
}

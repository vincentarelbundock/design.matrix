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

## Columns computed once, together (dmc_freeze() of their dmc_cbind(): one
## frame remotely, so that they combine with each other): a list of the
## columns. In memory, the columns themselves (no copy).
dmc_freeze_all <- function(..., table = FALSE) {
    a <- list(...)
    if (is.null(do.call(dm_rcols_frame, a))) return(a)
    m <- dmc_freeze(do.call(dmc_cbind, a), table = table)
    lapply(seq_along(a), function(j) dmc_col(m, j))
}

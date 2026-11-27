## GLM families as functions of columns (column_family()): the links,
## variances, deviance residuals and starting values of R's family objects,
## computed on numeric vectors in memory and in SQL remotely, with the same
## code.

xtx_glm_link_names <- c("identity", "log", "logit", "probit", "cloglog", "cauchit", "inverse", "sqrt", "1/mu^2")

## What column_family() supports (and so glm_xtx() of xtx.fit): the
## families of stats and quasi() with five variances, and the standard
## links. Returns the description used by the column functions and the
## checks.
xtx_glm_kind <- function(family) {
    fam <- family$family
    link <- family$link
    ok_fam <- c("gaussian", "binomial", "poisson", "Gamma", "inverse.gaussian",
                "quasibinomial", "quasipoisson", "quasi")
    if (!isTRUE(fam %in% ok_fam) || !isTRUE(link %in% xtx_glm_link_names))
        stop(sprintf("column_family() supports the families %s and the links %s (family `%s`, link `%s`).",
                     paste(ok_fam, collapse = ", "), paste(xtx_glm_link_names, collapse = ", "),
                     format(fam), format(link)), call. = FALSE)
    k <- switch(fam,
                gaussian = list(devk = "gaussian", init = "y", var = "1", validmu = "none"),
                binomial = , quasibinomial = list(devk = "binomial", init = "binomial", var = "mu(1-mu)",
                                                  validmu = "unit"),
                poisson = , quasipoisson = list(devk = "poisson", init = "poisson", var = "mu",
                                                validmu = "positive"),
                Gamma = list(devk = "gamma", init = "y", var = "mu^2", validmu = "positive"),
                inverse.gaussian = list(devk = "invgauss", init = "y", var = "mu^3", validmu = "none"),
                quasi = switch(family$varfun,
                               constant = list(devk = "gaussian", init = "y", var = "1", validmu = "none"),
                               "mu(1-mu)" = list(devk = "binomial", init = "clamp", var = "mu(1-mu)", validmu = "unit_nf"),
                               mu = list(devk = "poisson", init = "y01", var = "mu", validmu = "positive_nf"),
                               "mu^2" = list(devk = "qgamma", init = "y01", var = "mu^2", validmu = "positive_nf"),
                               "mu^3" = list(devk = "invgauss", init = "y01", var = "mu^3", validmu = "positive_nf"),
                               stop(sprintf("column_family() does not support quasi(variance = \"%s\").", family$varfun),
                                    call. = FALSE)))
    k$family <- fam
    k$link <- link
    k$binomial <- fam %in% c("binomial", "quasibinomial")
    k$aic <- if (fam %in% c("quasibinomial", "quasipoisson", "quasi")) "none" else fam
    k$fixed_disp <- fam %in% c("binomial", "poisson")
    k
}

## The family as functions of columns (vectors in memory, SQL remotely),
## with the arithmetic of R's family objects (so that, in memory, the
## values are those of glm.fit(), to the last bit).
xtx_glm_cf <- function(k) {
    eps <- .Machine$double.eps
    ## The clamps of the links use the fast pmax/pmin (xtx_fast_pmax(), plain
    ## greatest()/least() remotely, qpr8): eta is never NA or NaN where they
    ## are used. In the fit, eta = X beta + offset over the kept rows, with X,
    ## offset and the coefficients finite (checked: "NA/NaN/Inf in 'x'",
    ## xf_glm_checks() of xtx.fit, xf_glm_loop() of xtx.fit stops at non-finite coefficients, a
    ## non-finite `start` is refused by xtx_p_mult()); in fitted()
    ## and residuals() it is that of the fit's rows; predict() applies the
    ## inverse link only where eta is present (xf_glm_response() of xtx.fit). On finite
    ## (or infinite) eta, exp(), dnorm(), pnorm() and the like are never NA
    ## or NaN either.
    link <- switch(k$link,
        identity = list(fun = function(m) m, inv = function(e) e, me = function(e) xtx_c_const(e, 1)),
        log = list(fun = function(m) log(m), inv = function(e) xtx_fast_pmax(exp(e), eps),
                   me = function(e) xtx_fast_pmax(exp(e), eps)),
        logit = list(fun = function(m) log(m / (1 - m)), inv = xtx_c_logit_inv, me = xtx_c_logit_mu_eta),
        probit = {
            thr <- -stats::qnorm(eps)
            list(fun = function(m) xtx_c_qnorm(m), inv = function(e) xtx_c_pnorm(xtx_fast_pmin(xtx_fast_pmax(e, -thr), thr)),
                 me = function(e) xtx_fast_pmax(xtx_c_dnorm(e), eps))
        },
        cloglog = list(fun = function(m) log(-log(1 - m)),
                       inv = function(e) xtx_fast_pmax(xtx_fast_pmin(xtx_c_mexpm1m(exp(e)), 1 - eps), eps),
                       me = function(e) {
                           e <- xtx_fast_pmin(e, 700)
                           xtx_fast_pmax(exp(e) * exp(-exp(e)), eps)
                       }),
        cauchit = {
            thr <- -stats::qcauchy(eps)
            list(fun = function(m) xtx_c_qcauchy(m),
                 inv = function(e) xtx_c_pcauchy(xtx_fast_pmin(xtx_fast_pmax(e, -thr), thr)),
                 me = function(e) xtx_fast_pmax(1 / (pi * (1 + e * e)), eps))
        },
        inverse = list(fun = function(m) 1 / m, inv = function(e) 1 / e, me = function(e) -1 / (e * e),
                       valid = function(e) xtx_c_is_finite(e) & e != 0),
        sqrt = list(fun = function(m) sqrt(m), inv = function(e) e * e, me = function(e) 2 * e,
                    valid = function(e) xtx_c_is_finite(e) & e > 0),
        "1/mu^2" = list(fun = function(m) 1 / (m * m), inv = function(e) 1 / sqrt(e),
                        me = function(e) -1 / (2 * e^1.5), valid = function(e) xtx_c_is_finite(e) & e > 0))
    variance <- switch(k$var, "1" = function(m) xtx_c_const(m, 1), "mu(1-mu)" = function(m) m * (1 - m),
                       mu = function(m) m, "mu^2" = function(m) m * m, "mu^3" = function(m) m^3)
    validmu <- switch(k$validmu, none = NULL,
                      unit = function(m) xtx_c_is_finite(m) & m > 0 & m < 1,
                      unit_nf = function(m) xtx_c_present(m) & m > 0 & m < 1,
                      positive = function(m) xtx_c_is_finite(m) & m > 0,
                      positive_nf = function(m) xtx_c_present(m) & m > 0)
    xlogx <- function(a, b) xtx_c_ifelse(a != 0, a * log(a / b), 0)
    dev <- switch(k$devk,
        gaussian = function(y, m, w) w * ((y - m)^2),
        ## (fast ln where mu > 0 by construction, qpr8: the logit, probit,
        ## cloglog and cauchit links clamp mu into [eps, 1 - eps] (logit:
        ## |eta| <= 30), the log link to mu >= eps; y is finite on the
        ## fit's rows: xf_glm_checks() of xtx.fit)
        binomial = if (k$link %in% c("logit", "probit", "cloglog", "cauchit"))
                       function(y, m, w) xtx_c_binomial_dev(y, m, w, positive = TRUE) else xtx_c_binomial_dev,
        poisson = if (k$link == "log") function(y, m, w) xtx_c_poisson_dev(y, m, w, positive = TRUE)
                  else xtx_c_poisson_dev,
        gamma = function(y, m, w) -2 * w * (log(xtx_c_ifelse(y == 0, 1, y / m)) - (y - m) / m),
        qgamma = function(y, m, w) xtx_c_pmax(-2 * w * (log(xtx_c_ifelse(y == 0, 1, y) / m) - (y - m) / m), 0),
        invgauss = function(y, m, w) w * ((y - m)^2) / (y * m^2))
    ## dev(y, mu, w) = w (sat(y) + a part linear in y): the part of the
    ## deviance that does not depend on mu (for compression)
    sat <- switch(k$devk, gaussian = NULL,
                  binomial = function(y) 2 * (xlogx(y, 1) + xlogx(1 - y, 1)),
                  poisson = function(y) 2 * (xlogx(y, 1) - y),
                  gamma = , qgamma = function(y) -2 * log(xtx_c_ifelse(y == 0, 1, y)),
                  invgauss = function(y) 1 / y)
    mustart <- switch(k$init,
        y = function(y, w, n, cb) y,
        binomial = function(y, w, n, cb) if (cb) (n * y + 0.5) / (n + 1) else (w * y + 0.5) / (w + 1),
        poisson = function(y, w, n, cb) y + 0.1,
        y01 = function(y, w, n, cb) y + 0.1 * xtx_c_count(y == 0),
        ## (y is finite: xf_glm_checks() of xtx.fit; fast clamps, see the links above)
        clamp = function(y, w, n, cb) xtx_fast_pmax(0.001, xtx_fast_pmin(0.999, y)))
    ## (never NA: each check starts with is.finite() or !is.na(), FALSE where
    ## the value is missing; remotely flagged as never NULL, so that its
    ## sums need no count of NULL values)
    valid <- function(eta, mu) {
        a <- link$valid
        b <- validmu
        if (is.null(a) && is.null(b)) return(NULL)
        v <- if (is.null(a)) b(mu) else if (is.null(b)) a(eta) else a(eta) & b(mu)
        xtx_c_notnull(v)
    }
    ## fitted values numerically at the boundary, for glm.fit()'s warnings
    ## (binomial: below 10 eps or above 1 - 10 eps; poisson: below 10 eps),
    ## as a count (1 or 0, never NA); NULL for the other families
    beps <- 10 * .Machine$double.eps
    boundary <- switch(k$family,
                       binomial = function(mu) xtx_c_count(mu > 1 - beps | mu < beps),
                       poisson = function(mu) xtx_c_count(mu < beps),
                       NULL)
    list(linkfun = link$fun, linkinv = link$inv, mu.eta = link$me, variance = variance, valid = valid,
         boundary = boundary, dev = dev, sat = sat, mustart = mustart, me_is_mu = k$link == "log")
}

## d mu / d eta at eta, mu = linkinv(eta): for the log link mu itself (the
## same expression, not computed twice).
xtx_glm_mu_eta <- function(cf, eta, mu) {
    if (isTRUE(cf$me_is_mu)) return(mu)
    cf$mu.eta(eta)
}


#' GLM families as functions of columns
#'
#' `column_family()` describes a family object of [stats::family] by
#' functions of columns: numeric vectors in memory, [design_column()]s of a
#' remote design (SQL expressions, evaluated in the database). In memory the
#' values are those of R's family object, to the last bit; remotely the same
#' arithmetic in SQL. This is how `xtx.fit::glm_xtx()` runs one IRLS code on both
#' backends.
#'
#' The families of [stats::family] are supported: `binomial`, `poisson`,
#' `gaussian`, `Gamma`, `inverse.gaussian`, `quasibinomial`, `quasipoisson`
#' and `quasi` (with variance `"constant"`, `"mu(1-mu)"`, `"mu"`, `"mu^2"` or
#' `"mu^3"`), with the links `identity`, `log`, `logit`, `probit`,
#' `cloglog`, `cauchit`, `inverse`, `sqrt` and `1/mu^2`; other families are an
#' error.
#'
#' Names may change before the first release.
#'
#' @param family A family object, a family function or its name (as in
#'   [stats::glm()]).
#' @return An object of class `xtx_column_family`, a list with
#' \describe{
#'   \item{`family`, `link`}{The names of the family and the link.}
#'   \item{`binomial`}{Whether the family is binomial or quasibinomial.}
#'   \item{`aic`}{The family whose AIC applies (`"none"` for quasi families).}
#'   \item{`linkfun(mu)`, `linkinv(eta)`, `mu.eta(eta)`, `variance(mu)`}{As in
#'     the family object, on columns (the inverse links clamp as R does).}
#'   \item{`mu_eta(eta, mu)`}{`mu.eta(eta)` at `mu = linkinv(eta)` (for the
#'     log link, `mu` itself).}
#'   \item{`valid(eta, mu)`}{A logical column, never `NA`: where `eta` and
#'     `mu` are valid (`validmu()` and `valideta()`), or `NULL` when every
#'     value is.}
#'   \item{`boundary(mu)`}{For `binomial` and `poisson`, a column that counts
#'     (`TRUE`/1) the fitted values numerically at the boundary (`glm.fit()`'s
#'     "fitted probabilities numerically 0 or 1" and "fitted rates numerically
#'     0"), never `NA`; `NULL` for the other families.}
#'   \item{`dev(y, mu, w)`}{The deviance residuals, as `dev.resids()`.}
#'   \item{`sat(y)`}{The part of the deviance residuals that does not depend
#'     on `mu`, per unit of weight (`NULL` for the Gaussian): `dev(y, mu, w) =
#'     w (sat(y) + a part linear in y)`, which makes the deviance of rows that
#'     share `mu` that of their weighted mean response plus a constant.}
#'   \item{`mustart(y, w, n, cbind)`}{The starting values of `initialize`
#'     (`n` the trials of a binomial response, `cbind` whether it had two
#'     columns).}
#'   \item{`loglik(y, mu, w, m, disp)`}{The weighted log densities of the AIC
#'     (binomial: `m` the trials; Gamma: `disp` the dispersion; other
#'     families: zero).}
#' }
#' @examples
#' cf <- column_family(binomial("probit"))
#' eta <- c(-1, 0, 2)
#' all.equal(cf$linkinv(eta), binomial("probit")$linkinv(eta))
#' all.equal(cf$dev(c(0, 1, 1), cf$linkinv(eta), 1),
#'           binomial("probit")$dev.resids(c(0, 1, 1), cf$linkinv(eta), 1))
#' @export
column_family <- function(family) {
    if (is.character(family)) family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) family <- family()
    if (!is.list(family) || is.null(family$family)) stop("'family' not recognized", call. = FALSE)
    k <- xtx_glm_kind(family)
    cf <- xtx_glm_cf(k)
    cf$mu_eta <- function(eta, mu) xtx_glm_mu_eta(cf, eta, mu)
    cf$loglik <- function(y, mu, w, m = NULL, disp = NULL) {
        switch(k$aic,
            binomial = xtx_c_binomial_ll(m, y, w, mu),
            poisson = w * xtx_c_ifelse(y == 0, -mu, y * log(mu) - mu - xtx_c_lgamma(y + 1)),
            Gamma = {
                a <- 1 / disp
                w * (-lgamma(a) - a * log(mu * disp) + (a - 1) * log(y) - y / (mu * disp))
            },
            0 * y)
    }
    structure(c(k, cf), class = "xtx_column_family")
}

#' @rdname column_family
#' @param x A `xtx_column_family`.
#' @param ... Not used.
#' @export
print.xtx_column_family <- function(x, ...) {
    cat(sprintf("<column_family: %s, link %s>\n", x$family, x$link))
    invisible(x)
}

## The column flagged as never NULL (remotely; the caller knows it is not).
xtx_c_notnull <- function(x) {
    if (!xtx_is_column(x)) return(x)
    u <- unclass(x)
    u$notnull <- rep(TRUE, length(u$sql))
    structure(u, class = "xtx_column")
}

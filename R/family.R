## GLM families as functions of columns (column_family()): the links,
## variances, deviance residuals and starting values of R's family objects,
## computed on numeric vectors in memory and in SQL remotely, with the same
## code.

xtx_glm_link_names <- c("identity", "log", "logit", "probit", "cloglog", "cauchit", "inverse", "sqrt", "1/mu^2")

## What design_glm() supports: the families of stats and quasi() with five
## variances, and the standard links. Returns the description used by the
## column functions and the checks.
xtx_glm_kind <- function(family) {
    fam <- family$family
    link <- family$link
    ok_fam <- c("gaussian", "binomial", "poisson", "Gamma", "inverse.gaussian",
                "quasibinomial", "quasipoisson", "quasi")
    if (!isTRUE(fam %in% ok_fam) || !isTRUE(link %in% xtx_glm_link_names))
        stop(sprintf("design_glm() supports the families %s and the links %s (family `%s`, link `%s`).",
                     paste(ok_fam, collapse = ", "), paste(xtx_glm_link_names, collapse = ", "),
                     format(fam), format(link)), call. = FALSE)
    k <- switch(fam,
                gaussian = list(dev = "gaussian", init = "y", var = "1", validmu = "none"),
                binomial = , quasibinomial = list(dev = "binomial", init = "binomial", var = "mu(1-mu)",
                                                  validmu = "unit"),
                poisson = , quasipoisson = list(dev = "poisson", init = "poisson", var = "mu",
                                                validmu = "positive"),
                Gamma = list(dev = "gamma", init = "y", var = "mu^2", validmu = "positive"),
                inverse.gaussian = list(dev = "invgauss", init = "y", var = "mu^3", validmu = "none"),
                quasi = switch(family$varfun,
                               constant = list(dev = "gaussian", init = "y", var = "1", validmu = "none"),
                               "mu(1-mu)" = list(dev = "binomial", init = "clamp", var = "mu(1-mu)", validmu = "unit_nf"),
                               mu = list(dev = "poisson", init = "y01", var = "mu", validmu = "positive_nf"),
                               "mu^2" = list(dev = "qgamma", init = "y01", var = "mu^2", validmu = "positive_nf"),
                               "mu^3" = list(dev = "invgauss", init = "y01", var = "mu^3", validmu = "positive_nf"),
                               stop(sprintf("design_glm() does not support quasi(variance = \"%s\").", family$varfun),
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
    ## xtx_glm_checks(), xtx_glm_loop() stops at non-finite coefficients, a
    ## non-finite `start` is refused by xtx_p_mult()); in fitted()
    ## and residuals() it is that of the fit's rows; predict() applies the
    ## inverse link only where eta is present (xtx_glm_response()). On finite
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
    dev <- switch(k$dev,
        gaussian = function(y, m, w) w * ((y - m)^2),
        ## (fast ln where mu > 0 by construction, qpr8: the logit, probit,
        ## cloglog and cauchit links clamp mu into [eps, 1 - eps] (logit:
        ## |eta| <= 30), the log link to mu >= eps; y is finite on the
        ## fit's rows: xtx_glm_checks())
        binomial = if (k$link %in% c("logit", "probit", "cloglog", "cauchit"))
                       function(y, m, w) xtx_c_binomial_dev(y, m, w, positive = TRUE) else xtx_c_binomial_dev,
        poisson = if (k$link == "log") function(y, m, w) xtx_c_poisson_dev(y, m, w, positive = TRUE)
                  else xtx_c_poisson_dev,
        gamma = function(y, m, w) -2 * w * (log(xtx_c_ifelse(y == 0, 1, y / m)) - (y - m) / m),
        qgamma = function(y, m, w) xtx_c_pmax(-2 * w * (log(xtx_c_ifelse(y == 0, 1, y) / m) - (y - m) / m), 0),
        invgauss = function(y, m, w) w * ((y - m)^2) / (y * m^2))
    ## dev(y, mu, w) = w (sat(y) + a part linear in y): the part of the
    ## deviance that does not depend on mu (for compression)
    sat <- switch(k$dev, gaussian = NULL,
                  binomial = function(y) 2 * (xlogx(y, 1) + xlogx(1 - y, 1)),
                  poisson = function(y) 2 * (xlogx(y, 1) - y),
                  gamma = , qgamma = function(y) -2 * log(xtx_c_ifelse(y == 0, 1, y)),
                  invgauss = function(y) 1 / y)
    mustart <- switch(k$init,
        y = function(y, w, n, cb) y,
        binomial = function(y, w, n, cb) if (cb) (n * y + 0.5) / (n + 1) else (w * y + 0.5) / (w + 1),
        poisson = function(y, w, n, cb) y + 0.1,
        y01 = function(y, w, n, cb) y + 0.1 * xtx_c_count(y == 0),
        ## (y is finite: xtx_glm_checks(); fast clamps, see the links above)
        clamp = function(y, w, n, cb) xtx_fast_pmax(0.001, xtx_fast_pmin(0.999, y)))
    valid <- function(eta, mu) {
        a <- link$valid
        b <- validmu
        if (is.null(a) && is.null(b)) return(NULL)
        if (is.null(a)) return(b(mu))
        if (is.null(b)) return(a(eta))
        a(eta) & b(mu)
    }
    list(linkfun = link$fun, linkinv = link$inv, mu.eta = link$me, variance = variance, valid = valid,
         dev = dev, sat = sat, mustart = mustart, me_is_mu = k$link == "log")
}

## d mu / d eta at eta, mu = linkinv(eta): for the log link mu itself (the
## same expression, not computed twice).
xtx_glm_mu_eta <- function(cf, eta, mu) {
    if (isTRUE(cf$me_is_mu)) return(mu)
    cf$mu.eta(eta)
}


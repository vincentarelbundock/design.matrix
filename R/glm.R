## Generalized linear models fitted by IRLS on a design, in memory or remote,
## with one implementation on the primitives of primitives.R: the response,
## weights, linear predictor, fitted means and working weights are columns
## (numeric vectors in memory, SQL expressions remotely), the weighted least
## squares steps are qr() of the scaled design, and the sums over rows are
## dm_colsum().
##
## The iterations are those of glm.fit(), step for step: same starting
## values (family$initialize, or start), same weighted least squares problem
## at each iteration (lm()'s pivoting with glm.fit()'s tolerance), same step
## halving, the same convergence rule and the same warnings.
##
## Compression. Rows that share a design pattern share eta and mu (a fit
## with an offset is not compressed). For every supported family the
## deviance has the form dev(y, mu, w) = w (sat(y) + y h(mu) + k(mu)), so
## at fixed eta the working weights and responses of a pattern depend on
## its rows only through the sum of the prior weights W and the weighted
## mean response ybar: the iterations on the patterns, with weights W and
## response ybar, are those on the rows, and the deviance of the rows is
## that of the patterns plus C = sum_i w_i sat(y_i) - sum_g W_g sat(ybar_g)
## (for the Gaussian, the within-pattern sum of squares). The first
## iteration starts from per-row values, so it runs on the rows; the AIC,
## dispersion and null deviance are computed on the rows at the end.

#' Generalized linear models without the model matrix
#'
#' `design_glm()` fits a generalized linear model by iteratively reweighted
#' least squares (IRLS), like [stats::glm()], without forming the dense
#' model matrix. It takes a formula and a data frame, a [design_matrix()]
#' and a response `y` (like [stats::glm.fit()]), or remote data (see
#' *Remote data* below).
#'
#' **Method.** The iterations are those of [glm.fit()]: the same starting
#' values (from `family$initialize`, or `start`), the
#' same weighted least squares problem at each iteration, solved with
#' [qr()][design_qr] on the scaled design (columns pivoted, and aliased
#' columns detected, as in `glm.fit()`, with its tolerance
#' `min(1e-7, epsilon / 1000)`), the same step halving when the deviance is
#' not finite or the fitted values are invalid, the same convergence rule
#' `|dev - dev_old| / (|dev| + 0.1) < epsilon`, and the same warnings
#' (non-convergence, boundary, fitted probabilities numerically 0 or 1).
#' Results agree with those of `glm()` up to floating point rounding: the
#' number of iterations is the same, and the coefficients, covariance
#' matrix, deviances and AIC agree to many more digits than `epsilon`.
#'
#' **Compression.** Rows with the same design pattern and the same offset
#' have the same fitted value. For every family that `design_glm()`
#' supports, the deviance has the form
#' `w * (a(y) + y * h(mu) + k(mu))`, so at fixed coefficients a pattern
#' contributes to each iteration only through the sum of its prior weights
#' and its weighted mean response, and to the deviance through a constant
#' that is computed once. The fit then runs on the patterns and gives the
#' results of the fit on every row, with the same iterations (the first
#' iteration, which starts from per-row values, runs on the rows). The AIC,
#' the dispersion and the null deviance are computed from the rows. A fit
#' with an offset is not compressed. `compress` works as in [design_lm()].
#'
#' **Families.** The families of [stats::family]: `binomial`, `poisson`,
#' `gaussian`, `Gamma`, `inverse.gaussian`, `quasibinomial`, `quasipoisson`
#' and `quasi` (with variance `"constant"`, `"mu(1-mu)"`, `"mu"`, `"mu^2"`
#' or `"mu^3"`), with the links `identity`, `log`, `logit`, `probit`,
#' `cloglog`, `cauchit`, `inverse`, `sqrt` and `1/mu^2`. These are the
#' families and links that also have a SQL translation, so that every fit
#' works in memory and on remote data alike; other family objects are an
#' error.
#'
#' **Remote data.** `data` can be a DuckDB connection (with `table` or
#' `sql`) or Parquet paths, and the first argument a remote design (see
#' [design_remote]); `weights` and `offset` are then column names, `DBI::SQL()`
#' or R expressions over the columns.
#' The fit is the same code as in memory: its columns (response, weights,
#' linear predictor, working weights) are SQL expressions, and its products
#' with the design run in the database. See the section
#' *GLMs* of [design_remote].
#'
#' **Fixed effects.** A formula `y ~ x | fe1 + fe2` absorbs fixed effects
#' (FE), as [fixest::feglm()] and [fixest::fepois()]: each FE is a variable
#' or an interaction `a^b`. The iterations are those of `glm.fit()` (same
#' starting values, step halving, deviance convergence rule and warnings),
#' with the slopes and the FE values as parameters; the weighted least
#' squares step of each iteration is the estimator of linear models with
#' FE of [design_lm()] (sums over rows, the largest FE eliminated exactly,
#' the others solved densely or by conjugate gradients), so the dummies are
#' never formed. The same code fits data in memory and remote data.
#' * Rows: rows with a missing value in the model, FE or cluster variables,
#'   or a zero weight, are dropped; then, iteratively, rows whose FE level
#'   the FE fit perfectly, as fixest's `fixef.rm`: `"perfect_fit"`
#'   (default) drops singletons and, when the inverse link maps onto
#'   (0, Inf) or (0, 1) (fixest's rule: log, logit, probit, cloglog links),
#'   levels whose response is all 0 or all 1 (their FE would be infinite);
#'   `"singleton"`, `"infinite_coef"` and `"none"` drop only the first, only
#'   the second, or nothing. `fit$nsingletons` and `fit$nperfect` count the
#'   rows dropped.
#' * Coefficients: the intercept is absorbed; a column collinear with the
#'   FE or earlier columns is `NA` (fixest's rule, `collin.tol = 1e-9`).
#' * Standard errors: `vcov = "iid"` (default: the inverse information of
#'   the last weighted least squares step, times the dispersion) or
#'   `vcov = ~cl` (clustered by one variable), with fixest's default small
#'   sample corrections, `(n - 1) / (n - K)` with `K` the number of slopes
#'   and FE parameters, and `G / (G - 1)` for clusters (FE nested in the
#'   clusters not counted in `K`).
#' * Methods: those of a `design_glm` fit (class `c("design_glm_fe",
#'   "design_glm")`); `logLik()`, `AIC()` and `BIC()` count the FE
#'   parameters; `fixef()` (the generic of fixest or nlme, when loaded)
#'   returns the FE values, also in `fit$fixef`, normalized as for
#'   [design_lm()]; `predict()` gives the link or the response (`NA` for
#'   levels not in the fit); remotely, per-row results are tables.
#'   `fit$fe_iter` holds the iterations of conjugate gradients of each
#'   step (0 for a direct solve).
#' * Not supported with FE: `start`, `compress = TRUE` (compression is
#'   turned off), `na.action` other than `na.omit` and `na.exclude`,
#'   clustering by more than one variable, a design matrix as first
#'   argument.
#'
#' **What is kept.** As in [design_lm()]: the coefficients, the factor `R`
#' of the last weighted least squares problem (which gives the covariance
#' matrix, as in [summary.glm()]), the deviances, AIC, degrees of freedom,
#' the number of iterations, and, with `model = TRUE`, the model frame and
#' design, from which [residuals()], [fitted()] and [predict()] are
#' recomputed.
#'
#' **Methods.** [coef()], [vcov()] (`NA` for aliased coefficients, scaled
#' by the dispersion), [summary()] (coefficient table, dispersion, null and
#' residual deviance, AIC, iterations), [predict()] (`type = "link"` or
#' `"response"`), [residuals()] (types `"deviance"`,
#' `"pearson"`, `"working"` and `"response"`), [fitted()], [deviance()],
#' [df.residual()], [nobs()], [logLik()] (hence [AIC()] and [BIC()]),
#' [family()], [formula()], [terms()], [model.frame()] and `print()`.
#'
#' **Differences with `glm()`.** The result is not a `glm` object: it has
#' no `qr`, `effects`, `residuals`, `fitted.values` or `weights`
#' components, so functions that expect those (such as `anova.glm()` or
#' `plot.lm()`) do not apply. Every option works in memory and on remote
#' data alike; so `subset`, `etastart`, `mustart`, `method`, `x`,
#' `singular.ok = FALSE`, `na.action` other than `na.omit`, `na.exclude` and
#' `na.fail`, and `se.fit` and `type = "terms"` in `predict()` are not
#' supported. The null deviance of a
#' fit with an offset and an intercept is refitted, as `glm()` does (also
#' for a design matrix, where `glm.fit()` would not). Extra arguments do not
#' go to [glm.control()]: give `control = list(maxit = 50)`.
#'
#' @param formula A model formula, as in [stats::glm()], a [design_matrix()]
#'   (then give the response as `y`) or a remote design.
#' @param family A family object, function or name, as in [stats::glm()].
#' @param data,weights,na.action,start,offset,contrasts
#'   As in [stats::glm()] (`na.action`: `na.omit`, `na.exclude` or
#'   `na.fail`; with fixed effects `na.omit` or `na.exclude`). With a design matrix, `weights` and `offset` are numeric
#'   vectors with one value per row, and `data`, `na.action` and `contrasts`
#'   do not apply. For remote data, `weights` and `offset` are column
#'   names, `DBI::SQL()` expressions or R expressions over the columns
#'   (evaluated in batches when they have no SQL translation), as in
#'   [design_lm()].
#' @param control A list of arguments of [stats::glm.control()].
#' @param compress `"auto"`, `TRUE` or `FALSE`: fit on the design patterns
#'   (see Details).
#' @param model Keep the model frame and its design matrix, for
#'   [residuals()], [fitted()] and [predict()] without `newdata`.
#' @param y With a design matrix as first argument, the response: a numeric
#'   vector, a factor, or a two-column matrix of successes and failures (as
#'   in [glm.fit()]). With a remote design, the name of the response column,
#'   or the names of two columns (successes and failures).
#' @param vcov With fixed effects: `"iid"` or a one-sided formula naming a
#'   cluster variable (`~firm`). Without fixed effects only `"iid"`.
#' @param fixef.rm With fixed effects: `"perfect_fit"` (default),
#'   `"singleton"`, `"infinite_coef"` or `"none"`; see *Fixed effects*.
#' @param ... For remote data: `table`, `sql`, `engine` and `chunk_size`; see
#'   [design_remote]. Otherwise not used.
#' @return An object of class `"design_glm"`: a list with, among others,
#'   `coefficients`, `deviance`, `null.deviance`, `aic`, `iter`, `converged`,
#'   `boundary`, `df.residual`, `df.null`, `nobs`, `rank`, `pivot`, `R`,
#'   `family`, `compressed`, `patterns`, `rows`, `spec`, `call` and `terms`.
#' @examples
#' set.seed(1)
#' n <- 2000
#' d <- data.frame(x = rnorm(n), f = factor(sample(letters[1:4], n, TRUE)),
#'                 g = factor(sample(1:3, n, TRUE)))
#' d$y <- rbinom(n, 1, plogis(d$x + as.integer(d$f) / 4))
#' fit <- design_glm(y ~ x * f, binomial, d)
#' all.equal(coef(fit), coef(glm(y ~ x * f, binomial, d)))
#' all.equal(vcov(fit), vcov(glm(y ~ x * f, binomial, d)))
#' summary(fit)
#'
#' ## factors only: fitted on 12 patterns instead of 2000 rows
#' d$k <- rpois(n, exp(0.5 + as.integer(d$g) / 3))
#' fit <- design_glm(k ~ f * g, poisson, d)
#' fit$patterns
#' all.equal(deviance(fit), deviance(glm(k ~ f * g, poisson, d)))
#'
#' predict(fit, newdata = d[1:3, ], type = "response")
#'
#' ## fixed effects
#' panel <- data.frame(id = rep(1:100, each = 10), t = rep(1:10, 100), x = rnorm(1000))
#' panel$y <- rpois(1000, exp(0.5 * panel$x + panel$id / 100))
#' fe <- design_glm(y ~ x | id + t, poisson, panel, vcov = ~id)
#' summary(fe)
#' @export
design_glm <- function(formula, family = stats::gaussian, data, weights, na.action,
                       start = NULL, offset, control = list(),
                       contrasts = NULL, compress = "auto", model = TRUE, y = NULL,
                       vcov = "iid", fixef.rm = "perfect_fit", ...) {
    family <- dmglm_family(family, parent.frame())
    dmglm_kind(family)
    control <- do.call(stats::glm.control, as.list(control))
    cl <- match.call()
    fe <- dm_fe_has_bar(formula)
    if (inherits(formula, c("design_matrix", "dm_remote"))) {
        dmfit_check_compress(compress)
        checkmate::assert_flag(model)
        if (!identical(vcov, "iid"))
            stop("`vcov` other than \"iid\" requires fixed effects (y ~ x | fe).", call. = FALSE)
        if (!missing(data) || !missing(na.action) || !is.null(contrasts))
            stop("With a design matrix, give the response as `y`; `data`, ",
                 "`na.action` and `contrasts` do not apply.", call. = FALSE)
        if (...length()) stop("Unknown arguments: ", paste(names(list(...)), collapse = ", "), call. = FALSE)
        ## (remotely, `weights` and `offset` are columns or expressions)
        fit <- if (inherits(formula, "dm_remote")) {
            dm_remote_glm_design(formula, y, if (missing(weights)) NULL else substitute(weights),
                                 if (missing(offset)) NULL else substitute(offset), family, start, control,
                                 compress, parent.frame())
        } else {
            dmglm_design(formula, y, if (missing(weights)) NULL else weights,
                         if (missing(offset)) NULL else offset, start, family, control, compress)
        }
        fit$call <- cl
        return(fit)
    }
    if (!is.null(y)) stop("`y` applies only when the first argument is a design matrix.", call. = FALSE)
    if (!missing(na.action)) {
        if (fe) dm_fe_na_kind(na.action) else dm_na_kind(na.action, c("omit", "exclude", "fail"))
    }
    ## fixed effects (fe_glm.R)
    if (fe) {
        return(dm_fe_glm(cl, parent.frame(), formula, family, if (missing(data)) NULL else data, start, control,
                         contrasts, compress, model, vcov, fixef.rm, dm_dots(...)))
    }
    ## (the checks of design_lm(), in the same order)
    if (!identical(vcov, "iid"))
        stop("`vcov` other than \"iid\" requires fixed effects (y ~ x | fe).", call. = FALSE)
    dmfit_check_compress(compress)
    checkmate::assert_flag(model)
    ## the design and columns of a data frame or remote data (remote_fit.R)
    P <- design_problem(cl, parent.frame(), formula, if (missing(data)) NULL else data, contrasts, dm_dots(...),
                        family = dmglm_kind(family))
    fit <- dmglm_fit(P$x, P$cols, family, start, control, attr(P$terms, "intercept") > 0L, compress,
                     wmsg = "negative weights not allowed")
    fit$call <- cl
    fit$terms <- P$terms
    fit$formula <- P$formula
    fit$contrasts <- P$contrasts
    fit$xlevels <- P$xlevels
    fit$na.action <- P$na.action
    fit$offset_call <- cl$offset
    fit$env <- P$env
    ## (the data of a remote fit stay in the database: always kept)
    if (model) fit$model <- P$frame else if (!is.null(P$frame)) fit$data <- NULL
    fit
}

## The family, as glm() resolves it.
dmglm_family <- function(family, env) {
    if (is.character(family)) family <- get(family, mode = "function", envir = env)
    if (is.function(family)) family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized", call. = FALSE)
    }
    family
}

## ---- families on columns ----------------------------------------------------

dmglm_link_names <- c("identity", "log", "logit", "probit", "cloglog", "cauchit", "inverse", "sqrt", "1/mu^2")

## What design_glm() supports: the families of stats and quasi() with five
## variances, and the standard links. Returns the description used by the
## column functions and the checks.
dmglm_kind <- function(family) {
    fam <- family$family
    link <- family$link
    ok_fam <- c("gaussian", "binomial", "poisson", "Gamma", "inverse.gaussian",
                "quasibinomial", "quasipoisson", "quasi")
    if (!isTRUE(fam %in% ok_fam) || !isTRUE(link %in% dmglm_link_names))
        stop(sprintf("design_glm() supports the families %s and the links %s (family `%s`, link `%s`).",
                     paste(ok_fam, collapse = ", "), paste(dmglm_link_names, collapse = ", "),
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
dmglm_cf <- function(k) {
    eps <- .Machine$double.eps
    link <- switch(k$link,
        identity = list(fun = function(m) m, inv = function(e) e, me = function(e) dmc_const(e, 1)),
        log = list(fun = function(m) log(m), inv = function(e) dmc_pmax(exp(e), eps),
                   me = function(e) dmc_pmax(exp(e), eps)),
        logit = list(fun = function(m) log(m / (1 - m)), inv = dmc_logit_inv, me = dmc_logit_mu_eta),
        probit = {
            thr <- -stats::qnorm(eps)
            list(fun = function(m) dmc_qnorm(m), inv = function(e) dmc_pnorm(dmc_pmin(dmc_pmax(e, -thr), thr)),
                 me = function(e) dmc_pmax(dmc_dnorm(e), eps))
        },
        cloglog = list(fun = function(m) log(-log(1 - m)),
                       inv = function(e) dmc_pmax(dmc_pmin(dmc_mexpm1m(exp(e)), 1 - eps), eps),
                       me = function(e) {
                           e <- dmc_pmin(e, 700)
                           dmc_pmax(exp(e) * exp(-exp(e)), eps)
                       }),
        cauchit = {
            thr <- -stats::qcauchy(eps)
            list(fun = function(m) dmc_qcauchy(m),
                 inv = function(e) dmc_pcauchy(dmc_pmin(dmc_pmax(e, -thr), thr)),
                 me = function(e) dmc_pmax(1 / (pi * (1 + e * e)), eps))
        },
        inverse = list(fun = function(m) 1 / m, inv = function(e) 1 / e, me = function(e) -1 / (e * e),
                       valid = function(e) dmc_is_finite(e) & e != 0),
        sqrt = list(fun = function(m) sqrt(m), inv = function(e) e * e, me = function(e) 2 * e,
                    valid = function(e) dmc_is_finite(e) & e > 0),
        "1/mu^2" = list(fun = function(m) 1 / (m * m), inv = function(e) 1 / sqrt(e),
                        me = function(e) -1 / (2 * e^1.5), valid = function(e) dmc_is_finite(e) & e > 0))
    variance <- switch(k$var, "1" = function(m) dmc_const(m, 1), "mu(1-mu)" = function(m) m * (1 - m),
                       mu = function(m) m, "mu^2" = function(m) m * m, "mu^3" = function(m) m^3)
    validmu <- switch(k$validmu, none = NULL,
                      unit = function(m) dmc_is_finite(m) & m > 0 & m < 1,
                      unit_nf = function(m) dmc_notnan(m) & m > 0 & m < 1,
                      positive = function(m) dmc_is_finite(m) & m > 0,
                      positive_nf = function(m) dmc_notnan(m) & m > 0)
    xlogx <- function(a, b) dmc_ifelse(a != 0, a * log(a / b), 0)
    dev <- switch(k$dev,
        gaussian = function(y, m, w) w * ((y - m)^2),
        binomial = dmc_binomial_dev,
        poisson = function(y, m, w) {
            ## (in memory one pass in C)
            if (dmglm_mem_working(list(y, m, w), max(length(y), length(m), length(w))))
                return(.Call(C_dm_poisson_dev, y, m, w))
            dmc_ifelse(y > 0, 2 * w * (y * log(y / m) - (y - m)), 2 * m * w)
        },
        gamma = function(y, m, w) -2 * w * (log(dmc_ifelse(y == 0, 1, y / m)) - (y - m) / m),
        qgamma = function(y, m, w) dmc_pmax(-2 * w * (log(dmc_ifelse(y == 0, 1, y) / m) - (y - m) / m), 0),
        invgauss = function(y, m, w) w * ((y - m)^2) / (y * m^2))
    ## dev(y, mu, w) = w (sat(y) + a part linear in y): the part of the
    ## deviance that does not depend on mu (for compression)
    sat <- switch(k$dev, gaussian = NULL,
                  binomial = function(y) 2 * (xlogx(y, 1) + xlogx(1 - y, 1)),
                  poisson = function(y) 2 * (xlogx(y, 1) - y),
                  gamma = , qgamma = function(y) -2 * log(dmc_ifelse(y == 0, 1, y)),
                  invgauss = function(y) 1 / y)
    mustart <- switch(k$init,
        y = function(y, w, n, cb) y,
        binomial = function(y, w, n, cb) if (cb) (n * y + 0.5) / (n + 1) else (w * y + 0.5) / (w + 1),
        poisson = function(y, w, n, cb) y + 0.1,
        y01 = function(y, w, n, cb) y + 0.1 * dmc_count(y == 0),
        clamp = function(y, w, n, cb) dmc_pmax(0.001, dmc_pmin(0.999, y)))
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
## same expression; in memory, no second pass).
dmglm_mu_eta <- function(cf, eta, mu) {
    if (isTRUE(cf$me_is_mu) && !dm_is_rcols(eta) && !dm_is_rcols(mu)) return(mu)
    cf$mu.eta(eta)
}

## ---- columns of the response ------------------------------------------------

## The response, prior weights, trials and offset as glm()'s initialize sees
## them (columns: vectors in memory, SQL remotely): y in [0, 1] and weights
## times trials for a two-column binomial response (a successes, b
## failures), y = 0 where the weight is 0.
dmglm_cols <- function(y, a, b, w, off, k, cb, weighted = TRUE) {
    y0 <- y
    if (cb) {
        if (!k$binomial)
            stop("A two-column response (successes and failures) needs a binomial family.", call. = FALSE)
        nt <- a + b
        y <- dmc_ifelse(nt == 0, 0, a / nt)
        w <- w * nt
    } else {
        nt <- dmc_const(y, 1)
        if (k$binomial) y <- dmc_ifelse(w == 0, 0, y)
    }
    ## (trials weight the rows of a two-column response; y0 the response as
    ## given)
    list(y = y, w = w, n = nt, off = off, cbind = cb, a = a, b = b, weighted = weighted || cb, y0 = y0)
}

## The response of a model frame in memory (a vector, a factor or character
## vector, or a two-column matrix) as columns.
dmglm_mem_cols <- function(Y, w, off, k) {
    if (length(dim(Y)) == 1L) dim(Y) <- NULL
    if (is.matrix(Y) && ncol(Y) == 1L) Y <- Y[, 1L]
    n <- NROW(Y)
    w <- if (is.null(w)) rep.int(1, n) else as.double(w)
    off <- if (!is.null(off)) as.double(off)
    if (is.matrix(Y)) {
        if (ncol(Y) != 2L) stop("A two-column response needs two columns (successes and failures).", call. = FALSE)
        return(dmglm_cols(NULL, as.double(Y[, 1L]), as.double(Y[, 2L]), w, off, k, TRUE))
    }
    ## a character response is a factor, as remotely
    if (is.character(Y)) Y <- factor(Y)
    if (is.factor(Y)) {
        if (!k$binomial) stop("A factor response needs a binomial family.", call. = FALSE)
        Y <- Y != levels(Y)[1L]
    }
    if (!is.numeric(Y) && !is.logical(Y)) stop("The response must be numeric, logical or a factor.", call. = FALSE)
    y <- as.double(Y)
    attributes(y) <- NULL
    dmglm_cols(y, NULL, NULL, w, off, k, FALSE)
}

## design_glm() on a design matrix x and a response y.
dmglm_design <- function(x, y, weights, offset, start, family, control, compress) {
    checkmate::assert_class(x, "design_matrix")
    if (is.null(y)) stop("Give the response as `y`.", call. = FALSE)
    n <- nrow(x)
    if (NROW(y) != n) stop(sprintf("The response must have %d values.", n), call. = FALSE)
    if (!is.null(weights)) checkmate::assert_numeric(weights, len = n)
    if (!is.null(offset)) checkmate::assert_numeric(offset, len = n)
    tt <- x$training_terms
    intercept <- if (!is.null(tt)) attr(tt, "intercept") > 0L else "(Intercept)" %in% colnames(x)
    cols <- dmglm_mem_cols(y, weights, offset, dmglm_kind(family))
    dmglm_check_design(cols)
    dmglm_fit(x, cols, family, start, control, intercept, compress, wmsg = "negative weights not allowed")
}

## The checks of a fit on a design (no rows dropped): no missing response
## or weight, a finite offset.
dmglm_check_design <- function(cols) {
    y <- if (cols$cbind) cols$a + cols$b else cols$y0
    s <- dmc_sums(dmc_count(dmc_notnan(y), negate = TRUE),
                  if (cols$weighted && !cols$cbind) dmc_count(dmc_notnan(cols$w), negate = TRUE) else 0,
                  if (!is.null(cols$off)) dmc_count(dmc_is_finite(cols$off), negate = TRUE) else 0)
    bad <- is.na(s) | s > 0
    if (bad[1L]) stop("NA/NaN/Inf in 'y'", call. = FALSE)
    if (bad[2L]) stop("missing or negative weights not allowed", call. = FALSE)
    if (bad[3L]) stop("NA/NaN/Inf in 'offset'", call. = FALSE)
}

## ---- checks ------------------------------------------------------------------

## The checks of glm() and of the family's initialize, and the sums every
## fit needs, in one pass over the rows.
dmglm_checks <- function(cols, k, start, wmsg) {
    y <- cols$y
    w <- cols$w
    one <- 0 * y + 1
    zero <- 0 * y
    sums <- list(n = one, nz = dmc_count(w != 0), sw = w, swy = w * y,
                 maxn = if (k$binomial) dmc_count(cols$n > 1) else zero,
                 slw = if (k$aic == "gaussian") dmc_ifelse(w > 0, log(w), 0) else zero,
                 nzero = dmc_count(w == 0), neg = dmc_count(w < 0), infw = dmc_count(abs(w) == Inf),
                 ly = if (k$aic == "inverse.gaussian") dmc_ifelse(y > 0, w * log(y), 0) else zero,
                 nonint = if (k$aic == "poisson") dmc_count(y != dmc_round(y)) else zero)
    chk <- list()
    add <- function(cond, msg, type = "error") chk[[length(chk) + 1L]] <<- list(cond = cond, msg = msg, type = type)
    fam <- k$family
    if (k$binomial && !cols$cbind) add(y < 0 | y > 1, "y values must be 0 <= y <= 1")
    if (fam == "binomial" && !cols$cbind)
        add(abs(w * y - dmc_round(w * y)) > 0.001, "non-integer #successes in a binomial glm!", "warning")
    if (fam == "binomial" && cols$cbind)
        add(abs(cols$a - dmc_round(cols$a)) > 0.001 | abs(cols$b - dmc_round(cols$b)) > 0.001,
            "non-integer counts in a binomial glm!", "warning")
    if (fam == "poisson") add(y < 0, "negative values not allowed for the 'Poisson' family")
    if (fam == "quasipoisson") add(y < 0, "negative values not allowed for the 'quasiPoisson' family")
    if (fam == "Gamma") add(y <= 0, "non-positive values not allowed for the 'Gamma' family")
    if (fam == "inverse.gaussian") add(y <= 0, "positive values only are allowed for the 'inverse.gaussian' family")
    if (fam == "gaussian" && is.null(start) && k$link %in% c("inverse", "log"))
        add(if (k$link == "inverse") y == 0 else y <= 0, "cannot find valid starting values: please specify some")
    all_cols <- c(sums, lapply(chk, function(x) dmc_count(x$cond)))
    s <- do.call(dmc_sums, unname(all_cols))
    ck <- stats::setNames(as.list(s[seq_along(sums)]), names(sums))
    if (ck$n == 0) stop("0 (non-NA) cases", call. = FALSE)
    if (ck$neg > 0) stop(wmsg, call. = FALSE)
    if (ck$infw > 0) stop("NA/NaN/Inf in 'weights'", call. = FALSE)
    for (j in seq_along(chk)) {
        v <- s[length(sums) + j]
        if (is.na(v) || v > 0) {
            if (chk[[j]]$type == "error") stop(chk[[j]]$msg, call. = FALSE) else warning(chk[[j]]$msg, call. = FALSE)
        }
    }
    ck
}

## ---- the fit ------------------------------------------------------------------
##
## The iterations of glm.fit() (dmglm_loop()) run on "units": the rows of the
## design (columns: vectors in memory, SQL remotely), or the patterns of
## rows that share their design inputs (dmp_patterns(), vectors), with the
## response ybar (the weighted mean of y), the prior weights W (their sum)
## and the deviance plus C, the part of the deviance of the rows that does
## not depend on the coefficients. The first iteration starts from per-row
## values, so it always runs on the rows.

## `solver`: NULL (the weighted least squares steps are qr() of the scaled
## design), or the steps of another problem with the same loop: list(lin,
## wls, np), lin(theta, along) the linear predictor (without offset) at
## parameters theta (a vector, or a matrix with one column per parameter
## vector: several columns), wls(w2, z) the weighted least squares step
## (as dmglm_wls_design()), np the number of parameters (fixed effects:
## dm_fe_solver()). The fit is then not compressed and has no `start`.
dmglm_fit <- function(x, cols, family, start, control, intercept, compress, wmsg, solver = NULL) {
    k <- dmglm_kind(family)
    cf <- dmglm_cf(k)
    ck <- dmglm_checks(cols, k, start, wmsg)
    tol <- min(1e-07, control$epsilon / 1000)
    p <- ncol(x)
    if (!is.null(start) && length(start) != p)
        stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
                      p, paste(deparse(colnames(x)), collapse = ", ")), call. = FALSE)
    y <- cols$y
    w <- cols$w
    off <- cols$off %||% 0
    lin <- if (is.null(solver)) function(b, along = list()) dmp_mult(x, b, along = along) else solver$lin
    ## (with another solver, whose steps read the working columns many
    ## times, they are a table remotely: `materialize`)
    rows <- list(lin = lin, wls = if (is.null(solver)) dmglm_wls_design(x, tol) else solver$wls,
                 y = y, W = w, off = off, C = 0, materialize = !is.null(solver))
    np <- if (is.null(solver)) p else solver$np
    ## starting values, on the rows
    eta0 <- dmc_freeze(if (!is.null(start)) dmp_mult(x, start) + off
                         else cf$linkfun(cf$mustart(y, w, cols$n, cols$cbind)))
    mu0 <- dmc_freeze(cf$linkinv(eta0))
    v0 <- cf$valid(eta0, mu0)
    s0 <- dmc_sums(cf$dev(y, mu0, w), if (is.null(v0)) 0 * y else dmc_count(v0, negate = TRUE))
    if (isTRUE(is.na(s0[2L]) || s0[2L] > 0))
        stop("cannot find valid starting values: please specify some", call. = FALSE)
    devold <- s0[1L]

    pt <- if (p > 0L && is.null(solver) && !isFALSE(compress))
        dmp_patterns(x, y, if (isTRUE(cols$weighted)) w, compress, keys = cols$off)
    units <- rows
    if (!is.null(pt)) {
        W <- pt$w %||% rep(1, length(pt$z))
        xp <- pt$x
        C <- if (is.null(cf$sat)) pt$within else dm_colsum(w * cf$sat(y)) - sum(W * cf$sat(pt$z))
        if (is.finite(C)) {
            units <- list(lin = function(b) dmp_mult(xp, b), wls = dmglm_wls_design(xp, tol), y = pt$z, W = W,
                          off = if (is.null(pt$keys)) 0 else pt$keys[, 1L], C = C)
        } else {
            pt <- NULL
        }
    }
    B <- dmglm_units(units, cf, k, first = if (!is.null(pt)) c(rows, list(eta = eta0, mu = mu0)),
                     eta = eta0, mu = mu0)
    if (np == 0L) {
        res <- list(coef = numeric(0), fit = list(rank = 0L, pivot = integer(0), R = matrix(0, 0, 0)),
                    iter = 0L, conv = TRUE, boundary = TRUE, dev = devold)
    } else {
        res <- dmglm_loop(B, k, control, ck$n, devold, start)
    }
    ## the summaries, on the rows
    beta <- res$coef
    beta[is.na(beta)] <- 0
    old <- B$state$wls_beta
    fin <- dmglm_final(if (np == 0L) NULL else lin, cols, cf, k, ck, beta, old, eta0, res$dev, intercept)
    if (!is.null(cols$off) && intercept && np > 0L)
        fin$nulldev <- dmglm_null_refit(cols, cf, k, control, ck$n, fin$mu)
    nobs <- ck$nz
    if (nobs < .Machine$integer.max) nobs <- as.integer(nobs)
    rank <- res$fit$rank
    out <- list(coefficients = if (is.null(solver)) stats::setNames(res$coef, colnames(x)) else res$coef,
                rank = rank, pivot = res$fit$pivot,
                R = res$fit$R, tol = tol, deviance = res$dev, null.deviance = fin$nulldev,
                aic = dmglm_aic(k, ck, res$dev, fin$ll) + 2 * rank, iter = res$iter, converged = res$conv,
                boundary = res$boundary, df.residual = nobs - rank, df.null = nobs - as.integer(intercept),
                nobs = nobs, rows = ck$n, pearson = fin$pearson, zero_weights = fin$zero, family = family,
                intercept = intercept, compressed = !is.null(pt), patterns = if (is.null(pt)) ck$n else length(units$W),
                control = control, spec = dmfit_spec(x),
                method = if (!is.null(pt)) "patterns" else switch(dm_engine(x), sql = "SQL iterations", NULL),
                data = list(design = x, y = y, weights = w, offset = cols$off))
    ## (the last weighted least squares step of another solver)
    if (!is.null(solver)) out$solver_fit <- res$fit
    class(out) <- "design_glm"
    out
}

## Weighted least squares of z on a design (rows or patterns) with squared
## IRLS weights w2 (0 for rows left out): the least squares core of
## design_lm() (dmfit_ls(): qr() of the scaled design, lm()'s pivoting), at
## glm.fit()'s tolerance, as glm.fit()'s Cdqrls.
dmglm_wls_design <- function(x, tol) {
    force(x)
    function(w2, z) {
        f <- dmfit_ls(x, z, w2, tol, FALSE, sums = FALSE)
        b <- f$coefficients
        b[is.na(b)] <- 0
        list(coefficients = unname(b), rank = f$rank, pivot = f$pivot, R = unname(f$R))
    }
}

## The same for an intercept alone (the null model): w2 and z sums.
dmglm_wls_mean <- function(w2, z) {
    s <- dmc_sums(w2, w2 * z)
    if (!all(is.finite(s))) return(list(coefficients = NaN, rank = 0L, pivot = 1L, R = matrix(0, 1L, 1L)))
    list(coefficients = if (s[1L] > 0) s[2L] / s[1L] else 0, rank = if (s[1L] > 0) 1L else 0L,
         pivot = 1L, R = matrix(sqrt(s[1L]), 1L, 1L))
}

## The working weights and responses of IRLS at (eta, mu), with glm.fit()'s
## checks; NULL when no observation is informative.
dmglm_working <- function(u, eta, mu, cf) {
    W <- u$W
    v <- cf$variance(mu)
    me <- dmglm_mu_eta(cf, eta, mu)
    if (dmglm_mem_working(list(W, eta, mu, u$y, u$off, v, me), length(eta))) {
        ## in memory: the same values and counts in one pass (C)
        r <- .Call(C_dm_glm_working, W, eta, mu, u$y, u$off, v, me)
        return(dmglm_working_checks(r[[3L]], r[[1L]], r[[2L]]))
    }
    good <- (W > 0) & (me != 0)
    wz <- dmc_freeze(dmc_cbind(dmc_ifelse(good, W * me * me / v, 0),
                                 dmc_ifelse(good, (eta - u$off) + (u$y - mu) / me, 0)),
                     table = isTRUE(u$materialize))
    w2 <- dmc_col(wz, 1L)
    z <- dmc_col(wz, 2L)
    s <- dmc_sums(dmc_count((W > 0) & !dmc_notnan(v)), dmc_count((W > 0) & v == 0),
                               dmc_count((W > 0) & !dmc_notnan(me)), dmc_count(good),
                               dmc_count(dmc_is_finite(w2) & dmc_is_finite(z), negate = TRUE))
    dmglm_working_checks(s, w2, z)
}

## Can the working columns be computed by the kernel (in memory, plain
## double vectors of length 1 or n)?
dmglm_mem_working <- function(cols, n) {
    all(vapply(cols, function(v) is.double(v) && is.null(dim(v)) && length(v) %in% c(1L, n), NA)) && n > 0L
}

## The checks of dmglm_working() from its counts s.
dmglm_working_checks <- function(s, w2, z) {
    if (isTRUE(is.na(s[1L]) || s[1L] > 0)) stop("NAs in V(mu)", call. = FALSE)
    if (isTRUE(s[2L] > 0)) stop("0s in V(mu)", call. = FALSE)
    if (isTRUE(is.na(s[3L]) || s[3L] > 0)) stop("NAs in d(mu)/d(eta)", call. = FALSE)
    if (!isTRUE(s[4L] > 0)) return(NULL)
    list(w2 = w2, z = z, finite = isTRUE(s[5L] == 0))
}

## The backend of the loop on units (see dmglm_loop()): `first`, the rows
## and their starting values, for the first iteration of a fit on
## patterns; otherwise the units start at (eta, mu).
dmglm_units <- function(u, cf, k, first = NULL, eta = NULL, mu = NULL) {
    st <- new.env(parent = emptyenv())
    st$first <- first
    st$eta <- eta
    st$mu <- mu
    eps <- 10 * .Machine$double.eps
    list(
        state = st,
        wls = function() {
            if (!is.null(st$first)) {
                f <- st$first
                st$first <- NULL
                wk <- dmglm_working(f, f$eta, f$mu, cf)
                if (is.null(wk)) return(NULL)
                st$wls_beta <- NULL
                if (!wk$finite) return(list(coefficients = NaN, rank = 0L))
                return(f$wls(wk$w2, wk$z))
            }
            wk <- dmglm_working(u, st$eta, st$mu, cf)
            if (is.null(wk)) return(NULL)
            st$wls_beta <- st$beta
            if (!wk$finite) return(list(coefficients = NaN, rank = 0L))
            u$wls(wk$w2, wk$z)
        },
        eval = function(beta) {
            eta <- u$lin(beta)
            ## (in memory, no pass to add a zero offset)
            if (dm_is_rcols(eta) || !identical(u$off, 0)) eta <- eta + u$off
            mu <- dmc_freeze(cf$linkinv(eta))
            v <- cf$valid(eta, mu)
            ext <- if (k$family == "binomial") dmc_count(mu > 1 - eps | mu < eps)
                   else if (k$family == "poisson") dmc_count(mu < eps) else 0 * mu
            s <- dmc_sums(cf$dev(u$y, mu, u$W), if (is.null(v)) 0 * mu else dmc_count(v, negate = TRUE), ext)
            st$eta <- eta
            st$mu <- mu
            st$beta <- beta
            st$extreme <- s[3L]
            list(dev = s[1L] + u$C, valid = isTRUE(s[2L] == 0))
        },
        extreme = function() isTRUE(st$extreme > 0))
}

## The IRLS loop of glm.fit(), on a backend: wls() (the weighted least
## squares fit at the current eta, or NULL when no observation is
## informative), eval(beta) (move to beta: list(dev, valid)) and extreme()
## (fitted values numerically at the boundary, for the warnings).
dmglm_loop <- function(B, k, control, nobs, devold, coefold = NULL) {
    conv <- boundary <- FALSE
    coef <- fit <- NULL
    dev <- devold
    maxit <- control$maxit
    iter <- 0L
    for (iter in seq_len(maxit)) {
        f <- B$wls()
        if (is.null(f)) {
            warning(gettextf("no observations informative at iteration %d", iter), call. = FALSE)
            break
        }
        fit <- f
        if (any(!is.finite(f$coefficients))) {
            warning(gettextf("non-finite coefficients at iteration %d", iter), call. = FALSE)
            break
        }
        if (nobs < f$rank)
            stop(sprintf(ngettext(nobs, "X matrix has rank %d, but only %d observation",
                                  "X matrix has rank %d, but only %d observations"), f$rank, nobs),
                 call. = FALSE)
        start <- f$coefficients
        ev <- B$eval(start)
        dev <- ev$dev
        if (control$trace) cat("Deviance = ", dev, " Iterations - ", iter, "\n", sep = "")
        boundary <- FALSE
        if (!is.finite(dev)) {
            if (is.null(coefold))
                stop("no valid set of coefficients has been found: please supply starting values", call. = FALSE)
            warning("step size truncated due to divergence", call. = FALSE)
            ii <- 1
            while (!is.finite(dev)) {
                if (ii > maxit) stop("inner loop 1; cannot correct step size", call. = FALSE)
                ii <- ii + 1
                start <- (start + coefold) / 2
                ev <- B$eval(start)
                dev <- ev$dev
            }
            boundary <- TRUE
            if (control$trace) cat("Step halved: new deviance = ", dev, "\n", sep = "")
        }
        if (!ev$valid) {
            if (is.null(coefold))
                stop("no valid set of coefficients has been found: please supply starting values", call. = FALSE)
            warning("step size truncated: out of bounds", call. = FALSE)
            ii <- 1
            while (!ev$valid) {
                if (ii > maxit) stop("inner loop 2; cannot correct step size", call. = FALSE)
                ii <- ii + 1
                start <- (start + coefold) / 2
                ev <- B$eval(start)
            }
            boundary <- TRUE
            dev <- ev$dev
            if (control$trace) cat("Step halved: new deviance = ", dev, "\n", sep = "")
        }
        if (abs(dev - devold) / (0.1 + abs(dev)) < control$epsilon) {
            conv <- TRUE
            coef <- start
            break
        }
        devold <- dev
        coef <- coefold <- start
    }
    if (is.null(coef)) stop("The fit failed: no valid coefficients were found.", call. = FALSE)
    if (!conv) warning("glm.fit: algorithm did not converge", call. = FALSE)
    if (boundary) warning("glm.fit: algorithm stopped at boundary value", call. = FALSE)
    if (k$family == "binomial" && B$extreme())
        warning("glm.fit: fitted probabilities numerically 0 or 1 occurred", call. = FALSE)
    if (k$family == "poisson" && B$extreme())
        warning("glm.fit: fitted rates numerically 0 occurred", call. = FALSE)
    p <- length(coef)
    if (!is.null(fit$aliased)) {
        coef[fit$aliased] <- NA
    } else if (fit$rank < p) {
        coef[fit$pivot[seq.int(fit$rank + 1L, p)]] <- NA
    }
    list(coef = coef, fit = fit, iter = iter, conv = conv, boundary = boundary, dev = dev)
}

## The summaries of a fit, on the rows: Pearson sum for the dispersion
## (working weights of the last weighted least squares step, at `old`, or
## at the starting values eta0 when that step was the first; working
## residuals at beta), the null deviance (without the refit), the sum of
## weighted log densities of the AIC, and the fitted means. `lin`: the
## linear predictor at parameters (see dmglm_fit()), NULL without any.
dmglm_final <- function(lin, cols, cf, k, ck, beta, old, eta0, dev, intercept) {
    y <- cols$y
    w <- cols$w
    off <- cols$off %||% 0
    if (is.null(lin)) {
        eta <- 0 * y + off
        eta_old <- eta
    } else if (is.null(old)) {
        eta <- lin(beta, along = list(eta0)) + off
        eta_old <- eta0
    } else {
        P <- lin(cbind(beta, old))
        eta <- dmc_col(P, 1L) + off
        eta_old <- dmc_col(P, 2L) + off
    }
    ## (frozen together: columns of one frame; in memory, no need)
    if (dm_is_rcols(eta)) {
        mm <- dmc_freeze(dmc_cbind(cf$linkinv(eta), cf$linkinv(eta_old)))
        mu <- dmc_col(mm, 1L)
        mu0 <- dmc_col(mm, 2L)
    } else {
        mu <- cf$linkinv(eta)
        mu0 <- cf$linkinv(eta_old)
    }
    me <- dmglm_mu_eta(cf, eta, mu)
    me0 <- dmglm_mu_eta(cf, eta_old, mu0)
    good0 <- (w > 0) & (me0 != 0)
    w0 <- dmc_ifelse(good0, w * me0 * me0 / cf$variance(mu0), 0)
    r <- (y - mu) / me
    wtd <- if (intercept) ck$swy / ck$sw else cf$linkinv(0 * y + off)
    ll <- switch(k$aic,
        binomial = {
            m <- if (ck$maxn > 0) cols$n else w
            if (dmglm_mem_working(list(m, y, w, mu), max(length(m), length(y), length(w), length(mu)))) {
                ## (in memory one pass in C)
                .Call(C_dm_binomial_ll, m, y, w, mu)
            } else {
                kk <- dmc_round(m * y)
                mm <- dmc_round(m)
                dmc_ifelse(m > 0, w / m, 0) *
                    (dmc_lgamma(mm + 1) - dmc_lgamma(kk + 1) - dmc_lgamma(mm - kk + 1) +
                     dmc_ifelse(kk > 0, kk * log(mu), 0) + dmc_ifelse(mm - kk > 0, (mm - kk) * log(1 - mu), 0))
            }
        },
        poisson = w * dmc_ifelse(y == 0, -mu, y * log(mu) - mu - dmc_lgamma(y + 1)),
        Gamma = {
            disp <- dev / ck$sw
            a <- 1 / disp
            w * (-lgamma(a) - a * log(mu * disp) + (a - 1) * log(y) - y / (mu * disp))
        },
        0 * y)
    s <- dmc_sums(dmc_ifelse(good0, w0 * r * r, 0), dmc_count(good0, negate = TRUE), cf$dev(y, wtd, w), ll)
    list(pearson = s[1L], zero = isTRUE(s[2L] > 0), nulldev = s[3L], ll = s[4L], mu = mu)
}

## The AIC of glm.fit() (without 2 * rank) from sums over rows: `ll` is the
## sum of weighted log densities for binomial, poisson and Gamma.
dmglm_aic <- function(k, ck, dev, ll) {
    switch(k$aic,
           gaussian = {
               nobs <- ck$n
               if (ck$nzero > 0) Inf else nobs * (log(dev / nobs * 2 * pi) + 1) + 2 - ck$slw
           },
           binomial = -2 * ll,
           poisson = if (ck$nonint > 0) Inf else -2 * ll,
           Gamma = -2 * ll + 2,
           inverse.gaussian = ck$sw * (1 + log(dev / ck$sw * 2 * pi)) + 3 * ck$ly + 2,
           NA_real_)
}

## The null deviance of a fit with an offset and an intercept, refitted as
## glm() does: intercept only, starting from the fitted means of the model.
dmglm_null_refit <- function(cols, cf, k, control, nobs, mu) {
    u <- list(lin = function(b) b[1L], wls = dmglm_wls_mean, y = cols$y, W = cols$w, off = cols$off, C = 0)
    eta <- cf$linkfun(mu)
    mu1 <- cf$linkinv(eta)
    v <- cf$valid(eta, mu1)
    s <- dmc_sums(cf$dev(cols$y, mu1, cols$w), if (is.null(v)) 0 * mu1 else dmc_count(v, negate = TRUE))
    if (isTRUE(is.na(s[2L]) || s[2L] > 0)) stop("cannot find valid starting values: please specify some", call. = FALSE)
    B <- dmglm_units(u, cf, k, eta = eta, mu = mu1)
    res <- dmglm_loop(B, k, control, nobs, s[1L])
    if (!res$conv)
        warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?", call. = FALSE)
    res$dev
}

## The dispersion, as summary.glm() computes it (with its warning).
dmglm_dispersion <- function(object, warn = TRUE) {
    fam <- object$family
    if (!is.null(fam$dispersion) && !is.na(fam$dispersion)) return(list(value = fam$dispersion, est = FALSE))
    if (fam$family %in% c("poisson", "binomial")) return(list(value = 1, est = FALSE))
    df.r <- object$df.residual
    if (df.r > 0) {
        if (warn && isTRUE(object$zero_weights))
            warning("observations with zero weight not used for calculating dispersion", call. = FALSE)
        return(list(value = object$pearson / df.r, est = TRUE))
    }
    list(value = NaN, est = TRUE)
}

## ---- row-wise results ------------------------------------------------------

## The kept rows of a fit: design and columns.
dmglm_rows <- function(object, what) {
    if (is.null(object$data)) {
        stop(sprintf("%s needs the fit's data: refit with `model = TRUE`, or use `predict(fit, newdata)`.", what),
             call. = FALSE)
    }
    object$data
}

## Row-wise results: a vector in memory (the last column, named, padded by
## pad(na.action)), a remote table remotely (dmp_output()).
dmglm_output <- function(x, cols, na_act, pad = stats::napredict) {
    dmp_output(x, cols, na = na_act, memory = function(d) {
        v <- d[[ncol(d)]]
        names(v) <- dimnames(x)[[1L]]
        pad(na_act, v)
    })
}

## A linear predictor (with its NAs kept as NA by the inverse link).
dmglm_response <- function(cf, eta) dmc_ifelse(dmc_notnan(eta), cf$linkinv(eta), eta)

#' @export
fitted.design_glm <- function(object, ...) {
    d <- dmglm_rows(object, "fitted()")
    cf <- dmglm_cf(dmglm_kind(object$family))
    eta <- dmglm_eta(object, d$design, d$offset, d$codes)
    dmglm_output(d$design, list(.fitted = cf$linkinv(eta)), object$na.action)
}

#' @export
residuals.design_glm <- function(object, type = c("deviance", "pearson", "working", "response"), ...) {
    type <- match.arg(type)
    d <- dmglm_rows(object, "residuals()")
    cf <- dmglm_cf(dmglm_kind(object$family))
    eta <- dmglm_eta(object, d$design, d$offset, d$codes)
    mu <- cf$linkinv(eta)
    y <- d$y
    w <- d$weights
    res <- switch(type,
                  deviance = if (object$df.residual > 0) {
                      dd <- sqrt(dmc_pmax(cf$dev(y, mu, w), 0))
                      dmc_ifelse(y > mu, dd, -dd)
                  } else 0 * mu,
                  pearson = (y - mu) * sqrt(w) / sqrt(cf$variance(mu)),
                  working = (y - mu) / cf$mu.eta(eta),
                  response = y - mu)
    dmglm_output(d$design, list(.fitted = mu, .resid = res), object$na.action, pad = stats::naresid)
}

#' @param object A `"design_glm"` object.
#' @param newdata A data frame of new observations, or remote data (a DuckDB
#'   connection with `table` or `sql`, or Parquet paths). Omitted,
#'   predictions are those of the fit's rows (which need `model = TRUE`).
#' @param type `"link"` (the linear predictor) or `"response"` (the fitted
#'   means).
#' @rdname design_glm
#' @export
predict.design_glm <- function(object, newdata = NULL, type = c("link", "response"),
                               na.action = stats::na.pass, ...) {
    type <- match.arg(type)
    cf <- dmglm_cf(dmglm_kind(object$family))
    if (is.null(newdata)) {
        if (...length()) stop("Unknown arguments: ", paste(names(list(...)), collapse = ", "), call. = FALSE)
        d <- dmglm_rows(object, "predict() without `newdata`")
        x <- d$design
        off <- d$offset
        keys <- d$codes
        levels <- NULL
        na_act <- object$na.action
    } else {
        ## (with fixed effects, every row of newdata is kept)
        if (!is.null(object$fixef) && !identical(na.action, stats::na.pass))
            stop("With fixed effects, predict() keeps every row of `newdata` (na.action = na.pass).", call. = FALSE)
        nd <- dmfit_newdata(object, newdata, na.action, list(...))
        x <- nd$D
        off <- nd$off
        keys <- if (!is.null(object$fixef)) dm_fe_new_keys(object, x, newdata)
        levels <- object$fe_levels
        na_act <- NULL
        if (object$rank < length(object$coefficients))
            warning("prediction from rank-deficient fit; consider predict(., rankdeficient=\"NA\")", call. = FALSE)
    }
    eta <- dmglm_eta(object, x, off, keys, levels)
    dmglm_output(x, list(.fitted = if (type == "response") dmglm_response(cf, eta) else eta), na_act)
}

#' @export
print.design_glm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
    dmfit_print_rows(x)
    dmfit_print_coef(x, digits, none = "\nNo coefficients\n\n")
    cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ", x$df.residual, "Residual\n")
    if (nzchar(mess <- stats::naprint(x$na.action))) cat("  (", mess, ")\n", sep = "")
    cat("Null Deviance:\t   ", format(signif(x$null.deviance, digits)),
        "\nResidual Deviance:", format(signif(x$deviance, digits)), "\tAIC:",
        format(signif(x$aic, digits)))
    cat("\n")
    invisible(x)
}

#' @export
coef.design_glm <- function(object, ...) object$coefficients

#' @export
vcov.design_glm <- function(object, complete = TRUE, ...) {
    V <- dmfit_vcov(object, dmglm_dispersion(object)$value)
    if (!complete) {
        keep <- !is.na(object$coefficients)
        V <- V[keep, keep, drop = FALSE]
    }
    V
}

#' @export
nobs.design_glm <- function(object, ...) object$nobs

#' @export
deviance.design_glm <- function(object, ...) object$deviance

#' @export
df.residual.design_glm <- function(object, ...) object$df.residual

#' @export
family.design_glm <- function(object, ...) object$family

#' @export
formula.design_glm <- function(x, ...) {
    if (!is.null(x$formula)) return(x$formula)
    if (is.null(x$terms)) return(NULL)
    stats::formula(x$terms)
}

#' @export
terms.design_glm <- function(x, ...) {
    if (is.null(x$terms)) stop("The fit has no terms: it was not fitted from a formula.", call. = FALSE)
    x$terms
}

#' @export
model.frame.design_glm <- function(formula, ...) {
    if (is.null(formula$model))
        stop("The fit has no model frame: refit `design_glm()` with `model = TRUE`.", call. = FALSE)
    formula$model
}

#' @export
logLik.design_glm <- function(object, ...) {
    fam <- object$family$family
    dispersion <- object$family$dispersion
    p <- object$rank
    if (!is.null(dispersion)) {
        if (is.na(dispersion)) p <- p + 1L
    } else if (fam %in% c("gaussian", "Gamma", "inverse.gaussian")) {
        p <- p + 1L
    }
    structure(p - object$aic / 2, nobs = object$rows, df = p, class = "logLik")
}

#' @export
summary.design_glm <- function(object, dispersion = NULL, ...) {
    if (is.null(dispersion)) {
        d <- dmglm_dispersion(object)
        dispersion <- d$value
        est.disp <- d$est
    } else {
        est.disp <- FALSE
    }
    df.r <- object$df.residual
    aliased <- is.na(object$coefficients)
    p <- object$rank
    if (p > 0) {
        keep <- object$pivot[seq_len(p)]
        coef.p <- object$coefficients[keep]
        covmat.unscaled <- chol2inv(object$R[seq_len(p), seq_len(p), drop = FALSE])
        dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
        covmat <- dispersion * covmat.unscaled
        s.err <- sqrt(diag(covmat))
        tvalue <- coef.p / s.err
        dn <- c("Estimate", "Std. Error")
        if (!est.disp) {
            coef.table <- cbind(coef.p, s.err, tvalue, 2 * stats::pnorm(-abs(tvalue)))
            dimnames(coef.table) <- list(names(coef.p), c(dn, "z value", "Pr(>|z|)"))
        } else if (df.r > 0) {
            coef.table <- cbind(coef.p, s.err, tvalue, 2 * stats::pt(-abs(tvalue), df.r))
            dimnames(coef.table) <- list(names(coef.p), c(dn, "t value", "Pr(>|t|)"))
        } else {
            coef.table <- cbind(coef.p, NaN, NaN, NaN)
            dimnames(coef.table) <- list(names(coef.p), c(dn, "t value", "Pr(>|t|)"))
        }
    } else {
        coef.table <- matrix(NA_real_, 0L, 4L, dimnames = list(NULL, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")))
        covmat.unscaled <- covmat <- matrix(NA_real_, 0L, 0L)
    }
    ans <- c(object[c("call", "terms", "family", "deviance", "aic", "contrasts", "df.residual",
                      "null.deviance", "df.null", "iter", "na.action", "compressed", "patterns",
                      "rows", "method")],
             list(coefficients = coef.table, aliased = aliased, dispersion = dispersion,
                  df = c(object$rank, df.r, length(aliased)), cov.unscaled = covmat.unscaled,
                  cov.scaled = covmat))
    class(ans) <- "summary.design_glm"
    ans
}

#' @export
print.summary.design_glm <- function(x, digits = max(3L, getOption("digits") - 3L),
                                     signif.stars = getOption("show.signif.stars"), ...) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
    dmfit_print_rows(x)
    dmfit_print_coefmat(x, digits, signif.stars, ...)
    cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ", format(x$dispersion),
        ")\n\n", apply(cbind(paste(format(c("Null", "Residual"), justify = "right"), "deviance:"),
                             format(unlist(x[c("null.deviance", "deviance")]), digits = max(5L, digits + 1L)),
                             " on", format(unlist(x[c("df.null", "df.residual")])), " degrees of freedom\n"),
                       1L, paste, collapse = " "), sep = "")
    if (nzchar(mess <- stats::naprint(x$na.action))) cat("  (", mess, ")\n", sep = "")
    cat("AIC: ", format(x$aic, digits = max(4L, digits + 1L)), "\n\n",
        "Number of Fisher Scoring iterations: ", x$iter, "\n\n", sep = "")
    invisible(x)
}

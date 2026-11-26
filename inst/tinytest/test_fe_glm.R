# Fixed effects in design_glm() (y ~ x | fe) against fixest::feglm() and
# fepois(): coefficients, iid and clustered covariance matrices, number of
# observations (rows removed: singletons and levels whose FE is infinite),
# deviance, log-likelihood, predictions and residuals, with fixest's
# default small sample corrections. Every case runs in memory and on
# remote data (a DuckDB table, and Parquet in the full tier). fixest
# demeans iteratively; tight glm.tol and fixef.tol make its results exact
# enough for a 1e-8 comparison.

if (!requireNamespace("fixest", quietly = TRUE)) exit_file("fixest not installed")
source("helpers.R")
library(design.matrix)
remote <- requireNamespace("DBI", quietly = TRUE) && requireNamespace("duckdb", quietly = TRUE)

TOL <- 1e-8
CTRL <- list(epsilon = 1e-10)

fx_fit <- function(fml, d, family, weights = NULL, offset = NULL, ...) {
    args <- list(fml, d, family = family, glm.tol = 1e-10, fixef.tol = 1e-11, notes = FALSE, warn = FALSE, ...)
    if (!is.null(weights)) args$weights <- stats::as.formula(paste0("~", weights))
    if (!is.null(offset)) args$offset <- stats::as.formula(paste0("~", offset))
    do.call(fixest::feglm, args)
}

plain <- function(V) {
    V <- unclass(as.matrix(V))
    attributes(V) <- list(dim = dim(V))
    V
}

if (remote) {
    con <- DBI::dbConnect(duckdb::duckdb())
    pq_dir <- tempfile("feglm_pq")
    dir.create(pq_dir)
}

# Write a data set to DuckDB (table) and to a Parquet file (NULL without
# DuckDB).
put <- function(d, name) if (remote) write_fe_source(con, pq_dir, d, name)

# Fits of the same model in memory and on remote data (`weights` and
# `offset` name columns; `modes`: the remote modes).
fe_fits <- function(fml, family, d, src, vcov = "iid", weights = NULL, offset = NULL,
                    modes = remote_modes, ...) {
    marg <- list(fml, family, d, vcov = vcov, control = CTRL, ...)
    if (!is.null(weights)) marg$weights <- d[[weights]]
    if (!is.null(offset)) marg$offset <- log(d[[offset]])
    fits <- list(memory = do.call(design_glm, marg))
    if (remote && !is.null(src)) {
        rfit <- function(source, ...) {
            a <- list(fml, family, source, vcov = vcov, control = CTRL, ...)
            ## (as columns: `weights` a name, `offset` an SQL expression)
            if (!is.null(weights)) a$weights <- weights
            if (!is.null(offset)) a$offset <- DBI::SQL(sprintf("ln(%s)", offset))
            suppressMessages(do.call(design_glm, a))
        }
        if ("sql" %in% modes) fits$sql <- rfit(con, table = src$table, ...)
        if ("parquet" %in% modes) fits$parquet <- rfit(src$parquet, ...)
    }
    fits
}

# Compare the fits of a model with a feglm() fit.
# (ll = FALSE: fixest's log-likelihood of families with a dispersion
# differs from glm()'s, which design_glm() follows)
check_fe <- function(fml, family, d, src, cluster = NULL, weights = NULL, offset = NULL, ...,
                     modes = remote_modes, ll = TRUE, info = "") {
    vc <- if (is.null(cluster)) "iid" else cluster
    fits <- fe_fits(fml, family, d, src, vc, weights, offset, modes = modes, ...)
    b <- fx_fit(fml, d, family, weights = weights, offset = if (!is.null(offset)) sprintf("log(%s)", offset))
    fb <- unname(predict(b, d, type = "response"))
    for (nm in names(fits)) {
        a <- fits[[nm]]
        inf <- paste(info, nm)
        ok <- !is.na(coef(a))
        expect_equal(names(coef(a))[ok], names(coef(b)), info = paste(inf, "names"))
        expect_equal(unname(coef(a)[ok]), unname(coef(b)), tolerance = TOL, info = paste(inf, "coef"))
        expect_equal(nobs(a), nobs(b), info = paste(inf, "nobs"))
        expect_equal(plain(vcov(a)[ok, ok]),
                     plain(if (is.null(cluster)) vcov(b, vcov = "iid") else vcov(b, cluster = cluster)),
                     tolerance = TOL, info = paste(inf, "vcov"))
        expect_equal(deviance(a), deviance(b), tolerance = TOL, info = paste(inf, "deviance"))
        if (ll) {
            expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)), tolerance = TOL, info = paste(inf, "logLik"))
            expect_equal(BIC(a), BIC(b), tolerance = TOL, info = paste(inf, "BIC"))
        }
        if (nm == "memory") {
            expect_equal(unname(predict(a, d, type = "response")), fb, tolerance = TOL, info = paste(inf, "predict"))
            r <- resid(a, type = "response")
            expect_equal(unname(r[!is.na(r)]), unname(resid(b, type = "response")), tolerance = TOL,
                         info = paste(inf, "resid"))
        } else {
            p <- as.data.frame(predict(a, con, table = src$table, type = "response"))
            p <- p[order(p$.row), ]
            expect_equal(p$.fitted, fb, tolerance = TOL, info = paste(inf, "predict"))
        }
    }
    invisible(fits)
}

# remote fits of FE GLMs cost one to several seconds each: the default tier
# fits three cases remotely (Poisson with weights, an offset and clusters;
# logit with a factor; missing values), the full tier all of them (the
# families and offset() terms are also fitted remotely without FE in
# test_parity_glm.R)
full_modes <- if (full_tests) remote_modes else character()

set.seed(20261126)
n <- 1500
d <- data.frame(
    a = sample(150, n, TRUE),
    b = sample(15, n, TRUE),
    c = sample(5, n, TRUE),
    cl = sample(25, n, TRUE),
    x1 = rnorm(n), x2 = rnorm(n),
    f = factor(sample(letters[1:3], n, TRUE)),
    w = runif(n, 0.3, 2),
    e = runif(n, 1, 4))
eta <- 0.4 * d$x1 - 0.3 * d$x2 + d$a / 150 - d$b / 15
d$y <- rpois(n, exp(eta) * d$e)
d$yb <- rbinom(n, 1, plogis(eta + as.integer(d$f) / 3 - 0.5))
# levels whose FE is infinite: all-0 counts, all-0 and all-1 binary responses
d$y[d$a %% 19 == 0] <- 0
d$yb[d$a %% 23 == 0] <- 1
d$yb[d$a %% 29 == 0] <- 0
src <- put(d, "dglm")

# Poisson --------------------------------------------------------------------
check_fe(y ~ x1 + x2 | a, poisson, d, src, cluster = ~cl, modes = full_modes, info = "poisson one-way, cluster not nested")
check_fe(y ~ x1 + x2 | a + b, poisson, d, src, cluster = ~a, weights = "w", offset = "e",
         info = "poisson two-way, weights, offset, cluster nested")
check_fe(y ~ x1 + x2 + offset(log(e)) | a + b + c, poisson, d, src, cluster = ~cl, modes = full_modes,
         info = "poisson three-way, offset() term")

# Logit, probit --------------------------------------------------------------
check_fe(yb ~ x1 + x2 + f | a + b, binomial, d, src, cluster = ~cl, info = "logit two-way")
# (non-integer weights: glm()'s warning)
suppressWarnings(check_fe(yb ~ x1 + x2 | a + b + c, binomial, d, src, weights = "w", modes = full_modes, info = "logit three-way weighted"))
check_fe(yb ~ x1 + x2 | a + b, binomial("probit"), d, src, cluster = ~b, modes = full_modes, info = "probit two-way")
check_fe(yb ~ x1 * f | a^c + b, binomial("probit"), d, src, cluster = ~cl,
         modes = intersect(full_modes, "sql"), info = "probit interacted FE")

# Other families (fixest has the same estimator for any family) -------------
check_fe(y ~ x1 + x2 | a + b, quasipoisson, d, src, modes = intersect(full_modes, "sql"), info = "quasipoisson")
check_fe(y ~ x1 + x2 | a + b, binomial("cloglog"), within(d, y <- pmin(y, 1)), NULL, info = "cloglog")
dg <- d
dg$yg <- rgamma(n, 2, 2 / exp(0.2 * d$x1 + d$b / 15))
ga <- check_fe(yg ~ x1 + x2 | a + b, Gamma("log"), dg, NULL, cluster = ~cl, ll = FALSE, info = "Gamma log")
# its log-likelihood and AIC are those of glm() with the FE as dummies
gg <- glm(yg ~ x1 + x2 + factor(a) + factor(b), Gamma("log"), dg, control = CTRL)
expect_equal(as.numeric(logLik(ga$memory)), as.numeric(logLik(gg)), tolerance = 1e-7)
expect_equal(AIC(ga$memory), AIC(gg), tolerance = 1e-7)
expect_equal(deviance(ga$memory), deviance(gg), tolerance = 1e-7)

# Rows removed ---------------------------------------------------------------
# perfect fit (fixest's default fixef.rm = "perfect_fit"): rows of levels
# whose response is all 0 (Poisson) or all 0 / all 1 (binomial), and
# singletons, iteratively
pf <- design_glm(y ~ x1 | a + b, poisson, d)
bf <- fx_fit(y ~ x1 | a + b, d, poisson)
expect_equal(nobs(pf), nobs(bf))
expect_equal(pf$nperfect + pf$nsingletons, n - nobs(bf))
expect_true(pf$nperfect > 0)
pb <- design_glm(yb ~ x1 | a, binomial, d)
expect_equal(nobs(pb), nobs(fx_fit(yb ~ x1 | a, d, binomial())))
# the other values of fixef.rm, as fixest's
for (rm in c("singleton", "infinite_coef", "none")) {
    a <- design_glm(y ~ x1 + x2 | a + b, poisson, d, fixef.rm = rm, control = CTRL)
    b <- fx_fit(y ~ x1 + x2 | a + b, d, poisson, fixef.rm = if (rm == "singleton") "singletons" else rm)
    expect_equal(nobs(a), nobs(b), info = rm)
    expect_equal(unname(coef(a)), unname(coef(b)), tolerance = 1e-6, info = rm)
}
if (remote) {
    a <- design_glm(yb ~ x1 | a + c, binomial, con, table = "dglm")
    b <- design_glm(yb ~ x1 | a + c, binomial, d)
    expect_equal(c(a$nperfect, a$nsingletons, nobs(a)), c(b$nperfect, b$nsingletons, nobs(b)))
}
# a level that is both a singleton and all 0 counts as perfectly fit; a
# singleton with y > 0 as a singleton (the masks of each level count matrix
# are computed once, design.matrix#6y4d)
set.seed(11)
dp <- data.frame(a = c(1, 12, rep(2:11, each = 5)))
dp$x1 <- rnorm(nrow(dp))
dp$y <- c(0, 3, rpois(50, 2) + 1)
srcp <- put(dp, "dglm_prec")
for (rm in c("perfect_fit", "singleton", "infinite_coef")) {
    fp <- design_glm(y ~ x1 | a, poisson, dp, fixef.rm = rm)
    want <- switch(rm, perfect_fit = c(1, 1, 50), singleton = c(0, 2, 50), infinite_coef = c(1, 0, 51))
    expect_equal(c(fp$nperfect, fp$nsingletons, nobs(fp)), want, info = paste("precedence", rm))
    if (remote) {
        fr <- design_glm(y ~ x1 | a, poisson, con, table = srcp$table, fixef.rm = rm)
        expect_equal(c(fr$nperfect, fr$nsingletons, nobs(fr)), want, info = paste("precedence remote", rm))
        expect_equal(coef(fr), coef(fp), tolerance = TOL, info = paste("precedence remote coef", rm))
    }
}

# Missing values: na.omit, na.exclude ----------------------------------------
dn <- d
dn$x1[c(3, 50)] <- NA
dn$a[7] <- NA
dn$y[9] <- NA
srcn <- put(dn, "dglm_na")
check_fe(y ~ x1 + x2 | a + b, poisson, dn, srcn, cluster = ~cl, info = "missing values")
ax <- design_glm(y ~ x1 + x2 | a + b, poisson, dn, na.action = na.exclude)
rx <- residuals(ax, type = "response")
expect_equal(length(rx), n)
expect_true(all(is.na(rx[c(3, 7, 9, 50)])))
expect_equal(length(fitted(ax)), n)
if (remote) {
    rr <- design_glm(y ~ x1 + x2 | a + b, poisson, con, table = "dglm_na", na.action = na.exclude)
    fr <- as.data.frame(fitted(rr))
    expect_equal(nrow(fr), n)
    fr <- fr[order(fr$.row), ]
    expect_equal(fr$.fitted, unname(fitted(ax)), tolerance = TOL)
}

# Methods --------------------------------------------------------------------
a <- design_glm(y ~ x1 + x2 | a + b, poisson, d, vcov = ~cl, control = CTRL)
b <- fx_fit(y ~ x1 + x2 | a + b, d, poisson, vcov = ~cl)
expect_inherits(a, "design_glm_fe")
expect_equal(unname(summary(a)$coefficients), unname(fixest::coeftable(b)[, 1:4]), tolerance = TOL)
expect_equal(AIC(a), AIC(b), tolerance = TOL)
# fixef(): differences within a connected set are those of fixest
fa <- a$fixef$a
fb <- fixest::fixef(b)$a
expect_equal(unname(fa - fa[1]), unname(fb[names(fa)] - fb[names(fa)[1]]), tolerance = 1e-7)
# link predictions, new levels give NA
nd <- d[1:5, ]
nd$a[2] <- 9999
pl <- predict(a, nd)
expect_true(is.na(pl[2]))
expect_equal(unname(pl[-2]), unname(predict(b, nd[-2, ], type = "link")), tolerance = TOL)
expect_equal(unname(residuals(a, "pearson")), unname(resid(b, "pearson")), tolerance = TOL)
expect_equal(unname(residuals(a, "deviance")), unname(resid(b, "deviance")), tolerance = TOL)
expect_true(length(a$fe_iter) == a$iter)
expect_error(predict(a, nd, se.fit = TRUE))
expect_stdout(print(a), "Fixed effects")
expect_stdout(print(summary(a)), "clustered by cl")
# the Gaussian identity GLM is the linear model with FE
g <- design_glm(x2 ~ x1 | a + b, gaussian, d)
l <- design_lm(x2 ~ x1 | a + b, d)
expect_equal(coef(g), coef(l), tolerance = TOL)

# What is refused (in memory and remotely) -----------------------------------
expect_error(design_glm(y ~ x1 | a, poisson, d, compress = TRUE), "compress")
expect_error(design_glm(y ~ x1 | a, poisson, d, start = 0), "start")
expect_error(design_glm(y ~ x1 | a, poisson, d, na.action = na.fail), "na.omit or na.exclude")
expect_error(design_glm(y ~ x1, poisson, d, vcov = ~cl), "requires fixed effects")
expect_error(design_glm(y ~ x1 | a, poisson, d, vcov = ~a + b), "one-way")
expect_error(design_glm(y ~ x1 | a, poisson, d, fixef.rm = "x"), "fixef.rm")
if (remote) {
    expect_error(design_glm(y ~ x1 | a, poisson, con, table = "dglm", compress = TRUE), "compress")
    expect_error(design_glm(y ~ x1 | a, poisson, con, table = "dglm", start = 0), "start")
}

# Remote query size stays bounded (probit: DuckDB expands the qnorm and
# pnorm macros inline; frozen columns keep queries short; full tier) --------
if (remote && full_tests) {
    qlen <- integer()
    withCallingHandlers({
        op <- options(design.matrix.remote_debug = TRUE)
        fp <- design_glm(yb ~ x1 + x2 | a + b + c, binomial("probit"), con, table = "dglm")
        options(op)
    }, message = function(m) {
        qlen <<- c(qlen, nchar(conditionMessage(m)))
        invokeRestart("muffleMessage")
    })
    expect_true(max(qlen) < 2e4)
    expect_equal(coef(fp), coef(design_glm(yb ~ x1 + x2 | a + b + c, binomial("probit"), d)), tolerance = TOL)
}
if (remote) DBI::dbDisconnect(con)

# Cases shared by test_glm.R (in memory) and test_remote_glm.R (remote):
# the same models, each compared with glm() on the same rows.

glm_data <- function(n = 1500, seed = 42) {
    set.seed(seed)
    d <- data.frame(
        x = rnorm(n),
        xd = sample(1:4, n, TRUE),
        f = sample(c("b", "a", "c", "d"), n, TRUE),
        g = sample(c("u", "v", "w"), n, TRUE),
        w = runif(n, 0.5, 2),
        e = runif(n, 1, 3),
        id = seq_len(n))
    eta <- 0.3 * d$x + as.integer(factor(d$f)) / 4 - 0.5
    d$yb <- rbinom(n, 1, plogis(eta))
    d$yp <- rpois(n, d$e * exp(0.2 * d$xd + as.integer(factor(d$g)) / 5 - 1))
    d$yg <- rgamma(n, 2, scale = exp(0.1 * d$xd + as.integer(factor(d$f)) / 10) / 2)
    d$s <- rbinom(n, 5, plogis(0.2 * d$xd - 0.5))
    d$fl <- 5 - d$s
    d$prop <- d$s / 5
    d$five <- 5
    d$le <- log(d$e)
    d$xd2 <- 2 * d$xd   # aliased with xd
    d$fy <- ifelse(d$yb == 1, "yes", "no")
    d$sep <- as.integer(d$x > 0)  # complete separation on x
    d$yb[5] <- NA
    d$x[7] <- NA
    d$w[11] <- 0
    d
}

# name, formula, family, weights column, offset column
glm_cases <- list(
    list("logit x*f", yb ~ x * f, binomial()),
    list("logit factors", yb ~ xd * f + g, binomial()),
    list("logit weights", yb ~ xd + g, binomial(), "w"),
    list("probit", yb ~ xd + f, binomial("probit")),
    list("cloglog", yb ~ xd + f, binomial("cloglog")),
    list("cauchit", yb ~ x + g, binomial("cauchit")),
    list("binomial cbind", cbind(s, fl) ~ xd + g, binomial()),
    list("binomial proportions", prop ~ xd + g, binomial(), "five"),
    list("poisson offset()", yp ~ xd + g + offset(log(e)), poisson()),
    list("poisson offset arg", yp ~ x + g, poisson(), NULL, "le"),
    list("poisson no intercept", yp ~ 0 + g + xd, poisson()),
    list("gaussian", yg ~ xd + f, gaussian()),
    list("gaussian log", yg ~ x + f, gaussian("log"), "w"),
    list("Gamma inverse", yg ~ x + f, Gamma()),
    list("Gamma log", yg ~ xd + f, Gamma("log"), "w"),
    list("inverse.gaussian", yg ~ xd + f, inverse.gaussian()),
    list("quasibinomial", yb ~ x + f, quasibinomial()),
    list("quasipoisson", yp ~ xd + g, quasipoisson()),
    list("aliased", yb ~ xd + f + xd2, binomial()))

glm_ref <- function(cs, d) {
    eval(bquote(glm(.(cs[[2]]), family = .(cs[[3]]), data = d,
                    weights = .(if (length(cs) >= 4 && !is.null(cs[[4]])) as.name(cs[[4]])),
                    offset = .(if (length(cs) >= 5 && !is.null(cs[[5]])) as.name(cs[[5]])))))
}

glm_same <- function(a, b, info, tol = 1e-8) {
    same <- function(u, v) isTRUE(all.equal(unname(u), unname(v), tolerance = tol))
    expect_identical(is.na(coef(a)), is.na(coef(b)), info = paste(info, "aliased"))
    expect_true(same(coef(a), coef(b)), info = paste(info, "coef"))
    expect_true(same(vcov(a), vcov(b)), info = paste(info, "vcov"))
    expect_true(same(deviance(a), deviance(b)), info = paste(info, "deviance"))
    expect_true(same(a$null.deviance, b$null.deviance), info = paste(info, "null deviance"))
    expect_true(same(AIC(a), AIC(b)) || (is.na(AIC(a)) && is.na(AIC(b))), info = paste(info, "AIC"))
    expect_equal(df.residual(a), df.residual(b), info = paste(info, "df"))
    expect_equal(a$df.null, b$df.null, info = paste(info, "df null"))
    expect_equal(nobs(a), nobs(b), info = paste(info, "nobs"))
    expect_equal(a$iter, b$iter, info = paste(info, "iter"))
    expect_true(same(summary(a)$dispersion, summary(b)$dispersion), info = paste(info, "dispersion"))
}

# Warnings of an expression, and its value.
glm_warnings <- function(expr) {
    w <- character()
    v <- withCallingHandlers(expr, warning = function(cond) {
        w <<- c(w, conditionMessage(cond))
        invokeRestart("muffleWarning")
    })
    list(value = v, warnings = w)
}

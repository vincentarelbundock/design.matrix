source("helpers.R")
library(design.matrix)

set.seed(4)
n <- 200
d <- data.frame(y = rnorm(n), x = rnorm(n), z = runif(n, 1, 2),
                f = factor(sample(letters[1:3], n, TRUE)),
                ch = sample(c("u", "v", "w"), n, TRUE),
                g = factor(rep(1:20, length.out = n)))
d$yb <- stats::plogis(d$y)
d$yo <- cut(d$y, 3, ordered_result = TRUE)
d$yc <- rpois(n, 3)
# new data: some levels missing, rows reordered, no response
nd <- d[c(7, 3, 12, 40, 41, 55, 90, 91, 120, 150), c("x", "z", "f", "ch", "g")]
nd <- nd[nd$f != "c", ]
rhs <- ~ poly(x, 2) * f + log(z) + ch

reference <- function(plan, newdata) {
    mf <- stats::model.frame(plan$terms, newdata, xlev = plan$xlevels,
                             na.action = stats::na.pass)
    X <- stats::model.matrix(plan$terms, mf, contrasts.arg = plan$contrasts)
    if (isTRUE(plan$drop_intercept)) {
        att <- attributes(X)
        X <- X[, -1L, drop = FALSE]
        attr(X, "assign") <- att$assign[-1L]
        if (!is.null(att$contrasts)) attr(X, "contrasts") <- att$contrasts
    }
    X
}

check_fit <- function(fit, beta, pred, label, ...) {
    p <- design_matrix(fit, newdata = NULL, ...)
    expect_inherits(p, "design_plan", info = label)
    X <- as.matrix(design_matrix(p, nd))
    expect_true(strict_identical(X, reference(p, nd)), info = label)
    expect_equal(unname(drop(X %*% beta)), unname(as.numeric(pred)),
                 tolerance = 1e-10, info = label)
    # the second call takes the fast path prepared by the first
    D <- design_matrix(p, nd)
    expect_identical(D$path, "fast", info = label)
    expect_true(strict_identical(as.matrix(D), X), info = label)
}

# stats
fit <- lm(update(rhs, y ~ .), d)
check_fit(fit, coef(fit), predict(fit, nd), "lm")
fit <- glm(update(rhs, yc ~ .), d, family = poisson())
check_fit(fit, coef(fit), predict(fit, nd), "glm")
fit <- lm(update(rhs, y ~ .), d, contrasts = list(f = "contr.sum", ch = "contr.helmert"))
check_fit(fit, coef(fit), predict(fit, nd), "lm with contrasts")
fit <- lm(y ~ x, d)
expect_null(design.matrix:::design_parts(fit)$contrasts)
check_fit(fit, coef(fit), predict(fit, nd), "lm without factors")

if (requireNamespace("MASS", quietly = TRUE)) {
    fit <- MASS::rlm(update(rhs, y ~ .), d)
    check_fit(fit, coef(fit), predict(fit, nd), "rlm")
}

if (requireNamespace("lme4", quietly = TRUE)) {
    fit <- suppressMessages(lme4::lmer(y ~ poly(x, 2) * f + log(z) + ch + (1 | g), d))
    check_fit(fit, lme4::fixef(fit), predict(fit, nd, re.form = NA), "lmer")
}

if (requireNamespace("nlme", quietly = TRUE)) {
    # nlme keeps no levels for character variables: `ch` as a factor
    dn <- transform(d, ch = factor(ch))
    fit <- nlme::lme(update(rhs, y ~ .), random = ~ 1 | g, data = dn)
    check_fit(fit, nlme::fixef(fit), predict(fit, nd, level = 0), "lme")
    fit <- nlme::gls(update(rhs, y ~ .), data = dn)
    check_fit(fit, coef(fit), predict(fit, nd), "gls")
}

if (requireNamespace("glmmTMB", quietly = TRUE)) {
    fit <- glmmTMB::glmmTMB(y ~ poly(x, 2) * f + log(z) + ch + (1 | g), d)
    check_fit(fit, glmmTMB::fixef(fit)$cond, predict(fit, nd, re.form = NA), "glmmTMB")
}

if (requireNamespace("betareg", quietly = TRUE)) {
    fit <- betareg::betareg(yb ~ poly(x, 2) * f + log(z) + ch | x + f, d)
    check_fit(fit, coef(fit, model = "mean"), predict(fit, nd, type = "link"),
              "betareg mean")
    phi <- predict(fit, nd, type = "precision")
    check_fit(fit, coef(fit, model = "precision"), fit$link$precision$linkfun(phi),
              "betareg precision", component = "precision")
}

if (requireNamespace("ordinal", quietly = TRUE)) {
    # clm parameterises the location as -X beta, without an intercept
    fit <- ordinal::clm(update(rhs, yo ~ .), data = d)
    X <- as.matrix(design_matrix(design.matrix:::design_plan(fit), nd))
    expect_true(strict_identical(X, reference(design.matrix:::design_plan(fit), nd)))
    eta <- drop(X[, names(fit$beta), drop = FALSE] %*% fit$beta)
    lp <- predict(fit, nd, type = "linear.predictor")$eta1[, 1]
    expect_equal(unname(fit$alpha[1] - eta), unname(lp), tolerance = 1e-10)
}

# more training variables (after `nd` is drawn, so the draws above are unchanged)
d$zi <- d$x + rnorm(n)
d$time <- rexp(n, exp(0.3 * d$x))
d$status <- rbinom(n, 1, 0.8)
d$yz <- d$yc * rbinom(n, 1, 0.7)
nd$zi <- d[rownames(nd), "zi"]

if (requireNamespace("lme4", quietly = TRUE)) {
    d$y01 <- rbinom(n, 1, 0.5)
    fit <- suppressMessages(lme4::glmer(y01 ~ poly(x, 2) + f + ch + (1 | g), d, family = binomial()))
    check_fit(fit, lme4::fixef(fit), predict(fit, nd, re.form = NA), "glmer")
}

if (requireNamespace("glmmTMB", quietly = TRUE)) {
    fit <- suppressWarnings(glmmTMB::glmmTMB(yz ~ poly(x, 2) + f + (1 | g),
                                             ziformula = ~ poly(x, 2) + ch,
                                             family = poisson(), data = d))
    check_fit(fit, glmmTMB::fixef(fit)$zi, predict(fit, nd, type = "zlink", re.form = NA),
              "glmmTMB zi", component = "zi")
    expect_error(design.matrix:::design_plan(fit, component = "disp"), "no fixed effects")
    # contrasts the user did not supply are fixed by the plan
    expect_identical(design.matrix:::design_parts(fit)$contrasts, list(f = "contr.treatment"))
}

if (requireNamespace("MASS", quietly = TRUE)) {
    # polr: the intercept column is not used, the cut points take its place
    fit <- MASS::polr(update(rhs, yo ~ .), data = d)
    p <- design.matrix:::design_plan(fit)
    X <- as.matrix(design_matrix(p, nd))
    expect_true(strict_identical(X, reference(p, nd)))
    eta <- drop(X[, names(coef(fit))] %*% coef(fit))
    expect_equal(unname(stats::plogis(fit$zeta[1] - eta)),
                 unname(predict(fit, nd, type = "probs")[, 1]), tolerance = 1e-8)
    fit <- MASS::glm.nb(update(rhs, yc ~ .), d)
    check_fit(fit, coef(fit), predict(fit, nd), "glm.nb")
}

if (requireNamespace("nnet", quietly = TRUE)) {
    fit <- nnet::multinom(update(rhs, yo ~ .), d, trace = FALSE)
    p <- design.matrix:::design_plan(fit)
    X <- as.matrix(design_matrix(p, nd))
    expect_true(strict_identical(X, reference(p, nd)))
    eta <- X %*% t(coef(fit))
    pr <- cbind(1, exp(eta)) / rowSums(cbind(1, exp(eta)))
    expect_equal(unname(pr), unname(predict(fit, nd, type = "probs")), tolerance = 1e-6)
}

if (requireNamespace("survival", quietly = TRUE)) {
    Surv <- survival::Surv
    # coxph: the intercept column is dropped, the factors are coded as with it
    fit <- survival::coxph(Surv(time, status) ~ poly(x, 2) * f + ch, d)
    lp <- predict(fit, nd, type = "lp", reference = "zero")
    check_fit(fit, coef(fit), lp, "coxph")
    X <- as.matrix(design_matrix(design.matrix:::design_plan(fit), nd))
    expect_false("(Intercept)" %in% colnames(X))
    expect_true(strict_identical(X, stats::model.matrix(fit, data = nd)))
    # strata() and cluster() terms are removed
    fit <- survival::coxph(Surv(time, status) ~ poly(x, 2) + strata(f) + ch + cluster(g), d)
    check_fit(fit, coef(fit), predict(fit, nd, type = "lp", reference = "zero"),
              "coxph with strata and cluster")
    X <- as.matrix(design_matrix(design.matrix:::design_plan(fit), nd))
    expect_equal(unclass(X)[, ], unclass(stats::model.matrix(fit, data = nd))[, ],
                 check.attributes = FALSE)
    fit <- survival::coxph(Surv(time, status) ~ x * strata(f), d)
    expect_error(design.matrix:::design_plan(fit), "strata")

    fit <- survival::survreg(Surv(time, status) ~ poly(x, 2) * f + log(z) + ch, d)
    check_fit(fit, coef(fit), predict(fit, nd, type = "lp"), "survreg")
}

if (requireNamespace("estimatr", quietly = TRUE)) {
    fit <- estimatr::lm_robust(update(rhs, y ~ .), d)
    check_fit(fit, coef(fit), predict(fit, nd), "lm_robust")
    # terms() of iv_robust models folds the instruments into the formula.
    # predict() of estimatr and ivreg recomputes poly() on new data, so the
    # rows of `nd`, which are training rows, are checked against the fit.
    rows <- match(rownames(nd), rownames(d))
    fit <- estimatr::iv_robust(y ~ poly(x, 2) + f + ch | poly(zi, 2) + f + ch, d)
    check_fit(fit, coef(fit), fit$fitted.values[rows], "iv_robust")
}

if (requireNamespace("ivreg", quietly = TRUE)) {
    rows <- match(rownames(nd), rownames(d))
    fit <- ivreg::ivreg(y ~ poly(x, 2) + f + ch | poly(zi, 2) + f + ch, data = d)
    check_fit(fit, coef(fit), fitted(fit)[rows], "ivreg")
}

if (requireNamespace("pscl", quietly = TRUE)) {
    # the component terms of pscl models carry no predvars
    fit <- pscl::zeroinfl(yz ~ poly(x, 2) + f | poly(x, 2) + ch, d)
    check_fit(fit, coef(fit, model = "count"), log(predict(fit, nd, type = "count")),
              "zeroinfl count")
    check_fit(fit, coef(fit, model = "zero"), stats::qlogis(predict(fit, nd, type = "zero")),
              "zeroinfl zero", component = "zero")
    fit <- pscl::hurdle(yz ~ poly(x, 2) + f | x + ch, d)
    check_fit(fit, coef(fit, model = "count"), log(predict(fit, nd, type = "count")),
              "hurdle count")
}

# a plan without training data prepares its fast path on first use
# (plans are memoized on their pieces, so an equal model fitted above would
# hand back its warmed plan: start from an empty cache)
design.matrix:::dm_plan_cache_clear()
fit <- lm(update(rhs, y ~ .), d)
p <- design.matrix:::design_plan(fit)
expect_null(p$fast)
expect_stdout(print(p), "on first use")
expect_identical(design_matrix(p, nd)$path, "general")
expect_stdout(print(p), "fast path: yes")
expect_identical(design_matrix(p, d[1:3, ])$path, "fast")
# the same matrices as a plan built from the formula and training data
pf <- design.matrix:::design_plan(update(rhs, y ~ .), d)
expect_true(strict_identical(as.matrix(design_matrix(p, d[5:9, ])), as.matrix(design_matrix(pf, d[5:9, ]))))

# the lazily prepared fast path gives up on a relevelled factor and hands
# the variables it evaluated to the general path
fit <- lm(y ~ log(z) + f, d)
p <- design.matrix:::design_plan(fit)
invisible(design_matrix(p, d[1:5, ]))
nd2 <- d[1:6, ]
nd2$f <- factor(as.character(nd2$f), levels = rev(levels(d$f)))
D <- design_matrix(p, nd2)
expect_identical(D$path, "general")
expect_true(strict_identical(as.matrix(D), reference(p, nd2)))

# the response on the right-hand side is refused for fitted models too
expect_error(design.matrix:::design_plan(suppressWarnings(lm(y ~ y + x, d))), "right-hand side")

# data that warns does not prepare the fast path; the next data does
fit <- lm(y ~ log(x) + f, d[d$x > 0, ])
p <- design.matrix:::design_plan(fit)
expect_warning(design_matrix(p, d[1:20, ]))
expect_null(p$lazy$fast)
expect_identical(design_matrix(p, d[d$x > 0, ][1:5, ])$path, "general")
expect_identical(design_matrix(p, d[d$x > 0, ][6:9, ])$path, "fast")

# errors on new data reach the caller unchanged, and the plan still works
fit <- lm(y ~ x + ch, d)
p <- design.matrix:::design_plan(fit)
bad <- d[1:3, ]
bad$ch[1] <- "zzz"
expect_true(same_outcome(capture(function() as.matrix(design_matrix(p, bad))),
                         capture(function() reference(p, bad))))
expect_true(strict_identical(as.matrix(design_matrix(p, d[1:3, ])), reference(p, d[1:3, ])))

# a custom method; malformed parts are rejected
design_parts.my_model <- function(model, ...) model$parts
registerS3method("design_parts", "my_model", design_parts.my_model, asNamespace("design.matrix"))
obj <- structure(list(parts = list(terms = terms(y ~ x + f), xlevels = list(f = letters[1:3]),
                                   contrasts = list(f = "contr.sum"))),
                 class = "my_model")
p <- design.matrix:::design_plan(obj)
expect_true(strict_identical(as.matrix(design_matrix(p, nd)), reference(p, nd)))
obj$parts$terms <- y ~ x
expect_error(design.matrix:::design_plan(obj), "terms")

# without stored levels, a character variable is encoded from each
# batch's own values: the first batch must not fix its levels for later ones
obj <- structure(list(parts = list(terms = terms(~ x + ch + lg), xlevels = NULL,
                                   contrasts = NULL)), class = "my_model")
p <- design.matrix:::design_plan(obj)
b1 <- data.frame(x = 1:4, ch = c("a", "b", "c", "a"), lg = c(TRUE, FALSE, TRUE, TRUE))
b2 <- data.frame(x = 1:3, ch = c("a", "b", "a"), lg = c(TRUE, FALSE, TRUE))
b3 <- data.frame(x = 1:2, ch = c("a", "c"), lg = c(FALSE, TRUE))
for (b in list(b1, b2, b3, b1)) {
    expect_true(strict_identical(as.matrix(design_matrix(p, b)), reference(p, b)))
}

# an intercept-only model: one row per row of the new data
fit <- lm(y ~ 1, d)
p <- design.matrix:::design_plan(fit)
for (k in 1:2) {
    X <- as.matrix(design_matrix(p, nd))
    expect_true(strict_identical(X, reference(p, nd)))
    expect_identical(nrow(X), nrow(nd))
}

# contrasts the plan does not fix follow options("contrasts") at every call
obj <- structure(list(parts = list(terms = terms(~ x + f), xlevels = list(f = letters[1:3]),
                                   contrasts = NULL)), class = "my_model")
p <- design.matrix:::design_plan(obj)
invisible(as.matrix(design_matrix(p, nd)))
old_opts <- options(contrasts = c("contr.sum", "contr.poly"))
expect_true(strict_identical(as.matrix(design_matrix(p, nd)), reference(p, nd)))
options(old_opts)
expect_true(strict_identical(as.matrix(design_matrix(p, nd)), reference(p, nd)))

# on first use, each transformation is evaluated once
n_calls <- 0
count_id <- function(x) {
    n_calls <<- n_calls + 1
    x
}
fit <- lm(y ~ count_id(x) + f, d)
p <- design.matrix:::design_plan(fit)
n_calls <- 0
X <- as.matrix(design_matrix(p, nd))
expect_equal(n_calls, 1)
expect_false(is.null(p$lazy$fast))
expect_true(strict_identical(X, reference(p, nd)))
expect_error(design.matrix:::design_plan(structure(list(), class = "not_a_model")), "terms\\(\\) method")

# plans of fitted models are memoized on their pieces
design.matrix:::dm_plan_cache_clear()
fit1 <- lm(update(rhs, y ~ .), d)
p1 <- design.matrix:::design_plan(fit1)
invisible(design_matrix(p1, nd))
expect_identical(design_matrix(fit1, newdata = nd[1:2, ])$path, "fast")
expect_true(identical(design.matrix:::design_plan(fit1)$lazy, p1$lazy))
# different contrasts: a different plan
fit2 <- lm(update(rhs, y ~ .), d, contrasts = list(f = "contr.sum"))
expect_false(identical(design.matrix:::design_plan(fit2)$lazy, p1$lazy))
expect_true(strict_identical(as.matrix(design_matrix(fit2, newdata = nd)),
                             model.matrix(delete.response(terms(fit2)),
                                          model.frame(delete.response(terms(fit2)), nd, xlev = fit2$xlevels),
                                          contrasts.arg = fit2$contrasts)))

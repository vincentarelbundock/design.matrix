source("helpers.R")
library(xtx)

# design_lm() must agree with lm(): coefficients (with the same NA pattern
# for aliased columns), vcov, sigma, degrees of freedom, number of
# observations, and predictions on new data.

same_fit <- function(a, b, tol = 1e-7) {
    isTRUE(all.equal(coef(a), coef(b), tolerance = tol)) &&
        identical(is.na(coef(a)), is.na(coef(b))) &&
        ## vcov.lm() warns on a perfect fit
        isTRUE(all.equal(suppressWarnings(vcov(a)), vcov(b), tolerance = tol)) &&
        isTRUE(all.equal(sigma(a), sigma(b), tolerance = tol)) &&
        identical(as.numeric(df.residual(a)), as.numeric(df.residual(b))) &&
        identical(as.numeric(nobs(a)), as.numeric(nobs(b)))
}

# Fit with lm() and design_lm() from the same arguments (the call is
# rebuilt, so that `weights`, `subset` and `offset` are evaluated as in a
# direct call).
both <- function(...) {
    cl <- match.call()
    cl[[1L]] <- quote(lm)
    a <- eval(cl, parent.frame())
    cl[[1L]] <- quote(design_lm)
    b <- eval(cl, parent.frame())
    list(a = a, b = b)
}

# Small examples, every method ---------------------------------------------

set.seed(1)
n <- 200
d <- data.frame(y = rnorm(n), x = rnorm(n),
                f = factor(sample(letters[1:4], n, TRUE)),
                g = factor(sample(1:3, n, TRUE)),
                w = runif(n), o = rnorm(n))
d$x[3] <- NA
d$w[5] <- 0
d$f2 <- d$f # aliased with f

check_all <- function(a, b, nd = d[1:10, ]) {
    expect_true(same_fit(a, b))
    expect_equal(residuals(b), residuals(a))
    expect_equal(fitted(b), fitted(a))
    expect_equal(deviance(b), deviance(a))
    expect_equal(logLik(b), logLik(a))
    expect_equal(AIC(b), AIC(a))
    expect_equal(coef(summary(b)), coef(summary(a)))
    expect_equal(summary(b)$sigma, summary(a)$sigma)
    expect_equal(summary(b)$r.squared, summary(a)$r.squared)
    expect_equal(summary(b)$adj.r.squared, summary(a)$adj.r.squared)
    expect_equal(summary(b)$fstatistic, summary(a)$fstatistic)
    expect_equal(summary(b)$aliased, summary(a)$aliased)
    expect_equal(predict(b, nd), predict(a, nd))
    expect_equal(predict(b), predict(a))
    pa <- suppressWarnings(predict(a, nd, se.fit = TRUE, interval = "confidence"))
    pb <- predict(b, nd, se.fit = TRUE, interval = "confidence")
    expect_equal(pb$fit, pa$fit)
    ## predict.lm() drops the names of se.fit for an intercept-only model
    expect_equal(unname(pb$se.fit), unname(pa$se.fit))
    expect_equal(pb$df, pa$df)
    expect_equal(pb$residual.scale, pa$residual.scale)
    ## same values and warnings, on new data and on the training rows
    pa <- capture(function() predict(a, nd, interval = "prediction"))
    pb <- capture(function() predict(b, nd, interval = "prediction"))
    expect_equal(pb$value, pa$value)
    expect_equal(pb$warnings, pa$warnings)
    ## (with na.exclude and weights, predict.lm() recycles the weights of
    ## the complete rows against the padded predictions: no reference)
    if (!inherits(a$na.action, "exclude")) {
        pa <- capture(function() predict(a, interval = "prediction"))
        pb <- capture(function() predict(b, interval = "prediction"))
        expect_equal(pb$value, pa$value)
        expect_equal(pb$warnings, pa$warnings)
    }
    expect_equal(formula(b), formula(a))
}

fits <- list(
    both(y ~ x * f, d),
    both(y ~ x * f + f2, d, weights = w),
    both(y ~ f * g, d, weights = w),
    both(y ~ f * g + offset(o), d),
    both(y ~ x + f + g + f2, d, offset = o, weights = w, na.action = na.exclude),
    both(y ~ 0 + f:g, d),
    both(y ~ 1, d),
    both(y ~ log(abs(x)) + poly(o, 2) + f, d, subset = g != "2"),
    both(y ~ x + f, d, subset = 1:150, weights = w, offset = 0.5 * o))
for (ab in fits) check_all(ab$a, ab$b)

# factor-only designs are compressed by default, and equal uncompressed fits
b <- design_lm(y ~ f * g, d, weights = w)
expect_true(b$compressed)
expect_equal(b$patterns, 12L)
expect_false(design_lm(y ~ x * f, d)$compressed)
u <- design_lm(y ~ f * g, d, weights = w, compress = FALSE)
expect_false(u$compressed)
expect_true(same_fit(b, u, tol = 1e-10))
expect_equal(summary(b)$r.squared, summary(u)$r.squared)
expect_true(design_lm(y ~ f * g, d, compress = TRUE)$compressed)

# the empty model
a <- lm(y ~ 0, d)
b <- design_lm(y ~ 0, d)
expect_identical(coef(b), coef(a))
expect_equal(vcov(b), vcov(a))
expect_equal(sigma(b), sigma(a))
expect_equal(df.residual(b), df.residual(a))
expect_equal(predict(b, d[1:3, ]), predict(a, d[1:3, ]))

# a perfect fit with no residual degrees of freedom
d4 <- data.frame(y = c(1, 2, 3), f = factor(c("a", "b", "c")))
expect_true(same_fit(lm(y ~ f, d4), design_lm(y ~ f, d4)))

# na.exclude pads residuals, fitted and predictions without newdata
b <- design_lm(y ~ x + f, d, na.action = na.exclude)
expect_true(is.na(residuals(b)[3]))
expect_equal(length(fitted(b)), n)
expect_equal(predict(b), predict(lm(y ~ x + f, d, na.action = na.exclude)))

# NA in new data gives NA predictions, as predict.lm()
nd <- d[1:6, ]
nd$x[2] <- NA
expect_equal(predict(design_lm(y ~ x * f, d), nd), predict(lm(y ~ x * f, d), nd))

# predictions from a rank-deficient fit: non-estimable rows are flagged as
# by predict.lm()
nd <- d[1:8, ]
nd$f2 <- factor(c("a", "b", "c", "d", "b", "a", "d", "c"), levels = levels(d$f))
a <- lm(y ~ x + f + f2, d)
b <- design_lm(y ~ x + f + f2, d)
wa <- capture(function() predict(a, nd))
wb <- capture(function() predict(b, nd))
expect_equal(wb$value, wa$value)
expect_equal(wb$warnings, wa$warnings)
expect_true(length(attr(wb$value, "non-estim")) > 0)
for (rd in c("simple", "non-estim", "NA", "NAwarn")) {
    wa <- capture(function() predict(a, nd, rankdeficient = rd))
    wb <- capture(function() predict(b, nd, rankdeficient = rd))
    expect_equal(wb$value, wa$value)
    expect_equal(wb$warnings, wa$warnings)
}
expect_equal(predict(b, d[1:8, ], rankdeficient = "non-estim"),
             predict(a, d[1:8, ], rankdeficient = "non-estim"))

# model = FALSE keeps no data
b <- design_lm(y ~ x * f, d, model = FALSE)
expect_null(b$model)
expect_error(residuals(b), "model = TRUE")
expect_error(predict(b), "model = TRUE")
expect_equal(predict(b, d[1:5, ]), predict(lm(y ~ x * f, d), d[1:5, ]))
expect_equal(nrow(model.frame(design_lm(y ~ x * f, d))), n - 1L)

# print and summary print without error
expect_stdout(print(design_lm(y ~ f * g + f2, d)), "singular|NA")
expect_stdout(print(summary(design_lm(y ~ f * g + f2, d))), "not defined because of singularities")
expect_stdout(print(design_lm(y ~ f * g, d)), "design patterns")

# errors as in lm()
expect_error(design_lm(y ~ x, d, weights = -w), "negative weights")
d5 <- d
d5$y[1] <- Inf
expect_error(design_lm(y ~ x, d5), "NA/NaN/Inf in 'y'")
expect_error(design_lm(cbind(y, x) ~ f, d), "Multivariate")
expect_error(design_lm(y ~ x, d, compress = "yes"), "compress")
expect_error(design_lm(y ~ x, d[0, ]), "0 \\(non-NA\\) cases")

# design_lm() on a design matrix --------------------------------------------------------------

dc <- d[complete.cases(d), ]
D <- design_matrix(y ~ x * f + f2, dc)
b <- design_lm(D, y = dc$y, weights = dc$w, offset = dc$o)
a <- lm(y ~ x * f + f2, dc, weights = w, offset = o)
expect_true(same_fit(a, b))
expect_equal(residuals(b), residuals(a))
expect_equal(predict(b, dc[1:5, ]), predict(a, dc[1:5, ]) - dc$o[1:5])
expect_equal(summary(b)$r.squared, summary(a)$r.squared)
D0 <- design_matrix(y ~ 0 + f, dc)
expect_equal(summary(design_lm(D0, y = dc$y))$r.squared, summary(lm(y ~ 0 + f, dc))$r.squared)
expect_error(design_lm(D, y = dc$y[-1]))
expect_error(design_lm(D, dc, y = dc$y), "do not apply")
expect_error(design_lm(y ~ f, dc, y = dc$y), "applies only")

# Fuzz: random formulas, weights, offsets, missing values -------------------

set.seed(20261125)
compared <- 0L
for (i in 1:80) {
    nn <- sample(c(30, 80, 200), 1)
    dd <- fuzz_data(nn, special = FALSE, na = runif(1) < 0.3)
    dd$w <- if (runif(1) < 0.4) runif(nn) else 1
    if (runif(1) < 0.2) dd$w[sample(nn, 2)] <- 0
    dd$o <- rnorm(nn)
    fo <- fuzz_formula()
    ct <- fuzz_contrasts(dd)
    args <- list(fo, dd, contrasts = ct)
    if (runif(1) < 0.4) args$weights <- quote(w)
    if (runif(1) < 0.3) args$offset <- quote(o)
    a <- capture(function() do.call(lm, args))
    if (inherits(a$value, "captured_error")) next
    a <- a$value
    ## (X'WX)^-1 is only determined to about kappa^2 * eps, whatever the
    ## method (lm()'s own vcov too): on badly conditioned designs, the
    ## tolerance grows with the condition number
    kap <- if (a$rank > 0) kappa(qr.R(a$qr)[seq_len(a$rank), seq_len(a$rank), drop = FALSE], exact = TRUE) else 1
    tol <- max(1e-6, 1000 * .Machine$double.eps * kap^2)
    for (cmp in list("auto", FALSE)) {
        b <- capture(function() do.call(design_lm, c(args, list(compress = cmp))))
        ok <- !inherits(b$value, "captured_error") && same_fit(a, b$value, tol = tol)
        if (!ok) print(list(i = i, formula = fo, compress = cmp, args = names(args)))
        expect_true(ok)
    }
    b <- b$value
    if (!is.null(b$spec) && is.null(args$offset)) {
        nd <- fuzz_data(10, special = FALSE)
        for (lv in c("f1", "f2", "o1", "ch", "lg")) {
            ## new data with levels seen in training
            if (is.factor(dd[[lv]])) nd[[lv]] <- factor(sample(levels(droplevels(dd[[lv]])), 10, TRUE), levels = levels(dd[[lv]]), ordered = is.ordered(dd[[lv]]))
            else nd[[lv]] <- sample(stats::na.omit(unique(dd[[lv]])), 10, TRUE)
        }
        pa <- capture(function() suppressWarnings(predict(a, nd)))
        if (!inherits(pa$value, "captured_error")) {
            pb <- suppressWarnings(predict(b, nd))
            ## aliased coefficients count as zero in both
            expect_equal(pb, pa$value, tolerance = tol)
        }
    }
    compared <- compared + 1L
}
expect_true(compared > 50)

# Fuzz: factor-only designs, compressed against uncompressed and lm() --------

set.seed(7)
for (i in 1:25) {
    nn <- sample(c(50, 500, 3000), 1)
    k <- sample(2:6, 3, TRUE)
    dd <- data.frame(
        a = factor(sample(seq_len(k[1]), nn, TRUE)),
        b = factor(sample(seq_len(k[2]), nn, TRUE)),
        c = factor(sample(seq_len(k[3]), nn, TRUE)),
        i = sample(1:3, nn, TRUE),
        w = runif(nn), o = rnorm(nn))
    dd$y <- rnorm(nn, mean = as.integer(dd$a))
    dd$w[sample(nn, 3)] <- 0
    if (runif(1) < 0.3) dd$a[sample(nn, 2)] <- NA
    fo <- sample(c(y ~ a * b, y ~ a + b + c, y ~ a:b + c, y ~ 0 + a * b * c, y ~ a * i + b,
                   y ~ a + b + a:b:c), 1)[[1]]
    wt <- runif(1) < 0.5
    of <- runif(1) < 0.3
    args <- list(fo, dd)
    if (wt) args$weights <- quote(w)
    if (of) args$offset <- quote(o)
    a <- do.call(lm, args)
    bc <- do.call(design_lm, c(args, list(compress = TRUE)))
    bu <- do.call(design_lm, c(args, list(compress = FALSE)))
    expect_true(bc$compressed)
    expect_true(same_fit(a, bc))
    expect_true(same_fit(a, bu))
    expect_true(same_fit(bc, bu, tol = 1e-9))
    expect_equal(logLik(bc), logLik(a))
    expect_equal(summary(bc)$r.squared, suppressWarnings(summary(a))$r.squared)
    expect_equal(residuals(bc), residuals(a))
    nd <- dd[sample(nn, 20), ]
    expect_equal(predict(bc, nd, se.fit = TRUE)$se.fit, predict(a, nd, se.fit = TRUE)$se.fit)
}

# compression of rows with a numeric variable that repeats
dd <- data.frame(x = rep(1:5, 400), f = gl(4, 500), y = rnorm(2000))
b <- design_lm(y ~ x * f, dd)
expect_true(b$compressed)
expect_true(same_fit(lm(y ~ x * f, dd), b))

# compress = "auto" follows the option xtx.dedup
op <- options(xtx.dedup = FALSE)
expect_false(design_lm(y ~ x * f, transform(dd))$compressed)
options(op)

# Review findings (Roborev job 1273) ---------------------------------------------

# new data whose variables have another type than in training: the error of
# predict.lm(), from the plan (general and fast paths) and from predict()
dt <- data.frame(y = rnorm(20), x = rnorm(20), f = gl(2, 10),
                 ch = rep(c("u", "v"), 10), lg = rep(c(TRUE, FALSE), 10),
                 stringsAsFactors = FALSE)
nd_bad <- data.frame(x = factor(c("a", "b")), f = factor(c("1", "2")),
                     ch = c("u", "v"), lg = c(TRUE, FALSE))
msg <- tryCatch(predict(lm(y ~ x + f, dt), nd_bad), error = conditionMessage)
expect_true(grepl("fitted with type", msg))
spec <- design_matrix(design_matrix(y ~ x + f, dt), newdata = NULL)
expect_error(design_matrix(spec, nd_bad), msg, fixed = TRUE)
expect_error(design_matrix(design_matrix(y ~ x + f, dt), nd_bad), msg, fixed = TRUE)
expect_error(design_matrix(xtx:::design_plan(y ~ x + f, dt), nd_bad), msg, fixed = TRUE)
fit <- design_lm(y ~ x + f, dt)
expect_error(predict(fit, nd_bad), msg, fixed = TRUE)
# the plan's fast path is learned from good data, then gets bad data
p_lazy <- design_matrix(design_matrix(y ~ x + f, dt), newdata = NULL)
invisible(design_matrix(p_lazy, dt[1:3, ]))
expect_error(design_matrix(p_lazy, nd_bad), msg, fixed = TRUE)
# two variables of the wrong type
nd_bad2 <- transform(nd_bad, f = c(1, 2))
msg2 <- tryCatch(predict(lm(y ~ x + f, dt), nd_bad2), error = conditionMessage)
expect_error(predict(fit, nd_bad2), msg2, fixed = TRUE)
# substitutes that predict.lm() accepts: a factor for a character variable,
# integers for doubles, logicals for logicals
fit <- design_lm(y ~ x + ch + lg, dt)
nd_ok <- data.frame(x = 1:2, ch = factor(c("v", "u")), lg = c(FALSE, TRUE))
expect_equal(predict(fit, nd_ok), predict(lm(y ~ x + ch + lg, dt), nd_ok))

# offsets of new data with rows omitted by na.action: the offsets of the
# rows kept, in order
do <- data.frame(y = rnorm(30), x = rnorm(30), o = rnorm(30))
nd_na <- data.frame(x = c(NA, 1, 2, 3), o = c(10, 20, 30, 40))
for (fo in list(y ~ x + offset(o), y ~ x)) {
    fit <- if (length(attr(terms(fo), "offset"))) design_lm(fo, do) else design_lm(fo, do, offset = o)
    p <- predict(fit, nd_na, na.action = na.omit)
    expect_equal(names(p), c("2", "3", "4"))
    expect_equal(unname(p), coef(fit)[[1]] + coef(fit)[[2]] * c(1, 2, 3) + c(20, 30, 40))
    expect_equal(predict(fit, nd_na), c(NA, p), check.attributes = FALSE)
}

# integer weights whose sum in a pattern exceeds the integer range
dw <- data.frame(f = gl(3, 3), y = rnorm(9), w = rep(1000000000L, 9))
b_c <- design_lm(y ~ f, dw, weights = w, compress = TRUE)
b_u <- design_lm(y ~ f, dw, weights = w, compress = FALSE)
expect_true(b_c$compressed)
expect_true(same_fit(b_c, b_u, tol = 1e-10))
expect_true(same_fit(lm(y ~ f, dw, weights = w), b_c))

# a zero-weight row with a huge response does not enter the fit
dz <- data.frame(f = gl(3, 4), y = rnorm(12), w = 1)
dz$y[1] <- 1e200
dz$w[1] <- 0
b_c <- design_lm(y ~ f, dz, weights = w, compress = TRUE)
expect_true(b_c$compressed)
expect_true(is.finite(sigma(b_c)))
## (vcov(lm()) is NaN here: summary.lm() multiplies the row's infinite
## squared residual by its zero weight; lm() without the row is the reference)
expect_true(same_fit(lm(y ~ f, dz[-1, ]), b_c))
expect_true(same_fit(design_lm(y ~ f, dz, weights = w, compress = FALSE), b_c))

# a saturated fit: prediction standard errors are NaN, as in predict.lm()
d4 <- data.frame(y = c(1, 2, 3), f = factor(c("a", "b", "c")), f2 = factor(c("a", "b", "c")))
for (fo in list(y ~ f, y ~ f + f2)) {
    pa <- suppressWarnings(predict(lm(fo, d4), d4, se.fit = TRUE))
    pb <- suppressWarnings(predict(design_lm(fo, d4), d4, se.fit = TRUE))
    expect_true(all(is.nan(pb$se.fit)))
    expect_equal(unname(pb$se.fit), unname(pa$se.fit))
}

# residuals: (y - offset) - X b, which keeps the digits when y and the
# offset are large and close
dr <- data.frame(x = rnorm(20), o = 1e16 + 2 * (1:20))
dr$y <- dr$o + rnorm(20)
a <- lm(y ~ x, dr, offset = o)
b <- design_lm(y ~ x, dr, offset = o)
expect_equal(residuals(b), residuals(a))
expect_equal(sum(residuals(b)^2), deviance(b))

# predictions at the training rows are not checked for estimability, as in
# predict.lm(): a zero-weight row outside the estimable space keeps its
# fitted value with rankdeficient = "NA"
dz <- data.frame(y = rnorm(12), f = factor(rep(c("a", "b", "c"), each = 4)))
dz$g <- factor(ifelse(dz$f == "a", "u", "v"))
dz$g[12] <- "u"
w <- c(rep(1, 11), 0)
a <- lm(y ~ f + g, dz, weights = w)
b <- design_lm(y ~ f + g, dz, weights = w)
expect_true(anyNA(coef(b)))
expect_equal(predict(b, rankdeficient = "NA"), fitted(b))
expect_equal(predict(b, rankdeficient = "NA"), predict(a, rankdeficient = "NA"))

# compressed sums that overflow fall back to the fit on the rows
dov <- data.frame(x = rep(2^-300, 4), y = rep(2^450, 4))
fo <- design_lm(y ~ x - 1, dov, weights = rep(2^600, 4), compress = TRUE)
expect_false(fo$compressed)
expect_true(is.finite(coef(fo)))
expect_equal(unname(coef(fo)), 2^750)

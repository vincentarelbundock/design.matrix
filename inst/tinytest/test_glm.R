source("helpers.R")
source("glm_cases.R", local = TRUE)
library(design.matrix)

# design_glm() must agree with glm(): coefficients (with the NA pattern of
# aliased columns), vcov, deviances, AIC, degrees of freedom, dispersion and
# the number of iterations, compressed or not.
d <- glm_data()
d$f <- factor(d$f)
d$g <- factor(d$g)

for (cs in glm_cases) {
    b <- suppressWarnings(glm_ref(cs, d))
    for (cp in list("auto", FALSE)) {
        a <- suppressWarnings(eval(bquote(design_glm(.(cs[[2]]), family = .(cs[[3]]), data = d,
                     weights = .(if (length(cs) >= 4 && !is.null(cs[[4]])) as.name(cs[[4]])),
                     offset = .(if (length(cs) >= 5 && !is.null(cs[[5]])) as.name(cs[[5]])),
                     compress = .(cp)))))
        glm_same(a, b, paste(cs[[1]], "compress", cp))
    }
}

# compression is used where rows share patterns
fit <- design_glm(yb ~ xd * f + g, binomial, d)
expect_true(fit$compressed)
expect_equal(fit$patterns, 48L)
expect_false(design_glm(yb ~ xd * f + g, binomial, d, compress = FALSE)$compressed)
# a continuous offset splits patterns: not worth compressing
expect_false(design_glm(yp ~ xd + g + offset(log(e)), poisson, d)$compressed)
# an offset with few values: compressed by pattern and offset
fo <- design_glm(yp ~ g + offset(log(xd)), poisson, d)
expect_true(fo$compressed)
expect_equal(fo$patterns, 12L)
glm_same(fo, glm(yp ~ g + offset(log(xd)), poisson, d), "offset with few values")

# warnings: separation and non-convergence, as glm()
w1 <- glm_warnings(glm(sep ~ x, binomial, d))$warnings
expect_identical(glm_warnings(design_glm(sep ~ x, binomial, d))$warnings, w1)
w2 <- glm_warnings(glm(yb ~ x + f, binomial, d, control = list(maxit = 2)))$warnings
r2 <- glm_warnings(design_glm(yb ~ x + f, binomial, d, control = list(maxit = 2)))
expect_identical(r2$warnings, w2)
expect_false(r2$value$converged)
# one iteration only: the dispersion uses the working weights of the rows
b1 <- suppressWarnings(glm(yg ~ xd + f, Gamma, d, control = list(maxit = 1)))
a1 <- suppressWarnings(design_glm(yg ~ xd + f, Gamma, d, control = list(maxit = 1)))
glm_same(a1, b1, "maxit 1")
# initialize's checks and warnings
expect_error(design_glm(yg - 3 ~ x, Gamma, d), "non-positive values")
expect_error(design_glm(yp - 1 ~ x, poisson, d), "negative values")
expect_warning(design_glm(I(prop * 0.9) ~ xd, binomial, d, weights = five), "non-integer")
# zero weights: the dispersion warning of summary.glm()
expect_warning(summary(design_glm(yg ~ x, Gamma, d, weights = w)), "zero weight")

# NA rows: na.omit and na.exclude
fx <- design_glm(yp ~ x + g, poisson, d, na.action = na.exclude)
bx <- glm(yp ~ x + g, poisson, d, na.action = na.exclude)
expect_equal(residuals(fx), residuals(bx))
expect_equal(fitted(fx), fitted(bx))
expect_error(design_glm(yp ~ x + g, poisson, d, na.action = na.fail), "missing values")

# methods
b <- glm(yp ~ x + g + offset(log(e)), poisson, d)
a <- design_glm(yp ~ x + g + offset(log(e)), poisson, d)
for (ty in c("deviance", "pearson", "working", "response"))
    expect_equal(residuals(a, ty), residuals(b, ty), info = ty)
expect_equal(fitted(a), fitted(b))
expect_equal(predict(a), predict(b))
expect_equal(predict(a, d[1:20, ], type = "response"), predict(b, d[1:20, ], type = "response"))
expect_equal(logLik(a), logLik(b))
expect_equal(BIC(a), BIC(b))
expect_identical(family(a)$family, "poisson")
expect_equal(summary(a)$coefficients, summary(b)$coefficients)
expect_true(any(grepl("Fisher Scoring", capture.output(print(summary(a))))))
expect_true(any(grepl("Residual Deviance", capture.output(print(a)))))

# a factor response and a design matrix with y
d$ff <- factor(d$fy)
expect_equal(coef(design_glm(ff ~ x, binomial, d)), coef(glm(ff ~ x, binomial, d)))
dc <- d[complete.cases(d[c("x", "g", "yb")]), ]
D <- design_matrix(~ x + g, dc)
expect_equal(coef(design_glm(D, binomial, y = dc$yb)), coef(glm(yb ~ x + g, binomial, dc)))
expect_equal(coef(design_glm(D, binomial, y = cbind(dc$s, dc$fl))),
             coef(glm(cbind(s, fl) ~ x + g, binomial, dc)))

# start values
s0 <- coef(glm(yp ~ xd + g, poisson, d)) + 0.1
expect_equal(design_glm(yp ~ xd + g, poisson, d, start = s0)$iter, glm(yp ~ xd + g, poisson, d, start = s0)$iter)

# what design_glm() does not support
expect_error(design_glm(yp ~ x, poisson("identity") |> (\(f) { f$family <- "other"; f })(), d), "supports the families")
expect_error(design_glm(yp ~ x, quasi(variance = "mu", link = power(1/3)), d), "supports the families")

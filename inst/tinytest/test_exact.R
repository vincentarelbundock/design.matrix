source("helpers.R")
library(xtx)

old_threads <- design_threads()

# Fixed cases ------------------------------------------------------------
set.seed(1)
d <- fuzz_data(200)
fixed <- list(y ~ x1 * f1 + o1, ~ 0 + f1:x1 + poly(z, 3), y ~ (x1 + f1 + ch)^2 + lg,
              ~ log(z) + f1 %in% o1, y ~ splines::ns(z, 3) * f2, ~ 1, ~ 0 + f1,
              y ~ x1 + y, y ~ m * f1, ~ 0 + ch:lg + i1)
for (fo in fixed) {
    expect_true(same_outcome(capture(function() as.matrix(design_matrix(fo, d))),
                             capture(function() stats::model.matrix(fo, d))),
                info = deparse(fo))
}

# Signed zeros, NaN payloads, Inf: bitwise, not just ==
d0 <- data.frame(x = c(-1.5, 2, NA, Inf, -0, NaN), f = factor(c("a", "b", "a", "b", "a", "b")))
mf0 <- stats::model.frame(~ x:f, d0, na.action = stats::na.pass)
expect_true(strict_identical(as.matrix(design_matrix(~ 0 + x:f, mf0)), stats::model.matrix(~ 0 + x:f, mf0)))
X0 <- as.matrix(design_matrix(~ 0 + x:f, mf0))
expect_true(any(X0 == 0 & 1 / X0 < 0, na.rm = TRUE)) # negative zeros are there

# contrasts.arg, xlev, missing data argument
expect_true(strict_identical(
    as.matrix(design_matrix(y ~ f1 * o1, d, contrasts.arg = list(f1 = "contr.sum", o1 = stats::contr.treatment(3)))),
    stats::model.matrix(y ~ f1 * o1, d, contrasts.arg = list(f1 = "contr.sum", o1 = stats::contr.treatment(3)))))
y <- d$y; x1 <- d$x1; f1 <- d$f1
expect_true(strict_identical(as.matrix(design_matrix(y ~ x1 + f1)), stats::model.matrix(y ~ x1 + f1)))

# terms object and an existing model frame
mf <- stats::model.frame(y ~ x1 * f2 + ch, d)
expect_true(strict_identical(as.matrix(design_matrix(attr(mf, "terms"), mf)), stats::model.matrix(attr(mf, "terms"), mf)))

# Errors and warnings match stats
d1 <- d
d1$f1 <- factor(rep("a", nrow(d1)))
expect_true(same_outcome(capture(function() as.matrix(design_matrix(y ~ f1, d1))),
                         capture(function() stats::model.matrix(y ~ f1, d1))))
dc <- data.frame(y = 1:3, z = complex(real = 1:3, imaginary = 1))
expect_true(same_outcome(capture(function() as.matrix(design_matrix(y ~ z, dc))),
                         capture(function() stats::model.matrix(y ~ z, dc))))

# Non-formula objects fall back to stats
fit <- stats::lm(y ~ i1 + f1, data = d)
expect_true(strict_identical(as.matrix(design_matrix(fit)), stats::model.matrix(fit)))
expect_identical(design_matrix(fit)$path, "dense")

# Threads: same bits for several row chunks
big <- fuzz_data(20000)
for (fo in list(y ~ x1 * f1 + o1 + poly(z, 2), ~ 0 + f2:ch + m:lg)) {
    ref <- stats::model.matrix(fo, big)
    for (th in c(1L, 3L)) {
        design_threads(th)
        expect_true(strict_identical(as.matrix(design_matrix(fo, big)), ref), info = paste(deparse(fo), th))
    }
}
design_threads(old_threads)

# Randomized comparison --------------------------------------------------
set.seed(20261125)
n_cases <- 300
ok <- logical(n_cases)
for (k in seq_len(n_cases)) {
    n <- sample(c(1, 2, 5, 30, 200), 1)
    dk <- fuzz_data(n, special = TRUE, na = stats::runif(1) < 0.3)
    fo <- fuzz_formula(response = stats::runif(1) < 0.7)
    ca <- fuzz_contrasts(dk)
    design_threads(sample(1:3, 1))
    a <- capture(function() as.matrix(design_matrix(fo, dk, contrasts.arg = ca)))
    b <- capture(function() stats::model.matrix(fo, dk, contrasts.arg = ca))
    ok[k] <- same_outcome(a, b)
    if (!ok[k]) {
        expect_true(ok[k], info = sprintf("case %d: %s", k, deparse(fo)))
    }
}
design_threads(old_threads)
expect_true(all(ok), info = sprintf("%d of %d random cases differ", sum(!ok), n_cases))

# Zero rows
d0 <- fuzz_data(5)[0, ]
for (fo in list(y ~ x1 * f1, ~ 0 + f2 + ch)) {
    expect_true(same_outcome(capture(function() as.matrix(design_matrix(fo, d0))),
                             capture(function() stats::model.matrix(fo, d0))), info = deparse(fo))
}

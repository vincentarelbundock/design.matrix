# Grouping identical rows never changes results.
source("helpers.R")
library(design.matrix)
old <- options(design.matrix.dedup = "auto")

with_mode <- function(mode, f) {
    op <- options(design.matrix.dedup = mode)
    on.exit(options(op))
    f()
}

set.seed(31)
base <- fuzz_data(40, special = TRUE, na = TRUE)
d <- base[sample(nrow(base), 4000, replace = TRUE), ]  # many duplicates
d$x1[1:3] <- c(0, -0, NaN)                              # distinct bits, same value
p <- design_matrix(y ~ x1 * f1 + poly(z, 2) + ch + lg + o1, base[stats::complete.cases(base), ],
             contrasts.arg = list(f1 = "contr.sum"))

for (k in 1:2) { # fresh objects, then cached groups
    D_on <- design_matrix(p, d)
    D_off <- design_matrix(p, d)
    P <- ncol(D_on)
    b <- rnorm(P); B <- matrix(rnorm(2 * P), P); w <- rnorm(nrow(d))
    w[5] <- NA
    V <- crossprod(matrix(rnorm(P * P), P))
    g <- sample(c("a", "b", NA), nrow(d), TRUE)
    for (rep in 1:2) {
        expect_true(strict_identical(with_mode(TRUE, function() quadratic_diagonal(D_on, V)),
                                     with_mode(FALSE, function() quadratic_diagonal(D_off, V))))
        expect_true(strict_identical(with_mode(TRUE, function() D_on %*% b),
                                     with_mode(FALSE, function() D_off %*% b)))
        expect_true(strict_identical(with_mode(TRUE, function() D_on %*% B),
                                     with_mode(FALSE, function() D_off %*% B)))
        a <- with_mode(TRUE, function() crossprod(D_on, w))
        e <- with_mode(FALSE, function() crossprod(D_off, w))
        expect_identical(is.na(a), is.na(e))
        expect_equal(a, e, tolerance = 1e-10)
        a <- with_mode(TRUE, function() suppressWarnings(colSums(w * D_on, group = g)))
        e <- with_mode(FALSE, function() suppressWarnings(colSums(w * D_off, group = g)))
        expect_identical(is.na(a), is.na(e))
        expect_equal(a, e, tolerance = 1e-10)
    }
    expect_true(is.list(D_on$cache$groups))
    expect_true(length(D_on$cache$groups$first) < nrow(d) / 2)
}

# Jacobians: scales are part of each row's pattern
X0 <- design_matrix(p, d)
X1 <- design_matrix(p, transform(d, x1 = x1 + 1))
s0 <- stats::plogis(rnorm(40))[match(rownames(d), rownames(d))]
s0 <- rep_len(stats::runif(13), nrow(d))  # breaks some duplicate patterns
J_on <- 0.5 * X1 - s0 * X0
J_off <- 0.5 * X1 - s0 * X0
V <- crossprod(matrix(rnorm(ncol(X0)^2), ncol(X0)))
expect_true(strict_identical(with_mode(TRUE, function() quadratic_diagonal(J_on, V)),
                             with_mode(FALSE, function() quadratic_diagonal(J_off, V))))
expect_equal(with_mode(TRUE, function() crossprod(J_on, rep(1, nrow(d)))),
             with_mode(FALSE, function() crossprod(J_off, rep(1, nrow(d)))), tolerance = 1e-10)

# auto mode: no groups when every row is distinct (and the failure is cached)
du <- data.frame(x = rnorm(5000), f = gl(5, 1000))
Du <- design_matrix(~ x * f, du)
q <- quadratic_diagonal(Du, diag(ncol(Du)))
expect_false(is.list(Du$cache$groups)) # not attempted (cheap rows) or attempted and cached as FALSE
expect_true(strict_identical(q, with_mode(FALSE, function() quadratic_diagonal(design_matrix(~ x * f, du), diag(ncol(Du))))))

options(old)

# a costly kernel on distinct rows: attempted, given up, failure cached
op <- options(design.matrix.dedup = "auto")
Dd <- design_matrix(~ poly(x, 3) * f, du, contrasts.arg = list(f = "contr.sum"))
invisible(quadratic_diagonal(Dd, diag(ncol(Dd))))
expect_identical(Dd$cache$groups, FALSE)
options(op)

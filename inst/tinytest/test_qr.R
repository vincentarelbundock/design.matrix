source("helpers.R")
library(design.matrix)

# Least squares through qr() on a design matrix agrees with qr() on the
# dense matrix, up to rounding: same rank and pivot, coefficients (NA for
# aliased columns), fitted values, residuals, and R, Q and t(Q) y up to sign.
close_to <- function(a, b, tol = 1e-8) {
    a <- unname(as.matrix(a)); b <- unname(as.matrix(b))
    identical(dim(a), dim(b)) && identical(is.na(a), is.na(b)) &&
        isTRUE(all(abs(a[!is.na(a)] - b[!is.na(b)]) <= tol * pmax(1, abs(b[!is.na(b)]))))
}

check_qr <- function(D, y, info) {
    X <- as.matrix(D)
    q <- qr(D)
    q0 <- qr(X)
    expect_inherits(q, "design_qr", info = info)
    expect_equal(q$rank, q0$rank, info = paste(info, "rank"))
    expect_equal(q$pivot, q0$pivot, info = paste(info, "pivot"))
    expect_true(close_to(qr.coef(q, y), qr.coef(q0, y)), info = paste(info, "coef"))
    expect_identical(names(qr.coef(q, y)), names(qr.coef(q0, y)), info = paste(info, "coef names"))
    expect_true(close_to(qr.fitted(q, y), qr.fitted(q0, y)), info = paste(info, "fitted"))
    expect_true(close_to(qr.resid(q, y), qr.resid(q0, y)), info = paste(info, "resid"))
    k <- q$rank
    i <- seq_len(k)
    expect_true(close_to(abs(qr.R(q)[i, i]), abs(qr.R(q0)[i, i])), info = paste(info, "R"))
    expect_true(close_to(abs(qr.R(q)[i, -i]), abs(qr.R(q0)[i, -i])), info = paste(info, "R12"))
    expect_true(close_to(abs(qr.qty(q, y)), abs(qr.qty(q0, y)[i])), info = paste(info, "qty"))
    Q <- qr.Q(q)
    expect_true(close_to(abs(Q), abs(qr.Q(q0)[, i, drop = FALSE])), info = paste(info, "Q"))
    expect_true(close_to(qr.qy(q, seq_len(k)), Q %*% seq_len(k)), info = paste(info, "qy"))
    Y <- cbind(a = y, b = rev(y))
    expect_true(close_to(qr.coef(q, Y), qr.coef(q0, Y)), info = paste(info, "coef matrix"))
    expect_identical(dimnames(qr.coef(q, Y)), dimnames(qr.coef(q0, Y)), info = paste(info, "coef dimnames"))
}

set.seed(11)
d <- fuzz_data(300, special = FALSE)
forms <- list(y ~ x1 * f1 + f2, y ~ poly(z, 3) * ch + m:lg, y ~ splines::ns(z, 3):f2 + i1,
              y ~ 1, y ~ f1 * f2 * o1, y ~ x1 + x2 + I(x1 + x2) + f1, y ~ 0 + f1 + f2)
for (fo in forms) {
    D <- suppressWarnings(design_matrix(fo, d, contrasts.arg = list(f2 = "contr.sum", o1 = "contr.helmert")))
    check_qr(D, d$y + d$x1, deparse(fo))
}

# random formulas, contrasts and data
for (r in 1:40) {
    d <- fuzz_data(sample(c(5, 30, 200), 1), special = FALSE)
    fo <- fuzz_formula()
    D <- tryCatch(suppressWarnings(design_matrix(fo, d, contrasts.arg = fuzz_contrasts(d))),
                  error = function(e) NULL)
    if (is.null(D) || ncol(D) == 0L || !all(is.finite(as.matrix(D)))) next
    check_qr(D, d$x1 * 2 + rnorm(nrow(d)), paste(r, deparse(fo)))
}

# matches lm()
d <- fuzz_data(500, special = FALSE)
fit <- lm(y ~ x1 * f1 + poly(z, 2), d)
q <- qr(design_matrix(fit))
expect_equal(qr.coef(q, d$y), coef(fit))
expect_equal(qr.resid(q, d$y), unname(residuals(fit)))
expect_equal(solve(q, d$y), coef(fit))

# weighted least squares through scaled rows
w <- runif(500)
fitw <- lm(y ~ x1 * f1, d, weights = w)
Dw <- sqrt(w) * design_matrix(y ~ x1 * f1, d)
expect_equal(qr.coef(qr(Dw), sqrt(w) * d$y), coef(fitw))

# badly conditioned designs: refinement recovers the accuracy of a
# Householder QR, where the plain normal equations do not
set.seed(3)
n <- 5000
dc <- data.frame(x = runif(n) + 100, f = gl(5, n / 5))
dc$y <- 1 + dc$x + rnorm(n)
D <- design_matrix(y ~ f + x + I(x^2), dc)
X <- as.matrix(D)
expect_true(kappa(X, exact = TRUE) > 1e9)
expect_equal(qr(D)$rank, qr(X)$rank)
expect_equal(qr.coef(qr(D), dc$y), qr.coef(qr(X), dc$y), tolerance = 1e-8)

# rank detection follows `tol`
expect_equal(qr(D, tol = 1e-3)$rank, qr(X, tol = 1e-3)$rank)

# empty and degenerate designs
d1 <- data.frame(y = rnorm(5))
q <- qr(design_matrix(y ~ 0, d1))
expect_equal(q$rank, 0L)
expect_identical(qr.coef(q, d1$y), setNames(numeric(0), character(0)))
expect_equal(qr.fitted(q, d1$y), rep(0, 5))
expect_equal(qr.resid(q, d1$y), d1$y)
expect_equal(qr.coef(qr(design_matrix(y ~ 1, d1)), d1$y), c("(Intercept)" = mean(d1$y)))

# errors
d <- data.frame(y = 1:4, x = c(1, Inf, 2, 3))
expect_error(qr(design_matrix(y ~ x, d)), "NA/NaN/Inf")
D <- design_matrix(y ~ x, data.frame(y = rnorm(6), x = rnorm(6)))
q <- qr(D)
expect_error(qr.coef(q, c(1, NA, 1, 1, 1, 1)), "finite")
expect_error(qr.coef(q, 1:3), "6 rows")
expect_error(qr.fitted(q, rnorm(6), k = 1), "k = qr\\$rank")
expect_error(qr.Q(q, complete = TRUE), "complete")
expect_error(qr.qy(q, rnorm(6)), "2 rows")
expect_error(solve(q), "only square")
expect_error(qr(D, LAPACK = TRUE), "LAPACK")
dd <- data.frame(y = rnorm(6), x = rnorm(6)); dd$h <- 2 * dd$x
expect_error(solve(qr(design_matrix(y ~ x + h, dd)), dd$y), "singular")

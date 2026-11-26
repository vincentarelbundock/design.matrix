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

# Roborev 1271 regressions ----------------------------------------------

same_qr <- function(D, y, tol = 1e-07) {
    X <- as.matrix(D)
    q <- qr(D, tol = tol); q0 <- qr(X, tol = tol)
    identical(q$rank, q0$rank) && identical(q$pivot, q0$pivot) &&
        isTRUE(all.equal(qr.coef(q, y), qr.coef(q0, y)))
}

# 1. After a suspect column is accepted, later suspects are judged again:
#    with orthonormal a and z, a + 2e-5 z is a combination of a and a + 1e-5 z
set.seed(1)
Z <- qr.Q(qr(matrix(rnorm(300), 100)))
a <- Z[, 1]; z <- Z[, 2]
d1 <- data.frame(y = rnorm(100), x1 = a, x2 = a + 1e-5 * z, x3 = a + 2e-5 * z)
D1 <- design_matrix(y ~ 0 + x1 + x2 + x3, d1)
expect_equal(qr(D1)$rank, 2L)
expect_true(same_qr(D1, d1$y))

# 3. A column above the suspect cutoff but below `tol` is aliased
d3 <- data.frame(y = rnorm(100), x1 = a, x2 = a + 5e-4 * z)
D3 <- design_matrix(y ~ 0 + x1 + x2, d3)
expect_equal(qr(D3, tol = 1e-3)$rank, 1L)
expect_true(same_qr(D3, d3$y, tol = 1e-3))
expect_equal(qr(D3)$rank, 2L)

# 6. Columns whose squares overflow or underflow keep their rank (dense
#    fallback), and zero columns are still aliased
for (s in c(1e200, 1e-200, 1)) {
    d6 <- data.frame(y = rnorm(20), x = s * rnorm(20), w = rnorm(20), z0 = 0)
    D6 <- design_matrix(y ~ x + w + z0, d6)
    expect_equal(qr(D6)$rank, 3L, info = s)
    expect_true(same_qr(D6, d6$y), info = s)
}

# 8. solve() fills coefficients outside the basis with zero, as solve.qr()
set.seed(2)
dw <- data.frame(y = rnorm(3), x1 = rnorm(3), x2 = rnorm(3), x3 = rnorm(3), x4 = rnorm(3))
Dw <- design_matrix(y ~ x1 + x2 + x3 + x4, dw)
sol <- solve(qr(Dw), dw$y)
expect_false(anyNA(sol))
expect_equal(drop(Dw %*% sol), dw$y)
expect_equal(sol, solve.qr(qr(as.matrix(Dw)), dw$y), check.attributes = FALSE)

# 9. The methods dispatch through the package namespace
q9 <- qr(D1)
expect_equal(design.matrix::qr.coef(q9, d1$y), qr.coef(q9, d1$y))

# Roborev 1272 regressions -------------------------------------------------

# 4. The dense fallback computes Q products and fitted values from the
#    stored QR, without overflowing intermediates
D4 <- design_matrix(~ 0 + x, data.frame(x = 1e-200))
q4 <- qr(D4)
expect_false(is.null(q4$dense))
expect_true(is.finite(qr.qy(q4, 1e200)))
expect_equal(abs(qr.qy(q4, 1e200)), abs(drop(base::qr.qy(qr(as.matrix(D4)), 1e200))))
expect_equal(qr.fitted(q4, 3), 3)
expect_equal(qr.resid(q4, 3), 0)
expect_equal(abs(qr.Q(q4)), abs(base::qr.Q(qr(as.matrix(D4)))))

# 5. A residual ratio equal to tol keeps the column, as in qr()
D5 <- design_matrix(~ 0 + x, data.frame(x = c(1, 0, 0)))
expect_equal(qr(D5, tol = 1)$rank, qr(as.matrix(D5), tol = 1)$rank)
expect_equal(qr(D5, tol = 1)$rank, 1L)

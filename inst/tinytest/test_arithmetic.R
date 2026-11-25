source("helpers.R")
library(design.matrix)
design_row_groups <- design.matrix:::design_row_groups

set.seed(11)
d <- data.frame(x = rnorm(30), f = gl(3, 10))
p <- design_matrix(~ x * f, d)
X0 <- design_matrix(p, d)
X1 <- design_matrix(p, transform(d, x = x + 1))
M0 <- model.matrix(~ x * f, d)
M1 <- model.matrix(~ x * f, transform(d, x = x + 1))
s0 <- runif(30)
s1 <- runif(30)

# same values, bit for bit (the dense model matrices also carry `assign`)
same <- function(a, b) {
    a <- as.matrix(a)
    identical(dim(a), dim(b)) && strict_identical(as.vector(a), as.vector(b))
}

# one class for everything
for (D in list(X0, s0 * X0, X1 - X0, s1 * X1 - s0 * X0)) {
    expect_true(inherits(D, "design_matrix"))
    expect_identical(class(D), "design_matrix")
    expect_identical(dim(D), dim(M0))
    expect_identical(dimnames(D), dimnames(as.matrix(X0)))
}

# the dense matrix is the same arithmetic on the dense model matrices
expect_true(same(s0 * X0, s0 * M0))
expect_true(same(X0 * s0, s0 * M0))
expect_true(same(2 * X0, 2 * M0))
expect_true(same(-X0, -1 * M0))
expect_true(same(X0 / 4, (1 / 4) * M0))
expect_true(same(X0 / s0, (1 / s0) * M0))
expect_true(same(X1 - X0, M1 - M0))
expect_true(same(X1 + X0, M1 + M0))
expect_true(same(s1 * X1 - s0 * X0, s1 * M1 - s0 * M0))
expect_true(same(-(s1 * X1 - s0 * X0), -s1 * M1 + s0 * M0))
# scales are combined first
expect_true(same(s1 * (s0 * X0), (s1 * s0) * M0))

# unary plus and scaling nothing away keep the plain design
expect_identical(+X0, X0)

# products and conversions accept combinations
J <- s1 * X1 - s0 * X0
MJ <- s1 * M1 - s0 * M0
b <- rnorm(ncol(M0))
V <- crossprod(matrix(rnorm(36), 6)) / 10
expect_equal(drop(J %*% b), drop(MJ %*% b), check.attributes = FALSE)
expect_equal(drop(crossprod(J, s0)), drop(crossprod(MJ, s0)), check.attributes = FALSE)
expect_equal(quadratic_diagonal(J, V), rowSums((MJ %*% V) * MJ), check.attributes = FALSE)
expect_equal(as.matrix(as(J, "CsparseMatrix")), MJ, check.attributes = FALSE)

# printing
expect_stdout(print(X0), "design_matrix: 30 rows x 6 columns")
expect_stdout(print(s0 * X0), "scaled rows of one design")
expect_stdout(print(X1 - X0), "combination of two designs")

# errors
expect_error(X0 * X1, "between two design matrices")
expect_error(1 / X0, "design first")
expect_error(X0 + 1, "two design matrices")
expect_error(X0 ^ 2, "not supported")
expect_error(1:3 * X0, "length 30")
expect_error(X1 - X0 + X0, "at most two designs")
X2 <- design_matrix(~ x, d)
expect_error(X0 - X2, "same rows and columns")

# fuzz: the dense matrix of a combination is the dense arithmetic, bit for
# bit, with Inf, -0, NA and NaN in the data and in the scales
special_scale <- function(n) {
    s <- rnorm(n)
    s[sample(n, 6)] <- c(0, -0, Inf, -Inf, NA, NaN)
    s
}
set.seed(29)
for (k in 1:60) {
    fd <- fuzz_data(sample(c(8, 40, 200), 1), special = TRUE, na = k %% 2 == 0)
    fo <- fuzz_formula()
    ctr <- fuzz_contrasts(fd)
    fd2 <- fd
    fd2$x1 <- -fd2$x1
    fd2$f1 <- factor(sample(levels(fd$f1), nrow(fd), TRUE), levels = levels(fd$f1))
    D1 <- tryCatch(suppressWarnings(design_matrix(fo, fd, contrasts.arg = ctr)), error = function(e) NULL)
    D2 <- tryCatch(suppressWarnings(design_matrix(fo, fd2, contrasts.arg = ctr)), error = function(e) NULL)
    if (is.null(D1) || is.null(D2) || !identical(dim(D1), dim(D2)) ||
        !identical(dimnames(D1)[[2]], dimnames(D2)[[2]])) next
    M1 <- as.matrix(D1)
    M2 <- as.matrix(D2)
    s1 <- special_scale(nrow(D1))
    s2 <- special_scale(nrow(D1))
    label <- paste(deparse(fo), collapse = " ")
    expect_true(same(s1 * D1, s1 * M1), info = label)
    expect_true(same(s1 * D1 + s2 * D2, s1 * M1 + s2 * M2), info = label)
    expect_true(same(D1 - D2, M1 + (-1) * M2), info = label)
}

# Constant scales and previously expanded scales agree in all row kernels.
weights <- seq_len(nrow(X0)) / nrow(X0)
for (s in c(-2, -0, Inf, NA_real_, NaN)) {
    scalar <- s * X0
    expanded <- rep(s, nrow(X0)) * X0
    expect_identical(length(scalar$scale), 1L)
    expect_true(same(scalar, as.matrix(expanded)))
    expect_true(same(weights * scalar, as.matrix(weights * expanded)))
    expect_identical(design_row_groups(scalar), design_row_groups(expanded))
    expect_identical(scalar %*% b, expanded %*% b)
    expect_identical(crossprod(scalar, weights), crossprod(expanded, weights))
    expect_identical(quadratic_diagonal(scalar, V), quadratic_diagonal(expanded, V))
    expect_identical(as(scalar, "CsparseMatrix"), as(expanded, "CsparseMatrix"))
}

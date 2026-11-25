source("helpers.R")
library(design.matrix)

# x %*% B with many columns uses blocks of dense rows and the BLAS: the
# same values as the dense product up to rounding, NA and NaN in the same
# places, on more rows than one block, for designs and combinations.
K <- design.matrix:::dm_blas_columns + 5L
set.seed(303)
dense_mult <- function(x, B) unname(as.matrix(x) %*% B)
for (k in 1:12) {
    n <- sample(c(30, 700, 1300), 1)
    d <- fuzz_data(n, special = k %% 2 == 0, na = k %% 3 == 0)
    fo <- list(y ~ x1 * f1 + f2, y ~ poly(z, 2) * f1 + ch + lg, y ~ x1:f2 + m)[[1 + k %% 3]]
    D <- design_matrix(design_matrix(fo, d), d)
    B <- matrix(rnorm(ncol(D) * K), ncol = K)
    a <- D %*% B
    b <- dense_mult(D, B)
    expect_identical(is.na(a), is.na(b), info = paste("NA pattern, k =", k))
    expect_identical(is.nan(a), is.nan(b), info = paste("NaN pattern, k =", k))
    expect_equal(a, b, info = paste("values, k =", k))
    s <- runif(n)
    J <- s * D - 2 * design_matrix(D, transform(d, z = z + 1))
    expect_equal(J %*% B, dense_mult(J, B), info = paste("combination, k =", k))
}

# grouping identical rows does not apply here, so it cannot change results
d <- fuzz_data(2000, special = FALSE)
d <- d[sample(nrow(d), 6000, TRUE), ]
D <- design_matrix(design_matrix(y ~ f1 * f2 + o1, d), d)
B <- matrix(rnorm(ncol(D) * K), ncol = K)
with_mode <- function(mode) { op <- options(design.matrix.dedup = mode); on.exit(options(op)); D %*% B }
expect_true(strict_identical(with_mode(TRUE), with_mode(FALSE)))

# threads do not apply to this path either
old <- design_threads()
design_threads(4)
a4 <- D %*% B
design_threads(1)
expect_true(strict_identical(a4, D %*% B))
design_threads(old)

# no rows, and a design with no columns
expect_identical(dim(D[integer(0), ] %*% B), c(0L, K))
D0 <- design_matrix(~ 0, data.frame(x = 1:3))
expect_identical(D0 %*% matrix(0, 0, K), matrix(0, 3, K))

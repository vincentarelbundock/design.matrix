## quadratic_diagonal(D, V) is rowSums((X %*% V) * X) for any finite V: the
## remote design must agree with the dense computation when V is not
## positive semi-definite (it used to keep only the positive eigenvalues of
## V) and when V is not symmetric (only its symmetric part matters).
source("helpers.R")
library(xtx)

d <- data.frame(id = 1:2, y = c(0, 0), x = c(0, 2))
V <- diag(c(1, -1))
expect_equal(quadratic_diagonal(design_matrix(y ~ x, d), V), c(1, -3))

set.seed(7)
n <- 60
d2 <- data.frame(id = seq_len(n), y = 0, x = rnorm(n), f = gl(5, n / 5))
D2 <- design_matrix(y ~ x * f, d2)
X2 <- as.matrix(D2)
p <- ncol(X2)
A <- matrix(rnorm(p * p), p)
Vs <- list(indefinite = diag(rep(c(1, -2), length.out = p)), negative = -crossprod(A), asymmetric = A,
           semidefinite = tcrossprod(A[, 1:2]), zero = matrix(0, p, p), symmetric = A + t(A))
ref <- lapply(Vs, function(V) unname(rowSums((X2 %*% V) * X2)))
for (nm in names(Vs)) expect_equal(quadratic_diagonal(D2, Vs[[nm]]), ref[[nm]], info = nm)

if (requireNamespace("duckdb", quietly = TRUE)) {
    con <- test_con()
    DBI::dbWriteTable(con, "d", d)
    DBI::dbWriteTable(con, "d2", d2)
    rows <- function(col) {
        out <- as.data.frame(col, keep = "id")
        out[[2L]][order(out$id)]
    }
    R <- design_matrix(y ~ x, con, table = "d")
    expect_equal(rows(quadratic_diagonal(R, V)), c(1, -3))
    R2 <- design_matrix(y ~ x * f, con, table = "d2")
    for (nm in names(Vs)) expect_equal(rows(quadratic_diagonal(R2, Vs[[nm]])), ref[[nm]], tolerance = 1e-10, info = nm)
    ## scaled rows: the scale enters both factors of the form
    w <- design_column(R2, "x")
    expect_equal(rows(quadratic_diagonal(w * R2, Vs$indefinite)), d2$x^2 * ref$indefinite, tolerance = 1e-10)
    ## the result combines with other columns of the design
    expect_equal(rows(quadratic_diagonal(R2, Vs$indefinite) + R2 %*% rep(1, p)), ref$indefinite + unname(rowSums(X2)),
                 tolerance = 1e-10)
    test_con_done(con)
}

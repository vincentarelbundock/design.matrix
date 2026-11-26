source("helpers.R")
library(xtx)
if (!requireNamespace("Matrix", quietly = TRUE)) exit_file("Matrix not installed")

sparse_of <- function(x) methods::as(x, "CsparseMatrix")
same_sparse <- function(a, b) {
    identical(class(a), class(b)) && identical(a@Dim, b@Dim) &&
        identical(a@Dimnames, b@Dimnames) && identical(a@p, b@p) &&
        identical(a@i, b@i) && strict_identical(a@x, b@x)
}
# as(<matrix>, "CsparseMatrix") may return a triangular or symmetric class
# for square matrices; as(<design>, "CsparseMatrix") always returns a general dgCMatrix
dense_ref <- function(x) methods::as(methods::as(methods::as(as.matrix(x), "dMatrix"),
                                                  "generalMatrix"), "CsparseMatrix")

# fixed and random designs, including NA rows, Inf, -0, overflow
set.seed(21)
d <- fuzz_data(300, special = TRUE)
for (fo in list(y ~ x1 * f1 + f2, ~ 0 + f1:f2:o1, y ~ poly(z, 3) * ch + m:lg, ~ 1,
                ~ 0 + x1:x2:f1, y ~ f1 * f2)) {
    D <- suppressWarnings(design_matrix(fo, d, contrasts.arg = list(f2 = "contr.sum")))
    expect_true(same_sparse(sparse_of(D), dense_ref(D)), info = deparse(fo))
}
dn <- fuzz_data(200, special = TRUE, na = TRUE)
p <- design_matrix(y ~ x1 * f1 + x2:o1 + ch + lg, dn[stats::complete.cases(dn), ])
D <- design_matrix(p, dn)
expect_true(anyNA(as.matrix(D)))
expect_true(same_sparse(sparse_of(D), dense_ref(D)))

set.seed(22)
for (k in 1:80) {
    dk <- fuzz_data(sample(c(1, 5, 60), 1), special = TRUE, na = stats::runif(1) < 0.3)
    fk <- fuzz_formula(stats::runif(1) < 0.7)
    D <- tryCatch(suppressWarnings(design_matrix(fk, dk, contrasts.arg = fuzz_contrasts(dk))),
                  error = function(e) NULL)
    if (is.null(D)) next
    expect_true(same_sparse(sparse_of(D), dense_ref(D)), info = paste(k, deparse(fk)))
}

# Jacobians, including non-finite scales and cancellation
tr <- fuzz_data(100, special = FALSE)
p <- design_matrix(y ~ x1 * f1 + ch, tr)
X0 <- design_matrix(p, tr)
X1 <- design_matrix(p, transform(tr, x1 = x1 + 1))
s0 <- runif(100); s0[3] <- Inf
J <- runif(100) * X1 - s0 * X0
expect_true(same_sparse(sparse_of(J), dense_ref(J)))
expect_identical(Matrix::nnzero(sparse_of(rep(1, 100) * X0 - X0)), 0L)

# dense fallback and dimensions
fit <- stats::lm(y ~ i1 + f1, data = tr)
expect_true(same_sparse(sparse_of(design_matrix(fit)), dense_ref(design_matrix(fit))))
S <- sparse_of(design_matrix(y ~ f1 * f2, tr))
expect_identical(dim(S), dim(design_matrix(y ~ f1 * f2, tr)))

# both target classes give the same matrix
D <- design_matrix(y ~ f1 * f2 + x1, tr)
expect_true(same_sparse(methods::as(D, "dgCMatrix"), methods::as(D, "CsparseMatrix")))
J <- runif(nrow(tr)) * D
expect_true(same_sparse(methods::as(J, "dgCMatrix"), dense_ref(J)))
expect_false("xtx_sparse" %in% getNamespaceExports("xtx"))

library(Matrix)
library(MatrixModels)
base_sparse <- function(formula, data) {
    X <- model.matrix(formula, data)
    as(as(as(X, "dMatrix"), "generalMatrix"), "CsparseMatrix")
}
fast_sparse <- function(formula, data) {
    D <- design_matrix(formula, data)
    as(D, "CsparseMatrix")
}

sfo <- y ~ x + f1 * f2
reference <- fast_sparse(sfo, big)
same_entries <- function(x) {
    identical(dim(x), dim(reference)) && Matrix::nnzero(x - reference) == 0
}
stopifnot(
    identical(base_sparse(sfo, big), reference),
    same_entries(Matrix::sparse.model.matrix(sfo, big)),
    same_entries(MatrixModels::model.Matrix(sfo, big, sparse = TRUE)))
bench::mark(
    base_sparse(sfo, big),
    fast_sparse(sfo, big),
    Matrix::sparse.model.matrix(sfo, big),
    MatrixModels::model.Matrix(sfo, big, sparse = TRUE),
    iterations = 5, check = FALSE)

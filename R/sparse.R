#' Sparse model matrices
#'
#' `as(x, "CsparseMatrix")` (or `as(x, "dgCMatrix")`) builds a sparse
#' column-compressed matrix (a `Matrix::dgCMatrix`) from a design matrix,
#' without forming the dense matrix. The non-zero values are those
#' of [as.matrix()], bit for bit; exact zeros are dropped and `NA` and `NaN`
#' are kept. The result is identical to converting `as.matrix(x)` to a general
#' `dgCMatrix` (through `"dMatrix"`, `"generalMatrix"` and `"CsparseMatrix"`);
#' unlike a plain `as(<matrix>, "CsparseMatrix")`, it never returns a
#' triangular or symmetric class for square matrices. Dropping zeros loses
#' the sign of `-0`, which the dense matrix keeps.
#'
#' This is not the same matrix as [Matrix::sparse.model.matrix()], which has
#' its own rules for contrasts and column names; it follows
#' [stats::model.matrix()].
#'
#' @param from A `design_matrix`.
#' @return A `dgCMatrix` with the row and column names of `as.matrix(from)`.
#' @examples
#' d <- data.frame(x = rnorm(20), f = gl(4, 5), g = gl(2, 10))
#' D <- design_matrix(~ x + f * g, d)
#' as(D, "CsparseMatrix")
#' @importClassesFrom Matrix CsparseMatrix dgCMatrix
#' @importFrom methods as new setAs setOldClass
#' @importMethodsFrom methods coerce
#' @name sparse
#' @aliases coerce,design_matrix,CsparseMatrix-method coerce,design_matrix,dgCMatrix-method
NULL

## the conversion behind as(<design_matrix>, "CsparseMatrix")
xtx_sparse <- function(x) {
    checkmate::assert_class(x, "xtx_matrix")
    J <- as_jac(x)
    if (!jac_kernel_ok(J)) {
        X <- methods::as(as.matrix(J), "dMatrix")
        return(methods::as(methods::as(X, "generalMatrix"), "CsparseMatrix"))
    }
    parts <- do.call(.Call, c(list(C_xtx_jac_sparse), jac_args(J)))
    d <- dim(x)
    dn <- dimnames(x)
    if (is.null(dn[[1L]]) && is.null(dn[[2L]])) dn <- list(NULL, NULL)
    methods::new("dgCMatrix", i = parts[[2L]], p = parts[[1L]], x = parts[[3L]],
                 Dim = as.integer(d), Dimnames = dn)
}

## as(<design_matrix>, "CsparseMatrix") and as(..., "dgCMatrix")
#' @exportMethod coerce
NULL
setAs("xtx_matrix", "CsparseMatrix", function(from) xtx_sparse(from))
setAs("xtx_matrix", "dgCMatrix", function(from) xtx_sparse(from))

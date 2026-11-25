#' design.matrix: Exact and Fast Model Matrices
#'
#' Builds model matrices that are identical, bit for bit, to those of
#' [stats::model.matrix()], faster, and offers algebra on the design
#' without forming the dense matrix.
#'
#' @useDynLib design.matrix, .registration = TRUE, .fixes = "C_"
#' @import stats
#' @keywords internal
"_PACKAGE"

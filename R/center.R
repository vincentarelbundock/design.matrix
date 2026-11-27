#' Centered designs
#'
#' `design_center()` shifts the numeric variables of a design by their mean
#' (in a sample of rows) when that is an exact linear transform of the
#' design, `X~ = X T` with `T` unit upper triangular: sums of squares and
#' cross products of the centered design keep their accuracy when a variable
#' is far from zero. In memory the design is rebuilt from its model frame;
#' remotely the variables' SQL changes. Coefficients `b~` on the centered
#' design map back to the design as `b = T b~`, and their covariance `V` as
#' `T V T'` (as the fits with fixed effects of 'xtx.fit' do).
#'
#' Names may change before the first release.
#'
#' @param x A design: a `design_matrix` made from a formula and a data frame
#'   (one that keeps its model frame), or a remote design.
#' @param cols Columns of `x` to carry over to the centered design (a list,
#'   possibly nested, of columns and groupings): remotely re-made over the
#'   new design, in memory unchanged.
#' @return `NULL` when no exact centering applies, else a list with `x` (the
#'   centered design), `T` (the transform) and `cols`.
#' @examples
#' d <- data.frame(y = rnorm(100), x = rnorm(100) + 1e6, id = rep(1:10, 10))
#' D <- design_matrix(y ~ x, d)
#' ctr <- design_center(D)
#' colMeans(as.matrix(ctr$x))
#' @keywords internal
#' @export
design_center <- function(x, cols = list()) {
    checkmate::assert_multi_class(x, c("xtx_matrix", "xtx_remote"))
    ctr <- xtx_p_center(x)
    if (is.null(ctr)) return(NULL)
    xc <- ctr$x
    list(x = xc, T = ctr$T, cols = xtx_cols_map(cols, function(col) xtx_c_on(col, xc)))
}

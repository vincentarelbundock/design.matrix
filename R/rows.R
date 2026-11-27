#' Rows of a design, with the columns that go with them
#'
#' `design_subset()` keeps the rows of a design where a logical column is
#' `TRUE`, as `x[keep, ]` does, and carries other columns over the rows to
#' the new design: in memory their values on the rows kept; remotely the
#' same SQL over the new design. With `keep = NULL` no row is dropped, and
#' the design and the columns are only aligned: a remote design is moved to
#' the latest frame among columns computed from it (products, lookups,
#' columns computed once) and the columns are returned as columns of that
#' design, so that the design and the columns can be combined again; in
#' memory the arguments are returned as they are.
#'
#' @param x A design: a `design_matrix` or a remote design.
#' @param keep The rows to keep: a logical vector in memory (missing values
#'   count as `FALSE`), a logical [design_column()] of `x` remotely; or
#'   `NULL`, to keep every row and only align the design and the columns.
#' @param cols Columns over the rows of `x`: a list, possibly nested, of
#'   columns (vectors, matrices or data frames in memory, [design_column()]s
#'   remotely) and `NULL`s. With `keep = NULL` it may also hold bare group
#'   codes (see [key_codes()]), which are carried over; otherwise groupings
#'   are refused, since codes change with the rows.
#' @return A list with `x` (the design) and `cols` (the columns, in the
#'   shape of `cols`).
#' @examples
#' d <- data.frame(y = rnorm(10), x = rnorm(10))
#' D <- design_matrix(y ~ x, d)
#' s <- design_subset(D, d$x > 0, list(y = d$y))
#' nrow(s$x) == length(s$cols$y)
#' @keywords internal
#' @export
design_subset <- function(x, keep = NULL, cols = list()) {
    checkmate::assert_multi_class(x, c("xtx_matrix", "xtx_remote"))
    if (is.null(keep)) {
        x <- xtx_p_rebase(x, cols)
        return(list(x = x, cols = xtx_cols_map(cols, function(col) xtx_c_on(col, x))))
    }
    x2 <- xtx_p_subset(x, keep)
    walk <- function(v) {
        if (is.list(v) && is.null(oldClass(v))) lapply(v, walk) else xtx_c_subset(v, keep, x2)
    }
    list(x = x2, cols = walk(cols))
}

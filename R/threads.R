#' Number of threads used by design.matrix
#'
#' The dense fill and the row-wise kernels (`%*%`, [quadratic_diagonal()])
#' can run on several threads. Each row is computed independently, so the
#' results do not depend on the number of threads.
#'
#' Sums over rows (`crossprod()`, `colSums()` and `colMeans()`, with or
#' without groups) can also run on several threads: the rows are split into
#' fixed blocks, one per thread, whose sums are added in a fixed order.
#' Results are identical for a given number of threads; with different
#' numbers of threads they can differ in the last digits. Small problems
#' run on one thread.
#'
#' The default is one thread, or the value of the option
#' `design.matrix.threads`.
#'
#' @param n Number of threads. `NULL` returns the current setting.
#' @return The number of threads (invisibly when setting).
#' @examples
#' old <- design_threads()
#' design_threads(2)
#' design_threads(old)
#' @export
design_threads <- function(n = NULL) {
    checkmate::assert_count(n, positive = TRUE, null.ok = TRUE)
    if (is.null(n)) {
        n <- getOption("design.matrix.threads", 1L)
        if (n <= 1L) return(1L)
        return(as.integer(min(n, .Call(C_dm_omp_max_threads))))
    }
    options(design.matrix.threads = as.integer(n))
    invisible(as.integer(n))
}

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
#' The default follows fixest: half of the available threads (rounded up),
#' and one thread during `R CMD check`. The option `design.matrix.threads`,
#' which takes the same values as `n`, replaces that default. The available
#' threads are those OpenMP reports, capped by the environment variable
#' `OMP_THREAD_LIMIT`; without OpenMP, one.
#'
#' @param n Number of threads: an integer, `0` for all available threads,
#'   or a number strictly between 0 and 1 for that fraction of them (rounded
#'   up). More threads than available are capped, with a warning. `NULL`
#'   returns the current setting.
#' @return The number of threads (invisibly when setting).
#' @examples
#' old <- design_threads()
#' design_threads(2)
#' design_threads(0.5)
#' design_threads(old)
#' @export
design_threads <- function(n = NULL) {
    if (is.null(n)) {
        n <- getOption("design.matrix.threads")
        if (is.null(n)) n <- if (dm_is_r_check()) 1L else 0.5
        return(dm_resolve_threads(n, warn = FALSE))
    }
    n <- dm_resolve_threads(n, warn = TRUE)
    options(design.matrix.threads = n)
    invisible(n)
}

## The available threads: OpenMP's maximum, capped by OMP_THREAD_LIMIT.
dm_max_threads <- function() {
    lim <- suppressWarnings(as.integer(Sys.getenv("OMP_THREAD_LIMIT")))
    if (is.na(lim) || lim < 1L) lim <- .Machine$integer.max
    max(1L, min(.Call(C_dm_omp_max_threads), lim))
}

## n (a count >= 0, 0 for all, or a fraction in (0, 1)) as a number of
## threads in 1 .. dm_max_threads(), as fixest's check_set_nthreads().
dm_resolve_threads <- function(n, warn) {
    ok <- is.numeric(n) && length(n) == 1L && !is.na(n) &&
        ((n > 0 && n < 1) || (n >= 0 && n == round(n)))
    if (!ok)
        stop("`n` must be a number of threads (0 for all of them) or a fraction strictly between 0 and 1.",
             call. = FALSE)
    m <- dm_max_threads()
    if (n == 0) return(m)
    if (n < 1) return(max(1L, as.integer(ceiling(m * n))))
    if (n > m) {
        if (warn) warning(sprintf("Asked for %d threads while the maximum is %d. Set to %d threads instead.",
                                  as.integer(n), m, m), call. = FALSE)
        return(m)
    }
    as.integer(n)
}

## Under R CMD check, detected as fixest does (many _R_CHECK* variables).
dm_is_r_check <- function() {
    sum(grepl("_R_CHECK", names(Sys.getenv()), fixed = TRUE)) > 10
}

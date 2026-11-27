#' Number of threads used by xtx
#'
#' Resolves the option `xtx.threads` (see the section Options of
#' [xtx-package]) into a number of threads, or sets that option after
#' checking the value. Users set the option directly with
#' `options(xtx.threads = )`; this function is exported for the package
#' 'xtx.fit', whose own kernels need the resolved number.
#'
#' @param n Number of threads: an integer, `0` for all available threads,
#'   or a number strictly between 0 and 1 for that fraction of them (rounded
#'   up). More threads than available are capped, with a warning. `NULL`
#'   returns the current setting.
#' @return The number of threads (invisibly when setting).
#' @examples
#' design_threads()
#' old <- options(xtx.threads = 2)
#' design_threads()
#' options(old)
#' @keywords internal
#' @export
design_threads <- function(n = NULL) {
    if (is.null(n)) {
        n <- getOption("xtx.threads")
        if (is.null(n)) n <- if (xtx_is_r_check()) 1L else 0.5
        return(xtx_resolve_threads(n, warn = FALSE))
    }
    n <- xtx_resolve_threads(n, warn = TRUE)
    options(xtx.threads = n)
    invisible(n)
}

## The available threads: OpenMP's maximum, capped by OMP_THREAD_LIMIT.
xtx_max_threads <- function() {
    lim <- suppressWarnings(as.integer(Sys.getenv("OMP_THREAD_LIMIT")))
    if (is.na(lim) || lim < 1L) lim <- .Machine$integer.max
    max(1L, min(.Call(C_xtx_omp_max_threads), lim))
}

## n (a count >= 0, 0 for all, or a fraction in (0, 1)) as a number of
## threads in 1 .. xtx_max_threads(), as fixest's check_set_nthreads().
xtx_resolve_threads <- function(n, warn) {
    ok <- is.numeric(n) && length(n) == 1L && !is.na(n) &&
        ((n > 0 && n < 1) || (n >= 0 && n == round(n)))
    if (!ok)
        stop("`n` must be a number of threads (0 for all of them) or a fraction strictly between 0 and 1.",
             call. = FALSE)
    m <- xtx_max_threads()
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
xtx_is_r_check <- function() {
    sum(grepl("_R_CHECK", names(Sys.getenv()), fixed = TRUE)) > 10
}

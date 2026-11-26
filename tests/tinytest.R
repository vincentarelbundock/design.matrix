# Default (fast) tier: every feature in memory and on a DuckDB table. Set
# XTX_TESTS=full for the full tier (Parquet, every parity loop,
# the slow regression tests); see inst/tinytest/helpers.R.
#
# The test files run in parallel worker processes: XTX_TEST_NCPU
# of them if set, one (in this process, as tinytest does by default) when
# _R_CHECK_LIMIT_CORES_ is set (CRAN), else min(4, cores - 1). Each worker
# gets an equal share of the threads (option xtx.threads), so the
# workers do not oversubscribe the cores; files that set threads themselves
# (design_threads()) restore the worker's value.
if (requireNamespace("tinytest", quietly = TRUE)) {
    ncpu <- suppressWarnings(as.integer(Sys.getenv("XTX_TEST_NCPU")))
    cores <- parallel::detectCores()
    if (is.na(cores)) cores <- 1L
    if (is.na(ncpu)) {
        limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_")
        ncpu <- if (nzchar(limit) && !identical(toupper(limit), "FALSE")) 1L else min(4L, cores - 1L)
    }
    if (ncpu > 1L) {
        cl <- parallel::makeCluster(ncpu, outfile = "")
        parallel::clusterCall(cl, function(n) options(xtx.threads = n),
                              max(1L, cores %/% ncpu))
        tryCatch(tinytest::test_package("xtx", ncpu = cl),
                 finally = parallel::stopCluster(cl))
    } else {
        tinytest::test_package("xtx")
    }
}

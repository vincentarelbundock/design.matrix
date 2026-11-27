# Run the website benchmarks and save their results in results/.
#
#   Rscript website/benchmarks/run.R              # every benchmark
#   Rscript website/benchmarks/run.R crossprod qr # only these
#   BENCH_SIZE=quick Rscript website/benchmarks/run.R   # quick check
#
# Each benchmark runs in a fresh R process, on top of setup.R, and returns
# either a bench::mark() result, saved as its expression, min, median and
# mem_alloc columns, a data frame, saved as is, or NULL (nothing to save:
# the script saved its own files). The website reads these files and never
# runs a benchmark itself.
#
# BENCH_SIZE is "full" (the default: the published results) or "quick"
# (make bench-quick): fewer rows and iterations (setup.R) and a small
# fixed-effects grid (fixed_effects.R), to check that the scripts run. Quick
# results never go to results/: they go to $BENCH_OUT/website (default
# $TMPDIR/xtx-bench-quick, or /tmp/...), printed at the end. The
# workers get the folder as BENCH_RESULTS.
#
# The BLAS runs on one thread (OMP_NUM_THREADS=1, since OpenBLAS built with
# OpenMP ignores OPENBLAS_NUM_THREADS). The option xtx.threads is capped by the same
# OpenMP limit, so the benchmarks that time several threads run without it.
multithreaded <- c("threads")

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE))))
self <- file.path(here, "run.R")
setwd(here)
wanted <- commandArgs(trailingOnly = TRUE)
bench_size <- Sys.getenv("BENCH_SIZE", "full")
if (!bench_size %in% c("full", "quick")) stop("BENCH_SIZE must be \"full\" or \"quick\", not \"", bench_size, "\"")

summarize <- function(out) {
    if (!inherits(out, "bench_mark")) return(out)
    data.frame(
        expression = format(out$expression),
        min = format(out$min),
        median = format(out$median),
        mem_alloc = format(out$mem_alloc))
}

if (identical(wanted[1], "--one")) {
    # Worker: run one benchmark in this process.
    s <- wanted[2]
    shared <- new.env()
    sys.source("setup.R", envir = shared)
    set.seed(1)
    out <- eval(parse(paste0(s, ".R")), new.env(parent = shared))
    if (!is.null(out)) {
        write.csv(summarize(out), file.path(Sys.getenv("BENCH_RESULTS", "results"), paste0(s, ".csv")),
                  row.names = FALSE)
    }
    quit(save = "no")
}

if (!capabilities("profmem")) {
    warning("this R build lacks profmem support: the mem_alloc column will be NA", call. = FALSE)
}

scripts <- setdiff(sub("\\.R$", "", list.files(".", "\\.R$")), c("setup", "run"))
if (length(wanted) > 0) {
    unknown <- setdiff(wanted, scripts)
    if (length(unknown) > 0) stop("unknown benchmark(s): ", paste(unknown, collapse = ", "))
    scripts <- wanted
}
results <- if (bench_size == "full") file.path(here, "results") else {
    tmp <- Sys.getenv("TMPDIR")
    out <- Sys.getenv("BENCH_OUT")
    if (!nzchar(out)) out <- file.path(if (nzchar(tmp)) tmp else "/tmp", "xtx-bench-quick")
    file.path(out, "website")
}
dir.create(results, showWarnings = FALSE, recursive = TRUE)
results <- normalizePath(results)

rscript <- file.path(R.home("bin"), "Rscript")
for (s in scripts) {
    message("running ", s)
    env <- c("OPENBLAS_NUM_THREADS=1", if (!s %in% multithreaded) "OMP_NUM_THREADS=1",
             paste0("BENCH_SIZE=", bench_size), paste0("BENCH_RESULTS=", shQuote(results)))
    status <- system2(rscript, c(shQuote(self), "--one", s), env = env)
    if (status != 0) stop("benchmark ", s, " failed")
}

writeLines(c(
    R.version.string,
    paste("cores:", parallel::detectCores()),
    paste("BLAS:", basename(extSoftVersion()[["BLAS"]]), "(one thread)"),
    paste("run:", format(Sys.Date()), if (bench_size == "quick") "(quick: BENCH_SIZE=quick)")),
    file.path(results, "machine.txt"))
if (bench_size == "quick") message("quick results in ", results)

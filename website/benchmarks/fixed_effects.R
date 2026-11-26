# Fixed-effects OLS: design_lm() against other R packages, after the
# benchmark of pyfixest (itself adapted from Kyle Butts' fixest_benchmarks).
# The full grid is one command, pinned to one core:
#
#   taskset -c 2 Rscript website/benchmarks/run.R fixed_effects
#
# Each cell (package, data, size, fixed effects, regressors) runs in a fresh R
# process (fixed_effects/cell.R): one warmup call, then the median of three
# timed calls. Memory is the extra peak resident memory of the warmup call,
# so it counts allocations made in C, C++ and DuckDB. Once a package's fit
# takes more than `limit` seconds at some size, or fails after that long
# (dbreg's alternating projections can stop without converging), it is not
# run at larger sizes; a cell that runs longer than `kill` seconds is stopped.
# The feols cells also fit the reference coefficients at a tight tolerance,
# after the timed fits: if that reference fails or is stopped, the timing is
# kept, the `reference` column says what happened, and the skip still depends
# on the timed fit only.
#
# Environment variables shrink the grid for a quick check; the results are
# then marked provisional:
#
#   FE_BENCH_N="1e3,1e4,1e5" FE_BENCH_K="1,10" taskset -c 2 Rscript website/benchmarks/run.R fixed_effects
#
# FE_BENCH_PACKAGES selects packages (comma separated, names as in `full`),
# FE_BENCH_LIMIT and FE_BENCH_KILL set `limit` and `kill` in seconds (default
# 120 and 1800); both are saved in results/fixed_effects_machine.txt.
full <- list(
    n = c(1e3, 1e4, 1e5, 1e6, 1e7),
    k = c(1, 5, 10),
    package = c("design_lm", "design_lm (DuckDB)", "feols", "felm", "dbreg (DuckDB)"))
grid_env <- function(name, default, numeric = TRUE) {
    v <- Sys.getenv(name)
    if (v == "") return(default)
    v <- trimws(strsplit(v, ",")[[1]])
    if (numeric) as.numeric(v) else v
}
sizes <- grid_env("FE_BENCH_N", full$n)
ks <- grid_env("FE_BENCH_K", full$k)
packages <- grid_env("FE_BENCH_PACKAGES", full$package, numeric = FALSE)
limit <- as.numeric(Sys.getenv("FE_BENCH_LIMIT", "120"))
kill <- as.numeric(Sys.getenv("FE_BENCH_KILL", "1800"))
provisional <- !(setequal(sizes, full$n) && setequal(ks, full$k) && setequal(packages, full$package))

rscript <- file.path(R.home("bin"), "Rscript")
tmp <- tempfile(fileext = ".csv")
rows <- list()
for (dgp in c("simple", "difficult")) for (nfe in 2:3) for (k in ks) for (package in packages) {
    too_slow <- FALSE
    for (n in sort(sizes)) {
        cell <- data.frame(package = package, dgp = dgp, n = n, nfe = nfe, k = k,
                           seconds = NA_real_, reps = NA_integer_, mem_mb = NA_real_,
                           coef = NA_character_, ref_coef = NA_character_,
                           reference = NA_character_)
        if (too_slow) {
            cell$status <- "skipped"
        } else {
            message(sprintf("%s, %s, %d FE, k = %d, n = %g", package, dgp, nfe, k, n))
            unlink(tmp)
            wall <- Sys.time()
            status <- system2(rscript,
                c("fixed_effects/cell.R", shQuote(package), dgp, n, nfe, k, limit, tmp),
                env = c("OMP_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1"),
                timeout = kill)
            wall <- as.numeric(difftime(Sys.time(), wall, units = "secs"))
            killed <- identical(status, 124L)
            out <- if (file.exists(tmp)) {
                read.csv(tmp, colClasses = c(coef = "character", ref_coef = "character",
                                             fit_status = "character", ref_status = "character"))
            }
            if (!is.null(out) && out$fit_status == "ok") {
                # The timed fit finished: keep it, whatever the reference did
                cell[c("n", "seconds", "reps", "mem_mb", "coef", "ref_coef")] <-
                    out[c("n", "seconds", "reps", "mem_mb", "coef", "ref_coef")]
                cell$status <- if (out$reps == 0) "single" else "ok"
                cell$reference <- if (is.na(out$ref_status)) NA_character_
                    else if (out$ref_status == "pending") (if (killed) "killed" else "error")
                    else out$ref_status
                fit_seconds <- out$fit_seconds
            } else {
                if (!is.null(out)) cell$n <- out$n
                cell$status <- if (killed) "killed" else "error"
                # The time to failure of the fit, or of the whole process if
                # it died before saying so
                fit_seconds <- if (!is.null(out)) out$fit_seconds else wall
            }
            # A fit that failed (e.g. did not converge) after more than `limit`
            # seconds would only fail more slowly at larger sizes.
            too_slow <- cell$status %in% c("single", "killed") || isTRUE(cell$seconds > limit) ||
                (cell$status == "error" && isTRUE(fit_seconds > limit))
        }
        rows[[length(rows) + 1]] <- cell
    }
}
res <- do.call(rbind, rows)

# Largest absolute difference with feols's coefficients at a tight tolerance
key <- paste(res$dgp, res$n, res$nfe, res$k)
ref <- res$ref_coef[res$package == "feols"]
names(ref) <- key[res$package == "feols"]
split_num <- function(s) if (is.na(s)) NA_real_ else as.numeric(strsplit(s, ";")[[1]])
res$max_coef_diff <- mapply(function(b, r) {
    if (is.na(b) || is.na(r)) return(NA_real_)
    signif(max(abs(split_num(b) - split_num(r))), 3)
}, res$coef, ref[key])

res$seconds <- signif(res$seconds, 4)
res$provisional <- provisional

# The machine, the package versions and the thresholds, next to the results.
# Only the packages that ran are looked up, and a failed lookup never stops
# the results from being saved.
deps <- list("design_lm" = "design.matrix", "design_lm (DuckDB)" = c("design.matrix", "duckdb"),
             "feols" = "fixest", "felm" = "lfe", "dbreg (DuckDB)" = c("dbreg", "duckdb"))
used <- unique(unlist(deps[packages]))
versions <- vapply(used, function(p) tryCatch(as.character(packageVersion(p)),
                                              error = function(e) "not installed"), "")
cpu <- tryCatch(trimws(sub(".*:", "", grep("^model name", readLines("/proc/cpuinfo"), value = TRUE)[1])),
                error = function(e) "unknown")
tryCatch(writeLines(c(
    paste("CPU:", cpu, "(one core, one thread)"),
    R.version.string,
    paste("Packages:", paste(names(versions), versions, collapse = ", ")),
    sprintf("Limits: larger sizes skipped after %g s (FE_BENCH_LIMIT); cells stopped after %g s (FE_BENCH_KILL)",
            limit, kill),
    paste("Run:", format(Sys.Date()), if (provisional) "(provisional: reduced grid)")),
    file.path("results", "fixed_effects_machine.txt")),
    error = function(e) warning("could not write fixed_effects_machine.txt: ", conditionMessage(e)))

res[, c("package", "dgp", "n", "nfe", "k", "seconds", "reps", "mem_mb",
        "max_coef_diff", "status", "provisional", "reference")]

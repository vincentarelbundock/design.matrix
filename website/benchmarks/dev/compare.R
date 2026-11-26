# Development benchmark: before/after between two versions of design.matrix
# =========================================================================
#
# Compare two installed libraries (each containing design.matrix):
#
#   Rscript website/benchmarks/dev/compare.R --base <libpath> --new <libpath>
#
# or two git refs, each built once into a cache (by commit) with
# `git archive <ref> | tar -x; R CMD INSTALL -l <lib>`:
#
#   Rscript website/benchmarks/dev/compare.R --base-ref main --new-ref HEAD
#
# (refs are committed trees: commit before comparing; uncommitted changes are
# not included). With only --base / --base-ref, one version is timed.
#
# Options (defaults in brackets):
#   --cases fe,glm,...  groups or case names [all]. Groups: fe (FE OLS), glm
#                       (FE Poisson/logit), nofe (no FE, in memory), remote
#                       (DuckDB, with query counts), products (crossprod,
#                       colSums(group)). Case names: see cases.R, or --list.
#   --preset quick|standard|large   n = 1e5 / 1e6 / 1e7 [quick]
#   --n 1e6             rows (overrides the preset)
#   --reps 5            timed runs after one warmup (median) [quick 3, else 5]
#   --limit 60          seconds: a warmup slower than this is reported alone
#   --fixest-limit 5    same for the fixest timing (difficult 3FE: 30-50 s/call)
#   --no-fixest         skip fixest timing and references (d_ref then empty)
#   --threads 1         threads for design.matrix, fixest, data.table, DuckDB,
#                       OpenMP and BLAS [1]
#   --core 2            taskset CPU list ("2", "2-5"); "none" to not pin [2]
#   --csv path          CSV output [~/.cache/design.matrix-bench/results/<time>.csv]
#   --cache dir         where refs are built [~/.cache/design.matrix-bench]
#   --list              list the cases and exit
#
# Each case runs in a fresh R process per library (worker.R): simulate data,
# warmup, then reps timed runs. The base process also times fixest (default
# settings) and computes the reference (fixest with tight fixef.tol for FE
# models; lm/glm without FE; model.matrix() for products).
#
# Columns: base_s / new_s median seconds; ratio = new/base; fixest_s and
# vs_fixest = new/fixest; base_q / new_q DBI queries in one remote fit;
# d_base = max |coef new - coef base|; d_ref = max |coef new - reference|.
#
# Wall time (core 2, base + new): quick, all cases, about 3 minutes (mostly
# process startup and the tight fixest references); standard fe,glm took 26
# minutes before --fixest-limit (fixest on the difficult 3FE GLMs: 30-50 s per
# call; the tight reference is still computed once per case). At the quick preset many timings are a few ms: a ratio off
# by 10-40% on one case is noise; rerun that case (--cases name) with more
# --reps, or use the standard preset, before concluding.
#
# Cases are small functions in cases.R; add one there. The comparator for
# development is fixest only.
#
# Examples:
#   # a branch against main, quick, FE only
#   Rscript website/benchmarks/dev/compare.R --base-ref main --new-ref HEAD --cases fe
#   # standard size, FE OLS + GLM, two prebuilt libs
#   Rscript website/benchmarks/dev/compare.R --base /tmp/lib_a --new /tmp/lib_b \
#       --cases fe,glm --preset standard

args <- commandArgs(trailingOnly = TRUE)
# an environment, not a list: `$` must not partially match (--n vs --new)
opt <- list2env(list(cases = "all", preset = "quick", threads = "1", core = "2", limit = "60", fixest_limit = "5",
                     cache = file.path(Sys.getenv("HOME"), ".cache", "design.matrix-bench")))
i <- 1
while (i <= length(args)) {
    a <- args[i]
    if (!startsWith(a, "--")) stop("unexpected argument: ", a)
    key <- gsub("-", "_", substring(a, 3))
    if (key %in% c("list", "no_fixest")) {
        opt[[key]] <- TRUE
        i <- i + 1
        next
    }
    opt[[key]] <- args[i + 1]
    i <- i + 2
}

here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))))
source(file.path(here, "cases.R"))

if (isTRUE(opt$list)) {
    print(data.frame(case = names(case_group), group = unname(case_group)), right = FALSE)
    quit(save = "no")
}

presets <- list(quick = c(n = 1e5, reps = 3), standard = c(n = 1e6, reps = 5),
                large = c(n = 1e7, reps = 5))
if (!opt$preset %in% names(presets)) stop("unknown preset: ", opt$preset)
n <- as.numeric(opt$n %||% presets[[opt$preset]][["n"]])
reps <- as.integer(opt$reps %||% presets[[opt$preset]][["reps"]])
threads <- as.integer(opt$threads)

## Cases
sel <- strsplit(opt$cases, ",")[[1]]
run <- if (identical(sel, "all")) names(cases) else {
    bad <- setdiff(sel, c(names(cases), case_group))
    if (length(bad)) stop("unknown cases or groups: ", paste(bad, collapse = ", "))
    names(cases)[names(cases) %in% sel | case_group %in% sel]
}

## Libraries: given, or built from git refs into the cache
build_lib <- function(ref) {
    repo <- system2("git", c("-C", here, "rev-parse", "--show-toplevel"), stdout = TRUE)
    sha <- system2("git", c("-C", repo, "rev-parse", "--verify", paste0(ref, "^{commit}")), stdout = TRUE)
    lib <- file.path(opt$cache, "libs", substr(sha, 1, 12))
    if (file.exists(file.path(lib, "design.matrix", "DESCRIPTION"))) {
        message("using cached build of ", ref, " (", substr(sha, 1, 12), "): ", lib)
        return(lib)
    }
    message("building ", ref, " (", substr(sha, 1, 12), ") into ", lib)
    src <- tempfile("dmsrc")
    dir.create(src)
    dir.create(lib, recursive = TRUE, showWarnings = FALSE)
    status <- system(sprintf("git -C %s archive %s | tar -x -C %s", shQuote(repo), sha, shQuote(src)))
    if (status != 0) stop("git archive failed for ", ref)
    log <- file.path(lib, "..", paste0(basename(lib), ".install.log"))
    status <- system2("R", c("CMD", "INSTALL", "--no-test-load", "-l", shQuote(lib), shQuote(src)),
                      stdout = log, stderr = log)
    if (status != 0) {
        unlink(lib, recursive = TRUE)
        stop("R CMD INSTALL failed for ", ref, "; see ", log)
    }
    unlink(src, recursive = TRUE)
    lib
}
base <- if (!is.null(opt$base_ref)) build_lib(opt$base_ref) else opt$base
new <- if (!is.null(opt$new_ref)) build_lib(opt$new_ref) else opt$new
if (is.null(base)) stop("give --base <lib> or --base-ref <ref>")
for (l in c(base, new))
    if (!file.exists(file.path(l, "design.matrix", "DESCRIPTION")))
        stop("no design.matrix in library ", l)
base <- normalizePath(base)
if (!is.null(new)) new <- normalizePath(new)

## Run one case with one library, pinned, in a fresh process
thread_env <- sprintf(c("OMP_NUM_THREADS=%d", "OPENBLAS_NUM_THREADS=%d", "MKL_NUM_THREADS=%d",
                        "OMP_THREAD_LIMIT=%d", "R_DATATABLE_NUM_THREADS=%d"), threads)
run_one <- function(lib, case, fixest) {
    out <- tempfile(fileext = ".rds")
    cmd <- c(file.path(R.home("bin"), "Rscript"), file.path(here, "worker.R"),
             lib, case, n, reps, threads, opt$limit, as.integer(fixest), out, opt$fixest_limit)
    if (opt$core != "none") cmd <- c("taskset", "-c", opt$core, cmd)
    log <- tempfile(fileext = ".log")
    status <- system2("env", c(thread_env, shQuote(cmd)), stdout = log, stderr = log)
    if (status != 0 || !file.exists(out)) {
        message("  FAILED: ", case, " (", lib, ")\n", paste(tail(readLines(log), 15), collapse = "\n"))
        return(NULL)
    }
    readRDS(out)
}

maxdiff <- function(a, b) {
    if (is.null(a) || is.null(b)) return(NA_real_)
    if (!is.null(names(a)) && !is.null(names(b)) && all(names(a) %in% names(b))) b <- b[names(a)]
    if (length(a) != length(b)) return(NA_real_)
    max(abs(unname(a) - unname(b)))
}

message(sprintf("n = %g, reps = %d, threads = %d, core = %s\nbase: %s\nnew:  %s",
                n, reps, threads, opt$core, base, new %||% "-"))
rows <- list()
for (case in run) {
    message("* ", case)
    b <- run_one(base, case, fixest = !isTRUE(opt$no_fixest))
    w <- if (!is.null(new)) run_one(new, case, fixest = FALSE)
    if (is.null(b) && is.null(w)) next
    num <- function(x, f) if (is.null(x) || is.null(x[[f]])) NA else x[[f]]
    rows[[case]] <- data.frame(
        case = case, group = case_group[[case]], n = n,
        base_s = num(b, "seconds"), new_s = num(w, "seconds"),
        ratio = num(w, "seconds") / num(b, "seconds"),
        fixest_s = num(b, "fixest_seconds"),
        vs_fixest = (if (is.null(new)) num(b, "seconds") else num(w, "seconds")) / num(b, "fixest_seconds"),
        base_q = num(b, "queries"), new_q = num(w, "queries"),
        d_base = maxdiff(w$coef, b$coef),
        d_ref = maxdiff(if (is.null(new)) b$coef else w$coef, b$ref),
        ref = num(b, "ref_name"),
        reps = num(b, "reps"), fixest_reps = num(b, "fixest_reps"))
}
if (!length(rows)) stop("every case failed")
res <- do.call(rbind, rows)
rownames(res) <- NULL

## Output
show <- res
for (v in c("base_s", "new_s", "fixest_s")) show[[v]] <- sprintf("%.3f", show[[v]])
for (v in c("ratio", "vs_fixest")) show[[v]] <- sprintf("%.2f", show[[v]])
for (v in c("d_base", "d_ref")) show[[v]] <- sprintf("%.1e", show[[v]])
show$n <- format(show$n, scientific = TRUE)
show[show == "NA"] <- ""
show$base_q[is.na(show$base_q)] <- ""
show$new_q[is.na(show$new_q)] <- ""
if (is.null(new)) show$new_s <- show$ratio <- show$new_q <- show$d_base <- NULL
cat("\n")
options(width = 250)
print(show, right = TRUE, row.names = FALSE)

csv <- opt$csv %||% file.path(opt$cache, "results", paste0(format(Sys.time(), "%Y%m%d-%H%M%S"), ".csv"))
dir.create(dirname(csv), recursive = TRUE, showWarnings = FALSE)
res$base_lib <- base
res$new_lib <- new %||% NA
res$threads <- threads
write.csv(res, csv, row.names = FALSE)
cat("\nCSV:", csv, "\n")

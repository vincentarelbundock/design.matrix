# One case of the development benchmark in this (fresh) R process; called by
# compare.R, not meant to be run by hand:
#
#   Rscript worker.R <lib> <case> <n> <reps> <threads> <limit> <fixest 0|1> <out.rds> <fixest_limit>
#
# Loads xtx from <lib>, simulates the case's data, runs fit() once
# (warmup; queries counted for remote cases) and then <reps> timed times
# (only the warmup when it took more than <limit> seconds). With fixest = 1,
# also times the fixest call (warmup alone when slower than <fixest_limit>)
# and computes the reference. Saves a list.
args <- commandArgs(trailingOnly = TRUE)
lib <- args[1]
case <- args[2]
n <- as.numeric(args[3])
reps <- as.integer(args[4])
threads <- as.integer(args[5])
limit <- as.numeric(args[6])
do_fixest <- args[7] == "1"
out <- args[8]
fixest_limit <- as.numeric(args[9])

.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({
    library(xtx, lib.loc = lib)
    library(fixest)
})
design_threads(threads)
setFixest_nthreads(threads)
data.table::setDTthreads(threads)

here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))))
source(file.path(here, "cases.R"))
if (!case %in% names(cases)) stop("unknown case: ", case)

con <- NULL
if (grepl("remote", case) || grepl("remote", case_group[[case]])) {
    con <- DBI::dbConnect(duckdb::duckdb())
    invisible(DBI::dbExecute(con, sprintf("SET threads = %d", threads)))
}

cs <- cases[[case]](n)
if (!is.null(cs$setup)) cs$setup()

# Count DBI queries (the generics xtx calls) during one fit
queries <- new.env()
queries$n <- 0L
count_on <- function() {
    for (g in c("dbGetQuery", "dbExecute", "dbSendQuery", "dbSendStatement"))
        suppressMessages(trace(g, quote(queries$n <- queries$n + 1L), print = FALSE,
                               where = asNamespace("DBI")))
}
count_off <- function() {
    for (g in c("dbGetQuery", "dbExecute", "dbSendQuery", "dbSendStatement"))
        suppressMessages(untrace(g, where = asNamespace("DBI")))
}

# warmup (with query counting when `count`), then reps timed runs with
# bench::mark(): median seconds, R memory allocated per run (needs
# capabilities("profmem"); remote fits: R side only, not DuckDB's memory)
# and garbage collections per run
timed <- function(f, count = FALSE, limit = get("limit", globalenv())) {
    invisible(gc())
    if (count) count_on()
    t0 <- system.time(v <- f())[["elapsed"]]
    if (count) count_off()
    if (t0 > limit || reps < 1)
        return(list(value = v, seconds = t0, reps = 0L, mem = NA_real_, gc = NA_real_))
    bm <- bench::mark(f(), iterations = reps, check = FALSE, filter_gc = FALSE,
                      memory = capabilities("profmem"))
    list(value = v, seconds = as.numeric(bm$median), reps = reps,
         mem = if (capabilities("profmem")) as.numeric(bm$mem_alloc) else NA_real_,
         gc = bm$n_gc / reps)
}

res <- timed(cs$fit, count = isTRUE(cs$remote))

row <- list(
    case = case, group = case_group[[case]], n = n,
    version = as.character(packageVersion("xtx", lib.loc = lib)),
    seconds = res$seconds, reps = res$reps, mem = res$mem, gc = res$gc,
    queries = if (isTRUE(cs$remote)) queries$n else NA_integer_,
    coef = res$value, ref_name = cs$ref_name %||% NA_character_)

if (do_fixest) {
    if (!is.null(cs$fixest)) {
        fx <- timed(cs$fixest, limit = fixest_limit)
        row$fixest_seconds <- fx$seconds
        row$fixest_reps <- fx$reps
    }
    if (!is.null(cs$ref)) row$ref <- cs$ref()
}

if (!is.null(con)) DBI::dbDisconnect(con, shutdown = TRUE)
saveRDS(row, out)

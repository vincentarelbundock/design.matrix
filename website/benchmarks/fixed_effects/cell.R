# One cell of the fixed-effects benchmark (see ../fixed_effects.R), in a fresh
# R process:
#
#   Rscript cell.R <package> <dgp> <n> <nfe> <k> <limit> <out.csv>
#
# Simulates the data, fits the model once (warmup), then three timed times,
# and writes one row to <out.csv>: the median time in seconds, the extra peak
# resident memory of the warmup call in MB, and the coefficients. When the
# warmup takes more than <limit> seconds, its time is reported alone (reps =
# 0) instead of three more runs. If the warmup fails, the row says so
# (fit_status = "error") with the seconds it took to fail.
#
# For feols, the row is written before the reference fit at a tight
# tolerance (ref_status = "pending"), then rewritten with its coefficients
# (ref_status = "ok") or without (ref_status = "error"): a reference that
# fails, or is stopped by the driver, never loses the timed result.
args <- commandArgs(trailingOnly = TRUE)
package <- args[1]
dgp_type <- args[2]
n <- as.numeric(args[3])
nfe <- as.integer(args[4])
k <- as.integer(args[5])
limit <- as.numeric(args[6])
out <- args[7]

# The data generating process of Kyle Butts' fixest_benchmarks, as used by
# pyfixest: a panel of individuals over 10 years, about 23 individuals per
# firm. "simple": each row's firm is drawn at random. "difficult": firms are
# assigned round robin, so individuals move between firms every year and the
# fixed effects are poorly separated (slow convergence for demeaning).
dgp <- function(n, type, k, seed = 42) {
    set.seed(seed)
    nb_year <- 10
    nb_indiv <- round(n / nb_year)
    nb_firm <- round(nb_indiv / 23)
    n <- nb_indiv * nb_year
    indiv_id <- rep(seq_len(nb_indiv), each = nb_year)
    year <- rep(seq_len(nb_year), nb_indiv)
    firm_id <- if (type == "simple") sample.int(nb_firm, n, TRUE) else rep_len(seq_len(nb_firm), n)
    # x1..xk ~ N(0, 1), column by column (same draws as an n by k matrix)
    d <- data.frame(y = numeric(n))
    y <- numeric(n)
    for (j in seq_len(k)) {
        x <- rnorm(n)
        d[[paste0("x", j)]] <- x
        y <- y + x / j
    }
    firm_fe <- rnorm(nb_firm)
    unit_fe <- rnorm(nb_indiv)
    year_fe <- rnorm(nb_year)
    d$y <- y + firm_fe[firm_id] + unit_fe[indiv_id] + year_fe[year] + rnorm(n)
    d$indiv_id <- indiv_id
    d$year <- year
    d$firm_id <- firm_id
    d
}

d <- dgp(n, dgp_type, k)
xs <- paste0("x", seq_len(k))
fe <- c("indiv_id", "year", "firm_id")[seq_len(nfe)]
f <- as.formula(paste("y ~", paste(xs, collapse = " + "), "|", paste(fe, collapse = " + ")))
environment(f) <- globalenv()

# Every package on one thread
remote <- package %in% c("design_lm (DuckDB)", "dbreg (DuckDB)")
if (remote) {
    suppressMessages({
        con <- DBI::dbConnect(duckdb::duckdb())
        DBI::dbExecute(con, "SET threads = 1")
        DBI::dbWriteTable(con, "d", d)
    })
}
suppressPackageStartupMessages(switch(package,
    "design_lm" = ,
    "design_lm (DuckDB)" = {
        library(design.matrix)
        design_threads(1)
    },
    "feols" = {
        library(fixest)
        setFixest_nthreads(1)
        data.table::setDTthreads(1)
    },
    "felm" = {
        options(lfe.threads = 1)
        library(lfe)
    },
    "dbreg (DuckDB)" = {
        library(dbreg)
        # The default (100 iterations of alternating projections) stops with
        # an error on the difficult design with three fixed effects.
        options(dbreg.ap_max_iter = 1e4)
    },
    stop("unknown package: ", package)))

fit <- switch(package,
    "design_lm" = function() design_lm(f, d),
    "design_lm (DuckDB)" = function() design_lm(f, con, table = "d"),
    "feols" = function() feols(f, d, vcov = "iid", notes = FALSE),
    "felm" = function() felm(f, d),
    "dbreg (DuckDB)" = function() dbreg(f, conn = con, table = "d", vcov = "iid"))

# Peak resident memory of this process (Linux): reset the high-water mark,
# then read it after the call.
status_mb <- function(field) {
    line <- grep(paste0("^", field, ":"), readLines("/proc/self/status"), value = TRUE)
    as.numeric(gsub("\\D", "", line)) / 1024
}

# One row of results, written to a temporary file then renamed, so that a
# process stopped while writing never leaves a truncated <out.csv>.
save_row <- function(...) {
    row <- data.frame(package = package, dgp = dgp_type, n = nrow(d), nfe = nfe, k = k, ...)
    part <- paste0(out, ".part")
    write.csv(row, part, row.names = FALSE)
    invisible(file.rename(part, out))
}

invisible(gc())
writeLines("5", "/proc/self/clear_refs")
rss0 <- status_mb("VmRSS")
start <- Sys.time()
t0 <- tryCatch(system.time(model <- fit())[["elapsed"]], error = function(e) {
    message("fit failed: ", conditionMessage(e))
    NULL
})
if (is.null(t0)) {
    save_row(seconds = NA_real_, reps = NA_integer_, mem_mb = NA_real_,
             coef = NA_character_, ref_coef = NA_character_, fit_status = "error",
             fit_seconds = as.numeric(difftime(Sys.time(), start, units = "secs")),
             ref_status = NA_character_)
    quit(save = "no", status = 1)
}
mem_mb <- status_mb("VmHWM") - rss0
b <- coef(model)[xs]

if (t0 > limit) {
    seconds <- t0
    reps <- 0L
} else {
    rm(model)
    times <- replicate(3, system.time(fit())[["elapsed"]])
    seconds <- median(times)
    reps <- 3L
}

timed <- list(seconds = seconds, reps = reps, mem_mb = round(mem_mb, 1),
              coef = paste(sprintf("%.17g", b), collapse = ";"))
reference <- package == "feols"
do.call(save_row, c(timed, list(ref_coef = NA_character_, fit_status = "ok", fit_seconds = t0,
                                ref_status = if (reference) "pending" else NA_character_)))

# feols also gives the reference coefficients, with a tight tolerance
if (reference) {
    ref <- tryCatch(
        coef(feols(f, d, vcov = "iid", notes = FALSE, fixef.tol = 1e-11, fixef.iter = 1e5))[xs],
        error = function(e) {
            message("reference fit failed: ", conditionMessage(e))
            NULL
        })
    ok <- !is.null(ref) && length(ref) == length(xs) && all(is.finite(ref))
    do.call(save_row, c(timed, list(
        ref_coef = if (ok) paste(sprintf("%.17g", ref), collapse = ";") else NA_character_,
        fit_status = "ok", fit_seconds = t0, ref_status = if (ok) "ok" else "error")))
}

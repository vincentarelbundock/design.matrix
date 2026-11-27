## Remote (DuckDB) auto-compression: compressed vs uncompressed products and
## fits (kata design.matrix#8jp3). Base R + DBI + duckdb + bench + xtx.
##
## Run from the package root, with the package and xtx.fit (the fits)
## installed in a private library:
##   R CMD INSTALL -l .rlib .
##   R_LIBS=$PWD/.rlib R CMD INSTALL -l .rlib ../xtx.fit
##   R_LIBS=$PWD/.rlib Rscript bench/compression/run.R [grid]
## grid: "small" (n = 1e5, 1e6), "medium" (+ 5e6; default), "full" (+ 5e7).
## Environment overrides: BENCH_N (comma-separated), BENCH_RATIOS,
## BENCH_OPS, BENCH_KINDS, BENCH_REPS, BENCH_SIDE (0 = skip the side grids
## of ENUM keys and 1 thread), BENCH_MAXGB (skip a forced compressed fit
## whose dense pattern design would exceed this many GB, default 4),
## BENCH_OUT (CSV path).
##
## Data: one DuckDB table per (n, ratio, kind, key type), generated in SQL
## with deterministic hashes of the row number:
##   kind "factor":   factors f1..fk (text or ENUM), design ~ f1 + ... + fk,
##                    patterns ~ L^k (k, L chosen for the target ratio);
##   kind "discrete": f1 (5 levels), f2 (4 levels), xd1, xd2 DOUBLE with K
##                    values each, design ~ f1 + f2 + xd1 + xd2, patterns
##                    ~ 20 K^2;
##   kind "wide":     the data of "discrete", design ~ f1 + f2 +
##                    (xd1 + I(xd1^2) + I(xd1^3)) * (xd2 + I(xd2^2) + I(xd2^3))
##                    (15 numeric columns).
## The realized ratio (distinct design patterns / rows) is measured and
## recorded. y = signal + uniform noise, yb binary (logistic), w in 1:3,
## g with 10 groups (colSums(group = "g")).
##
## Modes: compressed (compress = TRUE and option
## xtx.remote_compress = TRUE) and uncompressed (FALSE / FALSE).
## The "auto" decision (xtx_compress_worth_ask(), the internal function the
## package calls; its result is cached per design, so it is timed apart
## on its own) is timed separately: decision_s and the decision it takes.
## FE fits are never compressed by the package; they are timed with
## compress = "auto" and FALSE to show that.
## Timing: bench::mark() per cell (both modes in one call, explicit
## iterations: max(5, BENCH_REPS + 2) for products, BENCH_REPS for fits;
## check = FALSE, the results are compared separately with all.equal(),
## column `agree`); median, mem_alloc (R allocations only: DuckDB's own
## memory is not counted) and n_gc. Peak RSS above the baseline (R and
## DuckDB, which runs in the process) from an untimed run of each mode.
## Each block records the 1-minute load average when it ran.

suppressMessages({
    library(DBI)
    library(duckdb)
    library(xtx)
    library(xtx.fit)
})

args <- commandArgs(trailingOnly = TRUE)
grid <- if (length(args)) args[1L] else "medium"
env_num <- function(name, default) {
    v <- Sys.getenv(name)
    if (nzchar(v)) as.numeric(strsplit(v, ",")[[1L]]) else default
}
env_chr <- function(name, default) {
    v <- Sys.getenv(name)
    if (nzchar(v)) strsplit(v, ",")[[1L]] else default
}
N <- env_num("BENCH_N", switch(grid, small = c(1e5, 1e6), medium = c(1e5, 1e6, 5e6),
                               full = c(1e5, 1e6, 5e6, 5e7), stop("grid: small, medium or full")))
RATIOS <- env_num("BENCH_RATIOS", c(1e-4, 1e-3, 1e-2, 0.1, 0.5))
OPS <- env_chr("BENCH_OPS", c("crossprod", "colSums_group", "lm", "wlm", "glm", "fe"))
KINDS <- env_chr("BENCH_KINDS", c("factor", "discrete", "wide"))
REPS <- env_num("BENCH_REPS", 3)
SIDE <- env_num("BENCH_SIDE", 1) > 0
MAXGB <- env_num("BENCH_MAXGB", 4)
OUT <- Sys.getenv("BENCH_OUT", file.path("bench", "compression", "results",
                                         sprintf("results-%s.csv", grid)))
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)

load1 <- function() as.numeric(strsplit(readLines("/proc/loadavg"), " ")[[1L]][1L])

## number of cells M such that n uniform draws hit about r * n distinct cells
cells_for <- function(n, r) {
    if (r >= 0.99) return(Inf)
    stats::uniroot(function(M) M * (1 - exp(-n / M)) - r * n, c(r * n * 0.999, n * 1e3))$root
}

## the data, its design formula and FE formula
make_data <- function(con, n, r, kind, keytype) {
    M <- cells_for(n, r)
    h <- function(salt, m) sprintf("CAST(hash(i * 7919 + %d) %% %d AS INTEGER)", salt, as.integer(m))
    u <- function(salt) sprintf("(CAST(hash(i * 7919 + %d) %% 1000000 AS DOUBLE) / 1e6)", salt)
    if (kind == "factor") {
        k <- max(1L, min(6L, floor(log2(M))))
        L <- max(2L, round(M^(1 / k)))
        fac <- sprintf("f%d", seq_len(k))
        cols <- c(sprintf("'l' || %s AS %s", vapply(seq_len(k), function(j) h(j, L), ""), fac))
        sig <- sprintf("CAST(%s AS DOUBLE) / %d", h(1, L), L)
        rhs <- paste(fac, collapse = " + ")
        fe <- if (k > 1L) sprintf("y ~ %s | f1", paste(fac[-1L], collapse = " + ")) else "y ~ 1 | f1"
        info <- sprintf("k=%d,L=%d", k, L)
    } else {
        K <- max(1L, round(sqrt(M / 20)))
        fac <- c("f1", "f2")
        cols <- c(sprintf("'l' || %s AS f1", h(1, 5)), sprintf("'l' || %s AS f2", h(2, 4)),
                  sprintf("CAST(%s AS DOUBLE) / 2 AS xd1", h(3, K)), sprintf("CAST(%s AS DOUBLE) / 2 AS xd2", h(4, K)))
        sig <- sprintf("CAST(%s AS DOUBLE) / %d", h(3, K), max(K, 1L))
        if (kind == "wide") {
            ## many numeric columns from two discrete inputs: 15 numeric
            ## columns, 120 sums of products per factor cell uncompressed
            num <- "(xd1 + I(xd1^2) + I(xd1^3)) * (xd2 + I(xd2^2) + I(xd2^3))"
            rhs <- paste("f1 + f2 +", num)
            fe <- paste("y ~ f2 +", num, "| f1")
        } else {
            rhs <- "f1 + f2 + xd1 + xd2"
            fe <- "y ~ xd1 + xd2 + f2 | f1"
        }
        info <- sprintf("K=%d", K)
    }
    dbExecute(con, "DROP TABLE IF EXISTS d")
    dbExecute(con, "DROP TABLE IF EXISTS d0")
    sql <- sprintf(paste("CREATE TABLE d0 AS SELECT %s, %s + %s - 0.5 AS y,",
                         "CASE WHEN %s < 1 / (1 + exp(-(2 * %s - 1))) THEN 1 ELSE 0 END AS yb,",
                         "1 + %s AS w, %s AS g FROM range(%.0f) t(i)"),
                   paste(cols, collapse = ", "), sig, u(101), u(102), sig, h(103, 3), h(104, 10), n)
    dbExecute(con, sql)
    if (keytype == "enum") {
        for (f in fac) {
            dbExecute(con, sprintf("DROP TYPE IF EXISTS %s_t", f))
            dbExecute(con, sprintf("CREATE TYPE %s_t AS ENUM (SELECT DISTINCT %s FROM d0 ORDER BY 1)", f, f))
        }
        sel <- paste(c(sprintf("CAST(%s AS %s_t) AS %s", fac, fac, fac), "* EXCLUDE (%s)"), collapse = ", ")
        dbExecute(con, sprintf(paste("CREATE TABLE d AS SELECT", sel, "FROM d0"),
                               paste(fac, collapse = ", ")))
        dbExecute(con, "DROP TABLE d0")
    } else {
        dbExecute(con, "ALTER TABLE d0 RENAME TO d")
    }
    keys <- if (kind == "factor") fac else c(fac, "xd1", "xd2")
    G <- dbGetQuery(con, sprintf("SELECT COUNT(*) AS G FROM (SELECT DISTINCT %s FROM d)",
                                 paste(keys, collapse = ", ")))$G
    list(rhs = rhs, fe = fe, G = as.numeric(G), info = info, keys = keys)
}

## time the "auto" decision: the internal query of the package, on a fresh
## design (no cache), with or without the group of colSums()
time_decision <- function(D, group = NULL) {
    ns <- asNamespace("xtx")
    keys <- ns$xtx_compress_keys(D)
    g <- character()
    if (!is.null(group)) {
        gc <- ns$xtx_remote_cols(D, group, "group", numeric = FALSE)
        D <- gc$D
        g <- gc$sql
    }
    n <- D$n
    dec <- ns$xtx_compress_worth_ask(D, keys, n, g)
    bm <- bench::mark(ns$xtx_compress_worth_ask(D, keys, n, g), iterations = max(5L, REPS + 2L),
                      check = FALSE, filter_gc = FALSE, memory = FALSE)
    list(time = as.numeric(bm$median), decision = isTRUE(dec), est = D$cache$compress_estimate %||% NA)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

## one operation in one mode, as a function of no arguments returning
## its result (a matrix, or a fit)
op_fun <- function(con, op, dat, D, mode) {
    fo <- stats::as.formula(paste("y ~", dat$rhs))
    fob <- stats::as.formula(paste("yb ~", dat$rhs))
    f <- switch(op,
        crossprod = function() crossprod(D),
        colSums_group = function() colSums(D, group = "g"),
        lm = function() lm_xtx(fo, con, table = "d", compress = mode),
        wlm = function() lm_xtx(fo, con, table = "d", weights = "w", compress = mode),
        glm = function() suppressWarnings(glm_xtx(fob, binomial, con, table = "d", compress = mode)),
        fe = function() lm_xtx(stats::as.formula(dat$fe), con, table = "d", compress = mode))
    ## (the products follow the option; FE fits are never compressed, so
    ## the option is "auto" for them as for the package default)
    opt <- if (op == "fe") "auto" else mode
    function() {
        oo <- options(xtx.remote_compress = opt)
        on.exit(options(oo))
        suppressMessages(f())
    }
}

## peak resident memory of the process (R and DuckDB, which runs in it)
## while f() runs, above the resident memory before, in MB: the peak is
## reset first (Linux, /proc/self/clear_refs). DuckDB keeps freed buffers,
## so this is a lower bound of what an operation needs.
status_kb <- function(field) {
    l <- grep(paste0("^", field, ":"), readLines("/proc/self/status"), value = TRUE)
    as.numeric(gsub("[^0-9]", "", l))
}
with_peak <- function(f) {
    gc(FALSE)
    ok <- tryCatch({ cat("5", file = "/proc/self/clear_refs"); TRUE }, error = function(e) FALSE)
    base <- status_kb("VmRSS")
    v <- f()
    list(value = v, peak_mb = if (ok) (status_kb("VmHWM") - base) / 1024 else NA)
}

## the numbers to compare between modes
result_of <- function(v) if (is.list(v) && !is.null(v$coefficients)) v$coefficients else v
agree <- function(a, b) isTRUE(all.equal(unname(result_of(a)), unname(result_of(b)), tolerance = 1e-6))

## the configurations: the main grid (text keys, default threads) and the
## side grids (ENUM keys; 1 thread) at the larger n and two ratios
configs <- expand.grid(n = N, ratio = RATIOS, kind = KINDS, keytype = "text", threads = "default",
                       stringsAsFactors = FALSE)
if (SIDE) {
    nside <- N[N >= 1e6 & N <= 5e6]
    if (!length(nside)) nside <- max(N)
    rs <- intersect(c(1e-3, 0.1), RATIOS)
    if (!length(rs)) rs <- RATIOS[1L]
    ks <- intersect(c("factor", "discrete"), KINDS)
    if (!length(ks)) ks <- KINDS
    configs <- rbind(configs,
                     expand.grid(n = nside, ratio = rs, kind = ks, keytype = "enum", threads = "default",
                                 stringsAsFactors = FALSE),
                     expand.grid(n = nside, ratio = rs, kind = ks, keytype = "text", threads = "1",
                                 stringsAsFactors = FALSE))
}
configs <- configs[order(configs$n, configs$kind, configs$keytype, configs$ratio), ]
cat(sprintf("%d configurations, ops: %s, reps: %d -> %s\n", nrow(configs), paste(OPS, collapse = ","), REPS, OUT))

results <- if (file.exists(OUT)) utils::read.csv(OUT, stringsAsFactors = FALSE) else NULL
done <- function(cfg, op) {
    if (is.null(results)) return(FALSE)
    any(results$n == cfg$n & results$ratio_target == cfg$ratio & results$kind == cfg$kind &
            results$keytype == cfg$keytype & results$threads == cfg$threads & results$op == op)
}

for (i in seq_len(nrow(configs))) {
    cfg <- configs[i, ]
    todo <- Filter(function(op) !done(cfg, op), OPS)
    if (!length(todo)) next
    con <- dbConnect(duckdb())
    if (cfg$threads != "default") dbExecute(con, sprintf("SET threads = %s", cfg$threads))
    nthreads <- dbGetQuery(con, "SELECT current_setting('threads') AS t")$t
    dat <- make_data(con, cfg$n, cfg$ratio, cfg$kind, cfg$keytype)
    D <- suppressMessages(design_matrix(stats::as.formula(paste("~", dat$rhs)), con, table = "d"))
    p <- ncol(D)
    dec <- time_decision(suppressMessages(design_matrix(stats::as.formula(paste("~", dat$rhs)), con, table = "d")))
    decg <- time_decision(suppressMessages(design_matrix(stats::as.formula(paste("~", dat$rhs)), con, table = "d")),
                          group = "g")
    for (op in todo) {
        ld <- load1()
        product <- op %in% c("crossprod", "colSums_group")
        fs <- list(comp = op_fun(con, op, dat, D, if (op == "fe") "auto" else TRUE),
                   unc = op_fun(con, op, dat, D, FALSE))
        if (op != "fe" && dat$G * p * 8 / 1e9 > MAXGB) fs$comp <- NULL
        ## untimed runs: the results, whether the fit compressed, the peak
        ## memory (and the warm-up of bench::mark())
        runs <- lapply(fs, function(f) tryCatch(with_peak(f), error = function(e) {
            message(sprintf("  %s: %s", op, conditionMessage(e)))
            NULL
        }))
        fs <- fs[!vapply(runs, is.null, NA)]
        runs <- runs[names(fs)]
        it <- if (product) max(5L, REPS + 2L) else REPS
        bm <- bench::mark(exprs = lapply(names(fs), function(m) bquote(fs[[.(m)]]())), env = environment(),
                          iterations = it, check = FALSE, filter_gc = FALSE, memory = TRUE)
        get <- function(m, col) {
            k <- match(m, names(fs))
            if (is.na(k)) return(NA)
            as.numeric(bm[[col]][k])
        }
        fitc <- if (!is.null(runs$comp) && is.list(runs$comp$value)) isTRUE(runs$comp$value$compressed) else NA
        row <- data.frame(n = cfg$n, ratio_target = cfg$ratio, kind = cfg$kind, keytype = cfg$keytype,
                          threads = cfg$threads, nthreads = nthreads, info = dat$info, p = p,
                          patterns = dat$G, ratio = dat$G / cfg$n, op = op, iterations = it,
                          decision_s = if (op == "colSums_group") decg$time else dec$time,
                          auto_compresses = if (op == "colSums_group") decg$decision else dec$decision,
                          compressed_s = get("comp", "median"), uncompressed_s = get("unc", "median"),
                          compressed_mem_mb = get("comp", "mem_alloc") / 2^20,
                          uncompressed_mem_mb = get("unc", "mem_alloc") / 2^20,
                          compressed_n_gc = get("comp", "n_gc"), uncompressed_n_gc = get("unc", "n_gc"),
                          compressed_peak_rss_mb = if (!is.null(runs$comp)) runs$comp$peak_mb else NA,
                          uncompressed_peak_rss_mb = if (!is.null(runs$unc)) runs$unc$peak_mb else NA,
                          fit_compressed = fitc,
                          agree = if (length(runs) == 2L) agree(runs$comp$value, runs$unc$value) else NA,
                          load = ld, stringsAsFactors = FALSE)
        cat(sprintf("n=%-8.0f %-8s %-4s thr=%-3s r=%.1e (G=%.0f p=%d) %-13s comp=%7.3f unc=%7.3f dec=%.3f mem=%.0f/%.0fMB rss=%.0f/%.0fMB agree=%s load=%.1f\n",
                    cfg$n, cfg$kind, cfg$keytype, nthreads, row$ratio, dat$G, p, op,
                    row$compressed_s, row$uncompressed_s, row$decision_s, row$compressed_mem_mb,
                    row$uncompressed_mem_mb, row$compressed_peak_rss_mb, row$uncompressed_peak_rss_mb,
                    row$agree, ld))
        results <- rbind(results, row)
        utils::write.csv(results, OUT, row.names = FALSE)
    }
    dbDisconnect(con, shutdown = TRUE)
}

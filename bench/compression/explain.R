## Where the time goes: every SQL statement that remote products and fits
## run, compressed vs uncompressed, with its wall time, and EXPLAIN ANALYZE
## of representative ones (kata design.matrix#8jp3). Base R + DBI + duckdb
## (the fits: xtx.fit, installed in .rlib with xtx).
##
##   R_LIBS=$PWD/.rlib Rscript bench/compression/explain.R [n]
##
## Wraps the package's internal xtx_query()/xtx_execute() in this R session
## only (assignInNamespace(); the package code is unchanged) to time each
## statement. Writes results/queries.csv (one row per statement) and
## results/explain/*.txt (EXPLAIN ANALYZE plans).

suppressMessages({
    library(DBI)
    library(duckdb)
    library(xtx)
    library(xtx.fit)
})
args <- commandArgs(trailingOnly = TRUE)
n <- if (length(args)) as.numeric(args[1L]) else 5e6
outdir <- file.path("bench", "compression", "results")
dir.create(file.path(outdir, "explain"), recursive = TRUE, showWarnings = FALSE)

con <- dbConnect(duckdb())
## the "discrete" data of run.R at a ratio of about 1e-3: f1 (5 levels),
## f2 (4), xd1, xd2 (K = 16 values each), 5120 patterns at n = 5e6
K <- max(1L, round(sqrt(n * 1e-3 / 20)))
h <- function(salt, m) sprintf("CAST(hash(i * 7919 + %d) %% %d AS INTEGER)", salt, as.integer(m))
u <- function(salt) sprintf("(CAST(hash(i * 7919 + %d) %% 1000000 AS DOUBLE) / 1e6)", salt)
sig <- sprintf("CAST(%s AS DOUBLE) / %d", h(3, K), K)
dbExecute(con, sprintf(paste("CREATE TABLE d AS SELECT 'l' || %s AS f1, 'l' || %s AS f2,",
                             "CAST(%s AS DOUBLE) / 2 AS xd1, CAST(%s AS DOUBLE) / 2 AS xd2,",
                             "%s + %s - 0.5 AS y, CASE WHEN %s < 1 / (1 + exp(-(2 * %s - 1))) THEN 1 ELSE 0 END AS yb,",
                             "1 + %s AS w, %s AS g FROM range(%.0f) t(i)"),
                       h(1, 5), h(2, 4), h(3, K), h(4, K), sig, u(101), u(102), sig, h(103, 3), h(104, 10), n))

ns <- asNamespace("xtx")
log <- new.env()
log$rows <- list()
log$explain <- NULL   # a function(sql) -> file name, or NULL
wrap <- function(fun) {
    force(fun)
    function(D, sql) {
        if (!is.null(log$explain) && grepl("^\\s*(SELECT|WITH)", sql)) {
            f <- log$explain(sql)
            if (!is.null(f)) {
                pl <- DBI::dbGetQuery(ns$xtx_con(D), paste("EXPLAIN ANALYZE", sql))
                writeLines(c(sql, "", pl[[2L]]), f)
            }
        }
        t <- system.time(out <- fun(D, sql), gcFirst = FALSE)[["elapsed"]]
        log$rows[[length(log$rows) + 1L]] <- data.frame(seconds = t, label = label(sql),
                                                        sql = substr(gsub("\\s+", " ", sql), 1, 160))
        out
    }
}
orig_q <- ns$xtx_query
orig_e <- ns$xtx_execute
utils::assignInNamespace("xtx_query", wrap(orig_q), "xtx")
utils::assignInNamespace("xtx_execute", wrap(orig_e), "xtx")

## a label for a statement
label <- function(sql) {
    if (grepl("approx_count_distinct", sql)) return("decision")
    if (grepl("^SELECT DISTINCT", sql)) return("levels")
    if (grepl("^CREATE", sql)) return("create view")
    if (grepl("^SELECT COUNT\\(\\*\\) AS n FROM", sql)) return("count")
    if (grepl("JOIN p ON", sql)) return("within pass (rows)")
    ## the scan of the SQL engine: GROUP BY the design cells (compressed:
    ## a LIMIT, the cap) or the factor cells with sums of numeric products
    if (grepl("^WITH b", sql) && grepl("GROUP BY ALL LIMIT", sql)) return("pattern scan")
    if (grepl("^WITH b", sql)) return("moments scan (rows)")
    if (grepl("LIMIT|TABLESAMPLE", sql)) return("sample")
    if (grepl("^SELECT SUM", sql)) return("sums (rows)")
    "other"
}

D <- suppressMessages(design_matrix(~ f1 + f2 + xd1 + xd2, con, table = "d"))
## 15 numeric columns from the two discrete inputs (run.R's kind "wide")
Dw <- suppressMessages(design_matrix(~ f1 + f2 + (xd1 + I(xd1^2) + I(xd1^3)) * (xd2 + I(xd2^2) + I(xd2^3)),
                                     con, table = "d"))
ops <- list(
    crossprod = function(m) crossprod(D),
    crossprod_wide = function(m) crossprod(Dw),
    colSums_group = function(m) colSums(D, group = "g"),
    lm = function(m) lm_xtx(y ~ f1 + f2 + xd1 + xd2, con, table = "d", compress = m),
    wlm = function(m) lm_xtx(y ~ f1 + f2 + xd1 + xd2, con, table = "d", weights = "w", compress = m),
    glm = function(m) glm_xtx(yb ~ f1 + f2 + xd1 + xd2, binomial, con, table = "d", compress = m),
    ## (FE fits are not compressed; compress = TRUE is an error: "auto")
    fe = function(m) lm_xtx(y ~ xd1 + xd2 + f2 | f1, con, table = "d", compress = if (isTRUE(m)) "auto" else m))
res <- list()
for (op in names(ops)) for (m in list(TRUE, FALSE, "auto")) {
    mode <- if (isTRUE(m)) "compressed" else if (isFALSE(m)) "uncompressed" else "auto"
    oo <- options(xtx.remote_compress = m)
    suppressMessages(ops[[op]](m))           # warm-up
    ## the plans of the first scans of a few kinds (one file per label)
    seen <- character()
    log$explain <- function(sql) {
        lb <- label(sql)
        if (!(lb %in% c("pattern scan", "moments scan (rows)", "sums (rows)", "decision"))) return(NULL)
        k <- sum(seen == lb) + 1L
        seen <<- c(seen, lb)
        if (k > if (op == "glm") 3L else 1L) return(NULL)
        file.path(outdir, "explain", sprintf("%s-%s-%s-%d.txt", op, mode, gsub("[^a-z]+", "_", lb), k))
    }
    suppressMessages(ops[[op]](m))
    log$explain <- NULL
    log$rows <- list()
    tt <- system.time(fit <- suppressMessages(ops[[op]](m)))[["elapsed"]]
    q <- do.call(rbind, log$rows)
    q$op <- op
    q$mode <- mode
    q$i <- seq_len(nrow(q))
    q$total <- tt
    q$iter <- if (is.list(fit) && !is.null(fit$iter)) fit$iter else NA
    res[[length(res) + 1L]] <- q
    options(oo)
    cat(sprintf("%-14s %-12s total %6.3f s, %2d statements, %6.3f s in SQL; by label: %s\n", op, mode, tt, nrow(q),
                sum(q$seconds), paste(sprintf("%s %d/%.2fs", names(tapply(q$seconds, q$label, sum)),
                                              as.vector(table(q$label)),
                                              tapply(q$seconds, q$label, sum)), collapse = "; ")))
}
q <- do.call(rbind, res)
utils::write.csv(q[c("op", "mode", "i", "label", "seconds", "total", "iter", "sql")],
                 file.path(outdir, sprintf("queries-%.0e.csv", n)), row.names = FALSE)

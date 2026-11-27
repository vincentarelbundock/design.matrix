## Compressed remote products are the scan of the SQL engine grouped by the
## cells of the design (factor codes and numeric variables): one query, the
## same results as the scan over the rows and as the dense computation.
if (!requireNamespace("duckdb", quietly = TRUE) || !requireNamespace("DBI", quietly = TRUE))
    exit_file("duckdb is not installed")
source("helpers.R")
library(xtx)

con <- test_con()
same <- function(a, b, tol = 1e-9) isTRUE(all.equal(unname(a), unname(b), tolerance = tol))
with_opts <- function(expr, ...) {
    op <- options(...)
    on.exit(options(op))
    expr
}
## the queries an expression runs
queries <- function(expr, ...) {
    op <- options(xtx.remote_debug = TRUE, ...)
    on.exit(options(op))
    q <- character()
    v <- withCallingHandlers(expr, message = function(m) {
        q <<- c(q, conditionMessage(m))
        invokeRestart("muffleMessage")
    })
    list(value = v, sql = q)
}

set.seed(7)
n <- 2000
d <- data.frame(xd = sample(0:4, n, TRUE) / 2, xe = sample(0:2, n, TRUE), f = sample(c("a", "b", "c"), n, TRUE),
                lg = sample(c(TRUE, FALSE), n, TRUE), w = sample(1:3, n, TRUE),
                g = sample(1:6, n, TRUE), xc = rnorm(n))
d$y <- d$xd + rnorm(n)
d$ym <- d$y
d$ym[5] <- NA
DBI::dbWriteTable(con, "d", d)
fo <- ~ xd * f + lg + log1p(xd) + I(xd^2)
D <- suppressMessages(design_matrix(fo, con, table = "d"))
M <- as.matrix(design_matrix(fo, d))

## TRUE: one query per product, the numeric variables as keys of the
## GROUP BY; the dense results
## (a numeric variable __xtx_v<j> among the keys of the SELECT of the cells)
grouped_scan <- function(sql) grepl("\\(SELECT (__xtx_[a-z0-9]+, )*__xtx_v[0-9]+,", sql) && grepl("GROUP BY ALL", sql, fixed = TRUE)
q <- queries(crossprod(D, "w" * D), xtx.remote_compress = TRUE)
expect_true(same(q$value, crossprod(M, d$w * M)))
expect_equal(length(q$sql), 1L)
expect_true(grouped_scan(q$sql))
q <- queries(colSums("w" * D, group = "g"), xtx.remote_compress = TRUE)
expect_true(same(q$value, rowsum(d$w * M, d$g)))
expect_true(grouped_scan(q$sql))

## the same values in every mode (compression changes the speed, never the
## results beyond rounding)
for (mode in list(TRUE, FALSE, "auto")) {
    info <- paste("mode", mode)
    expect_true(same(with_opts(crossprod(D, "w" * D), xtx.remote_compress = mode),
                     crossprod(M, d$w * M)), info = info)
    expect_true(same(with_opts(crossprod(D, "y"), xtx.remote_compress = mode), crossprod(M, d$y)),
                info = info)
    expect_true(same(with_opts(colSums("w" * D, group = "g"), xtx.remote_compress = mode),
                     rowsum(d$w * M, d$g)), info = info)
    expect_true(same(with_opts(colSums(D), xtx.remote_compress = mode), colSums(M)), info = info)
    expect_true(same(with_opts(colMeans(D), xtx.remote_compress = mode), colMeans(M)), info = info)
}

## "auto" by operation (xtx_compress_rule()): colSums() and colMeans() never
## compress, nor ask; products of a design with few numeric columns (here 3
## distinct products of numeric variables) neither
D <- suppressMessages(design_matrix(fo, con, table = "d"))
expect_equal(xtx:::xtx_compress_numeric(D), 3L)
for (e in list(quote(colSums(D, group = "g")), quote(colMeans(D)), quote(crossprod(D, "w" * D)),
               quote(crossprod(D, "y")))) {
    q <- queries(eval(e), xtx.remote_compress = "auto")
    expect_equal(length(q$sql), 1L, info = deparse(e))
    expect_false(grouped_scan(q$sql), info = deparse(e))
}
## a design with many numeric columns (11) and few patterns (15, under 1% of
## the rows): products compress, after the decision (asked once); colSums()
## still does not
fw <- ~ (xd + I(xd^2) + I(xd^3)) * (xe + I(xe^2))
DW <- suppressMessages(design_matrix(fw, con, table = "d"))
MW <- as.matrix(design_matrix(fw, d))
expect_equal(xtx:::xtx_compress_numeric(DW), 11L)
q <- queries(crossprod(DW, "w" * DW), xtx.remote_compress = "auto")
expect_true(same(q$value, crossprod(MW, d$w * MW)))
expect_true(length(q$sql) > 1L && grouped_scan(q$sql[length(q$sql)]))
q <- queries(crossprod(DW, "y"), xtx.remote_compress = "auto")
expect_true(same(q$value, crossprod(MW, d$y)))
expect_equal(length(q$sql), 1L)
expect_true(grouped_scan(q$sql))
q <- queries(colSums(DW, group = "g"), xtx.remote_compress = "auto")
expect_true(same(q$value, rowsum(MW, d$g)))
expect_false(grouped_scan(q$sql))
## more patterns than 1% of the rows (here with f, 45 for 2000 rows): not
## compressed, though a glm fit (ratio 0.2) would be; an unweighted lm (ratio
## 1e-3) would not, nor a weighted one
cap <- function(D, op) xtx:::xtx_compress_cap(D, "auto", op = op)
DW2 <- suppressMessages(design_matrix(update(fw, ~ . + f), con, table = "d"))
MW2 <- as.matrix(design_matrix(update(fw, ~ . + f), d))
q <- queries(crossprod(DW2), xtx.remote_compress = "auto")
expect_true(same(q$value, crossprod(MW2)))
expect_false(grouped_scan(q$sql[length(q$sql)]))
expect_null(cap(DW2, "products"))
expect_false(is.null(cap(DW2, "glm")))
expect_null(cap(DW2, "lm"))
expect_null(cap(DW2, "wlm"))
## the option xtx.remote_compress_ratio caps every ratio (never
## raises one)
rule <- xtx:::xtx_compress_rule
expect_null(rule("colsums"))
expect_equal(rule("products"), list(numeric = 8L, ratio = 1e-2))
expect_equal(rule("glm"), list(numeric = 0L, ratio = 0.2))
expect_equal(rule("lm"), list(numeric = 0L, ratio = 1e-3))
expect_equal(rule("wlm"), list(numeric = 8L, ratio = 1e-3))
expect_equal(with_opts(rule("glm")$ratio, xtx.remote_compress_ratio = 0.05), 0.05)
expect_equal(with_opts(rule("lm")$ratio, xtx.remote_compress_ratio = 0.05), 1e-3)
expect_equal(with_opts(rule("lm")$ratio, xtx.remote_compress_ratio = 1e-4), 1e-4)
expect_null(with_opts(cap(DW2, "glm"), xtx.remote_compress_ratio = 0.01))

## a missing response: compressed products put NA where the dense one does
xm <- with_opts(crossprod(D, "ym"), xtx.remote_compress = TRUE)
expect_true(all(is.na(xm)))

## too many cells: TRUE gives up after the scan (every row distinct), and
## the products are the scan over the rows
D2 <- design_matrix(~ xc + f, con, table = "d")
M2 <- as.matrix(design_matrix(~ xc + f, d))
expect_true(same(with_opts(crossprod(D2), xtx.remote_compress = TRUE), crossprod(M2)))
expect_null(xtx:::xtx_sql_moments(D2, cap = 10))

## ---- a trial that fails is discarded in the database ----
## The cells of a trial are counted where they are computed: when there are
## more than the cap, or as many as the rows, the trial returns no row, and
## the ordinary scan follows. What each mode fetches (record_queries()):
fetched <- function(expr, ...) {
    op <- options(...)
    on.exit(options(op))
    record_queries(expr)
}
trial <- function(r) grepl("c AS MATERIALIZED", r$log$sql, fixed = TRUE)
## forced on, every row its own cell (2000 cells of xc): nothing of the
## trial comes back; one more statement, the scan by the 3 levels of f
r <- fetched(crossprod(D2), xtx.remote_compress = TRUE)
expect_true(same(r$value, crossprod(M2)))
expect_identical(r$count, 2L)
expect_identical(trial(r), c(TRUE, FALSE))
expect_identical(r$log$rows, c(0L, 3L))
r <- fetched(crossprod(D2, "y"), xtx.remote_compress = TRUE)
expect_true(same(r$value, crossprod(M2, d$y)))
expect_identical(r$log$rows, c(0L, 3L))
r <- fetched(colSums(D2, group = "g"), xtx.remote_compress = TRUE)
expect_true(same(r$value, rowsum(M2, d$g)))
expect_identical(r$log$rows[trial(r)], 0L)
expect_true(all(r$log$rows <= 6L * 3L))
## forced on, few cells (30: 5 values of xd, 3 levels, 2 logicals): the
## cells and nothing else
r <- fetched(crossprod(D, "w" * D), xtx.remote_compress = TRUE)
expect_true(same(r$value, crossprod(M, d$w * M)))
## (the first statement reads the columns of the scaled source, no row)
expect_identical(trial(r), c(FALSE, TRUE))
expect_identical(r$log$rows, c(0L, nrow(unique(d[c("xd", "f", "lg")]))))
## forced off: no trial at all
for (e in list(quote(crossprod(D2)), quote(crossprod(D, "w" * D)), quote(colSums(D, group = "g")))) {
    r <- fetched(eval(e), xtx.remote_compress = FALSE)
    expect_false(any(trial(r)), info = deparse(e))
    expect_true(all(r$log$rows <= 6L * 3L * 2L), info = deparse(e))
}
## automatic, compressed: the decision (two statements of one row), then
## the 15 cells
DW <- suppressMessages(design_matrix(fw, con, table = "d"))
r <- fetched(crossprod(DW, "w" * DW), xtx.remote_compress = "auto")
expect_true(same(r$value, crossprod(MW, d$w * MW)))
expect_identical(trial(r), c(FALSE, FALSE, FALSE, TRUE))
expect_identical(r$log$rows, c(0L, 1L, 1L, nrow(unique(d[c("xd", "xe")]))))
## automatic, the estimate wrong: the cap (here 10 cells) is below the 15
## cells the estimate promised. The trial returns nothing and the ordinary
## scan (one row: the design has no factor) gives the products
DW <- suppressMessages(design_matrix(fw, con, table = "d"))
key <- paste(c(DW$filter, xtx:::xtx_compress_keys(DW)), collapse = "\r")
DW$cache$worth[[key]] <- c(est = 5, rows = n)
r <- fetched(crossprod(DW, "w" * DW), xtx.remote_compress = "auto", xtx.remote_compress_max = 10)
expect_true(same(r$value, crossprod(MW, d$w * MW)))
expect_identical(trial(r), c(FALSE, TRUE, FALSE))
expect_identical(r$log$rows, c(0L, 0L, 1L))
## and the fits' patterns (xtx_fit_suffstats() calls xtx_sql_patterns()):
## NULL, with nothing fetched, over the cap
pt <- fetched(xtx:::xtx_sql_patterns(DW, "CAST(s.y AS DOUBLE)", mode = "auto", fit = "glm"), xtx.remote_compress_max = 10)
expect_null(pt$value)
expect_identical(pt$log$rows[trial(pt)], 0L)
pt <- fetched(xtx:::xtx_sql_patterns(DW, "CAST(s.y AS DOUBLE)", mode = TRUE))
expect_identical(nrow(pt$value$x), 15L)
expect_identical(pt$log$rows[trial(pt)], 15L)
## weighted, every row its own cell: the second pass is never run
pt <- fetched(xtx:::xtx_sql_patterns(D2, "CAST(s.y AS DOUBLE)", "CAST(s.w AS DOUBLE)", mode = TRUE))
expect_null(pt$value)
expect_identical(pt$count, 1L)
expect_identical(pt$log$rows, 0L)

test_con_done(con)

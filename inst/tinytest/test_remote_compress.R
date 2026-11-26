## Compressed remote products are the scan of the SQL engine grouped by the
## cells of the design (factor codes and numeric variables): one query, the
## same results as the scan over the rows and as the dense computation.
if (!requireNamespace("duckdb", quietly = TRUE) || !requireNamespace("DBI", quietly = TRUE))
    exit_file("duckdb is not installed")
library(design.matrix)

con <- DBI::dbConnect(duckdb::duckdb())
same <- function(a, b, tol = 1e-9) isTRUE(all.equal(unname(a), unname(b), tolerance = tol))
with_opts <- function(expr, ...) {
    op <- options(...)
    on.exit(options(op))
    expr
}
## the queries an expression runs
queries <- function(expr, ...) {
    op <- options(design.matrix.remote_debug = TRUE, ...)
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
## (a numeric variable __dm_v<j> among the keys of the outer SELECT)
grouped_scan <- function(sql) grepl("\\) SELECT (__dm_[a-z0-9]+, )*__dm_v[0-9]+,", sql) && grepl("GROUP BY ALL", sql, fixed = TRUE)
q <- queries(crossprod(D, "w" * D), design.matrix.remote_compress = TRUE)
expect_true(same(q$value, crossprod(M, d$w * M)))
expect_equal(length(q$sql), 1L)
expect_true(grouped_scan(q$sql))
q <- queries(colSums("w" * D, group = "g"), design.matrix.remote_compress = TRUE)
expect_true(same(q$value, rowsum(d$w * M, d$g)))
expect_true(grouped_scan(q$sql))

## the same values in every mode (compression changes the speed, never the
## results beyond rounding)
for (mode in list(TRUE, FALSE, "auto")) {
    info <- paste("mode", mode)
    expect_true(same(with_opts(crossprod(D, "w" * D), design.matrix.remote_compress = mode),
                     crossprod(M, d$w * M)), info = info)
    expect_true(same(with_opts(crossprod(D, "y"), design.matrix.remote_compress = mode), crossprod(M, d$y)),
                info = info)
    expect_true(same(with_opts(colSums("w" * D, group = "g"), design.matrix.remote_compress = mode),
                     rowsum(d$w * M, d$g)), info = info)
    expect_true(same(with_opts(colSums(D), design.matrix.remote_compress = mode), colSums(M)), info = info)
    expect_true(same(with_opts(colMeans(D), design.matrix.remote_compress = mode), colMeans(M)), info = info)
}

## "auto" by operation (dm_compress_rule()): colSums() and colMeans() never
## compress, nor ask; products of a design with few numeric columns (here 3
## distinct products of numeric variables) neither
D <- suppressMessages(design_matrix(fo, con, table = "d"))
expect_equal(design.matrix:::dm_compress_numeric(D), 3L)
for (e in list(quote(colSums(D, group = "g")), quote(colMeans(D)), quote(crossprod(D, "w" * D)),
               quote(crossprod(D, "y")))) {
    q <- queries(eval(e), design.matrix.remote_compress = "auto")
    expect_equal(length(q$sql), 1L, info = deparse(e))
    expect_false(grouped_scan(q$sql), info = deparse(e))
}
## a design with many numeric columns (11) and few patterns (15, under 1% of
## the rows): products compress, after the decision (asked once); colSums()
## still does not
fw <- ~ (xd + I(xd^2) + I(xd^3)) * (xe + I(xe^2))
DW <- suppressMessages(design_matrix(fw, con, table = "d"))
MW <- as.matrix(design_matrix(fw, d))
expect_equal(design.matrix:::dm_compress_numeric(DW), 11L)
q <- queries(crossprod(DW, "w" * DW), design.matrix.remote_compress = "auto")
expect_true(same(q$value, crossprod(MW, d$w * MW)))
expect_true(length(q$sql) > 1L && grouped_scan(q$sql[length(q$sql)]))
q <- queries(crossprod(DW, "y"), design.matrix.remote_compress = "auto")
expect_true(same(q$value, crossprod(MW, d$y)))
expect_equal(length(q$sql), 1L)
expect_true(grouped_scan(q$sql))
q <- queries(colSums(DW, group = "g"), design.matrix.remote_compress = "auto")
expect_true(same(q$value, rowsum(MW, d$g)))
expect_false(grouped_scan(q$sql))
## more patterns than 1% of the rows (here with f, 45 for 2000 rows): not
## compressed, though fits (ratio 0.1) would be
DW2 <- suppressMessages(design_matrix(update(fw, ~ . + f), con, table = "d"))
MW2 <- as.matrix(design_matrix(update(fw, ~ . + f), d))
q <- queries(crossprod(DW2), design.matrix.remote_compress = "auto")
expect_true(same(q$value, crossprod(MW2)))
expect_false(grouped_scan(q$sql[length(q$sql)]))
expect_null(design.matrix:::dm_compress_cap(DW2, "auto", op = "products"))
expect_false(is.null(design.matrix:::dm_compress_cap(DW2, "auto", op = "fit")))

## a missing response: compressed products put NA where the dense one does
xm <- with_opts(crossprod(D, "ym"), design.matrix.remote_compress = TRUE)
expect_true(all(is.na(xm)))

## compressed and uncompressed fits: weighted (two passes for the sum of
## squares within the patterns) and with an offset in the keys (design_glm)
for (w in list(NULL, "w")) {
    a <- design_lm(y ~ xd * f + lg, con, table = "d", weights = w, compress = TRUE)
    b <- design_lm(y ~ xd * f + lg, con, table = "d", weights = w, compress = FALSE)
    l <- lm(y ~ xd * f + lg, d, weights = if (is.null(w)) NULL else d$w)
    expect_true(a$compressed)
    expect_false(b$compressed)
    expect_true(same(coef(a), coef(l)) && same(coef(b), coef(l)))
    expect_true(same(a$rss, b$rss) && same(a$rss, sum(weighted.residuals(l)^2)))
}
d$yp <- rpois(n, exp(0.2 * d$xd))
DBI::dbWriteTable(con, "d", d, overwrite = TRUE)
a <- design_glm(yp ~ xd + f, poisson, con, table = "d", offset = "g", compress = TRUE)
g0 <- glm(yp ~ xd + f, poisson, d, offset = g)
expect_true(a$compressed)
expect_true(same(coef(a), coef(g0), 1e-7) && same(deviance(a), deviance(g0), 1e-7))

## too many cells: TRUE gives up after the scan (every row distinct), and
## the products are the scan over the rows
D2 <- design_matrix(~ xc + f, con, table = "d")
M2 <- as.matrix(design_matrix(~ xc + f, d))
expect_true(same(with_opts(crossprod(D2), design.matrix.remote_compress = TRUE), crossprod(M2)))
expect_null(design.matrix:::dm_sql_moments(D2, cap = 10))

DBI::dbDisconnect(con, shutdown = TRUE)

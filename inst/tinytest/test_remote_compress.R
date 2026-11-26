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
d <- data.frame(xd = sample(0:4, n, TRUE) / 2, f = sample(c("a", "b", "c"), n, TRUE),
                lg = sample(c(TRUE, FALSE), n, TRUE), w = sample(1:3, n, TRUE),
                g = sample(1:6, n, TRUE), xc = rnorm(n))
d$y <- d$xd + rnorm(n)
d$ym <- d$y
d$ym[5] <- NA
DBI::dbWriteTable(con, "d", d)
fo <- ~ xd * f + lg + log1p(xd) + I(xd^2)
D <- suppressMessages(design_matrix(fo, con, table = "d"))
M <- as.matrix(design_matrix(fo, d))

## one query per product once the decision is cached; the dense results
for (mode in list(TRUE, "auto")) {
    info <- paste("mode", mode)
    q <- queries(crossprod(D, "w" * D), design.matrix.remote_compress = mode)
    expect_true(same(q$value, crossprod(M, d$w * M)), info = info)
    q <- queries(crossprod(D, "w" * D), design.matrix.remote_compress = mode)
    expect_equal(length(q$sql), 1L, info = info)
    ## (the numeric variables are keys of the GROUP BY)
    expect_true(grepl("SELECT __dm_v", q$sql, fixed = TRUE) && grepl("GROUP BY ALL", q$sql, fixed = TRUE),
                info = info)
    expect_true(same(with_opts(crossprod(D, "y"), design.matrix.remote_compress = mode), crossprod(M, d$y)),
                info = info)
    expect_true(same(with_opts(colSums("w" * D, group = "g"), design.matrix.remote_compress = mode),
                     rowsum(d$w * M, d$g)), info = info)
    expect_true(same(with_opts(colMeans(D), design.matrix.remote_compress = mode), colMeans(M)), info = info)
}

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

## too many cells: "auto" gives up after the scan, and the products are the
## scan over the rows
D2 <- design_matrix(~ xc + f, con, table = "d")
M2 <- as.matrix(design_matrix(~ xc + f, d))
expect_true(same(with_opts(crossprod(D2), design.matrix.remote_compress_ratio = 1,
                           design.matrix.remote_compress_max = 1e9), crossprod(M2)))
expect_null(design.matrix:::dm_sql_moments(D2, cap = 10))

DBI::dbDisconnect(con, shutdown = TRUE)

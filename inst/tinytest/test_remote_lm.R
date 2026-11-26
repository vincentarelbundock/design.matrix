## Compressed remote products and design_lm() on remote data: agreement with
## lm() and with the uncompressed computation.
if (!requireNamespace("duckdb", quietly = TRUE) || !requireNamespace("DBI", quietly = TRUE))
    exit_file("duckdb is not installed")
source("helpers.R")
library(design.matrix)

con <- test_con()
same <- function(a, b, tol = 1e-7) isTRUE(all.equal(unname(a), unname(b), tolerance = tol))
no_compress <- function(expr) {
    op <- options(design.matrix.remote_compress = FALSE)
    on.exit(options(op))
    expr
}

set.seed(42)
n <- 3000
d <- data.frame(
    xd = sample(1:5, n, TRUE),
    xc = rnorm(n),
    f1 = sample(c("b", "a", "c", "d"), n, TRUE),
    f2 = factor(sample(c("u", "v", "w"), n, TRUE), levels = c("w", "u", "v", "unused")),
    lg = sample(c(TRUE, FALSE), n, TRUE),
    w = runif(n, 0.5, 2),
    cl = sample(1:30, n, TRUE),
    id = seq_len(n))
d$y <- d$xd + d$xc + as.integer(factor(d$f1)) + rnorm(n)
d$xs <- 100 * d$xc + 1e7  # far from zero, not aliased by lm()
d$xt <- d$xc + 1e7        # so close to aliased that lm() drops it
d$xd7 <- d$xd + 1e7
d$f3 <- d$f1              # aliased with f1
d$y[c(3, 50)] <- NA
d$f1[7] <- NA
d$xc[9] <- NaN
d$w[11] <- 0
DBI::dbWriteTable(con, "d", d)
dd <- DBI::dbReadTable(con, "d")

check_fit <- function(f, l, info, tol = 1e-7) {
    expect_identical(is.na(coef(f)), is.na(coef(l)), info = paste(info, "aliased"))
    expect_true(same(coef(f), coef(l), tol), info = paste(info, "coef"))
    expect_true(same(vcov(f), vcov(l), tol), info = paste(info, "vcov"))
    expect_true(same(sigma(f), sigma(l), tol), info = paste(info, "sigma"))
    expect_identical(df.residual(f), df.residual(l), info = paste(info, "df"))
    expect_identical(as.numeric(nobs(f)), as.numeric(nobs(l)), info = paste(info, "nobs"))
    expect_true(same(summary(f)$r.squared, summary(l)$r.squared, tol), info = paste(info, "R2"))
    expect_true(same(logLik(f), logLik(l), tol), info = paste(info, "logLik"))
}

formulas <- list(
    list(fo = y ~ xd + f1 + f2, method = "patterns"),
    list(fo = y ~ xd * f2 + lg, method = "patterns"),
    list(fo = y ~ poly(xd, 2) + f1, method = "patterns"),
    list(fo = y ~ xd7 + f1, method = "patterns"),
    list(fo = y ~ f1 + f3 + xd, method = "patterns"),
    list(fo = y ~ xc * f1 + f2, method = "moments"),
    list(fo = y ~ xs + f1, method = "moments"),
    list(fo = y ~ xs * f2 + f1, method = "moments"),
    list(fo = y ~ xt + f1, method = "moments"),
    list(fo = y ~ xc + f1 + f3, method = "moments"),
    list(fo = y ~ log(xd) + I(xc^2) + f1, method = "moments"))
## ("patterns" designs are compressed with compress = TRUE: under "auto" an
## lm compresses only at most 0.1% patterns, a weighted one only with many
## numeric columns; "moments" designs are not compressed under "auto")
for (spec in formulas) {
    fo <- spec$fo
    cp <- if (spec$method == "patterns") TRUE else "auto"
    for (wt in c(FALSE, TRUE)) {
        l <- if (wt) lm(fo, dd, weights = w) else lm(fo, dd)
        for (engine in "auto") {
            info <- paste(deparse(fo), engine, if (wt) "weighted" else "")
            f <- suppressMessages(if (wt) {
                design_lm(fo, con, table = "d", engine = engine, weights = w, compress = cp)
            } else {
                design_lm(fo, con, table = "d", engine = engine, compress = cp)
            })
            expect_inherits(f, "design_lm")
            expect_identical(f$compressed, spec$method == "patterns", info = info)
            check_fit(f, l, info)
        }
    }
}

## weights as a string; compression forced off and on
l <- lm(y ~ xd + f1, dd, weights = w)
check_fit(design_lm(y ~ xd + f1, con, table = "d", weights = "w", compress = FALSE), l, "no compress")
expect_false(design_lm(y ~ xd + f1, con, table = "d", compress = FALSE)$compressed)
## (every row is a distinct pattern: nothing to compress)
expect_false(design_lm(y ~ xc + f1, con, table = "d", compress = TRUE)$compressed)
## with many patterns, "auto" gives up and TRUE compresses
op <- options(design.matrix.remote_compress_max = 10)
expect_false(design_lm(y ~ xd * f1, con, table = "d")$compressed)
expect_true(design_lm(y ~ xd * f1, con, table = "d", compress = TRUE)$compressed)
options(op)
check_fit(design_lm(y ~ xc + f1, con, table = "d", compress = TRUE), lm(y ~ xc + f1, dd), "compress = TRUE")

## the "auto" decision: the queries it runs, and its outcome
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
has <- function(q, pattern) any(grepl(pattern, q$sql, fixed = TRUE))
## a continuous DOUBLE key: skipped after a sample, without the estimate
q <- queries(design_lm(y ~ xc + f1, con, table = "d"))
expect_false(q$value$compressed)
expect_true(has(q, "TABLESAMPLE"))
expect_false(has(q, "approx_count_distinct(hash("))
check_fit(q$value, lm(y ~ xc + f1, dd), "auto: continuous")
## a discrete DOUBLE key: sample, estimate, and compression of a glm (about
## 25 patterns, 0.8% of the rows, at most 20%), not of an lm (at most 0.1%)
q <- queries(design_glm(y ~ xd7 + f1, gaussian, con, table = "d"))
expect_true(q$value$compressed)
expect_true(has(q, "TABLESAMPLE") && has(q, "approx_count_distinct(hash(") && has(q, "var_pop(__dm_y1)"))
q <- queries(design_lm(y ~ xd7 + f1, con, table = "d"))
expect_false(q$value$compressed)
expect_true(has(q, "TABLESAMPLE") && has(q, "approx_count_distinct(hash("))
expect_false(has(q, "var_pop(__dm_y1)"))
## integer keys: no sample, the estimate decides
q <- queries(design_glm(y ~ xd + f1, gaussian, con, table = "d"))
expect_true(q$value$compressed)
expect_false(has(q, "TABLESAMPLE"))
q <- queries(design_lm(y ~ id + f1, con, table = "d"))  # an integer id: one pattern per row
expect_false(q$value$compressed)
expect_true(has(q, "approx_count_distinct(hash("))
expect_false(has(q, "var_pop(__dm_y1)"))
check_fit(q$value, lm(y ~ id + f1, dd), "auto: integer id")
## the threshold on the estimate
q <- queries(design_lm(y ~ xd + f1, con, table = "d"), design.matrix.remote_compress_ratio = 1e-6)
expect_false(q$value$compressed)
expect_false(has(q, "var_pop(__dm_y1)"))
## the safety net after the GROUP BY: more than n / 2 patterns
op <- options(design.matrix.remote_compress_ratio = 1, design.matrix.remote_compress_max = 1e9)
D <- design_matrix(y ~ id + f1, con, table = "d")
expect_null(design.matrix:::dm_sql_patterns(D, "CAST(s.y AS DOUBLE)"))
options(op)

## printing
f <- design_lm(y ~ xd + f1, con, table = "d", compress = TRUE)
expect_stdout(print(f), "design patterns of 2,9")
expect_stdout(print(design_lm(y ~ xc + f1, con, table = "d")), "from SQL moments of")
expect_stdout(print(summary(f)), "Residual standard error")

## Parquet and sql= sources
pq <- tempfile(fileext = ".parquet")
DBI::dbExecute(con, sprintf("COPY (SELECT * FROM d) TO '%s' (FORMAT PARQUET)", pq))
dp <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s')", pq))
check_fit(design_lm(y ~ xd + f1, pq), lm(y ~ xd + f1, dp), "parquet")
check_fit(design_lm(y ~ xd + f1, con, sql = "SELECT * FROM d WHERE id <= 2000"),
          lm(y ~ xd + f1, dd[dd$id <= 2000, ]), "sql")

## design_lm() on a remote design
D <- design_matrix(y ~ xd + f1, con, table = "d")
check_fit(design_lm(D), lm(y ~ xd + f1, dd), "remote design")
check_fit(design_lm(D, weights = "w"), lm(y ~ xd + f1, dd, weights = w), "remote design, weighted")
D2 <- design_matrix(~ xd + f1, con, table = "d", na.action = stats::na.omit)
expect_error(design_lm(D2), "give `y`")
expect_error(design_lm(D, data = con), "With a design matrix")
expect_error(design_lm(y ~ xd, con, table = "d", y = "y"), "`y` applies only")
dk <- dd[stats::complete.cases(dd[c("xd", "f1")]), ]
check_fit(design_lm(D2, y = "w"), lm(w ~ xd + f1, dk), "remote design, y column")
## a response column with missing values: lm.fit()'s error, as in memory
expect_error(design_lm(D2, y = "xc"), "NA/NaN/Inf in 'y'")

## errors as lm()
DBI::dbWriteTable(con, "negw", transform(dd, w = -w))
expect_error(design_lm(y ~ xd, con, table = "negw", weights = "w"), "negative weights")

## predictions written to the database
f <- design_lm(y ~ xd + f1 + f2, con, table = "d")
l <- lm(y ~ xd + f1 + f2, dd)
p <- as.data.frame(predict(f, newdata = con, table = "d"))
p <- p[order(p$id), ]
ref <- predict(l, newdata = dd)
expect_true(same(p$.fitted, ref))
## an in-memory fit scored in the database
fm <- design_lm(y ~ xd + f1 + f2, dd)
p2 <- as.data.frame(predict(fm, newdata = con, table = "d"))
expect_true(same(p2$.fitted[order(p2$id)], ref))

## compressed products agree with the uncompressed ones
for (engine in "auto") {
    D <- suppressMessages(design_matrix(y ~ poly(xd, 2) * f2 + f1, con, table = "d", engine = engine))
    m1 <- design.matrix:::dm_remote_moments(D, ys = "s.y", w = "s.w")
    m0 <- no_compress(design.matrix:::dm_remote_moments(D, ys = "s.y", w = "s.w"))
    expect_true(same(m1$xtx, m0$xtx) && same(m1$xty, m0$xty) && same(m1$yty, m0$yty), info = engine)
    expect_equal(m1$n, m0$n, info = engine)
    expect_true(same(m1$sumw, m0$sumw), info = engine)
    expect_true(same(crossprod(D, "w" * D), no_compress(crossprod(D, "w" * D))), info = engine)
    expect_true(same(crossprod(D), no_compress(crossprod(D))), info = engine)
    expect_true(same(crossprod(D, "xc"), no_compress(crossprod(D, "xc"))), info = engine)
    expect_true(same(colSums("w" * D, group = "cl"), no_compress(colSums("w" * D, group = "cl"))),
                info = engine)
    expect_true(same(colMeans(D, group = "cl"), no_compress(colMeans(D, group = "cl"))), info = engine)
    expect_true(same(colMeans(D), no_compress(colMeans(D))), info = engine)
}
## a response or weights column with missing values: no compression, NA
D <- design_matrix(y ~ xd + f1, con, table = "d")
expect_true(all(is.na(crossprod(D, "xc"))))
## na.pass designs are not compressed (NA where the dense computation has it)
spec <- design_matrix(~ xd + f1, dd[1:100, ])
Dp <- design_matrix(spec, con, table = "d")
M <- design_matrix(spec, dd)
expect_true(same(crossprod(Dp), crossprod(M)))
expect_identical(is.na(crossprod(Dp)), is.na(crossprod(M)))

## residuals of a remote fit: sandwich covariances, as in memory
dk <- dd[stats::complete.cases(dd[c("y", "xc", "f1", "cl")]), ]
l <- lm(y ~ xc + f1, dk)
X <- model.matrix(l)
e <- unname(residuals(l))
for (engine in "auto") {
    f <- suppressMessages(design_lm(y ~ xc + f1, con, table = "d", engine = engine))
    R <- residuals(f)
    expect_inherits(R, "dm_remote")
    expect_true(same(crossprod(R, DBI::SQL('".resid" * ".resid"') * R), crossprod(X, e^2 * X)), info = engine)
    expect_true(same(colSums(".resid" * R, group = "cl"), rowsum(e * X, dk$cl)), info = engine)
}

test_con_done(con)
unlink(pq)

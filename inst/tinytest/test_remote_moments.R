## Compressed remote products: agreement with the uncompressed computation
## (the fits on remote data are tested in xtx.fit).
if (!requireNamespace("duckdb", quietly = TRUE) || !requireNamespace("DBI", quietly = TRUE))
    exit_file("duckdb is not installed")
source("helpers.R")
library(xtx)

con <- test_con()
same <- function(a, b, tol = 1e-7) isTRUE(all.equal(unname(a), unname(b), tolerance = tol))
no_compress <- function(expr) {
    op <- options(xtx.remote_compress = FALSE)
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

## the safety net after the GROUP BY: more than n / 2 patterns
op <- options(xtx.remote_compress_ratio = 1, xtx.remote_compress_max = 1e9)
D <- design_matrix(y ~ id + f1, con, table = "d")
expect_null(xtx:::xtx_sql_patterns(D, "CAST(s.y AS DOUBLE)"))
options(op)

## compressed products agree with the uncompressed ones
for (engine in "auto") {
    D <- suppressMessages(design_matrix(y ~ poly(xd, 2) * f2 + f1, con, table = "d", engine = engine))
    m1 <- xtx:::xtx_remote_moments(D, ys = "s.y", w = "s.w")
    m0 <- no_compress(xtx:::xtx_remote_moments(D, ys = "s.y", w = "s.w"))
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

test_con_done(con)

## Lookup columns (R/primitives_lookup.R): dmc_lookup(), dmc_cbind() and
## dmc_const(), in memory against base R, and remotely against memory.
library(design.matrix)
dmc_lookup <- design.matrix:::dmc_lookup
dmc_cbind <- design.matrix:::dmc_cbind
dmc_const <- design.matrix:::dmc_const
dmp_column <- design.matrix:::dmp_column
dmp_collect <- design.matrix:::dmp_collect
dmp_mult <- design.matrix:::dmp_mult
dmp_groupsum <- design.matrix:::dmp_groupsum
dm_colsum <- design.matrix:::dm_colsum

set.seed(11)
n <- 200
d <- data.frame(code = sample(0:6, n, TRUE), a = sample(c("u", "v", "w"), n, TRUE),
                b = sample(1:4, n, TRUE), x = rnorm(n), y = rnorm(n))
d$a[3] <- NA
vals <- c(10, 20, 30, 40, 50)          # codes 1..5; 0 and 6 have no value
lv <- data.frame(a = c("u", "u", "v", "w"), b = c(1L, 2L, 1L, 4L), stringsAsFactors = FALSE)
lvals <- c(1.5, 2.5, 3.5, 4.5)

## in memory
ref_codes <- ifelse(d$code >= 1 & d$code <= 5, vals[pmax(d$code, 1)], NA)
expect_equal(dmc_lookup(vals, d$code), ref_codes)
key <- paste(d$a, d$b)
ref_lv <- lvals[match(key, paste(lv$a, lv$b))]
ref_lv[is.na(d$a)] <- NA
expect_equal(dmc_lookup(lvals, d[c("a", "b")], lv), ref_lv)
dd <- d
dd$a <- factor(dd$a)
expect_equal(dmc_lookup(lvals, dd[c("a", "b")], lv), ref_lv)
expect_equal(dmc_lookup(c(7, 8), c("p", "q", "r"), data.frame(k = c("q", "p"))), c(8, 7, NA))
expect_equal(dmc_cbind(1:3, 4:6), cbind(1:3, 4:6))
expect_equal(dmc_const(1:4, 2), rep(2, 4))
expect_error(dmc_lookup(lvals, d["a"], lv), "numbers of columns")
## several columns of values: one lookup
vm <- cbind(vals, -vals)
expect_equal(dmc_lookup(vm, d$code), cbind(ref_codes, -ref_codes), check.attributes = FALSE)
expect_equal(dmc_lookup(cbind(lvals, 2 * lvals), d[c("a", "b")], lv), cbind(ref_lv, 2 * ref_lv),
             check.attributes = FALSE)
## dmp_materialize() and dmc_freeze(table = TRUE) are identities in memory
Dm <- design_matrix(~ x, d)
expect_identical(design.matrix:::dmp_materialize(Dm), Dm)
expect_identical(design.matrix:::dmc_freeze(d$x, table = TRUE), d$x)

if (requireNamespace("duckdb", quietly = TRUE)) {
    con <- DBI::dbConnect(duckdb::duckdb())
    d$.row <- seq_len(n)
    DBI::dbWriteTable(con, "d", d)
    D <- design_matrix(~ x, con, table = "d", na.action = na.pass)
    rc <- function(sql) dmp_column(D, sql)
    row <- rc("s.\".row\"")
    ## integer codes; values of rows in the database have no order: carry
    ## the row number through the same frame
    lk <- dmc_lookup(vals, rc("s.code"))
    got <- dmp_collect(dmc_cbind(row, lk))
    expect_equal(got[order(got[, 1]), 2], ref_codes)
    ## keys of several columns, character and integer
    lk2 <- dmc_lookup(lvals, dmc_cbind(rc("s.a"), rc("s.b")), lv)
    got <- dmp_collect(dmc_cbind(row, lk2))
    expect_equal(got[order(got[, 1]), 2], ref_lv)
    ## several columns of values
    lkm <- dmc_lookup(vm, rc("s.code"))
    got <- dmp_collect(dmc_cbind(row, lkm))
    expect_equal(unname(got[order(got[, 1]), 2:3]), unname(cbind(ref_codes, -ref_codes)))
    lkm2 <- dmc_lookup(cbind(lvals, 2 * lvals), dmc_cbind(rc("s.a"), rc("s.b")), lv)
    got <- dmp_collect(dmc_cbind(row, lkm2))
    expect_equal(unname(got[order(got[, 1]), 2:3]), unname(cbind(ref_lv, 2 * ref_lv)))
    ## a design and its columns over a table of its rows: the same values
    Ds <- design.matrix:::dmp_subset(D, DBI::SQL("s.code > 2"))
    Dt <- design.matrix:::dmp_materialize(Ds)
    expect_true(grepl("^SELECT \\* FROM", Dt$src$select))
    expect_equal(design.matrix:::dm_nrows(Dt), sum(d$code > 2))
    expect_equal(crossprod(Dt), crossprod(Ds))
    xt <- design.matrix:::dmc_on(rc("s.x * 2"), Dt)
    expect_equal(dm_colsum(xt), 2 * sum(d$x[d$code > 2]))
    ## a frozen column as a table: same sums, rows of the design only
    fz <- design.matrix:::dmc_freeze(design.matrix:::dmc_on(rc("s.x + 1"), Ds), table = TRUE)
    expect_equal(dm_colsum(fz), sum(d$x[d$code > 2] + 1))
    expect_true(grepl("^SELECT \\* FROM", unclass(fz)$frame$src$select))
    ## chained on a product: the lookup keeps the columns of the later frame
    b <- c(0.5, 2)
    xb <- dmp_mult(D, b)
    f <- xb + dmc_lookup(vals, rc("s.code"), along = list(xb))
    got <- dmp_collect(dmc_cbind(row, f))
    expect_equal(got[order(got[, 1]), 2], 0.5 + 2 * d$x + ref_codes)
    ## sums of lookup columns: NA as sum() gives it, grouped sums
    expect_true(is.na(dm_colsum(lk)))
    ok <- dmc_lookup(c(1, 1, 1, 1, 1, 1, 1), rc("s.code + 1"))
    expect_equal(dm_colsum(ok * rc("s.x")), sum(d$x))
    gs <- dmp_groupsum(ok * rc("s.y"), rc("s.b"))
    ref <- rowsum(d$y, d$b)
    expect_equal(unname(gs[, 1]), unname(ref[, 1]))
    expect_equal(dmc_const(lk, 3) |> dm_colsum(), 3 * n)
    DBI::dbDisconnect(con)
}

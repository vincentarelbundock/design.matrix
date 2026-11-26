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
    ## integer codes stay integer remotely
    colname <- function(v) gsub("^\\(s\\.|\"|\\)$", "", unclass(v)$sql)
    ci <- dmc_lookup(1:7, rc("s.code + 1"))
    expect_equal(unname(unclass(ci)$frame$src$types[colname(ci)]), "integer")
    ## dmc_groupsums(): several coded groupings and a total in one query,
    ## as dmp_groupsum() and dm_colsum() (NA for a group with a missing value)
    dmc_codes <- design.matrix:::dmc_codes
    dmc_groupsums <- design.matrix:::dmc_groupsums
    g1 <- dmc_codes(rc("s.code + 1"), 7L)
    g2 <- dmc_codes(rc("s.b"), 4L)
    m1 <- dmc_codes(d$code + 1L, 7L)
    m2 <- dmc_codes(d$b, 4L)
    yy <- dmc_cbind(rc("s.y"), rc("s.x * s.y"))
    r <- dmc_groupsums(yy, list(g1, g2), total = rc("s.x"))
    m <- dmc_groupsums(cbind(d$y, d$x * d$y), list(m1, m2), total = d$x)
    expect_equal(r, m)
    expect_equal(m$groups[[2]], unname(rowsum(cbind(d$y, d$x * d$y), d$b)), check.attributes = FALSE)
    expect_equal(m$total, sum(d$x))
    ## without a total; missing values
    r <- dmc_groupsums(dmc_lookup(vals, rc("s.code")), list(g1, g2))
    m <- dmc_groupsums(ref_codes, list(m1, m2))
    expect_equal(r, m)
    expect_true(is.na(dmc_groupsums(lk, list(g1), total = lk)$total))
    ## coded crosstabs: the remote one as in memory
    dmp_crosstab <- design.matrix:::dmp_crosstab
    expect_equal(dmp_crosstab(g1, g2, rc("s.x")), dmp_crosstab(m1, m2, d$x))
    expect_equal(dmp_crosstab(g2, g1), dmp_crosstab(m2, m1))
    ## dmp_materialize(x, cols): the data columns and the internal columns
    ## that the design or `cols` read
    lkc <- dmc_lookup(vals, rc("s.code"), along = list(lk2))
    Dl <- design.matrix:::dm_rebase(D, unclass(lkc)$frame)
    Dk <- design.matrix:::dmp_materialize(Dl, cols = list(lkc))
    nm <- names(Dk$src$types)
    expect_true(all(names(d) %in% nm))
    expect_true(colname(lkc) %in% nm)
    expect_false(colname(lk2) %in% nm)
    got <- dmp_collect(dmc_cbind(design.matrix:::dmc_on(row, Dk), design.matrix:::dmc_on(lkc, Dk)))
    expect_equal(got[order(got[, 1]), 2], ref_codes)
    expect_equal(crossprod(Dk), crossprod(D))
    ## frame lineage: columns built one on the other combine on the
    ## descendant frame; columns on two branches of the same design (two
    ## lookups, products or freezes from one frame) cannot be combined
    sib <- "two different branches of the same remote design"
    expect_error(lk + lk2, sib)
    expect_error(dmc_cbind(lk, lk2), sib)
    lk3 <- dmc_lookup(vals, rc("s.code"), along = list(lk2))
    expect_equal(dmp_collect(lk2 + lk3), ref_lv + ref_codes)
    expect_equal(dmp_collect(lk3 + lk2 + rc("s.x")), ref_codes + ref_lv + d$x)
    p1 <- dmp_mult(D, c(0, 2), along = list(lk))
    expect_equal(dmp_collect(p1 - 2 * rc("s.x") + lk), ref_codes)
    expect_error(p1 + lk2, sib)
    expect_error(dmp_mult(D, c(0, 1), along = list(lk, lk2)), sib)
    f1 <- design.matrix:::dmc_freeze(rc("s.x") + 1)
    f2 <- design.matrix:::dmc_freeze(rc("s.x") + 2)
    expect_error(f1 + f2, sib)
    f12 <- design.matrix:::dmc_freeze(dmc_cbind(rc("s.x") + 1, rc("s.x") + 2))
    expect_equal(dmp_collect(design.matrix:::dmc_col(f12, 2) - design.matrix:::dmc_col(f12, 1)), rep(1, n))
    ## (frames without a lineage, from older objects: the deeper one)
    nopath <- function(v) {
        v <- unclass(v)
        v$frame$path <- NULL
        structure(v, class = "dm_rcols")
    }
    expect_equal(design.matrix:::dm_frame_later(unclass(nopath(lk))$frame, unclass(nopath(lk3))$frame)$depth,
                 unclass(lk3)$frame$depth)
    expect_true(is.null(design.matrix:::dm_frame_later(unclass(nopath(lk))$frame, unclass(lk2)$frame)$path))
    DBI::dbDisconnect(con)
}

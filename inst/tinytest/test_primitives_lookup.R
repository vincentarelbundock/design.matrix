## Lookup columns (R/primitives_lookup.R): xtx_c_lookup(), xtx_c_cbind() and
## xtx_c_const(), in memory against base R, and remotely against memory.
source("helpers.R")
library(xtx)
xtx_c_lookup <- xtx:::xtx_c_lookup
xtx_c_cbind <- xtx:::xtx_c_cbind
xtx_c_const <- xtx:::xtx_c_const
xtx_p_column <- xtx:::xtx_p_column
xtx_p_collect <- xtx:::xtx_p_collect
xtx_p_mult <- xtx:::xtx_p_mult
xtx_p_groupsum <- xtx:::xtx_p_groupsum
xtx_colsum <- xtx:::xtx_colsum

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
expect_equal(xtx_c_lookup(vals, d$code), ref_codes)
key <- paste(d$a, d$b)
ref_lv <- lvals[match(key, paste(lv$a, lv$b))]
ref_lv[is.na(d$a)] <- NA
expect_equal(xtx_c_lookup(lvals, d[c("a", "b")], lv), ref_lv)
dd <- d
dd$a <- factor(dd$a)
expect_equal(xtx_c_lookup(lvals, dd[c("a", "b")], lv), ref_lv)
expect_equal(xtx_c_lookup(c(7, 8), c("p", "q", "r"), data.frame(k = c("q", "p"))), c(8, 7, NA))
expect_equal(xtx_c_cbind(1:3, 4:6), cbind(1:3, 4:6))
expect_equal(xtx_c_const(1:4, 2), rep(2, 4))
expect_error(xtx_c_lookup(lvals, d["a"], lv), "numbers of columns")
## several columns of values: one lookup
vm <- cbind(vals, -vals)
expect_equal(xtx_c_lookup(vm, d$code), cbind(ref_codes, -ref_codes), check.attributes = FALSE)
expect_equal(xtx_c_lookup(cbind(lvals, 2 * lvals), d[c("a", "b")], lv), cbind(ref_lv, 2 * ref_lv),
             check.attributes = FALSE)
## xtx_p_materialize() and xtx_c_freeze(table = TRUE) are identities in memory
Dm <- design_matrix(~ x, d)
expect_identical(xtx:::xtx_p_materialize(Dm), Dm)
expect_identical(xtx:::xtx_c_freeze(d$x, table = TRUE), d$x)

if (requireNamespace("duckdb", quietly = TRUE)) {
    con <- test_con()
    d$.row <- seq_len(n)
    DBI::dbWriteTable(con, "d", d)
    D <- design_matrix(~ x, con, table = "d", na.action = na.pass)
    rc <- function(sql) xtx_p_column(D, sql)
    row <- rc("s.\".row\"")
    ## integer codes; values of rows in the database have no order: carry
    ## the row number through the same frame
    lk <- xtx_c_lookup(vals, rc("s.code"))
    got <- xtx_p_collect(xtx_c_cbind(row, lk))
    expect_equal(got[order(got[, 1]), 2], ref_codes)
    ## keys of several columns, character and integer
    lk2 <- xtx_c_lookup(lvals, xtx_c_cbind(rc("s.a"), rc("s.b")), lv)
    got <- xtx_p_collect(xtx_c_cbind(row, lk2))
    expect_equal(got[order(got[, 1]), 2], ref_lv)
    ## several columns of values
    lkm <- xtx_c_lookup(vm, rc("s.code"))
    got <- xtx_p_collect(xtx_c_cbind(row, lkm))
    expect_equal(unname(got[order(got[, 1]), 2:3]), unname(cbind(ref_codes, -ref_codes)))
    lkm2 <- xtx_c_lookup(cbind(lvals, 2 * lvals), xtx_c_cbind(rc("s.a"), rc("s.b")), lv)
    got <- xtx_p_collect(xtx_c_cbind(row, lkm2))
    expect_equal(unname(got[order(got[, 1]), 2:3]), unname(cbind(ref_lv, 2 * ref_lv)))
    ## a design and its columns over a table of its rows: the same values
    Ds <- xtx:::xtx_p_subset(D, DBI::SQL("s.code > 2"))
    Dt <- xtx:::xtx_p_materialize(Ds)
    expect_true(grepl("^SELECT \\* FROM", Dt$src$select))
    expect_equal(xtx:::xtx_nrows(Dt), sum(d$code > 2))
    expect_equal(crossprod(Dt), crossprod(Ds))
    xt <- xtx:::xtx_c_on(rc("s.x * 2"), Dt)
    expect_equal(xtx_colsum(xt), 2 * sum(d$x[d$code > 2]))
    ## a frozen column as a table: same sums, rows of the design only
    fz <- xtx:::xtx_c_freeze(xtx:::xtx_c_on(rc("s.x + 1"), Ds), table = TRUE)
    expect_equal(xtx_colsum(fz), sum(d$x[d$code > 2] + 1))
    expect_true(grepl("^SELECT \\* FROM", unclass(fz)$frame$src$select))
    ## chained on a product: the lookup keeps the columns of the later frame
    ## (both pending: bound one on the other at the first use)
    b <- c(0.5, 2)
    xb <- xtx_p_mult(D, b)
    f <- xb + xtx_c_lookup(vals, rc("s.code"))
    got <- xtx_p_collect(xtx_c_cbind(row, f))
    expect_equal(got[order(got[, 1]), 2], 0.5 + 2 * d$x + ref_codes)
    ## sums of lookup columns: NA as sum() gives it, grouped sums
    expect_true(is.na(xtx_colsum(lk)))
    ok <- xtx_c_lookup(c(1, 1, 1, 1, 1, 1, 1), rc("s.code + 1"))
    expect_equal(xtx_colsum(ok * rc("s.x")), sum(d$x))
    gs <- xtx_p_groupsum(ok * rc("s.y"), rc("s.b"))
    ref <- rowsum(d$y, d$b)
    expect_equal(unname(gs[, 1]), unname(ref[, 1]))
    expect_equal(xtx_c_const(lk, 3) |> xtx_colsum(), 3 * n)
    ## integer codes stay integer remotely
    colname <- function(v) gsub("^\\(s\\.|\"|\\)$", "", unclass(v)$sql)
    ci <- xtx_c_lookup(1:7, rc("s.code + 1"))
    bound_frame <- function(v) xtx:::xtx_frame_bind(unclass(v)$frame)
    expect_equal(unname(bound_frame(ci)$src$types[colname(ci)]), "integer")
    ## xtx_c_groupsums(): several coded groupings and a total in one query,
    ## as xtx_p_groupsum() and xtx_colsum() (NA for a group with a missing value)
    xtx_c_codes <- xtx:::xtx_c_codes
    xtx_c_groupsums <- xtx:::xtx_c_groupsums
    g1 <- xtx_c_codes(rc("s.code + 1"), 7L)
    g2 <- xtx_c_codes(rc("s.b"), 4L)
    m1 <- xtx_c_codes(d$code + 1L, 7L)
    m2 <- xtx_c_codes(d$b, 4L)
    yy <- xtx_c_cbind(rc("s.y"), rc("s.x * s.y"))
    r <- xtx_c_groupsums(yy, list(g1, g2), total = rc("s.x"))
    m <- xtx_c_groupsums(cbind(d$y, d$x * d$y), list(m1, m2), total = d$x)
    expect_equal(r, m)
    expect_equal(m$groups[[2]], unname(rowsum(cbind(d$y, d$x * d$y), d$b)), check.attributes = FALSE)
    expect_equal(m$total, sum(d$x))
    ## without a total; missing values
    r <- xtx_c_groupsums(xtx_c_lookup(vals, rc("s.code")), list(g1, g2))
    m <- xtx_c_groupsums(ref_codes, list(m1, m2))
    expect_equal(r, m)
    expect_true(is.na(xtx_c_groupsums(lk, list(g1), total = lk)$total))
    ## coded crosstabs: the remote one as in memory
    xtx_p_crosstab <- xtx:::xtx_p_crosstab
    expect_equal(xtx_p_crosstab(g1, g2, rc("s.x")), xtx_p_crosstab(m1, m2, d$x))
    expect_equal(xtx_p_crosstab(g2, g1), xtx_p_crosstab(m2, m1))
    ## xtx_p_materialize(x, cols): the data columns and the internal columns
    ## that the design or `cols` read
    lkc <- xtx_c_lookup(vals, rc("s.code"))
    ## (lkc bound on top of lk2, already bound: the frame of their sum)
    Dl <- xtx:::xtx_rebase(D, unclass(lk2 + lkc)$frame)
    Dk <- xtx:::xtx_p_materialize(Dl, cols = list(lkc))
    nm <- names(Dk$src$types)
    expect_true(all(names(d) %in% nm))
    expect_true(colname(lkc) %in% nm)
    expect_false(colname(lk2) %in% nm)
    got <- xtx_p_collect(xtx_c_cbind(xtx:::xtx_c_on(row, Dk), xtx:::xtx_c_on(lkc, Dk)))
    expect_equal(got[order(got[, 1]), 2], ref_codes)
    expect_equal(crossprod(Dk), crossprod(D))
    ## frame lineage: columns built one on the other combine on the
    ## descendant frame; for columns settled on two branches of the same
    ## design (two lookups or products, each bound from one frame before
    ## they are combined: lk and lk2 were read above), the lookups and
    ## products of one branch are replayed on the other (epxj: no statement
    ## for a lookup, one view for a product); frozen columns on two branches
    ## cannot be combined
    sib <- "two different branches of the same remote design"
    nviews <- function() nrow(DBI::dbGetQuery(con, "SELECT * FROM duckdb_views() WHERE temporary AND NOT internal"))
    byrow <- function(v) {
        got <- xtx_p_collect(xtx_c_cbind(row, v))
        got[order(got[, 1]), -1]
    }
    v0 <- nviews()
    expect_equal(byrow(lk + lk2), ref_codes + ref_lv)
    expect_equal(unname(byrow(xtx_c_cbind(lk, lk2))), unname(cbind(ref_codes, ref_lv)))
    expect_equal(nviews(), v0)                   # lookups: no statement
    ## a lookup not read yet binds on top of the settled one
    lk3 <- xtx_c_lookup(vals, rc("s.code"))
    expect_equal(xtx_p_collect(lk2 + lk3), ref_lv + ref_codes)
    expect_equal(xtx_p_collect(lk3 + lk2 + rc("s.x")), ref_codes + ref_lv + d$x)
    p1 <- xtx_p_mult(D, c(0, 2))
    expect_equal(xtx_p_collect(p1 - 2 * rc("s.x") + lk), ref_codes)
    ## (p1 is now settled on top of lk: lk2's lookup is replayed on it)
    v1 <- nviews()
    expect_equal(byrow(p1 + lk2), 2 * d$x + ref_lv)
    expect_equal(nviews(), v1)
    expect_equal(byrow(xtx_p_mult(D, c(0, 1)) + lk + lk2), d$x + ref_codes + ref_lv)
    ## a product read on its own, combined with p1 (read on top of lk):
    ## one more view, the replayed product
    p5 <- xtx_p_mult(D, c(1, 0))
    expect_equal(xtx_colsum(p5), n)
    v2 <- nviews()
    expect_equal(byrow(p1 + p5 + lk), 2 * d$x + 1 + ref_codes)
    expect_equal(nviews(), v2 + 1L)
    ## the copies of a recipe add the columns of the original: a frame that
    ## has them (from an earlier replay) does not get them twice
    s15 <- p1 + p5
    expect_equal(byrow(s15 + p5), 2 * d$x + 2)
    f1 <- xtx:::xtx_c_freeze(rc("s.x") + 1)
    f2 <- xtx:::xtx_c_freeze(rc("s.x") + 2)
    expect_error(f1 + f2, sib)
    expect_error(f1 + f2, "stored \\(col_compute\\(\\)\\) or frozen")
    ## a frozen column and a lookup read on its own: the lookup is replayed
    expect_equal(byrow(f1 + lk2), d$x + 1 + ref_lv)
    f12 <- xtx:::xtx_c_freeze(xtx_c_cbind(rc("s.x") + 1, rc("s.x") + 2))
    expect_equal(xtx_p_collect(xtx:::xtx_c_col(f12, 2) - xtx:::xtx_c_col(f12, 1)), rep(1, n))
    ## (frames without a lineage, from older objects: the deeper one)
    nopath <- function(f) {
        f$path <- NULL
        f
    }
    expect_equal(xtx:::xtx_frame_later(nopath(bound_frame(lk)), nopath(bound_frame(lk3)))$depth,
                 bound_frame(lk3)$depth)
    expect_true(is.null(xtx:::xtx_frame_later(nopath(bound_frame(lk)), bound_frame(lk2))$path))

    ## ---- pending columns (nh9e): products and lookups bound at first use ----
    xtx_c_bind <- xtx:::xtx_c_bind
    xtx_c_freeze <- xtx:::xtx_c_freeze
    pending <- function(v) length(unclass(v)$frame$pend) > 0L
    ## a product is pending until read; arithmetic keeps it pending; it is
    ## bound once (later uses read the same view)
    pa <- xtx_p_mult(D, c(1, 1))
    expect_true(pending(pa))
    expect_true(pending(2 * pa + 1))
    nviews <- function() nrow(DBI::dbGetQuery(con, "SELECT * FROM duckdb_views() WHERE temporary AND NOT internal"))
    v0 <- nviews()
    expect_equal(xtx_colsum(pa), sum(1 + d$x))
    expect_equal(nviews(), v0 + 1L)
    expect_equal(xtx_colsum(pa * 2), 2 * sum(1 + d$x))
    expect_equal(nviews(), v0 + 1L)
    ## two products and a lookup, none read: one chain, in argument order
    pb <- xtx_p_mult(D, c(0, 3))
    pc <- xtx_p_mult(D, c(2, 0))
    lk4 <- xtx_c_lookup(vals, rc("s.code"))
    got <- xtx_p_collect(xtx_c_cbind(row, pb - pc + lk4))
    expect_equal(got[order(got[, 1]), 2], 3 * d$x - 2 + ref_codes)
    ## a lookup + product + freeze chain: the frozen column combines with the
    ## columns it was made from, and with a new product bound on top of it
    pe <- xtx_p_mult(D, c(0.5, -1))
    lk5 <- xtx_c_lookup(vals, rc("s.code"))
    fz5 <- xtx_c_freeze(pe * lk5)
    expect_false(pending(fz5))
    got <- xtx_p_collect(xtx_c_cbind(row, fz5 - pe * lk5))
    expect_equal(got[order(got[, 1]), 2], ifelse(is.na(ref_codes), NA_real_, 0))
    pf <- xtx_p_mult(D, c(1, 0))
    got <- xtx_p_collect(xtx_c_cbind(row, fz5 + pf + lk5))
    ref5 <- (0.5 - d$x) * ref_codes + 1 + ref_codes
    expect_equal(got[order(got[, 1]), 2], ref5)
    ## the same chain as a table (compute): a product made later binds on it
    fzt <- xtx_c_freeze(pe * lk5 + pf, table = TRUE)
    pg <- xtx_p_mult(D, c(0, 1))
    got <- xtx_p_collect(xtx_c_cbind(row, fzt + pg))
    expect_equal(got[order(got[, 1]), 2], (0.5 - d$x) * ref_codes + 1 + d$x)
    ## a lookup keyed by a pending lookup (codes of codes) binds after it
    lk6 <- xtx_c_lookup(c(5L, 4L, 3L, 2L, 1L), rc("s.code"))
    lk7 <- xtx_c_lookup(vals, lk6)
    got <- xtx_p_collect(xtx_c_cbind(row, lk7))
    ref7 <- ifelse(d$code >= 1 & d$code <= 5, vals[6L - pmin(pmax(d$code, 1), 5)], NA)
    expect_equal(got[order(got[, 1]), 2], ref7)
    ## xtx_c_bind(): columns used one by one, which still combine afterwards
    q1 <- xtx_c_lookup(vals, rc("s.code"))
    q2 <- xtx_p_mult(D, c(0, 1))
    qq <- xtx_c_bind(list(a = q1, b = list(q2, NULL)))
    expect_identical(names(qq), c("a", "b"))
    expect_false(pending(qq$a) || pending(qq$b[[1L]]))
    expect_true(is.null(qq$b[[2L]]))
    expect_true(is.na(xtx_colsum(qq$a)))
    expect_equal(xtx_colsum(qq$b[[1L]]), sum(d$x))
    got <- xtx_p_collect(xtx_c_cbind(row, qq$a + qq$b[[1L]]))
    expect_equal(got[order(got[, 1]), 2], ref_codes + d$x)
    expect_identical(xtx_c_bind(list(1:3, NULL)), list(1:3, NULL))
    ## without xtx_c_bind(), the same two columns read one by one sit on two
    ## branches: the lookup is replayed on the product's frame (no statement)
    r1 <- xtx_c_lookup(vals, rc("s.code"))
    r2 <- xtx_p_mult(D, c(0, 1))
    xtx_colsum(r2)
    expect_true(is.na(xtx_colsum(r1)))
    v3 <- nviews()
    got <- xtx_p_collect(xtx_c_cbind(row, r1 + r2))
    expect_equal(got[order(got[, 1]), 2], ref_codes + d$x)
    expect_equal(nviews(), v3)
    ## xtx_c_on() refuses a column whose recipes are not bound
    expect_error(xtx:::xtx_c_on(xtx_p_mult(D, c(1, 0)), D), "pending column")
    ## xtx_p_materialize(x, cols) binds the pending columns it is given
    lk8 <- xtx_c_lookup(vals, rc("s.code"))
    Dm8 <- xtx:::xtx_p_materialize(D, cols = list(lk8))
    expect_true(colname(lk8) %in% names(Dm8$src$types))
    got <- xtx_p_collect(xtx_c_cbind(xtx:::xtx_c_on(row, Dm8), xtx:::xtx_c_on(lk8, Dm8)))
    expect_equal(got[order(got[, 1]), 2], ref_codes)
    test_con_done(con)
}

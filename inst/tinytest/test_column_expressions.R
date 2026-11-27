## Remote columns as expressions that do not depend on what was read
## before: the rules for combining sources, storage, failure and retry, and
## the per-column flags (logical-valued, type). Per-row values are compared
## after matching on a key read in the same query (`id`); the data has
## repeated predictor values, and one fixture has no unique key at all.
source("helpers.R")
library(xtx)
if (!requireNamespace("duckdb", quietly = TRUE) || !requireNamespace("DBI", quietly = TRUE))
    exit_file("duckdb is not installed")

same <- function(a, b, tol = 1e-10) isTRUE(all.equal(as.vector(a), as.vector(b), tolerance = tol))
set.seed(7)
n <- 240
## (predictor values repeat: 6 distinct rows of predictors)
d <- data.frame(id = sample(n), x = rep(c(-1, 0, 2), 80), f = rep(c("a", "b"), each = 120),
                z = c(NA, NaN, rnorm(n - 2)), s = rep(c("u", "v", NA, "w"), 60), i = rep(1:5, 48),
                g = rep(c("k1", "k2", "k3", "k4"), 60), stringsAsFactors = FALSE)
d$y <- d$x + rnorm(n)
d <- d[order(d$id), ]
M <- design_matrix(y ~ x * f, d)
p <- ncol(M)
b1 <- rnorm(p)
b2 <- rnorm(p)

con <- test_con()
DBI::dbWriteTable(con, "d", d)
D <- design_matrix(y ~ x * f, con, table = "d")
by_id <- function(col) {
    r <- as.data.frame(col, keep = "id")
    r <- r[order(r$id), , drop = FALSE]
    out <- as.matrix(r[setdiff(names(r), "id")])
    dimnames(out) <- NULL
    if (ncol(out) == 1L) out[, 1L] else out
}
objects_now <- function() {
    list(views = DBI::dbGetQuery(con, "SELECT view_name FROM duckdb_views() WHERE NOT internal ORDER BY 1")$view_name,
         tables = DBI::dbGetQuery(con, "SELECT table_name FROM duckdb_tables() ORDER BY 1")$table_name)
}

## ---- products read in either order ----
## p, q, sum(p), collect q, collect p - q with a key; and the reverse
m1 <- drop(M %*% b1)
m2 <- drop(M %*% b2)
run <- function(order) {
    pp <- D %*% b1
    qq <- D %*% b2
    steps <- list(sum_p = function() sum(pp), q = function() by_id(qq), diff = function() by_id(pp - qq))
    out <- lapply(steps[order], function(f) f())
    names(out) <- order
    out
}
o0 <- objects_now()
fw <- run(c("sum_p", "q", "diff"))
bw <- run(c("diff", "q", "sum_p"))
expect_true(same(fw$diff, m1 - m2))
expect_true(same(fw$q, m2))
expect_true(same(fw$sum_p, sum(m1)))
expect_identical(fw$diff, bw$diff)
expect_identical(fw$q, bw$q)
expect_identical(fw$sum_p, bw$sum_p)
## reads created nothing on the connection
expect_identical(objects_now(), o0)
## a read is one query, and leaves the columns as they were
pp <- D %*% b1
qq <- D %*% b2
before <- list(pp, qq)
rq <- record_queries({
    sum(pp)
    as.vector(qq)
    as.data.frame(pp - qq, keep = "id")
})
expect_equal(rq$count, 3L)
expect_true(all(rq$log$fn == "dbGetQuery"))
expect_identical(list(pp, qq), before)
## the final read compiles both products over one relation: each once
sql <- rq$log$sql[3L]
node <- function(col) names(xtx:::xtx_frame_nodes(unclass(col)$frame))
count <- function(pattern, x) lengths(regmatches(x, gregexpr(pattern, x)))
expect_equal(count(paste0(" AS \"?", node(pp), "_1"), sql), 1L)
expect_equal(count(paste0(" AS \"?", node(qq), "_1"), sql), 1L)
expect_equal(count("FROM d\\b", sql), 1L)

## ---- p - p, p - q, exceptional values ----
expect_true(all(by_id(pp - pp) == 0))
expect_true(same(by_id(pp - qq), m1 - m2))
## (p - p shares the node, and keeps the subtraction: Inf - Inf is NaN)
big <- pp * 1e308 * 1e308
v <- by_id(big - big)
expect_true(all(is.nan(v) | v == 0))
expect_true(any(is.nan(v)))
sql <- record_queries(as.vector(pp - pp))$log$sql
expect_equal(count(paste0(" AS \"?", node(pp), "_1"), sql), 1L)

## ---- nested lookups, a lookup keyed by a computed column ----
g <- design_groups(D, "g")
vals <- c(10, 20, 30, 40)
l1 <- group_lookup(vals, g)
ref1 <- vals[match(d$g, c("k1", "k2", "k3", "k4"))]
expect_true(same(by_id(l1), ref1))
## a key computed from a product: groups of the sign of p
sgn <- (pp > 0) + 1
g2 <- design_groups(D, list(sg = sgn))
l2 <- group_lookup(c(-5, 5), g2)
ref2 <- c(-5, 5)[(m1 > 0) + 1]
expect_true(same(by_id(l2), ref2))
## nested: a lookup keyed by a lookup
l3 <- xtx:::xtx_c_lookup(c(100, 200, 300, 400), xtx:::xtx_c_lookup(4:1, g$codes))
expect_true(same(by_id(l3), c(400, 300, 200, 100)[match(d$g, c("k1", "k2", "k3", "k4"))]))
## all of them with products, in either order of reading
e1 <- l1 + l2 + l3 + pp - qq
expect_true(same(sum(l2), sum(ref2)))
expect_true(same(by_id(e1), ref1 + ref2 + by_id(l3) + m1 - m2))
## lookups neither multiply nor drop rows
expect_equal(length(as.vector(e1)), n)
expect_equal(nrow(as.data.frame(l1 + l3, keep = "id")), n)

## ---- cbind() and slicing ----
P <- D %*% cbind(b1, b2)
cb <- cbind(P[, 2], pp, l1)
expect_true(same(by_id(cb), cbind(m2, m1, ref1)))
expect_true(same(by_id(cb[, c(3, 1)]), cbind(ref1, m2)))
expect_true(same(by_id(P[, 1] - pp), rep(0, n)))
## a slice reads only the nodes it names
sql <- record_queries(as.vector(cb[, 3]))$log$sql
expect_false(grepl(node(pp), sql, fixed = TRUE))

## ---- two projections of one source combine (were refused before) ----
f1 <- col_compute(pp + 1, table = FALSE)
f2 <- col_compute(qq * 2, table = FALSE)
o1 <- objects_now()
expect_true(same(by_id(f1 + f2), m1 + 1 + 2 * m2))
expect_true(same(by_id(f2 - f1), 2 * m2 - m1 - 1))
expect_true(same(sum(f1), sum(m1 + 1)))
## a projection of projections, and with the columns they were made from
f3 <- col_compute(f1 * f2, table = FALSE)
expect_true(same(by_id(f3 - pp), (m1 + 1) * 2 * m2 - m1))
## projections are SQL only: no table, no view, no statement
expect_identical(objects_now(), o1)
expect_equal(record_queries(col_compute(pp * qq, table = FALSE))$count, 0L)
## a projection is evaluated once per row in a statement: a volatile
## expression projected gives one value per row wherever it is named
r <- col_compute(design_column(D, DBI::SQL("random()")), table = FALSE)
got <- as.matrix(cbind(r, r, r - r))
expect_identical(got[, 1], got[, 2])
expect_true(all(got[, 3] == 0))
expect_true(length(unique(got[, 1])) > 1L)

## ---- a stored ancestor and a lazy descendant ----
st <- col_compute(pp)
t1 <- objects_now()
expect_true(same(by_id(st + qq), m1 + m2))                # an original expression on the stored relation
expect_true(same(by_id(qq - st), m2 - m1))
expect_true(same(by_id(st - pp), rep(0, n)))              # the product it was made from
expect_true(same(by_id(st + l1 + f2), m1 + ref1 + 2 * m2))
expect_true(same(by_id(D %*% b2 + st), m1 + m2))          # a product made after the store
expect_identical(objects_now(), t1)
## the stored value is read from the table, not computed again
sql <- record_queries(as.vector(st))$log$sql
expect_false(grepl("FROM d\\b", sql))
## a store on a store: the deepest relation holds both values
st_b <- col_compute(st * 2 + qq)
expect_true(same(by_id(st_b - st), m1 + m2))
expect_true(same(by_id(st_b + st + pp), 4 * m1 + m2))
## a store captures its values: the source can change afterwards
DBI::dbWriteTable(con, "dm", d)
Dm <- design_matrix(y ~ x * f, con, table = "dm")
lazy <- Dm %*% b1
kept <- col_compute(lazy)
DBI::dbExecute(con, "UPDATE dm SET x = x + 1")
expect_true(same(by_id(kept), m1))
expect_false(same(by_id(lazy), m1))

## ---- two outputs stored together ----
st12 <- col_compute(cbind(pp, qq))
expect_true(same(by_id(st12[, 1] - st12[, 2]), m1 - m2))
expect_true(same(by_id(st12[, 2] + pp), m2 + m1))

## ---- two sibling stores: an error, before any query ----
sa <- col_compute(D %*% b1)
sb <- col_compute(D %*% b2)
rq <- record_queries(msg <- tryCatch(sa + sb, error = function(e) conditionMessage(e)))
expect_true(grepl("two different branches", msg))
expect_true(grepl("col_compute", msg))
expect_equal(rq$count, 0L)
expect_error(cbind(sa, sb), "two different branches")
expect_error(sum(sa, sb), "two different branches")
expect_error(crossprod(sa * D, sb * D), "two different branches")
## projections on sibling stores do not combine either
expect_error(col_compute(sa + 1, table = FALSE) + col_compute(sb + 1, table = FALSE), "two different branches")
## each sibling still works, alone and with lazy columns
expect_true(same(by_id(sa + qq), m1 + m2))
expect_true(same(by_id(sb + pp), m1 + m2))

## ---- another design, subset or connection ----
Ds <- D[DBI::SQL("TRUE"), ]                               # every row kept: still another design
expect_equal(nrow(Ds), n)
expect_error(pp + Ds %*% b1, "different remote designs")
expect_error(pp + design_column(Ds, "y"), "different remote designs")
expect_error(crossprod(Ds, pp), "different remote designs")
expect_error(group_sums(Ds %*% b1, g), "another remote design")
D2 <- design_matrix(y ~ x * f, con, table = "d")          # same data and formula
expect_error(pp + D2 %*% b1, "different remote designs")
con2 <- DBI::dbConnect(duckdb::duckdb())
DBI::dbWriteTable(con2, "d", d)
D3 <- design_matrix(y ~ x * f, con2, table = "d")
expect_error(pp + D3 %*% b1, "different remote designs")
## a subset translates the columns it is given onto its rows, explicitly
keep <- design_column(D, "x") >= 0
sub <- design_subset(D, keep, list(p = pp, l = l1, st = st, f = f1))
k <- d$x >= 0
by_id_sub <- function(col) {
    r <- as.data.frame(col, keep = "id")
    r[order(r$id), 2L]
}
expect_true(same(by_id_sub(sub$cols$p), m1[k]))
expect_true(same(by_id_sub(sub$cols$l + sub$cols$st + sub$cols$f), (ref1 + 2 * m1 + 1)[k]))
expect_true(same(by_id_sub(sub$cols$p - sub$x %*% b1), rep(0, sum(k))))
expect_error(sub$cols$p + pp, "different remote designs")
## and refuses the columns of a design it was not made from
expect_error(design_subset(sub$x, design_column(sub$x, "x") > 0, list(D2 %*% b1)), "another remote design")
## design_subset() without `keep`: the design and the columns as they are
al <- design_subset(D, cols = list(p = pp, st = st))
expect_identical(al$x, D)
expect_true(same(by_id(al$cols$p + al$cols$st), 2 * m1))

## ---- a disconnected handle ----
p3 <- D3 %*% b1
DBI::dbDisconnect(con2, shutdown = TRUE)
expect_error(as.vector(p3), "connection of this remote column is closed")
expect_error(col_compute(p3), "connection of this remote column is closed")
expect_true(any(grepl("closed", format(p3))))

## ---- duplicate observations, no unique key ----
dd <- d[rep(c(1, 1, 2, 3, 3, 3), 20), c("x", "f", "y", "g")]
DBI::dbWriteTable(con, "dd", dd)
Dd <- design_matrix(y ~ x * f, con, table = "dd")
Md <- design_matrix(y ~ x * f, dd)
pd <- Dd %*% b1
qd <- Dd %*% b2
expect_true(same(sum(pd), sum(Md %*% b1)))
expect_equal(length(as.vector(qd)), nrow(dd))
## (rows match on the values read in the same query)
got <- as.data.frame(cbind(pd - qd, col_compute(pd) + qd), keep = c("x", "f", "y"), nm = c("a", "b"))
ref <- cbind(dd[c("x", "f", "y")], a = drop(Md %*% (b1 - b2)), b = drop(Md %*% (b1 + b2)))
ord <- function(t) t[do.call(order, unname(as.list(t))), ]
expect_equal(ord(got), ord(ref), check.attributes = FALSE)
gd <- design_groups(Dd, "g")
expect_equal(nrow(as.data.frame(group_lookup(seq_len(gd$G), gd) + pd)), nrow(dd))

## ---- failure and retry ----
## a query error after staging inputs (the lookup's table of values): the
## columns still work, and only that table was added
t0 <- objects_now()
lk <- group_lookup(vals, g)
t1 <- objects_now()
expect_equal(length(setdiff(t1$tables, t0$tables)), 1L)
bad <- lk + design_column(D, DBI::SQL("CAST(s AS DOUBLE)"))
expect_error(as.vector(bad), "Could not convert")
expect_error(sum(bad), "Could not convert")
expect_identical(objects_now(), t1)
expect_true(same(by_id(lk + pp), ref1 + m1))               # retry with the same handles
## a compilation error (an invalid expression) is reported when the column
## is made, and nothing is left behind
expect_error(design_column(D, DBI::SQL("nope + 1")), "invalid SQL")
expect_identical(objects_now(), t1)
## a failed store: no table, no handle; the inputs still work
expect_error(col_compute(bad), "Could not convert")
expect_identical(objects_now(), t1)
expect_true(same(by_id(col_compute(lk + pp)), ref1 + m1))
t2 <- objects_now()
## a failed named store leaves an existing table of that name as it was
DBI::dbExecute(con, "CREATE TABLE keepme AS SELECT 42 AS v")
t3 <- objects_now()
expect_error(col_compute(bad, name = "keepme", overwrite = TRUE), "Could not convert")
expect_equal(DBI::dbGetQuery(con, "SELECT v FROM keepme")$v, 42)
expect_identical(objects_now(), t3)
## a name collision with overwrite off: an error before any statement that writes
rq <- record_queries(msg <- tryCatch(col_compute(pp, name = "keepme"), error = function(e) conditionMessage(e)))
expect_true(grepl("already exists", msg))
expect_false(any(grepl("CREATE", rq$log$sql)))
expect_equal(DBI::dbGetQuery(con, "SELECT v FROM keepme")$v, 42)
## names are checked before writing
expect_error(col_compute(pp, name = "t_new", nm = "x"), "already has a column")
expect_error(col_compute(pp, name = "t_new", nm = "__xtx_a"), "reserved")
expect_false("t_new" %in% objects_now()$tables)
## retry after the failure: overwrite, then the handle works with its design
pk <- col_compute(pp, name = "keepme", overwrite = TRUE, nm = "fit")
expect_true(same(DBI::dbGetQuery(con, "SELECT fit FROM keepme ORDER BY id")$fit, m1))
expect_true(same(by_id(pk + qq), m1 + m2))
## the connection is still usable and the handles made before are valid
expect_true(DBI::dbIsValid(con))
expect_true(same(by_id(pp - qq), m1 - m2))

## ---- two live designs on one connection ----
Da <- design_matrix(y ~ x * f, con, table = "d")
Db <- design_matrix(y ~ x, con, table = "d")
bb <- c(1, -2)
pa <- Da %*% b1
pb <- Db %*% bb
sta <- col_compute(pa + 1)
stb <- col_compute(pb + 1)
la <- group_lookup(vals, design_groups(Da, "g"))
lb <- group_lookup(-vals, design_groups(Db, "g"))
mb <- 1 - 2 * d$x
expect_true(same(by_id(sta + pa + la), 2 * m1 + 1 + ref1))
expect_true(same(by_id(stb + pb + lb), 2 * mb + 1 - ref1))
expect_true(same(by_id(sta - 1), m1))
## (distinct internal names: the nodes and tables of the two never meet)
expect_false(any(node(pa) %in% node(pb)))
expect_false(identical(xtx:::xtx_frame_root(unclass(sta)$frame)$select, xtx:::xtx_frame_root(unclass(stb)$frame)$select))
expect_error(sta + stb, "different remote designs")

## ---- the response of a QR is known by its identity, not its SQL text ----
y1 <- design_column(D, "y")
q <- qr(D)
qy <- xtx:::xtx_qr_for(D, 1e-7, y1)
expect_equal(xtx:::xtx_qr_xty(qy, D, y1), crossprod(D, y1), check.attributes = FALSE)
## the same SQL text over another stored relation is another response
DBI::dbWriteTable(con, "dz", transform(d, v = y))
Dz <- design_matrix(y ~ x * f, con, table = "dz")
vz <- design_column(Dz, "v")
qz <- xtx:::xtx_qr_for(Dz, 1e-7, vz)
out <- design_output(Dz, list(w = vz * 3), design = TRUE)   # a design whose data has another `v`? no: `w`
vw <- design_column(out, DBI::SQL("w"))
expect_false(identical(xtx:::xtx_column_key(vz), xtx:::xtx_column_key(vw)))
expect_identical(xtx:::xtx_column_key(vz), xtx:::xtx_column_key(design_column(Dz, "v")))
k1 <- xtx:::xtx_column_key(col_compute(vz, name = "zz1", nm = "val"))
k2 <- xtx:::xtx_column_key(col_compute(vz * 2, name = "zz2", nm = "val"))
expect_identical(k1$sql, k2$sql)                           # the same text
expect_false(identical(k1, k2))                            # two responses

## ---- flags per output column: logical values, types ----
x <- design_column(D, "x")
z <- design_column(D, "z")
lg <- cbind(x > 0, z > 0)
rx <- d$x > 0
rz <- d$z > 0
expect_equal(by_id(lg[, 1] + 1), rx + 1)
expect_equal(by_id(lg * 2), unname(cbind(rx, rz) * 2))
expect_equal(by_id(col_compute(x > 0) * 2), rx * 2)
expect_equal(by_id(col_compute(x > 0, table = FALSE) + 1), rx + 1)
expect_equal(by_id(design_subset(D, cols = list(b = x > 0))$cols$b + 1), rx + 1)
expect_equal(by_id_sub(design_subset(D, x > -1, cols = list(b = x > 0))$cols$b + 1), (rx + 1)[d$x > -1])
expect_equal(by_id(xtx:::xtx_c_col(lg, 2) + 1), rz + 1)
expect_equal(by_id(col_compute(lg)[, 2] + col_compute(lg, table = FALSE)[, 1]), rz + rx)
st_l <- col_compute(x > 0, name = "lgl_t", nm = "pos")
expect_equal(by_id(st_l + 1), rx + 1)
expect_true(is.logical(DBI::dbGetQuery(con, "SELECT pos FROM lgl_t")$pos))
## !!x of a number is x != 0 (NA stays NA); !is.na(x) keeps its short SQL
expect_equal(by_id(!x), as.numeric(!d$x))
expect_equal(by_id(!!x), as.numeric(!!d$x))
expect_equal(by_id((!!x) + 1), (!!d$x) + 1)
expect_equal(by_id(!z), as.numeric(!d$z))
expect_equal(by_id(!!z), as.numeric(!!d$z))
expect_equal(by_id(!!(x > 0)), as.numeric(rx))
expect_identical(unclass(!is.na(z))$sql, unclass(xtx:::xtx_c_present(z))$sql)
## is.na() by the type of the column, through slicing and translation
key <- design_key(D, c("s", "f", "i", "z"))
ref_na <- unname(as.matrix(is.na(d[c("s", "f", "i", "z")])) * 1)
expect_equal(by_id(is.na(key)), ref_na)
expect_equal(by_id(is.na(key[, 1])), ref_na[, 1])
expect_equal(by_id(is.na(key)[, 1] + 1), ref_na[, 1] + 1)
expect_equal(by_id(is.na(cbind(key[, 4], key[, 1]))), ref_na[, c(4, 1)])
ks <- design_subset(D, x > -1, cols = list(key))$cols[[1L]]
expect_equal(by_id_sub(is.na(ks[, 1])), ref_na[d$x > -1, 1])
expect_equal(by_id(is.na(col_compute(key[, 1]))), ref_na[, 1])
expect_equal(by_id(is.na(col_compute(key, table = FALSE))[, 1]), ref_na[, 1])
## text keys compare as text
expect_equal(by_id(key[, 1] == "u"), as.numeric(d$s == "u"))
expect_equal(by_id(key[, 1] == key[, 1]), as.numeric(d$s == d$s))
## the factor response of a remote frame is text
DBI::dbWriteTable(con, "dr", transform(d, r = factor(ifelse(i > 2, "yes", "no"))))
fr <- design_frame(r ~ x, con, table = "dr")
expect_inherits(fr$response$y, "xtx_column")
expect_equal(sum(is.na(fr$response$y)), 0)
expect_equal(sum(fr$response$y == "yes"), sum(d$i > 2))
## logical columns are 0/1 numbers wherever R takes them as such
expect_equal(by_id(exp(x > 0)), exp(rx))
expect_equal(by_id(sqrt(x > 0) + log1p(x > 0) + abs(x > 0) + floor(x > 0) + round(x > 0)), 4 * rx + log1p(rx))
expect_equal(by_id(log(x >= 0)), log(d$x >= 0))
expect_equal(by_id(rowSums(lg)), rowSums(cbind(rx, rz)))
expect_equal(by_id((x > 0) > 0.5), as.numeric(rx > 0.5))
expect_equal(by_id((x > 0) == 1), as.numeric(rx == 1))
expect_equal(by_id(0.5 <= (z > 0)), as.numeric(0.5 <= rz))
expect_equal(by_id((x > 0) == (z > 0)), as.numeric(rx == rz))
expect_equal(by_id((x > 0) < x), as.numeric(rx < d$x))
expect_equal(by_id((x > 0) == TRUE), as.numeric(rx == TRUE))
expect_equal(by_id(col_pmax(x > 0, 0.5)), pmax(rx, 0.5))
expect_equal(by_id(col_pmin(z > 0, 0.5)), pmin(rz, 0.5))
expect_equal(by_id(col_pmax(x > 0, x)), pmax(rx, d$x))
expect_equal(by_id(col_ifelse(z > 0, x > 1, 0)), ifelse(rz, d$x > 1, 0))
expect_equal(by_id(col_ifelse(z > 0, x > 1, 0) + 1), ifelse(rz, d$x > 1, 0) + 1)
expect_equal(by_id(col_ifelse(z > 0, x > 1, x < 0) + 1), ifelse(rz, d$x > 1, d$x < 0) + 1)
expect_equal(by_id(col_ifelse(x, 1, 2)), ifelse(d$x != 0, 1, 2))
expect_equal(by_id((x > 0) & x), as.numeric(rx & d$x))
expect_equal(by_id(x | (z > 0)), as.numeric(d$x | rz))
expect_equal(crossprod((x > 0) * D), crossprod(as.numeric(rx) * M), check.attributes = FALSE)
expect_equal(crossprod(D, x > 0), crossprod(M, as.numeric(rx)), check.attributes = FALSE)
expect_equal(colSums((x > 0) * D), colSums(as.numeric(rx) * as.matrix(M)))
## (%% and %/% stay refused for remote columns)
expect_error(x %% 2, "not supported")
expect_error(x %/% 2, "not supported")

test_con_done(con)

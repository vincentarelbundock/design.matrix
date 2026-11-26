## The public column type (Kata design.matrix#bvxy): design_column(),
## products of a remote design returning columns, as.vector(), as.matrix(),
## col_compute(), as.data.frame(), print/format and x[, j]. The same chain of public calls
## gives the same values in memory and on DuckDB. Only exported functions.
source("helpers.R")
library(xtx)

same <- function(a, b, tol = 1e-8) isTRUE(all.equal(as.vector(a), as.vector(b), tolerance = tol))

set.seed(41)
n <- 300
d <- data.frame(id = seq_len(n), y = rnorm(n), x = rnorm(n), z = runif(n, 1, 2), w = runif(n, 0.5, 2),
                f = sample(c("b", "a", "c"), n, TRUE), flag = sample(c(TRUE, FALSE), n, TRUE))
fo <- y ~ x * f + log(z)
M <- design_matrix(fo, d)
p <- ncol(M)
b <- rnorm(p)
B <- cbind(b1 = b, b2 = -b / 2)
V <- crossprod(matrix(rnorm(p^2), p))

## ---- in memory: a column is a plain numeric vector ----
expect_identical(design_column(M, "y"), as.double(d$y))
expect_identical(design_column(M, "log(z)"), log(d$z))
expect_identical(design_column(2 * M, "y"), as.double(d$y))
expect_error(design_column(M, "w"), "not a variable of the design")
expect_error(design_column(M, "f"), "not a numeric vector")
expect_error(design_column(as.matrix(M), "y"), "must be a design")
v <- design_column(M, "y")
expect_identical(as.vector(v), v)
expect_identical(col_compute(v), v)
expect_identical(col_compute(v, name = "never_made"), v)
Pm <- M %*% B
expect_identical(as.matrix(Pm), Pm)

## the chains of public calls compared below (the same code for both backends)
chain <- function(D) {
    eta <- D %*% b
    r <- design_column(D, "y") - exp(eta / 10)
    list(eta = eta, r = r, r2 = col_compute(r * r))
}
chain_cols <- function(D) {
    P <- D %*% B
    P[, 1] - P[, 2]
}
chain_quad <- function(D) design_column(D, "y") / sqrt(quadratic_diagonal(D, V))

if (!requireNamespace("duckdb", quietly = TRUE) || !requireNamespace("DBI", quietly = TRUE))
    exit_file("duckdb is not installed")

expect_error(design_column(M, DBI::SQL("y")), "only of a remote design")

con <- test_con()
DBI::dbWriteTable(con, "d", d)
D <- design_matrix(fo, con, table = "d")

## values of a remote column in the order of the rows of d (by the key `id`)
by_id <- function(col) {
    r <- as.data.frame(col, keep = "id")
    r <- r[order(r$id), , drop = FALSE]
    expect_identical(r$id, d$id)
    as.matrix(r[setdiff(names(r), "id")])
}

## ---- the class, printing ----
pr <- D %*% b
expect_inherits(pr, "xtx_column")
## dim() is c(n, k), length() n * k, as for the matrix in memory
expect_equal(dim(pr), c(n, 1))
expect_equal(dim(pr), dim(M %*% b))
expect_equal(nrow(pr), n)
expect_equal(ncol(pr), 1)
expect_equal(NROW(pr), n)
expect_equal(NCOL(pr), 1)
expect_equal(length(pr), n)
expect_equal(length(pr), length(M %*% b))
out <- capture.output(print(pr))
expect_true(grepl("^<design_column: 300 x 1, on", out[1L]))
expect_true(grepl("duckdb", out[1L]))
expect_identical(format(pr), out)
expect_identical(names(as.data.frame(pr, keep = "id")), c("id", ".fitted"))
expect_identical(names(as.data.frame(pr)), ".fitted")
expect_identical(names(as.data.frame(pr, nm = "eta")), "eta")
expect_true(same(by_id(pr), M %*% b))
expect_true(same(sort(as.vector(pr)), sort(drop(M %*% b))))
expect_true(is.double(as.vector(pr)) && is.null(attributes(as.vector(pr))))
expect_identical(dim(as.matrix(pr)), c(as.integer(n), 1L))
expect_null(dimnames(as.matrix(pr)))
expect_true(same(sort(as.matrix(pr)[, 1]), sort(drop(M %*% b))))

## ---- design_column() ----
y <- design_column(D, "y")
expect_inherits(y, "xtx_column")
expect_identical(names(as.data.frame(y, keep = "id")), c("id", "y"))
expect_true(same(by_id(y), d$y))
expect_true(same(by_id(design_column(D, "flag")), as.numeric(d$flag)))
yw <- design_column(D, DBI::SQL("y * w"))
expect_identical(names(as.data.frame(yw)), ".value")
expect_true(same(by_id(yw), d$y * d$w))
expect_true(same(by_id(design_column("w" * D, "y")), d$y))  # not scaled
expect_error(design_column(D, "nope"), "no column")
expect_error(design_column(D, "f"), "not numeric")
expect_error(design_column(D, DBI::SQL("nope + 1")), "invalid SQL")
expect_error(as.data.frame(y, keep = "nope"), "no column")
expect_error(as.data.frame(y, keep = "y"), "name of the values")
expect_identical(names(as.data.frame(y, keep = "y", nm = "y2")), c("y", "y2"))

## ---- chains: the same values in memory and remotely ----
cm <- chain(M)
cr <- chain(D)
expect_inherits(cr$r2, "xtx_column")
expect_true(same(by_id(cr$eta), cm$eta))
expect_true(same(by_id(cr$r), cm$r))
expect_true(same(by_id(cr$r2), cm$r2))
expect_true(same(sum(as.vector(cr$r2)), sum(as.vector(cm$r2))))
expect_true(same(sort(as.vector(cr$r)), sort(as.vector(cm$r))))
## a computed column combines with the columns it was made from
expect_true(same(by_id(cr$r2 - cr$r * cr$r), rep(0, n)))
expect_true(same(by_id(chain_cols(D)), chain_cols(M)))
expect_true(same(by_id(chain_quad(D)), chain_quad(M)))

## ---- rowSums, rowMeans, quadratic_diagonal, several columns ----
expect_inherits(rowSums(D), "xtx_column")
expect_true(same(by_id(rowSums(D)), rowSums(M)))
expect_true(same(by_id(rowMeans(D)), rowMeans(M)))
expect_true(same(by_id(rowSums("w" * D)), rowSums(d$w * M)))
expect_identical(names(as.data.frame(rowSums(D))), ".rowsum")
expect_identical(names(as.data.frame(rowMeans(D))), ".rowmean")
q <- quadratic_diagonal(D, V)
expect_inherits(q, "xtx_column")
expect_identical(names(as.data.frame(q)), ".quad")
expect_true(same(by_id(q), quadratic_diagonal(M, V)))
expect_true(same(by_id(sqrt(quadratic_diagonal("w" * D, V))), sqrt(quadratic_diagonal(d$w * M, V))))
P <- D %*% B
expect_equal(dim(P), c(n, 2))
expect_equal(length(P), 2 * n)
expect_equal(dim(P), dim(M %*% B))
expect_equal(ncol(P[, 2]), 1)
expect_true(grepl("^<design_column: 300 x 2, on", format(P)[1L]))
## (a subset: its own rows)
expect_equal(dim(D[DBI::SQL("x > 0"), ] %*% b), c(sum(d$x > 0), 1))
expect_identical(names(as.data.frame(P, keep = "id")), c("id", "b1", "b2"))
expect_identical(names(as.data.frame(D %*% unname(B))), c(".fitted1", ".fitted2"))
expect_true(same(by_id(P), M %*% B))
expect_equal(dim(as.matrix(P)), c(n, 2L))
expect_null(dimnames(as.matrix(P)))
expect_error(as.vector(P), "use as.matrix")
expect_true(same(as.matrix(P)[order(as.matrix(P)[, 1]), ], (M %*% B)[order((M %*% B)[, 1]), ]))
expect_true(same(by_id(P[, 2]), (M %*% B)[, 2]))
expect_true(same(by_id(P[, c(FALSE, TRUE)]), (M %*% B)[, 2]))
expect_error(P[1:3, ], "no order")
expect_error(P[, 3], "out of range")

## ---- separate products of one design combine (nh9e) ----
## (a product becomes a view at its first use, stacked on the products it
## is combined with: the same values as in memory)
b2 <- B[, 2]
two <- function(D) D %*% b + D %*% b2
ratio <- function(D) (D %*% b) / sqrt(quadratic_diagonal(D, V))
mixed <- function(D) rowSums(D) - D %*% b + design_column(D, "y") * quadratic_diagonal(D, V)
expect_true(same(by_id(two(D)), two(M)))
expect_true(same(by_id(ratio(D)), ratio(M)))
expect_true(same(by_id(mixed(D)), mixed(M)))
expect_true(same(by_id(D %*% b + rowSums(D)), M %*% b + rowSums(M)))
expect_true(same(sum(as.vector(two(D))), sum(two(M))))
## the products of one expression, kept and reused: each view is made once
e1 <- D %*% b
e2 <- quadratic_diagonal(D, V)
e <- e1 / sqrt(e2)
expect_true(same(by_id(e), ratio(M)))
expect_true(same(by_id(e1), M %*% b))
expect_true(same(by_id(e - e1 / sqrt(e2)), rep(0, n)))
## a computed (stored) column combines with a later product
st <- col_compute(D %*% b)
expect_true(same(by_id(st + D %*% b2), two(M)))
expect_true(same(by_id(D %*% b2 - st), M %*% b2 - M %*% b))
## two products each read before they are combined (epxj): their views are
## siblings, so the one read second is replayed on top of the other (one
## more view, made at the first use of the sum, once); the same values as
## in memory
nviews <- function() nrow(DBI::dbGetQuery(con, "SELECT * FROM duckdb_views() WHERE temporary AND NOT internal"))
p1 <- D %*% b
p2 <- D %*% b2
expect_equal(length(as.vector(p1)), n)
expect_equal(length(as.vector(p2)), n)
v0 <- nviews()
s12 <- p1 + p2
expect_equal(nviews(), v0)                       # combining runs nothing
expect_true(same(by_id(s12), two(M)))
expect_equal(nviews(), v0 + 1L)                  # one view: the replayed product
expect_true(same(sum(s12), sum(two(M))))
expect_equal(nviews(), v0 + 1L)                  # bound once
expect_true(same(by_id(p2 - p1), M %*% b2 - M %*% b))
expect_true(same(by_id(p1 * p2 + p1), (M %*% b) * (M %*% b2) + M %*% b))
## a product and a quadratic form, each read first
q1 <- quadratic_diagonal(D, V)
e1 <- D %*% b
expect_true(same(sum(q1), sum(quadratic_diagonal(M, V))))
expect_true(same(sum(e1), sum(M %*% b)))
expect_true(same(by_id(e1 / sqrt(q1)), ratio(M)))
## three products each read first
p3 <- D %*% (2 * b)
expect_true(same(sum(p3), 2 * sum(M %*% b)))
expect_true(same(by_id(p1 + p2 + p3), two(M) + 2 * (M %*% b)))
## a column stored by col_compute() combines with a product read on its own
## (the product is replayed on the stored table)
st1 <- col_compute(D %*% b)
p4 <- D %*% b2
expect_true(same(sum(p4), sum(M %*% b2)))
expect_true(same(by_id(st1 + p4), two(M)))
## two columns each stored separately cannot be replayed: an error
st2 <- col_compute(D %*% b2)
expect_error(st1 + st2, "two different branches")
expect_error(st1 + st2, "col_compute")
## (stored together, they combine)
st12 <- col_compute(D %*% cbind(b, b2))
expect_true(same(by_id(st12[, 1] + st12[, 2]), two(M)))
## columns of two designs (a subset is another design): an error
expect_error(D %*% b + D[DBI::SQL("x > 0"), ] %*% b, "different remote designs")
## (one product of both: always combines)
P12 <- D %*% cbind(b, b2)
expect_equal(length(as.vector(P12[, 1])), n)
expect_true(same(by_id(P12[, 1] + P12[, 2]), two(M)))

## ---- col_compute(name =): a permanent table, aligned with the design ----
r <- design_column(D, "y") - D %*% b
rt <- col_compute(r, name = "xtx_test_resid")
expect_inherits(rt, "xtx_column")
expect_true(DBI::dbExistsTable(con, "xtx_test_resid"))
tab <- DBI::dbReadTable(con, "xtx_test_resid")
expect_equal(nrow(tab), n)
expect_true(all(c("id", "y", "x", ".value") %in% names(tab)))
tab <- tab[order(tab$id), ]
expect_true(same(tab$.value, d$y - drop(M %*% b)))
expect_true(same(by_id(rt), d$y - drop(M %*% b)))
## it combines with the columns of its design, and with a later product
expect_true(same(by_id(rt - r), rep(0, n)))
expect_true(same(by_id(rt + D %*% b2), d$y - drop(M %*% b) + drop(M %*% b2)))
## an existing table: an error, unless overwrite = TRUE
expect_error(col_compute(r, name = "xtx_test_resid"), "already exists")
rt2 <- col_compute(2 * r, name = "xtx_test_resid", overwrite = TRUE)
expect_true(same(by_id(rt2), 2 * (d$y - drop(M %*% b))))
## value names: the label, `nm`, and a clash with a column of the data
expect_error(col_compute(design_column(D, "y"), name = "xtx_test_y"), "already has a column `y`")
expect_false(DBI::dbExistsTable(con, "xtx_test_y"))
yt <- col_compute(design_column(D, "y"), name = "xtx_test_y", nm = "y_copy")
expect_true("y_copy" %in% DBI::dbListFields(con, "xtx_test_y"))
expect_identical(names(as.data.frame(yt)), "y_copy")
pt <- col_compute(D %*% B, name = "xtx_test_p")
expect_true(all(c("b1", "b2") %in% DBI::dbListFields(con, "xtx_test_p")))
expect_true(same(by_id(pt), M %*% B))
expect_error(col_compute(r, name = "xtx_test_r2", nm = c("a", "b")), "nm")
expect_error(col_compute(r, name = c("a", "b")), "name")
## a subset design: only its rows
Dsub <- D[DBI::SQL("x > 0"), ]
ct <- col_compute(design_column(Dsub, "y") * 2, name = "xtx_test_sub")
expect_equal(nrow(DBI::dbReadTable(con, "xtx_test_sub")), sum(d$x > 0))
for (tn in c("xtx_test_resid", "xtx_test_y", "xtx_test_p", "xtx_test_sub")) DBI::dbRemoveTable(con, tn)

## ---- rows of the design ----
i <- d$x > 0
Ds <- D[DBI::SQL("x > 0"), ]
rs <- as.data.frame(design_column(Ds, "y") - Ds %*% b, keep = "id")
rs <- rs[order(rs$id), ]
expect_identical(rs$id, d$id[i])
expect_true(same(rs$.value, design_column(M[i, ], "y") - M[i, ] %*% b))
Dl <- D[design_column(D, DBI::SQL("x > 0")), ]
expect_equal(nrow(Dl), sum(i))

test_con_done(con)

## ---- lifetime: valid while the connection is open ----
con2 <- DBI::dbConnect(duckdb::duckdb())
DBI::dbWriteTable(con2, "d", d)
col <- design_matrix(fo, con2, table = "d") %*% b
expect_equal(length(as.vector(col)), n)
DBI::dbDisconnect(con2)
expect_true(grepl("(closed)", format(col)[1L], fixed = TRUE))
expect_error(as.vector(col), "connection of this remote column is closed")
expect_error(as.matrix(col), "connection of this remote column is closed")
expect_error(col_compute(col), "closed")
expect_error(col_compute(col, name = "t_closed"), "closed")
expect_error(as.data.frame(col), "closed")

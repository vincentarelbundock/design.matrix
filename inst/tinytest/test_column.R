## The public column type (Kata design.matrix#bvxy): design_column(),
## products of a remote design returning columns, collect(), compute(),
## as.data.frame(), print/format and x[, j]. The same chain of public calls
## gives the same values in memory and on DuckDB. Only exported functions.
source("helpers.R")
library(design.matrix)

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
expect_identical(collect(v), v)
expect_identical(compute(v), v)
Pm <- M %*% B
expect_identical(collect(Pm), Pm)

## the chains of public calls compared below (the same code for both backends)
chain <- function(D) {
    eta <- D %*% b
    r <- design_column(D, "y") - exp(eta / 10)
    list(eta = eta, r = r, r2 = compute(r * r))
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
expect_inherits(pr, "design_column")
expect_equal(length(pr), 1L)
out <- capture.output(print(pr))
expect_true(grepl("^<design_column: 1 column, 300 rows", out[1L]))
expect_true(grepl("duckdb", out[1L]))
expect_identical(format(pr), out)
expect_identical(names(as.data.frame(pr, keep = "id")), c("id", ".fitted"))
expect_identical(names(as.data.frame(pr)), ".fitted")
expect_identical(names(as.data.frame(pr, nm = "eta")), "eta")
expect_true(same(by_id(pr), M %*% b))
expect_true(same(sort(collect(pr)), sort(drop(M %*% b))))

## ---- design_column() ----
y <- design_column(D, "y")
expect_inherits(y, "design_column")
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
expect_inherits(cr$r2, "design_column")
expect_true(same(by_id(cr$eta), cm$eta))
expect_true(same(by_id(cr$r), cm$r))
expect_true(same(by_id(cr$r2), cm$r2))
expect_true(same(sum(collect(cr$r2)), sum(collect(cm$r2))))
expect_true(same(sort(collect(cr$r)), sort(collect(cm$r))))
## a computed column combines with the columns it was made from
expect_true(same(by_id(cr$r2 - cr$r * cr$r), rep(0, n)))
expect_true(same(by_id(chain_cols(D)), chain_cols(M)))
expect_true(same(by_id(chain_quad(D)), chain_quad(M)))

## ---- rowSums, rowMeans, quadratic_diagonal, several columns ----
expect_inherits(rowSums(D), "design_column")
expect_true(same(by_id(rowSums(D)), rowSums(M)))
expect_true(same(by_id(rowMeans(D)), rowMeans(M)))
expect_true(same(by_id(rowSums("w" * D)), rowSums(d$w * M)))
expect_identical(names(as.data.frame(rowSums(D))), ".rowsum")
expect_identical(names(as.data.frame(rowMeans(D))), ".rowmean")
q <- quadratic_diagonal(D, V)
expect_inherits(q, "design_column")
expect_identical(names(as.data.frame(q)), ".quad")
expect_true(same(by_id(q), quadratic_diagonal(M, V)))
expect_true(same(by_id(sqrt(quadratic_diagonal("w" * D, V))), sqrt(quadratic_diagonal(d$w * M, V))))
P <- D %*% B
expect_equal(length(P), 2L)
expect_identical(names(as.data.frame(P, keep = "id")), c("id", "b1", "b2"))
expect_identical(names(as.data.frame(D %*% unname(B))), c(".fitted1", ".fitted2"))
expect_true(same(by_id(P), M %*% B))
expect_equal(dim(collect(P)), c(n, 2L))
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
expect_true(same(sum(collect(two(D))), sum(two(M))))
## the products of one expression, kept and reused: each view is made once
e1 <- D %*% b
e2 <- quadratic_diagonal(D, V)
e <- e1 / sqrt(e2)
expect_true(same(by_id(e), ratio(M)))
expect_true(same(by_id(e1), M %*% b))
expect_true(same(by_id(e - e1 / sqrt(e2)), rep(0, n)))
## a computed (stored) column combines with a later product
st <- compute(D %*% b)
expect_true(same(by_id(st + D %*% b2), two(M)))
expect_true(same(by_id(D %*% b2 - st), M %*% b2 - M %*% b))
## two products each read before they are combined: refused, with the way out
p1 <- D %*% b
p2 <- D %*% b2
expect_equal(length(collect(p1)), n)
expect_equal(length(collect(p2)), n)
expect_error(p1 + p2, "two different branches")
expect_error(p1 + p2, "before reading either")
## (one product of both: always combines)
P12 <- D %*% cbind(b, b2)
expect_equal(length(collect(P12[, 1])), n)
expect_true(same(by_id(P12[, 1] + P12[, 2]), two(M)))

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
expect_equal(length(collect(col)), n)
DBI::dbDisconnect(con2)
expect_true(grepl("(closed)", format(col)[1L], fixed = TRUE))
expect_error(collect(col), "connection of this remote column is closed")
expect_error(compute(col), "closed")
expect_error(as.data.frame(col), "closed")

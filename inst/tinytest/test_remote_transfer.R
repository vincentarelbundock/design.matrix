## What crosses the connection. The database does everything that scales
## with the rows; R receives sums by cell and does the small dense algebra
## (the contrasts). So no operation of a remote design may fetch anything
## proportional to the number of rows, and the sums of a term (a pair of
## terms, for X'X) are grouped by the factors of that term, not by all the
## factors of the design: what is fetched is bounded by the cells of the
## terms involved. record_queries() (helpers.R) records every statement and
## the size of what it returns.
if (!requireNamespace("duckdb", quietly = TRUE) || !requireNamespace("DBI", quietly = TRUE))
    exit_file("duckdb is not installed")
library(xtx)
source("helpers.R")

con <- test_con()
same <- function(a, b, tol = 1e-9) isTRUE(all.equal(unname(a), unname(b), tolerance = tol))
set.seed(3)
n <- 30000
Lf <- 40L; Lg <- 30L; Lh <- 20L; G <- 7L
d <- data.frame(id = seq_len(n), y = rnorm(n), y2 = rnorm(n), x = rnorm(n), w = runif(n, 1, 2),
                grp = sample(seq_len(G), n, TRUE),
                f = sample(sprintf("f%02d", seq_len(Lf)), n, TRUE), g = sample(sprintf("g%02d", seq_len(Lg)), n, TRUE),
                h = sample(sprintf("h%02d", seq_len(Lh)), n, TRUE), stringsAsFactors = FALSE)
DBI::dbWriteTable(con, "d", d)
fo <- y ~ x + f + g + h
M <- model.matrix(fo, d)
p <- ncol(M)
## the fixture: many more cells of the three factors together than cells
## of the terms (one factor each), or of the pairs of terms
joint <- nrow(unique(d[c("f", "g", "h")]))
terms <- Lf + Lg + Lh
pairs <- Lf * Lg + Lf * Lh + Lg * Lh
expect_true(joint > 15000 && joint > 5 * pairs && pairs > 25 * terms)

op <- options(xtx.remote_compress = FALSE)

## preparation, recorded apart: two samples of at most 2000 rows (of the 5
## columns of the model), the levels, counts; nothing that grows with the
## rows
r <- record_queries(D <- design_matrix(fo, con, table = "d"))
expect_true(r$max_rows <= 2000L)
expect_true(r$fetched <= 2 * 2000 * 5 + 2 * terms + 20)
expect_true(all(r$log$up_rows == 0L))

## each operation on the prepared design: its statements, the rows of its
## largest result and the values fetched (rows x columns, over all its
## statements), within `rows` rows, `rows * cols` values and `count`
## statements
budget <- function(expr, rows, cols, count, ref = NULL, info) {
    r <- record_queries(expr)
    if (!is.null(ref)) expect_true(same(r$value, ref), info = paste(info, ": value"))
    expect_true(r$count <= count, info = sprintf("%s: %d statements (at most %d)", info, r$count, count))
    expect_true(r$max_rows <= rows, info = sprintf("%s: %d rows (at most %d)", info, r$max_rows, rows))
    expect_true(r$fetched <= rows * cols, info = sprintf("%s: %.0f values (at most %d)", info, r$fetched, rows * cols))
    expect_true(r$max_rows < n / 10, info = paste(info, ": far fewer rows than the data"))
    expect_true(all(r$log$up_rows == 0L), info = paste(info, ": nothing uploaded"))
    invisible(r)
}

## X'y, column sums: by term, one row per level of each factor (one
## statement, GROUP BY GROUPING SETS). Columns: the keys, the set, and a
## handful of sums.
r <- budget(crossprod(D, "y"), terms, 12L, 1L, crossprod(M, d$y), "X'y")
expect_true(grepl("GROUP BY GROUPING SETS", r$log$sql, fixed = TRUE))
expect_identical(r$log$rows, terms)
budget(crossprod(D, c("y", "y2")), terms, 16L, 1L, crossprod(M, cbind(d$y, d$y2)), "X'Y")
budget(crossprod("w" * D, "y"), terms, 12L, 2L, crossprod(d$w * M, d$y), "(wX)'y")
budget(colSums(D), terms, 8L, 1L, colSums(M), "colSums")
budget(colMeans("w" * D), terms, 9L, 2L, colMeans(d$w * M), "colMeans, scaled")
## sums within groups: by group and term
budget(colSums(D, group = "grp"), G * terms, 4L, 4L, rowsum(M, d$grp), "colSums by group")
budget(colMeans("w" * D, group = "grp"), G * terms, 4L, 5L,
       rowsum(d$w * M, d$grp) / as.vector(table(d$grp)), "colMeans by group, scaled")
## X'X: a block of two factors needs their crosstab, so the cells of the
## pairs of terms; never those of the three factors together
r <- budget(crossprod(D), pairs, 9L, 1L, crossprod(M), "X'X")
expect_true(r$max_rows > terms)
budget(crossprod(D, "w" * D), pairs, 10L, 2L, crossprod(M, d$w * M), "X'WX")
budget(crossprod("w" * D), pairs, 10L, 2L, crossprod(d$w * M), "(wX)'(wX)")

## rows of the design times coefficients, row sums, quadratic forms: lazy
## columns, computed in the database. Building one fetches nothing; a
## summary of one fetches the summary.
b <- rnorm(p)
V <- crossprod(matrix(rnorm(p * p), p)) / p
r <- budget(pr <- D %*% b, 0L, 0L, 1L, info = "D %*% b")
expect_identical(r$fetched, 0)
budget(sum(pr), 1L, 3L, 2L, sum(M %*% b), "sum(D %*% b)")
budget(sum(rowSums(D)), 1L, 3L, 2L, sum(M), "sum(rowSums(D))")
budget(sum(rowMeans(D)), 1L, 3L, 2L, sum(rowMeans(M)), "sum(rowMeans(D))")
budget(sum(quadratic_diagonal(D, V)), 1L, 3L, 2L, sum((M %*% V) * M), "sum(quadratic_diagonal(D, V))")
budget(sum((design_column(D, "y") - D %*% b)^2), 1L, 3L, 2L, sum((d$y - M %*% b)^2), "residual sum of squares")
## (collecting one is the caller's explicit request for its rows)
r <- record_queries(head(as.vector(pr)))
expect_true(r$max_rows <= n)

## summary(D, sql = TRUE) shows the statement that crossprod(D) sends, and
## sends nothing itself
r <- record_queries(s <- summary(D, sql = TRUE))
expect_identical(r$count, 0L)
shown <- sub("^-- crossprod\n", "", attr(s, "sql")[2L])
r <- record_queries(crossprod(D))
expect_identical(r$log$sql, shown)
expect_true(grepl("GROUP BY GROUPING SETS", shown, fixed = TRUE))

## the same budgets when the factors are ENUM and logical columns, with
## rows kept that have missing values (na.pass): X'y by term
d2 <- d[seq_len(5000), ]
d2$f <- factor(d2$f)
d2$l <- d2$grp > 3
d2$g[c(3, 50)] <- NA
d2$x[7] <- NA
DBI::dbWriteTable(con, "d2", d2)
m2 <- DBI::dbReadTable(con, "d2")
D2 <- design_matrix(D, newdata = con, table = "d2")
M2 <- model.matrix(fo, model.frame(fo, m2, na.action = na.pass, xlev = list(f = levels(m2$f))))
budget(crossprod(D2, "y"), terms + 1L, 12L, 2L, crossprod(M2, m2$y), "X'y, na.pass")
budget(colSums(D2), terms + 1L, 8L, 2L, colSums(M2), "colSums, na.pass")

## with one factor set inside another (f, f:g), the sums of the smaller
## come from the cells of the larger: one plain GROUP BY, as few rows as
## the cells of the term f:g
D3 <- design_matrix(y ~ x * f + f:g, con, table = "d")
M3 <- model.matrix(y ~ x * f + f:g, d)
r <- budget(crossprod(D3, "y"), Lf * Lg, 8L, 1L, crossprod(M3, d$y), "X'y, nested terms")
expect_false(grepl("GROUPING SETS", r$log$sql, fixed = TRUE))

options(op)
test_con_done(con)

## The primitive layer of the fitting core (R/primitives.R): in-memory
## methods against base R, and the remote methods against them.
source("helpers.R")
library(design.matrix)
dmp_mult <- design.matrix:::dmp_mult
dmp_scale <- design.matrix:::dmp_scale
dmp_crossprod <- design.matrix:::dmp_crossprod
dmp_colsums <- design.matrix:::dmp_colsums
dmp_groupsum <- design.matrix:::dmp_groupsum
dmp_crosstab <- design.matrix:::dmp_crosstab
dmp_quad <- design.matrix:::dmp_quad
dmp_patterns <- design.matrix:::dmp_patterns
dmp_output <- design.matrix:::dmp_output
dmp_column <- design.matrix:::dmp_column
dmp_collect <- design.matrix:::dmp_collect
dm_colsum <- design.matrix:::dm_colsum
dmp_colnorms <- design.matrix:::dmp_colnorms
same <- function(a, b, tol = 1e-8) isTRUE(all.equal(unname(a), unname(b), tolerance = tol))

set.seed(5)
n <- 300
d <- data.frame(y = rnorm(n), x = rnorm(n), xd = sample(1:3, n, TRUE),
                f = sample(c("a", "b", "c"), n, TRUE), g = sample(1:7, n, TRUE), w = runif(n))
D <- design_matrix(y ~ x * f + xd, d)
X <- as.matrix(D)
B <- matrix(rnorm(ncol(X) * 2), ncol(X))
V <- crossprod(matrix(rnorm(ncol(X)^2), ncol(X)))

## in memory: the primitives are the existing kernels
expect_true(same(dmp_mult(D, B), X %*% B))
expect_true(same(dmp_crossprod(D), crossprod(X)))
expect_true(same(dmp_crossprod(D, dmp_scale(D, d$w)), crossprod(X, d$w * X)))
expect_true(same(dmp_crossprod(D, cbind(d$y, d$w)), crossprod(X, cbind(d$y, d$w))))
expect_true(same(dmp_colsums(dmp_scale(D, d$w), d$g), rowsum(d$w * X, d$g)))
expect_true(same(dmp_colsums(D), colSums(X)))
expect_true(same(dmp_groupsum(cbind(d$y, d$w), d$g), rowsum(cbind(d$y, d$w), d$g)))
ct <- dmp_crosstab(d$f, d$g, d$w)
expect_inherits(ct, "dgCMatrix")
expect_true(same(as.vector(as.matrix(ct)), as.vector(xtabs(w ~ f + g, d))))
expect_true(same(as.matrix(dmp_crosstab(d$f, d$g)), unclass(table(d$f, d$g))))
expect_true(same(dmp_quad(D, V), rowSums((X %*% V) * X)))
expect_true(same(dm_colsum(d$y), sum(d$y)))
expect_true(same(dmp_colnorms(X %*% B), sqrt(colSums((X %*% B)^2))))
expect_identical(dmp_column(D, d$y), d$y)
expect_identical(dmp_collect(d$y), d$y)
out <- dmp_output(D, list(a = d$y, b = 2 * d$y))
expect_inherits(out, "data.frame")
expect_equal(out$b, 2 * d$y)
## compressed least squares: the fit on the patterns is the fit on the rows
Dd <- design_matrix(y ~ xd + f, d)
pt <- dmp_patterns(Dd, d$y, d$w, compress = TRUE)
expect_true(nrow(pt$x) < n)
l <- lm(y ~ xd + f, d, weights = w)
lp <- lm.wfit(as.matrix(pt$x), pt$z, pt$w)
expect_true(same(coef(lp), coef(l)))
expect_true(same(sum(lp$weights * lp$residuals^2) + pt$within, sum(l$weights * l$residuals^2)))
expect_null(dmp_patterns(D, d$y, NULL, compress = FALSE))

## columns: the same R code on numbers and on remote columns
if (requireNamespace("duckdb", quietly = TRUE)) {
    con <- DBI::dbConnect(duckdb::duckdb())
    DBI::dbWriteTable(con, "d", d)
    dd <- DBI::dbReadTable(con, "d")
    R <- design_matrix(y ~ x * f + xd, con, table = "d")
    Xr <- as.matrix(design_matrix(y ~ x * f + xd, dd))
    y <- dmp_column(R, 'CAST(s."y" AS DOUBLE)')
    w <- dmp_column(R, 'CAST(s."w" AS DOUBLE)')
    g <- dmp_column(R, 's."g"')
    f <- function(y, w) sum(dm_colsum(exp(y) * w - sqrt(abs(y)) / (1 + w^2)))
    expect_true(same(f(y, w), f(dd$y, dd$w)))
    expect_true(same(dmp_crossprod(R), crossprod(Xr)))
    expect_true(same(dmp_crossprod(R, dmp_scale(R, w)), crossprod(Xr, dd$w * Xr)))
    expect_true(same(dmp_crossprod(R, y * w), crossprod(Xr, dd$y * dd$w)))
    expect_true(same(dmp_colsums(dmp_scale(R, w), g), rowsum(dd$w * Xr, dd$g)))
    expect_true(same(dmp_groupsum(y * w, g), rowsum(dd$y * dd$w, dd$g)))
    expect_true(same(as.vector(as.matrix(dmp_crosstab(dmp_column(R, 's."f"'), g, w))), as.vector(xtabs(w ~ f + g, dd))))
    ## products are columns: residuals, refinement, quadratic forms
    b <- rnorm(ncol(Xr))
    r <- y - dmp_mult(R, b)
    expect_true(same(dm_colsum(r * r), sum((dd$y - Xr %*% b)^2)))
    expect_true(same(dmp_crossprod(R, r), crossprod(Xr, dd$y - Xr %*% b)))
    q <- dmp_quad(R, V, along = list(r))
    expect_true(same(dm_colsum(q * r), sum(rowSums((Xr %*% V) * Xr) * (dd$y - Xr %*% b))))
    expect_true(same(dmp_colnorms(dmp_mult(R, B)), sqrt(colSums((Xr %*% B)^2))))
    expect_true(same(sort(dmp_collect(r)), sort(drop(dd$y - Xr %*% b))))
    expect_true(same(dm_colsum(design.matrix:::dmc_pnorm(y)), sum(pnorm(dd$y)), 1e-13))
    expect_true(same(dm_colsum(design.matrix:::dmc_qnorm(w)), sum(qnorm(dd$w)), 1e-13))
    ## row-wise output: a table with the columns of the data
    tab <- dmp_output(R, list(.r = r, .q = q))
    o <- as.data.frame(tab)
    expect_true(all(c("y", "x", ".r", ".q") %in% names(o)))
    expect_true(same(sort(o$.r), sort(drop(dd$y - Xr %*% b))))
    ## NA semantics: sums with a NULL are NA
    DBI::dbWriteTable(con, "na", data.frame(y = c(1, NA, 3), x = 1:3))
    Rn <- design_matrix(~ x, con, table = "na")
    expect_true(is.na(dm_colsum(dmp_column(Rn, 'CAST(s."y" AS DOUBLE)'))))
    DBI::dbDisconnect(con)
}

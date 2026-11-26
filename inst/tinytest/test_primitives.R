## The primitive layer of the fitting core (R/primitives.R): in-memory
## methods against base R, and the remote methods against them.
source("helpers.R")
library(xtx)
xtx_p_mult <- xtx:::xtx_p_mult
xtx_p_scale <- xtx:::xtx_p_scale
xtx_p_crossprod <- xtx:::xtx_p_crossprod
xtx_p_colsums <- xtx:::xtx_p_colsums
xtx_p_groupsum <- xtx:::xtx_p_groupsum
xtx_p_crosstab <- xtx:::xtx_p_crosstab
xtx_p_quad <- xtx:::xtx_p_quad
xtx_p_patterns <- xtx:::xtx_p_patterns
xtx_p_output <- xtx:::xtx_p_output
xtx_p_column <- xtx:::xtx_p_column
xtx_p_collect <- xtx:::xtx_p_collect
xtx_colsum <- xtx:::xtx_colsum
xtx_p_colnorms <- xtx:::xtx_p_colnorms
same <- function(a, b, tol = 1e-8) isTRUE(all.equal(unname(a), unname(b), tolerance = tol))

set.seed(5)
n <- 300
d <- data.frame(y = rnorm(n), x = rnorm(n), xd = sample(1:3, n, TRUE),
                f = sample(c("a", "b", "c"), n, TRUE), g = sample(1:7, n, TRUE), w = runif(n))
D <- design_matrix(y ~ x * f + xd, d)
X <- as.matrix(D)
B <- matrix(rnorm(ncol(X) * 2), ncol(X))
V <- crossprod(matrix(rnorm(ncol(X)^2), ncol(X)))

## no function of the package takes `along` (nh9e: a remote product or
## lookup finds its frame at its first use)
ns <- asNamespace("xtx")
fns <- Filter(function(f) is.function(get(f, envir = ns)), ls(ns, all.names = TRUE))
expect_identical(fns[vapply(fns, function(f) "along" %in% names(formals(get(f, envir = ns))), NA)],
                 character(0))

## in memory: the primitives are the existing kernels
expect_true(same(xtx_p_mult(D, B), X %*% B))
expect_true(same(xtx_p_crossprod(D), crossprod(X)))
expect_true(same(xtx_p_crossprod(D, xtx_p_scale(D, d$w)), crossprod(X, d$w * X)))
expect_true(same(xtx_p_crossprod(D, cbind(d$y, d$w)), crossprod(X, cbind(d$y, d$w))))
expect_true(same(xtx_p_colsums(xtx_p_scale(D, d$w), d$g), rowsum(d$w * X, d$g)))
expect_true(same(xtx_p_colsums(D), colSums(X)))
expect_true(same(xtx_p_groupsum(cbind(d$y, d$w), d$g), rowsum(cbind(d$y, d$w), d$g)))
ct <- xtx_p_crosstab(d$f, d$g, d$w)
expect_inherits(ct, "dgCMatrix")
expect_true(same(as.vector(as.matrix(ct)), as.vector(xtabs(w ~ f + g, d))))
expect_true(same(as.matrix(xtx_p_crosstab(d$f, d$g)), unclass(table(d$f, d$g))))
expect_true(same(xtx_p_quad(D, V), rowSums((X %*% V) * X)))
expect_true(same(xtx_colsum(d$y), sum(d$y)))
expect_true(same(xtx_p_colnorms(X %*% B), sqrt(colSums((X %*% B)^2))))
expect_identical(xtx_p_column(D, d$y), d$y)
expect_identical(xtx_p_collect(d$y), d$y)
out <- xtx_p_output(D, list(a = d$y, b = 2 * d$y))
expect_inherits(out, "data.frame")
expect_equal(out$b, 2 * d$y)
## compressed least squares: the fit on the patterns is the fit on the rows
Dd <- design_matrix(y ~ xd + f, d)
pt <- xtx_p_patterns(Dd, d$y, d$w, compress = TRUE)
expect_true(nrow(pt$x) < n)
l <- lm(y ~ xd + f, d, weights = w)
lp <- lm.wfit(as.matrix(pt$x), pt$z, pt$w)
expect_true(same(coef(lp), coef(l)))
expect_true(same(sum(lp$weights * lp$residuals^2) + pt$within, sum(l$weights * l$residuals^2)))
expect_null(xtx_p_patterns(D, d$y, NULL, compress = FALSE))

## columns: the same R code on numbers and on remote columns
if (requireNamespace("duckdb", quietly = TRUE)) {
    con <- test_con()
    DBI::dbWriteTable(con, "d", d)
    dd <- DBI::dbReadTable(con, "d")
    R <- design_matrix(y ~ x * f + xd, con, table = "d")
    Xr <- as.matrix(design_matrix(y ~ x * f + xd, dd))
    y <- xtx_p_column(R, 'CAST(s."y" AS DOUBLE)')
    w <- xtx_p_column(R, 'CAST(s."w" AS DOUBLE)')
    g <- xtx_p_column(R, 's."g"')
    f <- function(y, w) sum(xtx_colsum(exp(y) * w - sqrt(abs(y)) / (1 + w^2)))
    expect_true(same(f(y, w), f(dd$y, dd$w)))
    expect_true(same(xtx_p_crossprod(R), crossprod(Xr)))
    expect_true(same(xtx_p_crossprod(R, xtx_p_scale(R, w)), crossprod(Xr, dd$w * Xr)))
    expect_true(same(xtx_p_crossprod(R, y * w), crossprod(Xr, dd$y * dd$w)))
    expect_true(same(xtx_p_colsums(xtx_p_scale(R, w), g), rowsum(dd$w * Xr, dd$g)))
    expect_true(same(xtx_p_groupsum(y * w, g), rowsum(dd$y * dd$w, dd$g)))
    expect_true(same(as.vector(as.matrix(xtx_p_crosstab(xtx_p_column(R, 's."f"'), g, w))), as.vector(xtabs(w ~ f + g, dd))))
    ## products are columns: residuals, refinement, quadratic forms
    b <- rnorm(ncol(Xr))
    r <- y - xtx_p_mult(R, b)
    expect_true(same(xtx_colsum(r * r), sum((dd$y - Xr %*% b)^2)))
    expect_true(same(xtx_p_crossprod(R, r), crossprod(Xr, dd$y - Xr %*% b)))
    ## (pending: binds on the frame of r at its first use with r)
    q <- xtx_p_quad(R, V)
    expect_true(same(xtx_colsum(q * r), sum(rowSums((Xr %*% V) * Xr) * (dd$y - Xr %*% b))))
    expect_true(same(xtx_p_colnorms(xtx_p_mult(R, B)), sqrt(colSums((Xr %*% B)^2))))
    expect_true(same(sort(xtx_p_collect(r)), sort(drop(dd$y - Xr %*% b))))
    expect_true(same(xtx_colsum(xtx:::xtx_c_pnorm(y)), sum(pnorm(dd$y)), 1e-13))
    expect_true(same(xtx_colsum(xtx:::xtx_c_qnorm(w)), sum(qnorm(dd$w)), 1e-13))
    ## row-wise output: a table with the columns of the data
    tab <- xtx_p_output(R, list(.r = r, .q = q))
    o <- as.data.frame(tab)
    expect_true(all(c("y", "x", ".r", ".q") %in% names(o)))
    expect_true(same(sort(o$.r), sort(drop(dd$y - Xr %*% b))))
    ## NA semantics: sums with a NULL are NA
    DBI::dbWriteTable(con, "na", data.frame(y = c(1, NA, 3), x = 1:3))
    Rn <- design_matrix(~ x, con, table = "na")
    expect_true(is.na(xtx_colsum(xtx_p_column(Rn, 'CAST(s."y" AS DOUBLE)'))))
    test_con_done(con)
}

## ---- contract extensions: coded groups, subsets, arguments ----
xtx_c_codes <- xtx:::xtx_c_codes
xtx_p_colsums_groups <- xtx:::xtx_p_colsums_groups
xtx_p_subset <- xtx:::xtx_p_subset
xtx_p_resolve <- xtx:::xtx_p_resolve
xtx_c_cbind <- xtx:::xtx_c_cbind
cg <- xtx_c_codes(d$g, 7L)
cf <- xtx_c_codes(match(d$f, c("a", "b", "c")), 3L)
expect_true(same(xtx_p_groupsum(cbind(d$y, d$w), cg), rowsum(cbind(d$y, d$w), d$g)))
expect_true(same(as.matrix(xtx_p_crosstab(cf, cg, d$w)), as.matrix(xtx_p_crosstab(d$f, d$g, d$w))))
gs <- xtx_p_colsums_groups(xtx_p_scale(D, d$w), list(cg, cf), extra = d$y)
expect_true(same(gs[[1]]$S, rowsum(d$w * X, d$g)))
expect_true(same(gs[[2]]$S, rowsum(d$w * X, d$f)))
expect_true(same(gs[[1]]$n, drop(rowsum(d$w, d$g))))
expect_true(same(gs[[2]]$extra, rowsum(d$w * d$y, d$f)))
gs0 <- xtx_p_colsums_groups(D, list(cg))
expect_true(same(gs0[[1]]$S, rowsum(X, d$g)) && same(gs0[[1]]$n, tabulate(d$g)))
keep <- d$x > 0
expect_true(same(c(as.matrix(xtx_p_subset(D, keep))), c(X[keep, ])))
expect_true(same(xtx_p_mult(D, B[, 1]), drop(X %*% B[, 1])))
expect_identical(xtx_p_resolve(D, list(w = quote(2 * w)), environment(), data = d)$cols$w, 2 * d$w)
## patterns split by keys (an offset, say): exact sums
kz <- round(d$x)
ptk <- xtx_p_patterns(Dd, d$y - kz, d$w, compress = TRUE, keys = kz)
lk <- lm(y ~ xd + f, d, weights = w, offset = kz)
lpk <- lm.wfit(as.matrix(ptk$x), ptk$z, ptk$w)
expect_true(same(coef(lpk), coef(lk)))
expect_equal(ncol(ptk$keys), 1L)

if (requireNamespace("duckdb", quietly = TRUE)) {
    con <- test_con()
    DBI::dbWriteTable(con, "d", d)
    dd <- DBI::dbReadTable(con, "d")
    R <- design_matrix(y ~ x * f + xd, con, table = "d")
    Xr <- as.matrix(design_matrix(y ~ x * f + xd, dd))
    y <- xtx_p_column(R, 'CAST(s."y" AS DOUBLE)')
    w <- xtx_p_column(R, 'CAST(s."w" AS DOUBLE)')
    gr <- xtx_c_codes(xtx_p_column(R, 's."g"'), 7L)
    fr <- xtx_c_codes(xtx_p_column(R, 'CASE s."f" WHEN \'a\' THEN 1 WHEN \'b\' THEN 2 ELSE 3 END'), 3L)
    expect_true(same(xtx_p_groupsum(xtx_c_cbind(y, w), gr), rowsum(cbind(dd$y, dd$w), dd$g)))
    expect_true(same(as.matrix(xtx_p_crosstab(fr, gr, w)), as.matrix(xtx_p_crosstab(dd$f, dd$g, dd$w))))
    gsr <- xtx_p_colsums_groups(xtx_p_scale(R, w), list(gr, fr), extra = y)
    expect_true(same(gsr[[1]]$S, rowsum(dd$w * Xr, dd$g)) && same(gsr[[2]]$extra, rowsum(dd$w * dd$y, dd$f)))
    expect_true(same(gsr[[1]]$n, drop(rowsum(dd$w, dd$g))))
    ## a group with a missing value sums to NA, as rowsum()
    yn <- xtx:::xtx_c_ifelse(xtx_p_column(R, 's."g"') == 3, xtx:::xtx_p_column(R, "CAST(NULL AS DOUBLE)"), y)
    expect_identical(unname(is.na(xtx_p_groupsum(yn, xtx_p_column(R, 's."g"'))[, 1])), sort(unique(dd$g)) == 3)
    ## crosstab labels as in memory
    ct <- xtx_p_crosstab(xtx_p_column(R, 's."f"'), xtx_p_column(R, 's."g"'))
    expect_identical(dimnames(ct), dimnames(xtx_p_crosstab(dd$f, dd$g)))
    ## subsets: a logical column or a SQL predicate
    Rs <- xtx_p_subset(R, xtx_p_column(R, 'CAST(s."x" AS DOUBLE)') > 0)
    expect_true(same(crossprod(Rs), crossprod(Xr[dd$x > 0, ])))
    expect_true(same(crossprod(xtx_p_subset(R, DBI::SQL("x > 0"))), crossprod(Xr[dd$x > 0, ])))
    expect_true(same(xtx_colsum(xtx:::xtx_c_on(y, Rs)), sum(dd$y[dd$x > 0])))
    ## arguments: SQL, translatable R, and R evaluated in batches
    myf <- function(v) v^2 + 1
    rr <- xtx_p_resolve(R, list(a = quote(2 * w), b = DBI::SQL("x + 1"), c = quote(myf(x)), s = quote(x > 0)),
                      environment(), logical = "s")
    expect_identical(names(rr$evaluated), "c")
    expect_true(same(xtx_colsum(rr$cols$a), sum(2 * dd$w)))
    expect_true(same(xtx_colsum(rr$cols$b), sum(dd$x + 1)))
    expect_true(same(xtx_colsum(rr$cols$c), sum(dd$x^2 + 1)))
    expect_true(same(crossprod(xtx_p_subset(rr$x, rr$cols$s)), crossprod(Xr[dd$x > 0, ])))
    ## patterns with keys, remotely
    Rd <- design_matrix(y ~ xd + f, con, table = "d")
    kzr <- xtx:::xtx_c_round(xtx_p_column(Rd, 'CAST(s."x" AS DOUBLE)'))
    ptr <- xtx_p_patterns(Rd, xtx_p_column(Rd, 'CAST(s."y" AS DOUBLE)') - kzr, xtx_p_column(Rd, 'CAST(s."w" AS DOUBLE)'),
                        compress = TRUE, keys = kzr)
    lpr <- lm.wfit(as.matrix(ptr$x), ptr$z, ptr$w)
    expect_true(same(coef(lpr), coef(lm(y ~ xd + f, dd, weights = w, offset = round(x)))))
    expect_true(same(xtx:::xtx_p_collect(xtx:::xtx_c_lgamma(w))[1:3], lgamma(dd$w[1:3])))
    test_con_done(con)
}

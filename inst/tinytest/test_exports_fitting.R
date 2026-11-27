## The exported functions that fitting code is written with (xtx.fit uses
## no internal function of xtx): col_math(), col_pmin() and col_pmax(checked
## = FALSE), design_moments(), qr(D, y =), pivoted_chol(), design_subset(),
## design_align(), design_key(), key_complete(), key_codes(), group_codes(),
## the group functions on bare codes, design_group_sums(), design_source()
## and design_remote_frame(). In memory and on DuckDB, against base R.
source("helpers.R")
library(xtx)

same <- function(a, b, tol = 1e-10) isTRUE(all.equal(as.vector(a), as.vector(b), tolerance = tol))

set.seed(31)
n <- 200
d <- data.frame(id = seq_len(n), y = rnorm(n), x = rnorm(n), w = runif(n, 0.5, 2),
                f = factor(sample(letters[1:4], n, TRUE)), a = sample(1:12, n, TRUE),
                b = sample(c("u", "v", "w"), n, TRUE), stringsAsFactors = FALSE)
d$a[c(3, 50)] <- NA
D <- design_matrix(y ~ x + f, d)

## ---- col_math(), col_pmin(), col_pmax() ----
u <- seq(-2, 2, length.out = 21)
expect_identical(col_math(u, "pnorm"), pnorm(u))
expect_identical(col_math(plogis(u), "qnorm"), qnorm(plogis(u)))
expect_identical(col_math(u, "dnorm"), dnorm(u))
expect_identical(col_math(u, "pcauchy"), pcauchy(u))
expect_identical(col_math(exp(u), "log_positive"), log(exp(u)))
expect_identical(col_math(exp(u), "one_minus_exp_neg"), -expm1(-exp(u)))
expect_error(col_math(u, "nope"))
expect_error(col_math(D, "pnorm"), "not a column")
expect_identical(col_pmax(u, 0, checked = FALSE), pmax(u, 0))
expect_identical(col_pmin(u, 0, checked = FALSE), pmin(u, 0))
expect_identical(col_pmax(c(1, NA), 0), pmax(c(1, NA), 0))
expect_error(col_pmax(u, 0, checked = NA))

## ---- design_moments() ----
m <- design_moments(D, d$y, d$w)
X <- as.matrix(D)
expect_true(same(m$xtx, crossprod(X, d$w * X)))
expect_true(same(m$xty, crossprod(X, d$w * d$y)))
expect_true(same(m$yty, sum(d$w * d$y^2)))
expect_true(same(m$tss, sum(d$w * (d$y - weighted.mean(d$y, d$w))^2)))
m0 <- design_moments(D, d$y)
expect_true(same(m0$tss, sum((d$y - mean(d$y))^2)))

## ---- pivoted_chol() ----
Z <- cbind(1, x = 1:6, x2 = 2 * (1:6), z = c(1, 0, 2, 1, 3, 0))
G <- crossprod(Z)
ch <- pivoted_chol(G, diag(G), 1e-14)
expect_identical(ch$rank, 3L)
expect_identical(ch$pivot, c(1L, 2L, 4L, 3L))
expect_true(same(crossprod(ch$R[1:3, 1:3]), G[c(1, 2, 4), c(1, 2, 4)]))

## ---- design_subset(), design_align() ----
keep <- d$x > 0
s <- design_subset(D, keep, list(y = d$y, more = list(w = d$w, none = NULL), key = d["a"]))
expect_identical(dim(s$x), c(sum(keep), ncol(D)))
expect_identical(s$cols$y, d$y[keep])
expect_identical(s$cols$more$w, d$w[keep])
expect_null(s$cols$more$none)
expect_identical(s$cols$key, d[keep, "a", drop = FALSE])
## missing values of `keep` count as FALSE
kn <- keep
kn[1:3] <- NA
expect_identical(nrow(design_subset(D, kn)$x), sum(keep[-(1:3)]))
al <- design_align(D, list(y = d$y))
expect_identical(al$cols$y, d$y)
expect_identical(dim(al$x), dim(D))

## ---- keys and bare codes ----
k <- design_key(D, "a", d)
expect_true(is.data.frame(k))
ok <- key_complete(k)
expect_identical(ok, !is.na(d$a))
expect_error(design_key(D, "nosuch", d), "No variable")
s <- design_subset(D, ok, list(k = k, y = d$y, w = d$w, b = d["b"]))
cd <- key_codes(s$cols$k)
ref <- sort(unique(d$a))
expect_identical(cd$levels$a, ref)
expect_identical(cd$codes$G, length(ref))
expect_identical(cd$codes$codes, match(d$a[ok], ref))
expect_equal(cd$n, as.vector(table(d$a[ok])))
g <- cd$codes
expect_inherits(g, "xtx_codes")
expect_identical(group_codes(g), g)
## the group functions on bare codes, against rowsum() and table()
expect_true(same(group_sums(s$cols$y, g), rowsum(d$y[ok], d$a[ok])))
gb <- key_codes(s$cols$b)$codes
ct <- group_crosstab(g, gb)
expect_true(same(as.matrix(ct), table(d$a[ok], d$b[ok])))
expect_true(same(as.matrix(group_crosstab(g, gb, s$cols$w)), xtabs(w ~ a + b, d[ok, ])))
expect_identical(group_lookup(seq_len(g$G) * 10, g), match(d$a[ok], ref) * 10)
expect_false(group_nested(g, gb))
expect_true(group_nested(g, g))
## several groupings and a total in one call
gs <- group_sums(s$cols$y, list(g, gb), total = s$cols$y^2)
expect_true(same(gs$groups[[1L]], rowsum(d$y[ok], d$a[ok])))
expect_true(same(gs$groups[[2L]], rowsum(d$y[ok], d$b[ok])))
expect_true(same(gs$total, sum(d$y[ok]^2)))
expect_error(group_sums(d$y, design_groups(D, d["b"]), total = d$y), "list of groupings")
## the bare codes of a grouping of design_groups()
gg <- design_groups(D, d["b"])
expect_identical(group_codes(gg)$codes, gg$codes)
expect_identical(group_codes(gg)$G, gg$G)
expect_error(group_codes(1), "design_groups")
## grouped sums of the design, for several groupings
Xs <- as.matrix(s$x)
dg <- design_group_sums(s$x, list(g, gb), w = s$cols$w, extra = s$cols$y)
expect_true(same(dg[[1L]]$S, rowsum(s$cols$w * Xs, d$a[ok])))
expect_true(same(dg[[2L]]$n, rowsum(s$cols$w, d$b[ok])))
expect_true(same(dg[[1L]]$extra, rowsum(s$cols$w * s$cols$y, d$a[ok])))
fs <- design_group_sums(s$x, list(g, gb), w = s$cols$w, extra = s$cols$y, pairs = list(c(1L, 2L)))
expect_true(same(fs$gs[[1L]]$S, dg[[1L]]$S))
expect_true(same(fs$gs[[2L]]$n, dg[[2L]]$n))
expect_identical(length(fs$ct), 1L)

## ---- remote: the same on DuckDB ----
if (requireNamespace("duckdb", quietly = TRUE)) {
    con <- test_con()
    DBI::dbWriteTable(con, "d", d)
    R <- design_matrix(y ~ x + f, con, table = "d")
    yr <- design_column(R, "y")
    wr <- design_column(R, "w")
    xr <- design_column(R, "x")
    ## col_math(), unchecked extremes
    expect_true(same(sum(col_math(xr, "pnorm")), sum(pnorm(d$x))))
    expect_true(same(sum(col_math(xr, "dnorm")), sum(dnorm(d$x))))
    expect_true(same(sum(col_math(exp(xr), "log_positive")), sum(d$x)))
    expect_true(same(sum(col_math(exp(xr), "one_minus_exp_neg")), sum(-expm1(-exp(d$x)))))
    expect_true(same(sum(col_pmax(xr, 0, checked = FALSE)), sum(pmax(d$x, 0))))
    expect_true(same(sum(col_pmin(xr, 0, checked = FALSE)), sum(pmin(d$x, 0))))
    ## design_moments()
    mr <- design_moments(R, yr, wr)
    expect_true(same(mr$xtx, m$xtx, 1e-8))
    expect_true(same(mr$xty, m$xty, 1e-8))
    expect_true(same(mr$tss, m$tss, 1e-8))
    ## qr(D, y =): X'y from the same scan, the same coefficients
    q <- qr(R, y = yr)
    expect_true(same(qr.coef(q, yr), coef(lm(y ~ x + f, d)), 1e-8))
    expect_true(same(qr.coef(qr(R), yr), coef(lm(y ~ x + f, d)), 1e-8))
    expect_error(qr(R, y = "y"), "column of the design")
    ## design_subset(): the rows, and the columns over them
    sr <- design_subset(R, xr > 0, list(y = yr, none = NULL))
    expect_identical(dim(sr$x)[1L], as.numeric(sum(keep)))
    expect_true(same(sum(sr$cols$y), sum(d$y[keep])))
    expect_null(sr$cols$none)
    ## keys with missing values, bare codes on the stored values
    kr <- design_key(R, "a")
    okr <- key_complete(kr)
    expect_equal(sum(okr), sum(ok))
    s2 <- design_subset(R, okr, list(k = kr, y = yr, w = wr, b = design_key(R, "b")))
    cdr <- key_codes(s2$cols$k)
    expect_identical(as.integer(cdr$levels[[1L]]), ref)
    expect_equal(as.numeric(cdr$n), as.vector(table(d$a[ok])))
    gr <- cdr$codes
    expect_true(same(group_sums(s2$cols$y, gr), rowsum(d$y[ok], d$a[ok])))
    gbr <- key_codes(s2$cols$b)$codes
    expect_true(same(as.matrix(group_crosstab(gr, gbr)), table(d$a[ok], d$b[ok])))
    gsr <- group_sums(s2$cols$y, list(gr, gbr), total = s2$cols$y^2)
    expect_true(same(gsr$groups[[2L]], rowsum(d$y[ok], d$b[ok])))
    expect_true(same(gsr$total, sum(d$y[ok]^2)))
    dgr <- design_group_sums(s2$x, list(gr, gbr), w = s2$cols$w)
    expect_true(same(dgr[[1L]]$S, dg[[1L]]$S, 1e-8))
    lk <- group_lookup(seq_len(gr$G) * 10, gr)
    expect_true(same(sum(lk), sum(match(d$a[ok], ref) * 10)))
    ## design_source(), design_remote_frame()
    src <- design_source(con, "d")
    expect_inherits(src, "xtx_src")
    expect_identical(unname(src$types[c("y", "b")]), c("double", "character"))
    expect_error(design_source(con), "exactly one")
    myf <- function(v) v * 3
    fr <- design_remote_frame(y ~ x + f, src, list(y = quote(y), weights = quote(2 * w), subset = quote(x > 0)),
                              environment(), logical = "subset")
    expect_inherits(fr$x, "xtx_remote")
    expect_identical(dim(fr$x), c(as.numeric(sum(keep)), ncol(D)))
    expect_true(same(sum(fr$cols$weights), 2 * sum(d$w[keep])))
    expect_true(same(sum(fr$cols$y), sum(d$y[keep])))
    expect_null(fr$off)
    ## an argument evaluated in R, and an offset term
    fr2 <- suppressMessages(design_remote_frame(y ~ x + offset(w), src, list(weights = quote(myf(w))), environment()))
    expect_true(same(sum(fr2$cols$weights), 3 * sum(d$w)))
    expect_true(same(sum(fr2$off), sum(d$w)))
    test_con_done(con)
}

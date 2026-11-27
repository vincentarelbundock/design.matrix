## The exported functions that fitting code is written with (xtx.fit uses
## no internal function of xtx): col_math(), col_pmin() and col_pmax(checked
## = FALSE), design_moments(), qr(D, y =), pivoted_chol(), design_subset()
## (with keep = NULL: the alignment), design_key(), is.na() of a key,
## key_codes(), the group functions on bare codes, design_group_sums() (bare
## codes and groupings), and the remote path of design_frame() with
## design_resolve(response =). In memory and on DuckDB, against base R.
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
expect_error(col_math(u, "nope"))
expect_identical(col_ifelse(c(TRUE, NA, FALSE), 1, 2, missing = 9), c(1, 9, 2))
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

## ---- design_subset() ----
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
## keep = NULL: no row is dropped, the design and the columns are aligned
al <- design_subset(D, cols = list(y = d$y))
expect_identical(al$cols$y, d$y)
expect_identical(dim(al$x), dim(D))
expect_identical(design_subset(D, NULL, list(y = d$y))$cols$y, d$y)

## ---- keys and bare codes ----
k <- design_key(D, "a", d)
expect_true(is.data.frame(k))
ok <- stats::complete.cases(k)
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
## groupings of design_groups(), alone or with bare codes
sD <- design_matrix(y ~ x + f, d[ok, ])
gg <- design_groups(sD, d[ok, "b", drop = FALSE])
dq <- design_group_sums(sD, list(g, gg), w = s$cols$w, extra = s$cols$y)
expect_true(same(dq[[1L]]$S, dg[[1L]]$S))
expect_true(same(dq[[2L]]$S, dg[[2L]]$S))
expect_true(same(dq[[2L]]$extra, dg[[2L]]$extra))
expect_error(design_group_sums(sD, list(1)), "design_groups")

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
    ## expm1() near zero, both signs, as R (not exp(x) - 1); a number times a logical column
    tiny <- xr * 1e-10
    expect_true(same(sum(expm1(tiny)), sum(expm1(d$x * 1e-10)), 1e-12))
    expect_true(same(sum(-expm1(-exp(xr))), sum(-expm1(-exp(d$x)))))
    expect_true(same(sum(2 * (xr > 0) + is.na(xr)), 2 * sum(d$x > 0)))
    expect_true(same(sum((xr > 0) + (xr < 1)), sum((d$x > 0) + (d$x < 1))))
    ## col_ifelse(missing =): a value where the test is missing
    tn <- col_ifelse(xr > 0, NA, xr) > -1
    expect_true(same(sum(col_ifelse(tn, 1, 0, missing = 0)), sum(d$x <= 0 & d$x > -1)))
    expect_true(same(sum(col_ifelse(tn, 1, 0, missing = 5)), sum(d$x <= 0 & d$x > -1) + 5 * sum(d$x > 0)))
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
    ## keys with missing values, bare codes on the stored values; is.na() of
    ## a key takes columns of any type, variable by variable
    kr <- design_key(R, "a")
    okr <- !is.na(kr)
    expect_equal(sum(okr), sum(ok))
    kb <- design_key(R, c("b", "a", "x"))
    nab <- is.na(kb)
    expect_identical(dim(nab)[2L], 3)
    expect_equal(sum(nab[, 1L]), 0)
    expect_equal(sum(nab[, 2L]), sum(is.na(d$a)))
    expect_equal(sum(nab[, 3L]), 0)
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
    ## (a grouping of design_groups() with bare codes)
    ggr <- design_groups(R, "b")
    dgg <- design_group_sums(R, list(ggr), w = wr)
    expect_true(same(dgg[[1L]]$n, rowsum(d$w, d$b)))
    lk <- group_lookup(seq_len(gr$G) * 10, gr)
    expect_true(same(sum(lk), sum(match(d$a[ok], ref) * 10)))
    ## keep = NULL remotely: the columns on the design returned
    alr <- design_subset(R, cols = list(y = yr, none = NULL))
    expect_true(same(sum(alr$cols$y), sum(d$y)))
    expect_null(alr$cols$none)

    ## ---- design_frame() on remote data ----
    d$lg <- d$x > 0
    d$dt <- as.Date("2020-01-01") + seq_len(n)
    d$s1 <- rbinom(n, 5, 0.4)
    d$s2 <- 5 - d$s1
    d$fz <- factor(sample(c("no", "yes", "never"), n, TRUE), levels = c("no", "yes", "never"))
    d$fz[d$fz == "never"] <- "no"
    DBI::dbWriteTable(con, "d2", d)
    myf <- function(v) v * 3
    fr <- design_frame(y ~ x + f, con, table = "d2", weights = 2 * w, subset = x > 0)
    expect_inherits(fr, "xtx_remote_frame")
    expect_inherits(fr$x, "xtx_remote")
    expect_identical(dim(fr$x), c(as.numeric(sum(keep)), ncol(D)))
    expect_true(same(sum(fr$weights), 2 * sum(d$w[keep])))
    expect_identical(fr$response$kind, "numeric")
    expect_identical(fr$response$ncol, 1L)
    expect_true(same(sum(fr$response$y), sum(d$y[keep])))
    expect_null(fr$offset)
    expect_identical(unname(fr$types[c("y", "b", "lg", "dt")]), c("double", "character", "logical", "other"))
    expect_identical(fr$na.action, "omit")
    expect_identical(fr$xlevels, list(f = levels(droplevels(d$f[keep]))))
    expect_identical(attr(fr$terms, "term.labels"), c("x", "f"))
    expect_true(any(grepl("remote model frame", capture.output(print(fr)))))
    ## the same model by the columns of the frame, against lm()
    sw <- sqrt(fr$weights)
    b <- solve(crossprod(sw * fr$x), crossprod(sw * fr$x, sw * fr$response$y))
    expect_true(same(b, coef(lm(y ~ x + f, d, weights = 2 * w, subset = x > 0)), 1e-8))
    ## `data` by name, in any position; no `data`, and data in memory: a model frame
    fr1 <- design_frame(table = "d2", data = con, formula = y ~ x)
    expect_inherits(fr1, "xtx_remote_frame")
    expect_identical(design_frame(y ~ x, d), model.frame(y ~ x, d))
    expect_identical(design_frame(data = d, y ~ x, weights = w), model.frame(data = d, y ~ x, weights = w))
    yy <- d$y; xx <- d$x
    expect_identical(design_frame(yy ~ xx), model.frame(yy ~ xx))
    ## an argument evaluated in R, an offset() term and an `offset`: one sum
    expect_message(fr2 <- design_frame(y ~ x + offset(w), con, table = "d2", weights = myf(w), offset = x),
                   "Evaluated in R: weights `myf(w)`", fixed = TRUE)
    expect_true(same(sum(fr2$weights), 3 * sum(d$w)))
    expect_true(same(sum(fr2$offset), sum(d$w + d$x)))
    expect_error(design_frame(y ~ x, con, table = "d2", weights = myf(w), engine = "sql"), "SQL engine cannot")
    ## weights and offset evaluated in R that give matrices: their numbers of columns
    fr3 <- suppressMessages(design_frame(y ~ x, con, table = "d2", weights = cbind(w, w), offset = cbind(x, x, x)))
    expect_identical(fr3$ncol, c(weights = 2L, offset = 3L))
    expect_identical(fr$ncol, c(weights = 1L, offset = 1L))
    ## extra columns that must be present (bare names); unknown arguments
    fr4 <- design_frame(y ~ x, con, table = "d2", a, b)
    expect_identical(dim(fr4$x)[1L], as.numeric(sum(!is.na(d$a))))
    expect_error(design_frame(y ~ x, con, table = "d2", nosuch), "No column `nosuch` in the remote data.", fixed = TRUE)
    expect_error(design_frame(y ~ x, con, table = "d2", foo = 1), "Unknown arguments: foo", fixed = TRUE)
    expect_error(design_frame(y ~ x, con), "exactly one")
    ## na.action
    expect_identical(design_frame(y ~ x + a, con, table = "d2", na.action = na.exclude)$na.action, "exclude")
    expect_error(design_frame(y ~ x + a, con, table = "d2", na.action = "na.fail"), "missing values in object")
    expect_error(design_frame(y ~ x + a, con, table = "d2", na.action = na.pass), "NA/NaN/Inf in 'x'")
    expect_error(design_frame(y ~ x, con, table = "d2", na.action = function(x) x), "`na.action` must be")
    ## columns of the wrong type as weights and offset: the errors of lm()
    expect_error(design_frame(y ~ x, con, table = "d2", weights = b), "'weights' must be a numeric vector", fixed = TRUE)
    expect_error(design_frame(y ~ x, con, table = "d2", weights = lg), "'weights' must be a numeric vector", fixed = TRUE)
    expect_error(design_frame(y ~ x, con, table = "d2", offset = dt), "'offset' must be numeric", fixed = TRUE)
    expect_error(design_frame(y ~ x + offset(b), con, table = "d2"), "non-numeric argument to binary operator",
                 fixed = TRUE)
    expect_error(design_frame(y ~ x, con, table = "d2", subset = x > 100), "0 (non-NA) cases", fixed = TRUE)
    ## the response, of any type
    r0 <- design_frame(~ x, con, table = "d2")$response
    expect_identical(r0$kind, "none")
    rl <- design_frame(lg ~ f, con, table = "d2")$response
    expect_identical(rl$kind, "logical")
    expect_equal(sum(rl$y), sum(d$lg))
    rf <- design_frame(fz ~ x, con, table = "d2", subset = x > 0)$response
    expect_identical(rf$kind, "factor")
    expect_identical(rf$levels, c("no", "yes"))
    expect_equal(sum(col_ifelse(rf$y == "yes", 1, 0)), sum(d$fz[keep] == "yes"))
    expect_equal(sum(col_ifelse(rf$y != "yes", 1, 0)), sum(d$fz[keep] != "yes"))
    expect_error(rf$y == c("a", "b"), "one string")
    rc <- design_frame(b ~ x, con, table = "d2")$response
    expect_identical(rc$kind, "factor")
    expect_identical(rc$levels, sort(unique(d$b)))
    ro <- design_frame(dt ~ x, con, table = "d2")$response
    expect_identical(ro[c("kind", "ncol")], list(kind = "other", ncol = 1L))
    expect_null(ro$y)
    rm2 <- design_frame(cbind(s1, s2) ~ x, con, table = "d2")$response
    expect_identical(rm2[c("kind", "ncol")], list(kind = "matrix", ncol = 2L))
    expect_equal(sum(rm2$cols[[1L]] + rm2$cols[[2L]]), 5 * n)
    rm3 <- design_frame(cbind(s1, s2, a) ~ x, con, table = "d2")
    expect_identical(rm3$response$ncol, 3L)
    expect_identical(length(rm3$response$cols), 2L)
    expect_identical(dim(rm3$x)[1L], as.numeric(sum(!is.na(d$a))))
    ## cbind() of one column is that column
    expect_identical(design_frame(cbind(lg) ~ x, con, table = "d2")$response$kind, "logical")
    expect_identical(design_frame(cbind(y) ~ x, con, table = "d2")$response$kind, "numeric")
    ## (but of a factor column: a matrix of one column, its codes, as in memory)
    rz <- suppressMessages(design_frame(cbind(fz) ~ x, con, table = "d2"))$response
    expect_identical(rz[c("kind", "ncol")], list(kind = "matrix", ncol = 1L))
    expect_equal(sum(rz$cols[[1L]]), sum(as.integer(d$fz)))
    expect_identical(design_frame(cbind() ~ x, con, table = "d2")$response$cols, list())
    ## a response evaluated in R that gives a matrix: its number of columns
    rk <- suppressMessages(design_frame(as.matrix(data.frame(y, x)) ~ f, con, table = "d2"))$response
    expect_identical(rk[c("kind", "ncol")], list(kind = "numeric", ncol = 2L))
    ## Parquet paths as data
    pq <- tempfile(fileext = ".parquet")
    DBI::dbExecute(con, sprintf("COPY d2 TO '%s' (FORMAT PARQUET)", pq))
    fp <- design_frame(y ~ x + f, pq, weights = w)
    expect_inherits(fp, "xtx_remote_frame")
    expect_true(same(sum(fp$weights), sum(d$w)))
    unlink(pq)

    ## ---- design_resolve(response =): the response on an existing design ----
    R2 <- design_matrix(~ x, con, table = "d2")
    rr <- design_resolve(R2, list(weights = quote(w)), environment(), response = quote(cbind(s1, s2)))
    expect_identical(rr$response$kind, "matrix")
    expect_equal(sum(rr$response$cols[[2L]]), sum(d$s2))
    expect_true(same(sum(rr$cols$weights), sum(d$w)))
    expect_identical(names(rr$cols), "weights")
    rr <- design_resolve(R2, list(), environment(), response = quote(fz))
    expect_identical(rr$response$levels, c("no", "yes"))
    expect_error(design_resolve(D, list(), response = quote(y)), "remote design")
    test_con_done(con)
}

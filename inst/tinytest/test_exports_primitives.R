## The exported functions that fitting code builds on (the fits themselves
## are in xtx.fit, which tests its own API): design_resolve(),
## is_remote_source(), design_patterns(fit =, check =), design_output(),
## col_const(), cbind() and rowSums() of columns, col_compute(table =
## FALSE), col_sums() with one na.rm per argument, design_groups(levels =).
## In memory and on DuckDB, against base R. Only exported functions.
source("helpers.R")
library(xtx)

same <- function(a, b, tol = 1e-8) isTRUE(all.equal(as.vector(a), as.vector(b), tolerance = tol))

set.seed(71)
n <- 300
d <- data.frame(id = seq_len(n), y = rnorm(n), x = rnorm(n), x2 = rnorm(n) + 5, w = runif(n, 0.5, 2),
                f = factor(sample(letters[1:4], n, TRUE)), a = sample(1:15, n, TRUE),
                b = sample(1:6, n, TRUE), stringsAsFactors = FALSE)
d$yb <- rbinom(n, 1, plogis(0.4 * d$x))
d$yp <- rpois(n, exp(0.3 * d$x))
## a singleton level of a
d$a[1] <- 99L

## ---- is_remote_source(), design_resolve() ----
expect_false(is_remote_source(d))
expect_false(is_remote_source(1))
M <- design_matrix(y ~ x + f, d)
expect_equal(design_resolve(M, list(w = quote(2 * w), z = NULL), environment(), data = d)$cols$w, 2 * d$w)

## ---- design_patterns(fit =, check = FALSE) ----
Mp <- design_matrix(y ~ f, d)
p1 <- design_patterns(Mp, d$y, compress = TRUE)
p2 <- design_patterns(Mp, d$y, compress = TRUE, check = FALSE)
p3 <- design_patterns(Mp, d$y, compress = TRUE, fit = "glm")
expect_equal(p1, p2)
expect_equal(p1, p3)
expect_error(design_patterns(Mp, d$y, fit = "x"))

## ---- columns: col_const(), cbind(), rowSums(), col_compute(table = FALSE) ----
expect_identical(col_const(M, 2), rep(2, n))
expect_identical(col_const(d$y, 0), rep(0, n))
expect_error(col_const(M, 1:2))
expect_identical(col_compute(d$y, table = FALSE), d$y)
expect_error(col_sums(d$y, na.rm = c(TRUE, FALSE)), "one value per argument")
## one na.rm per argument
yn <- d$y
yn[3] <- NA
expect_equal(unname(col_sums(yn, yn, na.rm = c(FALSE, TRUE))), c(NA, sum(yn, na.rm = TRUE)))
expect_equal(unname(col_sums(a = yn, b = d$y, na.rm = c(TRUE, FALSE))), c(sum(yn, na.rm = TRUE), sum(d$y)))

## ---- design_output() ----
ref <- lm(y ~ x + f, d, weights = w)
f <- drop(M %*% coef(ref))
out <- design_output(M, list(.fitted = f, .resid = d$y - f))
expect_true(is.data.frame(out))
expect_equal(out$.fitted, f)
expect_identical(design_output(M, list(.fitted = f), memory = function(dd) dd$.fitted), f)

## ---- design_groups(levels =): new rows coded against training levels ----
g <- design_groups(M, d["a"])
nd <- d[1:5, ]
nd$a <- c(2L, 3L, 1000L, 2L, NA)
Mn <- design_matrix(y ~ x + f, nd)
gn <- design_groups(Mn, nd["a"], levels = g$levels)
expect_identical(gn$codes, match(nd$a, g$levels$a))
expect_equal(group_lookup(seq_len(g$G) * 10, gn), match(nd$a, g$levels$a) * 10)
expect_error(group_sums(1, gn), "no level")

## ---- remote: the same on DuckDB ----
if (requireNamespace("duckdb", quietly = TRUE)) {
    con <- test_con()
    DBI::dbWriteTable(con, "d", d)
    expect_true(is_remote_source(con))
    R <- design_matrix(y ~ x + f, con, table = "d")
    yr <- design_column(R, "y")
    wr <- design_column(R, "w")
    ## columns
    expect_equal(sum(col_const(R, 2)), 2 * n)
    cb <- cbind(yr, 1, wr)
    expect_equal(dim(cb), c(n, 3L))
    expect_true(same(sum(rowSums(cb)), sum(d$y + 1 + d$w)))
    pc <- col_compute(yr * 2, table = FALSE)
    expect_true(same(sum(pc), 2 * sum(d$y)))
    expect_error(col_compute(yr, name = "t", table = FALSE), "only with")
    ynr <- col_ifelse(design_column(R, "id") == 3, NA, yr)
    expect_equal(unname(col_sums(ynr, ynr, yr, na.rm = c(FALSE, TRUE, FALSE))), c(NA, sum(yn, na.rm = TRUE), sum(d$y)))
    ## NaN without NA gives NaN, as in R; a count with NA rows (never NULL
    ## after !is.na()) counts them as FALSE
    expect_true(is.nan(col_sums(col_ifelse(design_column(R, "id") == 3, NaN, yr))))
    cond <- ynr > 0
    expect_equal(unname(col_sums(cond & !is.na(cond))), sum(yn > 0, na.rm = TRUE))
    expect_true(is.na(col_sums(cond)))
    ## design_patterns(check = FALSE)
    Rp <- design_matrix(y ~ f, con, table = "d")
    q1 <- design_patterns(Rp, "y", compress = TRUE)
    q2 <- design_patterns(Rp, design_column(Rp, "y"), compress = TRUE, check = FALSE)
    expect_true(same(sort(q1$z), sort(q2$z)))
    expect_true(same(sort(q1$z), sort(p1$z)))
    ## design_output(): a remote table
    fo <- R %*% coef(lm(y ~ x + f, d))
    tab <- design_output(R, list(.fitted = fo))
    expect_inherits(tab, "xtx_remote_table")
    got <- DBI::dbGetQuery(con, sprintf("SELECT id, \".fitted\" AS v FROM %s ORDER BY id", tab))
    expect_true(same(got$v, fitted(lm(y ~ x + f, d))))
    ## design_resolve(): SQL, and an expression evaluated in R
    myf <- function(v) v * 3
    rr <- design_resolve(R, list(a = quote(2 * w), b = quote(myf(w))), environment())
    expect_true(same(sum(rr$cols$a), 2 * sum(d$w)))
    expect_true(same(sum(rr$cols$b), 3 * sum(d$w)))
    expect_identical(names(rr$evaluated), "b")
    ## design_groups(levels =) remotely
    gr <- design_groups(R, "a")
    DBI::dbWriteTable(con, "nd", nd)
    Rnd <- design_matrix(y ~ x + f, con, table = "nd")
    gnr <- design_groups(Rnd, "a", levels = gr$levels)
    lk <- as.data.frame(group_lookup(seq_len(gr$G) * 10, gnr), keep = "id")
    lk <- lk[order(lk$id), ]
    expect_equal(lk[[2L]], match(nd$a, g$levels$a)[order(nd$id)] * 10)
    test_con_done(con)
}

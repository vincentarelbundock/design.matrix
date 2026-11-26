## Public groupings (Kata design.matrix#0ncx): design_groups(), group_sums(),
## group_crosstab(), group_lookup() and group_nested(), in memory and on
## DuckDB with the same results; multi-key groupings, weights, missing keys,
## errors, and a one-way fixed-effect slope (within transformation through
## group_sums()/group_lookup()) that reproduces lm() with factor dummies on
## both backends. Only exported functions.
source("helpers.R")
library(design.matrix)

same <- function(a, b, tol = 1e-8) isTRUE(all.equal(as.vector(a), as.vector(b), tolerance = tol))

set.seed(59)
n <- 400
d <- data.frame(id = seq_len(n), y = rnorm(n), x = rnorm(n), w = runif(n, 0.5, 2),
                ## mixed case: the levels are sorted in the C locale ("C" first)
                firm = sample(c("b", "a", "C", "d", "e"), n, TRUE),
                year = sample(2001:2004, n, TRUE),
                flag = sample(c(TRUE, FALSE), n, TRUE),
                stringsAsFactors = FALSE)
d$region <- ifelse(d$firm %in% c("a", "b"), "north", "south")
## a column and a weight with missing values (not in the formula: no row
## of the design is dropped)
d$yn <- d$y
d$yn[c(1, 7)] <- NA
d$wn <- d$w
d$wn[2] <- NA

fo <- y ~ x
M <- design_matrix(fo, d)
firms <- sort(unique(d$firm), method = "radix")

## ---- in memory ----
g <- design_groups(M, d["firm"])
expect_true(inherits(g, "design_groups"))
expect_identical(g$G, 5L)
expect_identical(firms[1L], "C")
expect_identical(g$levels$firm, firms)
expect_identical(names(g$levels), "firm")
expect_true(is.integer(g$codes))
expect_identical(g$codes, match(d$firm, firms))
expect_identical(g$n, nrow(M))
expect_true(grepl("design_groups: 5 groups of firm", format(g)[1L], fixed = TRUE))
expect_true(grepl("in memory", format(g)[1L], fixed = TRUE))

## other forms of the key: a list, a vector, a variable of the model frame,
## a scaled design
expect_identical(design_groups(M, list(firm = d$firm))$codes, g$codes)
expect_identical(design_groups(2 * M, d["firm"])$codes, g$codes)
gy <- design_groups(M, d$year)
expect_identical(names(gy$levels), "key")
expect_identical(gy$levels$key, sort(unique(d$year)))
Mf <- design_matrix(y ~ x + firm, d)
gf <- design_groups(Mf, "firm")
expect_identical(gf$codes, g$codes)
expect_identical(gf$levels, g$levels)
## logical keys: integer levels 0 and 1
gl <- design_groups(M, d["flag"])
expect_identical(gl$levels$flag, c(0L, 1L))
expect_identical(gl$codes, as.integer(d$flag) + 1L)
## factors in memory: in the order of their levels
gfa <- design_groups(M, data.frame(f = factor(d$firm, levels = rev(firms))))
expect_identical(gfa$levels$f, rev(firms))

## group_sums: rowsum() in code order
expect_true(same(group_sums(d$y, g), rowsum(d$y, d$firm)[firms, 1]))
expect_identical(dim(group_sums(d$y, g)), c(5L, 1L))
expect_null(dimnames(group_sums(d$y, g)))
expect_true(same(group_sums(1, g), as.vector(table(d$firm)[firms])))
expect_true(same(group_sums(cbind(d$y, d$x), g), rowsum(cbind(d$y, d$x), d$firm)[firms, ]))
expect_true(same(group_sums(d[c("y", "x")], g), group_sums(cbind(d$y, d$x), g)))
expect_true(same(group_sums(list(d$y, d$x), g), group_sums(cbind(d$y, d$x), g)))
expect_true(same(group_sums("y", g), group_sums(d$y, g)))
expect_true(same(group_sums(c("y", "x"), g), group_sums(cbind(d$y, d$x), g)))
expect_true(same(group_sums(d$flag, g), rowsum(as.double(d$flag), d$firm)[firms, 1]))
## a missing value: NA for its group only (as rowsum())
sn <- group_sums(d$yn, g)
bad_g <- unique(g$codes[c(1, 7)])
expect_true(all(is.na(sn[bad_g, 1])))
expect_true(all(!is.na(sn[-bad_g, 1])))

## group_lookup: values[codes]
vals <- c(10, 20, 30, 40, 50)
expect_identical(group_lookup(vals, g), vals[g$codes])
expect_identical(group_lookup(1:5, g), as.double(1:5)[g$codes])
V2 <- cbind(vals, -vals)
expect_identical(group_lookup(V2, g), unname(V2[g$codes, ]))
expect_identical(group_lookup(c(1, NA, 3, 4, 5), g), c(1, NA, 3, 4, 5)[g$codes])

## group_crosstab
g2 <- design_groups(M, d["year"])
ct <- group_crosstab(g, g2)
expect_true(inherits(ct, "dgCMatrix"))
expect_identical(dim(ct), c(5L, 4L))
tab <- table(factor(d$firm, levels = firms), d$year)
expect_true(same(as.matrix(ct), unclass(tab)))
wt <- tapply(d$w, list(factor(d$firm, levels = firms), d$year), sum)
wt[is.na(wt)] <- 0
expect_true(same(as.matrix(group_crosstab(g, g2, w = d$w)), wt))
expect_true(same(as.matrix(group_crosstab(g, g2, w = "x")), {
    xt <- tapply(d$x, list(factor(d$firm, levels = firms), d$year), sum)
    xt[is.na(xt)] <- 0
    xt
}))
## a missing weight: NA for its pair only
ctn <- as.matrix(group_crosstab(g, g2, w = d$wn))
expect_true(is.na(ctn[g$codes[2], g2$codes[2]]))
expect_identical(sum(is.na(ctn)), 1L)
## a grouping with itself: the counts on the diagonal
expect_true(same(as.matrix(group_crosstab(g, g)), diag(as.vector(table(d$firm)[firms]))))

## group_nested
greg <- design_groups(M, d["region"])
expect_true(group_nested(g, greg))
expect_false(group_nested(greg, g))
expect_true(group_nested(g, g))
expect_false(group_nested(g, g2))

## multi-key groupings: lexicographic, by the first variable then the next
gm <- design_groups(M, d[c("region", "firm")])
u <- unique(d[c("region", "firm")])
u <- u[order(u$region, u$firm, method = "radix"), ]
rownames(u) <- NULL
expect_identical(gm$G, nrow(u))
expect_equal(gm$levels, u)
expect_identical(gm$codes, match(paste(d$region, d$firm), paste(u$region, u$firm)))
expect_true(group_nested(g, gm) && group_nested(gm, g))
gmy <- design_groups(M, list(firm = d$firm, year = d$year))
expect_identical(names(gmy$levels), c("firm", "year"))
expect_identical(gmy$G, nrow(unique(d[c("firm", "year")])))
expect_true(group_nested(gmy, g))

## missing keys: an error
fk <- d$firm
fk[5] <- NA
expect_error(design_groups(M, data.frame(firm = fk)), "missing values")
expect_error(design_groups(M, c(NaN, d$x[-1])), "missing values")

## errors
expect_error(design_groups(as.matrix(M), d["firm"]), "must be a design")
expect_error(design_groups(M, "nope"), "not a variable of the design")
expect_error(design_groups(M, d$year[-1]), "one per row")
expect_error(design_groups(M, list(m = matrix(1, n, 2))), "must be a vector")
expect_error(group_sums(d$y[-1], g), "one value \\(row\\) per row")
expect_error(group_sums(letters[1:3], g), "not a variable")
expect_error(group_sums(d$y, list(codes = g$codes, G = 5L)), "made by design_groups")
expect_error(group_lookup(1:4, g), "one value \\(row\\) per group")
expect_error(group_lookup("a", g), "numeric")
expect_error(group_crosstab(g, g2, w = cbind(d$w, d$w)), "one column")
bad <- g
bad$codes[1] <- 6L
expect_error(group_sums(d$y, bad), "within 1..G")
expect_error(group_lookup(vals, bad), "within 1..G")
expect_error(group_crosstab(bad, g2), "within 1..G")
expect_error(group_nested(g2, bad), "within 1..G")
bad$codes[1] <- NA_integer_
expect_error(group_sums(d$y, bad), "within 1..G")
bad <- g
bad$codes <- as.double(g$codes)
expect_error(group_sums(d$y, bad), "integer codes")
bad <- g
bad$G <- 4L
expect_error(group_sums(d$y, bad), "not a valid grouping")
## groupings of different designs
Mh <- design_matrix(fo, d[1:200, ])
gh <- design_groups(Mh, d[1:200, "firm", drop = FALSE])
expect_error(group_crosstab(g, gh), "different designs")
expect_error(group_nested(g, gh), "different designs")
M3 <- design_matrix(fo, d)
g3 <- design_groups(M3, d["year"])
expect_error(group_crosstab(g, g3), "different designs")
## (a scaled design is the same design)
expect_true(same(as.matrix(group_crosstab(g, design_groups(2 * M, d["year"]))), as.matrix(ct)))

## one-way fixed effect: the within slope, the same code on both backends
within_slope <- function(D, g, w = NULL) {
    if (is.null(w)) {
        n_g <- group_sums(1, g)[, 1]
        demean <- function(v) v - group_lookup(group_sums(v, g)[, 1] / n_g, g)
        yt <- demean(design_column(D, "y"))
        xt <- demean(design_column(D, "x"))
        return(sum(as.vector(xt * yt)) / sum(as.vector(xt * xt)))
    }
    sw <- group_sums(w, g)[, 1]
    demean <- function(v) v - group_lookup(group_sums(w * v, g)[, 1] / sw, g)
    yt <- demean(design_column(D, "y"))
    xt <- demean(design_column(D, "x"))
    sum(as.vector(w * xt * yt)) / sum(as.vector(w * xt * xt))
}
ref1 <- coef(lm(y ~ x + firm, d))[["x"]]
ref2 <- coef(lm(y ~ x + interaction(firm, year), d))[["x"]]
refw <- coef(lm(y ~ x + firm, d, weights = w))[["x"]]
expect_equal(within_slope(M, g), ref1, tolerance = 1e-10)
expect_equal(within_slope(M, gmy), ref2, tolerance = 1e-10)
expect_equal(within_slope(M, g, d$w), refw, tolerance = 1e-10)

if (!requireNamespace("duckdb", quietly = TRUE) || !requireNamespace("DBI", quietly = TRUE))
    exit_file("duckdb is not installed")

expect_error(design_groups(M, DBI::SQL("firm")), "only of a remote design")

## ---- on DuckDB: the same results ----
con <- test_con()
DBI::dbWriteTable(con, "d", d)
R <- design_matrix(fo, con, table = "d")

## values of a remote column in the order of the rows of d (by the key `id`)
by_id <- function(col) {
    r <- as.data.frame(col, keep = "id")
    r <- r[order(r$id), , drop = FALSE]
    expect_identical(r$id, d$id)
    m <- as.matrix(r[setdiff(names(r), "id")])
    if (ncol(m) == 1L) unname(m[, 1L]) else unname(m)
}

gr <- design_groups(R, "firm")
expect_true(inherits(gr, "design_groups"))
expect_true(inherits(gr$codes, "design_column"))
expect_identical(gr$G, g$G)
expect_equal(gr$levels, g$levels)
expect_identical(as.integer(by_id(gr$codes)), g$codes)
expect_identical(as.integer(sort(as.vector(gr$codes))), sort(g$codes))
expect_true(grepl("design_groups: 5 groups of firm", format(gr)[1L], fixed = TRUE))
expect_true(grepl("duckdb", format(gr)[1L], fixed = TRUE))
expect_identical(design_groups(2 * R, "firm")$levels, gr$levels)

## group_sums
expect_true(same(group_sums("y", gr), group_sums(d$y, g)))
expect_true(same(group_sums(design_column(R, "y"), gr), group_sums(d$y, g)))
expect_true(same(group_sums(DBI::SQL("y * 2"), gr), 2 * group_sums(d$y, g)))
expect_identical(group_sums(1, gr), group_sums(1, g))
expect_true(same(group_sums(c("y", "x"), gr), group_sums(cbind(d$y, d$x), g)))
expect_true(same(group_sums(list(design_column(R, "y"), design_column(R, "x")), gr),
                 group_sums(cbind(d$y, d$x), g)))
expect_identical(dim(group_sums(c("y", "x"), gr)), c(5L, 2L))
expect_null(dimnames(group_sums("y", gr)))
expect_true(same(group_sums("yn", gr), group_sums(d$yn, g)))
expect_true(all(is.na(group_sums("yn", gr)[bad_g, 1])))
## a product of the design (pending), and a column computed first
b <- c(0.5, -1)
expect_true(same(group_sums(R %*% b, gr), group_sums(drop(M %*% b), g)))
expect_true(same(group_sums(compute(design_column(R, "y") - R %*% b), gr),
                 group_sums(d$y - drop(M %*% b), g)))

## group_lookup
lk <- group_lookup(vals, gr)
expect_true(inherits(lk, "design_column"))
expect_identical(by_id(lk), vals[g$codes])
expect_identical(by_id(group_lookup(V2, gr)), unname(V2[g$codes, ]))
expect_identical(by_id(group_lookup(c(1, NA, 3, 4, 5), gr)), c(1, NA, 3, 4, 5)[g$codes])

## group_crosstab
gr2 <- design_groups(R, "year")
expect_equal(gr2$levels, g2$levels)
ctr <- group_crosstab(gr, gr2)
expect_true(inherits(ctr, "dgCMatrix"))
expect_identical(dim(ctr), dim(ct))
expect_true(same(as.matrix(ctr), as.matrix(ct)))
expect_true(same(as.matrix(group_crosstab(gr, gr2, w = "w")), as.matrix(group_crosstab(g, g2, w = d$w))))
expect_true(same(as.matrix(group_crosstab(gr, gr2, w = design_column(R, "w"))), wt))
ctrn <- as.matrix(group_crosstab(gr, gr2, w = "wn"))
expect_true(same(ctrn, ctn))
expect_identical(sum(is.na(ctrn)), 1L)
expect_true(same(as.matrix(group_crosstab(gr, gr)), as.matrix(group_crosstab(g, g))))

## group_nested
grr <- design_groups(R, "region")
expect_true(group_nested(gr, grr))
expect_false(group_nested(grr, gr))
expect_true(group_nested(gr, gr))
expect_false(group_nested(gr, gr2))

## multi-key groupings, keys from SQL() and from columns
grm <- design_groups(R, c("region", "firm"))
expect_equal(grm$levels, gm$levels)
expect_identical(as.integer(by_id(grm$codes)), gm$codes)
grmy <- design_groups(R, list("firm", year = DBI::SQL("year")))
expect_identical(names(grmy$levels), c("firm", "year"))
expect_equal(grmy$levels, gmy$levels)
expect_identical(as.integer(by_id(grmy$codes)), gmy$codes)
grs <- design_groups(R, DBI::SQL("year % 2"))
gms <- design_groups(M, list(key = d$year %% 2L))
expect_true(same(grs$levels[[1L]], gms$levels[[1L]]))
expect_identical(as.integer(by_id(grs$codes)), gms$codes)
grc <- design_groups(R, design_column(R, "x") > 0)
gmc <- design_groups(M, d$x > 0)
expect_identical(grc$levels$key, c(0L, 1L))
expect_identical(gmc$levels$key, c(0L, 1L))
expect_identical(as.integer(by_id(grc$codes)), gmc$codes)
grl <- design_groups(R, "flag")
expect_identical(grl$levels$flag, c(0L, 1L))
expect_identical(as.integer(by_id(grl$codes)), gl$codes)

## groupings and columns used in any order (each use of a grouping joins
## its codes on top of the other columns of the use)
p <- R %*% b
pv <- as.vector(p)
expect_true(same(group_sums(p, gr), group_sums(drop(M %*% b), g)))
expect_true(same(as.matrix(group_crosstab(gr, gr2)), as.matrix(ct)))
mix <- group_lookup(vals, gr) + group_lookup(c(1, 2, 3, 4), gr2) * design_column(R, "x")
expect_true(same(by_id(mix), vals[g$codes] + c(1, 2, 3, 4)[g2$codes] * d$x))
expect_true(same(group_sums(mix, grm), group_sums(vals[g$codes] + c(1, 2, 3, 4)[g2$codes] * d$x, gm)))

## the within slope on DuckDB
expect_equal(within_slope(R, gr), ref1, tolerance = 1e-8)
expect_equal(within_slope(R, grmy), ref2, tolerance = 1e-8)
expect_equal(within_slope(R, gr, design_column(R, "w")), refw, tolerance = 1e-8)

## missing keys: an error, and a row selection avoids it
dn <- d
dn$firm[c(3, 9)] <- NA
dn$xk <- dn$x
dn$xk[4] <- NaN
DBI::dbWriteTable(con, "dn", dn)
Rn <- design_matrix(fo, con, table = "dn")
expect_error(design_groups(Rn, "firm"), "missing values")
expect_error(design_groups(Rn, list("region", "firm")), "variable `firm`")
expect_error(design_groups(Rn, DBI::SQL("CAST(xk AS DOUBLE)")), "missing values")
Rk <- Rn[DBI::SQL("firm IS NOT NULL"), ]
gk <- design_groups(Rk, "firm")
expect_identical(gk$G, 5L)
expect_true(same(group_sums(1, gk), as.vector(table(dn$firm[!is.na(dn$firm)])[firms])))

## errors
expect_error(design_groups(R, d$year), "remote rows have no order")
expect_error(design_groups(R, "nope"), "nope")
expect_error(group_sums(d$y, gr), "remote rows have no order")
expect_error(group_sums(design_column(R, "y"), g), "remote column, and the grouping is in memory")
expect_error(design_groups(M, design_column(R, "x")), "remote column")
expect_error(group_crosstab(g, gr), "in memory and the other remote")
expect_error(group_nested(gr, g), "in memory and the other remote")
expect_error(group_lookup(1:3, gr), "one value \\(row\\) per group")
expect_error(group_crosstab(gr, gr2, w = d$w), "remote rows have no order")
R2 <- design_matrix(fo, con, table = "d")
gr3 <- design_groups(R2, "year")
expect_error(group_crosstab(gr, gr3), "different designs")
expect_error(group_nested(gr3, gr), "different designs")
expect_error(group_sums(design_column(R2, "y"), gr), "another remote design")
expect_error(design_groups(R, design_column(R2, "x") > 0), "another remote design")
expect_error(group_sums(DBI::SQL("no_such_column + 1"), gr), "invalid SQL")

test_con_done(con)

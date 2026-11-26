## R's NA/NaN semantics on remote columns (design.matrix#p701): every
## operation on columns holding NA, NaN, Inf, -Inf, 0 and ordinary numbers
## gives, in DuckDB, exactly what it gives in memory (NA where R gives NA,
## NaN where R gives NaN). In DuckDB NaN = NaN and NaN > 0 are true, NULL is
## NA, and is.na() on a remote column used to be per list element.
if (!requireNamespace("duckdb", quietly = TRUE) || !requireNamespace("DBI", quietly = TRUE))
    exit_file("duckdb is not installed")
source("helpers.R")
library(design.matrix)

dmp_column <- design.matrix:::dmp_column
dmp_collect <- design.matrix:::dmp_collect
dmc_cbind <- design.matrix:::dmc_cbind
dmc_ifelse <- design.matrix:::dmc_ifelse
dmc_count <- design.matrix:::dmc_count
dmc_notnan <- design.matrix:::dmc_notnan
dmc_is_finite <- design.matrix:::dmc_is_finite
dmc_key <- design.matrix:::dmc_key
dmc_not_na <- design.matrix:::dmc_not_na
dm_colsum <- design.matrix:::dm_colsum

vals <- c(NA, NaN, Inf, -Inf, 0, 1, -2.5)
g <- expand.grid(a = vals, b = vals)
d <- data.frame(id = seq_len(nrow(g)), y = 1, a = g$a, b = g$b,
                f = ifelse(is.na(g$a) & !is.nan(g$a), NA, c("u", "v")[1 + (seq_len(nrow(g)) %% 2)]),
                i = ifelse(is.nan(g$b), NA, seq_len(nrow(g))), stringsAsFactors = FALSE)
con <- test_con()
DBI::dbWriteTable(con, "d", d)
R <- design_matrix(y ~ id, con, table = "d")
col <- function(v) dmp_column(R, sprintf('CAST(s."%s" AS DOUBLE)', v))
id <- col("id")
a <- col("a")
b <- col("b")

## a remote column in the order of d's rows, as doubles (logicals as 0/1)
get <- function(x) {
    m <- dmp_collect(dmc_cbind(id, x))
    m[order(m[, 1]), 2]
}
## the data went in with NA and NaN distinct
expect_identical(get(a), d$a)
expect_identical(get(b), d$b)

check <- function(remote, local, what) {
    expect_true(inherits(remote, "design_column"), info = what)
    expect_identical(get(remote), as.double(local), info = what)
}

## is.na() and friends: never NA; NaN is NA but NULL is not NaN
check(is.na(a), is.na(d$a), "is.na")
check(is.nan(a), is.nan(d$a), "is.nan")
check(is.finite(a), is.finite(d$a), "is.finite")
check(is.infinite(a), is.infinite(d$a), "is.infinite")
check(!is.na(a), !is.na(d$a), "!is.na")
check(dmc_notnan(a), !is.na(d$a), "dmc_notnan")
check(dmc_is_finite(a), is.finite(d$a), "dmc_is_finite")
## on a multi-column column, one column each
ab <- dmc_cbind(a, b)
m <- dmp_collect(dmc_cbind(id, is.na(ab), is.infinite(ab)))
m <- m[order(m[, 1]), -1]
expect_identical(unname(m), unname(1 * cbind(is.na(d$a), is.na(d$b), is.infinite(d$a), is.infinite(d$b))))

## comparisons: NA when either side is NA or NaN
for (op in c("==", "!=", "<", "<=", ">", ">=")) {
    f <- match.fun(op)
    check(f(a, b), f(d$a, d$b), op)
    check(f(a, 0), f(d$a, 0), paste(op, "0"))
    check(f(0, b), f(0, d$b), paste("0", op))
    check(f(a, Inf), f(d$a, Inf), paste(op, "Inf"))
    check(f(a, NaN), f(d$a, NaN), paste(op, "NaN"))
    check(f(a, NA), f(d$a, NA), paste(op, "NA"))
    check(f(a, a), f(d$a, d$a), paste(op, "self"))
}

## logical operators: three-valued, as R
pa <- a > 0
pb <- b <= 0
check(pa & pb, d$a > 0 & d$b <= 0, "&")
check(pa | pb, d$a > 0 | d$b <= 0, "|")
check(!pa, !(d$a > 0), "!")
check(pa & FALSE, d$a > 0 & FALSE, "& FALSE")
check(pa | TRUE, d$a > 0 | TRUE, "| TRUE")
check(is.na(a) | pb, is.na(d$a) | d$b <= 0, "is.na |")

## arithmetic: NaN and NA propagate (where both sides are missing, R's
## result, NA or NaN, depends on the platform: only missingness compared)
both <- is.na(d$a) & is.na(d$b)
for (op in c("+", "-", "*", "/")) {
    f <- match.fun(op)
    r <- get(f(a, b))
    l <- f(d$a, d$b)
    expect_identical(r[!both], l[!both], info = op)
    expect_true(all(is.na(r[both])), info = op)
    check(f(a, 2), f(d$a, 2), paste(op, "2"))
    check(f(a, NaN), f(d$a, NaN), paste(op, "NaN"))
    check(f(a, 0), f(d$a, 0), paste(op, "0"))
    check(-a, -d$a, "unary -")
}
check(a^2, d$a^2, "^2")

## pmax/pmin (design.matrix#2tr6): NA or NaN where either side is missing,
## y's where both are (R's pmax takes the later missing argument), with
## nothing platform-dependent
dmc_pmax <- design.matrix:::dmc_pmax
dmc_pmin <- design.matrix:::dmc_pmin
check(dmc_pmax(a, b), pmax(d$a, d$b), "pmax")
check(dmc_pmin(a, b), pmin(d$a, d$b), "pmin")
check(dmc_pmax(b, a), pmax(d$b, d$a), "pmax swapped")
check(dmc_pmin(b, a), pmin(d$b, d$a), "pmin swapped")
for (k in list(0, 1, -2.5, Inf, -Inf, NA, NaN)) {
    lab <- format(k)
    check(dmc_pmax(a, k), pmax(d$a, k), paste("pmax", lab))
    check(dmc_pmin(a, k), pmin(d$a, k), paste("pmin", lab))
    check(dmc_pmax(k, a), pmax(k, d$a), paste(lab, "pmax"))
    check(dmc_pmin(k, a), pmin(k, d$a), paste(lab, "pmin"))
}
## nested, as the glm links clamp
check(dmc_pmin(dmc_pmax(a, -1), 1), pmin(pmax(d$a, -1), 1), "clamp")
check(dmc_pmax(0.001, dmc_pmin(0.999, a)), pmax(0.001, pmin(0.999, d$a)), "clamp constant first")
## on a multi-column column, one column each
m <- dmp_collect(dmc_cbind(id, dmc_pmax(ab, 0), dmc_pmin(ab, 0)))
m <- m[order(m[, 1]), -1]
expect_identical(unname(m), unname(cbind(pmax(d$a, 0), pmax(d$b, 0), pmin(d$a, 0), pmin(d$b, 0))))

## sign: NaN for NaN, NA for NA
check(sign(a), sign(d$a), "sign")
check(sign(-a), sign(-d$a), "sign of -a")

## ^ as R's R_pow: x^0 = 1 and 1^y = 1 for every x and y (NA and NaN too); a
## negative base with an infinite power, and -Inf with a power that is not
## an integer, give NaN (C's pow() gives a number). Where both sides are
## missing, R's result (x + y) depends on the platform: only missingness
## compared there.
checkm <- function(remote, local, miss, what) {
    expect_true(inherits(remote, "design_column"), info = what)
    r <- get(remote)
    expect_identical(r[!miss], as.double(local)[!miss], info = what)
    expect_true(all(is.na(r[miss])), info = what)
}
pm <- is.na(d$a) & is.na(d$b)
checkm(a^b, d$a^d$b, pm, "^")
checkm(b^a, d$b^d$a, pm, "^ swapped")
check(a^0, d$a^0, "^0")
check(a^1, d$a^1, "^1")
check(a^3, d$a^3, "^3")
check(a^-1, d$a^-1, "^-1")
check(a^0.5, d$a^0.5, "^0.5")
check(a^-2.5, d$a^-2.5, "^-2.5")
check(a^Inf, d$a^Inf, "^Inf")
check(a^-Inf, d$a^-Inf, "^-Inf")
checkm(a^NA, d$a^NA, is.nan(d$a), "^NA")
checkm(a^NaN, d$a^NaN, is.na(d$a) & !is.nan(d$a), "^NaN")
for (k in list(1, 0, -1, -2.5, 2, Inf, -Inf)) check(k^b, k^d$b, paste(format(k), "^"))
checkm(NA^b, NA^d$b, is.nan(d$b), "NA ^")
checkm(NaN^b, NaN^d$b, is.na(d$b) & !is.nan(d$b), "NaN ^")

## ifelse and counts: a NA test gives NA (ifelse) and does not count
check(dmc_ifelse(a > 0, a, 7), dmc_ifelse(d$a > 0, d$a, 7), "ifelse")
check(dmc_ifelse(a == b, 1, 0), dmc_ifelse(d$a == d$b, 1, 0), "ifelse ==")
check(dmc_count(a > b), dmc_count(d$a > d$b), "count")
check(dmc_count(a > b, negate = TRUE), dmc_count(d$a > d$b, negate = TRUE), "count negate")
expect_equal(dm_colsum(dmc_count(a != 0)), sum(dmc_count(d$a != 0)))
expect_equal(dm_colsum(dmc_count(is.na(a))), sum(is.na(d$a)))

## complete keys: NaN is missing in both (complete.cases); a text or
## integer key only NA
for (vars in list("a", c("a", "b"), "f", c("f", "i"), c("f", "b"))) {
    kl <- dmc_key(NULL, vars, data = d)
    kr <- dmc_key(R, vars)
    check(dmc_not_na(kr), dmc_not_na(kl), paste("not_na", paste(vars, collapse = ",")))
}

test_con_done(con)

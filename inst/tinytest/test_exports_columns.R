## Column arithmetic for package developers (Kata design.matrix#tapx): the
## Summary group, mean(), round/signif/trunc/gamma/lgamma on remote columns,
## and the exported helpers col_ifelse(), col_pmin(), col_pmax(),
## col_sums(). Each result on a DuckDB table is compared with base R on the
## same values in memory (where a column is a plain vector), over NA, NaN,
## Inf and -Inf. Only exported functions (and base R, DBI) are called.
if (!requireNamespace("duckdb", quietly = TRUE) || !requireNamespace("DBI", quietly = TRUE))
    exit_file("duckdb is not installed")
source("helpers.R")
library(design.matrix)

expect_true(all(c("col_ifelse", "col_pmin", "col_pmax", "col_sums") %in% getNamespaceExports("design.matrix")))

## Values whose sums and products are exact in any order, so that remote
## results are identical() to R's.
vals <- c(NA, NaN, Inf, -Inf, 0, 1, -2.5, 3)
g <- expand.grid(a = vals, b = vals)
d <- data.frame(id = seq_len(nrow(g)), y = 1, a = g$a, b = g$b, k = rep(1:4, length.out = nrow(g)))
con <- test_con()
DBI::dbWriteTable(con, "d", d)
R <- design_matrix(y ~ id, con, table = "d")
a <- design_column(R, "a")
b <- design_column(R, "b")

## a remote column (one column) in the order of d's rows
get <- function(x, id = "id") {
    v <- as.data.frame(x, keep = id, nm = ".v")
    v$.v[order(v[[id]])]
}
## the data went in with NA and NaN distinct
expect_identical(get(a), d$a)
expect_identical(get(b), d$b)

## value and warnings of f(), for comparing remote and in-memory calls
run <- function(f) capture(f)
same_run <- function(r, l, info, loose = FALSE) {
    if (loose && length(l$value) && all(is.na(l$value))) {
        ## (R's result for a mix of NA and NaN depends on the platform)
        expect_true(all(is.na(r$value)), info = info)
    } else {
        expect_identical(r$value, l$value, info = info)
    }
    expect_identical(r$warnings, l$warnings, info = info)
}

## ---- Summary group and mean, over subsets of the NA/NaN/Inf grid ----

## the rows of each case, as a SQL predicate and in R
cases <- list(
    finite = list("isfinite(a)", function(v) is.finite(v)),
    inf = list("NOT isnan(a)", function(v) !is.na(v)),
    nan = list("a IS NOT NULL", function(v) !(is.na(v) & !is.nan(v))),
    na = list("a IS NULL OR NOT isnan(a)", function(v) !is.nan(v)),
    missing = list("a IS NULL OR isnan(a)", function(v) is.na(v)),
    all = list("TRUE", function(v) rep(TRUE, length(v))))

for (nm in names(cases)) {
    cs <- cases[[nm]]
    Rs <- R[DBI::SQL(cs[[1L]]), ]
    ar <- design_column(Rs, "a")
    br <- design_column(Rs, "b")
    keep <- cs[[2L]](d$a)
    al <- d$a[keep]
    bl <- d$b[keep]
    mixed <- anyNA(al) && any(is.nan(al)) || anyNA(bl) && any(is.nan(bl))
    for (na.rm in c(FALSE, TRUE)) {
        info <- function(f) sprintf("%s(), case %s, na.rm = %s", f, nm, na.rm)
        for (f in c("sum", "prod", "mean")) {
            fun <- match.fun(f)
            same_run(run(function() fun(ar, na.rm = na.rm)), run(function() fun(al, na.rm = na.rm)), info(f),
                     loose = mixed)
        }
        same_run(run(function() sum(ar, br, na.rm = na.rm)), run(function() sum(al, bl, na.rm = na.rm)),
                 info("sum(a, b)"), loose = mixed)
        same_run(run(function() sum(ar, 1, na.rm = na.rm)), run(function() sum(al, 1, na.rm = na.rm)),
                 info("sum(a, 1)"), loose = mixed)
        for (f in c("min", "max", "range")) {
            fun <- match.fun(f)
            same_run(run(function() fun(ar, na.rm = na.rm)), run(function() fun(al, na.rm = na.rm)), info(f))
            same_run(run(function() fun(ar, br, na.rm = na.rm)), run(function() fun(al, bl, na.rm = na.rm)),
                     info(paste0(f, "(a, b)")))
            same_run(run(function() fun(ar, -10, NA, na.rm = na.rm)),
                     run(function() fun(al, -10, NA, na.rm = na.rm)), info(paste0(f, "(a, -10, NA)")))
        }
        same_run(run(function() range(ar, na.rm = na.rm, finite = TRUE)),
                 run(function() range(al, na.rm = na.rm, finite = TRUE)), info("range(finite = TRUE)"))
        for (f in c("any", "all")) {
            fun <- match.fun(f)
            same_run(run(function() fun(ar > 0, na.rm = na.rm)), run(function() fun(al > 0, na.rm = na.rm)), info(f))
            same_run(run(function() fun(ar > 0, br < 1, na.rm = na.rm)),
                     run(function() fun(al > 0, bl < 1, na.rm = na.rm)), info(paste0(f, "(a > 0, b < 1)")))
            same_run(run(function() fun(is.na(ar), na.rm = na.rm)), run(function() fun(is.na(al), na.rm = na.rm)),
                     info(paste0(f, "(is.na(a))")))
        }
    }
}

## types: remote summaries are doubles (logical columns too), any/all logical
expect_identical(sum(a > 0, na.rm = TRUE), as.double(sum(d$a > 0, na.rm = TRUE)))
expect_true(is.double(sum(is.na(a))))
expect_true(is.logical(any(a > 0)))
expect_identical(sum(is.na(a)), as.double(sum(is.na(d$a))))
expect_identical(sum(is.nan(a)), as.double(sum(is.nan(d$a))))
expect_identical(sum(is.finite(a)), as.double(sum(is.finite(d$a))))

## one query per summary, whatever the number of arguments
nq <- function(f) {
    msgs <- character()
    op <- options(design.matrix.remote_debug = TRUE)
    on.exit(options(op))
    withCallingHandlers(f(), message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
    })
    sum(grepl("SELECT", msgs))
}
expect_identical(nq(function() sum(a, b, a * b, na.rm = TRUE)), 1L)
expect_identical(nq(function() range(a, b)), 1L)
expect_identical(nq(function() mean(a, na.rm = TRUE)), 1L)
expect_identical(nq(function() col_sums(a, b, n = is.na(a))), 1L)

## multi-column columns: every value of every column, as sum() of a matrix
M <- design_matrix(y ~ id, d)
B <- cbind(c(1, 0), c(0, 1), c(2, -1))
Pr <- R %*% B
Pm <- M %*% B
expect_identical(sum(Pr), sum(Pm))
expect_identical(range(Pr), range(Pm))
expect_identical(max(Pr, 1000), max(Pm, 1000))
expect_equal(mean(Pr), mean(Pm))
expect_identical(prod(Pr[, 1]), prod(Pm[, 1]))
expect_identical(all(Pr > -1000), all(Pm > -1000))
expect_identical(col_sums(P = Pr, a = 2), c(P = sum(Pm), a = 2))
## two separate products of one design, summarized together
expect_identical(sum(R %*% c(1, 0), R %*% c(0, 1)), sum(M %*% c(1, 0), M %*% c(0, 1)))

## ---- mean() with ordinary data: equal to rounding ----
set.seed(7)
n <- 400
dd <- data.frame(id = seq_len(n), y = rnorm(n), x = rnorm(n), g = sample(1:12, n, TRUE), w = runif(n, 0.5, 2))
DBI::dbWriteTable(con, "dd", dd)
Q <- design_matrix(y ~ x, con, table = "dd")
expect_equal(mean(design_column(Q, "y")), mean(dd$y))
expect_equal(sum(design_column(Q, "y")^2), sum(dd$y^2))
expect_equal(mean(exp(design_column(Q, "x"))), mean(exp(dd$x)))
expect_error(mean(design_column(Q, "y"), trim = 0.1), "trim")

## ---- Math: round, signif, trunc, gamma, lgamma ----
rv <- c(NA, NaN, Inf, -Inf, 0, -0.5, 0.5, 1.5, 2.5, -2.5, 0.15, 0.25, 0.35, -0.045, 1.005, 2.675, 123.456,
        -987.654321, 1234.5, 1e-300, 1e300, 5e-324, 1.23456789e-310, 0.1 + 0.2, 1 / 3, 2^52 + 0.5,
        12345678.9, -3, 7, 10.5, 170.5, 200)
dr <- data.frame(id = seq_along(rv), y = 1, v = rv)
DBI::dbWriteTable(con, "dr", dr)
Rr <- design_matrix(y ~ id, con, table = "dr")
v <- design_column(Rr, "v")
expect_identical(get(v), dr$v)
for (dig in c(0, 1, 2, 3, -1, -2, 0.4, 2.6, 15, 16, 309, 320, 400, -308, -400))
    expect_identical(get(round(v, dig)), round(dr$v, dig), info = sprintf("round(v, %s)", dig))
expect_identical(get(round(v)), round(dr$v))
expect_identical(get(round(v, digits = 2)), round(dr$v, digits = 2))
for (dig in c(1, 2, 3, 6, 10, 15, 17, 22, 23, 0, -1, 2.5, Inf, -Inf))
    expect_identical(get(signif(v, dig)), signif(dr$v, dig), info = sprintf("signif(v, %s)", dig))
expect_identical(get(signif(v)), signif(dr$v))
expect_identical(get(trunc(v)), trunc(dr$v))
## gamma and lgamma: the C library remotely, R's own code in memory
near <- function(r, l, info) {
    expect_identical(is.na(r), is.na(l), info = info)
    expect_identical(is.nan(r), is.nan(l), info = info)
    inf <- is.infinite(l)
    expect_identical(r[inf], l[inf], info = info)
    fin <- is.finite(l)
    expect_identical(is.finite(r), fin, info = info)
    expect_equal(r[fin], l[fin], tolerance = 1e-12, info = info)
}
near(get(lgamma(v)), suppressWarnings(lgamma(dr$v)), "lgamma")
near(get(gamma(v)), suppressWarnings(gamma(dr$v)), "gamma")
## Math functions on a logical column
expect_identical(get(round(v > 1, 1)), as.double(round(dr$v > 1, 1)))
## errors
expect_error(round(v, digits = v), "digits")
expect_error(round(v, digits = 1:2), "digits")
expect_error(cumsum(v), "no order")
expect_error(digamma(v), "not supported")

## ---- helpers: elementwise, in memory and remote ----

check <- function(remote, local, info) {
    expect_true(inherits(remote, "design_column"), info = info)
    expect_identical(get(remote), local, info = info)
}

## col_ifelse(): NA where the test is NA, double values, numbers recycle
check(col_ifelse(a > 0, a, b), col_ifelse(d$a > 0, d$a, d$b), "col_ifelse(a > 0, a, b)")
expect_identical(col_ifelse(d$a > 0, d$a, d$b), as.double(ifelse(d$a > 0, d$a, d$b)))
check(col_ifelse(is.na(a), 0, a), col_ifelse(is.na(d$a), 0, d$a), "col_ifelse(is.na(a), 0, a)")
check(col_ifelse(a > b, 1, NA), col_ifelse(d$a > d$b, 1, NA), "col_ifelse(a > b, 1, NA)")
check(col_ifelse(a >= 0, NaN, -1), col_ifelse(d$a >= 0, NaN, -1), "col_ifelse(a >= 0, NaN, -1)")
expect_identical(col_ifelse(d$a > 0, 1L, 0L), as.double(ifelse(d$a > 0, 1L, 0L)))
expect_identical(col_ifelse(c(TRUE, NA, FALSE), c(1, 2, 3), 0), c(1, NA, 0))
expect_identical(col_ifelse(logical(), 1, 0), numeric())

## col_pmin() and col_pmax(): pmin()/pmax() with R's missing values
for (f in c("col_pmin", "col_pmax")) {
    fun <- match.fun(f)
    base <- match.fun(sub("col_", "", f))
    check(fun(a, b), base(d$a, d$b), paste0(f, "(a, b)"))
    check(fun(b, a), base(d$b, d$a), paste0(f, "(b, a)"))
    for (cc in list(0, 1, -2.5, Inf, -Inf, NA, NaN)) {
        check(fun(a, cc), base(d$a, cc), sprintf("%s(a, %s)", f, cc))
        check(fun(cc, a), base(cc, d$a), sprintf("%s(%s, a)", f, cc))
    }
    expect_identical(fun(d$a, d$b), base(d$a, d$b))
    expect_identical(fun(2, 3), base(2, 3))
}
## nested clamps, and on a multi-column column
check(col_pmax(col_pmin(a, 2), -1), pmax(pmin(d$a, 2), -1), "clamp")
P2 <- as.data.frame(col_pmax(Pr[, 1:2], 3), keep = "id", nm = c("p1", "p2"))
P2 <- P2[order(P2$id), ]
expect_identical(unname(cbind(P2$p1, P2$p2)), pmax(Pm[, 1:2], 3))

## col_sums(): one value per argument, named after the arguments
## (no argument mixes NA and NaN, whose sum in R depends on the platform)
sums <- function(a, b) col_sums(n = is.finite(a), s = col_ifelse(is.finite(a), a, 0),
                                 a = col_ifelse(is.nan(a), 0, a), b = col_ifelse(is.na(b), -1, b), 2)
expect_identical(sums(a, b), sums(d$a, d$b))
expect_identical(names(sums(a, b)), c("n", "s", "a", "b", ""))
bare <- function(a, b) col_sums(a, b, 2)
expect_identical(names(bare(a, b)), c("a", "b", ""))
expect_identical(col_sums(a, b, na.rm = TRUE), col_sums(a = d$a, b = d$b, na.rm = TRUE))
expect_identical(col_sums(is.na(a), is.nan(a)), c(sum(is.na(d$a)), sum(is.nan(d$a))) + 0)
expect_identical(col_sums(), numeric())
expect_null(names(col_sums(a + 0, 1)))
## (a mix of NA and NaN: NA remotely; R's sum depends on the platform)
expect_true(is.na(col_sums(a)))
expect_true(is.na(col_sums(d$a)))
expect_identical(col_sums(x = col_ifelse(is.na(a), 0, a)), c(x = sum(ifelse(is.na(d$a), 0, d$a))))

## ---- errors: columns of two designs, vectors with remote columns, designs ----
R2 <- design_matrix(y ~ id, con, table = "d")
a2 <- design_column(R2, "a")
expect_error(col_pmin(a, a2), "different remote designs")
expect_error(col_pmax(a, a2), "different remote designs")
expect_error(col_ifelse(a > 0, a2, 0), "different remote designs")
expect_error(col_sums(a, a2), "different remote designs")
expect_error(sum(a, a2), "different remote designs")
expect_error(col_pmin(a, d$b), "combine only with numbers")
expect_error(col_ifelse(a > 0, d$b, 0), "combine only with numbers")
expect_error(sum(a, "x"), "combine only with numbers")
expect_error(col_sums(R), "is a design")
expect_error(col_pmin(R, 1), "is a design")
expect_error(col_pmin("x", 1), "must be a column or a number")
expect_error(col_ifelse("x", 1, 0), "logical")
expect_error(sum(a, foo = 1), "unknown argument")

## ---- rows scaled by a remote column: e * D, D / w ----
Mq <- design_matrix(y ~ x, dd)
bq <- drop(solve(crossprod(Mq), crossprod(Mq, dd$y)))
scores <- function(D, g) {
    e <- design_column(D, "y") - drop(D %*% bq)
    list(S = unname(colSums(e * D, group = g)), S2 = unname(colSums(D * e, group = g)),
         xwx = crossprod(D, design_column(D, "y")^2 * D),
         xdx = crossprod(D, D / (1 + design_column(D, "y")^2)))
}
sm <- scores(Mq, dd$g)
sr <- scores(Q, "g")
expect_equal(sr, sm)
expect_error(design_column(Q, "y") * R, "another remote design|different remote designs")

## ---- the estimators of the website page (website/primitives.typ) ----
ols_cluster <- function(D, cl) {
    y <- design_column(D, "y")
    XtX <- crossprod(D)
    b <- drop(solve(XtX, crossprod(D, y)))
    e <- y - drop(D %*% b)
    S <- colSums(e * D, group = cl)
    G <- nrow(S)
    bread <- solve(XtX)
    V <- G / (G - 1) * bread %*% crossprod(S) %*% bread
    fit <- cbind(estimate = b, std.error = sqrt(diag(V)))
    attr(fit, "rss") <- sum(e^2)
    fit
}
huber_step <- function(D, b, k = 1.345) {
    y <- design_column(D, "y")
    e <- col_compute(y - drop(D %*% b))
    s <- sqrt(sum(e^2) / (nrow(D) - ncol(D)))
    u <- abs(e) / s
    w <- col_pmin(1, k / u)
    b1 <- drop(solve(crossprod(D, w * D), crossprod(D, w * y)))
    list(coef = b1, scale = s, range = range(e),
         sums = col_sums(downweighted = u > k, weight = w, big = col_ifelse(abs(e) > 10, 1, 0)))
}
set.seed(48)
n <- 2000
dw <- data.frame(id = seq_len(n), x = rnorm(n), f = sample(c("a", "b", "c"), n, TRUE), cl = sample(1:40, n, TRUE))
dw$y <- 1 + dw$x + (dw$f == "b") + rnorm(40)[dw$cl] + rnorm(n)
out <- sample(n, 20)
dw$y[out] <- dw$y[out] + 20
DBI::dbWriteTable(con, "dw", dw)
Dm <- design_matrix(y ~ x + f, dw)
Dr <- design_matrix(y ~ x + f, con, table = "dw")
fm <- ols_cluster(Dm, dw$cl)
fr <- ols_cluster(Dr, "cl")
expect_equal(fr, fm)
## against the dense computation
X <- model.matrix(y ~ x + f, dw)
bx <- solve(crossprod(X), crossprod(X, dw$y))
ex <- drop(dw$y - X %*% bx)
Sx <- rowsum(ex * X, dw$cl)
Bx <- solve(crossprod(X))
expect_equal(unname(fr[, "std.error"]), unname(sqrt(diag(40 / 39 * Bx %*% crossprod(Sx) %*% Bx))))
expect_equal(unname(fr[, "estimate"]), unname(drop(bx)))
hm <- huber_step(Dm, fm[, "estimate"])
hr <- huber_step(Dr, fm[, "estimate"])
expect_equal(hr, hm)
expect_identical(hr$sums[["big"]], hm$sums[["big"]])
within_slope <- function(D, g) {
    n_g <- group_sums(1, g)[, 1]
    demean <- function(v) v - group_lookup(group_sums(v, g)[, 1] / n_g, g)
    yt <- demean(design_column(D, "y"))
    xt <- demean(design_column(D, "x"))
    col <- col_sums(xy = xt * yt, xx = xt^2)
    col[["xy"]] / col[["xx"]]
}
ws_m <- within_slope(Dm, design_groups(Dm, dw["cl"]))
ws_r <- within_slope(Dr, design_groups(Dr, "cl"))
expect_equal(ws_r, ws_m)
expect_equal(ws_m, coef(lm(y ~ x + factor(cl), dw))[["x"]])

test_con_done(con)

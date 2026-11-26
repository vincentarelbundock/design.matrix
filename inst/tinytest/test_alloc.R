## The fused in-memory kernels (src/fused.c, Kata design.matrix#bg3b) give
## exactly what the R expressions they replace give.
source("helpers.R")
library(design.matrix)
ns <- asNamespace("design.matrix")
cl <- function(nm) get(paste0("C_", nm), ns)

set.seed(11)
n <- 500L
z <- rnorm(n)
z[c(3, 9)] <- c(NA, NaN)
w <- runif(n)
w[c(5, 7)] <- c(0, NA)

## dmc_ifelse() (a missing test gives NA) and dmc_notnan() (a missing value
## gives FALSE, never NA)
test <- z > 0
r_ifelse <- function(test, yes, no) {
    out <- rep_len(as.double(no), length(test))
    i <- which(test)
    out[i] <- if (length(yes) == 1L) yes else yes[i]
    out[is.na(test)] <- NA
    out
}
expect_identical(ns$dmc_ifelse(test, w, 0), r_ifelse(test, w, 0))
expect_identical(ns$dmc_ifelse(test, 2L, w), r_ifelse(test, 2L, w))
expect_identical(ns$dmc_notnan(z), !is.na(z))
expect_identical(ns$dmc_notnan(c(1, NA, NaN, Inf, -Inf)), c(TRUE, FALSE, FALSE, TRUE, TRUE))
expect_identical(ns$dmc_notnan(c(1L, NA)), c(TRUE, FALSE))
expect_identical(ns$dmc_notnan(c(TRUE, NA)), c(TRUE, FALSE))
expect_identical(ns$dmc_notnan(matrix(c(1, NaN), 1L)), matrix(c(TRUE, FALSE), 1L))

## residuals and moments of dm_fe_core()
G <- c(7L, 4L)
codes <- list(sample.int(7L, n, TRUE), sample.int(4L, n, TRUE))
fe <- list(rnorm(7), rnorm(4))
zz <- rnorm(n)
xb <- rnorm(n)
ww <- runif(n)
for (wt in list(NULL, ww)) {
    r <- .Call(cl("dm_fe_resid"), zz, xb, wt, fe, codes, G)
    f <- xb + fe[[1L]][codes[[1L]]] + fe[[2L]][codes[[2L]]]
    expect_identical(f, .Call(cl("dm_fe_fitted_add"), xb, fe, codes))
    e <- zz - f
    we <- if (is.null(wt)) e else wt * e
    expect_identical(drop(r[[1L]]), we)
    expect_identical(r[[3L]], sum(we * e))
    for (k in 1:2) expect_identical(r[[2L]][[k]][, 1L],
                                    unname(drop(rowsum(we, codes[[k]], reorder = TRUE))[as.character(seq_len(G[k]))]))
    m <- .Call(cl("dm_fe_moments"), zz, wt)
    wz <- if (is.null(wt)) zz else wt * zz
    sw <- if (is.null(wt)) n else sum(wt)
    expect_identical(m[c(1L, 3L)], c(sum(wz * zz), sum(wz)))
    expect_identical(m[4L], if (is.null(wt)) sum((zz - sum(wz) / sw)^2) else sum(wt * (zz - sum(wz) / sw)^2))
}

## residual sums of dmfit_ls()
for (wt in list(NULL, w)) for (icpt in c(TRUE, FALSE)) {
    s <- .Call(cl("dm_ls_sums"), zz, xb, wt, icpt)
    r <- zz - xb
    rss <- if (is.null(wt)) sum(r^2) else sum(ifelse(wt > 0, wt * r^2, 0))
    ww1 <- if (is.null(wt)) 1 else wt
    fbar <- if (is.null(wt)) mean(xb) else sum(wt * xb) / sum(wt)
    mss <- if (icpt) sum(ww1 * (xb - fbar)^2) else sum(ww1 * xb^2)
    expect_identical(s, c(rss, mss))
}

## IRLS working columns, Poisson deviance, binomial log densities
mu <- runif(n, 0.05, 3)
eta <- log(mu)
y <- rpois(n, mu) + 0
W <- runif(n)
W[1:3] <- c(0, NA, -1)
v <- mu
me <- pmax(exp(eta), .Machine$double.eps)
k <- .Call(cl("dm_glm_working"), W, eta, mu, y, 0, v, me)
good <- (W > 0) & (me != 0)
expect_identical(k[[1L]], r_ifelse(good, W * me * me / v, 0))
expect_identical(k[[2L]], r_ifelse(good, (eta - 0) + (y - mu) / me, 0))
r_counts <- function(k, W, v, me, good) {
    wp <- !is.na(W) & W > 0
    as.double(c(sum(wp & is.na(v)), sum(wp & !is.na(v) & v == 0), sum(wp & is.na(me)), sum(good, na.rm = TRUE),
      sum(!(is.finite(k[[1L]]) & is.finite(k[[2L]])))))
}
expect_identical(k[[3L]], r_counts(k, W, v, me, good))
## missing variances and d mu / d eta are counted where W > 0 (as glm.fit's
## "NAs in V(mu)" and "NAs in d(mu)/d(eta)" checks)
v2 <- v
v2[c(1, 2, 5, 6)] <- c(NaN, NA, NaN, NA)
me2 <- me
me2[c(3, 7, 8)] <- c(NaN, NaN, NA)
k <- .Call(cl("dm_glm_working"), W, eta, mu, y, 0, v2, me2)
good <- (W > 0) & (me2 != 0)
expect_identical(k[[3L]], r_counts(k, W, v2, me2, good))
expect_identical(k[[3L]][c(1L, 3L)], c(2, 2))
expect_identical(.Call(cl("dm_poisson_dev"), y, mu, W),
                 r_ifelse(y > 0, 2 * W * (y * log(y / mu) - (y - mu)), 2 * mu * W))
p <- runif(n, 0.01, 0.99)
tr <- sample(1:5, n, TRUE) + 0
yb <- rbinom(n, tr, p) / tr
kk <- round(tr * yb)
mm <- round(tr)
expect_identical(.Call(cl("dm_binomial_ll"), tr, yb, W, p),
                 r_ifelse(tr > 0, W / tr, 0) * (lgamma(mm + 1) - lgamma(kk + 1) - lgamma(mm - kk + 1) +
                                                 r_ifelse(kk > 0, kk * log(p), 0) +
                                                 r_ifelse(mm - kk > 0, (mm - kk) * log(1 - p), 0)))

## The fused primitives (primitives_fused.R): the kernels on plain doubles,
## the generic code on other vectors (here named ones) and on remote
## columns, with the same values
nmd <- function(x) stats::setNames(x, seq_along(x))
for (wt in list(NULL, w)) for (icpt in c(TRUE, FALSE)) {
    s <- ns$dmc_ls_sums(zz, xb, wt, icpt)
    expect_identical(s, .Call(cl("dm_ls_sums"), zz, xb, wt, icpt))
    expect_identical(s, ns$dmc_ls_sums(nmd(zz), xb, wt, icpt))
}
expect_error(ns$dmc_ls_sums(zz, xb[-1L], NULL, TRUE), "one value per row")
k <- ns$dmc_irls_working(W, eta, mu, y, 0, v2, me2)
expect_identical(k, .Call(cl("dm_glm_working"), W, eta, mu, y, 0, v2, me2))
kg <- ns$dmc_irls_working(nmd(W), eta, mu, y, 0, v2, me2)
expect_identical(unname(kg[[1L]]), k[[1L]])
expect_identical(unname(kg[[2L]]), k[[2L]])
expect_identical(kg[[3L]], k[[3L]])
expect_error(ns$dmc_irls_working(W, eta, mu[-1L], y, 0, v, me), "lengths")
expect_identical(ns$dmc_poisson_dev(y, mu, W), .Call(cl("dm_poisson_dev"), y, mu, W))
expect_identical(unname(ns$dmc_poisson_dev(nmd(y), mu, W)), ns$dmc_poisson_dev(y, mu, W))
expect_identical(ns$dmc_binomial_ll(tr, yb, W, p), .Call(cl("dm_binomial_ll"), tr, yb, W, p))
expect_identical(unname(ns$dmc_binomial_ll(nmd(tr), yb, W, p)), ns$dmc_binomial_ll(tr, yb, W, p))
expect_identical(ns$dmc_freeze_all(y, mu), list(y, mu))
if (requireNamespace("duckdb", quietly = TRUE)) {
    ## (no missing values: SQL NULL is NA, not NaN, and the kernels count NaN)
    Wp <- abs(W)
    Wp[!is.finite(Wp)] <- 0
    wp <- w
    wp[is.na(wp)] <- 0
    df <- data.frame(zz = zz, xb = xb, w = wp, pw = Wp, eta = eta, mu = mu, y = y, v = v, me = me,
                     tr = tr, yb = yb, p = p)
    con <- DBI::dbConnect(duckdb::duckdb())
    DBI::dbWriteTable(con, "fz", df)
    R <- design_matrix(~ 1, con, table = "fz")
    rc <- function(nm) ns$dmp_column(R, sprintf('CAST(s."%s" AS DOUBLE)', nm))
    for (wt in list(NULL, "w")) for (icpt in c(TRUE, FALSE)) {
        expect_equal(ns$dmc_ls_sums(rc("zz"), rc("xb"), if (!is.null(wt)) rc(wt), icpt),
                     ns$dmc_ls_sums(zz, xb, if (!is.null(wt)) wp, icpt), tolerance = 1e-10)
    }
    km <- ns$dmc_irls_working(Wp, eta, mu, y, 0, v, me)
    for (tb in c(FALSE, TRUE)) {
        kr <- ns$dmc_irls_working(rc("pw"), rc("eta"), rc("mu"), rc("y"), 0, rc("v"), rc("me"), table = tb)
        expect_equal(ns$dmp_collect(kr[[1L]]), km[[1L]], tolerance = 1e-12)
        expect_equal(ns$dmp_collect(kr[[2L]]), km[[2L]], tolerance = 1e-12)
        expect_identical(kr[[3L]], km[[3L]])
    }
    expect_equal(ns$dmp_collect(ns$dmc_poisson_dev(rc("y"), rc("mu"), rc("pw"))),
                 ns$dmc_poisson_dev(y, mu, Wp), tolerance = 1e-12)
    expect_equal(ns$dmp_collect(ns$dmc_binomial_ll(rc("tr"), rc("yb"), rc("pw"), rc("p"))),
                 ns$dmc_binomial_ll(tr, yb, Wp, p), tolerance = 1e-12)
    fa <- ns$dmc_freeze_all(rc("y") + 1, rc("mu"))
    expect_equal(ns$dmp_collect(fa[[1L]] * fa[[2L]]), (y + 1) * mu, tolerance = 1e-14)
    DBI::dbDisconnect(con, shutdown = TRUE)
}

## unused factor levels, and ALTREP variables copied once
f <- factor(c("a", "b", NA), levels = c("a", "b", "c"))
expect_true(.Call(cl("dm_levels_unused"), f, 3L))
expect_false(.Call(cl("dm_levels_unused"), droplevels(f), 2L))
vars <- list(a = 1:10, b = structure(c(1L, 2L), levels = c("x", "y"), class = "factor"))
u <- .Call(cl("dm_unwrap_vars"), vars)
expect_identical(u, vars)
d <- data.frame(y = rnorm(20), g = factor(rep(c("p", "q"), 10)))
expect_equal(coef(design_lm(y ~ g, d)), coef(lm(y ~ g, d)), tolerance = 1e-12)

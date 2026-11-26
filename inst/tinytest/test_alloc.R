## The fused in-memory kernels (src/fused.c, Kata design.matrix#bg3b) give
## exactly what the R expressions they replace give.
source("helpers.R")
library(xtx)
ns <- asNamespace("xtx")
cl <- function(nm) get(paste0("C_", nm), ns)

set.seed(11)
n <- 500L
z <- rnorm(n)
z[c(3, 9)] <- c(NA, NaN)
w <- runif(n)
w[c(5, 7)] <- c(0, NA)

## xtx_c_ifelse() (a missing test gives NA) and xtx_c_present() (a missing value
## gives FALSE, never NA)
test <- z > 0
r_ifelse <- function(test, yes, no) {
    out <- rep_len(as.double(no), length(test))
    i <- which(test)
    out[i] <- if (length(yes) == 1L) yes else yes[i]
    out[is.na(test)] <- NA
    out
}
expect_identical(ns$xtx_c_ifelse(test, w, 0), r_ifelse(test, w, 0))
expect_identical(ns$xtx_c_ifelse(test, 2L, w), r_ifelse(test, 2L, w))
expect_identical(ns$xtx_c_present(z), !is.na(z))
expect_identical(ns$xtx_c_present(c(1, NA, NaN, Inf, -Inf)), c(TRUE, FALSE, FALSE, TRUE, TRUE))
expect_identical(ns$xtx_c_present(c(1L, NA)), c(TRUE, FALSE))
expect_identical(ns$xtx_c_present(c(TRUE, NA)), c(TRUE, FALSE))
expect_identical(ns$xtx_c_present(matrix(c(1, NaN), 1L)), matrix(c(TRUE, FALSE), 1L))

## residuals and moments of xtx_fe_core()
G <- c(7L, 4L)
codes <- list(sample.int(7L, n, TRUE), sample.int(4L, n, TRUE))
fe <- list(rnorm(7), rnorm(4))
zz <- rnorm(n)
xb <- rnorm(n)
ww <- runif(n)
for (wt in list(NULL, ww)) {
    r <- .Call(cl("xtx_fe_resid"), zz, xb, wt, fe, codes, G)
    f <- xb + fe[[1L]][codes[[1L]]] + fe[[2L]][codes[[2L]]]
    expect_identical(f, .Call(cl("xtx_fe_fitted_add"), xb, fe, codes))
    e <- zz - f
    we <- if (is.null(wt)) e else wt * e
    expect_identical(drop(r[[1L]]), we)
    expect_identical(r[[3L]], sum(we * e))
    for (k in 1:2) expect_identical(r[[2L]][[k]][, 1L],
                                    unname(drop(rowsum(we, codes[[k]], reorder = TRUE))[as.character(seq_len(G[k]))]))
    m <- .Call(cl("xtx_fe_moments"), zz, wt)
    wz <- if (is.null(wt)) zz else wt * zz
    sw <- if (is.null(wt)) n else sum(wt)
    expect_identical(m[c(1L, 3L)], c(sum(wz * zz), sum(wz)))
    expect_identical(m[4L], if (is.null(wt)) sum((zz - sum(wz) / sw)^2) else sum(wt * (zz - sum(wz) / sw)^2))
}

## residual sums of xtx_fit_ls()
for (wt in list(NULL, w)) for (icpt in c(TRUE, FALSE)) {
    s <- .Call(cl("xtx_ls_sums"), zz, xb, wt, icpt)
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
k <- .Call(cl("xtx_glm_working"), W, eta, mu, y, 0, v, me)
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
k <- .Call(cl("xtx_glm_working"), W, eta, mu, y, 0, v2, me2)
good <- (W > 0) & (me2 != 0)
expect_identical(k[[3L]], r_counts(k, W, v2, me2, good))
expect_identical(k[[3L]][c(1L, 3L)], c(2, 2))
expect_identical(.Call(cl("xtx_poisson_dev"), y, mu, W),
                 r_ifelse(y > 0, 2 * W * (y * log(y / mu) - (y - mu)), 2 * mu * W))
p <- runif(n, 0.01, 0.99)
tr <- sample(1:5, n, TRUE) + 0
yb <- rbinom(n, tr, p) / tr
kk <- round(tr * yb)
mm <- round(tr)
expect_identical(.Call(cl("xtx_binomial_ll"), tr, yb, W, p),
                 r_ifelse(tr > 0, W / tr, 0) * (lgamma(mm + 1) - lgamma(kk + 1) - lgamma(mm - kk + 1) +
                                                 r_ifelse(kk > 0, kk * log(p), 0) +
                                                 r_ifelse(mm - kk > 0, (mm - kk) * log(1 - p), 0)))

## The fused primitives (primitives_fused.R): the kernels on plain doubles,
## the generic code on other vectors (here named ones) and on remote
## columns, with the same values
nmd <- function(x) stats::setNames(x, seq_along(x))
for (wt in list(NULL, w)) for (icpt in c(TRUE, FALSE)) {
    s <- ns$xtx_c_ls_sums(zz, xb, wt, icpt)
    expect_identical(s, .Call(cl("xtx_ls_sums"), zz, xb, wt, icpt))
    expect_identical(s, ns$xtx_c_ls_sums(nmd(zz), xb, wt, icpt))
}
expect_error(ns$xtx_c_ls_sums(zz, xb[-1L], NULL, TRUE), "one value per row")
k <- ns$xtx_c_irls_working(W, eta, mu, y, 0, v2, me2)
expect_identical(k, .Call(cl("xtx_glm_working"), W, eta, mu, y, 0, v2, me2))
kg <- ns$xtx_c_irls_working(nmd(W), eta, mu, y, 0, v2, me2)
expect_identical(unname(kg[[1L]]), k[[1L]])
expect_identical(unname(kg[[2L]]), k[[2L]])
expect_identical(kg[[3L]], k[[3L]])
expect_error(ns$xtx_c_irls_working(W, eta, mu[-1L], y, 0, v, me), "lengths")
expect_identical(ns$xtx_c_poisson_dev(y, mu, W), .Call(cl("xtx_poisson_dev"), y, mu, W))
expect_identical(unname(ns$xtx_c_poisson_dev(nmd(y), mu, W)), ns$xtx_c_poisson_dev(y, mu, W))
expect_identical(ns$xtx_c_binomial_ll(tr, yb, W, p), .Call(cl("xtx_binomial_ll"), tr, yb, W, p))
expect_identical(unname(ns$xtx_c_binomial_ll(nmd(tr), yb, W, p)), ns$xtx_c_binomial_ll(tr, yb, W, p))
expect_identical(ns$xtx_c_freeze_all(y, mu), list(y, mu))
## the fused primitives of the FE fits: kernel on plain doubles, generic
## code on named vectors (and on codes with attributes), same values
gcodes <- Map(ns$xtx_c_codes, codes, G)
gnamed <- lapply(gcodes, function(g) { g$codes <- nmd(g$codes); g })
for (wt in list(NULL, ww)) {
    m <- ns$xtx_c_moments(zz, wt, n)
    expect_identical(m, .Call(cl("xtx_fe_moments"), zz, wt))
    expect_identical(m[2L], if (is.null(wt)) 0 else sum(wt))
    expect_equal(ns$xtx_c_moments(nmd(zz), wt, n), m, tolerance = 1e-15)
    r <- ns$xtx_c_fe_resid(zz, xb, wt, fe, gcodes)
    expect_identical(r, .Call(cl("xtx_fe_resid"), zz, xb, wt, fe, codes, G))
    rg <- ns$xtx_c_fe_resid(nmd(zz), xb, wt, fe, gcodes)
    expect_identical(unname(drop(rg[[1L]])), drop(r[[1L]]))
    expect_equal(rg[[2L]], r[[2L]], tolerance = 1e-14)
    expect_equal(rg[[3L]], r[[3L]], tolerance = 1e-14)
    expect_identical(ns$xtx_c_fe_resid(zz, xb, wt, fe, gnamed)[[1L]], drop(r[[1L]]))
}
fa <- ns$xtx_c_fe_add(xb, fe, gcodes)
expect_identical(fa, .Call(cl("xtx_fe_fitted_add"), xb, fe, codes))
expect_identical(unname(ns$xtx_c_fe_add(nmd(xb), fe, gcodes)), fa)
expect_identical(ns$xtx_c_fe_add(xb, fe, gnamed), fa)
expect_identical(ns$xtx_c_fe_add(xb, list(), list()), xb)
## (k columns: the generic code)
expect_identical(ns$xtx_c_fe_add(cbind(xb, xb), lapply(fe, function(v) cbind(v, 2 * v)), gcodes),
                 cbind(xb, xb) + cbind(fe[[1L]][codes[[1L]]], 2 * fe[[1L]][codes[[1L]]]) +
                 cbind(fe[[2L]][codes[[2L]]], 2 * fe[[2L]][codes[[2L]]]))
for (k in 1:2) {
    expect_identical(ns$xtx_p_groupcount(gcodes[[k]]), as.double(tabulate(codes[[k]], G[k])))
    expect_identical(ns$xtx_p_groupcount(gnamed[[k]]), ns$xtx_p_groupcount(gcodes[[k]]))
    expect_identical(ns$xtx_p_groupcount(gcodes[[k]], 1:3), c(1, 2, 3))
}
fl <- cbind(seq_len(7L) %% 2L == 0L, seq_len(7L) > 5L)
lf <- ns$xtx_c_lookup_flags(fl, gcodes[[1L]])
expect_identical(lf, list(fl[codes[[1L]], 1L], fl[codes[[1L]], 2L]))
expect_identical(ns$xtx_c_lookup_flags(fl, gnamed[[1L]]), lf)
## codes are checked before a kernel indexes by them (out of 1..G, NA)
for (bad in list(c(codes[[1L]][-1L], 8L), c(codes[[1L]][-1L], 0L), c(codes[[1L]][-1L], NA))) {
    expect_error(.Call(cl("xtx_fe_resid"), zz, xb, NULL, fe, list(bad, codes[[2L]]), G), "outside 1..G")
    expect_error(.Call(cl("xtx_fe_gsum"), bad, 7L, zz, 1L), "outside 1..G")
    expect_error(.Call(cl("xtx_fe_crosstabs"), bad, list(codes[[2L]]), NULL, 7L, 4L, 1L), "outside 1..G")
    expect_error(.Call(cl("xtx_fe_crosstabs"), codes[[2L]], list(bad), NULL, 4L, 7L, 1L), "outside 1..G")
    expect_error(.Call(cl("xtx_fe_nested"), codes[[2L]], bad, 7L), "outside 1..G")
    expect_error(.Call(cl("xtx_fe_recode"), bad, 7:1), "outside 1..G")
    ## (lookups give NA)
    expect_true(is.na(ns$xtx_c_fe_add(xb, fe, list(ns$xtx_c_codes(bad, 7L), gcodes[[2L]]))[n]))
}
expect_error(.Call(cl("xtx_fe_gsum"), codes[[1L]] + 0, 7L, zz, 1L), "integer")
if (requireNamespace("duckdb", quietly = TRUE)) {
    ## (no missing values: SQL NULL is NA, not NaN, and the kernels count NaN)
    Wp <- abs(W)
    Wp[!is.finite(Wp)] <- 0
    wp <- w
    wp[is.na(wp)] <- 0
    df <- data.frame(zz = zz, xb = xb, w = wp, pw = Wp, eta = eta, mu = mu, y = y, v = v, me = me,
                     tr = tr, yb = yb, p = p)
    con <- test_con()
    DBI::dbWriteTable(con, "fz", df)
    R <- design_matrix(~ 1, con, table = "fz")
    rc <- function(nm) ns$xtx_p_column(R, sprintf('CAST(s."%s" AS DOUBLE)', nm))
    for (wt in list(NULL, "w")) for (icpt in c(TRUE, FALSE)) {
        expect_equal(ns$xtx_c_ls_sums(rc("zz"), rc("xb"), if (!is.null(wt)) rc(wt), icpt),
                     ns$xtx_c_ls_sums(zz, xb, if (!is.null(wt)) wp, icpt), tolerance = 1e-10)
    }
    km <- ns$xtx_c_irls_working(Wp, eta, mu, y, 0, v, me)
    for (tb in c(FALSE, TRUE)) {
        kr <- ns$xtx_c_irls_working(rc("pw"), rc("eta"), rc("mu"), rc("y"), 0, rc("v"), rc("me"), table = tb)
        expect_equal(ns$xtx_p_collect(kr[[1L]]), km[[1L]], tolerance = 1e-12)
        expect_equal(ns$xtx_p_collect(kr[[2L]]), km[[2L]], tolerance = 1e-12)
        expect_identical(kr[[3L]], km[[3L]])
    }
    expect_equal(ns$xtx_p_collect(ns$xtx_c_poisson_dev(rc("y"), rc("mu"), rc("pw"))),
                 ns$xtx_c_poisson_dev(y, mu, Wp), tolerance = 1e-12)
    expect_equal(ns$xtx_p_collect(ns$xtx_c_binomial_ll(rc("tr"), rc("yb"), rc("pw"), rc("p"))),
                 ns$xtx_c_binomial_ll(tr, yb, Wp, p), tolerance = 1e-12)
    ## the fused primitives of the FE fits on remote columns
    DBI::dbWriteTable(con, "fzc", data.frame(zz = zz, xb = xb, w = ww, c1 = codes[[1L]], c2 = codes[[2L]]))
    RC <- design_matrix(~ 1, con, table = "fzc")
    rcc <- function(nm) ns$xtx_p_column(RC, sprintf('CAST(s."%s" AS DOUBLE)', nm))
    rcodes <- list(ns$xtx_c_codes(ns$xtx_p_column(RC, 's."c1"'), G[1L]),
                   ns$xtx_c_codes(ns$xtx_p_column(RC, 's."c2"'), G[2L]))
    for (wt in list(NULL, "w")) {
        wm <- if (!is.null(wt)) ww
        wr <- if (!is.null(wt)) rcc(wt)
        expect_equal(ns$xtx_c_moments(rcc("zz"), wr, n), ns$xtx_c_moments(zz, wm, n), tolerance = 1e-12)
        rr <- ns$xtx_c_fe_resid(rcc("zz"), rcc("xb"), wr, fe, rcodes)
        rm <- ns$xtx_c_fe_resid(zz, xb, wm, fe, gcodes)
        expect_equal(ns$xtx_p_collect(rr[[1L]]), drop(rm[[1L]]), tolerance = 1e-12)
        expect_equal(rr[[2L]], rm[[2L]], tolerance = 1e-12)
        expect_equal(rr[[3L]], rm[[3L]], tolerance = 1e-12)
    }
    expect_equal(ns$xtx_p_collect(ns$xtx_c_fe_add(rcc("xb"), fe, rcodes)), fa, tolerance = 1e-14)
    for (k in 1:2) expect_identical(ns$xtx_p_groupcount(rcodes[[k]]), ns$xtx_p_groupcount(gcodes[[k]]))
    lr <- ns$xtx_c_lookup_flags(fl, rcodes[[1L]])
    expect_identical(lapply(lr, function(v) as.logical(ns$xtx_p_collect(v))), lf)
    fa <- ns$xtx_c_freeze_all(rc("y") + 1, rc("mu"))
    expect_equal(ns$xtx_p_collect(fa[[1L]] * fa[[2L]]), (y + 1) * mu, tolerance = 1e-14)
    test_con_done(con)
}

## unused factor levels, and ALTREP variables copied once
f <- factor(c("a", "b", NA), levels = c("a", "b", "c"))
expect_true(.Call(cl("xtx_levels_unused"), f, 3L))
expect_false(.Call(cl("xtx_levels_unused"), droplevels(f), 2L))
vars <- list(a = 1:10, b = structure(c(1L, 2L), levels = c("x", "y"), class = "factor"))
u <- .Call(cl("xtx_unwrap_vars"), vars)
expect_identical(u, vars)
d <- data.frame(y = rnorm(20), g = factor(rep(c("p", "q"), 10)))
expect_equal(coef(design_lm(y ~ g, d)), coef(lm(y ~ g, d)), tolerance = 1e-12)

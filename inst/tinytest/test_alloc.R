## The fused in-memory kernels (src/fused.c) give exactly what the R
## expressions they replace give. (Those of the fits: xtx.fit, test_kernels.R.)
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

## the moments of a response (xtx_c_moments())
G <- c(7L, 4L)
codes <- list(sample.int(7L, n, TRUE), sample.int(4L, n, TRUE))
zz <- rnorm(n)
ww <- runif(n)
nmd <- function(x) stats::setNames(x, seq_along(x))
for (wt in list(NULL, ww)) {
    m <- .Call(cl("xtx_fe_moments"), zz, wt)
    wz <- if (is.null(wt)) zz else wt * zz
    sw <- if (is.null(wt)) n else sum(wt)
    expect_identical(m[c(1L, 3L)], c(sum(wz * zz), sum(wz)))
    expect_identical(m[4L], if (is.null(wt)) sum((zz - sum(wz) / sw)^2) else sum(wt * (zz - sum(wz) / sw)^2))
    ## the kernel on plain doubles, the generic code on other vectors
    expect_identical(ns$xtx_c_moments(zz, wt, n), m)
    expect_identical(m[2L], if (is.null(wt)) 0 else sum(wt))
    expect_equal(ns$xtx_c_moments(nmd(zz), wt, n), m, tolerance = 1e-15)
}

## codes are checked before a kernel indexes by them (out of 1..G, NA)
for (bad in list(c(codes[[1L]][-1L], 8L), c(codes[[1L]][-1L], 0L), c(codes[[1L]][-1L], NA))) {
    expect_error(.Call(cl("xtx_fe_gsum"), bad, 7L, zz, 1L), "outside 1..G")
    expect_error(.Call(cl("xtx_fe_crosstabs"), bad, list(codes[[2L]]), NULL, 7L, 4L, 1L), "outside 1..G")
    expect_error(.Call(cl("xtx_fe_crosstabs"), codes[[2L]], list(bad), NULL, 4L, 7L, 1L), "outside 1..G")
    expect_error(.Call(cl("xtx_fe_nested"), codes[[2L]], bad, 7L), "outside 1..G")
    expect_error(.Call(cl("xtx_fe_recode"), bad, 7:1), "outside 1..G")
}
expect_error(.Call(cl("xtx_fe_gsum"), codes[[1L]] + 0, 7L, zz, 1L), "integer")
if (requireNamespace("duckdb", quietly = TRUE)) {
    con <- test_con()
    DBI::dbWriteTable(con, "fzc", data.frame(zz = zz, w = ww, c1 = codes[[1L]], c2 = codes[[2L]]))
    RC <- design_matrix(~ 1, con, table = "fzc")
    rcc <- function(nm) ns$xtx_p_column(RC, sprintf('CAST(s."%s" AS DOUBLE)', nm))
    for (wt in list(NULL, "w")) {
        wm <- if (!is.null(wt)) ww
        wr <- if (!is.null(wt)) rcc(wt)
        expect_equal(ns$xtx_c_moments(rcc("zz"), wr, n), ns$xtx_c_moments(zz, wm, n), tolerance = 1e-12)
    }
    test_con_done(con)
}

## unused factor levels, and ALTREP variables copied once
f <- factor(c("a", "b", NA), levels = c("a", "b", "c"))
expect_true(.Call(cl("xtx_levels_unused"), f, 3L))
expect_false(.Call(cl("xtx_levels_unused"), droplevels(f), 2L))
vars <- list(a = 1:10, b = structure(c(1L, 2L), levels = c("x", "y"), class = "factor"))
u <- .Call(cl("xtx_unwrap_vars"), vars)
expect_identical(u, vars)

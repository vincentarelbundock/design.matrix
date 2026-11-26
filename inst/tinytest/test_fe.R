# Prototype fixed-effect absorption (R/fe.R, internal dm_fe_fit()) against
# fixest::feols(): coefficients, iid and clustered covariance matrices,
# residuals and number of observations, with feols()'s default small sample
# corrections and singleton removal. feols() demeans iteratively; a tight
# fixef.tol makes its results exact enough for a 1e-8 comparison.

if (!requireNamespace("fixest", quietly = TRUE)) exit_file("fixest not installed")
library(design.matrix)
fe_fit <- design.matrix:::dm_fe_fit

TOL <- 1e-8

fx_fit <- function(fml, d, weights = NULL, ...) {
    if (is.null(weights)) {
        fixest::feols(fml, d, fixef.tol = 1e-11, notes = FALSE, warn = FALSE, ...)
    } else {
        fixest::feols(fml, d, weights = weights, fixef.tol = 1e-11, notes = FALSE, warn = FALSE, ...)
    }
}

plain <- function(V) {
    V <- unclass(as.matrix(V))
    attributes(V) <- list(dim = dim(V))
    V
}

# Compare a dm_fe_fit() with a feols() fit on the same model.
check_fe <- function(fml, d, cluster = NULL, weights = NULL, ..., info = "") {
    a <- fe_fit(fml, d, cluster = cluster, weights = weights, ...)
    b <- fx_fit(fml, d, weights = weights)
    ok <- !is.na(a$coefficients)
    expect_equal(names(a$coefficients)[ok], names(coef(b)), info = paste(info, "names"))
    expect_equal(unname(a$coefficients[ok]), unname(coef(b)), tolerance = TOL, info = paste(info, "coef"))
    expect_equal(a$nobs, b$nobs, info = paste(info, "nobs"))
    expect_equal(a$residuals, unname(resid(b)), tolerance = TOL, info = paste(info, "resid"))
    expect_equal(plain(a$vcov_iid[ok, ok]), plain(vcov(b, vcov = "iid")),
                 tolerance = TOL, info = paste(info, "vcov iid"))
    if (!is.null(cluster))
        expect_equal(plain(a$vcov_cluster[ok, ok]), plain(vcov(b, cluster = cluster)),
                     tolerance = TOL, info = paste(info, "vcov cluster"))
    invisible(a)
}

set.seed(20261125)
n <- 4000
d <- data.frame(
    a = sample(600, n, TRUE),
    b = sample(40, n, TRUE),
    c = sample(7, n, TRUE),
    x1 = rnorm(n), x2 = rnorm(n),
    f = factor(sample(letters[1:4], n, TRUE)),
    h = factor(sample(c("u", "v", "w"), n, TRUE)),
    w = runif(n, 0.2, 3))
d$y <- d$x1 - 0.5 * d$x2 + d$a / 300 + d$b / 20 + as.integer(d$f) + rnorm(n)

# One-way ------------------------------------------------------------------
check_fe(y ~ x1 + x2 | a, d, info = "one-way")
check_fe(y ~ x1 + x2 | a, d, cluster = ~a, info = "one-way, cluster nested")
check_fe(y ~ x1 + x2 | a, d, cluster = ~b, info = "one-way, cluster not nested")
check_fe(y ~ x1 + x2 | b, d, cluster = ~b, weights = ~w, info = "one-way weighted")

# Two-way and three-way ------------------------------------------------------
check_fe(y ~ x1 + x2 | a + b, d, cluster = ~a, info = "two-way")
check_fe(y ~ x1 + x2 | a + b, d, cluster = ~b, info = "two-way, cluster b")
check_fe(y ~ x1 + x2 | a + b, d, cluster = ~c, info = "two-way, cluster not nested")
check_fe(y ~ x1 + x2 | a + b + c, d, cluster = ~a, info = "three-way")
check_fe(y ~ x1 + x2 | a + b + c, d, cluster = ~c, weights = ~w, info = "three-way weighted")
check_fe(y ~ x1 + x2 | a^c + b, d, cluster = ~b, info = "interacted FE")

# Every solver of the reduced system gives the same fit
for (s in c("dense", "chol", "pcg")) {
    a <- check_fe(y ~ x1 + x2 | a + b + c, d, cluster = ~b, solver = s, info = s)
    expect_equal(a$solver$solver, s, info = s)
}

# Factor regressors and interactions --------------------------------------
check_fe(y ~ x1 * f + x2 | a + b, d, cluster = ~a, info = "factor interaction")
check_fe(y ~ f:h + x1 | a + c, d, cluster = ~c, weights = ~w, info = "factor x factor")
check_fe(y ~ poly(x1, 2) + log(w) | b, d, cluster = ~b, info = "transformations")

# Unbalanced: a few large levels, many small ones ----------------------------
du <- d
du$a <- ifelse(runif(n) < 0.5, sample(3, n, TRUE), du$a)
du$b <- ifelse(runif(n) < 0.3, 1L, du$b)
check_fe(y ~ x1 + x2 | a + b, du, cluster = ~a, info = "unbalanced")

# Disconnected components of the FE graph -----------------------------------
dd <- d
half <- dd$a > 300
dd$b <- ifelse(half, dd$b + 100L, dd$b) # a in 1:300 only meets b in 1:40, a > 300 b > 100
a <- check_fe(y ~ x1 + x2 | a + b, dd, cluster = ~a, info = "disconnected")
expect_equal(a$solver$ncomp, 2L)
a <- check_fe(y ~ x1 + x2 | a + b, dd, cluster = ~a, solver = "chol", info = "disconnected chol")
a <- check_fe(y ~ x1 + x2 | a + b + c, dd, cluster = ~b, info = "disconnected three-way")

# Collinear regressors --------------------------------------------------------
dc <- d
dc$x3 <- dc$x1 + 2 * dc$x2                    # collinear with x1, x2
dc$xa <- (dc$a %% 13) / 7                     # constant within a: absorbed
dc$xb <- dc$b / 3 - (dc$a %% 5)              # sum of an a and a b effect
a <- check_fe(y ~ x1 + x2 + x3 | a + b, dc, cluster = ~a, info = "collinear x")
expect_true(is.na(a$coefficients[["x3"]]))
a <- check_fe(y ~ x1 + xa + x2 | a + b, dc, cluster = ~a, info = "collinear with FE")
expect_true(is.na(a$coefficients[["xa"]]))
a <- check_fe(y ~ xb + x1 | a + b, dc, cluster = ~b, info = "collinear with two FE")
expect_true(is.na(a$coefficients[["xb"]]))

# Singletons ----------------------------------------------------------------
ds <- d
ds$a[1:50] <- 10000L + seq_len(50)            # 50 singleton levels
ds$b[51] <- 999L; ds$a[51] <- 20000L          # a chain: 51 is alone in b ...
ds$a[52] <- 20000L                            # ... then 52 is alone in a
a <- check_fe(y ~ x1 + x2 | a + b, ds, cluster = ~b, info = "singletons")
# the random data has singletons of its own; nobs matched feols above
expect_true(a$nsingletons >= 52L)
expect_equal(a$nsingletons, nrow(ds) - fx_fit(y ~ x1 + x2 | a + b, ds)$nobs)
# kept singletons: feols(fixef.rm = "none")
a <- fe_fit(y ~ x1 + x2 | a + b, ds, singletons = "keep", cluster = ~b)
b <- fixest::feols(y ~ x1 + x2 | a + b, ds, fixef.rm = "none", fixef.tol = 1e-11, notes = FALSE)
expect_equal(a$nobs, nrow(ds))
expect_equal(unname(a$coefficients), unname(coef(b)), tolerance = TOL)
expect_equal(plain(a$vcov_iid), plain(vcov(b, vcov = "iid")), tolerance = TOL)
expect_equal(plain(a$vcov_cluster), plain(vcov(b, cluster = ~b)), tolerance = TOL)

# Missing values and zero weights ----------------------------------------
dm <- d
dm$x1[c(3, 9)] <- NA
dm$b[17] <- NA
dm$w[c(5, 6)] <- 0
check_fe(y ~ x1 + x2 | a + b, dm, cluster = ~a, weights = ~w, info = "NA and zero weights")

# Precision: regressors with a large mean relative to their spread ---------
# X'MX = X'WX - S'K^-S cancels like a one-pass variance: the relative error
# grows like eps * (mean / sd)^2. One step of iterative refinement from the
# row-wise residuals restores the coefficients (not the covariance).
dp <- d
dp$x1 <- 1e5 + dp$x1
dp$y <- dp$y + dp$x1
b <- fx_fit(y ~ x1 + x2 | a + b, dp)
a0 <- fe_fit(y ~ x1 + x2 | a + b, dp)
a1 <- fe_fit(y ~ x1 + x2 | a + b, dp, refine = 1L)
expect_true(max(abs(a0$coefficients / coef(b) - 1)) > 1e-9)
expect_equal(unname(a1$coefficients), unname(coef(b)), tolerance = TOL)
expect_equal(a1$residuals, unname(resid(b)), tolerance = TOL)
expect_equal(as.vector(sqrt(diag(a1$vcov_iid))), as.vector(fixest::se(b, vcov = "iid")), tolerance = 1e-5)

# Exact FE degrees of freedom (fixest's K.exact = TRUE) --------------------
a <- fe_fit(y ~ x1 + x2 | a + b, dd)
expect_equal(a$fe_dof_exact, sum(a$fe_sizes) - 2L) # two components
b <- fx_fit(y ~ x1 + x2 | a + b, dd)
expect_equal(plain(a$vcov_iid * (a$df.residual) / (a$nobs - a$rank - a$fe_dof_exact)),
             plain(vcov(b, vcov = "iid", ssc = fixest::ssc(K.exact = TRUE))), tolerance = TOL)
a <- fe_fit(y ~ x1 + x2 | a + b + c, d)
expect_equal(a$fe_dof_exact, a$fe_dof)

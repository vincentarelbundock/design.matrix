source("helpers.R")
library(design.matrix)

# In memory, qr() of a design with numeric variables far from zero factors
# the X'X of the centered design (as remote designs do), so that R, hence
# vcov(), keeps the accuracy of lm() on badly conditioned designs. The
# columns kept and the NA pattern still come from the uncentered X'X, as in
# lm(). The in-memory versions of the checks of test_remote_lm.R.
same <- function(a, b, tol = 1e-7) isTRUE(all.equal(unname(a), unname(b), tolerance = tol))

set.seed(42)
n <- 3000
d <- data.frame(
    xd = sample(1:5, n, TRUE),
    xc = rnorm(n),
    xe = rnorm(n),
    f1 = sample(c("b", "a", "c", "d"), n, TRUE),
    f2 = factor(sample(c("u", "v", "w"), n, TRUE), levels = c("w", "u", "v", "unused")),
    lg = sample(c(TRUE, FALSE), n, TRUE),
    w = runif(n, 0.5, 2))
d$y <- d$xd + d$xc + as.integer(factor(d$f1)) + rnorm(n)
d$xs <- 100 * d$xc + 1e7  # far from zero, not aliased by lm()
d$xt <- d$xc + 1e7        # so close to aliased that lm() drops it
d$xd7 <- d$xd + 1e7
d$xb <- 1e7 + 10 * d$xe   # ~1e7 with a small variance
d$f3 <- d$f1              # aliased with f1
d$y[c(3, 50)] <- NA
d$f1[7] <- NA
d$xc[9] <- NaN
d$w[11] <- 0

check_fit <- function(f, l, info, tol = 1e-7) {
    expect_identical(is.na(coef(f)), is.na(coef(l)), info = paste(info, "aliased"))
    expect_true(same(coef(f), coef(l), tol), info = paste(info, "coef"))
    expect_true(same(vcov(f), vcov(l), tol), info = paste(info, "vcov"))
    expect_true(same(sqrt(diag(vcov(f))), sqrt(diag(vcov(l))), tol), info = paste(info, "se"))
    expect_true(same(sigma(f), sigma(l), tol), info = paste(info, "sigma"))
    expect_identical(df.residual(f), df.residual(l), info = paste(info, "df"))
    expect_true(same(summary(f)$r.squared, summary(l)$r.squared, tol), info = paste(info, "R2"))
}

formulas <- list(
    y ~ xd + f1 + f2,
    y ~ xd7 + f1,
    y ~ f1 + f3 + xd,
    y ~ xc * f1 + f2,
    y ~ xs + f1,
    y ~ xs * f2 + f1,
    y ~ xt + f1,
    y ~ xc + f1 + f3,
    y ~ log(xd) + I(xc^2) + f1,
    y ~ xb + xs,
    y ~ xb * f2,
    y ~ xb * xs,
    y ~ xb * xs + xt + f1,
    y ~ xs + xt + xb + f3 + f1)
for (fo in formulas) {
    for (wt in c(FALSE, TRUE)) {
        info <- paste(deparse(fo), if (wt) "weighted" else "")
        l <- if (wt) lm(fo, d, weights = w) else lm(fo, d)
        f <- if (wt) design_lm(fo, d, weights = w) else design_lm(fo, d)
        check_fit(f, l, info)
    }
}

# qr() directly: rank, pivot and the factor of lm()'s QR, up to sign
d0 <- na.omit(d)
for (fo in list(y ~ xs + f1, y ~ xb * xs, y ~ xs * f2 + f1, y ~ xt + f1, y ~ xb * f2 + xt)) {
    info <- deparse(fo)
    D <- design_matrix(fo, d0)
    X <- as.matrix(D)
    q <- qr(D)
    q0 <- qr(X)
    expect_identical(q$rank, q0$rank, info = paste(info, "rank"))
    expect_identical(q$pivot, q0$pivot, info = paste(info, "pivot"))
    k <- q$rank
    i <- seq_len(k)
    expect_true(same(abs(qr.R(q)[i, i]), abs(qr.R(q0)[i, i])), info = paste(info, "R"))
    V <- chol2inv(qr.R(q)[i, i])
    V0 <- chol2inv(qr.R(q0)[i, i])
    expect_true(same(V, V0), info = paste(info, "(R'R)^-1"))
    expect_true(same(qr.coef(q, d0$y), qr.coef(q0, d0$y)), info = paste(info, "coef"))
    expect_true(same(qr.fitted(q, d0$y), qr.fitted(q0, d0$y)), info = paste(info, "fitted"))
}

# variables computed from a numeric variable far from zero (I(x^2) is a
# column of the model frame, not rebuilt from the shifted x)
set.seed(3)
dc <- data.frame(x = runif(5000) + 100, f = gl(5, 1000))
dc$y <- 1 + dc$x + rnorm(5000)
for (fo in list(y ~ f + x + I(x^2), y ~ f * log(x) + x)) {
    check_fit(design_lm(fo, dc), lm(fo, dc), deparse(fo))
    D <- design_matrix(fo, dc)
    expect_true(same(qr.coef(qr(D), dc$y), qr.coef(qr(as.matrix(D)), dc$y), 1e-8), info = deparse(fo))
}

# design_glm() (gaussian) takes the same factorization
l <- glm(y ~ xb * xs + f1, data = d)
g <- design_glm(y ~ xb * xs + f1, data = d)
expect_true(same(coef(g), coef(l)), info = "glm coef")
expect_true(same(vcov(g), vcov(l)), info = "glm vcov")

# centering leaves the user's random numbers alone
set.seed(1)
invisible(qr(design_matrix(y ~ xb * xs, d0)))
a <- runif(1)
set.seed(1)
expect_identical(runif(1), a)

# well-conditioned designs are not centered (no extra pass)
q <- qr(design_matrix(y ~ xc * f1 + xd, d0))
expect_null(q$centered)
q <- qr(design_matrix(y ~ xb * xs, d0))
expect_false(is.null(q$centered))

# in memory and remote agree (the remote fit centers too)
if (requireNamespace("duckdb", quietly = TRUE) && requireNamespace("DBI", quietly = TRUE)) {
    con <- DBI::dbConnect(duckdb::duckdb())
    DBI::dbWriteTable(con, "d", d)
    for (fo in list(y ~ xs + f1, y ~ xb * xs, y ~ xt + f1, y ~ xb * f2)) {
        info <- paste(deparse(fo), "remote")
        f <- design_lm(fo, d)
        r <- design_lm(fo, con, table = "d")
        expect_identical(is.na(coef(f)), is.na(coef(r)), info = paste(info, "aliased"))
        expect_true(same(coef(f), coef(r)), info = paste(info, "coef"))
        expect_true(same(vcov(f), vcov(r)), info = paste(info, "vcov"))
    }
    DBI::dbDisconnect(con, shutdown = TRUE)
}

# The pieces of the FE estimator (R/fe.R, design.matrix#6sq4): the
# pivoted Cholesky kernel shared with qr() (dm_pivchol(): per-column
# thresholds, tie mode, suspect band), the absorption operator
# (dm_absorb(): project, effects, scores, solve) against dense dummies, the
# refinement driver (dm_refine()), and the un-centering helpers.
library(design.matrix)
ns <- asNamespace("design.matrix")

# The kernel ----------------------------------------------------------------
# FE's rule (ties drop) is the sequential rule of the former R loop
# dm_fe_pivchol(), kept here as the reference
old_pivchol <- function(M, scale, tol) {
    keep <- integer(0)
    Rf <- matrix(0, 0, 0)
    for (j in seq_len(ncol(M))) {
        if (length(keep)) {
            r <- backsolve(Rf, M[keep, j], transpose = TRUE)
            d <- M[j, j] - sum(r^2)
        } else {
            r <- numeric(0)
            d <- M[j, j]
        }
        if (is.finite(d) && scale[j] > 0 && d > tol * scale[j]) {
            k <- length(keep)
            Rn <- matrix(0, k + 1L, k + 1L)
            if (k) {
                Rn[seq_len(k), seq_len(k)] <- Rf
                Rn[seq_len(k), k + 1L] <- r
            }
            Rn[k + 1L, k + 1L] <- sqrt(d)
            Rf <- Rn
            keep <- c(keep, j)
        }
    }
    list(keep = keep, R = Rf)
}
set.seed(6)
for (case in 1:20) {
    n <- 200
    p <- sample(1:8, 1)
    X <- matrix(rnorm(n * p), n, p)
    if (p > 2) X[, p] <- X[, 1] - 2 * X[, 2]            # aliased
    if (p > 3) X[, 3] <- X[, 3] * 1e-6 + X[, 1]          # nearly aliased
    if (p > 4) X[, 4] <- 0                               # zero column
    M <- crossprod(X)
    scale <- diag(M) * runif(p, 1, 3)
    for (tol in c(1e-10, 1e-7, 1e-3)) {
        a <- ns$dm_fe_pivchol(M, scale, tol)
        b <- old_pivchol(M, scale, tol)
        expect_identical(a$keep, b$keep, info = paste(case, tol))
        expect_equal(a$R, b$R, tolerance = 1e-13, info = paste(case, tol))
    }
}
# no column
z <- ns$dm_fe_pivchol(matrix(0, 0, 0), numeric(0), 1e-10)
expect_identical(z$keep, integer(0))
expect_identical(dim(z$R), c(0L, 0L))
# tie mode: d == tol * scale exactly (tol = 0.5, scale = 8, d = 4)
M <- diag(c(1, 4))
expect_identical(ns$dm_pivchol(M, c(1, 8), 0.5, "keep")$rank, 2L)
expect_identical(ns$dm_pivchol(M, c(1, 8), 0.5, "drop")$rank, 1L)
# suspect band: listed, judged by recheck() (which sees the factor of the
# columns taken before), and the factorization redone
seen <- NULL
ch <- ns$dm_pivchol(M, c(1, 8), 0.5, "drop", band = list(lo = c(-1, 3), hi = c(-1, 5)),
                    recheck = function(ch, cols) {
                        seen <<- c(seen, cols)
                        list(dres = 4.5, taken = TRUE)
                    })
expect_identical(seen, 2L)
expect_identical(ch$rank, 2L)
expect_equal(ch$R[2, 2], sqrt(4.5))
ch <- ns$dm_pivchol(M, c(1, 8), 0.5, "drop", band = list(lo = c(-1, 3), hi = c(-1, 5)),
                    recheck = function(ch, cols) list(dres = 3.5, taken = FALSE))
expect_identical(ch$rank, 1L)
expect_identical(ch$pivot, c(1L, 2L))
expect_null(ns$dm_pivchol(M, c(1, 8), 0.5, "drop", band = list(lo = c(-1, 3), hi = c(-1, 5)),
                          recheck = function(ch, cols) NULL))
# lm()'s rule (qr()): rank and pivot of qr() on aliased designs, a zero
# column moved to the end, lup = min(n, p)
d <- data.frame(x1 = rnorm(50), x2 = rnorm(50), f = gl(5, 10))
d$x3 <- d$x1 + d$x2
d$x0 <- 0
for (fo in list(~ x1 + x3 + x2 + f, ~ x0 + x1 + f, ~ x1 * f + x2 * f)) {
    D <- design_matrix(fo, d)
    q <- qr(D)
    q0 <- qr(as.matrix(D))
    expect_identical(q$rank, q0$rank, info = deparse(fo))
    expect_identical(q$pivot, q0$pivot, info = deparse(fo))
}
D <- design_matrix(~ x1 * f + x2 * f, d[1:7, ])
expect_identical(qr(D)$rank, qr(as.matrix(D))$rank)
expect_identical(qr(D)$pivot, qr(as.matrix(D))$pivot)

# The absorption operator against dense dummies -----------------------------
set.seed(7)
n <- 400
d <- data.frame(a = sample(30, n, TRUE), b = sample(8, n, TRUE), cl = sample(12, n, TRUE),
                x1 = rnorm(n), x2 = rnorm(n), w = runif(n, 0.5, 2))
d$y <- d$x1 - d$x2 + d$a / 10 + rnorm(n)
D <- design_matrix(~ x1 + x2, d)
xcols <- 2:3
codes <- lapply(d[c("a", "b")], function(v) ns$dmc_codes(match(v, sort(unique(v))), length(unique(v))))
cl <- ns$dmc_codes(d$cl, 12L)
Fd <- cbind(model.matrix(~ factor(a) - 1, d), model.matrix(~ factor(b) - 1, d))
X <- as.matrix(D)[, xcols]
for (w in list(NULL, d$w)) {
    info <- if (is.null(w)) "unweighted" else "weighted"
    ww <- if (is.null(w)) rep(1, n) else w
    A <- ns$dm_absorb(D, d$y, w, codes, xcols, ns$dm_fe_opts())
    expect_identical(unname(A$G), c(30L, 8L), info = info)
    Z <- cbind(X, d$y)
    M <- crossprod(Z, ww * Z)
    pj <- A$project(M)
    ## M[X z] by weighted least squares on the dummies
    sw <- sqrt(ww)
    MZ <- qr.resid(qr(sw * Fd), sw * Z)
    expect_equal(unname(pj$M), unname(crossprod(MZ)), tolerance = 1e-9, info = paste(info, "project"))
    ## the variances before absorption
    expect_equal(unname(A$variance(diag(M)[1:2])), unname(colSums(ww * sweep(X, 2, colSums(ww * X) / sum(ww))^2)),
                 tolerance = 1e-10, info = paste(info, "variance"))
    ## effects: the FE part of the fit at b
    b <- unname(qr.coef(qr(MZ[, 1:2]), MZ[, 3]))
    fe <- A$effects(b)
    fsum <- fe[[1L]][codes[[1L]]$codes] + fe[[2L]][codes[[2L]]$codes]
    cf <- qr.coef(qr(sw * Fd), sw * (d$y - X %*% b))
    cf[is.na(cf)] <- 0
    expect_equal(fsum, unname(drop(Fd %*% cf)), tolerance = 1e-9, info = paste(info, "effects"))
    ## scores: sum within clusters of w e (MX), MX the residual of X on the
    ## dummies (weighted projection), corrected for the FE not nested
    e <- d$y - X %*% b - fsum
    we <- drop(ww * e)
    MX <- MZ[, 1:2] / sw
    U <- rowsum(we * MX, d$cl, reorder = TRUE)
    sc <- A$scores(we, cl)
    expect_identical(sc$nested, c(FALSE, FALSE), info = paste(info, "nested"))
    expect_equal(unname(sc$U), unname(U), tolerance = 1e-8, info = paste(info, "scores"))
    ## solve: K^- B for any B in the range of K (here F'W v)
    v <- rnorm(n)
    Bv <- lapply(codes, function(g) matrix(rowsum(ww * v, g$codes, reorder = TRUE)[, 1L], ncol = 1L))
    sv <- A$solve(Bv)
    fv <- sv$A[[1L]][codes[[1L]]$codes, 1L] + sv$A[[2L]][codes[[2L]]$codes, 1L]
    expect_equal(fv, unname(drop(qr.fitted(qr(sw * Fd), sw * v) / sw)), tolerance = 1e-8, info = paste(info, "solve"))
}

# The refinement driver -----------------------------------------------------
# (residual, correct, the caller's stopping rule; r: the residual the last
# step started from)
calls <- 0
rf <- ns$dm_refine(1, function(s) {
    calls <<- calls + 1
    s
}, function(s, r) list(sol = s / 2), steps = 5, stop = function(step, r, it) it == 3)
expect_identical(rf$sol, 1 / 8)
expect_identical(rf$r, 1 / 4)
expect_true(rf$stopped)
expect_identical(calls, 3)
rf <- ns$dm_refine(1, identity, function(s, r) list(sol = s / 2), steps = 2)
expect_identical(rf$sol, 1 / 4)
expect_false(rf$stopped)
expect_null(ns$dm_refine(1, identity, function(s, r) list(sol = s), steps = 0)$r)

# fits with FE: the refinement options, as before the split
ref <- lm(y ~ x1 + x2 + factor(a) + factor(b), d, weights = w)
for (rf in list("auto", 0L, 1L, 2L)) {
    op <- options(design.matrix.fe_refine = rf)
    a <- design_lm(y ~ x1 + x2 | a + b, d, weights = w)
    options(op)
    expect_equal(coef(a), coef(ref)[c("x1", "x2")], tolerance = 1e-10, info = paste("refine", rf))
    expect_equal(sigma(a), sigma(ref), tolerance = 1e-10, info = paste("refine", rf))
    expect_true(a$fe_info$resid_pass %in% c("done", "skipped"), info = paste("refine", rf))
    if (!identical(rf, "auto")) expect_identical(a$fe_info$resid_pass, "done", info = paste("refine", rf))
}

# Un-centering helpers --------------------------------------------------------
Tm <- diag(3)
Tm[1, 2] <- 5
expect_true(ns$dm_uncenter_ok(Tm, c(1L, 2L), 3L))
expect_false(ns$dm_uncenter_ok(Tm, 2L, 1L))
expect_true(ns$dm_uncenter_ok(Tm, integer(0), 1:3))
V <- crossprod(matrix(rnorm(9), 3))
expect_equal(ns$dm_uncenter_cov(Tm, V), Tm %*% V %*% t(Tm))

# The suspect band of FE (design.matrix#6sq4, part b) -------------------------
# x2 = x1 + a e, e orthogonal to the FE and x1 (in the weighted inner
# product), a chosen so that the residual variance of x2 after the FE and
# x1 is m * tol times its variance. Near the threshold (m within about 2%
# of 1) the decision is taken from the data: it follows the exact ratio,
# computed here by a dense QR on the dummies, in memory and remotely.
remote <- requireNamespace("DBI", quietly = TRUE) && requireNamespace("duckdb", quietly = TRUE)
if (remote) con <- DBI::dbConnect(duckdb::duckdb())
set.seed(8)
n <- 2000
d <- data.frame(a = sample(50, n, TRUE), b = sample(10, n, TRUE), x1 = rnorm(n), e0 = rnorm(n),
                w = runif(n, 0.5, 2))
d$y <- d$x1 + d$a / 20 + rnorm(n)
Fd <- cbind(model.matrix(~ factor(a) - 1, d), model.matrix(~ factor(b) - 1, d))
tol <- 1e-10
for (wt in c(FALSE, TRUE)) {
    ww <- if (wt) d$w else rep(1, n)
    sw <- sqrt(ww)
    Q <- qr(sw * cbind(Fd, d$x1))
    e <- qr.resid(Q, sw * d$e0) / sw
    wvar <- function(v) sum(ww * (v - sum(ww * v) / sum(ww))^2)
    ratio <- function(v) sum(qr.resid(Q, sw * v)^2) / wvar(v)
    for (m in c(0.5, 0.98, 0.995, 0.999, 1.001, 1.005, 1.02, 2)) {
        t <- m * tol
        d$x2 <- d$x1 + sqrt(t * wvar(d$x1) / (sum(ww * e^2) * (1 - t))) * e
        r <- ratio(d$x2)
        info <- paste("band", m, if (wt) "weighted")
        expect_true(abs(r / t - 1) < 1e-4, info = paste(info, "construction"))
        fit <- if (wt) design_lm(y ~ x1 + x2 | a + b, d, weights = w) else design_lm(y ~ x1 + x2 | a + b, d)
        expect_identical(is.na(coef(fit)[["x2"]]), r <= tol, info = info)
        expect_false(is.na(coef(fit)[["x1"]]), info = info)
        ## the same decision for every solver of the reduced system
        for (s in c("dense", "pcg")) {
            op <- options(design.matrix.fe_solver = s)
            fs <- if (wt) design_lm(y ~ x1 + x2 | a + b, d, weights = w) else design_lm(y ~ x1 + x2 | a + b, d)
            options(op)
            expect_identical(is.na(coef(fs)), is.na(coef(fit)), info = paste(info, s))
        }
        ## remotely (other rounding of the sums): the same decision
        if (remote && m %in% c(0.995, 1.005)) {
            DBI::dbWriteTable(con, "dband", d, overwrite = TRUE)
            fr <- if (wt) design_lm(y ~ x1 + x2 | a + b, con, table = "dband", weights = w)
                  else design_lm(y ~ x1 + x2 | a + b, con, table = "dband")
            expect_identical(is.na(coef(fr)), is.na(coef(fit)), info = paste(info, "remote"))
        }
        ## clear cases: the NA pattern of feols()
        if (m %in% c(0.5, 2) && requireNamespace("fixest", quietly = TRUE)) {
            b <- fixest::feols(y ~ x1 + x2 | a + b, d, weights = if (wt) ~w, notes = FALSE, warn = FALSE)
            expect_identical(names(coef(fit))[!is.na(coef(fit))], names(coef(b)), info = paste(info, "feols"))
        }
    }
}
if (remote) DBI::dbDisconnect(con, shutdown = TRUE)

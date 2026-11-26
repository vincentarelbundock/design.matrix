source("helpers.R")
library(design.matrix)

# compare numbers to rounding, with NA/NaN in the same places
close_enough <- function(a, b, tol = 1e-10) {
    a <- unname(as.matrix(a)); b <- unname(as.matrix(b))
    if (!identical(dim(a), dim(b))) return(FALSE)
    if (!identical(is.na(a), is.na(b))) return(FALSE)
    ok <- !is.na(a)
    if (!identical(a[ok] == Inf, b[ok] == Inf) || !identical(a[ok] == -Inf, b[ok] == -Inf)) return(FALSE)
    fin <- ok & is.finite(a)
    isTRUE(all(abs(a[fin] - b[fin]) <= tol * pmax(1, abs(b[fin]))))
}

check_design <- function(D, info) {
    X <- as.matrix(D)
    p <- ncol(X); n <- nrow(X)
    b <- rnorm(p)
    B <- matrix(rnorm(2 * p), p)
    w <- rnorm(n)
    W <- matrix(rnorm(3 * n), n)
    V <- crossprod(matrix(rnorm(p * p), p))
    g <- sample(c("a", "b", NA), n, TRUE)
    expect_true(close_enough(D %*% b, X %*% b), info = paste(info, "mult"))
    expect_true(close_enough(D %*% B, X %*% B), info = paste(info, "mult matrix"))
    expect_true(close_enough(crossprod(D, w), crossprod(X, w)), info = paste(info, "crossprod"))
    expect_true(close_enough(crossprod(D, W), crossprod(X, W)), info = paste(info, "crossprod matrix"))
    expect_true(close_enough(suppressWarnings(colSums(w * D, group = g)), suppressWarnings(rowsum(w * X, g))),
                info = paste(info, "rowsum"))
    expect_true(close_enough(colSums(D), colSums(X)), info = paste(info, "colSums"))
    expect_true(close_enough(rowMeans(D), rowMeans(X)), info = paste(info, "rowMeans"))
    expect_true(close_enough(quadratic_diagonal(D, V), rowSums((X %*% V) * X)), info = paste(info, "quad"))
    expect_true(close_enough(crossprod(D), crossprod(X)), info = paste(info, "gram"))
    expect_true(close_enough(crossprod(D, w * D), crossprod(X, w * X)), info = paste(info, "weighted gram"))
    expect_true(close_enough(crossprod(w * D), crossprod(w * X)), info = paste(info, "scaled gram"))
}

# crossproducts of two different designs
set.seed(5)
d <- fuzz_data(300, special = TRUE)
for (fo in list(y ~ x1 * f1 + f2, y ~ poly(z, 3) * ch + m:lg, ~ 0 + x1:x2:f1)) {
    D0 <- design_matrix(fo, d)
    D1 <- design_matrix(D0, transform(d, x1 = x1 + 1, z = z * 1.1))
    X0 <- as.matrix(D0); X1 <- as.matrix(D1)
    w <- rnorm(nrow(d))
    expect_true(close_enough(crossprod(D0, D1), crossprod(X0, X1)), info = deparse(fo))
    expect_true(close_enough(crossprod(D1 - D0, w * D0), crossprod(X1 - X0, w * X0)), info = deparse(fo))
    expect_true(close_enough(crossprod(D1 - D0), crossprod(X1 - X0)), info = deparse(fo))
}

set.seed(4)
d <- fuzz_data(400, special = TRUE)
forms <- list(y ~ x1 * f1 + f2, ~ 0 + f1:f2:o1, y ~ poly(z, 3) * ch + m:lg,
              y ~ splines::ns(z, 3):f2 + i1, ~ 1, ~ 0 + x1:x2:f1,
              y ~ f1 * f2, y ~ x1 + y)
for (fo in forms) {
    D <- suppressWarnings(design_matrix(fo, d, contrasts.arg = list(f2 = "contr.sum", o1 = "contr.helmert")))
    check_design(D, deparse(fo))
}

# NA factor codes and non-finite numerics (plan with na.pass)
dn <- fuzz_data(300, special = TRUE, na = TRUE)
p <- design_matrix(y ~ x1 * f1 + x2:o1 + ch, dn[stats::complete.cases(dn), ])
D <- design_matrix(p, dn)
expect_true(anyNA(as.matrix(D)))
check_design(D, "NA rows")

# multithreaded kernels give the same values
old <- design_threads()
big <- fuzz_data(10000)
D <- design_matrix(y ~ x1 * f1 + poly(z, 2) + ch, big)
b <- rnorm(ncol(D)); V <- crossprod(matrix(rnorm(ncol(D)^2), ncol(D)))
design_threads(1); m1 <- D %*% b; q1 <- quadratic_diagonal(D, V)
design_threads(4); m4 <- D %*% b; q4 <- quadratic_diagonal(D, V)
design_threads(old)
expect_true(strict_identical(m1, m4))
expect_true(strict_identical(q1, q4))

# Jacobians
set.seed(5)
tr <- fuzz_data(500, special = FALSE)
p <- design_matrix(y ~ x1 * f1 + poly(z, 2) + ch, tr)
X0 <- design_matrix(p, tr)
X1 <- design_matrix(p, transform(tr, x1 = x1 + 1))
s0 <- runif(500); s1 <- runif(500)
J <- s1 * X1 - s0 * X0
Jd <- as.matrix(J)
expect_true(close_enough(Jd, s1 * as.matrix(X1) - s0 * as.matrix(X0)))
expect_identical(dim(J), dim(X0))
expect_identical(colnames(J), colnames(X0))
check_design(J, "jacobian")
J1 <- s0 * X0
check_design(J1, "jacobian one design")
expect_stdout(print(J), "combination of two designs")

# dense fallback: non-formula objects
fit <- stats::lm(y ~ i1 + f1, data = tr)
Dd <- design_matrix(fit)
expect_identical(Dd$path, "dense")
check_design(Dd, "dense fallback")
check_design(s0 * Dd, "dense jacobian")

# several columns of B: each column is bit for bit the product with that
# column alone (same order of the sums), on one or several threads
set.seed(5)
dm <- fuzz_data(300, special = TRUE)
Dm <- design_matrix(~ x1 * f1 + f2 + m, dm)
Bm <- matrix(rnorm(ncol(Dm) * 7), ncol(Dm))
by_column <- vapply(seq_len(ncol(Bm)), function(k) drop(Dm %*% Bm[, k]), numeric(nrow(Dm)))
for (nt in c(1, 4)) {
    design_threads(nt)
    expect_true(strict_identical(Dm %*% Bm, by_column), info = paste("threads", nt))
    Jm <- runif(300) * Dm - Dm
    expect_true(strict_identical(Jm %*% Bm, vapply(seq_len(ncol(Bm)), function(k) drop(Jm %*% Bm[, k]), numeric(300))))
}
design_threads(1)

# A plan can be saved after its contrast rows have been populated and reused
# with a data set whose level codes differ. Older setups lack that cache.
d <- data.frame(f = factor(rep(letters[1:3], 4)), x = seq_len(12))
plan <- design_matrix(~ f * x, d, contrasts.arg = list(f = "contr.sum"))
D <- design_matrix(plan, d)
beta <- seq_len(ncol(D))
invisible(D %*% beta)
restored <- unserialize(serialize(plan, NULL))
nd <- transform(d, f = rev(f))
D <- design_matrix(restored, nd)
ref <- unname(stats::model.matrix(~ f * x, nd, contrasts.arg = list(f = "contr.sum")) %*% beta)
expect_identical(D %*% beta, ref)
D$setup$contrast_rows <- NULL
expect_identical(D %*% beta, ref)

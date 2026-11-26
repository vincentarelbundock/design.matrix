source("helpers.R")
library(design.matrix)

# The column-structured kernels (used for finite designs) agree with the
# dense computation and with the row kernels, up to rounding.
close_to <- function(a, b, tol = 1e-10) {
    a <- unname(as.matrix(a)); b <- unname(as.matrix(b))
    identical(dim(a), dim(b)) &&
        isTRUE(all(abs(a - b) <= tol * pmax(1, abs(b))))
}
rows_only <- function(expr) {
    old <- options(design.matrix.structured = FALSE)
    on.exit(options(old))
    expr
}

set.seed(21)
used <- 0
for (r in 1:150) {
    d <- fuzz_data(sample(c(5, 40, 300), 1), special = FALSE)
    fo <- fuzz_formula()
    D <- tryCatch(suppressWarnings(design_matrix(fo, d, contrasts.arg = fuzz_contrasts(d))),
                  error = function(e) NULL)
    if (is.null(D) || ncol(D) == 0L) next
    X <- as.matrix(D)
    if (!all(is.finite(X)) || !design.matrix:::dm_struct_ok(D)) next
    used <- used + 1
    info <- paste(r, deparse(fo))
    n <- nrow(X); p <- ncol(X)
    w <- runif(n); b <- rnorm(p); B <- matrix(rnorm(3 * p), p); Y <- matrix(rnorm(2 * n), n)
    expect_true(close_to(crossprod(D), crossprod(X)), info = paste(info, "gram"))
    expect_true(close_to(crossprod(D, w * D), crossprod(X, w * X)), info = paste(info, "weighted gram"))
    expect_true(close_to(crossprod(sqrt(w) * D), crossprod(sqrt(w) * X)), info = paste(info, "scaled gram"))
    expect_true(close_to(crossprod(D, Y), crossprod(X, Y)), info = paste(info, "xty"))
    expect_true(close_to(crossprod(w * D, Y), crossprod(w * X, Y)), info = paste(info, "scaled xty"))
    expect_true(close_to(D %*% b, X %*% b), info = paste(info, "mult"))
    expect_true(close_to(D %*% B, X %*% B), info = paste(info, "mult matrix"))
    expect_true(close_to((w * D) %*% b, (w * X) %*% b), info = paste(info, "scaled mult"))
    expect_true(close_to(crossprod(D), rows_only(crossprod(D))), info = paste(info, "gram vs rows"))
    expect_true(close_to(D %*% b, rows_only(D %*% b)), info = paste(info, "mult vs rows"))
}
expect_true(used > 50)

# Designs the structured kernels decline: missing factor codes, non-finite
# or overflowing values, non-finite weights or scales. The row kernels then
# reproduce NA, NaN and Inf where the dense computation has them.
d <- data.frame(x = c(1, Inf, 3, 4), f = factor(c("a", "b", NA, "a")), y = 1:4)
p <- design_matrix(y ~ x + f, d[c(1, 4), ])
D <- design_matrix(p, d)
expect_false(design.matrix:::dm_struct_ok(D))
D2 <- design_matrix(y ~ x * f, data.frame(x = c(1e200, 2, 3), f = gl(3, 1), y = 1:3))
expect_false(design.matrix:::dm_struct_ok(D2))
D3 <- design_matrix(y ~ x * f, data.frame(x = c(1, 2, 3, 4), f = gl(2, 2), y = 1:4))
expect_true(design.matrix:::dm_struct_ok(D3))
X3 <- as.matrix(D3)
w <- c(1, NA, 2, 3)
expect_identical(is.na(crossprod(D3, w * D3)), is.na(crossprod(X3, w * X3)))
expect_identical(is.na(crossprod(D3, w)), is.na(crossprod(X3, w)))
s <- c(1, Inf, 1, 1)
expect_identical(unname(is.nan((s * D3) %*% c(1, 1, 1, 1))), unname(is.nan((s * X3) %*% c(1, 1, 1, 1))))

# the option turns the structured kernels off
expect_false(rows_only(design.matrix:::dm_struct_ok(D3)))

# The model frame skips na.omit()'s copy when nothing is missing; results
# stay identical to model.matrix() under every standard na.action, with and
# without missing values, and with a custom na.action.
mm_same <- function(fo, d, ...) {
    a <- capture(function() as.matrix(design_matrix(fo, d, ...)))
    b <- capture(function() model.matrix(fo, d, ...))
    same_outcome(a, b)
}
d0 <- data.frame(y = rnorm(6), x = rnorm(6), f = gl(3, 2), ch = letters[c(1, 2, 1, 2, 1, 2)])
d1 <- d0; d1$x[2] <- NA; d1$f[5] <- NA
d2 <- d0; rownames(d2) <- paste0("r", 1:6)
custom <- function(object, ...) object[-1, , drop = FALSE]
for (naa in list("na.omit", "na.exclude", "na.fail", "na.pass", custom)) {
    old <- options(na.action = naa)
    for (dd in list(d0, d1, d2)) {
        expect_true(mm_same(y ~ x * f + ch, dd), info = paste(deparse(naa)[1], nrow(dd)))
        expect_true(mm_same(~ log(abs(x)) + f, dd), info = paste(deparse(naa)[1], "no response"))
    }
    options(old)
}
d3 <- d1; attr(d3, "na.action") <- "na.fail"
expect_true(mm_same(y ~ x + f, d3))

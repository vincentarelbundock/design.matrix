# Regression tests for the findings of Roborev job 1262.
source("helpers.R")
library(design.matrix)

close_na <- function(a, b, tol = 1e-10) {
    a <- unname(as.matrix(a)); b <- unname(as.matrix(b))
    identical(dim(a), dim(b)) && identical(is.na(a), is.na(b)) &&
        identical(is.nan(a), is.nan(b)) &&
        isTRUE(all.equal(a[!is.na(a)], b[!is.na(b)], tolerance = tol))
}

# 1. A term with a zero-column matrix has no columns and writes nothing
d <- list(x = as.double(1:4), m = matrix(numeric(0), 4, 0))
attr(d, "row.names") <- c(NA_integer_, -4L)
class(d) <- "data.frame"
tt <- stats::terms(~ m:x + x)
attr(d, "terms") <- tt
# terms are reordered to x + m:x; m:x has no columns (R itself writes the
# x column of that term past the end of its matrix here)
X <- suppressWarnings(as.matrix(design_matrix(tt, d)))
expect_identical(dim(X), c(4L, 2L))
expect_identical(colnames(X), c("(Intercept)", "x"))
expect_identical(unname(X[, "x"]), as.double(1:4))
w <- capture(function() as.matrix(design_matrix(tt, d)))$warnings
expect_identical(w, "problem with term 2 in model.matrix: no columns are assigned")

# 2. The fast path evaluates the variables together, as model.frame() does
dd <- data.frame(x = rnorm(10), z = rnorm(10))
fo <- ~ I({ z <- x + 1; z }) + z
p <- design.matrix:::design_plan(fo, dd)
nd <- dd[3:6, ]
ref <- stats::model.matrix(p$terms, stats::model.frame(p$terms, nd, xlev = p$xlevels,
                           na.action = stats::na.pass), contrasts.arg = p$contrasts)
expect_true(strict_identical(as.matrix(design_matrix(p, nd)), ref))

# 3. Jacobian entries of the two designs are combined before products
d3 <- data.frame(a = c(2^54, 2^54), b = c(1, 1))
D3 <- design_matrix(~ 0 + a + b, d3)
J3 <- c(1, 1) * D3 - D3
expect_identical(drop(J3 %*% c(1, 1)), c(0, 0))
expect_identical(crossprod(J3, c(1, 1))[, 1], c(a = 0, b = 0))
expect_identical(quadratic_diagonal(J3, diag(2)), c(0, 0))

# 4. Non-finite weights reach every column, as 0 * NA does in the dense product
d4 <- data.frame(f = factor(c("a", "b", "a", "c")), x = c(1, 2, 3, 4))
D4 <- design_matrix(~ f + x, d4)
X4 <- as.matrix(D4)
for (bad in c(NA, Inf, NaN)) {
    w4 <- c(bad, 1, 2, 3)
    expect_true(close_na(crossprod(D4, w4), crossprod(X4, w4)), info = paste("weight", bad))
    expect_true(close_na(t(colSums(w4 * D4, group = c(1, 1, 2, 2))),
                         cbind(crossprod(X4[1:2, ], w4[1:2]), crossprod(X4[3:4, ], w4[3:4]))))
}

# 5. Non-finite scales give NaN where the dense Jacobian has Inf * 0
s5 <- c(Inf, 1, -Inf, 2)
J5 <- s5 * D4
expect_true(close_na(as.matrix(J5), s5 * X4))
b5 <- c(1, 2, 3, 4)
expect_true(close_na(J5 %*% b5, as.matrix(J5) %*% b5))
expect_true(close_na(crossprod(J5, rep(1, 4)), crossprod(as.matrix(J5), rep(1, 4))))
J5b <- rep(1, 4) * D4 + s5 * D4
expect_true(close_na(J5b %*% b5, as.matrix(J5b) %*% b5))

# 6. Products that overflow before meeting a zero contrast give NaN
d6 <- data.frame(x = c(1e200, 1, 2), z = c(1e200, 1, 3), f = factor(c("a", "b", "c")))
D6 <- design_matrix(~ x:z:f, d6)
X6 <- as.matrix(D6)
expect_true(any(is.nan(X6)))
b6 <- seq_len(ncol(X6))
expect_true(close_na(D6 %*% b6, X6 %*% b6))
expect_true(close_na(crossprod(D6, c(1, 2, 3)), crossprod(X6, c(1, 2, 3))))
d6b <- data.frame(x = c(1e300, 1), f = factor(c("a", "b")), g = factor(c("u", "v")))
D6b <- design_matrix(~ 0 + x:f:g, d6b, contrasts.arg = list(g = matrix(c(1e10, 0), 2, 1)))
X6b <- as.matrix(D6b)
expect_true(close_na(D6b %*% rep(1, ncol(X6b)), X6b %*% rep(1, ncol(X6b))))

# 7. The quadratic form gives NaN when X V overflows in an omitted column
d7 <- data.frame(f = factor(c("a", "b", "c")), x = c(1e300, 1, 1))
D7 <- design_matrix(~ 0 + f + x, d7)
X7 <- as.matrix(D7)
V7 <- matrix(1e10, 4, 4)
expect_true(close_na(quadratic_diagonal(D7, V7), rowSums((X7 %*% V7) * X7)))

# 8. Formula subclasses with their own method keep their dispatch, and
#    omitted arguments stay omitted
registerS3method("model.matrix", "dm_test_formula",
                 function(object) "custom method", envir = globalenv())
fo8 <- structure(y ~ x, class = c("dm_test_formula", "formula"))
expect_identical(as.matrix(design_matrix(fo8)), "custom method")
expect_identical(design_matrix(fo8)$path, "dense")
fit <- stats::lm(x ~ z, data = dd)
expect_true(strict_identical(as.matrix(design_matrix(fit)), stats::model.matrix(fit)))

# 9. Zero-column dense designs keep their dimensions
fit0 <- stats::lm(x ~ 0, data = dd)
D9 <- design_matrix(fit0)
expect_identical(dim(crossprod(D9, rep(1, 10))), c(0L, 1L))
expect_identical(dim(colSums(D9, group = rep(1:2, 5))), c(2L, 0L))

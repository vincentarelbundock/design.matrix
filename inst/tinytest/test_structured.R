source("helpers.R")
library(xtx)

# The column-structured kernels (used for finite designs) agree with the
# dense computation and with the row kernels, up to rounding.
close_to <- function(a, b, tol = 1e-10) {
    a <- unname(as.matrix(a)); b <- unname(as.matrix(b))
    identical(dim(a), dim(b)) &&
        isTRUE(all(abs(a - b) <= tol * pmax(1, abs(b))))
}
rows_only <- function(expr) {
    old <- options(xtx.structured = FALSE)
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
    if (!all(is.finite(X)) || !xtx:::xtx_struct_ok(D)) next
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
expect_false(xtx:::xtx_struct_ok(D))
D2 <- design_matrix(y ~ x * f, data.frame(x = c(1e200, 2, 3), f = gl(3, 1), y = 1:3))
expect_false(xtx:::xtx_struct_ok(D2))
D3 <- design_matrix(y ~ x * f, data.frame(x = c(1, 2, 3, 4), f = gl(2, 2), y = 1:4))
expect_true(xtx:::xtx_struct_ok(D3))
X3 <- as.matrix(D3)
w <- c(1, NA, 2, 3)
expect_identical(is.na(crossprod(D3, w * D3)), is.na(crossprod(X3, w * X3)))
expect_identical(is.na(crossprod(D3, w)), is.na(crossprod(X3, w)))
s <- c(1, Inf, 1, 1)
expect_identical(unname(is.nan((s * D3) %*% c(1, 1, 1, 1))), unname(is.nan((s * X3) %*% c(1, 1, 1, 1))))

# the option turns the structured kernels off
expect_false(rows_only(xtx:::xtx_struct_ok(D3)))

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

# Roborev 1271 regressions ----------------------------------------------

# 2. The model frame is evaluated once, with or without missing values, and
#    new levels on rows that na.omit() drops are not an error
count <- 0
tick <- function(x) { count <<- count + 1; x }
dd <- data.frame(y = rnorm(6), x = c(1, NA, 3, 4, 5, 6))
count <- 0; invisible(design_matrix(y ~ tick(x), dd)); expect_equal(count, 1)
dd$x[2] <- 2
count <- 0; invisible(design_matrix(y ~ tick(x), dd)); expect_equal(count, 1)
dx <- data.frame(y = rnorm(4), x = c(1, NA, 3, 4), f = factor(c("a", "c", "b", "a")))
expect_true(mm_same(y ~ x + f, dx, xlev = list(f = c("a", "b"))))

# 4-5. Row scales are not combined or moved when that would overflow
D45 <- design_matrix(~ 0 + x, data.frame(x = 1e-160))
X45 <- as.matrix(D45)
expect_equal(unname(crossprod(1e160 * D45)), unname(crossprod(1e160 * X45)))
expect_equal(unname(crossprod(1e160 * D45, 1e160)), unname(crossprod(1e160 * X45, 1e160)))
expect_equal(unname(crossprod(1e160 * D45, 1e160 * D45)), unname(crossprod(1e160 * X45)))
expect_equal(unname(drop((1e160 * D45) %*% 1e160)), unname(drop((1e160 * X45) %*% 1e160)))
expect_true(is.finite(crossprod(1e160 * D45)))

# 7. Combined factor levels beyond int range (1300^3 > 2^31): the structured
#    kernels decline instead of overflowing
lv <- as.character(1:1300) # 1300^3 combined levels > 2^31
d7 <- data.frame(f1 = factor(c("1", "2", "3"), levels = lv), f2 = factor(c("3", "1", "2"), levels = lv),
                 f3 = factor(c("2", "3", "1"), levels = lv))
cm <- matrix(seq_len(1300) / 1300, ncol = 1)
D7 <- design_matrix(~ f1 * f2 * f3, d7, contrasts.arg = list(f1 = cm, f2 = cm, f3 = cm))
expect_false(xtx:::xtx_struct_ok(D7))
X7 <- as.matrix(D7)
expect_equal(unname(crossprod(D7)), unname(crossprod(X7)))

# Roborev 1272 regressions -------------------------------------------------

# 1. When structured sums are not finite, the fallback uses the original
#    scaled operands
D1 <- design_matrix(~ 0 + z + x, data.frame(z = 1e150, x = 1e-150))
X1 <- as.matrix(D1)
expect_equal(unname(crossprod(1e100 * D1)), unname(crossprod(1e100 * X1)))
expect_true(is.finite(crossprod(1e100 * D1)[1, 2])) # the diagonal overflows densely too

# 2. No underflow from moving scales or regrouping products
D2 <- design_matrix(~ 0 + x, data.frame(x = 1e-150))
X2 <- as.matrix(D2)
expect_equal(unname(drop((1e100 * D2) %*% 1e-200)), unname(drop((1e100 * X2) %*% 1e-200)))
expect_true(drop((1e100 * D2) %*% 1e-200) != 0)
D2b <- design_matrix(~ 0 + x, data.frame(x = 1e150))
X2b <- as.matrix(D2b)
expect_equal(unname(crossprod(1e-100 * D2b, 1e-250)), unname(crossprod(1e-100 * X2b, 1e-250)))
expect_true(crossprod(1e-100 * D2b, 1e-250) != 0)

# 3. Tiny numerics times large contrasts: the structured kernels decline,
#    so X'X keeps entries whose squares would underflow if regrouped
d3 <- data.frame(y = rnorm(4), x = c(1, 2, 3, 4) * 1e-200, f = factor(c("a", "b", "a", "b")))
cm3 <- matrix(c(-1e100, 1e100), ncol = 1)
D3 <- design_matrix(y ~ x:f, d3, contrasts.arg = list(f = cm3))
expect_false(xtx:::xtx_struct_ok(D3))
X3 <- as.matrix(D3)
expect_equal(qr(D3)$rank, qr(X3)$rank)

# Factor x factor interactions (design.matrix#fj49): the structured kernels
# take each factor of a pair of terms in turn (selections in place, other
# contrasts by products), with moments by cell of the union of the two
# terms' factors. They agree with the row kernels and with model.matrix(),
# for 2- and 3-way interactions, with and without weights, with empty
# cells and unused levels, and for every kind of contrast.
set.seed(49)
ff_data <- function(n) {
    d <- data.frame(y = rnorm(n), x = rnorm(n), z = runif(n),
                    f1 = factor(sample(letters[1:6], n, TRUE), levels = letters[1:7]), # unused level
                    f2 = factor(sample(LETTERS[1:4], n, TRUE)),
                    f3 = factor(sample(c("u", "v", "w"), n, TRUE)),
                    lg = sample(c(TRUE, FALSE), n, TRUE))
    ## empty cells: f1 = "a" never meets f2 = "A", f3 = "u" never meets f2 = "B"
    d$f2[d$f1 == "a" & d$f2 == "A"] <- "C"
    d$f3[d$f3 == "u" & d$f2 == "B"] <- "v"
    d
}
sel_cm <- matrix(c(0, 2, 0, 0, 0, 0, 1, -1), 4, 2) # a selection: two levels in one column
ff_formulas <- list(y ~ f1 * f2, y ~ f1 * f2 * f3, y ~ f1:f2, y ~ 0 + f1:f2:f3, y ~ f1 * f3 + f2,
                    y ~ f1 + f1:f2, y ~ x * f1 * f2, y ~ x:f1:f2 + f3, y ~ poly(x, 2) * f1 * f3,
                    y ~ f1 * lg * f3, y ~ (f1 + f2 + f3)^2 + z)
ff_contrasts <- list(NULL,
                     list(f1 = "contr.sum", f2 = "contr.helmert", f3 = "contr.poly"),
                     list(f1 = "contr.SAS", f2 = "contr.sum"),
                     list(f2 = sel_cm))
for (n in c(30, 400)) {
    d <- ff_data(n)
    w <- runif(n)
    g <- sample(1:5, n, TRUE)
    for (fo in ff_formulas) for (ca in ff_contrasts) {
        info <- paste(n, deparse(fo), paste(names(ca), collapse = ","))
        D <- suppressWarnings(design_matrix(fo, d, contrasts.arg = ca))
        X <- suppressWarnings(model.matrix(fo, d, contrasts.arg = ca))
        expect_true(xtx:::xtx_struct_ok(D), info = info)
        expect_true(close_to(crossprod(D), crossprod(X)), info = paste(info, "gram"))
        expect_true(close_to(crossprod(D), rows_only(crossprod(D))), info = paste(info, "gram vs rows"))
        expect_true(isSymmetric(unname(crossprod(D))), info = paste(info, "symmetric"))
        expect_true(close_to(crossprod(D, w * D), crossprod(X, w * X)), info = paste(info, "weighted"))
        expect_true(close_to(crossprod(D, w * D), rows_only(crossprod(D, w * D))), info = paste(info, "weighted vs rows"))
        expect_true(close_to(crossprod(sqrt(w) * D), crossprod(sqrt(w) * X)), info = paste(info, "scaled"))
        Y <- cbind(d$y, w)
        expect_true(close_to(crossprod(D, Y), crossprod(X, Y)), info = paste(info, "xty"))
        B <- matrix(rnorm(2 * ncol(X)), ncol(X))
        expect_true(close_to(D %*% B, X %*% B), info = paste(info, "mult"))
        expect_true(close_to(colSums(D, group = g), rowsum(X, g)), info = paste(info, "group sums"))
        gs <- xtx:::xtx_p_colsums_groups(w * D, list(xtx:::xtx_c_codes(g, 5L)))[[1L]]
        expect_true(close_to(gs$S, rowsum(w * X, g)), info = paste(info, "structured group sums"))
        expect_true(close_to(quadratic_diagonal(D, diag(ncol(X))), rowSums(X^2)), info = paste(info, "quad"))
    }
}

source("helpers.R")
library(xtx)

ref <- function(x) unname(rowSums(!is.finite(as.matrix(x))) == 0)

# fuzzed data with Inf, -Inf, 1e300, NA and NaN, and NA factor codes
set.seed(202)
formulas <- list(
    y ~ x1 * f1 + f2,
    y ~ x1 * x2 + o1,           # 1e300 * 1e300 overflows
    y ~ poly(z, 2) * f1 + ch + lg,
    y ~ x1:f2 + m,
    y ~ 0 + f1:x2)
for (k in 1:40) {
    n <- sample(6:40, 1)
    d <- fuzz_data(n, special = TRUE, na = k %% 2 == 0)
    fo <- formulas[[1 + k %% length(formulas)]]
    D <- design_matrix(fo, d)      # NA rows dropped, as by model.matrix()
    expect_identical(xtx:::finite_rows(D), ref(D), info = paste(deparse(fo), "k =", k))
    D <- design_matrix(D, d)       # every row (na.pass)
    s <- rnorm(n)
    if (n >= 3) s[sample(n, 1)] <- sample(c(Inf, NaN, 1e300), 1)
    D1 <- design_matrix(D, transform(d, x2 = x2 * 2))
    for (J in list(s * D, D1 - D, s * D1 + 2 * D)) {
        expect_identical(xtx:::finite_rows(J), ref(J), info = paste(deparse(fo), "k =", k, "(combination)"))
    }
}

# a product of finite variables that overflows, and Inf times a zero
d <- data.frame(a = c(1, 1e200, 2), b = c(1, 1e200, 3), f = factor(c("u", "v", "u")))
expect_identical(xtx:::finite_rows(design_matrix(~ a:b, d)), c(TRUE, FALSE, TRUE))
d2 <- data.frame(x = c(1, Inf, 2), f = factor(c("u", "u", "v")))
expect_identical(xtx:::finite_rows(design_matrix(~ x:f, d2)), ref(design_matrix(~ x:f, d2)))

# threads do not change the answer
d <- fuzz_data(5000, special = TRUE, na = TRUE)
D <- design_matrix(y ~ x1 * f1 + x2 * f2, d)
old <- design_threads()
design_threads(4)
expect_identical(xtx:::finite_rows(D), ref(D))
design_threads(old)

expect_error(xtx:::finite_rows(matrix(1)), "xtx_matrix")

source("helpers.R")
library(design.matrix)
design_row_groups <- design.matrix:::design_row_groups
op <- options(design.matrix.dedup = "auto")

same_partition <- function(a, b) identical(a$map, b$map) && identical(a$first, b$first)

set.seed(41)
base <- fuzz_data(30, special = FALSE)
d <- base[sample(30, 3000, replace = TRUE), ]
p <- design_matrix(y ~ x1 * f1 + ch + poly(z, 2), base)

# design_row_groups(): map, first, counts
g <- design_row_groups(design_matrix(p, d))
expect_true(max(g$map) <= 30)
expect_identical(length(g$first), length(g$counts))
expect_identical(sum(g$counts), nrow(d))
expect_identical(g$map[g$first], seq_along(g$first))
key <- do.call(paste, d[c("x1", "f1", "ch", "z")])
expect_identical(length(g$first), length(unique(key)))

# Jacobian groups derived from the designs' cached groups
X0 <- design_matrix(p, d)
X1 <- design_matrix(p, transform(d, x1 = x1 + 1))
invisible(design_row_groups(X0)); invisible(design_row_groups(X1))
b <- rnorm(ncol(X0), sd = 0.3)
e0 <- drop(X0 %*% b); e1 <- drop(X1 %*% b)  # scales are functions of the row
J <- stats::dlogis(e1) * X1 - stats::dlogis(e0) * X0
derived <- design.matrix:::dm_groups_from_designs(J)
expect_false(is.null(derived))
J_hash <- stats::dlogis(e1) * X1 - stats::dlogis(e0) * X0
direct <- do.call(.Call, c(list(design.matrix:::C_dm_row_groups), design.matrix:::jac_args(J_hash), list(NULL, 1)))
names(direct) <- c("map", "first")
expect_true(same_partition(derived, direct))
V <- crossprod(matrix(rnorm(ncol(X0)^2), ncol(X0)))
q_seeded <- quadratic_diagonal(J, V)
op2 <- options(design.matrix.dedup = FALSE)
q_plain <- quadratic_diagonal(stats::dlogis(e1) * X1 - stats::dlogis(e0) * X0, V)
options(op2)
expect_true(strict_identical(q_seeded, q_plain))

# a scale that differs within a design group: not derived, hashed instead
s_noise <- stats::dlogis(e0) + seq_len(nrow(d)) * 1e-9
J2 <- s_noise * X0
expect_null(design.matrix:::dm_groups_from_designs(J2))
op2 <- options(design.matrix.dedup = FALSE)
ref2 <- quadratic_diagonal(s_noise * X0, V)
options(op2)
expect_true(strict_identical(quadratic_diagonal(J2, V), ref2))

# dense designs are refused
fit <- stats::lm(y ~ x1, data = base)
expect_error(design_row_groups(design_matrix(fit)), "formula or a plan")
expect_error(design_row_groups(1), "design")
options(op)

# Preserve first-occurrence ids as the hash table grows and rehashes.
set.seed(537)
keys <- sample.int(3000L, 20000L, replace = TRUE)
g <- design_row_groups(design_matrix(~ x, data.frame(x = as.double(keys))))
u <- unique(keys)
expect_identical(g$map, match(keys, u))
expect_identical(g$first, match(u, keys))
expect_identical(g$counts, tabulate(g$map))

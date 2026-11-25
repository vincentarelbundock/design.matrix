source("helpers.R")
library(design.matrix)

# Row subsetting: design_matrix(plan, nd)[i, ] is design_matrix(plan, nd[i, ]),
# in the strict sense, on fuzzed data and formulas with transformations,
# interactions, ordered and character factors, logicals and matrix columns.
set.seed(101)
formulas <- list(
    y ~ x1 * f1 + f2,
    y ~ poly(z, 2) * f1 + o1 + ch + lg,
    y ~ x1:f2 + m + log(z),
    y ~ f1 * f2 * ch,
    y ~ 0 + f1 + x2)
for (k in 1:30) {
    n <- sample(8:40, 1)
    d <- fuzz_data(n, special = TRUE, na = k %% 3 == 0)
    fo <- formulas[[1 + k %% length(formulas)]]
    p <- design_matrix(fo, d)
    D <- design_matrix(p, d)
    idx <- list(
        sample(n, sample(1:n, 1)),
        sample(n, n + 3, replace = TRUE),     # repeats
        -sample(n, 2),                        # negative
        sample(c(TRUE, FALSE), n, TRUE),      # logical
        integer(0),
        rev(seq_len(n)))
    for (i in idx) {
        a <- capture(function() as.matrix(design_matrix(p, d[i, , drop = FALSE])))
        b <- capture(function() as.matrix(D[i, ]))
        expect_true(same_outcome(a, b), info = paste(deparse(fo), "k =", k))
    }
}

# from a fitted model, and on the fast and general paths alike
d <- fuzz_data(30, special = FALSE)
fit <- lm(y ~ poly(z, 2) * f1 + ch + x1:f2, d)
D <- design_matrix(fit, newdata = d)
i <- c(5, 1, 30, 5)
expect_true(strict_identical(as.matrix(D[i, ]), as.matrix(design_matrix(fit, newdata = d[i, ]))))

# products and summaries on the subset agree with the subset of the products
b <- rnorm(ncol(D))
V <- crossprod(matrix(rnorm(ncol(D)^2), ncol(D))) / ncol(D)
X <- as.matrix(D)
expect_equal(drop(D[i, ] %*% b), unname(drop(X[i, ] %*% b)))
expect_equal(quadratic_diagonal(D[i, ], V), unname(rowSums((X[i, ] %*% V) * X[i, ])))
w <- rnorm(length(i))
expect_equal(drop(crossprod(D[i, ], w)), drop(crossprod(X[i, ], w)))

# combinations keep their scales for the selected rows
D1 <- design_matrix(fit, newdata = transform(d, z = z + 0.1))
s0 <- runif(30)
s1 <- runif(30)
J <- s1 * D1 - s0 * D
M <- s1 * as.matrix(D1) - s0 * X
expect_equal(unname(as.matrix(J[i, ])), unname(M[i, ]))
expect_equal(quadratic_diagonal(J[i, ], V), unname(rowSums((M[i, ] %*% V) * M[i, ])))
expect_equal(as.matrix((2 * D)[i, ]), 2 * as.matrix(D[i, ]), check.attributes = FALSE)

# grouped column sums and sparse matrices of a subset
g <- sample(c("a", "b", "c"), length(i), TRUE)
expect_equal(colSums(D[i, ], group = g), rowsum(X[i, ], g), check.attributes = FALSE)
expect_equal(colSums(J[i, ], group = g), rowsum(M[i, ], g), check.attributes = FALSE)
expect_true(strict_identical(as(D[i, ], "CsparseMatrix"), as(as.matrix(D[i, ]), "CsparseMatrix")))

# cached groups of identical rows carry over, renumbered, without changing results
dd <- fuzz_data(400, special = FALSE)
dd <- dd[sample(400, 3000, TRUE), ]
Dg <- design_matrix(design_matrix(y ~ f1 * f2 + o1, dd), dd)
Vg <- diag(ncol(Dg))
op <- options(design.matrix.dedup = TRUE)
invisible(quadratic_diagonal(Dg, Vg))
expect_true(is.list(Dg$cache$groups))
ig <- c(sample(3000, 500), 7, 7)
sub <- Dg[ig, ]
expect_true(is.list(sub$cache$groups))
expect_identical(sub$cache$groups$map, match(Dg$cache$groups$map[ig], unique(Dg$cache$groups$map[ig])))
q_carried <- quadratic_diagonal(sub, Vg)
options(design.matrix.dedup = FALSE)
expect_true(strict_identical(q_carried, quadratic_diagonal(Dg[ig, ], Vg)))
options(op)

# dense designs (non-default model.matrix methods) are subset too
Dd <- design.matrix:::dm_dense_design(model.matrix(~ x1 + f1, d))
expect_true(strict_identical(as.matrix(Dd[c(2, 4), ]), model.matrix(~ x1 + f1, d[c(2, 4), ])))

# the subset is a design of its own: reuse, and cached groups are not shared
expect_identical(dim(D[i, ]), c(length(i), ncol(D)))
expect_true(strict_identical(as.matrix(design_matrix(D[1:3, ], d[4:6, ])),
                             as.matrix(design_matrix(fit, newdata = d[4:6, ]))))

# errors
expect_error(D[, 1], "Only rows")
expect_error(D[1:2], "x\\[i, \\]")
expect_error(D[c(1, NA), ], "missing")
expect_error(D[c(1, -2), ], "mix")
expect_error(D[31, ], "out of bounds")
expect_error(D[c(TRUE, FALSE), ], "one per row")
expect_error(D["a", ], "integer or logical")

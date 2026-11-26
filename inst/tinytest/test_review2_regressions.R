# Regression tests for the findings of Roborev job 1265.
source("helpers.R")
library(xtx)

set.seed(11)
d <- fuzz_data(50, special = FALSE)

# 1. A plan refuses a formula whose response is also on the right-hand side
expect_error(xtx:::design_plan(y ~ y + x1, d), "right-hand side")
expect_inherits(xtx:::design_plan(y ~ x1 + y:f1, d), "xtx_plan") # interactions with y are fine

# 2. The stats fallback does not evaluate the caller's arguments twice
calls <- 0
make_data <- function() { calls <<- calls + 1; d }
X <- as.matrix(design_matrix(stats::lm(y ~ x1 + f1, data = make_data())))
expect_equal(calls, 1)
calls <- 0
D <- design_matrix(stats::lm(y ~ x1, data = make_data()))
expect_equal(calls, 1)
# omitted arguments stay omitted; supplied ones are passed on
registerS3method("model.matrix", "xtx_args_formula",
                 function(object, ...) names(list(...)), envir = globalenv())
fo <- structure(y ~ x1, class = c("xtx_args_formula", "formula"))
expect_null(as.matrix(design_matrix(fo)))
expect_identical(as.matrix(design_matrix(fo, data = d)), "data")
expect_identical(as.matrix(design_matrix(fo, xlev = list())), "xlev")

# 3. Fast path falling back after evaluating: transformations run once
counter <- new.env()
counter$n <- 0
tick <- function(x) { counter$n <- counter$n + 1; x }
p <- xtx:::design_plan(y ~ tick(x1) + f1, d)
nd <- d[1:5, ]
nd$f1 <- factor(as.character(nd$f1), levels = rev(levels(d$f1))) # forces fallback
counter$n <- 0
X3 <- as.matrix(design_matrix(p, nd))
expect_equal(counter$n, 1)
ref <- stats::model.matrix(p$terms, stats::model.frame(p$terms, nd, xlev = p$xlevels,
                           na.action = stats::na.pass), contrasts.arg = p$contrasts)
expect_true(strict_identical(X3, ref))
# fallback errors are the reference's errors
nd2 <- d[1:5, ]
nd2$f1 <- factor(c("a", "b", "q", "a", "b"))
expect_true(same_outcome(
    capture(function() as.matrix(design_matrix(p, nd2))),
    capture(function() stats::model.matrix(p$terms, stats::model.frame(
        p$terms, nd2, xlev = p$xlevels, na.action = stats::na.pass),
        contrasts.arg = p$contrasts))))
# a fallback with NA rows dropped by na.omit keeps the reference row names
nd3 <- transform(d[1:6, ], x1 = replace(x1, 2, NA))
nd3$f1 <- factor(as.character(nd3$f1), levels = rev(levels(d$f1)))
expect_true(strict_identical(
    as.matrix(design_matrix(p, nd3, na.action = stats::na.omit)),
    stats::model.matrix(p$terms, stats::model.frame(p$terms, nd3, xlev = p$xlevels,
                        na.action = stats::na.omit), contrasts.arg = p$contrasts)))

# 4. The quadratic form sums like rowSums(): no overflow from cancelling terms
d4 <- data.frame(a = 1, b = 1)
D4 <- design_matrix(~ 0 + a + b, d4)
V4 <- matrix(c(1e308, 1e308, 1e308, -1e308), 2)
X4 <- as.matrix(D4)
expect_identical(quadratic_diagonal(D4, V4), unname(rowSums((X4 %*% V4) * X4)))

# 5. A method returning NULL is returned as is
registerS3method("model.matrix", "xtx_null_formula",
                 function(object, ...) NULL, envir = globalenv())
fo5 <- structure(y ~ x1, class = c("xtx_null_formula", "formula"))
expect_null(as.matrix(design_matrix(fo5)))
expect_null(as.matrix(design_matrix(fo5)))

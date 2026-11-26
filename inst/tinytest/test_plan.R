source("helpers.R")
library(xtx)

reference <- function(plan, newdata, na.action = stats::na.pass) {
    mf <- stats::model.frame(plan$terms, newdata, xlev = plan$xlevels, na.action = na.action)
    stats::model.matrix(plan$terms, mf, contrasts.arg = plan$contrasts)
}

set.seed(2)
train <- fuzz_data(300, special = FALSE)
fo <- y ~ x1 * f1 + poly(z, 2) + ch + lg + o1 + m + splines::ns(x2, 3)
p <- xtx:::design_plan(fo, train)
expect_inherits(p, "xtx_plan")
expect_false(is.null(p$fast))

# fast path: subsets, reordered rows, a single row, NA values (na.pass)
nd_list <- list(train[1:5, ], train[c(9, 3, 3, 200), ], train[7, ],
                transform(train[1:10, ], x1 = replace(x1, 2, NA), ch = replace(ch, 3, NA),
                          f1 = replace(f1, 4, NA), lg = replace(lg, 5, NA)),
                transform(train[1:6, ], x1 = c(Inf, -0, NaN, 1, -Inf, 3)))
for (nd in nd_list) {
    D <- design_matrix(p, nd)
    expect_identical(D$path, "fast")
    expect_true(strict_identical(as.matrix(D), reference(p, nd)))
}

# a factor whose levels differ, a character given as factor: general path
nd <- train[1:8, ]
nd$f1 <- factor(as.character(nd$f1), levels = rev(levels(train$f1)))
D <- design_matrix(p, nd)
expect_identical(D$path, "general")
expect_true(strict_identical(as.matrix(D), reference(p, nd)))
nd <- train[1:8, ]
nd$ch <- factor(nd$ch)
expect_true(strict_identical(as.matrix(design_matrix(p, nd)), reference(p, nd)))

# new levels and missing columns: same errors as the reference
nd <- train[1:4, ]
nd$ch[1] <- "zzz"
expect_true(same_outcome(capture(function() as.matrix(design_matrix(p, nd))),
                         capture(function() reference(p, nd))))
nd <- train[1:4, ]
nd$f1 <- factor(c("a", "b", "q", "a"))
expect_true(same_outcome(capture(function() as.matrix(design_matrix(p, nd))),
                         capture(function() reference(p, nd))))
expect_error(as.matrix(design_matrix(p, train[1:4, setdiff(names(train), "x1")])))

# other na.action: general path, rows dropped as in the reference
nd <- transform(train[1:10, ], x1 = replace(x1, 2, NA))
D <- design_matrix(p, nd, na.action = stats::na.omit)
expect_identical(D$path, "general")
expect_true(strict_identical(as.matrix(D), reference(p, nd, stats::na.omit)))

# contrasts.arg in the plan; plan without response; numeric given as integer
p2 <- xtx:::design_plan(~ f2 * x1 + i1, train, contrasts.arg = list(f2 = "contr.sum"))
nd <- train[5:9, ]
nd$i1 <- as.double(nd$i1)
D <- design_matrix(p2, nd)
expect_identical(D$path, "fast")
expect_true(strict_identical(as.matrix(D), reference(p2, nd)))

# randomized plans
set.seed(3)
for (k in 1:60) {
    tr <- fuzz_data(sample(c(30, 200), 1), special = stats::runif(1) < 0.5)
    fk <- fuzz_formula(response = TRUE)
    # With the response also on the right-hand side, delete.response()
    # leaves a term without variables, and R's C routine returns that
    # column uninitialised: the reference itself is not deterministic.
    if ("y" %in% all.vars(fk[[3L]])) next
    pk <- tryCatch(xtx:::design_plan(fk, tr), error = function(e) NULL, warning = function(w) NULL)
    if (is.null(pk)) next
    rows <- sample(nrow(tr), sample(1:12, 1), replace = TRUE)
    nd <- tr[rows, ]
    a <- capture(function() as.matrix(design_matrix(pk, nd)))
    b <- capture(function() reference(pk, nd))
    expect_true(same_outcome(a, b), info = sprintf("plan case %d: %s", k, deparse(fk)))
}

# printing
expect_stdout(print(p), "design_plan")
expect_stdout(print(design_matrix(p, train[1:3, ])), "design_matrix")

# Rebinding uses training levels and the training polynomial basis on every batch.
train <- data.frame(y = 1:8, x = 1:8, f = factor(rep(c("a", "b"), 4)), l = c(TRUE, FALSE))
fo <- y ~ poly(x, 2) + f + l
p <- xtx:::design_plan(fo, train)
original <- design_matrix(fo, train)
a <- data.frame(x = c(10, 12), f = "a", l = TRUE)
b <- data.frame(x = c(-3, -1), f = "b", l = FALSE)
for (D in list(original, design_matrix(p, train), unserialize(serialize(original, NULL)))) {
    before <- as.matrix(D)
    A <- design_matrix(D, a)
    B <- design_matrix(A, b)
    expect_identical(as.matrix(B), as.matrix(design_matrix(p, b)))
    expect_identical(as.matrix(design_matrix(D, b)), as.matrix(B))
    expect_identical(as.matrix(A), as.matrix(design_matrix(p, a)))
    expect_identical(as.matrix(D), before)
}
expect_error(design_matrix(original), "data")
expect_error(design_matrix(-original, a), "formula or a plan")
expect_error(design_matrix(original, transform(a, f = "new")), "new level")

# Retain the fitted transformation even when model.matrix subsets a supplied frame.
mf <- model.frame(fo, train)
subdesign <- design_matrix(~ poly(x, 2), mf)
expect_identical(as.matrix(design_matrix(subdesign, b)), as.matrix(design_matrix(xtx:::design_plan(~ poly(x, 2), train), b)))

fit <- lm(fo, train)
trained <- design_matrix(fit, newdata = NULL)
deferred <- design_matrix(fo, newdata = NULL)
expect_identical(as.matrix(design_matrix(fit, newdata = b)), as.matrix(design_matrix(trained, b)))
expect_identical(as.matrix(design_matrix(design_matrix(deferred, train), newdata = b)), as.matrix(design_matrix(trained, b)))
expect_identical(as.matrix(design_matrix(fit)), model.matrix(fit))
## (exported for package developers since jf7f)
expect_true(all(c("design_plan", "design_parts") %in% getNamespaceExports("xtx")))

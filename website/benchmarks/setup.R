library(design.matrix)
design_threads(1)
# Build the model frame exactly as model.matrix() does, so that comparisons
# with functions that call it are fair, and compute every row.
options(design.matrix.skip_na_copy = FALSE, design.matrix.dedup = FALSE)

set.seed(1)
n <- 2e5
big <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    f1 = factor(sample(sprintf("a%02d", 1:50), n, TRUE)),
    f2 = factor(sample(letters[1:20], n, TRUE)),
    cl = sample(1:1000, n, TRUE))
fo <- y ~ x * f1 + f2

# Prebuilt dense and design matrices, and a fitted model
X <- model.matrix(fo, big)
DX <- design_matrix(fo, big)
delayedAssign("fit_big", lm(fo, big))
delayedAssign("spec_big", design_matrix(fit_big, newdata = NULL))

# Checks for exact and approximate results
same <- function(a, b) identical(a, b, num.eq = FALSE, single.NA = FALSE, attrib.as.set = FALSE)
flat <- function(v) unname(c(unlist(v)))
close <- function(a, b) isTRUE(all.equal(flat(a), flat(b)))

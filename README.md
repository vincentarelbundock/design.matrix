

# xtx

The `xtx` package builds model matrices that are identical, bit for bit, to those returned by `stats::model.matrix()`. Its main benefit comes when the dense model matrix is only an intermediate input to something else: least squares coefficients need `X'X` and `X'y`, predictions multiply it by a vector of coefficients, their standard errors need a quadratic form, and cluster-robust variances need weighted sums of its rows by group. `xtx` computes these results directly from the factor codes and numeric variables, without materializing the dense matrix, which often saves an order of magnitude or more in both time and memory.

**Documentation, benchmarks and recipes: <https://vincentarelbundock.github.io/xtx>**

## Installation

``` r
# install.packages("remotes")
remotes::install_github("vincentarelbundock/xtx")
```

## Example

Build a design with an intercept, a numeric predictor, and a two-level factor. Multiplying it by a coefficient vector gives one linear prediction per row, and the same design applies to new data.

``` r
library(xtx)

d <- data.frame(x = c(1, 2, 3), group = factor(c("a", "b", "a")))
D <- design_matrix(~ x + group, d)
beta <- c("(Intercept)" = 10, x = 2, groupb = -1)
drop(D %*% beta)
#> [1] 12 13 16

newdata <- data.frame(x = c(4, 5), group = "a")
D_new <- design_matrix(D, newdata)
drop(D_new %*% beta)
#> [1] 18 20

identical(as.matrix(D), model.matrix(~ x + group, d))
#> [1] TRUE
```

## Promises

1.  **Exact.** Model matrices and those from `model.matrix()` pass `identical(a, b, num.eq = FALSE, single.NA = FALSE, attrib.as.set = FALSE)`, the strictest comparison R offers, for any number of threads.
2.  **New data.** A design built from a fitted model uses the same columns, factor coding, and transformations that the model uses for prediction.
3.  **Algebra.** Products and sums computed without the dense matrix agree with the dense computation up to floating point rounding, with `NA` and `NaN` in the same places.
4.  **Sparse.** Sparse model matrices hold the same non-zero values as the dense ones, bit for bit.

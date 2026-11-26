

# design.matrix

## Introduction

The `design.matrix` package builds model matrices that are identical, bit for bit, to those returned by `stats::model.matrix()`, and computes with them without forming the dense matrix. Most uses of a model matrix only need something derived from it: predictions, their standard errors, or weighted and grouped sums. Computing those directly from the factor codes and numeric variables saves both time and memory, often by an order of magnitude or more.

Use `design.matrix` when you have a formula (or fitted model) and data, and you want to avoid the cost of allocating, filling, and storing a dense model matrix. If you already have the dense matrix `X`, use it directly.

``` r
# install.packages("remotes")
remotes::install_github("vincentarelbundock/design.matrix")
```

How much time can you save when predicting on a new batch? The plot compares two numeric predictors with a numeric predictor interacting with a 20-level factor. Both methods reuse the fitted model’s terms, factor coding, and coefficients. Timing starts after one-time setup and includes multiplying by the coefficients. The package uses one thread.

<details class="code-fold">
<summary>Show benchmark code</summary>

``` r
library(design.matrix)
library(ggplot2)
design_threads(1)
set.seed(27)

n <- 100000
d <- data.frame(x = rnorm(n), z = rnorm(n), region = factor(sample(1:20, n, replace = TRUE)))
train <- data.frame(x = rnorm(5000), z = rnorm(5000), region = factor(sample(1:20, 5000, replace = TRUE)))
train$y <- rnorm(nrow(train))
time_model <- function(formula) {
    fit <- lm(formula, train)
    tt <- delete.response(terms(fit))
    spec <- design_matrix(fit, newdata = NULL)
    beta <- coef(fit)

    base_predict <- function(batch) {
        mf <- model.frame(tt, batch, xlev = fit$xlevels)
        X <- model.matrix(tt, mf, contrasts.arg = fit$contrasts)
        drop(X %*% beta)
    }

    fast_predict <- function(batch) {
        D <- design_matrix(spec, batch)
        drop(D %*% beta)
    }

    time_batch <- function(n) {
        batch <- head(d, n)
        timings <- bench::mark(
            base_predict(batch), fast_predict(batch),
            iterations = 15, memory = FALSE, filter_gc = FALSE,
            check = function(a, b) isTRUE(all.equal(as.numeric(a), as.numeric(b))))
        data.frame(rows = n, ms = 1000 * as.numeric(timings$median),
                   method = c("model.matrix()", "design_matrix()"))
    }

    do.call(rbind, lapply(c(1000, 10000, 25000, 50000, 100000), time_batch))
}

times <- rbind(
    transform(time_model(y ~ x + z), model = "Two numeric predictors"),
    transform(time_model(y ~ x * region), model = "Numeric × 20-level factor"))
p <- ggplot(times, aes(rows, ms, colour = method)) +
    geom_line(linewidth = 1) + geom_point(size = 2.5) +
    facet_wrap(~ model) +
    expand_limits(y = 0) +
    scale_x_continuous(labels = scales::label_comma()) +
    scale_colour_manual(values = c("model.matrix()" = "#ad481e", "design_matrix()" = "#18756d")) +
    labs(title = "Predicting on new batches of rows",
         x = "Rows per batch", y = "Median prediction time (ms)", colour = NULL) +
    theme_minimal(base_size = 14) + theme(legend.position = "bottom")
ggsave("prediction-batch-timing.png", p, width = 10, height = 5, dpi = 160, bg = "white")
```

</details>

![Median prediction time by batch size for base R and design_matrix().](prediction-batch-timing.png)

The plot measures only time, and only for predictions. The savings extend to memory and to other common tasks. On data with 200,000 rows and 119 columns, on one thread:

| Task | Base R | `design.matrix` |
|----|----|----|
| Predictions on new data | 123ms, 388MB | 9.5ms, 3.12MB |
| Standard errors of predictions | 361ms, 940MB | 17ms, 9.67MB |
| Grouped weighted sums (cluster-robust variance) | 121ms, 397MB | 28.3ms, 37.2MB |
| Sparse model matrix (vs. `Matrix::sparse.model.matrix()`) | 447ms, 662MB | 28.3ms, 40.7MB |

These gains come from never forming the dense model matrix. Each task above uses that matrix only as an intermediate input, for example to multiply it by a vector of coefficients, so `design.matrix` computes the result directly from the factor codes and numeric variables. When the dense matrix itself is what you need, there is no step to skip: `design.matrix` returns the identical matrix at about the same speed on one thread and about twice as fast on four threads, but it takes the same memory. The benchmarks below give the details.

### Promises

1.  **Exact.** Model matrices and those from `model.matrix()` pass `identical(a, b, num.eq = FALSE, single.NA = FALSE, attrib.as.set = FALSE)`, the strictest comparison R offers: same bits (signed zeros and `NA`/`NaN` payloads included), attributes in the same order, and the same warnings and errors, for any number of threads.
2.  **New data.** A design built from a fitted model uses the same columns, factor coding, and transformations that the model uses for prediction.
3.  **Algebra.** Products and sums computed without the dense matrix agree with the dense computation up to floating point rounding, with `NA` and `NaN` in the same places.
4.  **Sparse.** Sparse model matrices hold the same non-zero values as the dense ones, bit for bit.

The test suite checks these on hundreds of random formulas and data sets.

### Performance

Under the hood, a design matrix stores only the factor codes, the numeric variables and the contrast matrices, which can take much less memory than the dense matrix. Products, sums and quadratic forms are computed by C kernels that read each row’s few non-zero entries straight from these inputs, optionally on several threads, and rows with identical inputs can be computed once. When the dense matrix is needed, it is filled by a multi-threaded port of R’s own C routine.

Each section shows a small example, checked with `identical()` when the result is exact and with `all.equal()` when sums run in a different order, followed by a benchmark on medium-sized data (200,000 rows, 119 columns). The relative costs depend on the size and structure of the design.

The benchmark tables report minimum and median elapsed time and `mem_alloc`, the total bytes allocated on the R heap during one evaluation. Reproducing the memory column requires an R build with `capabilities("profmem")` enabled.

### Simple example

Build a design with an intercept, a numeric predictor, and a two-level factor. Multiplying it by a coefficient vector gives one linear prediction per row.

``` r
library(design.matrix)

d <- data.frame(x = c(1, 2, 3), group = factor(c("a", "b", "a")))
D <- design_matrix(~ x + group, d)
beta <- c("(Intercept)" = 10, x = 2, groupb = -1)
drop(D %*% beta)
#> [1] 12 13 16

newdata <- data.frame(x = c(4, 5), group = "a")
D_new <- design_matrix(D, newdata)
drop(D_new %*% beta)
#> [1] 18 20
```

## Operations

### `design_matrix()`

A design matrix holds the prepared variables and the layout of the columns, but not the dense matrix. `as.matrix()` builds the latter.

``` r
library(design.matrix)
d <- data.frame(
    y = c(0.3, -1.2, 0.8, 1.5, -0.4, 0.9, -2.1, 0.2),
    x = c(-1.5, 0.2, 1.1, -0.3, 2.0, 0.7, -0.9, 1.4),
    f = factor(c("a", "b", "c", "a", "b", "c", "a", "b")),
    cl = c(1, 1, 2, 2, 3, 3, 4, 4))

D <- design_matrix(y ~ x * f, d)
D
#> <design_matrix: 8 rows x 6 columns>
#> columns: (Intercept), x, fb, fc, x:fb, x:fc

a <- model.matrix(y ~ x * f, d)
b <- as.matrix(D)
identical(a, b, num.eq = FALSE)
#> [1] TRUE
```

The following construction benchmark uses one thread. Building the design defers the cost of allocating and filling the dense matrix until `as.matrix()` is called; the next section times that step.

``` r
library(design.matrix)
design_threads(1)
options(design.matrix.dedup = FALSE)
set.seed(1)
n <- 2e5
big <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    f1 = factor(sample(sprintf("a%02d", 1:50), n, TRUE)),
    f2 = factor(sample(letters[1:20], n, TRUE)),
    cl = sample(1:1000, n, TRUE))
fo <- y ~ x * f1 + f2

bench::mark(
    "stats::model.matrix" = model.matrix(fo, big),
    "design_matrix"       = design_matrix(fo, big),
    iterations = 10, check = FALSE)[c("expression", "min", "median", "mem_alloc")]
#> # A tibble: 2 × 4
#>   expression               min   median mem_alloc
#>   <bch:expr>          <bch:tm> <bch:tm> <bch:byt>
#> 1 stats::model.matrix  52.26ms   52.6ms   211.1MB
#> 2 design_matrix         9.97ms   10.2ms    29.6MB
```

### Using multiple threads with `as.matrix()`

Writing the dense output accounts for much of the work in both implementations. Additional threads can speed things up (see `design_threads()`). The memory column includes both the dense output and temporary allocations made while building it.

``` r
library(design.matrix)
options(design.matrix.dedup = FALSE)
set.seed(1)
n <- 2e5
big <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    f1 = factor(sample(sprintf("a%02d", 1:50), n, TRUE)),
    f2 = factor(sample(letters[1:20], n, TRUE)),
    cl = sample(1:1000, n, TRUE))
fo <- y ~ x * f1 + f2
same <- function(a, b) identical(a, b, num.eq = FALSE, single.NA = FALSE, attrib.as.set = FALSE)

fast_model_matrix <- function(formula, data, threads) {
    design_threads(threads)
    as.matrix(design_matrix(formula, data))
}

bench::mark(
    model.matrix(fo, big),
    fast_model_matrix(fo, big, threads = 1),
    fast_model_matrix(fo, big, threads = 4),
    iterations = 10, check = same)[c("expression", "min", "median", "mem_alloc")]
#> # A tibble: 3 × 4
#>   expression                                   min   median mem_alloc
#>   <bch:expr>                              <bch:tm> <bch:tm> <bch:byt>
#> 1 model.matrix(fo, big)                     52.9ms   52.9ms     211MB
#> 2 fast_model_matrix(fo, big, threads = 1)   51.3ms   51.5ms     211MB
#> 3 fast_model_matrix(fo, big, threads = 4)     23ms   23.3ms     211MB
design_threads(1)
```

### Reusing a design on new data

To make predictions, new rows must be turned into columns exactly as the training rows were. Here, `poly(x, 2)` must use the training polynomial, and `f` must keep one column per training level even though the new data only contain `"b"`. Pass the training design and the new data to `design_matrix()`, then use `as.matrix()` to display the columns.

``` r
library(design.matrix)
d <- data.frame(
    y = c(0.3, -1.2, 0.8, 1.5, -0.4, 0.9, -2.1, 0.2),
    x = c(-1.5, 0.2, 1.1, -0.3, 2.0, 0.7, -0.9, 1.4),
    f = factor(c("a", "b", "c", "a", "b", "c", "a", "b")),
    cl = c(1, 1, 2, 2, 3, 3, 4, 4))

D <- design_matrix(y ~ poly(x, 2) + f, d)
nd <- data.frame(x = c(0, 1), f = "b")
as.matrix(design_matrix(D, nd))
#>   (Intercept) poly(x, 2)1 poly(x, 2)2 fb fc
#> 1           1  -0.1070552  -0.3717871  1  0
#> 2           1   0.2101455  -0.2062446  1  0
#> attr(,"assign")
#> [1] 0 1 1 2 2
#> attr(,"contrasts")
#> attr(,"contrasts")$f
#> [1] "contr.treatment"
```

A fitted model works too, and gives the same result:

``` r
fit <- lm(y ~ poly(x, 2) + f, d)
identical(as.matrix(design_matrix(D, nd)), as.matrix(design_matrix(fit, newdata = nd)), num.eq = FALSE)
#> [1] TRUE
```

When you have many batches, prepare once with `spec <- design_matrix(fit, newdata = NULL)`, then call `design_matrix(spec, batch)` for each batch. The benchmark below does this for a 100-row batch.

``` r
library(design.matrix)
options(design.matrix.dedup = FALSE)
set.seed(1)
n <- 2e5
big <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    f1 = factor(sample(sprintf("a%02d", 1:50), n, TRUE)),
    f2 = factor(sample(letters[1:20], n, TRUE)),
    cl = sample(1:1000, n, TRUE))
fo <- y ~ x * f1 + f2
same <- function(a, b) identical(a, b, num.eq = FALSE, single.NA = FALSE, attrib.as.set = FALSE)

fit <- lm(fo, big)
spec <- design_matrix(fit, newdata = NULL)
tt <- delete.response(terms(fit))
nd <- big[1:100, ]
base_model_matrix <- function(batch) {
    mf <- model.frame(tt, batch, xlev = fit$xlevels, na.action = na.pass)
    model.matrix(tt, mf, contrasts.arg = fit$contrasts)
}
fast_model_matrix <- function(batch) {
    as.matrix(design_matrix(spec, batch))
}

bench::mark(
    base_model_matrix(nd),
    fast_model_matrix(nd),
    iterations = 200, check = same)[c("expression", "min", "median", "mem_alloc")]
#> # A tibble: 2 × 4
#>   expression                 min   median mem_alloc
#>   <bch:expr>            <bch:tm> <bch:tm> <bch:byt>
#> 1 base_model_matrix(nd)  294.2µs  336.7µs   224.1KB
#> 2 fast_model_matrix(nd)   26.2µs   27.3µs    94.4KB
```

### Arithmetic

Multiply or divide a design matrix by a number to scale it, or by a vector to scale each row. You can also add or subtract two designs with matching rows and columns. The result is another design matrix.

``` r
library(design.matrix)
D <- design_matrix(~ x, data.frame(x = c(2, 3)))
as.matrix(c(1, 2) * D)
#>   (Intercept) x
#> 1           1 2
#> 2           2 6
```

### `%*%`

`D %*% B` computes `X B` from the factor codes and numeric variables. `B` can have several columns.

``` r
library(design.matrix)
d <- data.frame(
    y = c(0.3, -1.2, 0.8, 1.5, -0.4, 0.9, -2.1, 0.2),
    x = c(-1.5, 0.2, 1.1, -0.3, 2.0, 0.7, -0.9, 1.4),
    f = factor(c("a", "b", "c", "a", "b", "c", "a", "b")),
    cl = c(1, 1, 2, 2, 3, 3, 4, 4))
D <- design_matrix(y ~ x * f, d)

M <- model.matrix(y ~ x * f, d)
beta <- c(1, 0.5, -1, 2, 0.1, 0.3)

a <- M %*% beta
b <- D %*% beta
all.equal(a, b, check.attributes = FALSE)
#> [1] TRUE
```

The benchmarks in this section and the following ones time each operation in two situations. Starting from the data, both sides build their matrix and then use it once, as when predicting on new data. Starting from a prebuilt matrix (`X` or `DX`, built in the setup), only the operation is timed, as when the same matrix is reused; the prebuilt design typically uses less memory than the dense matrix.

``` r
library(design.matrix)
options(design.matrix.dedup = FALSE)
set.seed(1)
n <- 2e5
big <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    f1 = factor(sample(sprintf("a%02d", 1:50), n, TRUE)),
    f2 = factor(sample(letters[1:20], n, TRUE)),
    cl = sample(1:1000, n, TRUE))
fo <- y ~ x * f1 + f2
X <- model.matrix(fo, big)
DX <- design_matrix(fo, big)
flat <- function(v) unname(c(unlist(v)))
close <- function(a, b) isTRUE(all.equal(flat(a), flat(b)))

beta <- rnorm(ncol(X))
bench::mark(
    "model.matrix, then %*%"  = model.matrix(fo, big) %*% beta,
    "design_matrix, then %*%" = design_matrix(fo, big) %*% beta,
    "prebuilt X"              = X %*% beta,
    "prebuilt design matrix"  = DX %*% beta,
    iterations = 10, check = close)[c("expression", "min", "median", "mem_alloc")]
#> # A tibble: 4 × 4
#>   expression                   min   median mem_alloc
#>   <bch:expr>              <bch:tm> <bch:tm> <bch:byt>
#> 1 model.matrix, then %*%   67.62ms  68.41ms  212.67MB
#> 2 design_matrix, then %*%  19.16ms  19.38ms    31.1MB
#> 3 prebuilt X               14.64ms  14.81ms    1.53MB
#> 4 prebuilt design matrix    8.83ms   8.98ms    1.54MB
```

### `crossprod()`

`crossprod(D, W)` computes `X'W`, such as `X'y` in least squares or weighted sums of the columns of `X`.

``` r
library(design.matrix)
d <- data.frame(
    y = c(0.3, -1.2, 0.8, 1.5, -0.4, 0.9, -2.1, 0.2),
    x = c(-1.5, 0.2, 1.1, -0.3, 2.0, 0.7, -0.9, 1.4),
    f = factor(c("a", "b", "c", "a", "b", "c", "a", "b")),
    cl = c(1, 1, 2, 2, 3, 3, 4, 4))
M <- model.matrix(y ~ x * f, d)
D <- design_matrix(y ~ x * f, d)

w <- seq(-1, 1, length.out = 8)

a <- crossprod(M, w)
b <- crossprod(D, w)
all.equal(a, b, check.attributes = FALSE)
#> [1] TRUE
```

``` r
library(design.matrix)
options(design.matrix.dedup = FALSE)
set.seed(1)
n <- 2e5
big <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    f1 = factor(sample(sprintf("a%02d", 1:50), n, TRUE)),
    f2 = factor(sample(letters[1:20], n, TRUE)),
    cl = sample(1:1000, n, TRUE))
fo <- y ~ x * f1 + f2
X <- model.matrix(fo, big)
DX <- design_matrix(fo, big)
flat <- function(v) unname(c(unlist(v)))
close <- function(a, b) isTRUE(all.equal(flat(a), flat(b)))

w <- rnorm(n)
bench::mark(
    "model.matrix, then crossprod"  = crossprod(model.matrix(fo, big), w),
    "design_matrix, then crossprod" = crossprod(design_matrix(fo, big), w),
    "prebuilt X"                    = crossprod(X, w),
    "prebuilt design matrix"        = crossprod(DX, w),
    iterations = 10, check = close)[c("expression", "min", "median", "mem_alloc")]
#> # A tibble: 4 × 4
#>   expression                         min   median mem_alloc
#>   <bch:expr>                    <bch:tm> <bch:tm> <bch:byt>
#> 1 model.matrix, then crossprod   67.72ms  68.31ms  211.14MB
#> 2 design_matrix, then crossprod  19.24ms  19.59ms    31.1MB
#> 3 prebuilt X                     14.41ms  14.53ms     1000B
#> 4 prebuilt design matrix          9.21ms   9.42ms    1.54MB
```

### `colSums()`, `colMeans()`, `rowSums()` and `rowMeans()`

The row and column summaries of base R work on design matrices. Pass one group label per row to `colSums()` or `colMeans()` to get one row per group, as with `base::rowsum()`. Group sums are a building block of cluster-robust standard errors and subgroup averages.

``` r
library(design.matrix)
d <- data.frame(site = c("B", "A", "B", "B", "A"), x = c(2, 10, 4, 6, 20))
D <- design_matrix(~ x, d)

colMeans(D)
#> (Intercept)           x 
#>         1.0         8.4
colSums(D, group = d$site)
#>   (Intercept)  x
#> A           2 30
#> B           3 12
colMeans(D, group = d$site)
#>   (Intercept)  x
#> A           1 15
#> B           1  4

all.equal(colSums(D, group = d$site), rowsum(model.matrix(~ x, d), d$site))
#> [1] TRUE
```

The benchmark computes weighted group sums. To weight the rows, multiply the design by a vector first.

``` r
library(design.matrix)
options(design.matrix.dedup = FALSE)
set.seed(1)
n <- 2e5
big <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    f1 = factor(sample(sprintf("a%02d", 1:50), n, TRUE)),
    f2 = factor(sample(letters[1:20], n, TRUE)),
    cl = sample(1:1000, n, TRUE))
fo <- y ~ x * f1 + f2
X <- model.matrix(fo, big)
DX <- design_matrix(fo, big)
flat <- function(v) unname(c(unlist(v)))
close <- function(a, b) isTRUE(all.equal(flat(a), flat(b)))

w <- rnorm(n)
bench::mark(
    "model.matrix, then rowsum"   = rowsum(w * model.matrix(fo, big), big$f2),
    "design_matrix, then colSums" = colSums(w * design_matrix(fo, big), group = big$f2),
    "prebuilt X"                  = rowsum(w * X, big$f2),
    "prebuilt design matrix"      = colSums(w * DX, group = big$f2),
    iterations = 10, check = close)[c("expression", "min", "median", "mem_alloc")]
#> # A tibble: 4 × 4
#>   expression                       min   median mem_alloc
#>   <bch:expr>                  <bch:tm> <bch:tm> <bch:byt>
#> 1 model.matrix, then rowsum      117ms  117.5ms     395MB
#> 2 design_matrix, then colSums     22ms   22.8ms      35MB
#> 3 prebuilt X                    63.4ms   65.3ms   183.9MB
#> 4 prebuilt design matrix        11.6ms   11.7ms     5.4MB
```

### `quadratic_diagonal()`

A common task is to extract the diagonal of `X V X'`. For example, in a linear regression context, with `V` the variance-covariance matrix of the coefficients and `X` the design matrix, the diagonal holds the variance of each fitted value, and its square root is the standard error. With `V = solve(crossprod(X))` instead, the diagonal holds the leverage of each observation. `quadratic_diagonal()` computes this diagonal without forming `X` or the full `n` by `n` product.

``` r
library(design.matrix)
d <- data.frame(
    y = c(0.3, -1.2, 0.8, 1.5, -0.4, 0.9, -2.1, 0.2),
    x = c(-1.5, 0.2, 1.1, -0.3, 2.0, 0.7, -0.9, 1.4),
    f = factor(c("a", "b", "c", "a", "b", "c", "a", "b")),
    cl = c(1, 1, 2, 2, 3, 3, 4, 4))
M <- model.matrix(y ~ x * f, d)
D <- design_matrix(y ~ x * f, d)

V <- diag(6) / 10

a <- rowSums((M %*% V) * M)
b <- quadratic_diagonal(D, V)
all.equal(a, b, check.attributes = FALSE)
#> [1] TRUE
```

Base R builds an `n` by `p` intermediate; the kernel touches only the non-zero entries of each row.

``` r
library(design.matrix)
options(design.matrix.dedup = FALSE)
set.seed(1)
n <- 2e5
big <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    f1 = factor(sample(sprintf("a%02d", 1:50), n, TRUE)),
    f2 = factor(sample(letters[1:20], n, TRUE)),
    cl = sample(1:1000, n, TRUE))
fo <- y ~ x * f1 + f2
X <- model.matrix(fo, big)
DX <- design_matrix(fo, big)
flat <- function(v) unname(c(unlist(v)))
close <- function(a, b) isTRUE(all.equal(flat(a), flat(b)))

V <- crossprod(matrix(rnorm(ncol(X)^2), ncol(X))) / ncol(X)
base_diagonal <- function(X, V) rowSums((X %*% V) * X)

bench::mark(
    base_diagonal(model.matrix(fo, big), V),
    quadratic_diagonal(design_matrix(fo, big), V),
    base_diagonal(X, V),
    quadratic_diagonal(DX, V),
    iterations = 5, check = close)[c("expression", "min", "median", "mem_alloc")]
#> # A tibble: 4 × 4
#>   expression                                         min   median mem_alloc
#>   <bch:expr>                                    <bch:tm> <bch:tm> <bch:byt>
#> 1 base_diagonal(model.matrix(fo, big), V)        182.2ms  189.1ms   575.8MB
#> 2 quadratic_diagonal(design_matrix(fo, big), V)     21ms   21.2ms    31.3MB
#> 3 base_diagonal(X, V)                            135.2ms  136.7ms   364.7MB
#> 4 quadratic_diagonal(DX, V)                       10.5ms   10.6ms     1.7MB
```

### Sparse model matrix

`as(D, "CsparseMatrix")` builds a sparse matrix without forming the dense one. It is identical to converting the dense model matrix.

``` r
library(design.matrix)
library(Matrix)
d <- data.frame(
    y = c(0.3, -1.2, 0.8, 1.5, -0.4, 0.9, -2.1, 0.2),
    x = c(-1.5, 0.2, 1.1, -0.3, 2.0, 0.7, -0.9, 1.4),
    f = factor(c("a", "b", "c", "a", "b", "c", "a", "b")),
    cl = c(1, 1, 2, 2, 3, 3, 4, 4))
M <- model.matrix(y ~ x * f, d)
D <- design_matrix(y ~ x * f, d)

dense_to_sparse <- function(M) as(as(as(M, "dMatrix"), "generalMatrix"), "CsparseMatrix")
a <- dense_to_sparse(M)
b <- as(D, "CsparseMatrix")
identical(a, b)
#> [1] TRUE
```

`Matrix::sparse.model.matrix()` and `MatrixModels::model.Matrix()` are shown for reference: they have their own rules for contrasts and column names. This benchmark uses a factor interaction and checks that all four methods produce the same dimensions and entries.

``` r
library(design.matrix)
library(Matrix)
library(MatrixModels)
options(design.matrix.dedup = FALSE)
set.seed(1)
n <- 2e5
big <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    f1 = factor(sample(sprintf("a%02d", 1:50), n, TRUE)),
    f2 = factor(sample(letters[1:20], n, TRUE)),
    cl = sample(1:1000, n, TRUE))
base_sparse <- function(formula, data) {
    X <- model.matrix(formula, data)
    as(as(as(X, "dMatrix"), "generalMatrix"), "CsparseMatrix")
}
fast_sparse <- function(formula, data) {
    D <- design_matrix(formula, data)
    as(D, "CsparseMatrix")
}

sfo <- y ~ x + f1 * f2
reference <- fast_sparse(sfo, big)
same_entries <- function(x) {
    identical(dim(x), dim(reference)) && Matrix::nnzero(x - reference) == 0
}
stopifnot(
    identical(base_sparse(sfo, big), reference),
    same_entries(Matrix::sparse.model.matrix(sfo, big)),
    same_entries(MatrixModels::model.Matrix(sfo, big, sparse = TRUE)))
bench::mark(
    base_sparse(sfo, big),
    fast_sparse(sfo, big),
    Matrix::sparse.model.matrix(sfo, big),
    MatrixModels::model.Matrix(sfo, big, sparse = TRUE),
    iterations = 5, check = FALSE)[c("expression", "min", "median", "mem_alloc")]
#> # A tibble: 4 × 4
#>   expression                                               min   median mem_alloc
#>   <bch:expr>                                          <bch:tm> <bch:tm> <bch:byt>
#> 1 base_sparse(sfo, big)                                  1.14s     1.3s    3.02GB
#> 2 fast_sparse(sfo, big)                                28.04ms   28.4ms   40.75MB
#> 3 Matrix::sparse.model.matrix(sfo, big)               379.19ms  445.4ms  661.52MB
#> 4 MatrixModels::model.Matrix(sfo, big, sparse = TRUE)  350.4ms    360ms  661.53MB
```

### Repeated rows

When rows have identical predictor values, `design.matrix` can compute the result once per distinct row and reuse it for every matching row. Row-wise results stay identical bit for bit. The option `design.matrix.dedup` controls this: `"auto"` (default), `TRUE` (always) or `FALSE` (never).

Repeated rows are found once per design matrix and remembered, so the first call pays for the search. The first two timings include design construction. The cached timing reuses both the design and its groups. Here, 5,000 patterns repeat over 200,000 rows.

``` r
library(design.matrix)
options(design.matrix.dedup = FALSE)
set.seed(1)
n <- 2e5
dup <- data.frame(
    f1 = factor(sample.int(50, n, TRUE)),
    f3 = factor(sample.int(5, n, TRUE)),
    f2 = factor(sample.int(20, n, TRUE)))
pd <- design_matrix(~ f1 * f3 + f2, dup, contrasts.arg = list(f1 = "contr.sum", f3 = "contr.sum"))
Vd <- diag(ncol(design_matrix(pd, dup[1:2, ])))
same <- function(a, b) identical(a, b, num.eq = FALSE, single.NA = FALSE, attrib.as.set = FALSE)

quad <- function(mode, D = design_matrix(pd, dup)) {
    op <- options(design.matrix.dedup = mode)
    on.exit(options(op))
    quadratic_diagonal(D, Vd)
}
D_cached <- design_matrix(pd, dup)
invisible(quad("auto", D_cached))
bench::mark(
    "every row"              = quad(FALSE),
    "grouped, first call"    = quad("auto"),
    "grouped, groups cached" = quad("auto", D_cached),
    iterations = 5, check = same)[c("expression", "min", "median", "mem_alloc")]
#> # A tibble: 3 × 4
#>   expression                  min   median mem_alloc
#>   <bch:expr>             <bch:tm> <bch:tm> <bch:byt>
#> 1 every row               37.78ms  38.67ms    4.67MB
#> 2 grouped, first call      3.34ms   3.41ms    5.78MB
#> 3 grouped, groups cached   1.44ms   1.46ms    2.41MB
```

## Recipes

### Predictions

The predictions of a linear model are `X b`, the design matrix times the coefficients.

``` r
library(design.matrix)
d <- data.frame(
    y = c(0.3, -1.2, 0.8, 1.5, -0.4, 0.9, -2.1, 0.2),
    x = c(-1.5, 0.2, 1.1, -0.3, 2.0, 0.7, -0.9, 1.4),
    f = factor(c("a", "b", "c", "a", "b", "c", "a", "b")),
    cl = c(1, 1, 2, 2, 3, 3, 4, 4))

fast_predict <- function(fit, newdata) {
    D <- design_matrix(fit, newdata = newdata)
    drop(D %*% coef(fit))
}

fit <- lm(y ~ x + f, d)
a <- predict(fit, d)
b <- fast_predict(fit, d)
all.equal(a, b, check.attributes = FALSE)
#> [1] TRUE
```

``` r
library(design.matrix)
options(design.matrix.dedup = FALSE)
set.seed(1)
n <- 2e5
big <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    f1 = factor(sample(sprintf("a%02d", 1:50), n, TRUE)),
    f2 = factor(sample(letters[1:20], n, TRUE)),
    cl = sample(1:1000, n, TRUE))
fo <- y ~ x * f1 + f2
close <- function(a, b) isTRUE(all.equal(unname(a), unname(b)))

fast_predict <- function(fit, newdata, spec) {
    D <- design_matrix(spec, newdata)
    drop(D %*% coef(fit))
}

fit_big <- lm(fo, big)
spec_big <- design_matrix(fit_big, newdata = NULL)
bench::mark(
    "predict(fit, newdata)" = predict(fit_big, big),
    "fast_predict"          = fast_predict(fit_big, big, spec_big),
    iterations = 5, check = close)[c("expression", "min", "median", "mem_alloc")]
#> # A tibble: 2 × 4
#>   expression                 min   median mem_alloc
#>   <bch:expr>            <bch:tm> <bch:tm> <bch:byt>
#> 1 predict(fit, newdata)  56.98ms     87ms   387.7MB
#> 2 fast_predict            9.07ms    9.4ms    3.07MB
```

### Standard errors of predictions

The standard errors of the predictions are the square roots of the diagonal of `X V X'`, where `V` is the variance-covariance matrix of the coefficients.

``` r
library(design.matrix)
d <- data.frame(
    y = c(0.3, -1.2, 0.8, 1.5, -0.4, 0.9, -2.1, 0.2),
    x = c(-1.5, 0.2, 1.1, -0.3, 2.0, 0.7, -0.9, 1.4),
    f = factor(c("a", "b", "c", "a", "b", "c", "a", "b")),
    cl = c(1, 1, 2, 2, 3, 3, 4, 4))

base_predict_se <- function(fit, newdata) {
    predict(fit, newdata, se.fit = TRUE)$se.fit
}
fast_predict_se <- function(fit, newdata) {
    D <- design_matrix(fit, newdata = newdata)
    sqrt(quadratic_diagonal(D, vcov(fit)))
}

fit <- lm(y ~ x + f, d)
a <- base_predict_se(fit, d)
b <- fast_predict_se(fit, d)
all.equal(a, b, check.attributes = FALSE)
#> [1] TRUE
```

For this full-rank linear model, `predict(..., se.fit = TRUE)` returns predictions and standard errors. The helper computes only the standard errors and reuses a prepared design specification.

``` r
library(design.matrix)
options(design.matrix.dedup = FALSE)
set.seed(1)
n <- 2e5
big <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    f1 = factor(sample(sprintf("a%02d", 1:50), n, TRUE)),
    f2 = factor(sample(letters[1:20], n, TRUE)),
    cl = sample(1:1000, n, TRUE))
fo <- y ~ x * f1 + f2
close <- function(a, b) isTRUE(all.equal(unname(a), unname(b)))

base_predict_se <- function(fit, newdata) {
    predict(fit, newdata, se.fit = TRUE)$se.fit
}
fast_predict_se <- function(fit, newdata, spec) {
    D <- design_matrix(spec, newdata)
    sqrt(quadratic_diagonal(D, vcov(fit)))
}

fit_big <- lm(fo, big)
spec_big <- design_matrix(fit_big, newdata = NULL)
bench::mark(
    base_predict_se(fit_big, big),
    fast_predict_se(fit_big, big, spec_big),
    iterations = 5, check = close)[c("expression", "min", "median", "mem_alloc")]
#> # A tibble: 2 × 4
#>   expression                                   min   median mem_alloc
#>   <bch:expr>                              <bch:tm> <bch:tm> <bch:byt>
#> 1 base_predict_se(fit_big, big)            257.9ms    261ms   940.4MB
#> 2 fast_predict_se(fit_big, big, spec_big)   16.6ms     17ms    11.2MB
```

### Leverage

The hat values of a linear model are the diagonal of `X (X'X)^-1 X'`, and `(X'X)^-1` is `vcov(fit) / sigma(fit)^2`.

``` r
library(design.matrix)
d <- data.frame(
    y = c(0.3, -1.2, 0.8, 1.5, -0.4, 0.9, -2.1, 0.2),
    x = c(-1.5, 0.2, 1.1, -0.3, 2.0, 0.7, -0.9, 1.4),
    f = factor(c("a", "b", "c", "a", "b", "c", "a", "b")),
    cl = c(1, 1, 2, 2, 3, 3, 4, 4))

fit <- lm(y ~ x + f, d)
D <- design_matrix(fit, newdata = d)
a <- hatvalues(fit)
b <- quadratic_diagonal(D, vcov(fit) / sigma(fit)^2)
all.equal(a, b, check.attributes = FALSE)
#> [1] TRUE
```

`hatvalues()` works from the QR decomposition stored in the fit. Both manual calculations start with a prebuilt design and the same precomputed inverse crossproduct.

``` r
library(design.matrix)
options(design.matrix.dedup = FALSE)
set.seed(1)
n <- 2e5
big <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    f1 = factor(sample(sprintf("a%02d", 1:50), n, TRUE)),
    f2 = factor(sample(letters[1:20], n, TRUE)),
    cl = sample(1:1000, n, TRUE))
fo <- y ~ x * f1 + f2
fit_big <- lm(fo, big)
X <- model.matrix(fo, big)
D <- design_matrix(fo, big)
flat <- function(v) unname(c(unlist(v)))
close <- function(a, b) isTRUE(all.equal(flat(a), flat(b)))

XtXi <- vcov(fit_big) / sigma(fit_big)^2
base_leverage <- function(X, XtXi) {
    rowSums((X %*% XtXi) * X)
}
fast_leverage <- function(D, XtXi) {
    quadratic_diagonal(D, XtXi)
}

bench::mark(
    hatvalues(fit_big),
    base_leverage(X, XtXi),
    fast_leverage(D, XtXi),
    iterations = 5, check = close)[c("expression", "min", "median", "mem_alloc")]
#> # A tibble: 3 × 4
#>   expression                  min   median mem_alloc
#>   <bch:expr>             <bch:tm> <bch:tm> <bch:byt>
#> 1 hatvalues(fit_big)        2.22s    2.67s    15.3MB
#> 2 base_leverage(X, XtXi) 119.39ms 123.69ms   364.7MB
#> 3 fast_leverage(D, XtXi)  10.89ms  11.16ms     1.7MB
```

### Cluster-robust standard errors

The cluster-robust covariance matrix is the sandwich `B M B`, with bread `B = (X'X)^-1` and meat `M = S'S`, where row `c` of `S` sums `e_i x_i` over the observations in cluster `c` (without a small-sample correction). The same function works on the dense model matrix and on the design matrix.

``` r
library(design.matrix)
d <- data.frame(
    y = c(0.3, -1.2, 0.8, 1.5, -0.4, 0.9, -2.1, 0.2),
    x = c(-1.5, 0.2, 1.1, -0.3, 2.0, 0.7, -0.9, 1.4),
    f = factor(c("a", "b", "c", "a", "b", "c", "a", "b")),
    cl = c(1, 1, 2, 2, 3, 3, 4, 4))

base_cluster_vcov <- function(X, e, cluster, bread) {
    S <- rowsum(e * X, cluster)
    bread %*% crossprod(S) %*% bread
}
fast_cluster_vcov <- function(D, e, cluster, bread) {
    S <- colSums(e * D, group = cluster)
    bread %*% crossprod(S) %*% bread
}
fit <- lm(y ~ x + f, d)
D <- design_matrix(fit, newdata = d)
bread <- vcov(fit) / sigma(fit)^2
a <- base_cluster_vcov(model.matrix(fit), residuals(fit), d$cl, bread)
b <- fast_cluster_vcov(D, residuals(fit), d$cl, bread)
all.equal(a, b, check.attributes = FALSE)
#> [1] TRUE
```

``` r
library(design.matrix)
options(design.matrix.dedup = FALSE)
set.seed(1)
n <- 2e5
big <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    f1 = factor(sample(sprintf("a%02d", 1:50), n, TRUE)),
    f2 = factor(sample(letters[1:20], n, TRUE)),
    cl = sample(1:1000, n, TRUE))
fo <- y ~ x * f1 + f2
fit_big <- lm(fo, big)
X <- model.matrix(fo, big)
D <- design_matrix(fo, big)
XtXi <- vcov(fit_big) / sigma(fit_big)^2
base_cluster_vcov <- function(X, e, cluster, bread) {
    S <- rowsum(e * X, cluster)
    bread %*% crossprod(S) %*% bread
}
fast_cluster_vcov <- function(D, e, cluster, bread) {
    S <- colSums(e * D, group = cluster)
    bread %*% crossprod(S) %*% bread
}
flat <- function(v) unname(c(unlist(v)))
close <- function(a, b) isTRUE(all.equal(flat(a), flat(b)))

e <- residuals(fit_big)
bench::mark(
    base_cluster_vcov(X, e, big$cl, XtXi),
    fast_cluster_vcov(D, e, big$cl, XtXi),
    iterations = 5, check = close)[c("expression", "min", "median", "mem_alloc")]
#> # A tibble: 2 × 4
#>   expression                                 min   median mem_alloc
#>   <bch:expr>                            <bch:tm> <bch:tm> <bch:byt>
#> 1 base_cluster_vcov(X, e, big$cl, XtXi)   37.5ms   65.6ms     186MB
#> 2 fast_cluster_vcov(D, e, big$cl, XtXi)   12.6ms   12.7ms      11MB
```

## Setup

``` r
cat(R.version.string, "\n")
#> R version 4.5.3 (2026-03-11)
cat("cores:", parallel::detectCores(), "\n")
#> cores: 32
cat("BLAS:", basename(extSoftVersion()[["BLAS"]]), "\n")
#> BLAS: libblas.so.3
```

#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: [Operations])
#metadata((title: "Operations")) <website-metadata>

#import "/benchmarks/bench.typ": bench, bench-table, bench-code

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
options(width = 200)
```
]

#title()

Each section shows a small example, checked with `identical()` when the result is exact and with `all.equal()` when sums run in a different order, followed by a benchmark on medium-sized data (200,000 rows, 119 columns). The relative costs depend on the size and structure of the design.

= Data

The small examples use this data frame:

```r
library(xtx)
d <- data.frame(
    y = c(0.3, -1.2, 0.8, 1.5, -0.4, 0.9, -2.1, 0.2),
    x = c(-1.5, 0.2, 1.1, -0.3, 2.0, 0.7, -0.9, 1.4),
    f = factor(c("a", "b", "c", "a", "b", "c", "a", "b")),
    cl = c(1, 1, 2, 2, 3, 3, 4, 4))
```

The benchmarks use this one, with one thread unless stated otherwise. `X` and `DX` are prebuilt dense and design matrices, and `same()` and `close()` are the checks for exact and approximate results. Each benchmark runs in a fresh R process on top of this setup; the results shown were saved by `make bench` (see #link("https://github.com/vincentarelbundock/xtx/tree/main/website/benchmarks")[`website/benchmarks`]).

#bench-code("setup")

The benchmark tables report minimum and median elapsed time and `mem_alloc`, the total bytes allocated on the R heap during one evaluation. Reproducing the memory column requires an R build with `capabilities("profmem")` enabled. For fair comparisons with functions that call `model.matrix()`, the benchmarks set `options(xtx.skip_na_copy = FALSE)`, so that the model frame is built in exactly the same way (by default, `xtx` skips a copy of the frame that `na.omit()` makes even when no row has a missing value).

= Building designs

== `design_matrix()`

A design matrix holds the prepared variables and the layout of the columns, but not the dense matrix. `as.matrix()` builds the latter.

```r
D <- design_matrix(y ~ x * f, d)
D

a <- model.matrix(y ~ x * f, d)
b <- as.matrix(D)
identical(a, b, num.eq = FALSE)
```

Building the design defers the cost of allocating and filling the dense matrix until `as.matrix()` is called; the next section times that step.

#bench("construction")

== Using multiple threads with `as.matrix()`

Writing the dense output accounts for much of the work in both implementations. Additional threads can speed things up (see #link("reference/design_threads.html")[`design_threads()`]). The memory column includes both the dense output and temporary allocations made while building it.

#bench("threads")

== Reusing a design on new data

To make predictions, new rows must be turned into columns exactly as the training rows were. Here, `poly(x, 2)` must use the training polynomial, and `f` must keep one column per training level even though the new data only contain `"b"`. Pass the training design and the new data to `design_matrix()`, then use `as.matrix()` to display the columns.

```r
D <- design_matrix(y ~ poly(x, 2) + f, d)
nd <- data.frame(x = c(0, 1), f = "b")
as.matrix(design_matrix(D, nd))
```

A fitted model works too, and gives the same result:

```r
fit <- lm(y ~ poly(x, 2) + f, d)
identical(as.matrix(design_matrix(D, nd)), as.matrix(design_matrix(fit, newdata = nd)), num.eq = FALSE)
```

When you have many batches, prepare once with `spec <- design_matrix(fit, newdata = NULL)`, then call `design_matrix(spec, batch)` for each batch. The benchmark below does this for a 100-row batch.

#bench("reuse")

= Algebra

== Arithmetic

Multiply or divide a design matrix by a number to scale it, or by a vector to scale each row. You can also add or subtract two designs with matching rows and columns. The result is another design matrix.

```r
D <- design_matrix(~ x, data.frame(x = c(2, 3)))
as.matrix(c(1, 2) * D)
```

== `%*%`

`D %*% B` computes `X B` from the factor codes and numeric variables. `B` can have several columns.

```r
M <- model.matrix(y ~ x * f, d)
D <- design_matrix(y ~ x * f, d)
beta <- c(1, 0.5, -1, 2, 0.1, 0.3)

a <- M %*% beta
b <- D %*% beta
all.equal(a, b, check.attributes = FALSE)
```

The benchmarks in this section time each operation in two situations. Starting from the data, both sides build their matrix and then use it once, as when predicting on new data. Starting from a prebuilt matrix (`X` or `DX`), only the operation is timed, as when the same matrix is reused; the prebuilt design typically uses less memory than the dense matrix.

#bench("matmul")

== `crossprod()`

`crossprod(D, W)` computes `X'W`, such as `X'y` in least squares or weighted sums of the columns of `X`. `crossprod(D)` computes `X'X`, and `crossprod(D, w * D)` the weighted `X'WX`; both return ordinary matrices, which go straight to `chol()`, `chol2inv()`, `eigen()` or `solve()`.

```r
w <- seq(-1, 1, length.out = 8)

a <- crossprod(M, w)
b <- crossprod(D, w)
all.equal(a, b, check.attributes = FALSE)

all.equal(crossprod(M), crossprod(D))
all.equal(crossprod(M, w * M), crossprod(D, w * D))
```

#bench("crossprod_w")

For `X'X`, the design matrix adds up the products of each row's few non-zero entries, instead of all `p^2 / 2` products of the dense row.

#bench("crossprod")

== `qr()`

`qr(D)` decomposes a design matrix without forming `X` or `Q`, and the usual functions accept the result: `qr.coef()`, `qr.fitted()`, `qr.resid()`, `qr.qty()`, `qr.qy()`, `qr.Q()`, `qr.R()` and `solve()`. `R` is the Cholesky factor of `crossprod(D)`, with columns pivoted as in `qr()` and `lm()`, and `Q` is `X R^-1`, never formed. Coefficients are refined from the residuals `y - X b`, which recovers much of the accuracy that the normal equations lose on badly conditioned designs. Columns whose residual norm is close to the tolerance are checked against the data, so aliased columns are those of `qr()`.

Results agree with `qr(as.matrix(D))` up to floating point rounding, not bit for bit. Signs of the rows of `R` (and of the columns of `Q`) can differ, and `qr.Q()` and `qr.qty()` return only the first `rank` columns or rows, since the complete `Q` is never formed. See #link("reference/design_qr.html")[`design_qr`]. The #link("recipes/estimation.html")[estimation recipes] benchmark it against other least squares implementations.

```r
a <- qr(M)
b <- qr(D)
b
all.equal(qr.coef(a, d$y), qr.coef(b, d$y))
all.equal(qr.resid(a, d$y), qr.resid(b, d$y), check.attributes = FALSE)
all.equal(abs(qr.R(a)), abs(qr.R(b)), check.attributes = FALSE)
```

== `colSums()`, `colMeans()`, `rowSums()` and `rowMeans()`

The row and column summaries of base R work on design matrices. Pass one group label per row to `colSums()` or `colMeans()` to get one row per group, as with `base::rowsum()`. Group sums are a building block of cluster-robust standard errors and subgroup averages.

```r
sites <- data.frame(site = c("B", "A", "B", "B", "A"), x = c(2, 10, 4, 6, 20))
Ds <- design_matrix(~ x, sites)

colMeans(Ds)
colSums(Ds, group = sites$site)
colMeans(Ds, group = sites$site)

all.equal(colSums(Ds, group = sites$site), rowsum(model.matrix(~ x, sites), sites$site))
```

The benchmark computes weighted group sums. To weight the rows, multiply the design by a vector first.

#bench("colsums")

== `quadratic_diagonal()`

A common task is to extract the diagonal of `X V X'`. For example, in a linear regression context, with `V` the variance-covariance matrix of the coefficients and `X` the design matrix, the diagonal holds the variance of each fitted value, and its square root is the standard error. With `V = solve(crossprod(X))` instead, the diagonal holds the leverage of each observation. `quadratic_diagonal()` computes this diagonal without forming `X` or the full `n` by `n` product.

```r
V <- diag(6) / 10

a <- rowSums((M %*% V) * M)
b <- quadratic_diagonal(D, V)
all.equal(a, b, check.attributes = FALSE)
```

Base R builds an `n` by `p` intermediate; the kernel touches only the non-zero entries of each row.

#bench("quadratic")

= Sparse matrices

`as(D, "CsparseMatrix")` builds a sparse matrix without forming the dense one. It is identical to converting the dense model matrix.

```r
library(Matrix)
dense_to_sparse <- function(M) as(as(as(M, "dMatrix"), "generalMatrix"), "CsparseMatrix")
a <- dense_to_sparse(M)
b <- as(D, "CsparseMatrix")
identical(a, b)
```

`Matrix::sparse.model.matrix()` and `MatrixModels::model.Matrix()` are shown for reference: they have their own rules for contrasts and column names. This benchmark uses a factor interaction and checks that all four methods produce the same dimensions and entries.

#bench("sparse")

= Repeated rows

When rows have identical predictor values, `xtx` can compute the result once per distinct row and reuse it for every matching row. Row-wise results stay identical bit for bit. The option `xtx.dedup` controls this: `"auto"` (default), `TRUE` (always) or `FALSE` (never).

Repeated rows are found once per design matrix and remembered, so the first call pays for the search. The first two timings include design construction. The cached timing reuses both the design and its groups. Here, 5,000 patterns repeat over 200,000 rows.

#bench("dedup")

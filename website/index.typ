#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: [design.matrix])
#metadata((title: "Introduction")) <website-metadata>

#import "/benchmarks/bench.typ": bench-table

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)


#calepin.chunk(echo: false, results: "hide")[
```r
options(width = 200)
```
]

#title()

The `design.matrix` package builds model matrices that are identical, bit for bit, to those returned by `stats::model.matrix()`. Its main benefit comes when the dense model matrix is only an intermediate input to something else: least squares coefficients need `X'X` and `X'y`, predictions multiply it by a vector of coefficients, their standard errors need a quadratic form, and cluster-robust variances need weighted sums of its rows by group. `design.matrix` computes these results directly from the factor codes and numeric variables, without materializing the dense matrix, which often saves an order of magnitude or more in both time and memory. When the dense matrix itself is what you need, there is no step to skip, and the gains are much smaller.

Use `design.matrix` when you have a formula (or fitted model) and data, and you need the model matrix only on the way to something else, so that allocating, filling, and storing it would be wasted work. If you already have the dense matrix `X`, use it directly.

= Installation

```r
#| eval: false
# install.packages("remotes")
remotes::install_github("vincentarelbundock/design.matrix")
```

= Performance

How much time can you save when predicting on a new batch? The plot compares two numeric predictors with a numeric predictor interacting with a 20-level factor. Both methods reuse the fitted model's terms, factor coding, and coefficients. Timing starts after one-time setup and includes multiplying by the coefficients. The package uses one thread.

#calepin.chunk(
  echo: false,
  results: "render",
  fig-width: "100%",
  fig-device-width: 10,
  fig-device-height: 5,
  fig-alt-text: "Median prediction time by batch size for base R and design_matrix().",
)[
```r
library(ggplot2)
times <- read.csv("benchmarks/results/batch_timing.csv")
ggplot(times, aes(rows, ms, colour = method)) +
    geom_line(linewidth = 1) + geom_point(size = 2.5) +
    facet_wrap(~ model) +
    expand_limits(y = 0) +
    scale_x_continuous(labels = scales::label_comma()) +
    scale_colour_manual(values = c("model.matrix()" = "#ad481e", "design_matrix()" = "#18756d")) +
    labs(title = "Predicting on new batches of rows",
         x = "Rows per batch", y = "Median prediction time (ms)", colour = NULL) +
    theme_minimal(base_size = 14) + theme(legend.position = "bottom")
```
]

The plot measures only time, and only for predictions. The savings extend to memory and to other common tasks. On data with 200,000 rows and 119 columns, on one thread:

#bench-table("summary", code-column: false)

Each of these tasks uses the model matrix only as an intermediate input, so `design.matrix` never forms it. When you do need the dense matrix, `design.matrix` returns the identical matrix at about the same speed on one thread and about twice as fast on four threads, but it takes the same memory. The #link("operations.html")[Operations] page gives the details.

= Promises

+ *Exact.* Model matrices and those from `model.matrix()` pass `identical(a, b, num.eq = FALSE, single.NA = FALSE, attrib.as.set = FALSE)`, the strictest comparison R offers: same bits (signed zeros and `NA`/`NaN` payloads included), attributes in the same order, and the same warnings and errors, for any number of threads.
+ *New data.* A design built from a fitted model uses the same columns, factor coding, and transformations that the model uses for prediction.
+ *Algebra.* Products and sums computed without the dense matrix agree with the dense computation up to floating point rounding, with `NA` and `NaN` in the same places.
+ *Sparse.* Sparse model matrices hold the same non-zero values as the dense ones, bit for bit.

The test suite checks these on hundreds of random formulas and data sets.

= How it works

Under the hood, a design matrix stores only the factor codes, the numeric variables and the contrast matrices, which can take much less memory than the dense matrix. Products, sums and quadratic forms are computed by C kernels that read each row's few non-zero entries straight from these inputs, optionally on several threads, and rows with identical inputs can be computed once. When the dense matrix is needed, it is filled by a multi-threaded port of R's own C routine.

= Simple example

Build a design with an intercept, a numeric predictor, and a two-level factor. Multiplying it by a coefficient vector gives one linear prediction per row.

```r
library(design.matrix)

d <- data.frame(x = c(1, 2, 3), group = factor(c("a", "b", "a")))
D <- design_matrix(~ x + group, d)
beta <- c("(Intercept)" = 10, x = 2, groupb = -1)
drop(D %*% beta)

newdata <- data.frame(x = c(4, 5), group = "a")
D_new <- design_matrix(D, newdata)
drop(D_new %*% beta)
```

= Data in memory or in a database

The same code fits data frames, DuckDB tables and Parquet files: give a connection and a table instead of a data frame, and the products, `design_lm()`, `design_glm()` and fixed effects run in the database, with the same results as in memory.

```r
#| eval: false
con <- DBI::dbConnect(duckdb::duckdb())
DBI::dbWriteTable(con, "d", d)
design_lm(y ~ x + f, d)                        # a data frame
design_lm(y ~ x + f, con, table = "d")         # the same fit, in DuckDB
design_glm(k ~ x + f, poisson, con, table = "d")
design_lm(y ~ x | firm, con, table = "d", vcov = ~firm)
```

The #link("recipes/out-of-core.html")[Out-of-core data] recipe runs every call on both backends side by side (#link("recipes/out-of-core.html#one-code-two-backends")[One code, two backends]) and describes the remote designs.

= Setup

The benchmarks on this site ran on:

#raw(read("/benchmarks/results/machine.txt"), block: true)

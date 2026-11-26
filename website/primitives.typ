#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: [Primitives for package developers])
#metadata((title: "Primitives for package developers")) <website-metadata>

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
options(width = 200)
```
]

#title()

An estimator written with a few exported functions runs unchanged on a data frame and on a database table. The building blocks are the design and its products (`design_matrix()`, `%*%`, `crossprod()`, `colSums(group = )`, `quadratic_diagonal()`), *columns*, which hold one value per row, and base R arithmetic on columns. In memory a column is a numeric vector. For a remote design it is a `design_column`: a lazy SQL expression over the rows, which nothing evaluates until a sum, a product or `as.vector()` asks for its values. Only small results (numbers, $p times p$ matrices, one row per cluster) come back to R.

Remote examples need the DBI and duckdb packages.

#table(
  columns: 2,
  [*Task*], [*Functions*],
  [Designs and products], [`design_matrix()`, `D %*% b`, `crossprod(D)`, `crossprod(D, col)`, `crossprod(D, w * D)`, `colSums(e * D, group = )`, `quadratic_diagonal(D, V)`],
  [Columns], [`design_column(D, "y")`, arithmetic, comparisons, `exp()`, `log()`, `round()`, `lgamma()`, `is.na()`, ...],
  [Combining columns], [`col_ifelse()`, `col_pmin()`, `col_pmax()`],
  [Summaries], [`sum()`, `mean()`, `min()`, `max()`, `range()`, `any()`, `all()`, `col_sums()` (several sums, one query)],
  [Values], [`as.vector()`, `as.matrix()`, `as.data.frame(keep = )`, `col_compute()` (a temporary table, or a permanent one with `name =`)],
  [Groupings], [`design_groups()`, `group_sums()`, `group_crosstab()`, `group_lookup()`, `group_nested()`],
  [Compression], [`design_patterns()`],
)

= Data

The same data, as a data frame and as a DuckDB table: an outcome `y` with a few gross outliers, a numeric and a categorical regressor, and 40 clusters.

```r
library(design.matrix)
set.seed(48)
n <- 5000
d <- data.frame(
    id = seq_len(n),
    x = rnorm(n),
    f = sample(c("a", "b", "c"), n, TRUE),
    cl = sample(1:40, n, TRUE))
d$y <- 1 + d$x + (d$f == "b") + rnorm(40)[d$cl] + rnorm(n)
out <- sample(n, 50)
d$y[out] <- d$y[out] + 20

con <- DBI::dbConnect(duckdb::duckdb())
DBI::dbWriteTable(con, "d", d)

D_mem <- design_matrix(y ~ x + f, d)
D_db  <- design_matrix(y ~ x + f, con, table = "d")
```

= Least squares with clustered standard errors

The coefficients come from $X'X$ and $X'y$; the residuals are a column; the scores of each cluster are the column sums of the design scaled by the residuals, $sum_(i in g) e_i x_i$, one row per cluster. The sandwich $V = G/(G-1) (X'X)^(-1) S'S (X'X)^(-1)$ is then a small computation in R.

```r
ols_cluster <- function(D, cl) {
    y <- design_column(D, "y")
    XtX <- crossprod(D)
    b <- drop(solve(XtX, crossprod(D, y)))
    e <- y - drop(D %*% b)                  # a column
    S <- colSums(e * D, group = cl)         # G x p cluster scores
    G <- nrow(S)
    bread <- solve(XtX)
    V <- G / (G - 1) * bread %*% crossprod(S) %*% bread
    fit <- cbind(estimate = b, std.error = sqrt(diag(V)))
    attr(fit, "rss") <- sum(e^2)
    fit
}
```

The only difference between the two calls is how the clusters are named: in memory a vector, remotely a column of the table (a design in memory keeps only the variables of its formula).

```r
fit_mem <- ols_cluster(D_mem, d$cl)
fit_db  <- ols_cluster(D_db, "cl")
fit_db
all.equal(fit_mem, fit_db)
```

The same numbers in base R, from the dense model matrix:

```r
X <- model.matrix(y ~ x + f, d)
b <- solve(crossprod(X), crossprod(X, d$y))
e <- drop(d$y - X %*% b)
S <- rowsum(e * X, d$cl)
bread <- solve(crossprod(X))
V <- 40 / 39 * bread %*% crossprod(S) %*% bread
all.equal(unname(fit_db[, "std.error"]), unname(sqrt(diag(V))))
```

= A one-step robust estimator

One step of iteratively reweighted least squares for Huber's M-estimator, from the least squares fit: the residuals are scaled, the rows with large residuals are downweighted with `col_pmin()`, and a weighted least squares step gives the new coefficients. `col_compute()` stores the residuals in a temporary table, since they are used several times (in memory it returns the vector). `col_sums()` reads several sums in one query.

```r
huber_step <- function(D, b, k = 1.345) {
    y <- design_column(D, "y")
    e <- col_compute(y - drop(D %*% b))
    s <- sqrt(sum(e^2) / (nrow(D) - ncol(D)))
    u <- abs(e) / s
    w <- col_pmin(1, k / u)                 # Huber weights
    b1 <- drop(solve(crossprod(D, w * D), crossprod(D, w * y)))
    list(coef = b1,
         scale = s,
         range = range(e),
         sums = col_sums(downweighted = u > k, weight = w,
                         big = col_ifelse(abs(e) > 10, 1, 0)))
}

b0 <- fit_mem[, "estimate"]
step_mem <- huber_step(D_mem, b0)
step_db  <- huber_step(D_db, b0)
step_db
all.equal(step_mem, step_db)
```

The least squares slope is pulled by the outliers; one Huber step brings the coefficients back near those of the data without them:

```r
rbind(ols = b0, huber = step_db$coef,
      clean = coef(lm(y ~ x + f, d[-out, ])))
```

= Groupings

Fixed-effect estimators work from *groupings*: each row gets the number $1..G$ of its group, from the values of one or more key variables. Five functions cover them:

- `design_groups(D, key)` makes the grouping of the rows of `D` by `key`;
- `group_sums(col, g)`: the sums of a column (or several) within each group, a $G times k$ matrix, one query remotely; `group_sums(1, g)` counts the rows of each group;
- `group_crosstab(g1, g2, w)`: the (weighted) number of rows in each pair of groups, a sparse $G_1 times G_2$ matrix;
- `group_lookup(values, g)`: the column whose value in a row is the value of its group, `values[codes]`;
- `group_nested(g1, g2)`: whether each group of `g1` lies within one group of `g2`.

A grouping is a list with `codes` (the group of each row: an integer vector in memory, a `design_column` remotely), `G` (the number of groups), `levels` (a data frame of the key values of each group, in code order, sorted as the levels of fixed effects), `n` (the number of rows) and `design`. Results with one row per group are in code order: row `l` is `g$levels[l, ]`. Every row must have a group: a missing key (`NA`, `NaN`, SQL `NULL`) is an error on both backends, so select the rows with a key first (remotely `D[DBI::SQL("cl IS NOT NULL"), ]`).

The within (one-way fixed effect) slope of `y` on `x` demeans both columns by cluster, $tilde(y) = y - group_lookup(group_sums(y, g) / n_g, g)$, and regresses one on the other. The key is given as in the clustered example: in memory by its values, remotely by the name of a column.

```r
within_slope <- function(D, g) {
    n_g <- group_sums(1, g)[, 1]              # rows per cluster
    demean <- function(v) v - group_lookup(group_sums(v, g)[, 1] / n_g, g)
    yt <- demean(design_column(D, "y"))
    xt <- demean(design_column(D, "x"))
    col <- col_sums(xy = xt * yt, xx = xt^2)  # one query remotely
    col[["xy"]] / col[["xx"]]
}

g_mem <- design_groups(D_mem, d["cl"])
g_db  <- design_groups(D_db, "cl")
g_db$G
head(g_db$levels)
c(memory = within_slope(D_mem, g_mem), duckdb = within_slope(D_db, g_db),
  lm = coef(lm(y ~ x + factor(cl), d))[["x"]])
all.equal(within_slope(D_mem, g_mem), within_slope(D_db, g_db))
```

= Columns are lazy

Printing a remote column shows its dimensions ($n times k$, as `dim()`) and its SQL, without reading a row. Arithmetic builds a new expression, and the values are computed only when a summary, a product or `as.vector()` needs them:

```r
y <- design_column(D_db, "y")
r <- y - drop(D_db %*% b0)
r
z <- col_ifelse(abs(r) > 10, NA, r)
c(sum = sum(z), sum_na_rm = sum(z, na.rm = TRUE), mean = mean(z, na.rm = TRUE))
head(as.data.frame(round(r, 2), keep = "id", nm = "resid"))
```

Missing values follow R: SQL `NULL` is `NA` and `NaN` stays `NaN`, `is.na()` is true for both, and `sum()`, `mean()`, `min()` and the other summaries give `NA` where base R does, with `na.rm = TRUE` dropping both. Columns of one design combine with each other, also when they come from separate products; columns of two designs do not. A column is valid while its connection is open.

```r
DBI::dbDisconnect(con)
```

#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: [Out-of-core data])
#metadata((title: "Out-of-core data")) <website-metadata>

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
options(width = 200)
```
]

#title()

When the data do not fit in memory, `design_matrix()` and `design_lm()` take a DuckDB connection and a table (or a `sql` query), or the path of Parquet files, instead of a data frame. The rows stay in the database: products with the design run there as `GROUP BY` aggregates, or in batches, and only small results come back to R. These designs need the DBI and duckdb packages.

```r
library(design.matrix)
con <- DBI::dbConnect(duckdb::duckdb())
set.seed(1)
n <- 10000
d <- data.frame(
    x = rnorm(n),
    z = sample(1:5, n, TRUE),
    f = sample(c("a", "b", "c"), n, TRUE),
    w = runif(n, 0.5, 2),
    cl = sample(1:50, n, TRUE))
d$y <- d$x + d$z + rnorm(n) * (1 + abs(d$x))
DBI::dbWriteTable(con, "d", d)
```

= One code, two backends

A design built from a data frame and one built from a DuckDB table answer the same operations: products, sums, `qr()`, `design_lm()`, `predict()`. Results with one number per row come back as vectors in memory and as tables in the database; everything else is the same R object.

```r
fit_mem <- design_lm(y ~ x + f, d)
fit_db  <- design_lm(y ~ x + f, con, table = "d")
all.equal(coef(fit_mem), coef(fit_db))
all.equal(vcov(fit_mem), vcov(fit_db))

D_mem <- design_matrix(y ~ x + f, d)
D_db  <- design_matrix(y ~ x + f, con, table = "d")
all.equal(crossprod(D_mem, d$w * D_mem), crossprod(D_db, "w" * D_db))
all.equal(qr.coef(qr(D_mem), d$y), qr.coef(qr(D_db), "y"))
head(predict(fit_mem, newdata = d))
head(as.data.frame(predict(fit_db, newdata = con, table = "d")))
```

= Remote designs

`design_matrix(formula, con, table = )` returns a remote design. Factor levels are pinned in R, with R's sorting rules rather than the database collation, and rows with a missing value in a model variable are dropped, as `na.omit()` drops them.

```r
D <- design_matrix(y ~ x * f + log(z), con, table = "d")
D
```

The same works with `sql = "SELECT ..."`, or with Parquet files, which a private DuckDB connection reads:

```r
pq <- tempfile(fileext = ".parquet")
DBI::dbExecute(con, sprintf("COPY d TO '%s' (FORMAT PARQUET)", pq))
design_matrix(y ~ x * f + log(z), pq)
```

The SQL engine handles factors, character and logical columns, numeric columns, interactions, and arithmetic with `log()`, `exp()`, `sqrt()` and `abs()`. Any other formula, such as one with `poly()` or `splines::ns()`, uses the chunked engine: batches of rows go through `design_matrix(spec, batch)`, and the results are accumulated. `summary()` shows which terms the SQL engine computes, and `summary(D, sql = TRUE)` the queries.

```r
D2 <- design_matrix(y ~ splines::ns(x, 3) + f, con, table = "d")
summary(D2)
```

To fit in R and score in the database, pass a design or a fitted model instead of the formula: its levels, contrasts and transformations are used, and rows with missing values are kept, with `NULL` predictions, as in memory.

= Algebra

`crossprod()`, `colSums()`, `colMeans()`, `rowSums()`, `quadratic_diagonal()`, `qr()` and its methods, `as.matrix()` and sparse coercion (for small data), row selection with `D[DBI::SQL("x > 0"), ]` and row scaling work as in memory. Rows are scaled by a column name or a SQL expression (`DBI::SQL()`) instead of a vector, so `crossprod(D, "w" * D)` is `X'WX`:

```r
X <- model.matrix(y ~ x * f + log(z), d)
all.equal(crossprod(D), crossprod(X), check.attributes = FALSE)
all.equal(crossprod(D, "w" * D), crossprod(X, d$w * X), check.attributes = FALSE)
all.equal(crossprod("w" * D, "y"), crossprod(X, d$w * d$y), check.attributes = FALSE)
all.equal(colSums(D, group = "cl"), rowsum(X, d$cl), check.attributes = FALSE)
```

When the rows have few distinct design patterns (factors and discrete numeric variables), these products and `design_lm()` run one `GROUP BY` over the input columns of the design, and compute the result from the design of the patterns, in memory. The option `design.matrix.remote_compress` controls this.

= Predictions

`D %*% beta` writes the predictions, next to the columns of the remote data, to a temporary view, and returns its name. `as.data.frame()` reads it; a table keeps it.

```r
beta <- solve(crossprod(D), crossprod(D, "y"))
p <- D %*% beta
head(as.data.frame(p))
DBI::dbExecute(con, sprintf("CREATE TABLE preds AS SELECT * FROM %s", p))
```

= Linear models

`design_lm()` fits a linear model on remote data with the same algorithm as in memory: with few design patterns, on the patterns with their counts and sums of the response; otherwise, from `qr()` of the remote design, with `X'X` computed in the database. Either way, the coefficients, aliased columns, covariance and residual standard error are those of `lm()` on the same rows, and `weights` names a column.

```r
fit <- design_lm(y ~ x + f, con, table = "d", weights = "w")
all.equal(coef(fit), coef(lm(y ~ x + f, d, weights = w)))
all.equal(vcov(fit), vcov(lm(y ~ x + f, d, weights = w)))

fit_z <- design_lm(y ~ factor(z) * f, con, table = "d")
fit_z
```

`weights`, `offset` and `subset` are columns, `DBI::SQL()` expressions or R expressions with a SQL translation, `offset()` terms work, and `na.action` can be `na.omit`, `na.exclude`, `na.fail` or `na.pass`, as for `lm()`. Results with one value per row are tables in the database: `fitted(fit)` (with `.fitted` and `.resid`), and `predict(fit, newdata = con, table = )`, with `se.fit`, confidence and prediction intervals, and the non-estimable cases of a rank-deficient fit, as `predict.lm()`.

```r
fit_s <- design_lm(y ~ x + f, con, table = "d", subset = z > 2, offset = log(w))
p <- predict(fit_s, newdata = con, table = "d", interval = "confidence")
head(as.data.frame(p))
```

Expressions with no SQL translation, such as `weights = myfun(x)` or a response `my_transform(y)`, are evaluated in R, batch by batch, into a temporary table on which the fit runs; with `engine = "auto"`, a message names the expression.

= Generalized linear models

`design_glm()` takes remote data as `design_lm()` does, with `weights` and `offset` naming columns. With few design patterns, one `GROUP BY` returns the sums of each pattern (including the working weights of the first iteration and the parts of the deviance and AIC that do not depend on the coefficients), and the iterations run in memory. Otherwise, each iteration computes the linear predictor, the working weights and responses, `X'WX` and `X'Wz` in the database (the family and link are translated to SQL), or in batches with the chunked engine. Either way, the results are those of `glm()` on the same rows.

```r
d2 <- transform(d, k = rpois(n, exp(0.2 * z)), yb = as.integer(y > 5))
DBI::dbWriteTable(con, "d2", d2)
fit <- design_glm(k ~ x + f, poisson, con, table = "d2")
all.equal(coef(fit), coef(glm(k ~ x + f, poisson, d2)))
fit$method

fit_z <- design_glm(yb ~ factor(z) * f, binomial, con, table = "d2")
all.equal(coef(fit_z), coef(glm(yb ~ factor(z) * f, binomial, d2)))
fit_z$method
```

`predict(fit, newdata = con, table = , type = "response")` scores in the database, and `residuals(fit, type = )` returns the remote design with `.fitted` and `.resid` columns, as for linear models.

= Robust standard errors

`residuals()` of a remote fit returns the remote design with two more columns, `.fitted` and `.resid`, which can scale its rows. The heteroskedasticity-consistent (HC1) and clustered (CR1, as in `sandwich::vcovCL()`) covariance matrices follow:

```r
fit <- design_lm(y ~ x + f, con, table = "d")
R <- residuals(fit)
bread <- vcov(fit) / sigma(fit)^2
n <- nobs(fit)
k <- fit$rank

meat <- crossprod(R, DBI::SQL('".resid" * ".resid"') * R)
V_hc1 <- n / (n - k) * bread %*% meat %*% bread

scores <- colSums(".resid" * R, group = "cl")
G <- nrow(scores)
V_cl <- G / (G - 1) * (n - 1) / (n - k) * bread %*% crossprod(scores) %*% bread
sqrt(cbind(HC1 = diag(V_hc1), CR1 = diag(V_cl)))
```

The same computation on the dense model matrix:

```r
l <- lm(y ~ x + f, d)
X <- model.matrix(l)
e <- residuals(l)
B <- solve(crossprod(X))
all.equal(V_hc1, n / (n - k) * B %*% crossprod(X, e^2 * X) %*% B, check.attributes = FALSE)
S <- rowsum(e * X, d$cl)
all.equal(V_cl, G / (G - 1) * (n - 1) / (n - k) * B %*% crossprod(S) %*% B, check.attributes = FALSE)
```

= Differences from in-memory fits

- `model.frame()` is not available for a remote fit, and `model` is ignored: the data stay in the database (`fit$remote_design` is the remote design of the fit).
- Results with one value per row are tables in the database, not vectors: `fitted()` (`.fitted`, `.resid`), `predict()` (`.fitted`; `se.fit` as a column `.se.fit` with attributes `df` and `residual.scale`; `.nonestim` instead of the attribute `"non-estim"`), and `residuals()`, a remote design.
- These tables have no row order: rows are matched to the data through its columns, which every table keeps.
- `compress = "auto"` decides from other evidence than in memory; the results are the same either way.

= Limits

Remote designs keep only the algebra promise: results agree with the in-memory computation up to rounding. `as.matrix()` refuses to collect more than 100,000 rows unless `force = TRUE`. Only DuckDB connections are supported. With infinite values in the data, `NaN` and `Inf` can land in other places than in the dense computation. See `?design_remote` for the details.

#calepin.chunk(echo: false, results: "hide")[
```r
DBI::dbDisconnect(con)
```
]

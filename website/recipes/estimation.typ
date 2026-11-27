#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: [Estimation])
#metadata((title: "Estimation")) <website-metadata>

#import "/benchmarks/bench.typ": bench, bench-table, bench-code

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
options(width = 200)
```
]

#title()

These recipes fit a linear model and compute the quantities that follow from its design matrix: coefficients, leverage, and cluster-robust standard errors. Each one checks a small example against base R, then times both on 200,000 rows and 119 columns, on one thread.

```r
library(xtx)
library(xtx.fit)
d <- data.frame(
    y = c(0.3, -1.2, 0.8, 1.5, -0.4, 0.9, -2.1, 0.2),
    x = c(-1.5, 0.2, 1.1, -0.3, 2.0, 0.7, -0.9, 1.4),
    f = factor(c("a", "b", "c", "a", "b", "c", "a", "b")),
    cl = c(1, 1, 2, 2, 3, 3, 4, 4))
```

The benchmarks run on top of the #link("../operations.html#data")[benchmark setup], in a fresh R process each.

= Fitting linear models

`lm_xtx()`, from the package `xtx.fit`, fits a linear model like `lm()`, without the dense model matrix: the coefficients, covariance matrix, residual standard error and degrees of freedom agree with those of `lm()`, aliased columns included. `summary()`, `predict()`, `residuals()`, `fitted()`, `vcov()` and `logLik()` work as for `lm()`.

```r
fit <- lm_xtx(y ~ x + f, d)
all.equal(coef(fit), coef(lm(y ~ x + f, d)))
all.equal(vcov(fit), vcov(lm(y ~ x + f, d)))
summary(fit)
```

When many rows share a design pattern (factors and discrete numeric variables), the fit runs on the patterns, with their counts and the sums of the response, and gives the same results as the fit on every row. With a design matrix instead of a formula, give the response as `y`; with remote data (a DuckDB table or Parquet files), see #link("out-of-core.html")[out-of-core data].

```r
D <- design_matrix(~ x + f, d)
all.equal(coef(lm_xtx(D, y = d$y)), coef(fit))
```

= Fixed effects

A formula `y ~ x | fe1 + fe2` absorbs fixed effects, as `fixest::feols()` does: each fixed effect is a variable, or an interaction `a^b`. The fit never builds the dummies and never demeans the rows. It solves the normal equations from sums over rows: products with the design, sums within the levels of each fixed effect, and crosstabs of the fixed effects. The largest fixed effect is eliminated exactly, and the others are solved densely (up to 1,000 levels) or by conjugate gradients. When the reduced system is sparse (low mobility, such as workers who rarely change firms), it is formed and preconditioned by an incomplete Cholesky factor; when it would be near dense (high mobility), it is applied without being formed, with a Jacobi preconditioner. `options(xtx.fe_precond = "jacobi")` or `"ichol"` forces a preconditioner, and `xtx.fe_form_max` (`0` never forms the system, `Inf` always) moves the threshold. Singletons are dropped as `feols()` drops them, and `vcov = ~cl` clusters the standard errors with `feols()`'s small sample corrections. The same code fits remote data (a DuckDB table or Parquet files; see #link("out-of-core.html")[out-of-core data]).

```r
set.seed(3)
panel <- data.frame(id = rep(1:200, each = 5), t = rep(1:5, 200),
                    x = rnorm(1000), f = factor(sample(c("a", "b"), 1000, TRUE)))
panel$y <- panel$x + (panel$f == "b") + panel$id / 100 + panel$t / 10 + rnorm(1000)
fe <- lm_xtx(y ~ x + f | id + t, panel, vcov = ~id)
summary(fe)
if (requireNamespace("fixest", quietly = TRUE)) {
    fx <- fixest::feols(y ~ x + f | id + t, panel, cluster = ~id, fixef.tol = 1e-11, notes = FALSE)
    print(all.equal(coef(fe), coef(fx)))
    print(all.equal(unclass(vcov(fe)), unclass(vcov(fx)), check.attributes = FALSE))
}
```

`fe$fixef` (or `fixef(fe)` with fixest or nlme loaded) holds the values of the fixed effects, one named vector per fixed effect, and `predict(fe, newdata)` adds them to `X b`: levels absent from the fit give `NA`. Offsets (`offset =` or `offset()` terms) are subtracted from the response, as in `feols()`. `subset`, `se.fit` and intervals in `predict()`, and clustering by more than one variable are not supported with fixed effects.

```r
nd <- data.frame(id = c(1, 2, 999), t = c(1, 5, 2), x = 0, f = factor("a", levels = c("a", "b")))
predict(fe, nd)
```

= Generalized linear models

`glm_xtx()` fits a generalized linear model like `glm()`, by iteratively reweighted least squares on the design matrix: each iteration needs `X'WX`, `X'Wz` and `X beta`, which the design computes without the dense matrix. The iterations are those of `glm.fit()`, step for step (starting values, step halving, convergence rule, warnings), so the coefficients, covariance matrix, deviances, AIC and number of iterations agree with those of `glm()`.

```r
set.seed(2)
dg <- data.frame(x = rnorm(500), f = factor(sample(c("a", "b", "c"), 500, TRUE)),
                 e = sample(1:3, 500, TRUE))
dg$y <- rbinom(500, 1, plogis(dg$x))
dg$k <- rpois(500, dg$e * exp(0.3 * as.integer(dg$f)))
fit <- glm_xtx(y ~ x + f, binomial, dg)
all.equal(coef(fit), coef(glm(y ~ x + f, binomial, dg)))
all.equal(vcov(fit), vcov(glm(y ~ x + f, binomial, dg)))
summary(fit)
```

When rows share a design pattern and an offset, they share their fitted value, and the fit runs on the patterns: for every supported family, a pattern enters each iteration only through the sum of its prior weights and its weighted mean response, and the deviance of the rows is that of the patterns plus a constant. The first iteration, which starts from per-row values, is summed from the rows, so even the iterations are those of the fit on every row.

```r
fit <- glm_xtx(k ~ f + offset(log(e)), poisson, dg)
fit$patterns
all.equal(deviance(fit), deviance(glm(k ~ f + offset(log(e)), poisson, dg)))
```

Fixed effects work in `glm_xtx()` as in `lm_xtx()`, with the results of `fixest::feglm()` and `fixest::fepois()`. The iterations are those of `glm.fit()`, with the fixed effects among the parameters: the weighted least squares problem of each iteration is solved by the fixed-effect estimator of linear models, from sums over rows. As fixest does by default, rows are dropped (iteratively) when their level of some fixed effect is a singleton or has an all-zero count (Poisson) or all-zero or all-one outcome (binomial), since its fixed effect would be infinite; `fit$nperfect` and `fit$nsingletons` count them. Offsets, weights and `vcov = ~cl` work as above. Both fits below are run to a tight convergence tolerance, so that they agree to many digits (the covariance matrix comes from the weights of the last iteration, which is one step behind the coefficients).

```r
panel$k <- rpois(1000, exp(0.3 * panel$x + panel$id / 200))
panel$k[panel$id %% 25 == 0] <- 0
fp <- glm_xtx(k ~ x + f | id + t, poisson, panel, vcov = ~id, control = list(epsilon = 1e-10))
summary(fp)
if (requireNamespace("fixest", quietly = TRUE)) {
    fx <- fixest::fepois(k ~ x + f | id + t, panel, cluster = ~id, glm.tol = 1e-10, fixef.tol = 1e-11, notes = FALSE)
    print(c(nobs(fp), nobs(fx)))
    print(all.equal(coef(fp), coef(fx)))
    print(all.equal(unclass(vcov(fp)), unclass(vcov(fx)), check.attributes = FALSE))
}
```

The families are those of `stats` (binomial, Poisson, Gaussian, Gamma, inverse Gaussian and the quasi families) with their standard links. Every option of `glm_xtx()` also works on remote data (see #link("out-of-core.html")[out-of-core data]), which is why `subset`, `etastart`, `mustart` and standard errors of predictions are not offered.

= Least squares

Least squares coefficients only need `X'X` and `X'y`, so `qr()` on a design matrix gives them without forming the dense matrix. With weights, use `qr(sqrt(w) * D)` and `sqrt(w) * y`.

```r
fit <- lm(y ~ x + f, d)
D <- design_matrix(y ~ x + f, d)
all.equal(coef(fit), qr.coef(qr(D), d$y))
```

#calepin.chunk(echo: false, results: "typst")[
```r
## (the number of rows, from the results: the run sets it, OLS_BENCH_N)
ols_n <- unique(read.csv(file.path(normalizePath(file.path(Sys.getenv("XTX_ROOT", "../.."), "bench", "ols")),
                                   "results.csv"))$n)[1]
cat("The benchmark below compares `qr()` on a design matrix with eight other implementations of least squares, on ",
    format(ols_n, big.mark = ",", scientific = FALSE),
    " rows and designs from 3 to 401 columns. It runs outside this page, because a fair comparison needs a fresh R process for each method and design: `bench/ols/run.sh` runs it and writes `bench/ols/results.csv`, which the tables below read. The protocol:\n", sep = "")
```
]

- Two comparisons. _From formula_: every method starts from the formula and the data, as a user would. _Prebuilt_: the model matrix (or design matrix) is built before timing, and only the fit is timed.
- One thread everywhere, including the BLAS. The model frame is built exactly as `model.matrix()` builds it (option `xtx.skip_na_copy = FALSE`).
- Memory is the peak resident memory of the first call in a fresh process, so it counts every allocation, including those made in C++.
- Every result is compared with the coefficients of `lm()`, including which columns are aliased (`NA`). A method that gives a different answer is marked "no".

#calepin.chunk(echo: false, results: "hide")[
```r
bench_dir <- normalizePath(file.path(Sys.getenv("XTX_ROOT", "../.."), "bench", "ols"))
res <- read.csv(file.path(bench_dir, "results.csv"))
res$p <- as.integer(res$p)
ols_table <- function(mode) {
    x <- res[res$mode == mode, ]
    x$cell <- ifelse(is.na(x$ms), "error",
                     sprintf("%s ms, %s MB%s", format(x$ms, big.mark = ",", trim = TRUE),
                             x$MB, ifelse(x$matches_lm %in% TRUE, "", " (no)")))
    tab <- reshape(x[, c("method", "p", "cell")], idvar = "method", timevar = "p", direction = "wide")
    names(tab) <- c("Method", paste(sort(unique(x$p)), "columns"))
    tab <- tab[order(tab$Method != "xtx qr()", tab$Method), ]
    tinytable::tt(tab)
}
```
]

*From formula*

#calepin.chunk(echo: false, results: "typst")[
```r
ols_table("formula")
```
]

*Prebuilt*

#calepin.chunk(echo: false, results: "typst")[
```r
ols_table("prebuilt")
```
]

#calepin.chunk(echo: false, results: "verbatim")[
```r
cat(readLines(file.path(bench_dir, "machine.txt")), sep = "\n")
```
]

From 21 columns up, `qr()` on a design matrix is the fastest method that matches `lm()` in both comparisons. From formula, with 21 to 119 columns, it is 2 to 11 times faster than the next method, and it uses the least memory, because it never builds the dense matrix. Prebuilt, it is 2 to 7 times faster, though some methods need even less memory, since they read the prebuilt matrix in place. At 401 columns, 21 of them aliased, the only other method that matches `lm()` is `.lm.fit()` (the same routine that `lm()` uses), 45 times slower from formula and 140 times slower prebuilt; the faster alternatives stop with an error, drop the aliased columns, or return a different solution.

With 3 to 5 numeric columns, the dense matrix is small and the gains disappear. From formula, `qr()` ties for fastest with 3 columns and is about 4 ms slower than `RcppEigen`'s Cholesky decomposition with 5; `fixest`, which does not build a model frame, uses the least memory. Prebuilt, `RcppEigen` is faster by 1 to 3 ms.

= Leverage

The hat values of a linear model are the diagonal of `X (X'X)^-1 X'`, and `(X'X)^-1` is `vcov(fit) / sigma(fit)^2`.

```r
fit <- lm(y ~ x + f, d)
D <- design_matrix(fit, newdata = d)
a <- hatvalues(fit)
b <- quadratic_diagonal(D, vcov(fit) / sigma(fit)^2)
all.equal(a, b, check.attributes = FALSE)
```

`hatvalues()` works from the QR decomposition stored in the fit. Both manual calculations start with a prebuilt design and the same precomputed inverse crossproduct.

#bench("leverage")

= Cluster-robust standard errors

The cluster-robust covariance matrix is the sandwich `B M B`, with bread `B = (X'X)^-1` and meat `M = S'S`, where row `c` of `S` sums `e_i x_i` over the observations in cluster `c` (without a small-sample correction). The same function works on the dense model matrix and on the design matrix.

```r
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
```

#bench("cluster")

## The fit-specific side of remote data: the documentation of remote fits,
## the arguments of remote data a fit takes in `...`, and the offsets of new
## remote data for predictions. The input stage itself (design_problem(),
## materialization of R expressions) is in problem.R.

#' @name design_remote
#' @section Linear models:
#' [design_lm()] works on remote data as on a data frame:
#' `design_lm(y ~ x + f, con, table = "flights", weights = w, subset = year > 2000)`,
#' or on a remote design: `design_lm(D, y = "col", weights = "w")`.
#' `weights`, `offset`, `subset` and `y` are column names (strings or bare
#' names), `DBI::SQL()` expressions, or R expressions over the columns, as
#' are the response and `offset()` terms of the formula. R expressions with
#' a SQL translation (arithmetic as for the SQL engine; for `subset`, also
#' comparisons, `&`, `|`, `!`, `%in%` with constants and `is.na()`) run in
#' the database. Any other R expression, such as `weights = myfun(x)`, is
#' evaluated in R, batch by batch, and its values are written with the
#' columns of the remote data to a temporary table, on which the fit runs
#' (one more pass over the data, with a message naming the expression;
#' `engine = "sql"` refuses). Compression stays exact, since the computed values are columns
#' of that table. `na.action` can be [na.omit()] (default), [na.exclude()],
#' [na.fail()] or [na.pass()], with the results or errors of [lm()] (with
#' fixed effects, `na.omit` or `na.exclude`; other na.action functions are
#' refused, in memory too). The
#' rows used are those with no missing model variable, response, weight or
#' offset, and factor levels are those of these rows, as in memory.
#'
#' The fit is the same algorithm as in memory, on remote primitives: with
#' few design patterns it is compressed (see *Compression*); otherwise it
#' runs [qr()][design_qr] of the remote design (`X'WX` in the database,
#' columns judged with [lm()]'s rule, suspect columns from products with
#' the data), with the factor's values computed with numeric variables
#' shifted when that is an exact transform of the design (so that `x + 1e7`
#' keeps its accuracy). Coefficients, aliased columns, covariance, residual
#' standard error, degrees of freedom, `logLik()`, `AIC()`, `BIC()` and
#' `summary()` are those of [lm()] on the same rows, up to rounding.
#'
#' Results that have one value per row are tables in the database (names of
#' temporary views, of class `xtx_remote_table`, with all the columns of the
#' remote data): `fitted(fit)` (columns `.fitted` and `.resid`; with
#' `na.exclude`, the rows with missing values too, with `NULL`), and
#' `predict(fit)` or `predict(fit, newdata = con, table = ...)` (column
#' `.fitted`; with `se.fit = TRUE`, `.se.fit`, and the attributes `df` and
#' `residual.scale`; with intervals, `.lwr` and `.upr`; for a rank-deficient
#' fit, `.nonestim`, the non-estimable cases flagged by `rankdeficient` as in
#' [stats::predict.lm()], with its warnings), computed by the same code as
#' in memory, in SQL.
#' `residuals(fit)` returns the remote design of the fit with two more
#' columns, `.fitted` and `.resid`, which can scale its rows: for
#' heteroskedasticity-consistent and clustered standard errors,
#' ```
#' R <- residuals(fit)
#' meat_hc1 <- crossprod(R, DBI::SQL('".resid" * ".resid"') * R)
#' scores <- colSums(".resid" * R, group = "cluster")
#' meat_cl <- crossprod(scores)
#' ```
#' (with `bread <- vcov(fit) / sigma(fit)^2`, the HC1 covariance is
#' `n / (n - k) * bread %*% meat_hc1 %*% bread`).
#'
#' **Differences from in-memory fits.**
#' * `model.frame()` is not available for a remote fit, and `model` is
#'   ignored: the data stay in the database (`fit$data$design` is the
#'   remote design of the fit).
#' * Results with one value per row are tables in the database, not
#'   vectors: `fitted()` (`.fitted`, `.resid`), `predict()` (`.fitted`;
#'   `se.fit` as a column `.se.fit` with attributes `df` and
#'   `residual.scale`, instead of a list; `.nonestim` instead of the
#'   attribute `"non-estim"`), and `residuals()`, a remote design.
#' * These tables have no row order: rows are matched to the data through
#'   its columns, which every table keeps.
#' * `compress = "auto"` decides from other evidence than in memory (see
#'   *Compression*); the results are the same either way.
#' @section GLMs:
#' [design_glm()] accepts remote data like [design_lm()]:
#' `design_glm(y ~ x + f, binomial, con, table = "d", weights = w, offset = log(e))`,
#' with `weights` and `offset` column names, `DBI::SQL()` or R expressions
#' over the columns, `offset()` terms in the formula, and a remote design
#' with `y` naming the response column (two columns for successes and
#' failures). The fit is the same code as in memory, on the remote methods
#' of the products: with few design patterns (see *Compression*), the
#' iterations run on the patterns; otherwise each iteration computes its
#' products in the database. `predict(fit, newdata = con,
#' table = ..., type = "response")`, `fitted()` and `residuals(fit, type =
#' )` return a temporary view with a `.fitted` (and `.resid`) column next to
#' the columns of the remote data. Fixed effects (`y ~ x | fe`, see
#' [design_glm()]) work remotely as in memory: each iteration's sums over
#' rows are queries, the working weights and responses a temporary table.
NULL

## The arguments of remote data that a fit on a formula takes in `...`
## (those design_problem() accepts in `dots`).
xtx_fit_remote_keys <- c("table", "sql", "engine", "chunk_size")

## The arguments in `...` of a fit on a formula, as a named list: those of
## remote data (xtx_fit_remote_keys) evaluated, the others NULL,
## unevaluated, for the check of unknown arguments (design_problem(); so
## that design_glm(subset = x > 0) reports `subset`, not that `x` is not
## found).
xtx_dots <- function(...) {
    n <- ...length()
    nm <- ...names() %||% rep("", n)
    out <- vector("list", n)
    names(out) <- nm
    for (i in which(nm %in% xtx_fit_remote_keys)) out[i] <- list(...elt(i))
    out
}

## The offsets of new remote data D for a fit: offset() terms of the design
## and the `offset` argument of the fit, as columns of D (design_resolve():
## SQL, or evaluated in R into new columns of the data of D), or NULL.
## Returns list(D, off).
xtx_remote_new_offset <- function(object, D) {
    tt <- D$plan$terms
    vars <- attr(tt, "variables")
    exprs <- lapply(attr(tt, "offset"), function(i) vars[[i + 1L]][[2L]])
    names(exprs) <- if (length(exprs)) paste0("offset_term", seq_along(exprs))
    if (!is.null(object$offset_call)) exprs$offset <- object$offset_call
    if (!length(exprs)) return(list(D = D, off = NULL))
    env <- object$env %||% environment(object$terms) %||% globalenv()
    r <- design_resolve(D, exprs, env)
    list(D = r$x, off = Reduce(`+`, r$cols))
}

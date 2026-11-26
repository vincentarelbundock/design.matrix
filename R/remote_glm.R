## design_glm() on remote data: the input side only. The fit is that of
## glm.R (dmglm_fit()) on the primitives of primitives.R. Here the response,
## prior weights and offset become columns over the rows of the remote
## design (dmp_resolve(): column names, DBI::SQL() or R expressions, those
## without a SQL translation evaluated in batches), and the rows where one
## is missing are dropped (dmp_subset()), as na.omit() drops them.

#' @name design_remote
#' @section GLMs:
#' [design_glm()] accepts remote data like [design_lm()]:
#' `design_glm(y ~ x + f, binomial, con, table = "d", weights = w, offset = log(e))`,
#' with `weights` and `offset` column names, `DBI::SQL()` or R expressions
#' over the columns, `offset()` terms in the formula, and a remote design
#' with `y` naming the response column (two columns for successes and
#' failures). The fit is the same code as in memory, on the remote methods
#' of the products: with few design patterns (see *Compression*), the
#' iterations run on the patterns; otherwise each iteration computes its
#' products in the database or in batches. `predict(fit, newdata = con,
#' table = ..., type = "response")`, `fitted()` and `residuals(fit, type =
#' )` return a temporary view with a `.fitted` (and `.resid`) column next to
#' the columns of the remote data.
NULL

## design_glm(formula, family, <remote>, ...): `weights` and `offset` are
## unevaluated arguments (or NULL), evaluated in env.
dm_remote_glm <- function(formula, source, family, weights, offset, start, contrasts, compress, control,
                          call, dots, na_fail, env) {
    dm_need_duckdb()
    k <- dmglm_kind(family)
    src <- dm_remote_source(source, dots$table, dots$sql)
    dots$table <- dots$sql <- NULL
    tt <- stats::terms(formula, data = src$proto)
    if (attr(tt, "response") == 0L) stop("The formula has no response.", call. = FALSE)
    lhs <- attr(tt, "variables")[[attr(tt, "response") + 1L]]
    ## the design of the right-hand side (offset() terms included)
    rhs <- stats::delete.response(tt)
    attr(rhs, "predvars") <- NULL
    rhs <- stats::formula(rhs)
    environment(rhs) <- environment(formula)
    D <- dm_remote(rhs, src, contrasts, NULL, dots, drop_unused = TRUE, env = env)
    D$src$holder <- src$holder
    if (na_fail) {
        tot <- as.numeric(dm_query(D, sprintf("SELECT COUNT(*) AS n FROM (%s) AS s", src$select))$n)
        if (tot != D$n) stop("missing values in object", call. = FALSE)
    }
    pr <- dmglm_remote_prepare(D, lhs, weights, offset, k, env, offset_terms = TRUE)
    if (na_fail && dm_nrows(pr$D) != D$n) stop("missing values in object", call. = FALSE)
    fit <- dmglm_fit(pr$D, pr$cols, family, start, control, attr(tt, "intercept") > 0L, compress,
                     wmsg = "negative weights not allowed")
    fit$call <- call
    fit$terms <- tt
    fit$formula <- formula
    fit$contrasts <- D$plan$contrasts
    fit$xlevels <- D$plan$xlevels
    fit$remote_offset_arg <- offset
    fit$remote_env <- env
    fit
}

## design_glm(<remote design>, y = , weights = , offset = ).
dm_remote_glm_design <- function(D, y, weights, offset, family, start, control, compress, env) {
    k <- dmglm_kind(family)
    if (!is.null(D$scale)) stop("Fit an unscaled remote design, with `weights` naming a column.", call. = FALSE)
    if (!identical(D$na, "omit"))
        stop("Fit a remote design built with `na.action = na.omit` (or from a formula).", call. = FALSE)
    lhs <- if (!is.null(y)) {
        dm_remote_column(D, y, "y", several = TRUE)
        if (length(y) == 2L) call("cbind", as.name(y[1L]), as.name(y[2L]))
        else if (length(y) == 1L) as.name(y)
        else stop("`y` must name one column, or two (successes and failures).", call. = FALSE)
    } else if (!is.null(D$response)) {
        D$vars[[D$response]]$expr
    } else {
        stop("With a remote design, give `y`, the name of the response column (or use a design with a response).",
             call. = FALSE)
    }
    pr <- dmglm_remote_prepare(D, lhs, weights, offset, k, env, offset_terms = FALSE)
    intercept <- attr(D$plan$terms, "intercept") > 0L && !isTRUE(D$plan$drop_intercept)
    fit <- dmglm_fit(pr$D, pr$cols, family, start, control, intercept, compress,
                     wmsg = "missing or negative weights not allowed")
    fit$remote_offset_arg <- offset
    fit$remote_env <- env
    fit
}

## The columns of the fit over the rows of D where the response, weights and
## offset are not missing: list(D, cols) with cols as dmglm_cols() makes
## them. The response `lhs`: an expression, cbind(successes, failures), or
## a character or factor column (success: any level but the first).
dmglm_remote_prepare <- function(D, lhs, weights, offset, k, env, offset_terms) {
    types <- D$src$types
    con <- dm_con(D)
    exprs <- list()
    cb <- is.call(lhs) && identical(lhs[[1L]], as.name("cbind"))
    if (cb) {
        if (length(lhs) != 3L) stop("A two-column response needs two columns (successes and failures).", call. = FALSE)
        exprs$a <- lhs[[2L]]
        exprs$b <- lhs[[3L]]
    } else if (is.name(lhs) && isTRUE(types[as.character(lhs)] %in% c("character", "factor"))) {
        col <- as.character(lhs)
        q <- DBI::dbQuoteIdentifier(con, col)
        vals <- dm_query(D, sprintf("SELECT DISTINCT s.%s AS v FROM (%s) AS s WHERE %s AND s.%s IS NOT NULL",
                                    q, D$src$select, D$filter, q))$v
        f <- if (is.factor(D$src$proto[[col]])) factor(as.character(vals), levels = levels(D$src$proto[[col]]))
             else factor(vals)
        exprs$y <- DBI::SQL(sprintf("CASE WHEN s.%s IS NULL THEN NULL WHEN CAST(s.%s AS VARCHAR) = %s THEN 0 ELSE 1 END",
                                    q, q, DBI::dbQuoteString(con, levels(droplevels(f))[1L])))
    } else if (is.name(lhs) && isTRUE(types[as.character(lhs)] == "logical")) {
        exprs$y <- DBI::SQL(sprintf("CAST(s.%s AS DOUBLE)", DBI::dbQuoteIdentifier(con, as.character(lhs))))
    } else {
        exprs$y <- lhs
    }
    exprs$weights <- weights
    exprs$offset <- offset
    r <- dmp_resolve(D, exprs, env)
    D <- r$x
    cols <- r$cols
    ## offset() terms of the formula
    off <- cols$offset
    if (offset_terms) {
        ot <- dm_remote_offset_sql(D)
        if (!is.null(ot)) off <- if (is.null(off)) dmp_column(D, ot) else off + dmp_column(D, ot)
    }
    ## rows where every column is present (NULL and NaN are missing)
    present <- function(v) dmc_is_finite(v) | (abs(v) == Inf)
    used <- Filter(Negate(is.null), c(cols[c("y", "a", "b", "weights")], list(off)))
    keep <- Reduce(`&`, lapply(used, present))
    D <- dmp_subset(D, keep)
    on <- function(v) if (is.null(v)) NULL else dmc_on(v, D)
    list(D = D, cols = dmglm_cols(on(cols$y), on(cols$a), on(cols$b),
                                  if (is.null(cols$weights)) dmc_const(D, 1) else on(cols$weights),
                                  on(off), k, cb, weighted = !is.null(cols$weights)))
}

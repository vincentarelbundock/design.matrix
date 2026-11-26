## The input side of fits on remote data (design_lm(), design_glm() and
## fixed effects; see remote.R): the design of the rows used, and the
## response, weights, offset and subset as columns over them (in SQL, or
## evaluated in R into a temporary table), with one code path for the three
## fits. The fits themselves are the code of lm.R, glm.R and fe.R on the
## primitives of primitives.R.

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
#' [na.fail()] or [na.pass()], with the results or errors of [lm()]. The
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
#' temporary views, of class `dm_remote_table`, with all the columns of the
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

## na.action of a fit, as one of `allowed` ("omit", "exclude", "fail",
## "pass").
dm_na_kind <- function(na.action, allowed = c("omit", "exclude", "fail", "pass"), what = "") {
    if (is.null(na.action)) return("omit")
    for (nm in allowed) {
        f <- get(paste0("na.", nm), envir = asNamespace("stats"))
        if (identical(na.action, f) || identical(na.action, paste0("na.", nm))) return(nm)
    }
    stop(sprintf("%s`na.action` must be %s.", what,
                 paste(paste0("na.", allowed), collapse = ", ")), call. = FALSE)
}

## The remote design of a fit and the columns of its arguments (one code
## path for design_lm(), design_glm() and fixed effects).
##   formula: the model formula; the design is that of its right-hand side
##     (offset() terms included; their values are columns too), and its
##     response is given by the caller in `args`.
##   src: an opened source (dm_remote_source()).
##   args: named list of unevaluated arguments (NULL: absent), resolved in
##     env: column names, DBI::SQL() or R expressions over the columns.
##     Numeric, except those named in `logical` (`subset`, which restricts
##     the rows).
##   present: columns that must not be missing (FE and cluster variables, a
##     factor response).
##   na: "omit" (or "exclude": the same rows), "fail" or "pass" (lm()'s
##     errors when a row has a missing value).
## R expressions with no SQL translation are evaluated in R into a
## temporary table (`engine = "sql"` refuses). The
## rows used are those of the subset where no model variable, argument or
## `present` column is missing, and factor levels are those of these rows,
## as model.frame(drop.unused.levels = TRUE) gives them in memory.
## Returns list(x = the design, cols = named list of columns (NULL for
## absent arguments), off = the offset (offset() terms plus the `offset`
## argument) or NULL).
dm_remote_fit <- function(formula, src, args, env, dots, contrasts = NULL, logical = character(),
                          present = character(), na = "omit") {
    miss <- setdiff(present, names(src$types))
    if (length(miss)) stop(sprintf("No column `%s` in the remote data.", miss[1L]), call. = FALSE)
    tt <- stats::terms(formula, data = src$proto)
    tv <- dm_terms_vars(tt)
    rhs <- if (attr(tt, "response") > 0L) stats::delete.response(tt) else tt
    attr(rhs, "predvars") <- NULL
    rhs <- stats::formula(rhs)
    environment(rhs) <- environment(formula)
    ## offset() terms are arguments too: their SQL, or the column computed
    ## in R, replaces the variable of the design (`override`)
    offt <- list()
    for (k in attr(tt, "offset")) offt[[tv$names[[k]]]] <- tv$exprs[[k]][[2L]]
    all <- c(Filter(Negate(is.null), args), offt)
    label <- function(nm) if (nm %in% names(offt)) "offset term" else if (nm %in% c("y", "a", "b")) "response" else nm
    con <- src$holder$con
    q <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
    res <- list()
    items <- list()
    exprs <- list()
    for (nm in names(all)) {
        a <- dm_resolve_or_eval(all[[nm]], env, src$types, con, label(nm), bool = nm %in% logical)
        if (is.null(a$expr)) {
            res[[nm]] <- a
            next
        }
        col <- paste0("__dm_eval", length(exprs) + 1L)
        items[[label(nm)]] <- a
        exprs[[col]] <- a$expr
        res[[nm]] <- list(sql = if (nm %in% logical) sprintf("s.%s", q(col)) else sprintf("CAST(s.%s AS DOUBLE)", q(col)),
                          inputs = col, logical = nm %in% logical)
    }
    if (length(exprs)) {
        dots <- dm_untranslated_engine(items, dots)
        lg <- unlist(lapply(res, function(r) if (isTRUE(r$logical)) r$inputs))
        src <- dm_materialize(src, exprs, env, logical = lg,
                              chunk_size = dots$chunk_size %||% getOption("design.matrix.chunk_size", 250000))
    }
    dots$table <- dots$sql <- NULL
    nonlog <- setdiff(names(res), c(logical, names(offt)))
    extra <- stats::setNames(lapply(res[nonlog], function(r) list(sql = r$sql, inputs = r$inputs, kind = "numeric")),
                             sprintf("(%s)", nonlog))
    for (col in present) extra[[paste0("(present:", col, ")")]] <- list(inputs = col, kind = "present")
    override <- lapply(res[names(offt)], function(r) list(sql = r$sql, inputs = r$inputs))
    sub <- if (any(logical %in% names(res))) paste(sprintf("(%s)", vapply(res[intersect(logical, names(res))], `[[`, "", "sql")), collapse = " AND ")
    D <- dm_remote(rhs, src, contrasts, NULL, dots, drop_unused = TRUE, subset = sub, extra = extra,
                   override = override)
    if (na %in% c("fail", "pass")) dm_remote_fit_na(D, res, na, nonlog, names(offt))
    cols <- lapply(names(args), function(nm) if (!is.null(res[[nm]])) dm_rcols(D, res[[nm]]$sql))
    names(cols) <- names(args)
    off <- c(res[["offset"]]$sql, vapply(res[names(offt)], `[[`, "", "sql"))
    list(x = D, cols = cols, off = if (length(off)) dm_rcols(D, paste(sprintf("(%s)", off), collapse = " + ")))
}

## The source of a fit on remote data, the terms of its formula and the
## response expression: list(src, terms, lhs).
dm_remote_model <- function(formula, source, dots) {
    dm_need_duckdb()
    src <- dm_remote_source(source, dots$table, dots$sql)
    tt <- stats::terms(formula, data = src$proto)
    if (attr(tt, "response") == 0L) stop("The formula has no response.", call. = FALSE)
    list(src = src, terms = tt, lhs = attr(tt, "variables")[[attr(tt, "response") + 1L]])
}

## na.fail and na.pass: the errors of lm() when rows have missing values.
dm_remote_fit_na <- function(D, res, na, args, offt) {
    total <- as.numeric(dm_query(D, sprintf("SELECT COUNT(*) AS n FROM (%s) AS s", D$src$select))$n)
    if (total == D$n) return(invisible())
    if (na == "fail") stop("missing values in object", call. = FALSE)
    bad <- function(nm) {
        r <- res[[nm]]
        if (is.null(r)) return(FALSE)
        as.numeric(dm_query(D, sprintf("SELECT COUNT(*) AS n FROM (%s) AS s WHERE (%s) IS NULL OR isnan(%s)",
                                       D$src$select, r$sql, r$sql))$n) > 0
    }
    if (any(vapply(c("y", "a", "b", "offset", offt), bad, NA))) stop("NA/NaN/Inf in 'y'", call. = FALSE)
    if (bad("weights")) stop("missing or negative weights not allowed", call. = FALSE)
    stop("NA/NaN/Inf in 'x'", call. = FALSE)
}

## design_lm(formula, <remote data>, ...). `weights`, `offset` and `subset`
## are unevaluated arguments, resolved in `env`.
dm_remote_lm <- function(formula, source, weights, offset, subset, na.action, contrasts,
                         compress, tol, call, dots, env) {
    na <- dm_na_kind(na.action, what = "For remote data, ")
    m <- dm_remote_model(formula, source, dots)
    src <- m$src
    tt <- m$terms
    r <- dm_remote_fit(formula, src, list(y = m$lhs, weights = weights, offset = offset, subset = subset),
                       env, dots, contrasts, logical = "subset", na = na)
    fit <- dm_remote_lm_fit(r$x, r$cols$y, r$cols$weights, r$off, compress, tol)
    if (na == "exclude") fit$na.action <- dm_remote_exclude(r$x)
    fit$offset_call <- offset
    fit$env <- env
    fit$call <- call
    fit$terms <- tt
    fit
}

## design_lm(<remote design>, y = , weights = , offset = ): unevaluated
## arguments resolved in `env` against the columns of the design's data.
dm_remote_lm_design <- function(D, y, weights, offset, compress, tol, env) {
    if (is.null(y)) {
        if (is.null(D$response))
            stop("With a remote design, give `y`, the name of one column (or use a design with a response).",
                 call. = FALSE)
        v <- D$vars[[D$response]]
        if (v$kind != "numeric")
            stop("The response of a remote fit must be a numeric variable with a SQL translation.", call. = FALSE)
        y <- DBI::SQL(v$sql)
    }
    r <- dmp_resolve(D, list(y = y, weights = weights, offset = offset), env)
    if (length(r$evaluated)) message(dm_untranslated_msg(r$evaluated), ": evaluated in R, in batches of rows.")
    fit <- dm_remote_lm_fit(r$x, r$cols$y, r$cols$weights, r$cols$offset, compress, tol)
    fit$offset_call <- offset
    fit$env <- env
    fit
}

## The fit of a remote design: the fitting core of lm.R on the remote
## primitives (primitives.R).
dm_remote_lm_fit <- function(D, y, w, off, compress = "auto", tol = 1e-07) {
    dmfit_check_compress(compress)
    intercept <- attr(D$plan$terms, "intercept") > 0L && !isTRUE(D$plan$drop_intercept)
    out <- dmfit_lm(D, y, w, off, compress, tol, intercept)
    out$method <- if (out$compressed) "patterns" else "SQL moments"
    out$spec <- D$plan
    out$contrasts <- D$plan$contrasts
    out$xlevels <- D$plan$xlevels
    out$data <- list(design = D, y = y, weights = w, offset = off)
    D$cache$n <- out$rows
    out
}

## design_glm(formula, family, <remote>, ...): `weights` and `offset` are
## unevaluated arguments (or NULL), evaluated in env. The response: an
## expression, cbind(successes, failures), or a character, factor or
## logical column (success: any level but the first).
dm_remote_glm <- function(formula, source, family, weights, offset, start, contrasts, compress, control,
                          call, dots, na, env) {
    k <- dmglm_kind(family)
    m <- dm_remote_model(formula, source, dots)
    src <- m$src
    tt <- m$terms
    resp <- dmglm_remote_response(m$lhs, src$types, src$holder$con)
    r <- dm_remote_fit(formula, src, c(resp$args, list(weights = weights, offset = offset)), env, dots,
                       contrasts, present = resp$present, na = na)
    D <- r$x
    y <- if (length(resp$present)) dmglm_remote_factor(D, resp$present) else r$cols$y
    cols <- dmglm_cols(y, r$cols$a, r$cols$b, if (is.null(r$cols$weights)) dmc_const(D, 1) else r$cols$weights,
                       r$off, k, resp$cbind, weighted = !is.null(r$cols$weights))
    fit <- dmglm_fit(D, cols, family, start, control, attr(tt, "intercept") > 0L, compress,
                     wmsg = "negative weights not allowed")
    if (na == "exclude") fit$na.action <- dm_remote_exclude(D)
    fit$call <- call
    fit$terms <- tt
    fit$formula <- formula
    fit$contrasts <- D$plan$contrasts
    fit$xlevels <- D$plan$xlevels
    fit$offset_call <- offset
    fit$env <- env
    fit
}

## The response of a remote GLM as arguments: list(args, present, cbind).
dmglm_remote_response <- function(lhs, types, con) {
    if (is.call(lhs) && identical(lhs[[1L]], as.name("cbind"))) {
        if (length(lhs) != 3L) stop("A two-column response needs two columns (successes and failures).", call. = FALSE)
        return(list(args = list(a = lhs[[2L]], b = lhs[[3L]]), present = character(), cbind = TRUE))
    }
    if (is.name(lhs) && isTRUE(types[as.character(lhs)] %in% c("character", "factor")))
        return(list(args = list(), present = as.character(lhs), cbind = FALSE))
    if (is.name(lhs) && isTRUE(types[as.character(lhs)] == "logical"))
        lhs <- DBI::SQL(sprintf("CAST(s.%s AS DOUBLE)", DBI::dbQuoteIdentifier(con, as.character(lhs))))
    list(args = list(y = lhs), present = character(), cbind = FALSE)
}

## A character or factor response as a column: 0 for the first level of the
## rows used, 1 for the others.
dmglm_remote_factor <- function(D, col) {
    q <- dm_qi(D, col)
    vals <- dm_query(D, sprintf("SELECT DISTINCT s.%s AS v FROM (%s) AS s WHERE %s", q, D$src$select, D$filter))$v
    proto <- D$src$proto[[col]]
    f <- if (is.factor(proto)) factor(as.character(vals), levels = levels(proto)) else factor(vals)
    first <- DBI::dbQuoteString(dm_con(D), levels(droplevels(f))[1L])
    dm_rcols(D, sprintf("CASE WHEN CAST(s.%s AS VARCHAR) = %s THEN 0 ELSE 1 END", q, first))
}

## design_glm(<remote design>, y = , weights = , offset = ): `y` names one
## column, or two (successes and failures).
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
    resp <- dmglm_remote_response(lhs, D$src$types, dm_con(D))
    r <- dmp_resolve(D, c(resp$args, list(weights = weights, offset = offset)), env)
    D <- r$x
    y <- if (length(resp$present)) dmglm_remote_factor(D, resp$present) else r$cols$y
    cols <- dmglm_cols(y, r$cols$a, r$cols$b, if (is.null(r$cols$weights)) dmc_const(D, 1) else r$cols$weights,
                       r$cols$offset, k, resp$cbind, weighted = !is.null(r$cols$weights))
    dmglm_check_design(cols)
    intercept <- attr(D$plan$terms, "intercept") > 0L && !isTRUE(D$plan$drop_intercept)
    fit <- dmglm_fit(D, cols, family, start, control, intercept, compress,
                     wmsg = "missing or negative weights not allowed")
    fit$offset_call <- offset
    fit$env <- env
    fit
}

## na.exclude for a remote fit: results on the rows used, padded with the
## rows of the data that the fit dropped for missing values (NULL results),
## as napredict() pads in memory (see napredict.dm_remote_exclude()).
dm_remote_exclude <- function(D) structure(list(design = D), class = "dm_remote_exclude")

#' @exportS3Method stats::napredict
napredict.dm_remote_exclude <- function(omit, x, ...) {
    if (!inherits(x, "dm_remote_table")) return(x)
    D0 <- omit$design
    con <- attr(x, "connection")
    sel <- sprintf("SELECT * FROM %s UNION ALL BY NAME SELECT s.* FROM (%s) AS s WHERE NOT COALESCE(%s, FALSE)",
                   DBI::dbQuoteIdentifier(con, unclass(x)[1L]), D0$src$select, D0$filter)
    name <- dm_tmp_name("dm_out")
    DBI::dbExecute(con, sprintf("CREATE TEMPORARY VIEW %s AS %s", DBI::dbQuoteIdentifier(con, name), sel))
    structure(name, class = "dm_remote_table", connection = con, holder = attr(x, "holder"))
}

#' @exportS3Method stats::naresid
naresid.dm_remote_exclude <- function(omit, x, ...) napredict.dm_remote_exclude(omit, x, ...)

## The offsets of new remote data D for a fit: offset() terms of the plan
## and the `offset` argument of the fit, resolved against the columns of D
## (SQL), or NULL. Returns list(D, off): offsets with no SQL translation are
## evaluated in R into new columns of the data of D (dm_materialize()).
dm_remote_new_offset <- function(object, D) {
    con <- dm_con(D)
    tt <- D$plan$terms
    nms <- dm_terms_vars(tt)$names[attr(tt, "offset")]
    pieces <- character()
    exprs <- list()
    add <- function(e) {
        col <- paste0("__dm_offn", length(exprs) + 1L)
        exprs[[col]] <<- e
        pieces <<- c(pieces, sprintf("CAST(s.%s AS DOUBLE)", DBI::dbQuoteIdentifier(con, col)))
    }
    for (nm in nms) {
        v <- D$vars[[nm]]
        if (v$kind == "numeric") pieces <- c(pieces, v$sql) else add(v$expr[[2L]])
    }
    arg <- object$offset_call
    env <- object$env %||% environment(object$terms) %||% globalenv()
    if (!is.null(arg)) {
        a <- dm_resolve_or_eval(arg, env, D$src$types, con, "offset")
        if (is.null(a$expr)) pieces <- c(pieces, a$sql) else add(a$expr)
    }
    if (length(exprs)) D$src <- dm_materialize(D$src, exprs, env, chunk_size = D$chunk_size)
    list(D = D, off = if (length(pieces)) dm_rcols(D, paste(sprintf("(%s)", pieces), collapse = " + ")))
}

## R expressions over the remote data with no SQL translation (a response,
## an offset() term, the `weights`, `offset` or `subset` of a fit): they are
## evaluated in R, batch by batch, and written with all the columns of the
## source to a temporary table, whose new columns the fit then reads like
## any other. One pass over the data; the rows, levels, missing values and
## compressed sums are then exact, as for translated expressions.
##   exprs: named list of R expressions (the names of the new columns);
##   logical: names of the expressions to keep as logical (subset);
##   model: for the expressions that are model variables (named as in
##     exprs), list(inputs, factor, label), see dm_eval_model(). Their
##     columns and kinds are returned in the attribute "model" of the
##     result: list(kind, cols) per expression.
## Returns the new source (dm_remote_source() of the temporary table),
## which keeps the original one alive.
dm_materialize <- function(src, exprs, env, logical = character(), chunk_size = 250000, model = list()) {
    info <- list()
    con <- src$holder$con
    q <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
    from <- src$base_table
    raw <- NULL
    if (is.null(from)) {
        raw <- q(dm_tmp_name("dm_raw"))
        DBI::dbExecute(con, sprintf("CREATE TEMPORARY TABLE %s AS %s", raw, src$select))
        on.exit(DBI::dbExecute(con, sprintf("DROP TABLE IF EXISTS %s", raw)), add = TRUE)
        from <- raw
    }
    r <- DBI::dbGetQuery(con, sprintf("SELECT MIN(rowid) AS lo, MAX(rowid) AS hi FROM %s", from))
    out <- dm_tmp_name("dm_eval")
    created <- FALSE
    write <- function(df) {
        if (!created) {
            DBI::dbWriteTable(con, out, df, temporary = TRUE)
            created <<- TRUE
        } else {
            DBI::dbAppendTable(con, out, df)
        }
    }
    evaluate <- function(chunk) {
        for (nm in names(exprs)) {
            if (!is.null(model[[nm]])) {
                r <- dm_eval_model(exprs[[nm]], chunk, env, model[[nm]], nm, info[[nm]])
                info[[nm]] <<- r$info
                for (cn in names(r$cols)) chunk[[cn]] <- r$cols[[cn]]
                next
            }
            v <- eval(exprs[[nm]], chunk, env)
            if (length(v) == 1L && nrow(chunk) != 1L) v <- rep(v, nrow(chunk))
            if (length(v) != nrow(chunk) || !is.null(dim(v)))
                stop(sprintf("`%s` gives %d values for %d rows.", deparse2(exprs[[nm]]), length(v), nrow(chunk)),
                     call. = FALSE)
            if (nm %in% logical) {
                if (!is.logical(v)) stop(sprintf("`%s` is not logical.", deparse2(exprs[[nm]])), call. = FALSE)
            } else {
                if (!is.numeric(v) && !is.logical(v))
                    stop(sprintf("`%s` is not numeric.", deparse2(exprs[[nm]])), call. = FALSE)
                v <- as.double(v)
            }
            chunk[[nm]] <- v
        }
        chunk
    }
    lo <- as.numeric(r$lo)
    if (!is.na(lo)) {
        hi <- as.numeric(r$hi)
        a <- lo
        while (a <= hi) {
            chunk <- DBI::dbGetQuery(con, sprintf("SELECT * FROM %s WHERE rowid >= %.0f AND rowid < %.0f",
                                                  from, a, a + chunk_size))
            if (nrow(chunk)) write(evaluate(chunk))
            a <- a + chunk_size
        }
    }
    if (!created) write(evaluate(src$proto))
    new <- dm_remote_source(con, out, NULL)
    new$parent <- src$holder
    new$label <- paste(src$label, "(with columns computed in R)")
    attr(new, "model") <- info
    new
}

## The value of a model variable (the R expression `expr`) on a batch of
## rows, as columns named after `nm`: one numeric column, or one per column
## of a numeric matrix (nm_1, nm_2, ...: poly(), ns()); a text column for a
## factor or character value, or when the specification holds levels for
## the variable (m$factor); a logical column. The expression sees only the
## rows where its input columns (m$inputs) are present; the others get
## missing values. m$class: the class of the variable in training (one
## named element of the plan's data classes), checked on each batch.
## `prev`: the kind and columns found on earlier batches,
## which every batch must give. Returns list(cols, info = list(kind, cols)).
dm_eval_model <- function(expr, chunk, env, m, nm, prev = NULL) {
    n <- nrow(chunk)
    ok <- if (length(m$inputs)) stats::complete.cases(chunk[m$inputs]) else rep(TRUE, n)
    bad <- function(msg) stop(sprintf("`%s` %s.", m$label, msg), call. = FALSE)
    if (!any(ok) && !is.null(prev)) {
        v <- NULL
        info <- prev
    } else {
        v <- eval(expr, chunk[ok, , drop = FALSE], env)
        ## a variable of another type than in training: the error of
        ## predict.lm() (a character variable with levels is a factor, as
        ## model.frame() makes it)
        if (!is.null(m$class))
            stats::.checkMFClasses(m$class, stats::setNames(list(if (is.character(v) && isTRUE(m$factor)) factor(v) else v),
                                                            names(m$class)))
        k <- if (is.matrix(v)) ncol(v) else 1L
        if (!is.matrix(v) && length(v) == 1L && sum(ok) != 1L) v <- rep(v, sum(ok))
        if (NROW(v) != sum(ok)) bad(sprintf("gives %d values for %d rows", NROW(v), sum(ok)))
        kind <- if (is.matrix(v)) {
            if (!is.numeric(v) && !is.logical(v)) bad("gives a matrix that is not numeric")
            "numeric"
        } else if (isTRUE(m$factor) || is.factor(v) || is.character(v)) {
            "factor"
        } else if (is.logical(v)) {
            "logical"
        } else if (is.numeric(v)) {
            "numeric"
        } else {
            bad(sprintf("gives values of class %s, which a remote design cannot use", class(v)[1L]))
        }
        info <- list(kind = kind, cols = if (k == 1L) nm else paste0(nm, "_", seq_len(k)))
        if (!is.null(prev) && !identical(prev, info)) bad("gives values of different kinds on different rows")
    }
    fill <- function(x, na) {
        out <- rep(na, n)
        if (!is.null(v)) out[ok] <- x
        out
    }
    cols <- lapply(seq_along(info$cols), function(u) {
        switch(info$kind,
               numeric = fill(if (is.matrix(v)) as.double(v[, u]) else as.double(v), NA_real_),
               factor = fill(as.character(v), NA_character_),
               logical = fill(v, NA))
    })
    list(cols = stats::setNames(cols, info$cols), info = info)
}

## Resolve an argument (dm_resolve_arg()), or, when its R expression has no
## SQL translation, return it for evaluation in R: list(sql, inputs) or
## list(expr, reason).
dm_resolve_or_eval <- function(x, env, types, con, what, bool = FALSE) {
    if (is.null(x)) return(NULL)
    if (is.call(x) && !deparse2(x[[1L]]) %in% c("SQL", "DBI::SQL")) {
        tr <- if (bool) dm_sql_bool(x, types, con) else dm_sql_num(x, types, con)
        if (is.null(tr$sql)) return(list(expr = x, reason = tr$reason))
    }
    dm_resolve_arg(x, env, types, con, what, bool)
}

## The expressions of a fit that SQL cannot compute, as one message:
## "weights `myfun(x)` cannot be translated to SQL (myfun() has no SQL
## translation)".
dm_untranslated_msg <- function(items) {
    nm <- names(items)[1L]
    sprintf("%s `%s` cannot be translated to SQL (%s)", nm, deparse2(items[[1L]]$expr), items[[1L]]$reason)
}

## A fit with expressions evaluated in R: a message names the first one;
## `engine = "sql"` refuses.
dm_untranslated_engine <- function(items, dots) {
    if (!length(items)) return(dots)
    why <- dm_untranslated_msg(items)
    if (identical(dots$engine, "sql")) stop("The SQL engine cannot compute this fit: ", why, ".", call. = FALSE)
    message("Evaluated in R: ", why, ".")
    dots
}

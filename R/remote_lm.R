## design_lm() on remote data (see remote.R, remote_compress.R, lm.R): the
## arguments of a remote fit (expressions in SQL or evaluated in R,
## na.action, offsets) and the tables of its row-wise results. The fit
## itself is the core of lm.R (dmfit_lm()) on the primitives of
## primitives.R: compressed patterns when they pay, otherwise qr() of the
## remote design (X'X in the database, suspect columns judged from the
## data) with coefficients refined from the residuals, as in memory.

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
#' (one more pass over the data; with `engine = "auto"`, the chunked engine
#' is used, with a message naming the expression, and `engine = "sql"`
#' refuses). Compression stays exact, since the computed values are columns
#' of that table. `na.action` can be [na.omit()] (default), [na.exclude()],
#' [na.fail()] or [na.pass()], with the results or errors of [lm()].
#'
#' With few design patterns, the fit is compressed (see *Compression*);
#' otherwise it is computed from `X'WX`, `X'Wy` and `y'Wy` with [lm()]'s
#' column selection, numeric variables centered when that is an exact
#' transform of the design. Coefficients, covariance, residual standard
#' error, degrees of freedom, `logLik()`, `AIC()`, `BIC()` and `summary()`
#' are those of [lm()] on the same rows, up to rounding.
#'
#' Results that have one value per row are tables in the database (names of
#' temporary views, of class `dm_remote_table`, with all the columns of the
#' remote data): `fitted(fit)` (columns `.fitted` and `.resid`; with
#' `na.exclude`, the rows with missing values too, with `NULL`), and
#' `predict(fit)` or `predict(fit, newdata = con, table = ...)` (column
#' `.fitted`; with `se.fit = TRUE`, `.se.fit`, and the attributes `df` and
#' `residual.scale`; with intervals, `.lwr` and `.upr`; for a rank-deficient
#' fit, `.nonestim`, the non-estimable cases flagged by `rankdeficient` as in
#' [stats::predict.lm()], with its warnings). Standard errors, intervals and
#' non-estimability are computed in batches of rows; the fit alone in SQL.
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
#'   ignored: the data stay in the database (`fit$remote_design` is the
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
#' * Without compression, aliased columns are detected from the centered
#'   `X'WX`, with [lm()]'s rule but without a refinement from the data: the
#'   same columns as [lm()] unless a column's residual is within rounding
#'   of `tol`.
NULL

## na.action of a remote fit: "omit", "exclude", "fail" or "pass".
dm_remote_lm_na <- function(na.action) {
    if (is.null(na.action)) return("omit")
    if (is.character(na.action)) na.action <- get(na.action, mode = "function")
    for (nm in c("omit", "exclude", "fail", "pass"))
        if (identical(na.action, get(paste0("na.", nm), envir = asNamespace("stats")))) return(nm)
    stop("For remote data, `na.action` must be na.omit, na.exclude, na.fail or na.pass.", call. = FALSE)
}

## SQL of the offsets of a remote design: offset() terms of its plan, and
## the `offset` argument (the "(offset)" variable). NULL without offsets.
dm_remote_offset_sql <- function(D) {
    tt <- D$plan$terms
    nms <- dm_terms_vars(tt)$names[attr(tt, "offset")]
    nms <- c(nms, intersect("(offset)", names(D$vars)))
    if (!length(nms)) return(NULL)
    sqls <- vapply(D$vars[nms], function(v) {
        if (v$kind != "numeric")
            stop(sprintf("The offset `%s` has no SQL translation.", v$name), call. = FALSE)
        v$sql
    }, "")
    paste(sprintf("(%s)", sqls), collapse = " + ")
}

## design_lm(formula, <remote data>, ...). `weights`, `offset` and `subset`
## are unevaluated arguments (see dm_resolve_arg()), resolved in `env`.
dm_remote_lm <- function(formula, source, weights, offset, subset, na.action, contrasts,
                         compress, tol, call, dots, env) {
    na <- dm_remote_lm_na(na.action)
    offset_arg <- offset
    dm_need_duckdb()
    src <- dm_remote_source(source, dots$table, dots$sql)
    dots$table <- dots$sql <- NULL
    types <- src$types
    con <- src$holder$con
    ## expressions with no SQL translation: evaluated in R into new columns
    tt <- stats::terms(formula, data = src$proto)
    tv <- dm_terms_vars(tt)
    items <- list()
    exprs <- list()
    override <- list()
    col_sql <- function(col) sprintf("CAST(s.%s AS DOUBLE)", DBI::dbQuoteIdentifier(con, col))
    resp <- attr(tt, "response")
    if (resp > 0L) {
        e <- tv$exprs[[resp]]
        tr <- dm_sql_num(e, types, con)
        if (is.null(tr$sql) && !(is.name(e) && isTRUE(types[as.character(e)] %in% c("character", "factor", "logical")))) {
            items[["response"]] <- list(expr = e, reason = tr$reason)
            exprs[["__dm_y"]] <- e
            override[[tv$names[[resp]]]] <- list(sql = col_sql("__dm_y"), inputs = "__dm_y")
        }
    }
    for (k in attr(tt, "offset")) {
        e <- tv$exprs[[k]][[2L]]
        tr <- dm_sql_num(e, types, con)
        if (is.null(tr$sql)) {
            col <- paste0("__dm_offt", k)
            items[["offset term"]] <- list(expr = tv$exprs[[k]], reason = tr$reason)
            exprs[[col]] <- e
            override[[tv$names[[k]]]] <- list(sql = col_sql(col), inputs = col)
        }
    }
    a <- dm_resolve_or_eval(weights, env, types, con, "weights")
    if (!is.null(a$expr)) {
        items[["weights"]] <- a
        exprs[["__dm_w"]] <- a$expr
        weights <- "__dm_w"
    }
    a <- dm_resolve_or_eval(offset, env, types, con, "offset")
    if (!is.null(a$expr)) {
        items[["offset"]] <- a
        exprs[["__dm_o"]] <- a$expr
        offset <- "__dm_o"
    }
    a <- dm_resolve_or_eval(subset, env, types, con, "subset", bool = TRUE)
    if (!is.null(a$expr)) {
        items[["subset"]] <- a
        exprs[["__dm_subset"]] <- a$expr
        subset <- "__dm_subset"
    }
    if (length(exprs)) {
        dots <- dm_untranslated_engine(items, dots)
        src <- dm_materialize(src, exprs, env, logical = "__dm_subset",
                              chunk_size = dots$chunk_size %||% getOption("design.matrix.chunk_size", 250000))
    }
    D <- dm_remote(formula, src, contrasts, NULL, dots, weights = weights, drop_unused = TRUE,
                   subset = subset, offset = offset, env = env, override = override)
    if (is.null(D$response)) stop("The formula has no response.", call. = FALSE)
    ## rows dropped for missing values: na.fail and na.pass refuse them
    if (na %in% c("fail", "pass")) {
        total <- as.numeric(dm_query(D, sprintf("SELECT COUNT(*) AS n FROM (%s) AS s", D$src$select))$n)
        if (total > D$n) {
            if (na == "fail") stop("missing values in object", call. = FALSE)
            ## na.pass: the error of the fit on rows with missing values
            bad <- function(v) {
                if (is.null(v) || v$kind != "numeric") return(FALSE)
                as.numeric(dm_query(D, sprintf("SELECT COUNT(*) AS n FROM (%s) AS s WHERE NOT COALESCE(%s, FALSE)",
                                               D$src$select, dm_remote_var_ok(D, v)))$n) > 0
            }
            offs <- dm_terms_vars(D$plan$terms)$names[attr(D$plan$terms, "offset")]
            if (bad(D$vars[[D$response]]) || any(vapply(D$vars[c(offs, intersect("(offset)", names(D$vars)))], bad, NA)))
                stop("NA/NaN/Inf in 'y'", call. = FALSE)
            if (bad(D$vars[["(weights)"]])) stop("missing or negative weights not allowed", call. = FALSE)
            stop("NA/NaN/Inf in 'x'", call. = FALSE)
        }
    }
    v <- D$vars[[D$response]]
    if (v$kind != "numeric")
        stop("The response of a remote fit must be a numeric variable with a SQL translation.", call. = FALSE)
    off <- dm_remote_offset_sql(D)
    z <- if (is.null(off)) v$sql else sprintf("(%s) - (%s)", v$sql, off)
    w <- D$vars[["(weights)"]]$sql
    fit <- dm_remote_lm_fit(D, z, w, compress, tol)
    fit$remote_na <- na
    fit$remote_offset <- off
    fit$remote_offset_arg <- offset_arg
    fit$remote_env <- env
    fit$call <- call
    fit$terms <- D$terms_full
    fit$contrasts <- D$plan$contrasts
    fit$xlevels <- D$plan$xlevels
    fit
}

## design_lm(<remote design>, y = , weights = , offset = ): unevaluated
## arguments resolved in `env` against the columns of the design's data.
dm_remote_lm_design <- function(D, y, weights, offset, compress, tol, env) {
    types <- D$src$types
    con <- dm_con(D)
    ## expressions with no SQL translation: evaluated in R into new columns
    ## of the design's data
    items <- list()
    exprs <- list()
    for (arg in c("y", "weights", "offset")) {
        a <- dm_resolve_or_eval(get(arg), env, types, con, arg)
        if (!is.null(a$expr)) {
            items[[arg]] <- a
            exprs[[paste0("__dm_", arg)]] <- a$expr
            assign(arg, paste0("__dm_", arg))
        }
    }
    if (length(exprs)) {
        message(dm_untranslated_msg(items), ": evaluated in R, in batches of rows.")
        D$src <- dm_materialize(D$src, exprs, env, chunk_size = D$chunk_size)
        types <- D$src$types
    }
    if (is.null(y)) {
        if (is.null(D$response))
            stop("With a remote design, give `y`, the name of one column (or use a design with a response).",
                 call. = FALSE)
        v <- D$vars[[D$response]]
        if (v$kind != "numeric")
            stop("The response of a remote fit must be a numeric variable with a SQL translation.", call. = FALSE)
        ysql <- v$sql
    } else {
        ysql <- dm_resolve_arg(y, env, types, con, "y")$sql
    }
    o <- dm_resolve_arg(offset, env, types, con, "offset")
    w <- dm_resolve_arg(weights, env, types, con, "weights")
    z <- if (is.null(o)) ysql else sprintf("(%s) - (%s)", ysql, o$sql)
    fit <- dm_remote_lm_fit(D, z, w$sql, compress, tol)
    fit$remote_na <- "omit"
    fit$remote_offset <- o$sql
    fit
}

## The fit of a remote design on z (the response minus offsets, SQL over
## `s`) with weights w (SQL over `s`, or NULL): the fitting core of lm.R on
## the remote primitives (primitives.R).
dm_remote_lm_fit <- function(D, z, w, compress = "auto", tol = 1e-07) {
    dmfit_check_compress(compress)
    intercept <- attr(D$plan$terms, "intercept") > 0L
    out <- dmfit_lm(D, dm_rcols(D, z), if (!is.null(w)) dm_rcols(D, w), NULL, compress, tol, intercept)
    out$method <- if (out$compressed) "patterns" else if (D$engine == "sql") "SQL moments"
                  else "moments computed in batches"
    out$spec <- D$plan
    out$remote <- TRUE
    out$remote_design <- D
    out$remote_z <- z
    out$remote_w <- w
    D$cache$n <- out$rows
    out
}







## The offsets of new remote data D for a remote or in-memory fit: offset()
## terms of the plan and the `offset` argument of the fit, resolved
## against the columns of D (SQL), or NULL.
## Returns list(D, off): offsets with no SQL translation are evaluated in R
## into new columns of the data of D (dm_materialize()).
dm_remote_new_offset <- function(object, D) {
    con <- dm_con(D)
    tt <- D$plan$terms
    nms <- dm_terms_vars(tt)$names[attr(tt, "offset")]
    pieces <- character()
    exprs <- list()
    for (nm in nms) {
        v <- D$vars[[nm]]
        if (v$kind == "numeric") {
            pieces <- c(pieces, v$sql)
        } else {
            col <- paste0("__dm_offn", length(exprs) + 1L)
            exprs[[col]] <- v$expr[[2L]]
            pieces <- c(pieces, sprintf("CAST(s.%s AS DOUBLE)", DBI::dbQuoteIdentifier(con, col)))
        }
    }
    arg <- if (isTRUE(object$remote)) object$remote_offset_arg else object$offset_call
    env <- (if (isTRUE(object$remote)) object$remote_env else environment(object$terms)) %||% globalenv()
    if (!is.null(arg)) {
        a <- dm_resolve_or_eval(arg, env, D$src$types, con, "offset")
        if (is.null(a$expr)) {
            pieces <- c(pieces, a$sql)
        } else {
            col <- paste0("__dm_offn", length(exprs) + 1L)
            exprs[[col]] <- a$expr
            pieces <- c(pieces, sprintf("CAST(s.%s AS DOUBLE)", DBI::dbQuoteIdentifier(con, col)))
        }
    }
    if (length(exprs)) D$src <- dm_materialize(D$src, exprs, env, chunk_size = D$chunk_size)
    list(D = D, off = if (length(pieces)) paste(sprintf("(%s)", pieces), collapse = " + ") else NULL)
}

## Fitted values and residuals of a remote fit, or predictions on remote
## data, as a table: the columns of the remote data, `.fitted` (with the
## offset), and `.resid` (training rows, with `resid = TRUE`), `.se.fit`,
## `.lwr`, `.upr` and `.nonestim` when asked. With na.exclude, the training
## rows with missing values are kept, with NULL. Standard errors, intervals
## and non-estimability are computed in batches (the chunked engine's
## reader), the rest in SQL when the design's engine is SQL.
dm_remote_lm_table <- function(object, newdata = NULL, na.action = stats::na.pass, dots = list(),
                               resid = FALSE, se.fit = FALSE, interval = "none", level = 0.95,
                               rankdeficient = "warnif", tol = 1e-06) {
    beta <- dmfit_beta(object)
    training <- is.null(newdata)
    if (training) {
        if (!isTRUE(object$remote)) stop("This fit was not computed on remote data.", call. = FALSE)
        D <- object$remote_design
        off <- object$remote_offset
        z <- object$remote_z
        w <- object$remote_w
    } else {
        if (is.null(object$spec)) stop("This design cannot be applied to new data.", call. = FALSE)
        D <- do.call(design_matrix, c(list(object$spec, newdata, na.action = na.action), dots))
        no <- dm_remote_new_offset(object, D)
        D <- no$D
        off <- no$off
        z <- w <- NULL
    }
    D <- dm_remote_unscaled(D)
    srccols <- names(D$src$types)
    helpers <- c(`__dm_off` = off, `__dm_z` = z, `__dm_w` = w)
    if (length(helpers)) D <- dm_remote_extend(D, helpers, "offset")
    rankdef <- !training && object$rank < length(object$coefficients)
    if (rankdef && rankdeficient == "simple")
        warning("prediction from rank-deficient fit; consider predict(., rankdeficient=\"NA\")", call. = FALSE)
    batch <- se.fit || interval != "none" || (rankdef && rankdeficient != "simple")
    has <- function(h) h %in% names(helpers)
    if (!batch) {
        v <- dm_remote_predict(D, beta)
        items <- c(sprintf("%s AS \".fitted\"", if (has("__dm_off")) "\".fitted\" + \"__dm_off\"" else "\".fitted\""),
                   if (resid) "\"__dm_z\" - \".fitted\" AS \".resid\"")
        sel <- sprintf("SELECT %s FROM %s", paste(c(dm_qi(D, srccols), items), collapse = ", "), dm_qi(D, unclass(v)[1L]))
    } else {
        res_var <- object$sigma^2
        V <- dmfit_vcov(object, 1)
        V[is.na(V)] <- 0
        tq <- stats::qt((1 - level) / 2, object$df.residual)
        nonest_n <- 0
        if (interval == "prediction") {
            if (training) {
                warning("predictions on current data refer to _future_ responses\n", call. = FALSE)
                if (!is.null(w))
                    warning("assuming prediction variance inversely proportional to weights used for fitting\n",
                            call. = FALSE)
            } else if (object$weighted) {
                warning("Assuming constant prediction variance even though model fit is weighted\n", call. = FALSE)
            }
        }
        out <- dm_tmp_name("dm_pred")
        created <- FALSE
        con <- dm_con(D)
        dm_chunks(D, unique(c(srccols, names(helpers), dm_remote_inputs(D))), function(chunk) {
            cd <- dm_chunk_design(D, chunk)
            X <- cd$design
            k <- cd$keep
            xb <- dmfit_linpred(object, X)
            fit <- if (has("__dm_off")) xb + chunk$`__dm_off`[k] else xb
            res <- chunk[k, srccols, drop = FALSE]
            if (rankdef && rankdeficient != "simple") {
                ne <- seq_along(k) %in% dmfit_nonestimable(object, X, tol)
                nonest_n <<- nonest_n + sum(ne)
                if (startsWith(rankdeficient, "NA")) fit[ne] <- NA else res$.nonestim <- ne
            }
            res$.fitted <- fit
            if (resid) res$.resid <- chunk$`__dm_z`[k] - xb
            if (se.fit || interval != "none") {
                se <- sqrt((if (ncol(X)) quadratic_diagonal(X, V) else rep(0, length(k))) * res_var)
                if (se.fit) res$.se.fit <- se
                if (interval != "none") {
                    pv <- if (interval == "prediction" && training && !is.null(w)) res_var / chunk$`__dm_w`[k] else res_var
                    hw <- tq * if (interval == "confidence") se else sqrt(se^2 + pv)
                    res$.lwr <- fit + hw
                    res$.upr <- fit - hw
                }
            }
            attr(res, "row.names") <- c(NA_integer_, -length(k))
            if (!created) {
                DBI::dbWriteTable(con, out, res, temporary = TRUE)
                created <<- TRUE
            } else {
                DBI::dbAppendTable(con, out, res)
            }
        })
        if (!created) {
            res <- D$src$proto[srccols]
            res$.fitted <- numeric()
            DBI::dbWriteTable(con, out, res, temporary = TRUE)
        }
        if (rankdef && nonest_n > 0) {
            if (rankdeficient == "warnif")
                warning("prediction from rank-deficient fit; attr(*, \"non-estim\") has doubtful cases", call. = FALSE)
            if (rankdeficient == "NAwarn")
                warning("prediction from rank-deficient fit: NAs produced for non-estimable cases", call. = FALSE)
        }
        sel <- sprintf("SELECT * FROM %s", dm_qi(D, out))
    }
    ## na.exclude: the training rows with missing values, with NULL
    if (training && identical(object$remote_na, "exclude")) {
        D0 <- object$remote_design
        sel <- sprintf("%s UNION ALL BY NAME SELECT s.* FROM (%s) AS s WHERE NOT COALESCE(%s, FALSE)",
                       sel, D0$src$select, D0$filter)
    }
    name <- dm_tmp_name("dm_pred")
    dm_execute(D, sprintf("CREATE TEMPORARY VIEW %s AS %s", dm_qi(D, name), sel))
    out <- structure(name, class = "dm_remote_table", connection = dm_con(D), holder = D$src$holder)
    if (se.fit) {
        attr(out, "df") <- object$df.residual
        attr(out, "residual.scale") <- object$sigma
    }
    out
}

## R expressions over the remote data with no SQL translation (a response,
## an offset() term, the `weights`, `offset` or `subset` of a fit): they are
## evaluated in R, batch by batch, and written with all the columns of the
## source to a temporary table, whose new columns the fit then reads like
## any other. One pass over the data; the rows, levels, missing values and
## compressed sums are then exact, as for translated expressions.
##   exprs: named list of R expressions (the names of the new columns);
##   logical: names of the expressions to keep as logical (subset).
## Returns the new source (dm_remote_source() of the temporary table),
## which keeps the original one alive.
dm_materialize <- function(src, exprs, env, logical = character(), chunk_size = 250000) {
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
    new
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

## The engine for a fit with expressions evaluated in R: "auto" switches to
## the chunked engine with a message, "sql" refuses.
dm_untranslated_engine <- function(items, dots) {
    if (!length(items)) return(dots)
    why <- dm_untranslated_msg(items)
    engine <- dots$engine %||% "auto"
    if (identical(engine, "sql")) stop("The SQL engine cannot compute this fit: ", why, ".", call. = FALSE)
    if (identical(engine, "auto")) message("Using the chunked engine: ", why, ".")
    dots$engine <- "chunked"
    dots
}

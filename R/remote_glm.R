## design_glm() on remote data: the input side only. The fit is that of
## glm.R (dmglm_fit()) on the primitives of primitives.R; what is remote
## here is how the response, prior weights, trials and offset become
## columns: SQL over the source, computed as glm()'s initialize sees them
## (y in [0, 1] and weights times trials for a two-column binomial response,
## y = 0 where the weight is 0), with the rows where one is missing dropped,
## as na.omit() drops them.

#' @name design_remote
#' @section GLMs:
#' [design_glm()] accepts remote data like [design_lm()]:
#' `design_glm(y ~ x + f, binomial, con, table = "d", weights = "w", offset = "o")`,
#' with `weights` and `offset` naming columns, `offset()` terms in the
#' formula, and a remote design with `y` naming the response column (two
#' columns for successes and failures). The response, weights and offset
#' must have SQL translations. The fit is the same code as in memory, on
#' the remote methods of the products: with few design patterns (see
#' *Compression*), the iterations run on the patterns; otherwise each
#' iteration computes its products in the database or in batches.
#' `predict(fit, newdata = con, table = ..., type = "response")`, `fitted()`
#' and `residuals(fit, type = )` return the name of a temporary view with a
#' `.fitted` (and `.resid`) column next to the columns of the remote data.
NULL

## design_glm(formula, family, <remote>, ...): `weights` and `offset` are
## column names or NULL.
dm_remote_glm <- function(formula, source, family, weights, offset, start, contrasts, compress, control,
                          call, dots, na_fail = FALSE) {
    dm_need_duckdb()
    k <- dmglm_kind(family)
    nms <- names(dots) %||% rep("", length(dots))
    src <- dm_remote_source(source, dots$table, dots$sql)
    con <- src$holder$con
    tt <- stats::terms(formula, data = src$proto)
    if (attr(tt, "response") == 0L) stop("The formula has no response.", call. = FALSE)
    vars <- attr(tt, "variables")
    lhs <- vars[[attr(tt, "response") + 1L]]
    off_exprs <- lapply(attr(tt, "offset"), function(i) vars[[i + 1L]][[2L]])
    labs <- attr(tt, "term.labels")
    int <- attr(tt, "intercept") > 0L
    rhs <- if (length(labs)) stats::reformulate(labs, intercept = int) else if (int) ~1 else ~0
    environment(rhs) <- environment(formula)
    rc <- dmglm_remote_sql(con, src$types, src$proto, src$select, lhs, weights, off_exprs, offset, k)
    ## the rows: those where the response, weights and offset are not missing
    sel <- sprintf("SELECT s.*, %s FROM (%s) AS s WHERE %s",
                   paste(sprintf("%s AS %s", rc$exprs, DBI::dbQuoteIdentifier(con, names(rc$exprs))), collapse = ", "),
                   src$select, rc$notnull)
    D <- dm_remote(rhs, con, contrasts, NULL, c(list(sql = sel), dots[!nms %in% c("table", "sql")]),
                   drop_unused = TRUE)
    ## keep a private connection (Parquet) open as long as the design
    D$src$holder <- src$holder
    if (na_fail) {
        tot <- as.numeric(dm_query(D, sprintf("SELECT COUNT(*) AS n FROM (%s) AS s", src$select))$n)
        if (tot != D$n) stop("missing values in object", call. = FALSE)
    }
    fit <- dmglm_fit(D, dmglm_remote_cols(D, rc), family, start, control, int, compress,
                     wmsg = "negative weights not allowed")
    fit$call <- call
    fit$terms <- tt
    fit$formula <- formula
    fit$contrasts <- D$plan$contrasts
    fit$xlevels <- D$plan$xlevels
    if (!is.null(offset)) fit$offset_call <- as.name(offset)
    fit
}

## design_glm(<remote design>, y = , weights = , offset = ).
dm_remote_glm_design <- function(D, y, weights, offset, family, start, control, compress) {
    k <- dmglm_kind(family)
    if (!is.null(D$scale)) stop("Fit an unscaled remote design, with `weights` naming a column.", call. = FALSE)
    if (!identical(D$na, "omit"))
        stop("Fit a remote design built with `na.action = na.omit` (or from a formula).", call. = FALSE)
    if (!is.null(weights)) dm_remote_column(D, weights, "weights")
    if (!is.null(offset)) dm_remote_column(D, offset, "offset")
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
    rc <- dmglm_remote_sql(dm_con(D), D$src$types, D$src$proto, D$src$select, lhs, weights, list(), offset, k)
    D <- dm_remote_extend(D, rc$exprs, "response")
    D$filter <- sprintf("(%s) AND %s", D$filter, rc$notnull)
    D$cache <- new.env(parent = emptyenv())
    D$n <- as.numeric(dm_query(D, sprintf("SELECT COUNT(*) AS n FROM (%s) AS s WHERE %s", D$src$select, D$filter))$n)
    intercept <- attr(D$plan$terms, "intercept") > 0L && !isTRUE(D$plan$drop_intercept)
    fit <- dmglm_fit(D, dmglm_remote_cols(D, rc), family, start, control, intercept, compress,
                     wmsg = "missing or negative weights not allowed")
    if (!is.null(offset)) fit$offset_call <- as.name(offset)
    fit
}

## The columns of the fit (dmglm_fit()) from the source columns of
## dmglm_remote_sql().
dmglm_remote_cols <- function(D, rc) {
    col <- function(nm) dmp_column(D, sprintf("s.%s", dm_qi(D, nm)))
    list(y = col("__dm_y"), w = col("__dm_w"), n = col("__dm_n"), off = if (rc$offset) col("__dm_off"),
         cbind = rc$cbind, a = if (rc$cbind) col("__dm_ca"), b = if (rc$cbind) col("__dm_cb"),
         weighted = rc$weighted)
}

## SQL of the response, weights, trials and offset (over the columns of
## `s`), and the condition that none is missing.
dmglm_remote_sql <- function(con, types, proto, select, lhs, weights, off_exprs, off_col, k) {
    tr <- function(e, what) {
        r <- dm_sql_num(e, types, con)
        if (!is.null(r$reason))
            stop(sprintf("The %s `%s` has no SQL translation (%s).", what, deparse2(e), r$reason), call. = FALSE)
        r$sql
    }
    nn <- function(x) sprintf("((%s) IS NOT NULL AND NOT isnan(%s))", x, x)
    notnull <- character()
    w <- "CAST(1 AS DOUBLE)"
    if (!is.null(weights)) {
        if (!isTRUE(types[weights] %in% c("double", "integer", "logical")))
            stop("'weights' must be a numeric vector", call. = FALSE)
        w <- sprintf("CAST(s.%s AS DOUBLE)", DBI::dbQuoteIdentifier(con, weights))
        notnull <- c(notnull, nn(w))
    }
    cb <- is.call(lhs) && identical(lhs[[1L]], as.name("cbind"))
    n <- "CAST(1 AS DOUBLE)"
    raw <- list()
    if (cb) {
        if (length(lhs) != 3L || !k$binomial)
            stop("A two-column response (successes and failures) needs a binomial family.", call. = FALSE)
        a <- tr(lhs[[2L]], "response")
        b <- tr(lhs[[3L]], "response")
        notnull <- c(notnull, nn(a), nn(b))
        n <- sprintf("(%s + %s)", a, b)
        y <- sprintf("(CASE WHEN %s = 0 THEN 0 ELSE %s / %s END)", n, a, n)
        w <- sprintf("(%s * %s)", w, n)
        raw <- list(a = a, b = b)
    } else if (is.name(lhs) && isTRUE(types[as.character(lhs)] %in% c("character", "factor"))) {
        ## a factor response: success is any level but the first, as in binomial()
        if (!k$binomial) stop("A factor response needs a binomial family.", call. = FALSE)
        col <- as.character(lhs)
        q <- DBI::dbQuoteIdentifier(con, col)
        vals <- DBI::dbGetQuery(con, sprintf("SELECT DISTINCT s.%s AS v FROM (%s) AS s WHERE s.%s IS NOT NULL",
                                             q, select, q))$v
        f <- if (is.factor(proto[[col]])) factor(as.character(vals), levels = levels(proto[[col]])) else factor(vals)
        first <- levels(droplevels(f))[1L]
        y <- sprintf("(CASE WHEN CAST(s.%s AS VARCHAR) = %s THEN 0.0 ELSE 1.0 END)", q, DBI::dbQuoteString(con, first))
        notnull <- c(notnull, sprintf("s.%s IS NOT NULL", q))
    } else {
        y <- tr(lhs, "response")
        notnull <- c(notnull, nn(y))
    }
    if (k$binomial && !cb) y <- sprintf("(CASE WHEN %s = 0 THEN 0 ELSE %s END)", w, y)
    off <- dmglm_offset_sql(con, types, off_exprs, off_col)
    if (!is.null(off)) notnull <- c(notnull, nn(off))
    exprs <- c(`__dm_y` = y, `__dm_w` = w, `__dm_n` = n)
    if (!is.null(off)) exprs <- c(exprs, `__dm_off` = off)
    if (cb) exprs <- c(exprs, `__dm_ca` = raw$a, `__dm_cb` = raw$b)
    list(exprs = exprs, notnull = if (length(notnull)) paste(notnull, collapse = " AND ") else "TRUE",
         cbind = cb, offset = !is.null(off), weighted = !is.null(weights))
}

## The offset as SQL: offset() terms of the formula and an offset column.
dmglm_offset_sql <- function(con, types, off_exprs, off_col) {
    pieces <- vapply(off_exprs, function(e) {
        r <- dm_sql_num(e, types, con)
        if (!is.null(r$reason))
            stop(sprintf("The offset `%s` has no SQL translation (%s).", deparse2(e), r$reason), call. = FALSE)
        r$sql
    }, "")
    if (!is.null(off_col)) {
        if (!isTRUE(types[off_col] %in% c("double", "integer", "logical")))
            stop("The offset must be a numeric column.", call. = FALSE)
        pieces <- c(pieces, sprintf("CAST(s.%s AS DOUBLE)", DBI::dbQuoteIdentifier(con, off_col)))
    }
    if (!length(pieces)) return(NULL)
    paste0("(", paste(pieces, collapse = " + "), ")")
}

## The offset of new remote data for predict(), as a column of its design.
dm_remote_glm_offset <- function(object, x) {
    tt <- object$terms
    off_exprs <- if (is.null(tt)) list() else {
        vars <- attr(tt, "variables")
        lapply(attr(tt, "offset"), function(i) vars[[i + 1L]][[2L]])
    }
    col <- if (!is.null(object$offset_call)) as.character(object$offset_call)
    off <- dmglm_offset_sql(dm_con(x), x$src$types, off_exprs, col)
    if (is.null(off)) NULL else dmp_column(x, off)
}

## Out-of-core designs: the rows of the design live in a DuckDB table, a SQL
## query or Parquet files, and never all at once in R.
##
## A remote design holds the source (connection and SQL), a design_plan
## (terms with fixed transformations, levels, contrasts), a description of
## each model variable (how to compute it, or its missingness, in SQL) and
## the engine that computes products:
##
##  * "sql": each term t has X_t[i, ] = N_t(i) C_t[F_t(i), ], with F_t(i)
##    the combined level of its factors and N_t(i) the product of its
##    numeric variables (see structured.R). Sums over rows of products of
##    such columns are GROUP BY aggregates over factor codes; the small
##    dense products with the C_t are done in R (remote_sql.R).
##  * "chunked": batches of rows are read through DBI and handed to
##    design_matrix(plan, batch); the in-memory results are accumulated
##    (remote_chunked.R).

#' Out-of-core designs on DuckDB
#'
#' With a DBI connection to DuckDB, or the path of Parquet files, as `data`,
#' [design_matrix()] returns a *remote design*: the rows stay in the
#' database, and products with the design run there (or in batches), so that
#' data larger than memory can be used.
#'
#' ```
#' design_matrix(formula, con, table = "flights")
#' design_matrix(formula, con, sql = "SELECT * FROM flights WHERE year > 2000")
#' design_matrix(formula, "path/to/*.parquet")
#' design_matrix(spec, con, table = "new_flights")   # fit in R, score in the database
#' ```
#'
#' @section Sources:
#' `con` must be a [DBI][DBI::dbConnect()] connection to DuckDB (packages
#' DBI and duckdb). Give either `table`, the name of a table or view (a
#' character string or a [DBI::Id()]), or `sql`, a `SELECT` query. A
#' character vector of Parquet paths (globs allowed) opens a private
#' in-memory DuckDB connection that reads them with `read_parquet()`; it is
#' closed when the design is garbage collected.
#'
#' @section Training and new data:
#' With a formula, the design is learned from the remote data, as
#' `design_matrix(formula, data)` learns it from a data frame:
#' * Levels: those of an `ENUM` column (a factor in R) are its declared
#'   levels; those of a character column are its distinct values in the rows
#'   used, sorted as [factor()] sorts them in R (not by SQL collation).
#'   Logical columns are factors with levels `FALSE` and `TRUE`.
#' * Rows: rows with a missing value (SQL `NULL`, or `NaN`) in any model
#'   variable, response included, are dropped, as [stats::na.omit()] drops
#'   them.
#' * Transformations whose parameters depend on the data, such as
#'   [stats::poly()], [splines::ns()] or [scale()], are fitted as
#'   [stats::model.frame()] fits them, on the whole of their input columns:
#'   those columns (only those) are read into R once. To avoid this, fit
#'   the specification in R, on a sample or on other data, and pass it
#'   instead of the formula.
#' * Variables computed by an R expression that yields a factor, such as
#'   `factor(cyl)` or `interaction(a, b)`, take the levels that the
#'   expression gives on the distinct values of its input columns.
#'
#' With a specification (a design, a plan, a fitted model or another remote
#' design), the levels, contrasts and transformation parameters come from
#' it, and the remote data is only scored. Levels absent from the
#' specification raise the error that [stats::model.frame()] raises, and
#' columns of another type than in training (a numeric variable stored as
#' text, say) the error of [stats::predict.lm()]. As in memory, the rows
#' with missing values are kept by default (`na.action = na.pass`): their
#' predictions are `NULL`, and products involving their columns are `NA`.
#' With `na.action = na.omit`, they are dropped.
#'
#' @section Engines:
#' `engine = "sql"` computes everything in DuckDB. It covers factors,
#' character and logical columns, numeric columns, interactions, and
#' numeric expressions built from `+`, `-`, `*`, `/`, `^`, [I()], [log()]
#' (with base 2, 10 or a number), [log2()], [log10()], [exp()], [sqrt()] and
#' [abs()]. `engine = "chunked"` handles any formula: rows are read in
#' batches of `chunk_size` rows (from the table itself when it is a DuckDB
#' base table, otherwise from a temporary table holding the needed columns of
#' the rows used) and each batch goes through `design_matrix(spec, batch)`.
#' `engine = "auto"` (default) uses SQL when every variable can be
#' translated, and otherwise the chunked engine, with a message naming the
#' first term that cannot be translated. [design_explain()] shows the
#' decision for each term, and the generated SQL.
#'
#' In the chunked engine, a variable evaluated in R is taken as missing,
#' when selecting rows in SQL, if one of its input columns is missing; rows
#' where the R expression itself gives `NA` are then dropped in R. The
#' number of rows printed is that of the SQL selection until a pass over the
#' data has counted the rows kept.
#'
#' @section Operations:
#' * `crossprod(D)`: `X'X`.
#' * `crossprod(D, "y")`: `X'y` for a column (or several columns) of the
#'   remote data.
#' * [design_moments()]: `X'X`, `X'y`, `y'y` and the number of rows in
#'   one call, from the response of the formula by default, optionally
#'   weighted by a column (`weights = "w"`): what (weighted) least squares
#'   needs.
#' * `colSums(D)`, `colMeans(D)`, and with `group = "g"`, the sums or means
#'   within the groups of a column (sorted as [base::rowsum()] sorts them);
#'   `weights = "w"` weights the rows.
#' * `D %*% beta` or [design_predict()]: the linear predictor, written to a
#'   temporary view (SQL engine) or table (chunked engine) whose name is
#'   returned.
#' * `as.matrix(D)` collects the rows and builds the dense matrix, only when
#'   there are at most `getOption("design.matrix.remote_max_rows", 1e5)` rows
#'   or with `force = TRUE`.
#' * `design_matrix(D, newdata)` applies the design to a data frame, and
#'   `design_matrix(D, con, table = ...)` to other remote data.
#'
#' Remote designs keep only the *Algebra* promise of the package: results
#' agree with the in-memory computation up to floating point rounding. With
#' infinite values in the data, the placement of `NaN` and `Inf` in sums can
#' differ from the dense computation, and `NaN` may appear where the dense
#' result has `NA`.
#'
#' @param table,sql The remote rows: a table name, or a `SELECT` query.
#' @param engine `"auto"`, `"sql"` or `"chunked"`.
#' @param chunk_size Rows per batch for the chunked engine (default
#'   `getOption("design.matrix.chunk_size", 250000)`).
#' @examples
#' if (requireNamespace("duckdb", quietly = TRUE)) {
#'     con <- DBI::dbConnect(duckdb::duckdb())
#'     d <- data.frame(y = rnorm(100), x = rnorm(100), f = gl(4, 25))
#'     DBI::dbWriteTable(con, "d", d)
#'     D <- design_matrix(y ~ x * f, con, table = "d")
#'     D
#'     crossprod(D)
#'     m <- design_moments(D)
#'     solve(m$xtx, m$xty)
#'     p <- D %*% solve(m$xtx, m$xty)
#'     head(as.data.frame(p))
#'     DBI::dbDisconnect(con)
#' }
#' @name design_remote
NULL

## Is `x` remote data for design_matrix()? (A character vector is never
## valid data for model.frame(), so it can only mean Parquet paths.)
dm_is_remote_source <- function(x) {
    inherits(x, "DBIConnection") || (is.character(x) && length(x) > 0L && is.null(dim(x)))
}

dm_need_duckdb <- function() {
    for (p in c("DBI", "duckdb")) {
        if (!requireNamespace(p, quietly = TRUE))
            stop(sprintf("Remote designs need the '%s' package.", p), call. = FALSE)
    }
}


## Entry point from design_matrix().
dm_remote <- function(object, source, contrasts.arg, xlev, dots, na.action = stats::na.pass) {
    dm_need_duckdb()
    keys <- c("table", "sql", "engine", "chunk_size")
    nms <- names(dots)
    if (is.null(nms)) nms <- rep("", length(dots))
    ra <- dots[nms %in% keys]
    rest <- dots[!nms %in% keys]
    engine <- ra$engine %||% "auto"
    checkmate::assert_choice(engine, c("auto", "sql", "chunked"), .var.name = "engine")
    chunk_size <- ra$chunk_size %||% getOption("design.matrix.chunk_size", 250000)
    checkmate::assert_number(chunk_size, lower = 1, finite = TRUE, .var.name = "chunk_size")
    src <- dm_remote_source(source, ra$table, ra$sql)

    if (inherits(object, "design_formula")) {
        contrasts.arg <- object$contrasts
        xlev <- object$xlevels
        object <- object$formula
    }
    if (inherits(object, "formula")) {
        if (length(rest)) stop("Unknown arguments: ", paste(names(rest), collapse = ", "), call. = FALSE)
        D <- dm_remote_formula(object, src, contrasts.arg, xlev)
    } else {
        if (!is.null(contrasts.arg) || !is.null(xlev))
            stop("`contrasts.arg` and `xlev` apply only to a formula.", call. = FALSE)
        plan <- if (inherits(object, "design_plan")) {
            object
        } else if (inherits(object, "design_matrix")) {
            dm_plan_for_design(object)
        } else {
            do.call(design_plan, c(list(object), rest))
        }
        D <- dm_remote_spec(plan, src, dm_remote_na(na.action))
    }
    D$chunk_size <- as.integer(chunk_size)
    dm_remote_finish(D, engine)
}

## na.action of a specification applied to remote data: "pass" (rows with
## missing values kept, as by default in memory) or "omit".
dm_remote_na <- function(na.action) {
    if (identical(na.action, stats::na.pass) || identical(na.action, "na.pass")) return("pass")
    if (identical(na.action, stats::na.omit) || identical(na.action, "na.omit") ||
        identical(na.action, stats::na.exclude) || identical(na.action, "na.exclude")) return("omit")
    stop("For remote data, `na.action` must be na.pass or na.omit.", call. = FALSE)
}

## The source: connection (in a holder environment, so that a private
## connection is closed once nothing refers to it) and a SELECT statement.
dm_remote_source <- function(source, table, sql) {
    holder <- new.env(parent = emptyenv())
    base_table <- NULL
    if (is.character(source)) {
        if (!is.null(table) || !is.null(sql))
            stop("Give `table` or `sql` only with a connection.", call. = FALSE)
        holder$con <- DBI::dbConnect(duckdb::duckdb())
        holder$own <- TRUE
        reg.finalizer(holder, function(e) {
            try(DBI::dbDisconnect(e$con), silent = TRUE)
        }, onexit = TRUE)
        con <- holder$con
        paths <- paste(DBI::dbQuoteString(con, source), collapse = ", ")
        select <- sprintf("SELECT * FROM read_parquet([%s])", paths)
        label <- sprintf("parquet %s", paste(source, collapse = ", "))
    } else {
        if (!inherits(source, "duckdb_connection"))
            stop("Remote designs support DuckDB connections only.", call. = FALSE)
        holder$con <- source
        holder$own <- FALSE
        con <- source
        if (is.null(table) == is.null(sql))
            stop("Give exactly one of `table` and `sql`.", call. = FALSE)
        if (!is.null(table)) {
            qt <- as.character(DBI::dbQuoteIdentifier(con, table))
            select <- sprintf("SELECT * FROM %s", qt)
            label <- sprintf("table %s", qt)
            ## a base table can be read in batches by rowid, without a copy
            ok <- tryCatch({
                DBI::dbGetQuery(con, sprintf("SELECT rowid FROM %s LIMIT 0", qt))
                TRUE
            }, error = function(e) FALSE)
            if (ok) base_table <- qt
        } else {
            checkmate::assert_string(sql, .var.name = "sql")
            select <- sub(";[[:space:]]*$", "", sql)
            label <- "SQL query"
        }
    }
    proto <- DBI::dbGetQuery(con, sprintf("SELECT * FROM (%s) AS s LIMIT 0", select))
    types <- vapply(proto, dm_col_type, "")
    list(holder = holder, select = select, base_table = base_table, label = label,
         proto = proto, types = types)
}

dm_con <- function(D) D$src$holder$con

dm_col_type <- function(x) {
    if (is.factor(x)) return("factor")
    if (is.character(x)) return("character")
    if (is.logical(x)) return("logical")
    if (is.integer(x) && !is.object(x)) return("integer")
    if (is.double(x) && !is.object(x)) return("double")
    "other"
}

dm_qi <- function(D, x) as.character(DBI::dbQuoteIdentifier(dm_con(D), x))

## Run a query (with the option design.matrix.remote_debug, print it).
dm_query <- function(D, sql) {
    if (isTRUE(getOption("design.matrix.remote_debug"))) message(sql)
    DBI::dbGetQuery(dm_con(D), sql)
}

dm_execute <- function(D, sql) {
    if (isTRUE(getOption("design.matrix.remote_debug"))) message(sql)
    DBI::dbExecute(dm_con(D), sql)
}

dm_tmp_name <- function(prefix) {
    paste0(prefix, "_", paste(sample(c(letters, 0:9), 12, TRUE), collapse = ""))
}

## Describe one model variable: how SQL can compute it, or test whether it
## is missing.
##   kind: "factor" (a character or ENUM column), "logical" (a logical
##         column), "numeric" (an expression with a SQL translation) or
##         "other" (evaluated in R).
dm_remote_var <- function(name, expr, types, D, xlevels, formula_case) {
    v <- list(name = name, expr = expr, inputs = intersect(all.vars(expr), names(types)),
              kind = "other", reason = NULL, sql = NULL, col = NULL, levels = NULL)
    in_xlev <- name %in% names(xlevels)
    if (is.name(expr) && as.character(expr) %in% names(types)) {
        col <- as.character(expr)
        ty <- types[[col]]
        if (ty %in% c("character", "factor")) {
            if (formula_case || in_xlev) {
                v$kind <- "factor"
                v$col <- col
            } else {
                v$reason <- sprintf("`%s` is a character column but not a factor in the specification", col)
            }
            return(v)
        }
        if (ty == "logical") {
            if (formula_case || !in_xlev) {
                v$kind <- "logical"
                v$col <- col
            } else {
                v$reason <- sprintf("`%s` is a logical column but a factor in the specification", col)
            }
            return(v)
        }
    }
    if (in_xlev) {
        v$reason <- sprintf("`%s` is a factor in the specification", name)
        return(v)
    }
    tr <- dm_sql_num(expr, types, dm_con(D))
    if (is.null(tr$sql)) {
        v$reason <- tr$reason
    } else {
        v$kind <- "numeric"
        v$sql <- tr$sql
    }
    v
}

## SQL condition: is variable `v` present (not NULL, not NaN) in a row?
## Exact for variables with a SQL translation; a variable evaluated in R
## counts as present when all its input columns are.
dm_remote_var_ok <- function(D, v) {
    types <- D$src$types
    cond <- switch(v$kind,
        factor = , logical = sprintf("s.%s IS NOT NULL", dm_qi(D, v$col)),
        numeric = sprintf("(%s) IS NOT NULL AND NOT isnan(%s)", v$sql, v$sql),
        other = vapply(v$inputs, function(col) {
            q <- dm_qi(D, col)
            if (types[[col]] == "double") sprintf("s.%s IS NOT NULL AND NOT isnan(s.%s)", q, q)
            else sprintf("s.%s IS NOT NULL", q)
        }, ""))
    cond <- unique(cond)
    if (!length(cond)) "TRUE" else paste(sprintf("(%s)", cond), collapse = " AND ")
}

## WHERE clause: the rows used. With na.omit, the rows with no missing model
## variable; with na.pass (a specification scoring remote data), all rows.
dm_remote_filter <- function(D) {
    if (identical(D$na, "pass")) return("TRUE")
    cond <- unique(vapply(D$vars, function(v) dm_remote_var_ok(D, v), ""))
    cond <- setdiff(cond, "TRUE")
    if (!length(cond)) "TRUE" else paste(sprintf("(%s)", cond), collapse = " AND ")
}

## Variables of a terms object, with their predvars.
dm_terms_vars <- function(tt) {
    vars <- as.list(attr(tt, "variables"))[-1L]
    pv <- attr(tt, "predvars")
    exprs <- if (is.null(pv)) vars else as.list(pv)[-1L]
    list(names = vapply(vars, deparse2, ""), exprs = exprs, calls = vars)
}

## A remote design from a formula: learn levels and transformations from the
## remote data.
dm_remote_formula <- function(formula, src, contrasts.arg, xlev) {
    D <- list(src = src, na = "omit")
    tt <- stats::terms(formula, data = src$proto)
    tv <- dm_terms_vars(tt)
    types <- src$types
    D$vars <- Map(function(nm, e) dm_remote_var(nm, e, types, D, NULL, TRUE), tv$names, tv$exprs)
    names(D$vars) <- tv$names
    D$filter <- dm_remote_filter(D)

    ## a sample of the rows used, to see what the R expressions return
    sample <- dm_remote_sample(D)
    env <- environment(tt)
    predvars <- tv$exprs
    xlevels <- list()
    for (j in seq_along(D$vars)) {
        v <- D$vars[[j]]
        if (v$kind == "factor") {
            x <- if (types[[v$col]] == "factor") src$proto[[v$col]] else dm_remote_distinct(D, v$col)[[1L]]
            xlevels[[v$name]] <- levels(as.factor(x))
        } else if (v$kind == "other") {
            val <- tryCatch(eval(v$expr, sample, env), error = function(e) e)
            if (is.factor(val) || is.character(val)) {
                ## levels from the distinct values of the inputs
                inp <- dm_remote_distinct(D, v$inputs)
                lv <- levels(as.factor(eval(v$expr, inp, env)))
                xlevels[[v$name]] <- lv
            } else if (!is.logical(val)) {
                ## does model.frame() fix parameters from the data?
                changed <- inherits(val, "error") ||
                    !identical(stats::makepredictcall(val, tv$calls[[j]]), tv$calls[[j]])
                if (changed) {
                    full <- dm_query(D, sprintf("SELECT %s FROM (%s) AS s",
                                                paste0("s.", dm_qi(D, v$inputs), collapse = ", "),
                                                src$select))
                    predvars[[j]] <- stats::makepredictcall(eval(v$expr, full, env), tv$calls[[j]])
                }
            }
        }
    }
    for (nm in names(xlev)) xlevels[[nm]] <- xlev[[nm]]
    attr(tt, "predvars") <- as.call(c(list(as.name("list")), predvars))
    for (j in seq_along(D$vars)) D$vars[[j]]$expr <- predvars[[j]]
    if (!nrow(sample)) stop("The remote data has no complete rows.", call. = FALSE)
    D$plan <- design_plan.formula(tt, sample, contrasts.arg = contrasts.arg,
                                  xlev = if (length(xlevels)) xlevels else NULL)
    D$terms_full <- tt
    resp <- attr(tt, "response")
    D$response <- if (resp > 0L) tv$names[[resp]] else NULL
    D
}

## A remote design from a specification: score the remote data.
dm_remote_spec <- function(plan, src, na = "pass") {
    D <- list(src = src, plan = plan, terms_full = plan$terms, response = NULL, na = na)
    tv <- dm_terms_vars(plan$terms)
    xl <- plan$xlevels
    dm_remote_check_classes(plan, src, tv, xl)
    D$vars <- Map(function(nm, e) dm_remote_var(nm, e, src$types, D, xl, FALSE), tv$names, tv$exprs)
    names(D$vars) <- tv$names
    D$filter <- dm_remote_filter(D)
    D
}

## Columns of another type than in training (a numeric variable stored as
## text, say): the error of predict.lm(), before anything is built. Only
## plain columns are checked here; variables computed by R expressions are
## checked on each batch (dm_chunk_design()).
dm_remote_check_classes <- function(plan, src, tv, xl) {
    cl <- plan$data_classes
    if (is.null(cl)) return(invisible())
    m <- list()
    for (j in seq_along(tv$exprs)) {
        e <- tv$exprs[[j]]
        nm <- tv$names[[j]]
        if (!is.name(e) || !as.character(e) %in% names(src$proto) || !nm %in% names(cl)) next
        x <- src$proto[[as.character(e)]]
        ## model.frame() turns a character variable with levels into a factor
        if (is.character(x) && nm %in% names(xl)) x <- factor(x)
        m[[nm]] <- x
    }
    if (length(m)) stats::.checkMFClasses(cl, m)
    invisible()
}

## Levels, columns, engine.
dm_remote_finish <- function(D, engine) {
    plan <- D$plan
    xl <- plan$xlevels
    ## levels of the factors, and new levels (as model.frame() reports them)
    for (j in seq_along(D$vars)) {
        v <- D$vars[[j]]
        if (v$kind == "logical") D$vars[[j]]$levels <- c("FALSE", "TRUE")
        if (v$kind != "factor") next
        lv <- xl[[v$name]]
        D$vars[[j]]$levels <- lv
        x <- dm_remote_distinct(D, v$col)[[1L]]
        xi <- if (is.character(x)) as.factor(x) else x[, drop = TRUE]
        nxl <- levels(xi)
        m <- is.na(match(nxl, lv))
        if (any(m))
            stop(sprintf(ngettext(length(m), "factor %s has new level %s", "factor %s has new levels %s"),
                         v$name, paste(nxl[m], collapse = ", ")), call. = FALSE)
    }
    sample <- dm_remote_sample(D)
    proto <- dm_chunk_design(D, sample)$design
    D$assign <- if (is.null(proto$setup)) attr(proto$dense, "assign") else proto$setup$assign
    D$colnames <- colnames(proto)
    D$p <- length(D$colnames)
    D$n <- as.numeric(dm_query(D, sprintf("SELECT COUNT(*) AS n FROM (%s) AS s WHERE %s",
                                           D$src$select, D$filter))$n)
    D$cache <- new.env(parent = emptyenv())

    ## which terms SQL can compute
    D$terms <- dm_remote_terms(D)
    bad <- Filter(function(t) !t$sql_ok, D$terms)
    resp_bad <- !is.null(D$response) && D$vars[[D$response]]$kind == "other"
    ## variables in no term (an offset, say) decide which rows are used
    used <- unique(unlist(lapply(D$terms, function(t) t$vars)))
    loose <- setdiff(names(D$vars)[vapply(D$vars, function(v) v$kind == "other", NA)],
                     c(names(D$vars)[used], D$response))
    short <- if (length(bad)) {
        sprintf("term `%s` is evaluated in R", bad[[1L]]$label)
    } else if (resp_bad) {
        sprintf("response `%s` is evaluated in R", D$response)
    } else if (length(loose)) {
        sprintf("variable `%s` is evaluated in R", loose[1L])
    }
    why <- if (length(bad)) {
        sprintf("term `%s` cannot be translated to SQL (%s)", bad[[1L]]$label, bad[[1L]]$reason)
    } else if (resp_bad) {
        sprintf("response `%s` cannot be translated to SQL (%s)", D$response, D$vars[[D$response]]$reason)
    } else if (length(loose)) {
        sprintf("variable `%s` cannot be translated to SQL (%s)", loose[1L], D$vars[[loose[1L]]]$reason)
    } else {
        NULL
    }
    if (engine == "sql" && !is.null(why))
        stop("The SQL engine cannot compute this design: ", why, ".", call. = FALSE)
    if (engine == "auto") {
        engine <- if (is.null(why)) "sql" else "chunked"
        if (!is.null(why)) message("Using the chunked engine: ", why, ".")
    }
    D$engine <- engine
    D$why <- why
    D$why_short <- if (engine == "chunked") short %||% "engine = \"chunked\"" else NULL
    if (engine == "sql") D <- dm_sql_prepare(D)
    class(D) <- "dm_remote"
    D
}

## Distinct values of columns, in the rows used.
dm_remote_distinct <- function(D, cols) {
    dm_query(D, sprintf("SELECT DISTINCT %s FROM (%s) AS s WHERE %s",
                        paste0("s.", dm_qi(D, cols), collapse = ", "), D$src$select, D$filter))
}

## The columns that the model variables read.
dm_remote_inputs <- function(D) {
    unique(unlist(lapply(D$vars, function(v) v$inputs)))
}

dm_remote_sample <- function(D, n = 2000L) {
    cols <- dm_remote_inputs(D)
    sel <- if (length(cols)) paste0("s.", dm_qi(D, cols), collapse = ", ") else "1 AS dm_one"
    dm_query(D, sprintf("SELECT %s FROM (%s) AS s WHERE %s LIMIT %d", sel, D$src$select, D$filter, n))
}

## Terms of the design: their variables, and whether SQL computes them.
dm_remote_terms <- function(D) {
    tt <- D$plan$terms
    fac <- attr(tt, "factors")
    labels <- attr(tt, "term.labels")
    ## rows of `factors` are the variables, in order
    pv <- dm_terms_vars(tt)$names
    out <- list()
    if (attr(tt, "intercept") == 1L && !isTRUE(D$plan$drop_intercept))
        out[[1L]] <- list(index = 0L, label = "(Intercept)", vars = integer(), sql_ok = TRUE, reason = NULL)
    for (t in seq_along(labels)) {
        vnames <- pv[fac[, t] > 0]
        kinds <- vapply(D$vars[vnames], function(v) v$kind, "")
        reason <- NULL
        if (any(kinds == "other")) {
            k <- which(kinds == "other")[1L]
            reason <- D$vars[[vnames[k]]]$reason
        } else {
            L <- prod(vapply(D$vars[vnames], function(v) max(1L, length(v$levels)), 1L))
            if (L * sum(D$assign == t) > dm_sql_max_cells)
                reason <- sprintf("%s combined levels: too many for the SQL engine", format(L, big.mark = ","))
        }
        out[[length(out) + 1L]] <- list(index = t, label = labels[t], vars = match(vnames, names(D$vars)),
                                        sql_ok = is.null(reason), reason = reason)
    }
    out
}

#' @export
print.dm_remote <- function(x, ...) {
    cat(sprintf("<remote design_matrix: %s rows x %d columns>\n", dm_fmt_count(x$cache$n %||% x$n), x$p))
    cat("source: ", x$src$label, "\n", sep = "")
    if (x$engine == "sql") {
        cat("engine: sql (all terms pushed down)\n")
    } else {
        cat(sprintf("engine: chunked, %s rows per batch (%s)\n", dm_fmt_count(x$chunk_size), x$why_short))
    }
    cn <- x$colnames
    if (length(cn)) {
        shown <- utils::head(cn, 8L)
        more <- if (length(cn) > 8L) sprintf(" ... (%d more)", length(cn) - 8L) else ""
        cat("columns:", paste(shown, collapse = ", "), more, "\n")
    }
    invisible(x)
}

dm_fmt_count <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)

#' @export
dim.dm_remote <- function(x) c(x$cache$n %||% x$n, x$p)

#' @export
dimnames.dm_remote <- function(x) list(NULL, x$colnames)

#' @export
as.matrix.dm_remote <- function(x, force = FALSE, ...) {
    max_rows <- getOption("design.matrix.remote_max_rows", 1e5)
    if (!isTRUE(force) && nrow(x) > max_rows)
        stop(sprintf("The remote design has %s rows (more than %s): use `force = TRUE` to build the dense matrix in memory.",
                     dm_fmt_count(nrow(x)), dm_fmt_count(max_rows)), call. = FALSE)
    chunk <- dm_query(x, dm_chunk_sql(x, dm_remote_inputs(x)))
    as.matrix(dm_chunk_design(x, chunk)$design)
}

## A column name argument (y, weights, group): check and quote.
dm_remote_column <- function(D, x, what, several = FALSE) {
    if (several) checkmate::assert_character(x, min.len = 1L, any.missing = FALSE, .var.name = what)
    else checkmate::assert_string(x, .var.name = what)
    miss <- setdiff(x, names(D$src$types))
    if (length(miss))
        stop(sprintf("`%s`: no column %s in the remote data.", what, paste0("`", miss, "`", collapse = ", ")),
             call. = FALSE)
    x
}

## The response column of a formula design (for design_moments()).
dm_remote_y <- function(D, y) {
    if (!is.null(y)) return(list(cols = dm_remote_column(D, y, "y", several = TRUE), response = FALSE))
    if (is.null(D$response)) return(NULL)
    list(cols = D$response, response = TRUE)
}

#' Moments for least squares
#'
#' `design_moments()` returns what least squares needs from a design:
#' `X'X`, `X'y`, `y'y` and the number of rows, optionally weighted. For a
#' remote design (see [design_remote]), everything is computed in the
#' database (or in one pass of the chunked engine).
#'
#' @param x A design matrix, or a remote design.
#' @param y The response. For a remote design: `NULL` (default) for the
#'   response of the formula (if any; without one, only `X'X` is computed),
#'   or the names of one or more numeric columns.
#'   For a design in memory: a numeric vector or matrix with one row per row
#'   of `x`.
#' @param weights Optional weights: the name of a numeric column for a
#'   remote design, a numeric vector for a design in memory.
#' @return A list with `xtx` (`X'WX`), `xty` (`X'Wy`, one column per
#'   response, or `NULL`), `yty` (`y'Wy`, a matrix, or `NULL`), `n` (rows)
#'   and `sumw` (sum of the weights, `n` without weights).
#' @examples
#' d <- data.frame(y = rnorm(20), x = rnorm(20), f = gl(2, 10))
#' m <- design_moments(design_matrix(y ~ x + f, d), d$y)
#' solve(m$xtx, m$xty)
#' @export
design_moments <- function(x, y = NULL, weights = NULL) {
    if (inherits(x, "dm_remote")) {
        yy <- dm_remote_y(x, y)
        if (!is.null(weights)) dm_remote_column(x, weights, "weights")
        out <- if (x$engine == "sql") dm_sql_moments(x, yy, weights) else dm_chunked_moments(x, yy, weights)
        x$cache$n <- out$n
        return(out)
    }
    checkmate::assert_class(x, "design_matrix")
    if (is.null(y)) stop("Give `y` for a design in memory.", call. = FALSE)
    Y <- if (is.matrix(y)) y else matrix(y, ncol = 1L, dimnames = list(NULL, "y"))
    checkmate::assert_matrix(Y, mode = "numeric", nrows = nrow(x), .var.name = "y")
    if (is.null(weights)) {
        xtx <- crossprod(x)
        xty <- crossprod(x, Y)
        yty <- crossprod(Y)
        sumw <- nrow(x)
    } else {
        checkmate::assert_numeric(weights, len = nrow(x))
        xtx <- crossprod(x, weights * x)
        xty <- crossprod(x, weights * Y)
        yty <- crossprod(Y, weights * Y)
        sumw <- sum(weights)
    }
    list(xtx = xtx, xty = xty, yty = yty, n = nrow(x), sumw = sumw)
}

#' Linear predictions written to the database
#'
#' `design_predict(D, beta)` computes `X %*% beta` for a remote design (see
#' [design_remote]) and writes it, with columns of the remote data, to a new
#' view or table. `D %*% beta` is `design_predict(D, beta)`.
#'
#' With the SQL engine and no `name`, the result is a temporary view, which
#' computes the predictions each time it is read; with `name`, a table is
#' created and filled. The chunked engine writes a temporary table, or the
#' table `name`.
#'
#' @param x A remote design.
#' @param beta Coefficients: a numeric vector of length `ncol(x)`, or a
#'   matrix with `ncol(x)` rows (one output column per column).
#' @param name Name of a table to create and keep, or `NULL` for a temporary
#'   view or table that disappears with the connection.
#' @param column Name of the prediction column (with several columns in
#'   `beta`, their column names, or `column` followed by 1, 2, ...).
#' @param keep Columns of the remote data to copy next to the predictions,
#'   to join them back. `NULL` (default) keeps all the columns.
#' @param overwrite Replace an existing table `name`?
#' @return The name of the view or table, a character string of class
#'   `dm_remote_table`; [as.data.frame()] reads it into R, and its
#'   `"connection"` attribute is the connection where it lives (useful for
#'   a Parquet source, which has a private connection).
#' @export
design_predict <- function(x, beta, name = NULL, column = ".fitted", keep = NULL, overwrite = FALSE) {
    checkmate::assert_class(x, "dm_remote")
    p <- x$p
    B <- if (is.matrix(beta)) beta else matrix(beta, ncol = 1L)
    checkmate::assert_matrix(B, mode = "numeric", nrows = p, any.missing = FALSE, .var.name = "beta")
    checkmate::assert_true(all(is.finite(B)), .var.name = "all(is.finite(beta))")
    storage.mode(B) <- "double"
    checkmate::assert_string(name, null.ok = TRUE)
    checkmate::assert_string(column)
    checkmate::assert_flag(overwrite)
    keep <- if (is.null(keep)) names(x$src$types) else dm_remote_column(x, keep, "keep", several = TRUE)
    outcols <- if (ncol(B) == 1L) column else colnames(B) %||% paste0(column, seq_len(ncol(B)))
    if (any(outcols %in% keep))
        stop("The prediction column clashes with a kept column; change `column` or `keep`.", call. = FALSE)
    out <- if (x$engine == "sql") {
        dm_sql_predict(x, B, name, outcols, keep, overwrite)
    } else {
        dm_chunked_predict(x, B, name, outcols, keep, overwrite)
    }
    structure(out, class = "dm_remote_table", connection = dm_con(x), holder = x$src$holder)
}

#' @export
print.dm_remote_table <- function(x, ...) {
    cat(sprintf("<remote table \"%s\">: read it with as.data.frame(), or query it on attr(, \"connection\")\n",
                unclass(x)[1L]))
    invisible(x)
}

#' @export
as.data.frame.dm_remote_table <- function(x, ...) {
    con <- attr(x, "connection")
    DBI::dbGetQuery(con, sprintf("SELECT * FROM %s", DBI::dbQuoteIdentifier(con, unclass(x)[1L])))
}

#' @export
`%*%.dm_remote` <- function(x, y) {
    if (!inherits(x, "dm_remote"))
        stop("Only `design %*% beta` is supported, with the design first.", call. = FALSE)
    design_predict(x, y)
}

#' @method crossprod dm_remote
#' @export
crossprod.dm_remote <- function(x, y = NULL, ...) {
    if (!inherits(x, "dm_remote"))
        stop("Only `crossprod(design)` and `crossprod(design, \"column\")` are supported, with a design first.",
             call. = FALSE)
    ## crossprod() is a primitive: its dispatch drops argument names, so
    ## that `crossprod(D, weights = "w")` would arrive as `y = "w"`
    nm <- names(sys.call())
    if (length(list(...)) || (!is.null(nm) && any(!nm %in% c("", "x", "y"))))
        stop("`crossprod()` of a remote design takes only `y`; for weights, use `design_moments(x, weights = )`.",
             call. = FALSE)
    if (is.null(y)) {
        return(if (x$engine == "sql") dm_sql_gram(x) else dm_chunked_moments(x, NULL, NULL)$xtx)
    }
    if (!is.character(y))
        stop("For a remote design, `y` must name columns of the remote data.", call. = FALSE)
    dm_remote_column(x, y, "y", several = TRUE)
    yy <- list(cols = y, response = FALSE)
    if (x$engine == "sql") dm_sql_xty(x, yy) else dm_chunked_moments(x, yy, NULL, gram = FALSE)$xty
}

dm_remote_colsums <- function(x, group, weights, reorder, mean) {
    if (!is.null(weights)) dm_remote_column(x, weights, "weights")
    if (!isTRUE(reorder))
        stop("Remote designs support only `reorder = TRUE`.", call. = FALSE)
    if (!is.null(group)) dm_remote_column(x, group, "group")
    if (x$engine == "sql") dm_sql_colsums(x, group, weights, mean)
    else dm_chunked_colsums(x, group, weights, mean)
}

setOldClass("dm_remote")

#' @rdname design_remote
#' @param x A remote design.
#' @param na.rm,dims Only the defaults are supported.
#' @param ... Not used.
#' @param group For `colSums()` and `colMeans()`, the name of a column of the
#'   remote data whose values define the groups.
#' @param weights The name of a numeric column of weights.
#' @param reorder Only `TRUE` (groups sorted) is supported.
#' @export
setMethod("colSums", "dm_remote", function(x, na.rm = FALSE, dims = 1, ..., group = NULL,
                                           weights = NULL, reorder = TRUE) {
    dm_defaults_only(na.rm, dims)
    dm_remote_colsums(x, group, weights, reorder, mean = FALSE)
})

#' @rdname design_remote
#' @export
setMethod("colMeans", "dm_remote", function(x, na.rm = FALSE, dims = 1, ..., group = NULL,
                                            weights = NULL, reorder = TRUE) {
    dm_defaults_only(na.rm, dims)
    dm_remote_colsums(x, group, weights, reorder, mean = TRUE)
})

#' Explain how a remote design is computed
#'
#' Lists, for each term of a remote design (see [design_remote]), whether
#' the SQL engine can compute it and why not, and optionally prints the SQL
#' that the engine runs.
#'
#' @param x A remote design.
#' @param sql Print the generated SQL: the query that computes the model
#'   variables and, for the SQL engine, the aggregates of `crossprod(x)`;
#'   for the chunked engine, the query that reads a batch.
#' @return A data frame of class `dm_explain` with one row per term
#'   (`term`, `sql`, `reason`), printed with the engine of the design (and
#'   the SQL, with `sql = TRUE`).
#' @export
design_explain <- function(x, sql = FALSE) {
    checkmate::assert_class(x, "dm_remote")
    checkmate::assert_flag(sql)
    tab <- data.frame(
        term = vapply(x$terms, function(t) t$label, ""),
        sql = vapply(x$terms, function(t) t$sql_ok, NA),
        reason = vapply(x$terms, function(t) t$reason %||% "", ""),
        stringsAsFactors = FALSE)
    attr(tab, "engine") <- if (x$engine == "sql") "sql" else sprintf("chunked (%s)", x$why_short)
    if (sql) {
        attr(tab, "sql") <- if (x$engine == "sql") {
            c(paste0("-- model variables\n", dm_sql_base(x)),
              paste0("-- crossprod\n", vapply(dm_sql_gram_queries(x, NULL), function(q) q$sql, "")))
        } else {
            paste0("-- batch\n", dm_chunk_sql(x, dm_remote_inputs(x), range = c(0, x$chunk_size)))
        }
    }
    class(tab) <- c("dm_explain", class(tab))
    tab
}

#' @export
print.dm_explain <- function(x, ...) {
    cat("engine: ", attr(x, "engine"), "\n", sep = "")
    sql <- attr(x, "sql")
    tab <- x
    attr(tab, "engine") <- attr(tab, "sql") <- NULL
    class(tab) <- "data.frame"
    print(tab, row.names = FALSE)
    for (q in sql) cat("\n", q, "\n", sep = "")
    invisible(x)
}

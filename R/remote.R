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
#' first term that cannot be translated. `summary(D)` shows the decision
#' for each term, and `summary(D, sql = TRUE)` the generated SQL.
#'
#' In the chunked engine, a variable evaluated in R is taken as missing,
#' when selecting rows in SQL, if one of its input columns is missing; rows
#' where the R expression itself gives `NA` are then dropped in R. The
#' number of rows printed is that of the SQL selection until a pass over the
#' data has counted the rows kept.
#'
#' @section Operations:
#' The operations of a design in memory (see [design_algebra] and
#' [design_arithmetic]), computed in the database:
#' * `crossprod(D)`: `X'X`; `crossprod(D, "y")`: `X'y` for one or more
#'   columns of the remote data (or a `DBI::SQL()` expression).
#' * Scaled rows: `"w" * D` (a column name), `DBI::SQL("1 / v") * D` (an
#'   expression over the columns), `2 * D`, `D / "w"` and `-D` scale each
#'   row, as `w * D` with a numeric vector in memory. So `crossprod(D, "w" *
#'   D)` is `X'WX`, `crossprod("w" * D, "y")` is `X'Wy` and `colSums("w" * D,
#'   group = "g")` the weighted sums within groups. Only one design is
#'   involved: `crossprod(x, y)` needs `y` to be `x` with other scales.
#' * `colSums(D)`, `colMeans(D)`, and with `group = "g"` (a column, or a
#'   `DBI::SQL()` expression), the sums or means within groups (sorted as
#'   [base::rowsum()] sorts them).
#' * `D %*% beta`: the linear predictor, with all the columns of the remote
#'   data (to join it back), written to a temporary view (SQL engine) or
#'   table (chunked engine) whose name is returned, as a string of class
#'   `dm_remote_table`; [as.data.frame()] reads it, and its `"connection"`
#'   attribute is the connection where it lives (for a Parquet source, a
#'   private one). To keep it, create a table from it:
#'   `DBI::dbExecute(con, sprintf("CREATE TABLE preds AS SELECT * FROM %s", p))`.
#' * `as.matrix(D)` collects the rows and builds the dense matrix, only when
#'   there are at most `getOption("design.matrix.remote_max_rows", 1e5)` rows
#'   or with `force = TRUE`.
#' * `design_matrix(D, newdata)` applies the design to a data frame, and
#'   `design_matrix(D, con, table = ...)` to other remote data.
#' * `summary(D)`: the engine and, for each term, whether SQL computes it.
#'
#' @section Compression:
#' When the rows have few distinct *design patterns* (the same values in
#' every input column of the design, as with factors and discrete numeric
#' variables), `crossprod()`, `colSums()`, `colMeans()` and [design_lm()]
#' run one `GROUP BY` over those columns, which returns
#' each pattern with its number of rows and the sums of the weights and
#' responses, and compute the result from the design of the patterns, in
#' memory. This works with both engines (a variable evaluated in R is
#' evaluated on the patterns). The option `design.matrix.remote_compress`
#' controls it. With `"auto"` (default), a design with a continuous
#' numeric column (a `DOUBLE` column with many distinct values in a sample
#' of the rows) is not compressed, at almost no cost; otherwise the number
#' of patterns is estimated over all the rows (`approx_count_distinct()`),
#' and the products are compressed when it is at most
#' `getOption("design.matrix.remote_compress_ratio", 0.1)` times the number
#' of rows and at most `getOption("design.matrix.remote_compress_max", 1e6)`.
#' `TRUE` always compresses (unless every row is distinct); `FALSE` never
#' does. Designs that keep rows with missing values
#' (`na.action = na.pass`) are not compressed, nor are products with a
#' column or row scale that has missing values.
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
#'     beta <- solve(crossprod(D), crossprod(D, "y"))
#'     p <- D %*% beta
#'     head(as.data.frame(p))
#'     summary(D)
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
## `weights`, `offset`, `subset` (see dm_resolve_arg(), resolved in `env`)
## and `drop_unused` are for design_lm(): the rows are those of `subset`;
## rows with a missing weight or offset are dropped; factor levels absent
## from the rows used are dropped, as lm() builds its model frame.
## `source` can also be an opened source (dm_remote_source()); `override`
## gives SQL for model variables whose values were computed in R
## (dm_materialize()).
dm_remote <- function(object, source, contrasts.arg, xlev, dots, na.action = stats::na.pass,
                      weights = NULL, drop_unused = FALSE, subset = NULL, offset = NULL,
                      env = parent.frame(), override = list()) {
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
    src <- if (inherits(source, "dm_src")) source else dm_remote_source(source, ra$table, ra$sql)
    con <- src$holder$con
    if (!is.null(subset)) {
        sub <- dm_resolve_arg(subset, env, src$types, con, "subset", bool = TRUE)
        src$select <- sprintf("SELECT s.* FROM (%s) AS s WHERE %s", src$select, sub$sql)
        src$base_table <- NULL
        src$label <- paste(src$label, "(subset)")
    }
    extra <- list()
    w <- dm_resolve_arg(weights, env, src$types, con, "weights")
    if (!is.null(w)) extra[["(weights)"]] <- w
    o <- dm_resolve_arg(offset, env, src$types, con, "offset")
    if (!is.null(o)) extra[["(offset)"]] <- o

    if (inherits(object, "design_formula")) {
        contrasts.arg <- object$contrasts
        xlev <- object$xlevels
        object <- object$formula
    }
    if (inherits(object, "formula")) {
        if (length(rest)) stop("Unknown arguments: ", paste(names(rest), collapse = ", "), call. = FALSE)
        D <- dm_remote_formula(object, src, contrasts.arg, xlev, extra, drop_unused, override)
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
    structure(list(holder = holder, select = select, base_table = base_table, label = label,
         proto = proto, types = types), class = "dm_src")
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
    ## offset(e): a numeric variable in no term
    e <- if (is.call(expr) && identical(expr[[1L]], as.name("offset")) && length(expr) == 2L) expr[[2L]] else expr
    tr <- dm_sql_num(e, types, dm_con(D))
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
dm_remote_formula <- function(formula, src, contrasts.arg, xlev, extra = list(), drop_unused = FALSE,
                              override = list()) {
    D <- list(src = src, na = "omit")
    tt <- stats::terms(formula, data = src$proto)
    tv <- dm_terms_vars(tt)
    types <- src$types
    D$vars <- Map(function(nm, e) dm_remote_var(nm, e, types, D, NULL, TRUE), tv$names, tv$exprs)
    names(D$vars) <- tv$names
    ## variables computed in R into columns of the source
    for (nm in names(override)) {
        D$vars[[nm]]$kind <- "numeric"
        D$vars[[nm]]$sql <- override[[nm]]$sql
        D$vars[[nm]]$reason <- NULL
        D$vars[[nm]]$inputs <- unique(c(D$vars[[nm]]$inputs, override[[nm]]$inputs))
    }
    ## weights and offset argument: variables in no term, which select rows
    ## (and feed the fit)
    for (nm in names(extra)) {
        D$vars[[nm]] <- list(name = nm, expr = NULL, inputs = extra[[nm]]$inputs, kind = "numeric",
                             reason = NULL, col = NULL, levels = NULL, sql = extra[[nm]]$sql)
    }
    D$filter <- dm_remote_filter(D)

    ## a sample of the rows used, to see what the R expressions return
    sample <- dm_remote_sample(D)
    env <- environment(tt)
    predvars <- tv$exprs
    xlevels <- list()
    for (j in seq_along(D$vars)) {
        v <- D$vars[[j]]
        if (v$kind == "factor") {
            x <- if (types[[v$col]] == "factor" && !drop_unused) src$proto[[v$col]]
                 else dm_remote_distinct(D, v$col)[[1L]]
            x <- as.factor(x)
            if (drop_unused) x <- droplevels(x)
            xlevels[[v$name]] <- levels(x)
        } else if (v$kind == "other") {
            val <- tryCatch(eval(v$expr, sample, env), error = function(e) e)
            if (is.factor(val) || is.character(val)) {
                ## levels from the distinct values of the inputs
                inp <- dm_remote_distinct(D, v$inputs)
                fx <- as.factor(eval(v$expr, inp, env))
                lv <- levels(if (drop_unused) droplevels(fx) else fx)
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
    for (j in seq_along(predvars)) D$vars[[j]]$expr <- predvars[[j]]
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
    ## identifies the design through copies with scaled rows
    D$id <- dm_tmp_name("dm")

    ## which terms SQL can compute
    D$terms <- dm_remote_terms(D)
    bad <- Filter(function(t) !t$sql_ok, D$terms)
    resp_bad <- !is.null(D$response) && D$vars[[D$response]]$kind == "other"
    ## variables in no term (an offset, say) decide which rows are used
    used <- unique(unlist(lapply(D$terms, function(t) t$vars)))
    ## (offset() terms are not in the design: a fit or a prediction
    ## evaluates them, in R if need be)
    offs <- dm_terms_vars(D$plan$terms)$names[attr(D$plan$terms, "offset")]
    loose <- setdiff(names(D$vars)[vapply(D$vars, function(v) v$kind == "other", NA)],
                     c(names(D$vars)[used], D$response, offs))
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
    if (!is.null(x$scale)) cat("rows scaled by: ", x$scale, "\n", sep = "")
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

## A `weights` or `group` argument: a column name, or a SQL expression
## (`DBI::SQL()`) over the columns of the remote data. An expression is
## computed as an extra column of the source, so that both engines (and
## the handling of missing values) treat it as any other column.
## Returns list(D = the design to use, col = the column name, or NULL).
dm_remote_arg <- function(D, x, what) {
    if (is.null(x)) return(list(D = D, col = NULL))
    if (!inherits(x, "SQL")) return(list(D = D, col = dm_remote_column(D, x, what)))
    checkmate::assert_string(as.character(x), .var.name = what)
    col <- paste0("__dm_arg_", what)
    list(D = dm_remote_extend(D, stats::setNames(as.character(x), col), what), col = col)
}

## The design over the same rows, with extra columns `exprs` (named SQL
## expressions over the columns of the source) added to the source.
dm_remote_extend <- function(D, exprs, what = "expression") {
    con <- dm_con(D)
    items <- sprintf("(%s) AS %s", exprs, dm_qi(D, names(exprs)))
    select <- sprintf("SELECT s.*, %s FROM (%s) AS s", paste(items, collapse = ", "), D$src$select)
    proto <- tryCatch(DBI::dbGetQuery(con, sprintf("SELECT * FROM (%s) AS s LIMIT 0", select)),
                      error = function(e) stop(sprintf("`%s`: invalid SQL expression: %s", what,
                                                       conditionMessage(e)), call. = FALSE))
    D$src$select <- select
    D$src$proto <- proto
    D$src$types <- vapply(proto, dm_col_type, "")
    ## the extra columns are not in the base table: batches read a copy
    D$src$base_table <- NULL
    D
}

## Residuals of a linear fit: the remote design over the same rows, with two
## more columns in its remote data, `.fitted` (X beta) and `.resid` (the
## response minus `.fitted`), usable as any column, for instance to scale
## the design (`".resid" * R`) for sandwich covariance matrices. The
## predictions live in a temporary view (SQL engine) or table (chunked).
dm_remote_residuals <- function(x, beta, y = NULL, z = NULL, off = NULL) {
    checkmate::assert_class(x, "dm_remote")
    checkmate::assert_numeric(beta, len = x$p, any.missing = FALSE)
    x <- dm_remote_unscaled(x)
    if (is.null(z)) {
        yy <- dm_remote_y(x, y)
        if (is.null(yy)) stop("Give `y`: the design has no response.", call. = FALSE)
        if (length(yy$cols) != 1L) stop("`y` must be a single column.", call. = FALSE)
        z <- if (yy$response) x$vars[[yy$cols]]$sql else sprintf("CAST(s.%s AS DOUBLE)", dm_qi(x, yy$cols))
        if (is.null(z))
            stop("The response has no SQL translation; give `y`, a column of the remote data.", call. = FALSE)
    }
    taken <- intersect(c(".fitted", ".resid"), names(x$src$types))
    if (length(taken))
        stop(sprintf("The remote data already has a column `%s`.", taken[1L]), call. = FALSE)
    view <- dm_remote_predict(x, beta, column = ".fitted")
    out <- x
    ## (z and the offset are expressions over the columns of the data, which
    ## the view keeps; .fitted there is X beta)
    out$src$select <- sprintf("SELECT s.* EXCLUDE (\".fitted\"), %s AS \".fitted\", (%s) - s.\".fitted\" AS \".resid\" FROM %s AS s",
                              if (is.null(off)) "s.\".fitted\"" else sprintf("s.\".fitted\" + (%s)", off),
                              z, dm_qi(x, unclass(view)[1L]))
    out$src$base_table <- NULL
    proto <- DBI::dbGetQuery(dm_con(x), sprintf("SELECT * FROM (%s) AS s LIMIT 0", out$src$select))
    out$src$proto <- proto
    out$src$types <- vapply(proto, dm_col_type, "")
    out$src$view <- view
    out
}

## The response column of a formula design.
dm_remote_y <- function(D, y) {
    if (!is.null(y)) return(list(cols = dm_remote_column(D, y, "y", several = TRUE), response = FALSE))
    if (is.null(D$response)) return(NULL)
    list(cols = D$response, response = TRUE)
}

## X'WX, X'Wy, y'Wy, n and the sum of the weights of a remote design, in
## one pass where possible (compressed, SQL or chunked). `weights`: a column
## name, a DBI::SQL() expression or NULL; `y`: columns, or NULL for the
## response of the formula.
dm_remote_moments <- function(x, y = NULL, weights = NULL) {
    yy <- dm_remote_y(x, y)
    wa <- dm_remote_arg(x, weights, "weights")
    x <- wa$D
    weights <- wa$col
    out <- dm_compressed_moments(x, yy, weights)
    if (is.null(out))
        out <- if (x$engine == "sql") dm_sql_moments(x, yy, weights) else dm_chunked_moments(x, yy, weights)
    x$cache$n <- out$n
    out
}

## X %*% beta written to a temporary view (SQL engine) or table (chunked
## engine), or to the table `name`, with the columns `keep` (default: all)
## of the remote data. Returns the name, of class dm_remote_table.
dm_remote_predict <- function(x, beta, name = NULL, column = ".fitted", keep = NULL, overwrite = FALSE) {
    checkmate::assert_class(x, "dm_remote")
    p <- x$p
    B <- if (is.matrix(beta)) beta else matrix(beta, ncol = 1L)
    checkmate::assert_matrix(B, mode = "numeric", nrows = p, any.missing = FALSE, .var.name = "beta")
    checkmate::assert_true(all(is.finite(B)), .var.name = "all(is.finite(beta))")
    storage.mode(B) <- "double"
    sc <- x$scale
    x <- dm_remote_unscaled(x)
    keep <- if (is.null(keep)) names(x$src$types) else dm_remote_column(x, keep, "keep", several = TRUE)
    outcols <- if (ncol(B) == 1L) column else colnames(B) %||% paste0(column, seq_len(ncol(B)))
    if (any(outcols %in% keep))
        stop("The prediction column clashes with a column of the remote data.", call. = FALSE)
    if (!is.null(sc)) {
        ## rows scaled: the predictions times the scale
        x <- dm_remote_extend(x, c(`__dm_scale` = sc), "scale")
        tmp <- if (x$engine == "sql") dm_sql_predict(x, B, NULL, outcols, c(keep, "__dm_scale"), FALSE)
               else dm_chunked_predict(x, B, NULL, outcols, c(keep, "__dm_scale"), FALSE)
        sel <- sprintf("SELECT %s FROM %s", paste(c(dm_qi(x, keep),
                       sprintf("%s * \"__dm_scale\" AS %s", dm_qi(x, outcols), dm_qi(x, outcols))), collapse = ", "),
                       dm_qi(x, tmp))
        if (is.null(name)) {
            name <- dm_tmp_name("dm_pred")
            dm_execute(x, sprintf("CREATE TEMPORARY VIEW %s AS %s", dm_qi(x, name), sel))
        } else {
            dm_execute(x, sprintf("CREATE %sTABLE %s AS %s", if (overwrite) "OR REPLACE " else "", dm_qi(x, name), sel))
        }
        out <- name
    } else {
        out <- if (x$engine == "sql") dm_sql_predict(x, B, name, outcols, keep, overwrite)
               else dm_chunked_predict(x, B, name, outcols, keep, overwrite)
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
    dm_remote_predict(x, y)
}

## Scaled rows, as `w * D` for a design in memory: the scale is a SQL
## expression over the columns of the remote data (a number, a column
## name, or a DBI::SQL() expression), kept unevaluated in `x$scale`.
dm_remote_scale_sql <- function(D, s, op) {
    if (inherits(s, "SQL")) {
        checkmate::assert_string(as.character(s), .var.name = "scale")
        e <- sprintf("(%s)", as.character(s))
    } else if (is.character(s)) {
        dm_remote_column(D, s, "scale")
        if (!D$src$types[[s]] %in% c("double", "integer", "logical"))
            stop(sprintf("Column `%s` is not numeric.", s), call. = FALSE)
        e <- sprintf("CAST(%s AS DOUBLE)", dm_qi(D, s))
    } else if (is.numeric(s) && length(s) == 1L && is.null(dim(s))) {
        e <- dm_sql_double(s)
    } else {
        stop("A remote design can only be scaled by a number, the name of a column or a DBI::SQL() expression.",
             call. = FALSE)
    }
    if (op == "/") e <- sprintf("(CAST(1 AS DOUBLE) / %s)", e)
    out <- if (is.null(D$scale)) e else sprintf("%s * %s", D$scale, e)
    ## check the expression now
    dm_remote_extend(D, c(`__dm_scale` = out), "scale")
    out
}

#' @export
Ops.dm_remote <- function(e1, e2) {
    op <- .Generic
    if (nargs() == 1L) {
        if (op == "+") return(e1)
        if (op == "-") {
            e1$scale <- dm_remote_scale_sql(e1, -1, "*")
            return(e1)
        }
    } else if (op %in% c("*", "/")) {
        d1 <- inherits(e1, "dm_remote")
        d2 <- inherits(e2, "dm_remote")
        if (d1 && d2) stop(sprintf("`%s` between two remote designs is not supported.", op), call. = FALSE)
        if (op == "/" && d2) stop("Only `design / scale` is supported, with the design first.", call. = FALSE)
        D <- if (d1) e1 else e2
        D$scale <- dm_remote_scale_sql(D, if (d1) e2 else e1, op)
        return(D)
    }
    stop(sprintf("`%s` is not supported for remote designs.", op), call. = FALSE)
}

## The design without its row scale.
dm_remote_unscaled <- function(D) {
    D$scale <- NULL
    D
}

## A remote design and its scale as a weights argument (DBI::SQL() or NULL).
dm_remote_weights <- function(D, extra = NULL) {
    w <- c(D$scale, extra)
    list(D = dm_remote_unscaled(D),
         w = if (length(w)) DBI::SQL(paste(w, collapse = " * ")) else NULL)
}

## Is `b` the same remote design as `a` (up to the row scale)?
dm_remote_same <- function(a, b) identical(a$id, b$id)

#' @export
as.matrix.dm_remote <- function(x, force = FALSE, ...) {
    max_rows <- getOption("design.matrix.remote_max_rows", 1e5)
    if (!isTRUE(force) && nrow(x) > max_rows)
        stop(sprintf("The remote design has %s rows (more than %s): use `force = TRUE` to build the dense matrix in memory.",
                     dm_fmt_count(nrow(x)), dm_fmt_count(max_rows)), call. = FALSE)
    sc <- x$scale
    x <- dm_remote_unscaled(x)
    cols <- dm_remote_inputs(x)
    if (!is.null(sc)) {
        x <- dm_remote_extend(x, c(`__dm_scale` = sc), "scale")
        cols <- c(cols, "__dm_scale")
    }
    chunk <- dm_query(x, dm_chunk_sql(x, cols))
    cd <- dm_chunk_design(x, chunk)
    X <- as.matrix(cd$design)
    if (!is.null(sc)) {
        a <- attributes(X)
        X <- chunk$`__dm_scale`[cd$keep] * X
        attributes(X) <- a
    }
    X
}

#' @method crossprod dm_remote
#' @export
crossprod.dm_remote <- function(x, y = NULL, ...) {
    if (!inherits(x, "dm_remote"))
        stop("Only `crossprod(design)`, `crossprod(design, design)` and `crossprod(design, \"column\")` are supported, with a design first.",
             call. = FALSE)
    ## crossprod() is a primitive: its dispatch drops argument names, so
    ## that `crossprod(D, weights = "w")` would arrive as `y = "w"`
    nm <- names(sys.call())
    if (...length() || (!is.null(nm) && any(!nm %in% c("", "x", "y"))))
        stop("`crossprod()` of a remote design takes only `x` and `y`; for weights, scale the rows: `crossprod(D, \"w\" * D)`.",
             call. = FALSE)
    if (is.null(y) || inherits(y, "dm_remote")) {
        ## X' diag(s_x s_y) X
        if (!is.null(y) && !dm_remote_same(x, y))
            stop("`crossprod(x, y)` of remote designs needs the same design, possibly with scaled rows.", call. = FALSE)
        sy <- if (is.null(y)) x$scale else y$scale
        dw <- dm_remote_weights(x, sy)
        if (is.null(dw$w)) {
            D <- dw$D
            out <- dm_compressed_moments(D, NULL, NULL)
            if (!is.null(out)) return(out$xtx)
            return(if (D$engine == "sql") dm_sql_gram(D) else dm_chunked_moments(D, NULL, NULL)$xtx)
        }
        return(dm_remote_moments(dw$D, NULL, dw$w)$xtx)
    }
    if (inherits(y, "SQL")) {
        D <- dm_remote_extend(dm_remote_unscaled(x), c(`__dm_y` = as.character(y)), "y")
        D$scale <- x$scale
        out <- crossprod(D, "__dm_y")
        colnames(out) <- as.character(y)
        return(out)
    }
    if (!is.character(y))
        stop("For a remote design, `y` must be another remote design, column names or a DBI::SQL() expression.",
             call. = FALSE)
    dm_remote_column(x, y, "y", several = TRUE)
    dw <- dm_remote_weights(x)
    if (!is.null(dw$w)) {
        wa <- dm_remote_arg(dw$D, dw$w, "weights")
        D <- wa$D
        yy <- list(cols = y, response = FALSE)
        out <- dm_compressed_moments(D, yy, wa$col, gram = FALSE)
        if (is.null(out))
            out <- if (D$engine == "sql") list(xty = dm_sql_xty(D, yy, wa$col)) else dm_chunked_moments(D, yy, wa$col, gram = FALSE)
        return(out$xty)
    }
    D <- dw$D
    yy <- list(cols = y, response = FALSE)
    out <- dm_compressed_moments(D, yy, NULL, gram = FALSE)
    if (!is.null(out)) return(out$xty)
    if (D$engine == "sql") dm_sql_xty(D, yy) else dm_chunked_moments(D, yy, NULL, gram = FALSE)$xty
}

dm_remote_colsums <- function(x, group, reorder, mean) {
    if (!isTRUE(reorder))
        stop("Remote designs support only `reorder = TRUE`.", call. = FALSE)
    dw <- dm_remote_weights(x)
    wa <- dm_remote_arg(dw$D, dw$w, "weights")
    ga <- dm_remote_arg(wa$D, group, "group")
    x <- ga$D
    weights <- wa$col
    group <- ga$col
    out <- dm_compressed_colsums(x, group, weights, mean)
    if (!is.null(out)) return(out)
    if (x$engine == "sql") dm_sql_colsums(x, group, weights, mean)
    else dm_chunked_colsums(x, group, weights, mean)
}

setOldClass("dm_remote")

#' @rdname design_remote
#' @param x,object A remote design.
#' @param na.rm,dims Only the defaults are supported.
#' @param ... Not used.
#' @param group For `colSums()` and `colMeans()`, the name of a column of the
#'   remote data whose values define the groups, or a SQL expression over
#'   its columns (`DBI::SQL()`).
#' @param reorder Only `TRUE` (groups sorted) is supported.
#' @export
setMethod("colSums", "dm_remote", function(x, na.rm = FALSE, dims = 1, ..., group = NULL, reorder = TRUE) {
    dm_defaults_only(na.rm, dims)
    dm_remote_colsums(x, group, reorder, mean = FALSE)
})

#' @rdname design_remote
#' @export
setMethod("colMeans", "dm_remote", function(x, na.rm = FALSE, dims = 1, ..., group = NULL, reorder = TRUE) {
    dm_defaults_only(na.rm, dims)
    dm_remote_colsums(x, group, reorder, mean = TRUE)
})

#' @rdname design_remote
#' @param sql For `summary()`: also print the generated SQL (the query that
#'   computes the model variables and, for the SQL engine, the aggregates of
#'   `crossprod()`; for the chunked engine, the query that reads a batch).
#' @export
summary.dm_remote <- function(object, sql = FALSE, ...) {
    x <- object
    checkmate::assert_flag(sql)
    tab <- data.frame(
        term = vapply(x$terms, function(t) t$label, ""),
        sql = vapply(x$terms, function(t) t$sql_ok, NA),
        reason = vapply(x$terms, function(t) t$reason %||% "", ""),
        stringsAsFactors = FALSE)
    attr(tab, "engine") <- if (x$engine == "sql") "sql" else sprintf("chunked (%s)", x$why_short)
    attr(tab, "scale") <- x$scale
    if (sql) {
        attr(tab, "sql") <- if (x$engine == "sql") {
            c(paste0("-- model variables\n", dm_sql_base(x)),
              paste0("-- crossprod\n", vapply(dm_sql_gram_queries(x, NULL), function(q) q$sql, "")))
        } else {
            paste0("-- batch\n", dm_chunk_sql(x, dm_remote_inputs(x), range = c(0, x$chunk_size)))
        }
    }
    class(tab) <- c("summary.dm_remote", class(tab))
    tab
}

#' @export
print.summary.dm_remote <- function(x, ...) {
    cat("engine: ", attr(x, "engine"), "\n", sep = "")
    if (!is.null(attr(x, "scale"))) cat("rows scaled by: ", attr(x, "scale"), "\n", sep = "")
    sql <- attr(x, "sql")
    tab <- x
    attr(tab, "engine") <- attr(tab, "sql") <- attr(tab, "scale") <- NULL
    class(tab) <- "data.frame"
    print(tab, row.names = FALSE)
    for (q in sql) cat("\n", q, "\n", sep = "")
    invisible(x)
}

## The public column type (Kata design.matrix#bvxy): the remote column of the
## fitting core (primitives.R, class design_column, made by dm_column()) is
## also what users get from a remote design: design_column(D, "y"), the
## products D %*% b, rowSums(D), rowMeans(D) and quadratic_diagonal(D, V)
## return it, arithmetic composes it lazily, and collect(), as.data.frame()
## and compute() evaluate or store it. In memory a column is a plain numeric
## vector, and the same calls work on it (collect() and compute() return it).
## Frames, their depth and lineage stay internal.

#' Columns of a design
#'
#' A *column* holds one value per row of a design. In memory it is a plain
#' numeric vector. For a remote design (see [design_remote]), it is a
#' `design_column`: a lazy SQL expression over the rows of the design, which
#' nothing evaluates until its values are asked for. Products with a remote
#' design return one, arithmetic on it builds a new one, and [collect()],
#' [as.data.frame()] or [compute()] evaluate it. So the same code runs on a
#' data frame and on a DuckDB table:
#'
#' ```
#' eta <- D %*% beta                 # vector in memory, design_column remotely
#' mu  <- exp(eta)
#' r   <- design_column(D, "y") - mu
#' collect(r)                        # the values, in R
#' ```
#'
#' @section Making columns:
#' * `design_column(D, "y")`: the column `y` of the remote data, as a double
#'   (numeric and logical columns; a logical column gives 0 and 1).
#'   `design_column(D, DBI::SQL("y / w"))`: a SQL expression over the
#'   columns of the remote data, kept as SQL computes it (a logical
#'   expression stays logical, and can select rows: `D[col, ]`). The column
#'   is over the rows of `D` (its filter, for a design with
#'   `na.action = na.omit` or a subset), not its row scale.
#' * In memory, `design_column(D, "y")` returns the variable `y` of the
#'   model frame of the design (the response and the variables of the
#'   formula, by their names in the model frame, such as `"log(z)"`) as a
#'   double vector. A design in memory does not keep the other columns of
#'   its data, so any other name is an error: use the vector itself. A
#'   `DBI::SQL()` expression is an error in memory (there is no database to
#'   evaluate it).
#' * `D %*% b`, `rowSums(D)`, `rowMeans(D)` and `quadratic_diagonal(D, V)`
#'   of a remote design return a column (`b` with `k` columns gives a column
#'   with `k` columns).
#' * Arithmetic: `+`, `-`, `*`, `/`, `^` with numbers and other columns of
#'   the same design, comparisons, `&`, `|`, `!`, [is.na()], [is.nan()],
#'   [is.finite()], [is.infinite()] and the functions of the `Math` group
#'   (`exp()`, `log()` with one argument, `sqrt()`, `abs()`, ...), with R's
#'   missing values (SQL `NULL` is `NA`, and `NaN` stays `NaN`).
#'   `col_ifelse()`, `col_pmin()` and `col_pmax()` combine columns, and
#'   `sum()`, `mean()`, `min()`, `range()`, `any()`, ... and `col_sums()`
#'   summarize them, in one query: see [column_arithmetic].
#' * `e * D` (or `D / w`), with `e` a column of a remote design `D`: the
#'   design with rows scaled by the column, as for a vector in memory (for
#'   instance `colSums(e * D, group = "cl")`, or `crossprod(D, w * D)`).
#'
#' @section Getting values:
#' * `collect(x)`: the values in R, a double vector (a matrix for several
#'   columns). Remotely this reads every row: only for results that fit in
#'   memory. In memory it returns `x`.
#' * `as.data.frame(x, keep = "id")`: a data frame with one column per
#'   column of `x` (named `.fitted` for `D %*% b`, `.rowsum`, `.rowmean`,
#'   `.quad`, the column name for `design_column(D, "y")`, else `.value`;
#'   or `nm`), after the columns `keep` of the remote data, read in the same
#'   query. Remote rows have no order: to match the values with other data,
#'   keep a key column. Two calls to `collect()` need not return the rows
#'   in the same order.
#' * `compute(x)`: the column computed once into a temporary table, and a
#'   column that reads it (later queries scan the table instead of
#'   recomputing the expression and the views below it). In memory it
#'   returns `x`.
#' * `x[, j]`: the columns `j` of a column with several columns (as for the
#'   matrix in memory), for instance of `D %*% cbind(b1, b2)`. Rows cannot
#'   be selected: remote rows have no order.
#' * `print()` and `format()` show the number of columns and rows, the SQL
#'   and the connection, without reading any row. `length(x)` is the number
#'   of columns of `x` (not of rows: remotely, use `nrow(D)`).
#'
#' @section Lifetime:
#' A remote column is valid while its connection is open. The views and
#' temporary tables it reads (from products, and from `compute()`) live in
#' the session of that connection, and disappear when it is closed; for a
#' Parquet source, the column keeps its private connection open. To keep
#' values beyond the session, read them with `as.data.frame()` and write
#' them where they belong. Columns of one design combine with each other,
#' including the results of separate products (`D %*% b1 + D %*% b2`, or
#' `(D %*% b) / sqrt(quadratic_diagonal(D, V))`): a product becomes a view
#' only when its values are first needed, on top of the other products it
#' is combined with. Columns of different designs do not combine. Nor do
#' two products that were each read (by `collect()`, `as.data.frame()`,
#' `compute()` or a sum) before being combined, since their views then
#' stand side by side and the database has no key to align their rows:
#' combine the columns first, compute them in one product (`D %*% cbind(b1,
#' b2)`, then `[, 1]` and `[, 2]`), or read them with `as.data.frame(keep =
#' "id")` and match them in R. A column made by `design_column()` combines
#' with any product of its design.
#'
#' `collect()` and `compute()` are generics of this package. If another
#' package that defines functions of the same names is attached after
#' design.matrix, call `design.matrix::collect()` and
#' `design.matrix::compute()`.
#'
#' @param x For `design_column()`, a design: a `design_matrix` or a remote
#'   design. For the other functions, a column.
#' @param column The name of a column (a string), or for a remote design a
#'   SQL expression over its columns (`DBI::SQL()`).
#' @param row.names,optional Not used.
#' @param keep Names of columns of the remote data to read with the values
#'   (for instance a key, to match the rows with other data).
#' @param nm Names of the value columns of the data frame (one per column of
#'   `x`).
#' @param ... Not used.
#' @return `design_column()`: a double vector in memory, a `design_column`
#'   remotely. `collect()`: a double vector or matrix. `as.data.frame()`: a
#'   data frame. `compute()`: a `design_column` (in memory, `x`).
#' @examples
#' d <- data.frame(y = rnorm(100), x = rnorm(100), f = gl(4, 25))
#' D <- design_matrix(y ~ x * f, d)
#' b <- solve(crossprod(D), crossprod(D, d$y))
#' r <- design_column(D, "y") - D %*% b
#' head(collect(r))
#'
#' if (requireNamespace("duckdb", quietly = TRUE)) {
#'     con <- DBI::dbConnect(duckdb::duckdb())
#'     d$id <- seq_len(nrow(d))
#'     DBI::dbWriteTable(con, "d", d)
#'     Dr <- design_matrix(y ~ x * f, con, table = "d")
#'     r <- design_column(Dr, "y") - Dr %*% b
#'     r
#'     head(as.data.frame(r, keep = "id", nm = "resid"))
#'     r2 <- compute(r^2)
#'     sum(collect(r2))
#'     DBI::dbDisconnect(con)
#' }
#' @name design_column
#' @export
design_column <- function(x, column) {
    if (inherits(x, "dm_remote")) {
        D <- dm_remote_unscaled(x)
        if (inherits(column, "SQL")) {
            cc <- dm_remote_cols(D, column, "column")
            return(dm_column(D, cc$sql))
        }
        checkmate::assert_string(column, .var.name = "column")
        dm_numeric_cols(D, dm_remote_column(D, column, "column"))
        return(dm_column_label(dm_column(D, sprintf("CAST(s.%s AS DOUBLE)", dm_qi(D, column))), column))
    }
    if (!inherits(x, "design_matrix"))
        stop("`x` must be a design: a `design_matrix` or a remote design.", call. = FALSE)
    if (inherits(column, "SQL"))
        stop("A DBI::SQL() expression makes a column only of a remote design; in memory, use the vector itself.",
             call. = FALSE)
    checkmate::assert_string(column, .var.name = "column")
    vars <- (if (dm_is_comb(x)) x$design else x)$vars
    if (is.null(vars) || !column %in% names(vars))
        stop(sprintf(paste("`%s` is not a variable of the design: in memory, a design keeps only the",
                           "variables of its model frame; use the vector itself."), column), call. = FALSE)
    v <- vars[[column]]
    if (!is.null(dim(v)) || !(is.numeric(v) || is.logical(v)))
        stop(sprintf("Variable `%s` is not a numeric vector.", column), call. = FALSE)
    as.double(v)
}

## The names of the values of a column (for as.data.frame()), kept by the
## public calls that make it; arithmetic drops them.
dm_column_label <- function(x, label) {
    y <- unclass(x)
    y$label <- label
    structure(y, class = "design_column")
}

## Stop when the connection of a remote column is closed.
dm_column_valid <- function(x) {
    con <- dm_con(unclass(x)$frame)
    if (!DBI::dbIsValid(con))
        stop("The connection of this remote column is closed: a column is valid only while its connection is open.",
             call. = FALSE)
    invisible(con)
}

#' @rdname design_column
#' @export
collect <- function(x, ...) UseMethod("collect")

#' @rdname design_column
#' @export
collect.default <- function(x, ...) x

#' @rdname design_column
#' @export
collect.design_column <- function(x, ...) {
    dm_column_valid(x)
    out <- dmp_collect(x)
    if (is.matrix(out)) dimnames(out) <- NULL
    out
}

#' @rdname design_column
#' @export
compute <- function(x, ...) UseMethod("compute")

#' @rdname design_column
#' @export
compute.default <- function(x, ...) x

#' @rdname design_column
#' @export
compute.design_column <- function(x, ...) {
    dm_column_valid(x)
    out <- dmc_freeze(x, table = TRUE)
    lab <- unclass(x)$label
    if (is.null(lab)) out else dm_column_label(out, lab)
}

#' @rdname design_column
#' @export
as.data.frame.design_column <- function(x, row.names = NULL, optional = FALSE, ..., keep = NULL, nm = NULL) {
    dm_column_valid(x)
    fr <- dm_frame_bind(unclass(x)$frame)
    sql <- unclass(x)$sql
    k <- length(sql)
    nm <- nm %||% unclass(x)$label %||% (if (k == 1L) ".value" else paste0(".value", seq_len(k)))
    checkmate::assert_character(nm, len = k, any.missing = FALSE, unique = TRUE, .var.name = "nm")
    if (!is.null(keep)) {
        dm_remote_column(fr, keep, "keep", several = TRUE)
        if (any(startsWith(keep, "__dm_")))
            stop("`keep`: only columns of the remote data can be kept.", call. = FALSE)
        if (any(nm %in% keep))
            stop(sprintf("`keep`: column `%s` has the name of the values; rename them with `nm`.",
                         intersect(nm, keep)[1L]), call. = FALSE)
    }
    items <- c(if (length(keep)) paste0("s.", dm_qi(fr, keep)),
               sprintf("CAST(%s AS DOUBLE) AS %s", sql, dm_qi(fr, nm)))
    dm_query(fr, sprintf("SELECT %s FROM (%s) AS s WHERE %s", paste(items, collapse = ", "),
                         fr$src$select, fr$filter))
}

#' @rdname design_column
#' @param i Not supported for a remote column (remote rows have no order).
#' @param j Indices (numeric or logical) of columns of `x`.
#' @param drop Not used.
#' @export
`[.design_column` <- function(x, i, j, ..., drop = TRUE) {
    if (!missing(i))
        stop("Remote rows have no order: a remote column selects only columns, with `x[, j]`; ",
             "select rows of the design instead, with `D[i, ]`.", call. = FALSE)
    if (missing(j)) return(x)
    if (!is.numeric(j) && !is.logical(j))
        stop("Select the columns of a remote column by number or with a logical vector.", call. = FALSE)
    idx <- seq_len(length(x))[j]
    if (!length(idx) || anyNA(idx)) stop("Column index out of range.", call. = FALSE)
    u <- unclass(x)
    out <- dm_column(u$frame, u$sql[idx], notnull = u$notnull[idx])
    if (is.null(u$label)) out else dm_column_label(out, u$label[idx])
}

#' @rdname design_column
#' @export
format.design_column <- function(x, ...) {
    fr <- unclass(x)$frame
    sql <- unclass(x)$sql
    con <- dm_con(fr)
    open <- tryCatch(DBI::dbIsValid(con), error = function(e) FALSE)
    n <- dm_nrows(fr)
    top <- sprintf("<design_column: %d column%s, %s rows, on <%s>%s>", length(sql),
                    if (length(sql) == 1L) "" else "s",
                    if (is.null(n)) "?" else dm_fmt_count(n), class(con)[1L],
                    if (open) "" else " (closed)")
    shown <- utils::head(sql, 3L)
    lab <- unclass(x)$label
    pre <- if (is.null(lab)) "" else paste0(utils::head(lab, 3L), ": ")
    body <- paste0("  ", pre, ifelse(nchar(shown) > 70L, paste0(substr(shown, 1L, 67L), "..."), shown))
    if (length(sql) > 3L) body <- c(body, sprintf("  ... (%d more)", length(sql) - 3L))
    c(top, body)
}

#' @rdname design_column
#' @export
print.design_column <- function(x, ...) {
    cat(format(x), sep = "\n")
    invisible(x)
}

## The public column type (Kata design.matrix#bvxy, epxj): the remote column
## of the fitting core (primitives.R, class design_column, made by
## dm_column()) is also what users get from a remote design:
## design_column(D, "y"), the products D %*% b, rowSums(D), rowMeans(D) and
## quadratic_diagonal(D, V) return it, arithmetic composes it lazily, and
## the base generics as.vector(), as.matrix() and as.data.frame() evaluate
## it; col_compute() stores it in a table (temporary, or permanent with a
## name). In memory a column is a plain numeric vector (or matrix), base R
## applies to it, and col_compute() returns it.
## Frames, their depth and lineage stay internal.

#' Columns of a design
#'
#' A *column* holds one value per row of a design. In memory it is a plain
#' numeric vector. For a remote design (see [design_remote]), it is a
#' `design_column`: a lazy SQL expression over the rows of the design, which
#' nothing evaluates until its values are asked for. Products with a remote
#' design return one, arithmetic on it builds a new one, and the base
#' generics [as.vector()], [as.matrix()] and [as.data.frame()] evaluate it,
#' as they would a vector in memory. So the same code runs on a data frame
#' and on a DuckDB table:
#'
#' ```
#' eta <- D %*% beta                 # vector in memory, design_column remotely
#' mu  <- exp(eta)
#' r   <- design_column(D, "y") - mu
#' as.vector(r)                      # the values, in R
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
#' * `as.vector(x)`: the values in R, a double vector, for a column with one
#'   column (a column with several columns is an error: use `as.matrix()`).
#'   `as.matrix(x)`: an `n x k` double matrix, without dimnames (as the
#'   products in memory), for one column or several. Remotely both read
#'   every row: only for results that fit in memory. In memory they are
#'   base R's functions on the vector or matrix.
#' * `as.data.frame(x, keep = "id")`: a data frame with one column per
#'   column of `x` (named `.fitted` for `D %*% b`, `.rowsum`, `.rowmean`,
#'   `.quad`, the column name for `design_column(D, "y")`, else `.value`;
#'   or `nm`), after the columns `keep` of the remote data, read in the same
#'   query. Remote rows have no order: to match the values with other data,
#'   keep a key column. Two calls to `as.vector()` need not return the rows
#'   in the same order.
#' * `col_compute(x)`: the column computed once into a temporary table, and
#'   a column that reads it (later queries scan the table instead of
#'   recomputing the expression and the views below it).
#'   `col_compute(x, name = "preds")`: the same, into a permanent table
#'   `preds` of the connection (`CREATE TABLE ... AS`; an error if the
#'   table exists, unless `overwrite = TRUE`, which replaces it), for use
#'   outside the package. The table holds the rows of the design (after its
#'   filter: `na.action`, subsets), each with every column of the frame the
#'   values were computed on: the columns of the data, the values (named
#'   `nm`: by default as `as.data.frame()` names them) and internal columns
#'   of the package (names starting with `__dm_`, the products and lookups
#'   under the column), which can be ignored. The column returned reads
#'   that table and stays aligned with the rows of its design, as with a
#'   temporary table: the rows keep, next to the values, the columns that
#'   later columns of the design read, so a product made later is a view
#'   over the table. To match the values with other data outside the
#'   package, join on a key column of the data (such as `id`). In memory
#'   `col_compute()` returns `x` (there is no table).
#' * `x[, j]`: the columns `j` of a column with several columns (as for the
#'   matrix in memory), for instance of `D %*% cbind(b1, b2)`. Rows cannot
#'   be selected: remote rows have no order.
#' * `print()` and `format()` show the number of columns and rows, the SQL
#'   and the connection, without reading any row. `length(x)` is the number
#'   of columns of `x` (not of rows: remotely, use `nrow(D)`).
#'
#' @section Lifetime:
#' A remote column is valid while its connection is open. The views and
#' temporary tables it reads (from products, and from `col_compute()`)
#' live in the session of that connection, and disappear when it is
#' closed; for a Parquet source, the column keeps its private connection
#' open. To keep values beyond the session, store them in a permanent table
#' with `col_compute(x, name = )` (on the connection of the design: for a
#' Parquet source, a private in-memory connection, so read them with
#' `as.data.frame()` instead), or read them with `as.data.frame()` and
#' write them where they belong. Columns of one design combine with each other,
#' including the results of separate products (`D %*% b1 + D %*% b2`, or
#' `(D %*% b) / sqrt(quadratic_diagonal(D, V))`): a product becomes a view
#' only when its values are first needed, on top of the other products it
#' is combined with. Columns of different designs do not combine. Nor do
#' two products that were each read (by `as.vector()`, `as.data.frame()`,
#' `col_compute()` or a sum) before being combined, since their views then
#' stand side by side and the database has no key to align their rows:
#' combine the columns first, compute them in one product (`D %*% cbind(b1,
#' b2)`, then `[, 1]` and `[, 2]`), or read them with `as.data.frame(keep =
#' "id")` and match them in R. A column made by `design_column()` combines
#' with any product of its design.
#'
#' @param x For `design_column()`, a design: a `design_matrix` or a remote
#'   design. For the other functions, a column.
#' @param column The name of a column (a string), or for a remote design a
#'   SQL expression over its columns (`DBI::SQL()`).
#' @param row.names,optional Not used.
#' @param keep Names of columns of the remote data to read with the values
#'   (for instance a key, to match the rows with other data).
#' @param nm Names of the value columns (one per column of `x`): of the data
#'   frame for `as.data.frame()`, of the table for `col_compute()` with a
#'   `name`.
#' @param name For `col_compute()`: `NULL` (default) for a temporary table,
#'   or the name of a permanent table to create on the connection (a
#'   string, quoted as an identifier).
#' @param overwrite For `col_compute()` with a `name`: replace a table of
#'   that name if it exists (otherwise that is an error).
#' @param mode As in [as.vector()].
#' @param ... Not used.
#' @return `design_column()`: a double vector in memory, a `design_column`
#'   remotely. `as.vector()`: a double vector. `as.matrix()`: a double
#'   matrix. `as.data.frame()`: a data frame. `col_compute()`: a
#'   `design_column` that reads the table (in memory, `x`).
#' @examples
#' d <- data.frame(y = rnorm(100), x = rnorm(100), f = gl(4, 25))
#' D <- design_matrix(y ~ x * f, d)
#' b <- solve(crossprod(D), crossprod(D, d$y))
#' r <- design_column(D, "y") - D %*% b
#' head(as.vector(r))
#'
#' if (requireNamespace("duckdb", quietly = TRUE)) {
#'     con <- DBI::dbConnect(duckdb::duckdb())
#'     d$id <- seq_len(nrow(d))
#'     DBI::dbWriteTable(con, "d", d)
#'     Dr <- design_matrix(y ~ x * f, con, table = "d")
#'     r <- design_column(Dr, "y") - Dr %*% b
#'     r
#'     head(as.data.frame(r, keep = "id", nm = "resid"))
#'     r2 <- col_compute(r^2)
#'     sum(r2)
#'     rt <- col_compute(r, name = "resid", nm = "resid")
#'     head(DBI::dbGetQuery(con, "SELECT id, resid FROM resid ORDER BY id"))
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
as.vector.design_column <- function(x, mode = "any") {
    k <- length(unclass(x)$sql)
    if (k != 1L)
        stop(sprintf(paste("This remote column has %d columns: as.vector() reads a column with one column;",
                           "use as.matrix() for several, or x[, j] for one of them."), k), call. = FALSE)
    dm_column_valid(x)
    as.vector(dmp_collect(x), mode)
}

#' @rdname design_column
#' @export
as.matrix.design_column <- function(x, ...) {
    dm_column_valid(x)
    out <- dmp_collect(x)
    if (is.matrix(out)) {
        dimnames(out) <- NULL
        out
    } else {
        matrix(out, ncol = 1L)
    }
}

#' @rdname design_column
#' @export
col_compute <- function(x, name = NULL, overwrite = FALSE, nm = NULL) {
    ## (in memory a column is the vector itself: nothing to store)
    if (!dm_is_column(x)) return(x)
    checkmate::assert_flag(overwrite, .var.name = "overwrite")
    con <- dm_column_valid(x)
    lab <- unclass(x)$label
    if (is.null(name)) {
        out <- dmc_freeze(x, table = TRUE)
        return(if (is.null(lab)) out else dm_column_label(out, lab))
    }
    checkmate::assert_string(name, min.chars = 1L, .var.name = "name")
    if (!overwrite && DBI::dbExistsTable(con, name))
        stop(sprintf("A table `%s` already exists on the connection: choose another `name`, or use `overwrite = TRUE`.",
                     name), call. = FALSE)
    sql <- unclass(x)$sql
    k <- length(sql)
    nm <- nm %||% lab %||% (if (k == 1L) ".value" else paste0(".value", seq_len(k)))
    checkmate::assert_character(nm, len = k, any.missing = FALSE, min.chars = 1L, unique = TRUE, .var.name = "nm")
    if (any(startsWith(nm, "__dm_"))) stop("`nm`: names starting with `__dm_` are reserved.", call. = FALSE)
    ## The table: the rows of the column's frame (the design's filter), with
    ## every column of that frame and the values, as dmc_freeze(table =
    ## TRUE) makes the temporary one; so the column returned (on a frame
    ## one level deeper, in the lineage of the column's) combines with the
    ## columns of its design.
    fr <- dm_frame_bind(unclass(x)$frame)
    taken <- intersect(nm, names(fr$src$types))
    if (length(taken))
        stop(sprintf("`nm`: the table already has a column `%s` (from the data); name the values with `nm`.",
                     taken[1L]), call. = FALSE)
    tab <- dm_remote_extend(dm_remote_unscaled(fr), stats::setNames(sql, nm), "column")
    tab <- dm_remote_to_table(tab, name = name, overwrite = overwrite)
    tab$depth <- dm_frame_depth(fr) + 1L
    tab <- dm_frame_child(tab, fr)
    dm_column_label(dm_column(tab, sprintf("s.%s", dm_qi(tab, nm)), notnull = unclass(x)$notnull), nm)
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

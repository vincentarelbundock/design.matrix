## Public groupings (Kata design.matrix#0ncx, tier T2 of pg3k): documented
## wrappers over the coded groupings of the fitting core (xtx_c_dense_codes(),
## xtx_c_codes(), xtx_p_groupsum_codes(), xtx_p_crosstab_codes(), xtx_c_lookup(),
## xtx_p_nested(); see primitives.R and primitives_lookup.R). A grouping
## (class design_groups) holds integer codes 1..G of the rows of a design
## (an integer vector in memory, a pending lookup column remotely), G, and
## the levels of the key per code, sorted as the levels of fixed effects.
##
## Remotely, every use of a grouping makes a fresh lookup of its codes (a
## new recipe on the same lookup table, written once by design_groups()):
## the join then binds on top of the other columns of that use, whatever
## frames earlier uses of the grouping settled on, so that a grouping
## combines with any column of its design and with other groupings, used in
## any order. Re-making the recipe runs no statement (the join is part of
## the SQL of the query that reads it).

#' Groupings of the rows of a design
#'
#' A *grouping* gives each row of a design the number of its group, from the
#' values of one or more key variables: integer codes `1..G`, one per level
#' of the key. Grouped sums, crosstabs, lookups of a value per group and
#' nesting tests then work from the codes, in memory and on a remote design
#' (see [design_remote]) with the same code and the same results. These are
#' the operations fixed-effect estimators are built from: for instance, the
#' within transformation of a column is `y - group_lookup(group_sums(y, g) /
#' n_g, g)`.
#'
#' * `design_groups(x, key)` makes the grouping of the rows of design `x` by
#'   `key`.
#' * `group_sums(col, g)`: the sums of a column (or several) within each
#'   group, a `G x k` matrix (one query remotely).
#' * `group_crosstab(g1, g2, w)`: the (weighted) number of rows in each pair
#'   of groups, a sparse `G1 x G2` [Matrix::dgCMatrix-class].
#' * `group_lookup(values, g)`: the column whose value in a row is the value
#'   of its group, `values[codes]`.
#' * `group_nested(g1, g2)`: is `g1` nested in `g2` (does each group of `g1`
#'   lie within one group of `g2`)?
#'
#' @section The grouping object:
#' A list of class `xtx_groups` with elements
#' \describe{
#'   \item{`codes`}{The group of each row, integers `1..G`: in memory an
#'     integer vector with one value per row of the design; for a remote
#'     design, a [design_column()] (lazy: a join of the rows with a table of
#'     the levels, made when its values are first needed). `codes` is a
#'     column like any other: `as.vector(g$codes)` reads it.}
#'   \item{`G`}{The number of groups (an integer).}
#'   \item{`levels`}{A data frame with `G` rows, one column per key variable:
#'     the values of the key in each group, in code order (row `l` is group
#'     `l`).}
#'   \item{`n`}{The number of rows of the design (remotely, `NULL` when not
#'     yet known).}
#'   \item{`design`}{The design.}
#' }
#' Treat it as read-only: the functions check that the codes are within
#' `1..G` before they use them.
#'
#' @section Levels and their order:
#' The groups are the combinations of key values that occur in the rows of
#' the design, sorted as the levels of fixed effects are sorted (as
#' `xtx.fit::lm_xtx()` numbers them): by the first key variable, then the second,
#' and so on; numbers in increasing order, strings in the C locale (as
#' `sort(method = "radix")`: upper case before lower case), `FALSE` before
#' `TRUE`. Logical keys give integer levels `0` and `1`. Factors in memory
#' are sorted in the order of their levels; a remote data source has no
#' factors (strings, and DuckDB `ENUM` columns, are sorted as strings): for
#' the same groups on both backends, use character keys, or factors whose
#' levels are sorted.
#'
#' @section Missing keys:
#' Every row of the design must have a group: a key with a missing value
#' (`NA` or `NaN` in memory, SQL `NULL` or `NaN` remotely) in a row of the
#' design is an error, on both backends. Fits with fixed effects drop such
#' rows; a grouping has no rows of its own, so select the rows first: for a
#' remote design, `D[DBI::SQL("k IS NOT NULL"), ]`; in memory, make the
#' design from the rows with a key. Missing values in the *columns* are
#' kept: as [rowsum()] does, `group_sums()` gives `NA` for a group with a
#' missing value, and `group_crosstab()` an `NA` count for a pair with a
#' missing weight.
#'
#' @section Rows, designs and backends:
#' In memory, columns are plain vectors (or matrices) with one value per row
#' of the design, in the order of its rows. For a remote design, columns are
#' [design_column()]s of that design (or, where a column is expected, the
#' name of a column of the remote data or a `DBI::SQL()` expression): remote
#' rows have no order, so a vector in R cannot be matched with them, and
#' results with one value per row (`group_lookup()`) are columns. A
#' grouping combines only with columns and groupings of the design it was
#' made from (a row selection `D[i, ]` is another design; in memory, where
#' columns are plain vectors, only their length is checked); an in-memory
#' grouping does not combine with remote columns, nor the converse. Results
#' with one row per group (`group_sums()`, `group_crosstab()`) are in code
#' order, in R, without dimension names: row `l` is `g$levels[l, ]`. Both
#' backends give the same values, up to rounding of the sums.
#'
#' A remote grouping is valid while its connection is open (see
#' [design_column()], *Lifetime*): `design_groups()` reads the levels (one
#' `SELECT DISTINCT`) and writes them to a temporary table, and each use
#' joins that table to the rows, in the query that reads them.
#'
#' @param x For `design_groups()`, a design: a `design_matrix` or a remote
#'   design (a row scale, `s * D`, is ignored). For `format()` and
#'   `print()`, a grouping.
#' @param key The key: one or more variables whose combinations of values
#'   are the groups. A character vector is a vector of variable *names*: of
#'   the remote data for a remote design (columns of any type); in memory,
#'   of the variables of the model frame of the design (the response and the
#'   variables of the formula), the only data a design in memory keeps. In
#'   memory, a key can also be given by its values: a vector (numeric,
#'   logical, factor) with one value per row of `x`, or a data frame (such as
#'   `d["firm"]`) or a list of such vectors (strings as values must be in a
#'   data frame or a list, since a character vector is a vector of names).
#'   For a remote design, also a `DBI::SQL()` expression over the columns of
#'   the remote data, or a [design_column()] of `x` (each of its columns a
#'   key variable). A list combines several of these.
#' @param col A column: in memory a numeric (or logical) vector with one
#'   value per row, or a matrix or data frame of such columns; for a remote
#'   design a [design_column()] (with one or more columns). A character
#'   vector of names, or a `DBI::SQL()` expression, is taken as
#'   `design_column(g$design, .)` (in memory, variables of the model frame).
#'   A list of columns is bound into one with several columns. A single
#'   number is a constant column: `group_sums(1, g)` counts the rows of each
#'   group.
#' @param g,g1,g2 Groupings made by `design_groups()`, of the same design,
#'   or bare codes (see [key_codes()]), which are used as given, without
#'   the checks. For `group_sums()`, `g` can also be a list of groupings:
#'   the sums within each of them, from one pass over the rows.
#' @param total For `group_sums()` with a list of groupings: `NULL`, or a
#'   column whose sum over all rows is computed in the same pass.
#' @param w Optional weights, a column in the same form as `col` (one
#'   column). Without `w`, the counts of rows.
#' @param values The value of each group: a numeric vector of length `G`, or
#'   a matrix (or data frame) with `G` rows for several columns, in code
#'   order. Missing values give missing values.
#' @param levels `NULL`, or the levels of a grouping made before (its
#'   `levels`, a data frame with one column per key variable): the rows of
#'   `x` are then coded against those levels, as new data against the levels
#'   of a fit, instead of the levels that occur in `x`. Rows whose key is
#'   missing or not among `levels` get the code `NA` (no error): such a
#'   grouping works with `group_lookup()` (`NA` for those rows) and
#'   the predictions of `xtx.fit::lm_xtx()` with fixed effects, not with
#'   `group_sums()`, `group_crosstab()` or `group_nested()`. Names may change before the first release.
#' @param ... Not used.
#' @return `design_groups()`: a grouping (see *The grouping object*).
#'   `group_sums()` with a list of groupings: a list with `groups` (one
#'   matrix of sums per grouping) and `total`.
#'   `group_sums()`: a `G x k` numeric matrix, `k` the number of columns of
#'   `col`. `group_crosstab()`: a `G1 x G2` [Matrix::dgCMatrix-class] of the
#'   sums of `w` (counts without `w`) in each pair of groups; pairs with no
#'   row are not stored. `group_lookup()`: a column (in memory a double
#'   vector, or a matrix for a matrix of `values`; remotely a
#'   [design_column()]). `group_nested()`: `TRUE` or `FALSE`.
#' @examples
#' set.seed(1)
#' d <- data.frame(y = rnorm(200), x = rnorm(200),
#'                 firm = sample(c("a", "b", "c", "d"), 200, TRUE))
#' d$region <- ifelse(d$firm %in% c("a", "b"), "north", "south")
#' D <- design_matrix(y ~ x, d)
#'
#' # the same code for both backends: the within (one-way fixed effect)
#' # slope of y on x
#' within_slope <- function(D, g) {
#'     n_g <- group_sums(1, g)[, 1]      # rows per group
#'     demean <- function(v) v - group_lookup(group_sums(v, g)[, 1] / n_g, g)
#'     yt <- demean(design_column(D, "y"))
#'     xt <- demean(design_column(D, "x"))
#'     sum(xt * yt) / sum(xt * xt)
#' }
#'
#' g <- design_groups(D, d["firm"])
#' g
#' group_sums(design_column(D, "y"), g)
#' within_slope(D, g)
#' coef(lm(y ~ x + firm, d))[["x"]]
#' group_nested(g, design_groups(D, d["region"]))
#'
#' if (requireNamespace("duckdb", quietly = TRUE)) {
#'     con <- DBI::dbConnect(duckdb::duckdb())
#'     DBI::dbWriteTable(con, "d", d)
#'     R <- design_matrix(y ~ x, con, table = "d")
#'     gr <- design_groups(R, "firm")
#'     gr$levels
#'     group_sums("y", gr)
#'     within_slope(R, gr)
#'     group_nested(gr, design_groups(R, "region"))
#'     group_crosstab(gr, design_groups(R, c("region", "firm")))
#'     DBI::dbDisconnect(con)
#' }
#' @name design_groups
#' @export
design_groups <- function(x, key, levels = NULL) {
    remote <- inherits(x, "xtx_remote")
    if (!remote && !inherits(x, "xtx_matrix"))
        stop("`x` must be a design: a `design_matrix` or a remote design.", call. = FALSE)
    if (remote) x <- xtx_remote_unscaled(x)
    parts <- xtx_groups_parts(key)
    if (!length(parts)) stop("`key` must name at least one variable.", call. = FALSE)
    nms <- xtx_groups_names(parts)
    if (!is.null(levels)) return(xtx_groups_matched(x, parts, nms, levels, remote))
    if (remote) {
        kc <- xtx_groups_key_remote(x, parts)
        cd <- xtx_c_dense_codes(kc)
        lv <- cd$levels
        names(lv) <- nms
        miss <- vapply(lv, anyNA, NA)
        if (any(miss)) stop(xtx_groups_na_msg(nms[miss][1L]), call. = FALSE)
        ## (as in memory, where logical keys are coded as integers)
        for (j in seq_along(lv)) if (is.logical(lv[[j]])) lv[[j]] <- as.integer(lv[[j]])
        codes <- cd$codes$codes
        return(structure(list(codes = codes, G = as.integer(cd$codes$G), levels = lv, n = xtx_nrows(x),
                              design = x),
                         class = "xtx_groups"))
    }
    kd <- xtx_groups_key_memory(x, parts, nms)
    cd <- xtx_c_dense_codes(kd)
    lv <- cd$levels
    names(lv) <- nms
    rownames(lv) <- NULL
    structure(list(codes = cd$codes$codes, G = as.integer(cd$codes$G), levels = lv, n = nrow(x), design = x),
              class = "xtx_groups")
}

#' @rdname design_groups
#' @export
group_sums <- function(col, g, total = NULL) {
    ## bare codes (key_codes()): unchecked; a list of them,
    ## every grouping (and the sum of `total`) in one pass over the rows
    if (xtx_is_codes(g)) return(xtx_p_groupsum(col, g))
    if (is.list(g) && is.null(oldClass(g)))
        return(xtx_c_groupsums(col, lapply(g, xtx_group_codes), total = total))
    if (!is.null(total)) stop("`total` applies only with a list of groupings.", call. = FALSE)
    xtx_groups_check(g, "g")
    col <- xtx_groups_col(col, g, "col")
    out <- xtx_p_groupsum_codes(col, xtx_groups_codes(g))
    dimnames(out) <- NULL
    out
}

#' @rdname design_groups
#' @export
group_crosstab <- function(g1, g2, w = NULL) {
    if (xtx_is_codes(g1) && xtx_is_codes(g2)) return(xtx_p_crosstab(g1, g2, w))
    xtx_groups_check(g1, "g1")
    xtx_groups_check(g2, "g2")
    xtx_groups_same(g1, g2)
    if (!is.null(w)) {
        w <- xtx_groups_col(w, g1, "w")
        if (xtx_ncols(w) != 1L) stop("`w` must be one column.", call. = FALSE)
        if (is.matrix(w)) w <- w[, 1L]
    }
    a <- xtx_groups_codes(g1)
    b <- xtx_groups_codes(g2)
    if (xtx_is_column(a$codes) && !is.null(w)) {
        ## (as xtx_p_crosstab_codes(), with the count of missing weights: a
        ## pair with one is NA, as the sum in memory)
        res <- xtx_group_query(list(a$codes, b$codes), list(w), nulls = TRUE)
        v <- as.double(res$v1)
        v[as.numeric(res$n1) > 0] <- NA_real_
        ct <- .Call(C_xtx_fe_crosstabs, as.integer(res$g2), list(as.integer(res$g1)), v, b$G, a$G,
                    design_threads())[[1L]]
        return(xtx_ct_matrix(ct, a$G, b$G))
    }
    xtx_p_crosstab_codes(a, b, w)
}

#' @rdname design_groups
#' @export
group_lookup <- function(values, g) {
    if (xtx_is_codes(g)) return(xtx_c_lookup(values, g))
    xtx_groups_check(g, "g", na_ok = TRUE)
    if (is.data.frame(values)) values <- as.matrix(values)
    if (!(is.numeric(values) || is.logical(values)) || !(is.null(dim(values)) || length(dim(values)) == 2L))
        stop("`values` must be a numeric vector or matrix.", call. = FALSE)
    if (NROW(values) != g$G)
        stop(sprintf("`values` must have one value (row) per group: %d, not %d.", g$G, NROW(values)), call. = FALSE)
    if (is.matrix(values)) {
        dimnames(values) <- NULL
        storage.mode(values) <- "double"
    } else {
        values <- as.double(values)
    }
    xtx_c_lookup(values, xtx_groups_codes(g))
}

#' @rdname design_groups
#' @export
group_nested <- function(g1, g2) {
    if (xtx_is_codes(g1) && xtx_is_codes(g2)) return(xtx_p_nested(g2, g1))
    xtx_groups_check(g1, "g1")
    xtx_groups_check(g2, "g2")
    xtx_groups_same(g1, g2)
    if (identical(g1, g2)) return(TRUE)
    ## (xtx_p_nested(a, b): is each level of b within one level of a?)
    xtx_p_nested(xtx_groups_codes(g2), xtx_groups_codes(g1))
}

#' @rdname design_groups
#' @export
format.xtx_groups <- function(x, ...) {
    remote <- xtx_is_column(x$codes)
    where <- if (remote) {
        con <- xtx_con(unclass(x$codes)$frame)
        open <- tryCatch(DBI::dbIsValid(con), error = function(e) FALSE)
        sprintf("%s rows, on <%s>%s", if (is.null(x$n)) "?" else xtx_fmt_count(x$n), class(con)[1L],
                if (open) "" else " (closed)")
    } else {
        sprintf("%s rows, in memory", xtx_fmt_count(x$n))
    }
    top <- sprintf("<design_groups: %s group%s of %s; %s>", xtx_fmt_count(x$G), if (identical(x$G, 1L)) "" else "s",
                   paste(names(x$levels), collapse = ", "), where)
    shown <- utils::head(x$levels, 6L)
    body <- if (nrow(shown)) paste0("  ", utils::capture.output(print(shown, row.names = TRUE)))
    if (x$G > 6L) body <- c(body, sprintf("  ... (%s more)", xtx_fmt_count(x$G - 6L)))
    c(top, body)
}

#' @rdname design_groups
#' @export
print.xtx_groups <- function(x, ...) {
    cat(format(x), sep = "\n")
    invisible(x)
}

## ---- internals ----

## The parts of a key: a list of list(v = value, name = its name or NULL,
## is_name = v is the name of a variable). Character vectors at the top
## level are names; the columns of a data frame are values; in a list, a
## single string is a name and other vectors are values.
xtx_groups_parts <- function(key) {
    one <- function(v, nm = NULL, in_df = FALSE) {
        if (xtx_is_column(v)) {
            lab <- unclass(v)$label
            return(lapply(seq_len(xtx_ncols(v)), function(j) list(v = v[, j], name = if (!is.null(nm) && nzchar(nm)) {
                if (xtx_ncols(v) == 1L) nm else paste0(nm, j)
            } else lab[j], is_name = FALSE)))
        }
        if (inherits(v, "SQL"))
            return(lapply(seq_along(v), function(j) list(v = DBI::SQL(as.character(v)[j]),
                                                          name = if (!is.null(nm) && nzchar(nm)) nm else as.character(v)[j],
                                                          is_name = FALSE)))
        if (!in_df && is.character(v) && is.null(dim(v)) && length(v) == 1L && !is.na(v))
            return(list(list(v = v, name = if (!is.null(nm) && nzchar(nm)) nm else v, is_name = TRUE)))
        list(list(v = v, name = if (!is.null(nm) && nzchar(nm)) nm, is_name = FALSE))
    }
    if (xtx_is_column(key) || inherits(key, "SQL")) return(one(key))
    if (is.character(key) && is.null(dim(key))) {
        if (anyNA(key) || any(!nzchar(key))) stop("`key`: variable names must not be missing or empty.", call. = FALSE)
        return(unlist(lapply(key, one), recursive = FALSE))
    }
    if (is.data.frame(key)) {
        nm <- names(key)
        return(unlist(lapply(seq_along(key), function(j) one(key[[j]], nm[j], in_df = TRUE)), recursive = FALSE))
    }
    if (is.list(key) && is.null(oldClass(key))) {
        nm <- names(key) %||% rep("", length(key))
        return(unlist(lapply(seq_along(key), function(j) one(key[[j]], nm[j])), recursive = FALSE))
    }
    one(key)
}

## The names of the key variables (the columns of the levels): given, else
## "key" (one part) or "key1", "key2", ...; unique.
xtx_groups_names <- function(parts) {
    nms <- vapply(parts, function(p) if (is.null(p$name) || is.na(p$name)) "" else as.character(p$name), "")
    blank <- !nzchar(nms)
    nms[blank] <- if (length(parts) == 1L) "key" else paste0("key", seq_along(parts))[blank]
    make.unique(nms)
}

xtx_groups_na_msg <- function(nm) {
    sprintf(paste("`key` has missing values (variable `%s`): every row of the design needs a group.",
                  "Select the rows with a key first (for a remote design, `D[DBI::SQL(\"%s IS NOT NULL\"), ]`;",
                  "in memory, make the design from those rows)."), nm, nm)
}

## The key of a remote design as one column (one SQL expression per part).
xtx_groups_key_remote <- function(D, parts) {
    sqls <- character()
    for (p in parts) {
        v <- p$v
        kc <- if (xtx_is_column(v)) {
            if (!identical(unclass(v)$frame$id, D$id))
                stop("`key`: a column of another remote design.", call. = FALSE)
            xtx_remote_cols(D, v, "key")
        } else if (inherits(v, "SQL")) {
            xtx_remote_cols(D, v, "key")
        } else if (p$is_name) {
            xtx_remote_cols(D, v, "key", numeric = FALSE)
        } else {
            stop(paste("`key`: for a remote design, give names of columns of the remote data,",
                       "DBI::SQL() expressions or design_column()s of the design (remote rows have no",
                       "order, so values in R cannot be matched with them)."), call. = FALSE)
        }
        ## (a key column may move the design to a later frame, which keeps
        ## the columns read by the earlier parts)
        D <- kc$D
        sqls <- c(sqls, kc$sql)
    }
    xtx_column(D, sqls)
}

## The key of a design in memory as a data frame (one column per part).
xtx_groups_key_memory <- function(x, parts, nms) {
    n <- nrow(x)
    vars <- (if (xtx_is_comb(x)) x$design else x)$vars
    cols <- lapply(seq_along(parts), function(j) {
        p <- parts[[j]]
        v <- p$v
        if (xtx_is_column(v)) stop("`key`: a remote column, and the design is in memory.", call. = FALSE)
        if (inherits(v, "SQL"))
            stop("`key`: a DBI::SQL() expression is a key only of a remote design; in memory, pass the values.",
                 call. = FALSE)
        if (p$is_name) {
            if (is.null(vars) || !v %in% names(vars))
                stop(sprintf(paste("`%s` is not a variable of the design: in memory, a design keeps only the",
                                   "variables of its model frame; pass the values instead (for instance",
                                   "`d[\"%s\"]`, one value per row of the design)."), v, v), call. = FALSE)
            cls <- unname(attr(attr(vars, "terms"), "dataClasses")[v])
            v <- vars[[v]]
            ## (a character variable of the data may be kept as a factor whose
            ## levels are sorted in the locale: its key sorts as strings, as
            ## the same values given as a vector, or a remote key)
            if (is.factor(v) && identical(cls, "character")) v <- as.character(v)
        }
        if (!is.atomic(v) || !is.null(dim(v)) || !(is.numeric(v) || is.logical(v) || is.character(v) || is.factor(v)))
            stop(sprintf("`key`: variable `%s` must be a vector (numeric, logical, character or factor).", nms[j]),
                 call. = FALSE)
        if (length(v) != n)
            stop(sprintf("`key`: variable `%s` has %d values, not one per row of the design (%d).", nms[j],
                         length(v), n), call. = FALSE)
        if (anyNA(v)) stop(xtx_groups_na_msg(nms[j]), call. = FALSE)
        ## (other attributes, such as names, would change the coding path)
        if (!is.factor(v)) attributes(v) <- NULL
        v
    })
    names(cols) <- nms
    structure(cols, class = "data.frame", row.names = c(NA_integer_, -n))
}

## Stop unless g is a grouping whose codes can be used (in memory: integer
## codes within 1..G, one per row: the kernels index by them).
xtx_groups_check <- function(g, what, na_ok = FALSE) {
    if (!inherits(g, "xtx_groups"))
        stop(sprintf("`%s` must be a grouping made by design_groups().", what), call. = FALSE)
    G <- g$G
    if (!is.integer(G) || length(G) != 1L || is.na(G) || G < 0L || !is.data.frame(g$levels) || nrow(g$levels) != G)
        stop(sprintf("`%s` is not a valid grouping (G and levels).", what), call. = FALSE)
    cd <- g$codes
    if (xtx_is_column(cd)) {
        if (!inherits(g$design, "xtx_remote") || xtx_ncols(cd) != 1L)
            stop(sprintf("`%s` is not a valid grouping (remote codes).", what), call. = FALSE)
        xtx_column_valid(cd)
        return(invisible(g))
    }
    if (!is.integer(cd) || !is.null(dim(cd)) || length(cd) != g$n)
        stop(sprintf("`%s$codes` must be integer codes, one per row of the design.", what), call. = FALSE)
    if (length(cd)) {
        if (isTRUE(g$matched) && anyNA(cd)) {
            if (!na_ok)
                stop(sprintf("`%s` has rows with no level (codes NA, from `levels`): use it with group_lookup().",
                             what), call. = FALSE)
            cd <- cd[!is.na(cd)]
        }
        r <- if (length(cd)) range(cd) else c(1L, 1L)
        if (anyNA(r) || r[1L] < 1L || r[2L] > G)
            stop(sprintf("`%s$codes` must be within 1..G (%d), without missing values.", what, G), call. = FALSE)
    }
    invisible(g)
}

## Stop unless two groupings are of the same design.
xtx_groups_same <- function(g1, g2) {
    r1 <- xtx_is_column(g1$codes)
    r2 <- xtx_is_column(g2$codes)
    if (r1 != r2) stop("One grouping is in memory and the other remote.", call. = FALSE)
    same <- if (r1) {
        identical(g1$design$id, g2$design$id)
    } else {
        t1 <- xtx_groups_token(g1$design)
        t2 <- xtx_groups_token(g2$design)
        isTRUE(g1$n == g2$n) && (is.null(t1) || is.null(t2) || identical(t1, t2))
    }
    if (!same) stop("The groupings are of different designs.", call. = FALSE)
    invisible(TRUE)
}

## (in memory, a design made from a model frame keeps an environment of
## its own, which its scaled versions share; NULL for a dense design)
xtx_groups_token <- function(x) {
    e <- (if (xtx_is_comb(x)) x$design else x)$cache
    if (is.environment(e)) e
}

## The codes of a grouping as xtx_c_codes(). Remotely the lookup of the
## codes that design_groups() made: a node, which each read joins to its
## source with its other operands.
xtx_groups_codes <- function(g) xtx_c_codes(g$codes, g$G)

## A column for grouping g (see the documentation of `col`): in memory a
## double vector or matrix with one row per row of the design; remotely a
## design_column of the design of g.
xtx_groups_col <- function(col, g, what) {
    remote <- xtx_is_column(g$codes)
    ## (a number: a constant column)
    if (is.numeric(col) && length(col) == 1L && is.null(dim(col))) {
        if (!is.finite(col)) stop(sprintf("`%s`: a constant column must be finite.", what), call. = FALSE)
        return(if (remote) xtx_c_const(g$design, col) else rep(as.double(col), g$n))
    }
    if (is.character(col) || inherits(col, "SQL")) {
        if (!length(col)) stop(sprintf("`%s` must name at least one column.", what), call. = FALSE)
        cols <- if (inherits(col, "SQL")) lapply(as.character(col), DBI::SQL) else as.list(col)
        col <- lapply(cols, function(cn) design_column(g$design, cn))
        if (length(col) == 1L) col <- col[[1L]]
    }
    if (is.list(col) && !is.data.frame(col) && !xtx_is_column(col)) {
        if (!length(col)) stop(sprintf("`%s` has no columns.", what), call. = FALSE)
        col <- if (any(vapply(col, xtx_is_column, NA))) do.call(xtx_c_cbind, unname(col)) else do.call(cbind, unname(col))
    }
    if (remote) {
        if (!xtx_is_column(col))
            stop(sprintf(paste("`%s`: for a remote grouping, a column of its design (a design_column(), a",
                               "column name or DBI::SQL()); remote rows have no order, so values in R cannot",
                               "be matched with them."), what), call. = FALSE)
        if (!identical(unclass(col)$frame$id, g$design$id))
            stop(sprintf("`%s` is a column of another remote design than the grouping.", what), call. = FALSE)
        xtx_column_valid(col)
        return(col)
    }
    if (xtx_is_column(col)) stop(sprintf("`%s` is a remote column, and the grouping is in memory.", what), call. = FALSE)
    if (is.data.frame(col)) col <- as.matrix(col)
    if (!(is.numeric(col) || is.logical(col)) || !(is.null(dim(col)) || length(dim(col)) == 2L))
        stop(sprintf("`%s` must be a numeric vector or matrix.", what), call. = FALSE)
    if (NROW(col) != g$n)
        stop(sprintf("`%s` must have one value (row) per row of the design: %d, not %d.", what, g$n, NROW(col)),
             call. = FALSE)
    if (is.matrix(col)) {
        dimnames(col) <- NULL
        if (!is.double(col)) storage.mode(col) <- "double"
    } else if (!is.double(col) || !is.null(attributes(col))) {
        col <- as.double(col)
    }
    col
}

## design_groups(x, key, levels): the rows of x coded against `levels` (a
## data frame of the levels of an earlier grouping, one column per key
## part), NA for keys that are missing or not among them; as
## xtx_c_lookup() matches keys to levels (in memory: the values as text,
## doubles to 17 digits; remotely a join on the key columns, text levels
## compared as text).
xtx_groups_matched <- function(x, parts, nms, levels, remote) {
    if (!is.data.frame(levels) || ncol(levels) != length(parts))
        stop("`levels` must be a data frame with one column per key variable.", call. = FALSE)
    lv <- levels
    rownames(lv) <- NULL
    G <- nrow(lv)
    if (remote) {
        kc <- xtx_groups_key_remote(x, parts)
        codes <- xtx_c_lookup(seq_len(G), kc, lv)
        return(structure(list(codes = codes, G = as.integer(G), levels = lv, n = xtx_nrows(x), design = x,
                              matched = TRUE),
                         class = "xtx_groups"))
    }
    ## (the values as given, or the variables of the design's model frame)
    vars <- (if (xtx_is_comb(x)) x$design else x)$vars
    kd <- lapply(seq_along(parts), function(j) {
        p <- parts[[j]]
        if (xtx_is_column(p$v) || inherits(p$v, "SQL"))
            stop("`key`: in memory, give the values of the key (or names of variables of the design).", call. = FALSE)
        v <- if (p$is_name) {
            if (is.null(vars) || !p$v %in% names(vars)) stop(sprintf("No variable `%s` in the data.", p$v), call. = FALSE)
            vars[[p$v]]
        } else p$v
        if (length(v) != nrow(x))
            stop(sprintf("`key`: variable `%s` has %d values, not one per row of the design (%d).", nms[j],
                         length(v), nrow(x)), call. = FALSE)
        if (is.factor(v)) as.character(v) else v
    })
    kd <- as.data.frame(stats::setNames(kd, nms), stringsAsFactors = FALSE)
    codes <- match(xtx_lookup_keys(kd), xtx_lookup_keys(lv))
    structure(list(codes = as.integer(codes), G = as.integer(G), levels = lv, n = nrow(x), design = x, matched = TRUE),
              class = "xtx_groups")
}

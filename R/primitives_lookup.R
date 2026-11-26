## Lookup columns, keys and codes: primitives of the fitting core (see
## primitives.R for the column conventions), used by fixed effects (fe.R).
## Each has an in-memory and a remote method; code that calls them does not
## branch on the backend.
##
##   dmc_lookup(values, key, levels = NULL, along = list())
##       the column whose value in a row is values[l], l the level of the
##       row's key: a lookup table joined on the key. Without `levels`, the
##       key is one column of integer codes 1..length(values) (or a
##       dmc_codes() grouping); with `levels` (a data frame, one row per
##       value, one column per key part), the key has one column per part
##       and a row's level is the row of `levels` with the same values. Keys
##       with no level (codes out of range, new values, NA) give NA (SQL
##       NULL).
##       In memory: values[codes] (C), or values[match(keys)] (key: a vector,
##       or a list or data frame of vectors, one per part).
##       Remotely: the lookup table is written to a temporary table and
##       joined (LEFT JOIN) to the frame of the key (or the later frame of
##       `along`, as for dmp_mult()), which gives a new frame keeping every
##       column of the old one.
##   dmc_key(x, vars, data = NULL)
##       the values of variables (columns of the data, of any type) as a key
##       over the rows of design x: in memory a data frame (from `data`, or
##       from the model frame the design keeps); remotely one SQL column per
##       variable.
##   dmc_not_na(key)
##       a logical column: no part of the key is missing.
##   dmc_dense_codes(key, along = list())
##       list(codes = dmc_codes(), levels = data frame): the levels of the
##       key that occur, sorted, and each row's level as an integer code. In
##       memory from the values (integer coding in C); remotely a SELECT
##       DISTINCT of the levels and a lookup of their codes.
##   dmc_subset(col, keep, x2)
##       a column of x restricted to the rows of x2 = dmp_subset(x, keep): in
##       memory col[keep]; remotely the same SQL over x2 (dmc_on()).
##   dmp_nested(a, b)
##       for two dmc_codes() groupings: is each level of b in one level of
##       a? (In memory one pass in C; remotely from the crosstab.)
##   dmp_center(x)
##       the design with numeric variables shifted by their mean (in a
##       sample of rows), when that design is X T with T unit upper
##       triangular: list(x = the shifted design, T), or NULL. In memory the
##       design is rebuilt from its model frame; remotely the variables'
##       SQL (or the chunks' shift) changes, as in dm_center(). The shifted
##       design is a new design: re-make columns with dmc_on().

dmc_lookup <- function(values, key, levels = NULL, along = list()) {
    values <- as.double(values)
    if (dm_is_codes(key)) key <- key$codes
    if (!dm_is_rcols(key)) {
        if (is.null(levels)) {
            if (!is.integer(key)) key <- as.integer(key)
            return(.Call(C_dm_fe_lookup, values, key))
        }
        kd <- if (is.data.frame(key) || is.list(key)) as.data.frame(key, stringsAsFactors = FALSE)
              else data.frame(key, stringsAsFactors = FALSE)
        if (ncol(kd) != ncol(levels)) stop("`key` and `levels` have different numbers of columns.", call. = FALSE)
        for (j in seq_along(kd)) if (is.factor(kd[[j]])) kd[[j]] <- as.character(kd[[j]])
        return(values[match(dm_lookup_keys(kd), dm_lookup_keys(levels))])
    }
    frame <- do.call(dm_rcols_frame, c(list(key), as.list(along)))
    con <- dm_con(frame)
    q <- function(x) dm_qi(frame, x)
    ksql <- unclass(key)$sql
    if (is.null(levels)) {
        tab <- data.frame(code = seq_along(values), v = values)
        on <- sprintf("CAST(%s AS BIGINT) = l.code", ksql[1L])
    } else {
        if (length(ksql) != ncol(levels)) stop("`key` and `levels` have different numbers of columns.", call. = FALSE)
        lcols <- paste0("k", seq_along(ksql))
        tab <- stats::setNames(as.data.frame(levels, stringsAsFactors = FALSE), lcols)
        for (j in seq_along(tab)) if (is.factor(tab[[j]])) tab[[j]] <- as.character(tab[[j]])
        tab$v <- values
        on <- vapply(seq_along(ksql), function(j) {
            if (is.character(tab[[j]])) sprintf("CAST(%s AS VARCHAR) = l.%s", ksql[j], lcols[j])
            else sprintf("%s = l.%s", ksql[j], lcols[j])
        }, "")
        on <- paste(on, collapse = " AND ")
    }
    tn <- dm_tmp_name("dm_lookup")
    DBI::dbWriteTable(con, tn, tab, temporary = TRUE)
    col <- paste0("__dm_lk", substr(tn, 11L, 18L))
    fr <- dm_remote_unscaled(frame)
    fr$src$select <- sprintf("SELECT s.*, l.v AS %s FROM (%s) AS s LEFT JOIN %s AS l ON %s",
                             q(col), frame$src$select, q(tn), on)
    fr$src$base_table <- NULL
    fr$src$proto <- DBI::dbGetQuery(con, sprintf("SELECT * FROM (%s) AS s LIMIT 0", fr$src$select))
    fr$src$types <- vapply(fr$src$proto, dm_col_type, "")
    fr$depth <- dm_frame_depth(frame) + 1L
    dm_rcols(fr, sprintf("s.%s", q(col)))
}

## One string per row of a data frame of key values (doubles to 17 digits).
dm_lookup_keys <- function(d) {
    v <- lapply(d, function(x) if (is.double(x)) format(x, digits = 17, trim = TRUE) else as.character(x))
    if (length(v) == 1L) return(v[[1L]])
    do.call(paste, c(v, sep = "\r"))
}

dmc_key <- function(x, vars, data = NULL) {
    if (inherits(x, "dm_remote")) {
        miss <- setdiff(vars, names(x$src$types))
        if (length(miss)) stop(sprintf("No column `%s` in the remote data.", miss[1L]), call. = FALSE)
        return(dm_rcols(x, sprintf("s.%s", dm_qi(x, vars))))
    }
    d <- if (is.null(data)) x$vars else data
    miss <- setdiff(vars, names(d))
    if (length(miss)) stop(sprintf("No variable `%s` in the data.", miss[1L]), call. = FALSE)
    out <- as.data.frame(lapply(vars, function(v) d[[v]]), stringsAsFactors = FALSE)
    names(out) <- vars
    out
}

dmc_not_na <- function(key) {
    if (!dm_is_rcols(key)) return(if (anyNA(key)) stats::complete.cases(key) else rep(TRUE, NROW(key)))
    sql <- unclass(key)$sql
    dm_rcols(unclass(key)$frame, sprintf("(%s)", paste(sprintf("(%s) IS NOT NULL", sql), collapse = " AND ")))
}

dmc_dense_codes <- function(key, along = list()) {
    if (!dm_is_rcols(key)) {
        cd <- dm_fe_code_term(key, names(key))
        return(list(codes = dmc_codes(cd$codes, cd$G), levels = cd$levels))
    }
    frame <- unclass(key)$frame
    sql <- unclass(key)$sql
    lv <- dm_query(frame, sprintf("SELECT DISTINCT %s FROM (%s) AS s WHERE %s",
                                  paste(sprintf("%s AS k%d", sql, seq_along(sql)), collapse = ", "),
                                  frame$src$select, frame$filter))
    for (j in seq_along(lv)) if (is.factor(lv[[j]])) lv[[j]] <- as.character(lv[[j]])
    lv <- lv[do.call(order, unname(as.list(lv))), , drop = FALSE]
    rownames(lv) <- NULL
    codes <- dmc_lookup(seq_len(nrow(lv)), key, lv, along)
    list(codes = dmc_codes(codes, nrow(lv)), levels = lv)
}

dmc_subset <- function(col, keep, x2) {
    if (is.null(col)) return(NULL)
    if (dm_is_codes(col)) stop("Codes change with the rows: compute them again.", call. = FALSE)
    if (dm_is_rcols(col)) return(dmc_on(col, x2))
    k <- !is.na(keep) & keep
    if (is.data.frame(col) || is.matrix(col)) col[k, , drop = FALSE] else col[k]
}

dmp_nested <- function(a, b) {
    if (identical(a, b)) return(TRUE)
    if (!dm_is_rcols(a$codes) && !dm_is_rcols(b$codes))
        return(.Call(C_dm_fe_nested, a$codes, b$codes, b$G))
    ct <- dmp_crosstab(a, b) # Ga x Gb
    all(diff(ct@p) == 1L)
}

dmp_center <- function(x) {
    if (inherits(x, "dm_remote")) {
        if (!isTRUE(getOption("design.matrix.remote_center", TRUE))) return(NULL)
        ctr <- dm_center(x)
        if (is.null(ctr)) return(NULL)
        x2 <- x
        x2$cache <- new.env(parent = emptyenv())
        x2$id <- dm_tmp_name("dm")
        if (x$engine == "sql") {
            for (nm in names(ctr$shift))
                x2$vars[[nm]]$sql <- sprintf("(%s - %s)", x$vars[[nm]]$sql, dm_sql_double(ctr$shift[[nm]]))
        } else {
            x2$shift <- ctr$shift
        }
        return(list(x = x2, T = ctr$T))
    }
    mt <- x$training_terms
    mf <- x$vars
    if (is.null(mt) || is.null(mf)) return(NULL)
    ctr <- dm_fe_center(mt, mf, x$contrasts)
    if (is.null(ctr)) return(NULL)
    for (nm in names(ctr$shift)) mf[[nm]] <- mf[[nm]] - ctr$shift[[nm]]
    list(x = design_matrix(mt, mf, contrasts.arg = x$contrasts), T = ctr$T)
}

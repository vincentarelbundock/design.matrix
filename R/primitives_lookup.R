## Lookup columns: a primitive of the fitting core (see primitives.R for
## the column conventions), used by fixed effects (fe.R).
##
##   dmc_lookup(values, key, levels = NULL, along = list())
##       the column whose value in a row is values[l], l the level of the
##       row's key: a lookup table joined on the key. Without `levels`, the
##       key is one column of integer codes 1..length(values); with
##       `levels` (a data frame, one row per value, one column per key
##       part), the key has one column per part and a row's level is the
##       row of `levels` with the same values. Keys with no level (codes
##       out of range, new values, NA) give NA (SQL NULL).
##       In memory: values[codes], or values[match(keys)] (key: a vector,
##       or a list or data frame of vectors, one per part).
##       Remotely: the lookup table is written to a temporary table and
##       joined (LEFT JOIN) to the frame of the key (or the later frame of
##       `along`, as for dmp_mult()), which gives a new frame keeping every
##       column of the old one.
##   dmc_cbind(...)
##       columns side by side (k columns: the sum of theirs); remotely they
##       must come from the same design (the later frame is taken).
##   dmc_const(like, value)
##       a column with the same rows as the column `like`, every value
##       `value`.

dmc_lookup <- function(values, key, levels = NULL, along = list()) {
    values <- as.double(values)
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

dmc_cbind <- function(...) {
    cols <- list(...)
    if (!any(vapply(cols, dm_is_rcols, NA))) return(do.call(cbind, cols))
    frame <- do.call(dm_rcols_frame, cols)
    dm_rcols(frame, unlist(lapply(cols, dm_rcols_operand)))
}

dmc_const <- function(like, value) {
    if (!dm_is_rcols(like)) return(rep(as.double(value), NROW(like)))
    dm_rcols(unclass(like)$frame, dm_sql_double(value))
}

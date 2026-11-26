## Lookup columns, keys and codes: primitives of the fitting core (see
## primitives.R for the column conventions), used by fixed effects (fe.R).
## Each has an in-memory and a remote method; code that calls them does not
## branch on the backend.
##
##   dmc_lookup(values, key, levels = NULL)
##       the column whose value in a row is values[l], l the level of the
##       row's key: a lookup table joined on the key. Without `levels`, the
##       key is one column of integer codes 1..length(values) (or a
##       dmc_codes() grouping); with `levels` (a data frame, one row per
##       value, one column per key part), the key has one column per part
##       and a row's level is the row of `levels` with the same values. Keys
##       with no level (codes out of range, new values, NA) give NA (SQL
##       NULL). `values` can be a G x k matrix: k columns, from one lookup
##       table (one join remotely).
##       In memory: values[codes] (C), or values[match(keys)] (key: a vector,
##       or a list or data frame of vectors, one per part), as doubles.
##       Remotely: the lookup table is written to a temporary table at
##       once, and the column is pending (see primitives.R): at its first
##       use the table is joined (LEFT JOIN) to the latest frame among the
##       operands of that use (at least the frame of the key), which gives
##       a new frame keeping every column of the old one; integer `values`
##       (codes) give an integer column.
##   dmc_lookup_flags(flags, key)
##       for a G x m logical matrix of flags of the levels and key a
##       dmc_codes() grouping (codes 1..G): the m logical columns of the
##       rows' flags. In memory flags[codes, j]; remotely one dmc_lookup()
##       of the flags as 0/1.
##   dmc_key(x, vars, data = NULL)
##       the values of variables (columns of the data, of any type) as a key
##       over the rows of design x: in memory a data frame (from `data`, or
##       from the model frame the design keeps); remotely one SQL column per
##       variable.
##   dmc_complete(key)
##       a logical column: no part of the key is missing.
##   dmc_dense_codes(key)
##       list(codes = dmc_codes(), levels = data frame): the levels of the
##       key that occur, sorted, and each row's level as an integer code. In
##       memory from the values (integer coding in C); remotely a SELECT
##       DISTINCT of the levels and a lookup of their codes (pending: codes
##       of several keys that are used one by one but must combine later
##       are bound together with dmc_bind()).
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
##       design is rebuilt from its model frame, its columns as stored (a
##       column I(x^2) or log(x) is shifted on its own, not rebuilt from a
##       shifted x); remotely the variables' SQL (or the chunks' shift)
##       changes, as in dm_center(), each variable on its own. The shifted
##       design is a new design: re-make columns with dmc_on().

dmc_lookup <- function(values, key, levels = NULL) {
    mat <- is.matrix(values)
    ## (remotely, integer values stay integer: codes, which group and join
    ## faster and arrive in R as integers)
    int <- is.integer(values) && !mat
    values <- if (mat) {
        storage.mode(values) <- "double"
        values
    } else {
        as.double(values)
    }
    if (dm_is_codes(key)) {
        levels <- key$levels %||% levels
        key <- key$codes
    }
    if (!dm_is_column(key)) {
        if (is.null(levels)) {
            if (!is.integer(key)) key <- as.integer(key)
            ## (a matrix of values gives an n x k matrix, in one call)
            return(.Call(C_dm_fe_lookup, values, key))
        }
        kd <- if (is.data.frame(key) || is.list(key)) as.data.frame(key, stringsAsFactors = FALSE)
              else data.frame(key, stringsAsFactors = FALSE)
        if (ncol(kd) != ncol(levels)) stop("`key` and `levels` have different numbers of columns.", call. = FALSE)
        for (j in seq_along(kd)) if (is.factor(kd[[j]])) kd[[j]] <- as.character(kd[[j]])
        i <- match(dm_lookup_keys(kd), dm_lookup_keys(levels))
        return(if (mat) values[i, , drop = FALSE] else values[i])
    }
    vm <- if (mat) values else matrix(values, ncol = 1L)
    vcols <- paste0("v", seq_len(ncol(vm)))
    ## (the key's frame, possibly pending: the lookup binds on top of the
    ## key's recipes)
    frame <- unclass(key)$frame
    con <- dm_con(frame)
    q <- function(x) dm_qi(frame, x)
    ksql <- unclass(key)$sql
    if (is.null(levels)) {
        tab <- data.frame(code = seq_len(nrow(vm)))
        on <- sprintf("CAST(%s AS BIGINT) = l.code", ksql[1L])
    } else {
        if (length(ksql) != ncol(levels)) stop("`key` and `levels` have different numbers of columns.", call. = FALSE)
        lcols <- paste0("k", seq_along(ksql))
        tab <- stats::setNames(as.data.frame(levels, stringsAsFactors = FALSE), lcols)
        for (j in seq_along(tab)) if (is.factor(tab[[j]])) tab[[j]] <- as.character(tab[[j]])
        on <- vapply(seq_along(ksql), function(j) {
            if (is.character(tab[[j]])) sprintf("CAST(%s AS VARCHAR) = l.%s", ksql[j], lcols[j])
            else sprintf("%s = l.%s", ksql[j], lcols[j])
        }, "")
        on <- paste(on, collapse = " AND ")
    }
    for (j in seq_along(vcols)) tab[[vcols[j]]] <- if (int) as.integer(vm[, j]) else vm[, j]
    tn <- dm_tmp_name("dm_lookup")
    dm_write_temp(con, tn, tab)
    col <- paste0("__dm_lk", substr(tn, 11L, 18L), if (mat) paste0("_", seq_along(vcols)))
    ## (the join is made at the first use: dm_bind_lookup())
    r <- dm_recipe("lookup", deps = frame$pend %||% list(),
                   items = paste(sprintf("l.%s AS %s", vcols, q(col)), collapse = ", "),
                   table = q(tn), on = on, col = col, type = if (int) integer(0) else double(0))
    fr <- dm_remote_unscaled(frame)
    fr$pend <- c(frame$pend, list(r))
    dm_column(fr, sprintf("s.%s", q(col)))
}

dmc_lookup_flags <- function(flags, key) {
    if (dm_plain_codes(key)) return(lapply(seq_len(ncol(flags)), function(j) flags[, j][key$codes]))
    storage.mode(flags) <- "double"
    lk <- dmc_lookup(flags, key)
    lapply(seq_len(ncol(flags)), function(j) dmc_col(lk, j) == 1)
}

## The columns of the frame of a remote design (or column) plus the new
## columns `cols` (all of the prototype `type`): from the frame's
## prototype when it has one, else from the database (query `sql`).
dm_proto_add <- function(frame, cols, type, sql) {
    proto <- frame$src$proto
    if (is.null(proto) || any(cols %in% names(proto))) return(DBI::dbGetQuery(dm_con(frame), sql))
    for (cn in cols) proto[[cn]] <- type
    proto
}

## A data frame written to the temporary table `name`: registered with
## DuckDB and copied by one statement (DBI::dbWriteTable() first probes
## for the table, which costs much more than the copy of a lookup table).
dm_write_temp <- function(con, name, df) {
    reg <- paste0(name, "_r")
    duckdb::duckdb_register(con, reg, df)
    on.exit(duckdb::duckdb_unregister(con, reg))
    DBI::dbExecute(con, sprintf("CREATE TEMPORARY TABLE %s AS SELECT * FROM %s",
                                DBI::dbQuoteIdentifier(con, name), DBI::dbQuoteIdentifier(con, reg)))
    invisible(name)
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
        return(dm_column(x, sprintf("s.%s", dm_qi(x, vars))))
    }
    d <- if (is.null(data)) x$vars else data
    miss <- setdiff(vars, names(d))
    if (length(miss)) stop(sprintf("No variable `%s` in the data.", miss[1L]), call. = FALSE)
    out <- as.data.frame(lapply(vars, function(v) d[[v]]), stringsAsFactors = FALSE)
    names(out) <- vars
    out
}

## TRUE where no key column is missing: NA or NaN (as complete.cases() in
## memory; was dmc_not_na, renamed in tapx: dmc_present() is the same rule
## for one numeric column). Remotely, a key column of the data (s."v") that is not double
## cannot be NaN: only IS NOT NULL (isnan() would not take text); a double
## column, or any other expression, is NaN-checked as well.
dmc_complete <- function(key) {
    if (!dm_is_column(key)) return(if (anyNA(key)) stats::complete.cases(key) else rep(TRUE, NROW(key)))
    frame <- unclass(key)$frame
    sql <- unclass(key)$sql
    types <- frame$src$types
    cols <- names(types)[match(sql, sprintf("(s.%s)", dm_qi(frame, names(types))))]
    plain <- !is.na(cols) & !types[cols] %in% "double"
    one <- ifelse(plain, sprintf("(%s) IS NOT NULL", sql),
                  sprintf("COALESCE(NOT isnan(CAST(%s AS DOUBLE)), FALSE)", sql))
    dm_column(frame, sprintf("(%s)", paste(one, collapse = " AND ")))
}

dmc_dense_codes <- function(key, raw = FALSE) {
    if (!dm_is_column(key)) {
        cd <- dm_fe_code_term(key, names(key))
        ## (in memory, also the rows per level: n)
        return(list(codes = dmc_codes(cd$codes, cd$G), levels = cd$levels, n = cd$n))
    }
    ## (the levels are read now: the key on its settled frame)
    key <- dm_column_bind(key)
    frame <- unclass(key)$frame
    sql <- unclass(key)$sql
    fields <- sprintf("%s AS k%d", sql, seq_along(sql))
    composite <- raw && length(sql) > 1L
    if (composite) {
        ## A lossless scalar key for GROUPING SETS and joins on a tuple.
        ## Keep the original columns for R's level order and fixef labels.
        parts <- paste(sprintf("%s, 'ASC NULLS LAST'", sql), collapse = ", ")
        tuple <- sprintf("hex(create_sort_key(%s))", parts)
        fields <- c(fields, sprintf("%s AS tuple", tuple))
    }
    select <- if (raw) "SELECT" else "SELECT DISTINCT"
    if (raw) fields <- c(fields, "COUNT(*) AS n")
    lv <- dm_query(frame, paste(select, paste(fields, collapse = ", "), "FROM (",
                                frame$src$select, ") AS s WHERE", frame$filter,
                                if (raw) "GROUP BY ALL" else ""))
    counts <- lv$n
    tuple_values <- lv$tuple
    lv$n <- lv$tuple <- NULL
    for (j in seq_along(lv)) if (is.factor(lv[[j]])) lv[[j]] <- as.character(lv[[j]])
    ord <- do.call(order, c(unname(as.list(lv)), method = "radix"))
    lv <- lv[ord, , drop = FALSE]
    rownames(lv) <- NULL
    if (!raw) {
        codes <- dmc_lookup(seq_len(nrow(lv)), key, lv)
        return(list(codes = dmc_codes(codes, nrow(lv)), levels = lv))
    }
    g <- dmc_codes(if (composite) dm_column(frame, tuple) else key, nrow(lv))
    g$levels <- if (composite) data.frame(tuple = tuple_values[ord]) else lv
    list(codes = g, levels = lv, n = counts[ord])
}

## ---- codes in memory (dmc_dense_codes()) ----

## Integer codes 1..G of a vector without missing values, its levels (the
## values that occur, sorted, or in level order for a factor) and the rows
## per level (n), in one pass in C (dm_fe_code() in src/fe.c: direct
## indexing for integer values with a compact range, hashing otherwise;
## the levels are then sorted as sort(unique(v), method = "radix") sorts
## them). Other types, and strings that occur in several encodings, as
## sort(unique()) and match().
dm_fe_code <- function(v) {
    if (is.factor(v)) {
        r <- .Call(C_dm_fe_code, as.integer(v), nlevels(v) + 1)
        if (!is.null(r)) return(list(codes = r[[1L]], G = r[[2L]], levels = levels(v)[r[[3L]]], n = r[[4L]]))
        cd <- dm_fe_compact(as.integer(v), nlevels(v))
        return(list(codes = cd$codes, G = cd$G, levels = levels(v)[cd$used]))
    }
    if (is.logical(v)) v <- as.integer(v)
    if (is.numeric(v) || (is.character(v) && is.null(attributes(v)))) {
        r <- .Call(C_dm_fe_code, v, max(4 * length(v), 1e6))
        if (!is.null(r)) {
            if (length(r) == 4L) return(list(codes = r[[1L]], G = r[[2L]], levels = r[[3L]], n = r[[4L]]))
            ## hashed: levels in order of first appearance, sorted here
            ## (as plain values: not for a class with its own order)
            u <- r[[3L]]
            if (is.null(oldClass(v)) && (!is.character(u) || !anyDuplicated(u))) {
                o <- order(u, method = "radix")
                rank <- integer(length(o))
                rank[o] <- seq_along(o)
                return(list(codes = .Call(C_dm_fe_recode, r[[1L]], rank), G = r[[2L]], levels = u[o],
                            n = r[[4L]][o]))
            }
        }
    }
    u <- sort(unique(v), method = "radix")
    list(codes = match(v, u), G = length(u), levels = u)
}

## Renumber codes so that the levels that occur are 1..G.
dm_fe_compact <- function(cd, G) {
    used <- tabulate(cd, nbins = G) > 0L
    if (all(used)) return(list(codes = cd, G = G, used = used))
    map <- cumsum(used)
    list(codes = map[cd], G = sum(used), used = used)
}

## Codes of an FE (one or several variables, columns of `d`), with its
## levels as a data frame with one column per variable.
dm_fe_code_term <- function(d, vars) {
    if (length(vars) == 1L) {
        cd <- dm_fe_code(d[[vars]])
        lv <- list(cd$levels)
        names(lv) <- vars
        return(list(codes = cd$codes, G = cd$G, levels = as.data.frame(lv, stringsAsFactors = FALSE), n = cd$n))
    }
    parts <- lapply(vars, function(v) dm_fe_code(d[[v]]))
    cd <- parts[[1L]]$codes
    for (k in seq_along(parts)[-1L]) cd <- (cd - 1) * parts[[k]]$G + parts[[k]]$codes
    comb <- dm_fe_code(cd)
    ## the components of each combined level
    code <- comb$levels - 1
    lv <- vector("list", length(vars))
    for (k in rev(seq_along(parts))) {
        Gk <- parts[[k]]$G
        lv[[k]] <- parts[[k]]$levels[code %% Gk + 1]
        code <- code %/% Gk
    }
    names(lv) <- vars
    list(codes = comb$codes, G = comb$G, levels = as.data.frame(lv, stringsAsFactors = FALSE), n = comb$n)
}

dmc_subset <- function(col, keep, x2) {
    if (is.null(col)) return(NULL)
    if (dm_is_codes(col)) stop("Codes change with the rows: compute them again.", call. = FALSE)
    if (dm_is_column(col)) return(dmc_on(col, x2))
    k <- if (is.logical(keep) && !anyNA(keep)) keep else !is.na(keep) & keep
    if (is.data.frame(col) || is.matrix(col)) col[k, , drop = FALSE] else col[k]
}

dmp_nested <- function(a, b) {
    if (identical(a, b)) return(TRUE)
    if (!dm_is_column(a$codes) && !dm_is_column(b$codes))
        return(.Call(C_dm_fe_nested, a$codes, b$codes, b$G))
    ct <- dmp_crosstab(a, b) # Ga x Gb
    all(diff(ct@p) == 1L)
}

dmp_center <- function(x) {
    if (inherits(x, "dm_remote")) {
        if (!isTRUE(getOption("design.matrix.remote_center", TRUE))) return(NULL)
        ctr <- dm_center(x)
        if (is.null(ctr)) return(NULL)
        x2 <- dm_shifted(x, ctr$shift)
        x2$cache <- new.env(parent = emptyenv())
        x2$id <- dm_tmp_name("dm")
        return(list(x = x2, T = ctr$T))
    }
    mt <- x$training_terms
    mf <- x$vars
    if (is.null(mt) || is.null(mf)) return(NULL)
    ctr <- dm_center_memory(mt, mf, x$contrasts)
    if (is.null(ctr)) return(NULL)
    ## (the columns as stored, as T was found: design.matrix#r6g7)
    attr(mf, "terms") <- mt
    for (nm in names(ctr$shift)) mf[[nm]] <- mf[[nm]] - ctr$shift[[nm]]
    list(x = design_matrix(mt, mf, contrasts.arg = x$contrasts), T = ctr$T)
}

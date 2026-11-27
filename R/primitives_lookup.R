## Lookup columns, keys and codes: primitives of the fitting core (see
## primitives.R for the column conventions), used by fits with fixed effects (absorb.R of xtx.fit).
## Each has an in-memory and a remote method; code that calls them does not
## branch on the backend.
##
##   xtx_c_lookup(values, key, levels = NULL)
##       the column whose value in a row is values[l], l the level of the
##       row's key: a lookup table joined on the key. Without `levels`, the
##       key is one column of integer codes 1..length(values) (or a
##       xtx_c_codes() grouping); with `levels` (a data frame, one row per
##       value, one column per key part), the key has one column per part
##       and a row's level is the row of `levels` with the same values. Keys
##       with no level (codes out of range, new values, NA) give NA (SQL
##       NULL). `values` can be a G x k matrix: k columns, from one lookup
##       table (one join remotely).
##       In memory: values[codes] (C), or values[match(keys)] (key: a vector,
##       or a list or data frame of vectors, one per part), as doubles.
##       Remotely: the lookup table is written to a temporary table at
##       once, and the column is a node (see primitives.R): each read that
##       names it joins the table (LEFT JOIN) to the source relation of the
##       read, once, after the nodes its key reads; integer `values`
##       (codes) give an integer column.
##   xtx_c_key(x, vars, data = NULL)
##       the values of variables (columns of the data, of any type) as a key
##       over the rows of design x: in memory a data frame (from `data`, or
##       from the model frame the design keeps); remotely one SQL column per
##       variable.
##   (is.na(key): xtx_c_present(), which takes the columns of a key of
##       any type)
##   xtx_c_dense_codes(key)
##       list(codes = xtx_c_codes(), levels = data frame): the levels of the
##       key that occur, sorted, and each row's level as an integer code. In
##       memory from the values (integer coding in C); remotely a SELECT
##       DISTINCT of the levels and a lookup of their codes (a node).
##   xtx_c_subset(col, keep, x2)
##       a column of x restricted to the rows of x2 = xtx_p_subset(x, keep): in
##       memory col[keep]; remotely the column carried over to x2 (xtx_c_on()).
##   xtx_p_nested(a, b)
##       for two xtx_c_codes() groupings: is each level of b in one level of
##       a? (In memory one pass in C; remotely from the crosstab.)
##   xtx_p_center(x)
##       the design with numeric variables shifted by their mean (in a
##       sample of rows), when that design is X T with T unit upper
##       triangular: list(x = the shifted design, T), or NULL. In memory the
##       design is rebuilt from its model frame, its columns as stored (a
##       column I(x^2) or log(x) is shifted on its own, not rebuilt from a
##       shifted x); remotely the variables' SQL (or the chunks' shift)
##       changes, as in xtx_center(), each variable on its own. The shifted
##       design is a new design: carry columns over with xtx_c_on().

xtx_c_lookup <- function(values, key, levels = NULL) {
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
    if (xtx_is_codes(key)) {
        levels <- key$levels %||% levels
        key <- key$codes
    }
    if (!xtx_is_column(key)) {
        if (is.null(levels)) {
            if (!is.integer(key)) key <- as.integer(key)
            ## (a matrix of values gives an n x k matrix, in one call)
            return(.Call(C_xtx_fe_lookup, values, key))
        }
        kd <- if (is.data.frame(key) || is.list(key)) as.data.frame(key, stringsAsFactors = FALSE)
              else data.frame(key, stringsAsFactors = FALSE)
        if (ncol(kd) != ncol(levels)) stop("`key` and `levels` have different numbers of columns.", call. = FALSE)
        for (j in seq_along(kd)) if (is.factor(kd[[j]])) kd[[j]] <- as.character(kd[[j]])
        i <- match(xtx_lookup_keys(kd), xtx_lookup_keys(levels))
        return(if (mat) values[i, , drop = FALSE] else values[i])
    }
    vm <- if (mat) values else matrix(values, ncol = 1L)
    vcols <- paste0("v", seq_len(ncol(vm)))
    ## (the key's frame, with its nodes: the lookup reads them)
    frame <- unclass(key)$frame
    con <- xtx_con(frame)
    q <- function(x) xtx_qi(frame, x)
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
    ## (the table of values belongs to the session of the connection, as
    ## any temporary table: it lives until the connection is closed)
    tn <- xtx_tmp_name("xtx_lookup")
    xtx_write_temp(con, tn, tab)
    stem <- xtx_node_stem("lk")
    col <- paste0(stem, if (mat) paste0("_", seq_along(vcols)))
    ## (a node: each read that names the column joins the table to its
    ## source, once)
    node <- xtx_node("lookup", stem, col, xtx_frame_nodes(frame), sprintf("l.%s AS %s", vcols, q(col)),
                     sprintf(" LEFT JOIN %s AS l ON %s", q(tn), on),
                     ref = list(table = tn, key = ksql, type = if (int) "integer" else "double"))
    xtx_column(xtx_frame_add(xtx_remote_unscaled(frame), node), sprintf("s.%s", q(col)))
}

## A data frame written to the temporary table `name`: registered with
## DuckDB and copied by one statement (DBI::dbWriteTable() first probes
## for the table, which costs much more than the copy of a lookup table).
xtx_write_temp <- function(con, name, df) {
    reg <- paste0(name, "_r")
    duckdb::duckdb_register(con, reg, df)
    on.exit(duckdb::duckdb_unregister(con, reg))
    DBI::dbExecute(con, sprintf("CREATE TEMPORARY TABLE %s AS SELECT * FROM %s",
                                DBI::dbQuoteIdentifier(con, name), DBI::dbQuoteIdentifier(con, reg)))
    invisible(name)
}

## One string per row of a data frame of key values (doubles to 17 digits).
xtx_lookup_keys <- function(d) {
    v <- lapply(d, function(x) if (is.double(x)) format(x, digits = 17, trim = TRUE) else as.character(x))
    if (length(v) == 1L) return(v[[1L]])
    do.call(paste, c(v, sep = "\r"))
}

xtx_c_key <- function(x, vars, data = NULL) {
    if (inherits(x, "xtx_remote")) {
        miss <- setdiff(vars, names(x$src$types))
        if (length(miss)) stop(sprintf("No column `%s` in the remote data.", miss[1L]), call. = FALSE)
        return(xtx_column(x, sprintf("s.%s", xtx_qi(x, vars)), type = unname(x$src$types[vars])))
    }
    d <- if (is.null(data)) x$vars else data
    miss <- setdiff(vars, names(d))
    if (length(miss)) stop(sprintf("No variable `%s` in the data.", miss[1L]), call. = FALSE)
    out <- as.data.frame(lapply(vars, function(v) d[[v]]), stringsAsFactors = FALSE)
    names(out) <- vars
    out
}

xtx_c_dense_codes <- function(key, raw = FALSE) {
    if (!xtx_is_column(key)) {
        cd <- xtx_fe_code_term(key, names(key))
        ## (in memory, also the rows per level: n)
        return(list(codes = xtx_c_codes(cd$codes, cd$G), levels = cd$levels, n = cd$n))
    }
    ## (the levels are read now)
    sql <- unclass(key)$sql
    frame <- xtx_frame_ready(unclass(key)$frame, sql)
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
    lv <- xtx_query(frame, paste(select, paste(fields, collapse = ", "), "FROM (",
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
        codes <- xtx_c_lookup(seq_len(nrow(lv)), key, lv)
        return(list(codes = xtx_c_codes(codes, nrow(lv)), levels = lv))
    }
    g <- xtx_c_codes(if (composite) xtx_column(unclass(key)$frame, tuple) else key, nrow(lv))
    g$levels <- if (composite) data.frame(tuple = tuple_values[ord]) else lv
    list(codes = g, levels = lv, n = counts[ord])
}

## ---- codes in memory (xtx_c_dense_codes()) ----

## Integer codes 1..G of a vector without missing values, its levels (the
## values that occur, sorted, or in level order for a factor) and the rows
## per level (n), in one pass in C (xtx_fe_code() in src/fe.c: direct
## indexing for integer values with a compact range, hashing otherwise;
## the levels are then sorted as sort(unique(v), method = "radix") sorts
## them). Other types, and strings that occur in several encodings, as
## sort(unique()) and match().
xtx_fe_code <- function(v) {
    if (is.factor(v)) {
        r <- .Call(C_xtx_fe_code, as.integer(v), nlevels(v) + 1)
        if (!is.null(r)) return(list(codes = r[[1L]], G = r[[2L]], levels = levels(v)[r[[3L]]], n = r[[4L]]))
        cd <- xtx_fe_compact(as.integer(v), nlevels(v))
        return(list(codes = cd$codes, G = cd$G, levels = levels(v)[cd$used]))
    }
    if (is.logical(v)) v <- as.integer(v)
    if (is.numeric(v) || (is.character(v) && is.null(attributes(v)))) {
        r <- .Call(C_xtx_fe_code, v, max(4 * length(v), 1e6))
        if (!is.null(r)) {
            if (length(r) == 4L) return(list(codes = r[[1L]], G = r[[2L]], levels = r[[3L]], n = r[[4L]]))
            ## hashed: levels in order of first appearance, sorted here
            ## (as plain values: not for a class with its own order)
            u <- r[[3L]]
            if (is.null(oldClass(v)) && (!is.character(u) || !anyDuplicated(u))) {
                o <- order(u, method = "radix")
                rank <- integer(length(o))
                rank[o] <- seq_along(o)
                return(list(codes = .Call(C_xtx_fe_recode, r[[1L]], rank), G = r[[2L]], levels = u[o],
                            n = r[[4L]][o]))
            }
        }
    }
    u <- sort(unique(v), method = "radix")
    list(codes = match(v, u), G = length(u), levels = u)
}

## Renumber codes so that the levels that occur are 1..G.
xtx_fe_compact <- function(cd, G) {
    used <- tabulate(cd, nbins = G) > 0L
    if (all(used)) return(list(codes = cd, G = G, used = used))
    map <- cumsum(used)
    list(codes = map[cd], G = sum(used), used = used)
}

## Codes of an FE (one or several variables, columns of `d`), with its
## levels as a data frame with one column per variable.
xtx_fe_code_term <- function(d, vars) {
    if (length(vars) == 1L) {
        cd <- xtx_fe_code(d[[vars]])
        lv <- list(cd$levels)
        names(lv) <- vars
        return(list(codes = cd$codes, G = cd$G, levels = as.data.frame(lv, stringsAsFactors = FALSE), n = cd$n))
    }
    parts <- lapply(vars, function(v) xtx_fe_code(d[[v]]))
    cd <- parts[[1L]]$codes
    for (k in seq_along(parts)[-1L]) cd <- (cd - 1) * parts[[k]]$G + parts[[k]]$codes
    comb <- xtx_fe_code(cd)
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

xtx_c_subset <- function(col, keep, x2) {
    if (is.null(col)) return(NULL)
    if (xtx_is_codes(col)) stop("Codes change with the rows: compute them again.", call. = FALSE)
    if (xtx_is_column(col)) return(xtx_c_on(col, x2))
    k <- if (is.logical(keep) && !anyNA(keep)) keep else !is.na(keep) & keep
    if (is.data.frame(col) || is.matrix(col)) col[k, , drop = FALSE] else col[k]
}

xtx_p_nested <- function(a, b) {
    if (identical(a, b)) return(TRUE)
    if (!xtx_is_column(a$codes) && !xtx_is_column(b$codes))
        return(.Call(C_xtx_fe_nested, a$codes, b$codes, b$G))
    ct <- xtx_p_crosstab(a, b) # Ga x Gb
    all(diff(ct@p) == 1L)
}

xtx_p_center <- function(x) {
    if (inherits(x, "xtx_remote")) {
        if (!isTRUE(getOption("xtx.remote_center", TRUE))) return(NULL)
        ctr <- xtx_center(x)
        if (is.null(ctr)) return(NULL)
        x2 <- xtx_shifted(x, ctr$shift)
        x2$cache <- new.env(parent = emptyenv())
        x2$from <- c(x$id, x$from)
        x2$id <- xtx_tmp_name("xtx")
        return(list(x = x2, T = ctr$T))
    }
    mt <- x$training_terms
    mf <- x$vars
    if (is.null(mt) || is.null(mf)) return(NULL)
    ctr <- xtx_center_memory(mt, mf, x$contrasts)
    if (is.null(ctr)) return(NULL)
    ## (the columns as stored, as T was found: design.matrix#r6g7)
    attr(mf, "terms") <- mt
    for (nm in names(ctr$shift)) mf[[nm]] <- mf[[nm]] - ctr$shift[[nm]]
    list(x = design_matrix(mt, mf, contrasts.arg = x$contrasts), T = ctr$T)
}

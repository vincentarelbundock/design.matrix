## The chunked engine of remote designs (see remote.R): batches of rows are
## read through DBI and handed to the in-memory design.
##
## DuckDB's R client materializes a whole result before dbFetch() slices it,
## so batches are separate queries over rowid ranges: of the table itself
## when the source is a base table, otherwise of a temporary table holding
## the needed columns of the rows used (which DuckDB can spill to disk).

## The query reading columns of the rows used, optionally in a rowid range
## of `from` (a table).
dm_chunk_sql <- function(D, cols, range = NULL, from = NULL) {
    sel <- if (length(cols)) paste0("s.", dm_qi(D, cols), collapse = ", ") else "1 AS dm_one"
    src <- from %||% sprintf("(%s)", D$src$select)
    where <- D$filter
    if (!is.null(range))
        where <- sprintf("(%s) AND s.rowid >= %.0f AND s.rowid < %.0f", where, range[1L], range[2L])
    sprintf("SELECT %s FROM %s AS s WHERE %s", sel, src, where)
}

## Call FUN(chunk) on each non-empty batch of the rows used.
dm_chunks <- function(D, cols, FUN) {
    cs <- D$chunk_size
    if (!is.null(D$src$base_table)) {
        from <- D$src$base_table
        r <- dm_query(D, sprintf("SELECT MIN(rowid) AS lo, MAX(rowid) AS hi FROM %s", from))
        filtered <- TRUE
    } else {
        from <- dm_qi(D, dm_tmp_name("dm_rows"))
        dm_execute(D, sprintf("CREATE TEMPORARY TABLE %s AS %s", from, dm_chunk_sql(D, cols)))
        on.exit(dm_execute(D, sprintf("DROP TABLE IF EXISTS %s", from)), add = TRUE)
        r <- dm_query(D, sprintf("SELECT MIN(rowid) AS lo, MAX(rowid) AS hi FROM %s", from))
        filtered <- FALSE
    }
    lo <- as.numeric(r$lo)
    hi <- as.numeric(r$hi)
    if (is.na(lo)) return(invisible(NULL))
    a <- lo
    while (a <= hi) {
        rng <- c(a, a + cs)
        sql <- if (filtered) {
            dm_chunk_sql(D, cols, range = rng, from = from)
        } else {
            sprintf("SELECT * FROM %s AS s WHERE s.rowid >= %.0f AND s.rowid < %.0f", from, rng[1L], rng[2L])
        }
        chunk <- dm_query(D, sql)
        if (nrow(chunk)) FUN(chunk)
        a <- a + cs
    }
    invisible(NULL)
}

## The in-memory design of a batch: rows with a missing model variable
## (response included) are dropped, as na.omit() drops them, unless the
## design keeps them (na.pass).
## Returns list(design, keep = rows of the chunk kept, y = response).
dm_chunk_design <- function(D, chunk) {
    plan <- D$plan
    na <- if (identical(D$na, "pass")) stats::na.pass else stats::na.omit
    mf <- design_frame(D$terms_full, chunk, xlev = plan$xlevels, na.action = na)
    ## variables of another type than in training: the error of predict.lm()
    if (!is.null(plan$data_classes)) stats::.checkMFClasses(plan$data_classes, mf)
    ## numeric variables shifted for a better conditioned fit (qr.R)
    for (nm in names(D$shift)) mf[[nm]] <- mf[[nm]] - D$shift[[nm]]
    design <- design_matrix(plan$terms, mf, contrasts.arg = plan$contrasts)
    if (isTRUE(plan$drop_intercept)) design <- dm_drop_intercept(design)
    omit <- attr(mf, "na.action")
    keep <- seq_len(nrow(chunk))
    if (length(omit)) keep <- keep[-omit]
    y <- if (is.null(D$response)) NULL else unname(stats::model.response(mf))
    list(design = design, keep = keep, y = y)
}

dm_chunk_y <- function(D, yy, chunk, cd) {
    if (yy$response) {
        if (!is.numeric(cd$y)) stop("The response must be numeric.", call. = FALSE)
        return(matrix(as.double(cd$y), ncol = 1L))
    }
    Y <- vapply(yy$cols, function(col) {
        x <- chunk[[col]][cd$keep]
        if (!is.numeric(x) && !is.logical(x)) stop(sprintf("Column `%s` is not numeric.", col), call. = FALSE)
        as.double(x)
    }, numeric(length(cd$keep)))
    matrix(Y, ncol = length(yy$cols))
}

dm_chunked_moments <- function(D, yy, weights, gram = TRUE) {
    p <- D$p
    K <- if (is.null(yy)) 0L else length(yy$cols)
    ycols <- if (K && !yy$response) yy$cols
    xtx <- matrix(0, p, p)
    xty <- matrix(0, p, K)
    yty <- matrix(0, K, K)
    n <- 0
    sumw <- 0
    dm_chunks(D, unique(c(dm_remote_inputs(D), ycols, weights)), function(chunk) {
        cd <- dm_chunk_design(D, chunk)
        X <- cd$design
        w <- if (is.null(weights)) NULL else as.double(chunk[[weights]][cd$keep])
        if (gram) xtx <<- xtx + if (is.null(w)) crossprod(X) else crossprod(X, w * X)
        if (K) {
            Y <- dm_chunk_y(D, yy, chunk, cd)
            WY <- if (is.null(w)) Y else w * Y
            xty <<- xty + crossprod(X, WY)
            yty <<- yty + crossprod(Y, WY)
        }
        n <<- n + length(cd$keep)
        sumw <<- sumw + if (is.null(w)) length(cd$keep) else sum(w)
    })
    D$cache$n <- n
    yn <- if (K) yy$cols
    dimnames(xtx) <- list(D$colnames, D$colnames)
    dimnames(xty) <- list(D$colnames, yn)
    dimnames(yty) <- list(yn, yn)
    if (!K) xty <- yty <- NULL
    list(xtx = xtx, xty = xty, yty = yty, n = n, sumw = sumw)
}

dm_chunked_colsums <- function(D, group, weights, mean) {
    p <- D$p
    gvals <- NULL
    S <- matrix(0, 0L, p)
    counts <- numeric()
    total <- numeric(p)
    n <- 0
    dm_chunks(D, unique(c(dm_remote_inputs(D), group, weights)), function(chunk) {
        cd <- dm_chunk_design(D, chunk)
        X <- cd$design
        if (!is.null(weights)) X <- as.double(chunk[[weights]][cd$keep]) * X
        n <<- n + length(cd$keep)
        if (is.null(group)) {
            total <<- total + colSums(X)
            return()
        }
        g <- chunk[[group]][cd$keep]
        u <- unique(g)
        new <- u[is.na(match(u, gvals))]
        if (length(new)) {
            gvals <<- if (is.null(gvals)) new else c(gvals, new)
            S <<- rbind(S, matrix(0, length(new), p))
            counts <<- c(counts, numeric(length(new)))
        }
        codes <- match(g, gvals)
        s <- suppressWarnings(colSums(X, group = codes))
        idx <- as.integer(rownames(s))
        S[idx, ] <<- S[idx, , drop = FALSE] + s
        counts <<- counts + tabulate(codes, nbins = length(gvals))
    })
    D$cache$n <- n
    if (is.null(group)) {
        names(total) <- D$colnames
        return(if (mean) total / n else total)
    }
    if (anyNA(gvals)) warning("missing values for 'group'", call. = FALSE)
    ord <- order(gvals, na.last = TRUE)
    out <- S[ord, , drop = FALSE]
    if (mean) out <- out / counts[ord]
    dimnames(out) <- list(as.character(gvals[ord]), D$colnames)
    out
}

dm_chunked_predict <- function(D, B, name, outcols, keep, overwrite) {
    con <- dm_con(D)
    temporary <- is.null(name)
    if (temporary) name <- dm_tmp_name("dm_pred")
    if (!temporary && overwrite && DBI::dbExistsTable(con, name)) DBI::dbRemoveTable(con, name)
    created <- FALSE
    write <- function(df) {
        if (!created) {
            DBI::dbWriteTable(con, name, df, temporary = temporary)
            created <<- TRUE
        } else {
            DBI::dbAppendTable(con, name, df)
        }
    }
    dm_chunks(D, unique(c(dm_remote_inputs(D), keep)), function(chunk) {
        cd <- dm_chunk_design(D, chunk)
        P <- cd$design %*% B
        out <- chunk[cd$keep, keep, drop = FALSE]
        for (j in seq_along(outcols)) out[[outcols[j]]] <- P[, j]
        attr(out, "row.names") <- c(NA_integer_, -length(cd$keep))
        write(out)
    })
    if (!created) {
        out <- D$src$proto[keep]
        for (j in seq_along(outcols)) out[[outcols[j]]] <- numeric()
        write(out)
    }
    name
}

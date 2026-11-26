## Compressed remote products: one GROUP BY over every input column of the
## design (see remote.R).
##
## Rows with the same values in every input column of the design give the
## same design row. When the remote data has few such patterns (a design of
## factors and discrete numeric variables), one query returns, for each
## pattern, its input values, its number of rows and the sums of the weights
## and of the (weighted) responses; the design of the patterns is built in
## memory through the plan, and
##
##     X'WX = Xp' diag(sum w) Xp,   X'Wy = Xp' (sum w y),   y'Wy = sum w y^2,
##
## exactly (up to rounding). This is the compression of design_lm() (see
## dmfit_suffstats() in lm.R), done in the database.
##
## The rows used are those of the design (na.omit only: with na.pass, rows
## with missing values would put NA where the dense computation does, which
## the patterns cannot reproduce, so compression is skipped). A variable
## evaluated in R that gives NA on a pattern drops the whole pattern, as
## na.omit() drops each of its rows.
##
## Option design.matrix.remote_compress: "auto" (default) tries compression
## when the first rows have few distinct patterns, and gives up when the
## patterns exceed min(n / 2, getOption("design.matrix.remote_compress_max",
## 1e6)); TRUE always compresses (unless every row is distinct); FALSE never.

dm_compress_mode <- function(mode = getOption("design.matrix.remote_compress", "auto")) {
    if (!(isTRUE(mode) || isFALSE(mode) || identical(mode, "auto")))
        stop("Option `design.matrix.remote_compress` must be \"auto\", TRUE or FALSE.", call. = FALSE)
    mode
}

## Input columns of the design (not those of the response alone).
dm_compress_keys <- function(D) {
    pv <- dm_terms_vars(D$plan$terms)$names
    unique(unlist(lapply(D$vars[pv], function(v) v$inputs)))
}

## Per-pattern sums, or NULL when compression does not pay.
##   ys: SQL expressions over `s` for the responses (named), or empty;
##   w: SQL expression over `s` for the weights, or NULL;
##   group: a column added to the patterns;
##   within: also the weighted sum of squares of the first response around
##     its pattern means (for design_lm()).
## Returns list(P = input values, n, sw, sy (patterns x K), yy (K x K),
## within, rows).
dm_remote_patterns <- function(D, ys = character(), w = NULL, group = NULL, within = FALSE,
                               mode = dm_compress_mode()) {
    if (isFALSE(mode) || !identical(D$na, "omit")) return(NULL)
    keys <- unique(c(dm_compress_keys(D), group))
    n <- D$cache$n %||% D$n
    if (n < 2) return(NULL)
    src <- sprintf("(%s) AS s", D$src$select)
    qk <- if (length(keys)) paste0("s.", dm_qi(D, keys)) else character()
    if (identical(mode, "auto") && length(keys)) {
        ## few distinct patterns among the first rows?
        probe <- 100000
        pr <- dm_query(D, sprintf("SELECT COUNT(*) AS g, SUM(r) AS r FROM (SELECT COUNT(*) AS r FROM (SELECT %s FROM %s WHERE %s LIMIT %d) AS s GROUP BY %s)",
                                  paste(qk, collapse = ", "), src, D$filter, probe,
                                  paste(dm_qi(D, keys), collapse = ", ")))
        if (as.numeric(pr$g) > 0.1 * max(as.numeric(pr$r), 1000)) return(NULL)
    }
    cap <- if (isTRUE(mode)) n else min(n / 2, getOption("design.matrix.remote_compress_max", 1e6))
    K <- length(ys)
    wf <- if (is.null(w)) "" else " * __dm_w"
    inner <- c(qk, if (K) sprintf("%s AS __dm_y%d", ys, seq_len(K)), if (!is.null(w)) sprintf("%s AS __dm_w", w))
    if (!length(inner)) inner <- "1 AS __dm_one"
    idx <- expand.grid(k = seq_len(K), l = seq_len(K))
    idx <- idx[idx$k <= idx$l, , drop = FALSE]
    aggs <- c("COUNT(*) AS __dm_n",
              sprintf("SUM(CAST(1 AS DOUBLE)%s) AS __dm_sw", wf),
              if (K) sprintf("SUM(__dm_y%d%s) AS __dm_sy%d", seq_len(K), wf, seq_len(K)),
              if (K) sprintf("SUM(__dm_y%d * __dm_y%d%s) AS __dm_yy%d_%d", idx$k, idx$l, wf, idx$k, idx$l),
              if (within && is.null(w)) "var_pop(__dm_y1) * COUNT(*) AS __dm_within",
              if (!is.null(w)) c("SUM(CASE WHEN __dm_w <> 0 THEN 1 ELSE 0 END) AS __dm_nz",
                                 "SUM(CASE WHEN __dm_w <> 0 THEN ln(__dm_w) ELSE 0 END) AS __dm_slw"))
    kq <- if (length(keys)) dm_qi(D, keys) else character()
    sql <- sprintf("SELECT %s FROM (SELECT %s FROM %s WHERE %s) AS b%s%s",
                   paste(c(kq, aggs), collapse = ", "), paste(inner, collapse = ", "), src, D$filter,
                   if (length(keys)) " GROUP BY ALL" else "",
                   if (length(keys)) sprintf(" LIMIT %.0f", cap + 1) else "")
    res <- dm_query(D, sql)
    G <- nrow(res)
    if (G > cap || (G == n && length(keys))) return(NULL)
    num <- function(x) {
        x <- as.double(x)
        x[is.na(x)] <- 0
        x
    }
    out <- list(P = res[keys], n = num(res$`__dm_n`), sw = num(res$`__dm_sw`), rows = sum(num(res$`__dm_n`)))
    out$nz <- if (is.null(w)) out$n else num(res$`__dm_nz`)
    out$slw <- if (is.null(w)) numeric(G) else num(res$`__dm_slw`)
    out$sy <- matrix(vapply(seq_len(K), function(k) num(res[[paste0("__dm_sy", k)]]), numeric(G)), G, K)
    out$yy <- array(0, c(G, K, K))
    for (r in seq_len(nrow(idx))) {
        k <- idx$k[r]
        l <- idx$l[r]
        v <- num(res[[sprintf("__dm_yy%d_%d", k, l)]])
        out$yy[, k, l] <- v
        out$yy[, l, k] <- v
    }
    if (within) {
        out$within <- if (is.null(w)) {
            num(res$`__dm_within`)
        } else {
            ## weighted: a second pass, around the pattern means
            on <- if (length(keys)) {
                paste(sprintf("b.%s IS NOT DISTINCT FROM p.%s", kq, kq), collapse = " AND ")
            } else "TRUE"
            q <- sprintf(paste("WITH b AS (SELECT %s FROM %s WHERE %s),",
                               "p AS (SELECT %s FROM b%s)",
                               "SELECT %s FROM b JOIN p ON %s WHERE b.__dm_w <> 0%s"),
                         paste(inner, collapse = ", "), src, D$filter,
                         paste(c(kq, "SUM(__dm_w * __dm_y1) / SUM(__dm_w) AS __dm_m"), collapse = ", "),
                         if (length(keys)) " GROUP BY ALL" else "",
                         paste(c(paste0("b.", kq), "SUM(b.__dm_w * pow(b.__dm_y1 - p.__dm_m, 2)) AS __dm_within"),
                               collapse = ", "),
                         on, if (length(keys)) " GROUP BY ALL" else "")
            wr <- dm_query(D, q)
            ## in the order of the patterns
            mm <- if (length(keys)) dm_match_rows(res[keys], wr[keys]) else 1L
            v <- num(wr$`__dm_within`)[mm]
            v[is.na(v)] <- 0
            v
        }
    }
    out
}

## For each row of `a`, the matching row of `b` (same values, NA matching NA).
dm_match_rows <- function(a, b) {
    key <- function(d) do.call(paste, c(lapply(d, function(x) {
        if (is.double(x)) sprintf("%a", x) else as.character(x)
    }), sep = "\r"))
    match(key(a), key(b))
}

## The in-memory design of the patterns: list(design, keep = patterns kept).
dm_pattern_design <- function(D, P) {
    plan <- D$plan
    if (!ncol(P)) P <- structure(list(), row.names = c(NA_integer_, -max(1L, nrow(P))), class = "data.frame")
    mf <- design_frame(plan$terms, P, xlev = plan$xlevels, na.action = stats::na.omit)
    if (!is.null(plan$data_classes)) stats::.checkMFClasses(plan$data_classes, mf)
    design <- design_matrix(plan$terms, mf, contrasts.arg = plan$contrasts)
    if (isTRUE(plan$drop_intercept)) design <- dm_drop_intercept(design)
    omit <- attr(mf, "na.action")
    keep <- seq_len(nrow(P))
    if (length(omit)) keep <- keep[-omit]
    list(design = design, keep = keep)
}

## Responses as SQL over `s`: the response of the formula, or columns.
dm_compress_ys <- function(D, yy) {
    if (is.null(yy)) return(character())
    if (yy$response) {
        v <- D$vars[[yy$cols]]
        if (v$kind != "numeric") return(NULL)
        return(v$sql)
    }
    for (col in yy$cols) {
        if (!D$src$types[[col]] %in% c("double", "integer", "logical"))
            stop(sprintf("Column `%s` is not numeric.", col), call. = FALSE)
    }
    sprintf("CAST(s.%s AS DOUBLE)", dm_qi(D, yy$cols))
}

## design_moments() from patterns, or NULL.
dm_compressed_moments <- function(D, yy, weights, gram = TRUE) {
    ys <- dm_compress_ys(D, yy)
    if (is.null(ys)) return(NULL)
    ## missing responses or weights give NA (the uncompressed path does that)
    na_cols <- c(if (!is.null(yy) && !yy$response) yy$cols, weights)
    if (length(na_cols) && any(dm_sql_count_na(D, na_cols) > 0)) return(NULL)
    w <- if (is.null(weights)) NULL else sprintf("CAST(s.%s AS DOUBLE)", dm_qi(D, weights))
    pt <- dm_remote_patterns(D, ys, w)
    if (is.null(pt)) return(NULL)
    pd <- dm_pattern_design(D, pt$P[dm_compress_keys(D)])
    k <- pd$keep
    X <- pd$design
    sw <- pt$sw[k]
    xtx <- if (gram) crossprod(X, sw * X) else NULL
    K <- length(ys)
    yn <- if (K) yy$cols
    xty <- yty <- NULL
    if (K) {
        xty <- crossprod(X, pt$sy[k, , drop = FALSE])
        yty <- matrix(colSums(matrix(pt$yy[k, , , drop = FALSE], length(k), K * K)), K, K)
        dimnames(xty) <- list(D$colnames, yn)
        dimnames(yty) <- list(yn, yn)
    }
    if (gram) dimnames(xtx) <- list(D$colnames, D$colnames)
    n <- sum(pt$n[k])
    D$cache$n <- n
    D$cache$patterns <- length(k)
    list(xtx = xtx, xty = xty, yty = yty, n = n, sumw = sum(sw))
}

## colSums()/colMeans() from patterns (with the group in the key), or NULL.
dm_compressed_colsums <- function(D, group, weights, mean) {
    if (!is.null(weights) && dm_sql_count_na(D, weights) > 0) return(NULL)
    w <- if (is.null(weights)) NULL else sprintf("CAST(s.%s AS DOUBLE)", dm_qi(D, weights))
    pt <- dm_remote_patterns(D, character(), w, group = group)
    if (is.null(pt)) return(NULL)
    pd <- dm_pattern_design(D, pt$P[dm_compress_keys(D)])
    k <- pd$keep
    X <- pt$sw[k] * pd$design
    cnt <- pt$n[k]
    if (is.null(group)) {
        out <- colSums(X)
        names(out) <- D$colnames
        return(if (mean) out / sum(cnt) else out)
    }
    g <- pt$P[[group]][k]
    out <- colSums(X, group = g)
    if (mean) out <- out / drop(rowsum(cnt, g, reorder = TRUE))
    dimnames(out) <- list(rownames(out), D$colnames)
    out
}

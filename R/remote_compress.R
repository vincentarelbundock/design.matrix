## Compressed remote products: the scan of the SQL engine (remote_sql.R)
## with the numeric variables of the design in its GROUP BY.
##
## Rows with the same factor codes and the same values of the numeric
## variables give the same design row. When the remote data has few such
## cells (a design of factors and discrete numeric variables), the one scan
## (dm_sql_scan()) groups by every factor of the design and every numeric
## variable, and sums, per cell, the rows, the weights and the (weighted)
## responses; the design of the cells is X[c, cols_t] = N_t(c) C_t[F_t(c), ]
## with the engine's own C_t (dm_sql_cell_design()), and
##
##     X'WX = Xc' diag(sum w) Xc,   X'Wy = Xc' (sum w y),   y'Wy = sum w y^2,
##
## exactly (up to rounding). This is the compression of design_lm() (see
## dmfit_suffstats() in lm.R), done in the database.
##
## Only designs with na.omit are compressed (with na.pass, rows with
## missing values put NA where the dense computation does).
##
## Option design.matrix.remote_compress: "auto" (default) compresses by the
## rule of the operation (dm_compress_rule()) and dm_compress_worth(), and
## gives up when the scan returns more than
## min(n / 2, getOption("design.matrix.remote_compress_max", 1e6)) cells;
## TRUE always compresses (unless every row is distinct); FALSE never.

dm_compress_mode <- function(mode = getOption("design.matrix.remote_compress", "auto")) {
    if (!(isTRUE(mode) || isFALSE(mode) || identical(mode, "auto")))
        stop("Option `design.matrix.remote_compress` must be \"auto\", TRUE or FALSE.", call. = FALSE)
    mode
}

## Is compression worth it ("auto")? Two steps, before any GROUP BY:
##  1. continuous keys: a DOUBLE key column whose approximate number of
##     distinct values, in a sample of rows, exceeds half the sampled rows
##     makes most rows distinct patterns: skip. (The sample is a system
##     sample, which reads only some blocks, of 0.1% of the rows, or a
##     reservoir sample of 10,000 rows for at most a million rows.) Other
##     types (integers, text, ENUM, logical, dates) go to step 2;
##  2. an estimate of the number of patterns over all the rows used, by
##     approx_count_distinct() of a hash of the key columns (missing
##     values count as values): compress when it is at most `ratio` times
##     the number of rows and at most
##     getOption("design.matrix.remote_compress_max", 1e6).
## The real GROUP BY gives up as well when it returns more patterns than
## that (dm_sql_cell_design()). The estimate is asked once for the rows of
## a design and a set of keys; the operations differ only in `ratio`.
dm_compress_worth <- function(D, keys, n, group = character(),
                              ratio = getOption("design.matrix.remote_compress_ratio", 0.1)) {
    ck <- paste(c(D$filter, keys, group), collapse = "\r")
    if (is.null(D$cache$worth[[ck]])) D$cache$worth[[ck]] <- dm_compress_estimate(D, keys, n, group)
    dm_compress_ok(D$cache$worth[[ck]], ratio)
}

dm_compress_ok <- function(e, ratio) {
    e[["est"]] <= ratio * e[["rows"]] && e[["est"]] <= getOption("design.matrix.remote_compress_max", 1e6)
}

## (uncached; bench/compression/run.R times it)
dm_compress_worth_ask <- function(D, keys, n, group = character(),
                                  ratio = getOption("design.matrix.remote_compress_ratio", 0.1)) {
    dm_compress_ok(dm_compress_estimate(D, keys, n, group), ratio)
}

## c(est, rows): the estimated number of patterns (Inf for continuous keys,
## step 1) and the number of rows.
dm_compress_estimate <- function(D, keys, n, group = character()) {
    qk <- c(if (length(keys)) paste0("s.", dm_qi(D, keys)), sprintf("(%s)", group))
    dbl <- keys[D$src$types[keys] == "double"]
    if (length(dbl)) {
        from <- if (!is.null(D$src$base_table)) sprintf("%s AS s", D$src$base_table)
                else sprintf("(%s) AS s", D$src$select)
        smp <- if (n > 1e6) "TABLESAMPLE 0.1% (system)" else "TABLESAMPLE 10000 ROWS"
        qd <- paste0("s.", dm_qi(D, dbl))
        res <- dm_query(D, sprintf("SELECT COUNT(*) AS n, %s FROM %s %s WHERE %s",
                                   paste(sprintf("approx_count_distinct(%s) AS d%d", qd, seq_along(qd)), collapse = ", "),
                                   from, smp, D$filter))
        m <- as.numeric(res$n)
        dd <- vapply(seq_along(dbl), function(j) as.numeric(res[[paste0("d", j)]]), 0)
        if (m > 0 && any(dd > 0.5 * m)) return(c(est = Inf, rows = n))
    }
    res <- dm_query(D, sprintf("SELECT COUNT(*) AS n, approx_count_distinct(hash(%s)) AS est FROM (%s) AS s WHERE %s",
                               paste(qk, collapse = ", "), D$src$select, D$filter))
    est <- as.numeric(res$est)
    D$cache$compress_estimate <- est
    c(est = est, rows = as.numeric(res$n))
}

## Input columns of the design (not those of the response alone).
dm_compress_keys <- function(D) {
    tt <- D$plan$terms
    pv <- dm_terms_vars(tt)$names
    ## offset() terms are not in the design (the fit subtracts them from
    ## the response)
    off <- attr(tt, "offset")
    if (length(off)) pv <- pv[-off]
    unique(unlist(lapply(D$vars[pv], function(v) v$inputs)))
}

## The rule of "auto" for each operation (bench/compression/ANALYSIS.md,
## design.matrix#etg8): NULL, never compress; else list(numeric, ratio):
## compress designs with at least `numeric` numeric columns
## (dm_compress_numeric(); fewer skip the decision queries entirely) whose
## estimated patterns are at most `ratio` times the rows. The option
## design.matrix.remote_compress_ratio, when set, caps every ratio (the
## smaller of the two).
##  - "colsums" (colSums(), colMeans(), grouped or not): never. Compressed,
##    the GROUP BY splits the cells (and the groups) by every value of the
##    numeric variables, and there are no products to save: it lost almost
##    everywhere at n >= 1e6 (up to 50x).
##  - "products" (crossprod(), X'y): the uncompressed scan already groups by
##    the factor cells; compression only saves the per-row products of the
##    numeric columns, about m (m + 1) / 2 of them. It won 2.7-10x for 15
##    numeric columns up to ratio 1e-2 and lost for 2. The cutoff of 8 is
##    PROVISIONAL (between 2 and 15 is unmeasured): retune it on a quiet
##    machine with bench/compression/run.R.
##  - the fits (dmp_patterns(), dm_sql_patterns()):
##    - "glm" (design_glm()): IRLS reads the rows once per iteration
##      uncompressed, once in all compressed: it gained 2-4x up to ratio
##      1e-2, still 1.1-3x at 0.1, and broke even or lost at 0.5. Ratio 0.2;
##    - "lm" (unweighted design_lm()): one pass either way; compression
##      gained a modest 1.2x (2-9x wide) up to 1e-3 and lost at 1e-2 for
##      discrete designs at 5e6. Ratio 1e-3;
##    - "wlm" (weighted design_lm()): the weighted sum of squares within the
##      cells takes a second pass, and it lost except for designs with many
##      numeric columns (wide: 1.9-7.8x up to 1e-3). Only with at least 8
##      numeric columns (the cutoff of "products") at ratio 1e-3.
##    These ratios are PROVISIONAL, from a noisy machine: confirm them on a
##    quiet one with bench/compression/run.R.
dm_compress_rule <- function(op) {
    user <- getOption("design.matrix.remote_compress_ratio")
    cap <- function(r) if (is.null(user)) r else min(user, r)
    switch(op,
           colsums = NULL,
           products = list(numeric = 8L, ratio = cap(1e-2)),
           glm = list(numeric = 0L, ratio = cap(0.2)),
           lm = list(numeric = 0L, ratio = cap(1e-3)),
           wlm = list(numeric = 8L, ratio = cap(1e-3)),
           stop("unknown operation: ", op))
}

## The number of distinct numeric columns of the products: the distinct
## products of numeric variables in the SQL terms (x, x:f and x:g count
## once; f alone not at all).
dm_compress_numeric <- function(D) {
    st <- Filter(function(t) length(t$cols) > 0L && length(t$M) > 0L, D$sql_terms)
    length(unique(vapply(st, function(t) paste(sort(t$M), collapse = "*"), "")))
}

## The cap on the number of cells of compressed products of D, or NULL when
## they are not compressed. `group`: SQL expressions that also split the
## cells (the groups of colSums(), the keys of dmp_patterns()). `op`: the
## operation, for the rule of "auto" (dm_compress_rule()).
dm_compress_cap <- function(D, mode = dm_compress_mode(), group = character(),
                            op = c("lm", "wlm", "glm", "products", "colsums")) {
    op <- match.arg(op)
    if (isFALSE(mode) || !identical(D$na, "omit")) return(NULL)
    n <- D$cache$n %||% D$n
    if (n < 2) return(NULL)
    if (isTRUE(mode)) return(n)
    rule <- dm_compress_rule(op)
    if (is.null(rule) || dm_compress_numeric(D) < rule$numeric) return(NULL)
    cols <- dm_compress_keys(D)
    if (length(c(cols, group)) && !dm_compress_worth(D, cols, n, group, rule$ratio)) return(NULL)
    min(n / 2, getOption("design.matrix.remote_compress_max", 1e6))
}

## The keys of the cells of the scan: S, the factors of the design (every
## sum is within their combined levels), and `by`, the numeric columns of
## the base query (dm_sql_base()). Nothing without `cap` (not compressed).
dm_sql_cells <- function(D, cap) {
    if (is.null(cap)) return(list(S = integer(), by = character()))
    st <- Filter(function(t) length(t$cols) > 0L, D$sql_terms)
    list(S = sort(unique(as.integer(unlist(lapply(st, function(t) t$S))))),
         by = unique(unlist(lapply(st, function(t) t$M))))
}

## The design of the cells (the rows of `res`, a result of the scan), or
## NULL when they are more than `cap`, or as many as the rows.
dm_sql_cell_design <- function(D, res, cap) {
    G <- nrow(res)
    if (G > cap || (G == (D$cache$n %||% D$n) && G > 1L)) return(NULL)
    X <- matrix(0, G, D$p, dimnames = list(NULL, D$colnames))
    for (Tt in D$sql_terms) {
        if (!length(Tt$cols)) next
        Ct <- Tt$C[dm_sql_level(res, Tt$S, Tt$strides) + 1, , drop = FALSE]
        for (m in sort(Tt$M)) Ct <- as.double(res[[m]]) * Ct
        X[, Tt$cols] <- X[, Tt$cols] + Ct
    }
    attr(X, "assign") <- D$assign
    X
}

## Compressed least squares for dmp_patterns(): one scan of the cells
## (with the SQL expressions `keys` in the key, columns __dm_g1, ...),
## summing the weights `w` and the weighted response `z`, and the weighted
## sum of squares of z around the means of the cells. Returns list(x = the
## design of the cells, n, sw, sy, within, keys), or NULL when compression
## does not pay. `fit`: "lm" or "glm", for the rule of "auto" (a weighted
## "lm" is "wlm").
dm_sql_patterns <- function(D, z, w = NULL, keys = NULL, mode = dm_compress_mode(), fit = "lm") {
    op <- if (identical(fit, "glm")) "glm" else if (is.null(w)) "lm" else "wlm"
    cap <- dm_compress_cap(D, mode, keys, op = op)
    if (is.null(cap)) return(NULL)
    ck <- dm_sql_cells(D, cap)
    gn <- if (length(keys)) paste0("__dm_g", seq_along(keys)) else character()
    base <- dm_sql_base(D, c(`__dm_y1` = z, `__dm_w` = w, if (length(gn)) stats::setNames(keys, gn)), asis = gn)
    sc <- dm_sql_scan(D, base, vapply(D$vars, function(v) max(1L, length(v$levels)), 1L),
                      by = c(ck$by, gn), limit = cap)
    wf <- if (is.null(w)) "" else " * __dm_w"
    id <- list(n = sc$add(ck$S, "1"), sw = sc$add(ck$S, paste0("CAST(1 AS DOUBLE)", wf)),
               sy = sc$add(ck$S, paste0("__dm_y1", wf)),
               within = if (is.null(w)) sc$add(ck$S, "var_pop(__dm_y1) * COUNT(*)", agg = FALSE))
    sc$run()
    res <- sc$get(id$n)$res
    X <- dm_sql_cell_design(D, res, cap)
    if (is.null(X)) return(NULL)
    kc <- c(if (length(ck$S)) paste0("__dm_c", ck$S), ck$by, gn)
    within <- if (is.null(w)) {
        sc$get(id$within)$m
    } else {
        ## weighted: a second pass, around the means of the cells
        kq <- dm_qi(D, kc)
        on <- if (length(kc)) paste(sprintf("b.%s IS NOT DISTINCT FROM p.%s", kq, kq), collapse = " AND ") else "TRUE"
        q <- sprintf(paste("WITH b AS (%s),",
                           "p AS (SELECT %s FROM b%s)",
                           "SELECT %s FROM b JOIN p ON %s WHERE b.__dm_w <> 0%s"),
                     base, paste(c(kq, "SUM(__dm_w * __dm_y1) / SUM(__dm_w) AS __dm_m"), collapse = ", "),
                     if (length(kc)) " GROUP BY ALL" else "",
                     paste(c(paste0("b.", kq), "SUM(b.__dm_w * pow(b.__dm_y1 - p.__dm_m, 2)) AS __dm_within"),
                           collapse = ", "),
                     on, if (length(kc)) " GROUP BY ALL" else "")
        wr <- dm_sql_decode(D, dm_query(D, q))
        ## in the order of the cells
        v <- as.double(wr$`__dm_within`)[if (length(kc)) dm_match_rows(res[kc], wr[kc]) else 1L]
        v[is.na(v)] <- 0
        v
    }
    list(x = X, n = sc$get(id$n)$m, sw = sc$get(id$sw)$m, sy = sc$get(id$sy)$m, within = within,
         keys = if (length(gn)) as.matrix(res[gn]))
}

## For each row of `a`, the matching row of `b` (same values, NA matching NA).
dm_match_rows <- function(a, b) {
    key <- function(d) do.call(paste, c(lapply(d, function(x) {
        if (is.double(x)) sprintf("%a", x) else as.character(x)
    }), sep = "\r"))
    match(key(a), key(b))
}

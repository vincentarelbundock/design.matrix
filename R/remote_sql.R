## The SQL engine of remote designs (see remote.R).
##
## Term t of the design has X_t[i, ] = N_t(i) C_t[F_t(i), ], where F_t(i)
## is the combined level of the factors S_t of the term (the first factor
## varying fastest), N_t(i) the product of its numeric variables, and C_t a
## small matrix with one row per combined level. Then
##
##     (X'WX)[cols_s, cols_t] = C_s' M_st C_t,
##     M_st[a, b] = sum over rows with F_s = a, F_t = b of w N_s N_t,
##
## and M_st is a GROUP BY over the factor codes of S_s and S_t. X'Wy and
## grouped column sums are the same with one term. Predictions join each
## term's table of level -> C_t beta_t to the rows.
##
## C_t is read off the in-memory design of a few prototype rows (every
## combined level of S_t, numeric variables set to 1), so that columns,
## contrasts and their order are exactly those of model.matrix().

## Largest C_t (levels x columns) the SQL engine builds.
dm_sql_max_cells <- 2^23

## A numeric constant, exactly (the string is parsed with correct rounding).
dm_sql_double <- function(x) {
    x <- as.double(x)
    s <- ifelse(is.nan(x), "NaN", ifelse(x == Inf, "Infinity", ifelse(x == -Inf, "-Infinity",
                sprintf("%.17g", x))))
    sprintf("CAST('%s' AS DOUBLE)", s)
}

## Translate a numeric R expression into SQL over the columns of `s`.
## Returns list(sql = ) or list(reason = ).
dm_sql_num <- function(expr, types, con) {
    fail <- function(msg) {
        stop(structure(class = c("dm_untranslatable", "error", "condition"),
                       list(message = msg, call = NULL)))
    }
    unary <- c(exp = "exp", sqrt = "sqrt", abs = "abs", log2 = "log2", log10 = "log10")
    tr <- function(e) {
        if (is.name(e)) {
            nm <- as.character(e)
            ty <- types[nm]
            if (is.na(ty)) fail(sprintf("`%s` is not a column of the remote data", nm))
            if (!ty %in% c("double", "integer", "logical"))
                fail(sprintf("column `%s` is not numeric", nm))
            return(sprintf("CAST(s.%s AS DOUBLE)", DBI::dbQuoteIdentifier(con, nm)))
        }
        if (is.numeric(e) && length(e) == 1L && is.null(attributes(e))) {
            if (is.na(e)) fail("a missing constant")
            return(dm_sql_double(e))
        }
        if (is.call(e) && !is.name(e[[1L]])) fail(sprintf("%s() has no SQL translation", deparse2(e[[1L]])))
        if (!is.call(e)) fail("an unsupported expression")
        fn <- as.character(e[[1L]])
        args <- as.list(e)[-1L]
        an <- names(args)
        nargs <- length(args)
        if (fn == "log") {
            if (nargs < 1L || nargs > 2L) fail("log() with unsupported arguments")
            if (is.null(an)) an <- rep("", nargs)
            if (!all(an %in% c("", "x", "base"))) fail("log() with unsupported arguments")
            ## match the arguments as R would: names first, then positions
            ix <- match("x", an)
            ib <- match("base", an)
            pos <- which(an == "")
            if (is.na(ix)) { ix <- pos[1L]; pos <- pos[-1L] }
            if (is.na(ib) && nargs == 2L) ib <- pos[1L]
            x <- args[[ix]]
            if (nargs == 1L) return(sprintf("ln(%s)", tr(x)))
            b <- args[[ib]]
            if (!is.numeric(b) || length(b) != 1L) fail("log() with a base that is not a number")
            if (b == 10) return(sprintf("log10(%s)", tr(x)))
            if (b == 2) return(sprintf("log2(%s)", tr(x)))
            return(sprintf("(ln(%s) / ln(%s))", tr(x), dm_sql_double(b)))
        }
        if (!is.null(an) && any(an != "")) fail(sprintf("%s() with named arguments", fn))
        if (fn %in% c("(", "I") && nargs == 1L) return(sprintf("(%s)", tr(args[[1L]])))
        if (fn %in% names(unary) && nargs == 1L) return(sprintf("%s(%s)", unary[[fn]], tr(args[[1L]])))
        if (fn %in% c("+", "-") && nargs == 1L) return(sprintf("(%s%s)", fn, tr(args[[1L]])))
        if (fn %in% c("+", "-", "*", "/") && nargs == 2L)
            return(sprintf("(%s %s %s)", tr(args[[1L]]), fn, tr(args[[2L]])))
        if (fn == "^" && nargs == 2L) return(sprintf("pow(%s, %s)", tr(args[[1L]]), tr(args[[2L]])))
        fail(sprintf("%s() has no SQL translation", fn))
    }
    tryCatch(list(sql = tr(expr)), dm_untranslatable = function(e) list(reason = conditionMessage(e)))
}

## Term structures (S, M, L, strides, cols, C) for the SQL engine.
dm_sql_prepare <- function(D) {
    plan <- D$plan
    vi <- match(dm_terms_vars(plan$terms)$names, names(D$vars))
    kinds <- vapply(D$vars, function(v) v$kind, "")
    out <- vector("list", length(D$terms))
    for (k in seq_along(D$terms)) {
        t <- D$terms[[k]]
        S <- t$vars[kinds[t$vars] %in% c("factor", "logical")]
        M <- t$vars[kinds[t$vars] == "numeric"]
        L <- vapply(S, function(j) length(D$vars[[j]]$levels), 1L)
        Lt <- as.integer(prod(L))
        if (t$index == 0L) {
            C <- matrix(1, 1L, 1L)
        } else {
            grid <- expand.grid(lapply(L, seq_len), KEEP.OUT.ATTRS = FALSE)
            vals <- lapply(vi, function(j) {
                v <- D$vars[[j]]
                s <- match(j, S)
                codes <- if (is.na(s)) rep(1L, Lt) else grid[[s]]
                switch(v$kind,
                       factor = factor(v$levels[codes], levels = v$levels),
                       logical = c(FALSE, TRUE)[codes],
                       rep(1, Lt))
            })
            nd <- structure(list(), row.names = c(NA_integer_, -Lt), class = "data.frame")
            X <- as.matrix(dm_apply_general(plan, nd, stats::na.pass, vals = vals))
            C <- unname(X[, attr(X, "assign") == t$index, drop = FALSE])
        }
        out[[k]] <- list(S = S, M = M, L = L, Lt = Lt,
                         strides = cumprod(c(1, L))[seq_along(L)],
                         cols = which(D$assign == t$index), C = C)
    }
    D$sql_terms <- out
    D
}

## The query computing, for each row used, the factor codes (__dm_c<j>) and
## numeric variables (__dm_v<j>), plus `extra` select items.
dm_sql_base <- function(D, extra = character()) {
    con <- dm_con(D)
    items <- extra
    joins <- character()
    for (j in seq_along(D$vars)) {
        v <- D$vars[[j]]
        if (v$kind == "factor") {
            lv <- v$levels
            vals <- paste0("(", DBI::dbQuoteString(con, lv), ", ", seq_along(lv), ")", collapse = ", ")
            joins <- c(joins, sprintf("LEFT JOIN (VALUES %s) AS dm_l%d(value, code) ON CAST(s.%s AS VARCHAR) = dm_l%d.value",
                                      vals, j, dm_qi(D, v$col), j))
            items <- c(items, sprintf("dm_l%d.code AS __dm_c%d", j, j))
        } else if (v$kind == "logical") {
            items <- c(items, sprintf("CASE WHEN s.%s THEN 2 ELSE 1 END AS __dm_c%d", dm_qi(D, v$col), j))
        } else if (v$kind == "numeric") {
            items <- c(items, sprintf("%s AS __dm_v%d", v$sql, j))
        }
    }
    if (!length(items)) items <- "1 AS __dm_one"
    sprintf("SELECT %s FROM (%s) AS s %s WHERE %s", paste(items, collapse = ", "),
            D$src$select, paste(joins, collapse = " "), D$filter)
}

## Product of numeric variables (and weights), as SQL over the base query.
dm_sql_prod <- function(M, w = FALSE, y = NULL) {
    f <- c(if (length(M)) paste0("__dm_v", sort(M)), y, if (w) "__dm_w")
    if (!length(f)) "CAST(1 AS DOUBLE)" else paste(f, collapse = " * ")
}

## Combined level (0-based) of the factors S from their code columns.
dm_sql_level <- function(res, S, strides) {
    out <- numeric(nrow(res))
    for (k in seq_along(S)) out <- out + (res[[paste0("__dm_c", S[k])]] - 1) * strides[k]
    out
}

dm_sql_extra_w <- function(D, w) {
    if (is.null(w)) character() else sprintf("CAST(s.%s AS DOUBLE) AS __dm_w", dm_qi(D, w))
}

## Queries for X'WX: one per set of factors U = S_s + S_t.
dm_sql_gram_queries <- function(D, w) {
    st <- D$sql_terms
    base <- dm_sql_base(D, dm_sql_extra_w(D, w))
    groups <- list()
    for (t in seq_along(st)) for (s in seq_len(t)) {
        U <- sort(union(st[[s]]$S, st[[t]]$S))
        key <- paste(U, collapse = ",")
        prod <- dm_sql_prod(c(st[[s]]$M, st[[t]]$M), !is.null(w))
        g <- groups[[key]] %||% list(U = U, aggs = character(), pairs = list())
        a <- match(prod, g$aggs)
        if (is.na(a)) {
            g$aggs <- c(g$aggs, prod)
            a <- length(g$aggs)
        }
        g$pairs[[length(g$pairs) + 1L]] <- c(s, t, a)
        groups[[key]] <- g
    }
    lapply(groups, function(g) {
        keys <- if (length(g$U)) paste0("__dm_c", g$U) else character()
        aggs <- sprintf("SUM(%s) AS __dm_m%d", g$aggs, seq_along(g$aggs))
        g$sql <- sprintf("WITH b AS (%s) SELECT %s FROM b%s", base, paste(c(keys, aggs), collapse = ", "),
                         if (length(keys)) " GROUP BY ALL" else "")
        g
    })
}

## Number of rows used where a column is missing.
dm_sql_count_na <- function(D, cols) {
    if (!length(cols)) return(integer())
    items <- vapply(cols, function(col) {
        q <- dm_qi(D, col)
        nan <- if (D$src$types[[col]] == "double") sprintf(" OR isnan(s.%s)", q) else ""
        sprintf("SUM(CASE WHEN s.%s IS NULL%s THEN 1 ELSE 0 END)", q, nan)
    }, "")
    res <- dm_query(D, sprintf("SELECT %s FROM (%s) AS s WHERE %s",
                               paste(items, collapse = ", "), D$src$select, D$filter))
    out <- vapply(res, function(x) as.numeric(x %||% 0), 0)
    out[is.na(out)] <- 0
    names(out) <- cols
    out
}

dm_sql_gram <- function(D, w = NULL) {
    st <- D$sql_terms
    p <- D$p
    G <- matrix(0, p, p, dimnames = list(D$colnames, D$colnames))
    if (!is.null(w) && dm_sql_count_na(D, w) > 0) {
        G[] <- NA_real_
        return(G)
    }
    for (g in dm_sql_gram_queries(D, w)) {
        res <- dm_query(D, g$sql)
        for (pr in g$pairs) {
            s <- pr[1L]
            t <- pr[2L]
            m <- as.double(res[[paste0("__dm_m", pr[3L])]])
            m[is.na(m)] <- 0
            Ts <- st[[s]]
            Tt <- st[[t]]
            if (!length(Ts$cols) || !length(Tt$cols)) next
            a <- dm_sql_level(res, Ts$S, Ts$strides)
            if (s == t) {
                mm <- numeric(Ts$Lt)
                mm[a + 1] <- m
                G[Ts$cols, Ts$cols] <- crossprod(Ts$C, mm * Ts$C)
            } else {
                b <- dm_sql_level(res, Tt$S, Tt$strides)
                Mst <- Matrix::sparseMatrix(i = a + 1, j = b + 1, x = m, dims = c(Ts$Lt, Tt$Lt))
                blk <- as.matrix(crossprod(Ts$C, as.matrix(Mst %*% Tt$C)))
                G[Ts$cols, Tt$cols] <- blk
                G[Tt$cols, Ts$cols] <- t(blk)
            }
        }
    }
    G
}

## X'W Y for SQL expressions `ycols` over the base query (with `extra`
## select items defining them): p x K.
dm_sql_xty_core <- function(D, ycols, extra, w) {
    st <- D$sql_terms
    K <- length(ycols)
    out <- matrix(0, D$p, K)
    base <- dm_sql_base(D, c(extra, dm_sql_extra_w(D, w)))
    keys <- vapply(st, function(t) paste(t$S, collapse = ","), "")
    for (key in unique(keys)) {
        ts <- which(keys == key)
        S <- st[[ts[1L]]]$S
        aggs <- character()
        for (t in ts) for (k in seq_len(K))
            aggs <- c(aggs, sprintf("SUM(%s) AS __dm_m%d_%d", dm_sql_prod(st[[t]]$M, !is.null(w), ycols[k]), t, k))
        kc <- if (length(S)) paste0("__dm_c", S) else character()
        res <- dm_query(D, sprintf("WITH b AS (%s) SELECT %s FROM b%s", base, paste(c(kc, aggs), collapse = ", "),
                                   if (length(kc)) " GROUP BY ALL" else ""))
        a <- dm_sql_level(res, S, st[[ts[1L]]]$strides)
        for (t in ts) {
            Tt <- st[[t]]
            for (k in seq_len(K)) {
                m <- as.double(res[[sprintf("__dm_m%d_%d", t, k)]])
                m[is.na(m)] <- 0
                mm <- numeric(Tt$Lt)
                mm[a + 1] <- m
                out[Tt$cols, k] <- crossprod(Tt$C, mm)
            }
        }
    }
    out
}

## y columns as SQL over the base query: list(exprs, extra items).
dm_sql_ycols <- function(D, yy) {
    if (yy$response) {
        v <- D$vars[[yy$cols]]
        if (v$kind != "numeric") stop("The response must be numeric.", call. = FALSE)
        return(list(exprs = paste0("__dm_v", match(yy$cols, names(D$vars))), extra = character()))
    }
    for (col in yy$cols) {
        if (!D$src$types[[col]] %in% c("double", "integer", "logical"))
            stop(sprintf("Column `%s` is not numeric.", col), call. = FALSE)
    }
    K <- length(yy$cols)
    list(exprs = paste0("__dm_y", seq_len(K)),
         extra = sprintf("CAST(s.%s AS DOUBLE) AS __dm_y%d", dm_qi(D, yy$cols), seq_len(K)))
}

dm_sql_xty <- function(D, yy, w = NULL) {
    yc <- dm_sql_ycols(D, yy)
    out <- dm_sql_xty_core(D, yc$exprs, yc$extra, w)
    dimnames(out) <- list(D$colnames, yy$cols)
    na <- dm_sql_count_na(D, c(if (!yy$response) yy$cols, w))
    if (!is.null(w) && na[[w]] > 0) out[] <- NA_real_
    if (!yy$response) out[, na[yy$cols] > 0] <- NA_real_
    out
}

dm_sql_moments <- function(D, yy, w = NULL) {
    if (is.null(yy)) {
        xtx <- dm_sql_gram(D, w)
        base <- dm_sql_base(D, dm_sql_extra_w(D, w))
        res <- dm_query(D, sprintf("WITH b AS (%s) SELECT COUNT(*) AS __dm_n, SUM(%s) AS __dm_sw FROM b",
                                   base, if (is.null(w)) "CAST(1 AS DOUBLE)" else "__dm_w"))
        sumw <- as.double(res$`__dm_sw`)
        if (is.na(sumw) || (!is.null(w) && dm_sql_count_na(D, w) > 0)) sumw <- if (is.null(w)) 0 else NA_real_
        return(list(xtx = xtx, xty = NULL, yty = NULL, n = as.numeric(res$`__dm_n`), sumw = sumw))
    }
    yc <- dm_sql_ycols(D, yy)
    K <- length(yc$exprs)
    xtx <- dm_sql_gram(D, w)
    xty <- dm_sql_xty(D, yy, w)
    base <- dm_sql_base(D, c(yc$extra, dm_sql_extra_w(D, w)))
    wf <- if (is.null(w)) "" else " * __dm_w"
    idx <- expand.grid(k = seq_len(K), l = seq_len(K))
    aggs <- c("COUNT(*) AS __dm_n", sprintf("SUM(CAST(1 AS DOUBLE)%s) AS __dm_sw", wf),
              sprintf("SUM(%s * %s%s) AS __dm_yy%d", yc$exprs[idx$k], yc$exprs[idx$l], wf, seq_len(nrow(idx))))
    res <- dm_query(D, sprintf("WITH b AS (%s) SELECT %s FROM b", base, paste(aggs, collapse = ", ")))
    yty <- matrix(vapply(seq_len(nrow(idx)), function(r) as.double(res[[paste0("__dm_yy", r)]]), 0), K, K,
                  dimnames = list(yy$cols, yy$cols))
    yty[is.na(yty)] <- 0
    n <- as.numeric(res$`__dm_n`)
    sumw <- as.double(res$`__dm_sw`)
    if (is.na(sumw)) sumw <- 0
    na <- dm_sql_count_na(D, c(if (!yy$response) yy$cols, w))
    if (!is.null(w) && na[[w]] > 0) {
        yty[] <- NA_real_
        sumw <- NA_real_
    }
    if (!yy$response) {
        bad <- na[yy$cols] > 0
        yty[bad, ] <- NA_real_
        yty[, bad] <- NA_real_
    }
    list(xtx = xtx, xty = xty, yty = yty, n = n, sumw = sumw)
}

dm_sql_colsums <- function(D, group, w, mean) {
    if (is.null(group)) {
        out <- dm_sql_xty_core(D, "CAST(1 AS DOUBLE)", character(), w)[, 1L]
        if (!is.null(w) && dm_sql_count_na(D, w) > 0) out[] <- NA_real_
        names(out) <- D$colnames
        return(if (mean) out / D$n else out)
    }
    st <- D$sql_terms
    extra <- c(sprintf("s.%s AS __dm_g", dm_qi(D, group)), dm_sql_extra_w(D, w))
    base <- dm_sql_base(D, extra)
    cnt <- dm_query(D, sprintf("WITH b AS (%s) SELECT __dm_g, COUNT(*) AS __dm_n FROM b GROUP BY ALL", base))
    ug <- cnt$`__dm_g`
    if (anyNA(ug)) warning("missing values for 'group'", call. = FALSE)
    ug <- sort(ug, na.last = TRUE, method = "quick")
    G <- length(ug)
    out <- matrix(0, G, D$p, dimnames = list(as.character(ug), D$colnames))
    keys <- vapply(st, function(t) paste(t$S, collapse = ","), "")
    for (key in unique(keys)) {
        ts <- which(keys == key)
        S <- st[[ts[1L]]]$S
        aggs <- sprintf("SUM(%s) AS __dm_m%d", vapply(ts, function(t) dm_sql_prod(st[[t]]$M, !is.null(w)), ""), ts)
        kc <- c("__dm_g", if (length(S)) paste0("__dm_c", S))
        res <- dm_query(D, sprintf("WITH b AS (%s) SELECT %s FROM b GROUP BY ALL", base,
                                   paste(c(kc, aggs), collapse = ", ")))
        gi <- match(res$`__dm_g`, ug)
        a <- dm_sql_level(res, S, st[[ts[1L]]]$strides)
        for (t in ts) {
            Tt <- st[[t]]
            m <- as.double(res[[paste0("__dm_m", t)]])
            m[is.na(m)] <- 0
            A <- Matrix::sparseMatrix(i = gi, j = a + 1, x = m, dims = c(G, Tt$Lt))
            out[, Tt$cols] <- as.matrix(A %*% Tt$C)
        }
    }
    if (!is.null(w) && dm_sql_count_na(D, w) > 0) out[] <- NA_real_
    if (mean) out <- out / cnt$`__dm_n`[match(ug, cnt$`__dm_g`)]
    out
}

dm_sql_predict <- function(D, B, name, outcols, keep, overwrite) {
    con <- dm_con(D)
    st <- D$sql_terms
    K <- ncol(B)
    base <- dm_sql_base(D, paste0("s.", dm_qi(D, keep)))
    joins <- character()
    pieces <- rep(list(character()), K)
    for (k in seq_along(st)) {
        Tt <- st[[k]]
        if (!length(Tt$cols)) next
        V <- Tt$C %*% B[Tt$cols, , drop = FALSE]
        N <- if (length(Tt$M)) paste0(paste0("b.__dm_v", sort(Tt$M), collapse = " * "), " * ") else ""
        if (Tt$Lt == 1L) {
            for (j in seq_len(K)) pieces[[j]] <- c(pieces[[j]], paste0(N, dm_sql_double(V[1L, j])))
        } else {
            alias <- sprintf("dm_p%d", k)
            vals <- do.call(paste, c(list(seq_len(Tt$Lt) - 1L), lapply(seq_len(K), function(j) dm_sql_double(V[, j])),
                                     sep = ", "))
            keyexpr <- paste(sprintf("(b.__dm_c%d - 1) * %.0f", Tt$S, Tt$strides), collapse = " + ")
            joins <- c(joins, sprintf("LEFT JOIN (VALUES %s) AS %s(key, %s) ON %s = %s.key",
                                      paste0("(", vals, ")", collapse = ", "), alias,
                                      paste0("v", seq_len(K), collapse = ", "), keyexpr, alias))
            for (j in seq_len(K)) pieces[[j]] <- c(pieces[[j]], sprintf("%s%s.v%d", N, alias, j))
        }
    }
    preds <- vapply(pieces, function(x) if (length(x)) paste(x, collapse = " + ") else "CAST(0 AS DOUBLE)", "")
    select <- sprintf("WITH b AS (%s) SELECT %s FROM b %s", base,
                      paste(c(paste0("b.", dm_qi(D, keep)), sprintf("%s AS %s", preds, dm_qi(D, outcols))),
                            collapse = ", "),
                      paste(joins, collapse = " "))
    if (is.null(name)) {
        name <- dm_tmp_name("dm_pred")
        dm_execute(D, sprintf("CREATE TEMPORARY VIEW %s AS %s", dm_qi(D, name), select))
    } else {
        dm_execute(D, sprintf("CREATE %sTABLE %s AS %s", if (overwrite) "OR REPLACE " else "",
                              dm_qi(D, name), select))
    }
    name
}

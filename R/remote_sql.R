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

## Term structures (S, M, L, strides, cols, C) of a remote design. A term
## with a matrix variable (a variable computed in R with several columns,
## such as poly(x, 2)) is the sum of one part per column of the matrix (per
## combination of columns, with several such variables): each part has the
## columns of the term, its own C and one scalar numeric per variable, so
## the parts of a design may share columns, and their products add up.
## `term`: the index of the part's term in D$terms; `M`: the names of its
## numeric columns in the base query (dm_sql_base()).
dm_sql_prepare <- function(D) {
    plan <- D$plan
    vi <- match(dm_terms_vars(plan$terms)$names, names(D$vars))
    kinds <- vapply(D$vars, function(v) v$kind, "")
    nk <- vapply(D$vars, function(v) as.integer(v[["k"]] %||% 1L), 1L)
    out <- list()
    for (k in seq_along(D$terms)) {
        t <- D$terms[[k]]
        S <- t$vars[kinds[t$vars] %in% c("factor", "logical")]
        Mv <- t$vars[kinds[t$vars] == "numeric"]
        L <- vapply(S, function(j) length(D$vars[[j]]$levels), 1L)
        Lt <- as.integer(prod(L))
        parts <- if (length(Mv)) as.matrix(expand.grid(lapply(nk[Mv], seq_len), KEEP.OUT.ATTRS = FALSE))
                 else matrix(0L, 1L, 0L)
        for (r in seq_len(nrow(parts))) {
            u <- parts[r, ]
            if (t$index == 0L) {
                C <- matrix(1, 1L, 1L)
            } else {
                grid <- expand.grid(lapply(L, seq_len), KEEP.OUT.ATTRS = FALSE)
                vals <- lapply(vi, function(j) {
                    v <- D$vars[[j]]
                    s <- match(j, S)
                    codes <- if (is.na(s)) rep(1L, Lt) else grid[[s]]
                    if (v$kind == "factor") return(factor(v$levels[codes], levels = v$levels))
                    if (v$kind == "logical") return(c(FALSE, TRUE)[codes])
                    if (nk[j] == 1L) return(rep(1, Lt))
                    ## a matrix variable: the column of this part
                    m <- matrix(1, Lt, nk[j])
                    if (j %in% Mv) {
                        m[] <- 0
                        m[, u[match(j, Mv)]] <- 1
                    }
                    m
                })
                nd <- structure(list(), row.names = c(NA_integer_, -Lt), class = "data.frame")
                X <- as.matrix(dm_apply_general(plan, nd, stats::na.pass, vals = vals))
                C <- unname(X[, attr(X, "assign") == t$index, drop = FALSE])
            }
            out[[length(out) + 1L]] <- list(term = k, S = S, L = L, Lt = Lt,
                                            M = if (length(Mv)) dm_sql_vnames(Mv, nk[Mv], u) else character(),
                                            strides = cumprod(c(1, L))[seq_along(L)],
                                            cols = which(D$assign == t$index), C = C)
        }
    }
    D$sql_terms <- out
    D
}

## Names, in the base query, of the numeric columns of variables j (with
## k columns each): column u of a matrix variable is __dm_v<j>_<u>.
dm_sql_vnames <- function(j, k, u) {
    n <- max(length(j), length(u))
    j <- rep_len(j, n)
    k <- rep_len(k, n)
    ifelse(k == 1L, paste0("__dm_v", j), paste0("__dm_v", j, "_", rep_len(u, n)))
}

## The query computing, for each row used, the factor codes (__dm_c<j>) and
## numeric variables (__dm_v<j>), plus `extra`: named SQL expressions over
## the source `s` (the responses, weights and groups of a product), each
## computed as a column of that name (doubles, except those in `asis`).
dm_sql_base <- function(D, extra = character(), asis = character()) {
    con <- dm_con(D)
    items <- if (length(extra)) sprintf(ifelse(names(extra) %in% asis, "(%s) AS %s", "CAST((%s) AS DOUBLE) AS %s"),
                                        extra, names(extra)) else character()
    joins <- character()
    for (j in seq_along(D$vars)) {
        v <- D$vars[[j]]
        if (v$kind == "factor") {
            lv <- v$levels
            vals <- paste0("(", DBI::dbQuoteString(con, lv), ", ", seq_along(lv), ")", collapse = ", ")
            ## (names that no expression of the user over `s` can capture)
            joins <- c(joins, sprintf("LEFT JOIN (VALUES %s) AS dm_l%d(__dm_lv, __dm_lc) ON CAST(s.%s AS VARCHAR) = dm_l%d.__dm_lv",
                                      vals, j, dm_qi(D, v$col), j))
            items <- c(items, sprintf("dm_l%d.__dm_lc AS __dm_c%d", j, j))
        } else if (v$kind == "logical") {
            q <- dm_qi(D, v$col)
            items <- c(items, sprintf("CASE WHEN s.%s IS NULL THEN NULL WHEN s.%s THEN 2 ELSE 1 END AS __dm_c%d",
                                      q, q, j))
        } else if (v$kind == "numeric") {
            items <- c(items, sprintf("%s AS %s", v$sql, dm_sql_vnames(j, length(v$sql), seq_along(v$sql))))
        }
    }
    if (!length(items)) items <- "1 AS __dm_one"
    sprintf("SELECT %s FROM (%s) AS s %s WHERE %s", paste(items, collapse = ", "),
            D$src$select, paste(joins, collapse = " "), D$filter)
}

## Product of numeric columns (and weights), as SQL over the base query.
dm_sql_prod <- function(M, w = FALSE, y = NULL) {
    f <- c(sort(M), y, if (w) "__dm_w")
    if (!length(f)) "CAST(1 AS DOUBLE)" else paste(f, collapse = " * ")
}

## Combined level (0-based) of the factors S from their code columns.
dm_sql_level <- function(res, S, strides) {
    out <- numeric(nrow(res))
    for (k in seq_along(S)) out <- out + (res[[paste0("__dm_c", S[k])]] - 1) * strides[k]
    out
}

## One scan of the rows used for many sums at once. Each item is a SQL
## expression over the base query `b`, summed within the combined levels of
## a set U of factors: add(U, e) returns its id, run() executes, and
## get(id) gives list(res, m): the cells (their code columns __dm_c<j>, at
## least those of U) and the sums. All items share one GROUP BY over every
## factor involved, so that the data are read once; a cell finer than U is
## summed again by the caller (dm_sql_fold(), or a sparse matrix, which adds
## duplicates). When that joint grouping would return too many cells, items
## are grouped by their own U instead, one query per set. `by`: columns of
## the base query that split every sum (the groups of colSums()); their
## number is unknown, so the joint grouping is then used only for few cells.
dm_sql_scan <- function(D, base, levels, by = character()) {
    cap <- if (length(by)) 256 else dm_sql_max_cells
    U <- list()
    e <- character()
    out <- NULL
    plan <- function() {
        all <- sort(unique(unlist(U)))
        ue <- unique(e)
        joint <- prod(as.numeric(levels[all])) * length(ue) <= cap
        key <- if (joint) rep("", length(e)) else vapply(U, paste, "", collapse = ",")
        lapply(unique(key), function(k) {
            i <- which(key == k)
            keys <- if (joint) all else U[[i[1L]]]
            ge <- unique(e[i])
            kc <- c(by, if (length(keys)) paste0("__dm_c", keys))
            list(items = i, col = match(e[i], ge),
                 sql = sprintf("WITH b AS (%s) SELECT %s FROM b%s", base,
                               paste(c(kc, sprintf("SUM(%s) AS __dm_m%d", ge, seq_along(ge))), collapse = ", "),
                               if (length(kc)) " GROUP BY ALL" else ""))
        })
    }
    list(
        add = function(u, expr) {
            U[[length(U) + 1L]] <<- sort(as.integer(u))
            e <<- c(e, expr)
            length(e)
        },
        sql = function() vapply(plan(), function(g) g$sql, ""),
        run = function() {
            out <<- vector("list", length(e))
            for (g in plan()) {
                res <- dm_query(D, g$sql)
                for (k in seq_along(g$items)) {
                    m <- as.double(res[[paste0("__dm_m", g$col[k])]])
                    m[is.na(m)] <- 0
                    out[[g$items[k]]] <<- list(res = res, m = m)
                }
            }
            invisible()
        },
        get = function(id) out[[id]],
        total = function(id) sum(out[[id]]$m))
}

## 1 where a (double) expression of the base query is missing, to sum.
dm_sql_missing <- function(e) sprintf("CASE WHEN %s IS NULL OR isnan(%s) THEN 1 ELSE 0 END", e, e)

## Sums of x within the 0-based levels `a` (of L levels).
dm_sql_fold <- function(a, x, L) {
    out <- numeric(L)
    ok <- !is.na(a)
    if (!any(ok)) return(out)
    r <- rowsum(x[ok], a[ok], reorder = FALSE)
    out[as.numeric(rownames(r)) + 1] <- r[, 1L]
    out
}

## With na.pass, rows with a missing variable are kept: the columns of a
## term with a missing value in some row are NA in every product, as in the
## dense computation. SQL condition "no variable of term k is missing".
dm_sql_term_ok <- function(D, k) {
    vars <- D$terms[[k]]$vars
    if (!length(vars)) return("TRUE")
    paste(sprintf("(%s)", vapply(D$vars[vars], function(v) dm_remote_var_ok(D, v), "")), collapse = " AND ")
}

## Which terms have a missing value in a row used? (Cached.)
dm_sql_na_terms <- function(D) {
    nt <- length(D$terms)
    if (!identical(D$na, "pass") || !nt) return(rep(FALSE, nt))
    if (!is.null(D$cache$na_terms)) return(D$cache$na_terms)
    items <- sprintf("SUM(CASE WHEN COALESCE(%s, FALSE) THEN 0 ELSE 1 END) AS __dm_bad%d",
                     vapply(seq_len(nt), function(k) dm_sql_term_ok(D, k), ""), seq_len(nt))
    res <- dm_query(D, sprintf("SELECT %s FROM (%s) AS s", paste(items, collapse = ", "), D$src$select))
    out <- vapply(seq_len(nt), function(k) isTRUE(as.numeric(res[[k]]) > 0), NA)
    D$cache$na_terms <- out
    out
}

## The same for each part of D$sql_terms.
dm_sql_na_parts <- function(D) dm_sql_na_terms(D)[vapply(D$sql_terms, function(t) t$term, 1L)]

dm_sql_na_cols <- function(D) {
    unique(unlist(lapply(D$sql_terms[dm_sql_na_parts(D)], function(t) t$cols)))
}

## X'WX, X'WY, Y'WY, the number of rows and the sum of the weights, in one
## scan (dm_sql_scan()). Columns are SQL expressions over the source `s`
## (see dm_remote_cols()): `ys` the columns Y (a sum with a missing value is
## NA), `w` the weights, `wy` other weights for the sums with Y (X'[wy]Y;
## Y'Y is then not computed). With sql_only, the queries instead.
dm_sql_moments <- function(D, ys = character(), w = NULL, gram = TRUE, sql_only = FALSE, wy = NULL) {
    st <- D$sql_terms
    p <- D$p
    K <- length(ys)
    yc <- list(exprs = if (K) paste0("__dm_y", seq_len(K)) else character())
    hasw <- !is.null(w)
    wf <- if (hasw) " * __dm_w" else ""
    bad <- dm_sql_na_parts(D)
    live <- which(!bad & vapply(st, function(t) length(t$cols) > 0L, NA))
    sc <- dm_sql_scan(D, dm_sql_base(D, c(stats::setNames(ys, yc$exprs), `__dm_w` = w, `__dm_wy` = wy)),
                      vapply(D$vars, function(v) max(1L, length(v$levels)), 1L))
    yprod <- function(M, y) if (is.null(wy)) dm_sql_prod(M, hasw, y) else dm_sql_prod(M, FALSE, c(y, "__dm_wy"))
    id <- list(n = sc$add(integer(), "1"), sw = sc$add(integer(), paste0("CAST(1 AS DOUBLE)", wf)),
               naw = if (hasw) sc$add(integer(), dm_sql_missing("__dm_w")),
               nay = vapply(yc$exprs, function(e) sc$add(integer(), dm_sql_missing(e)), 1L))
    gid <- list()
    if (gram) for (t in live) for (s in live[live <= t])
        gid[[length(gid) + 1L]] <- c(s, t, sc$add(c(st[[s]]$S, st[[t]]$S), dm_sql_prod(c(st[[s]]$M, st[[t]]$M), hasw)))
    hid <- matrix(0L, length(st), K)
    for (t in live) for (k in seq_len(K))
        hid[t, k] <- sc$add(st[[t]]$S, yprod(st[[t]]$M, yc$exprs[k]))
    yid <- matrix(0L, K, K)
    if (is.null(wy)) for (k in seq_len(K)) for (l in seq_len(k))
        yid[k, l] <- yid[l, k] <- sc$add(integer(), sprintf("%s * %s%s", yc$exprs[l], yc$exprs[k], wf))
    if (sql_only) return(sc$sql())
    sc$run()

    xtx <- xty <- yty <- NULL
    if (gram) {
        xtx <- matrix(0, p, p, dimnames = list(D$colnames, D$colnames))
        for (pr in gid) {
            Ts <- st[[pr[1L]]]
            Tt <- st[[pr[2L]]]
            g <- sc$get(pr[3L])
            a <- dm_sql_level(g$res, Ts$S, Ts$strides)
            ## (added: the parts of a term share its columns)
            if (pr[1L] == pr[2L]) {
                xtx[Ts$cols, Ts$cols] <- xtx[Ts$cols, Ts$cols] + crossprod(Ts$C, dm_sql_fold(a, g$m, Ts$Lt) * Ts$C)
            } else {
                b <- dm_sql_level(g$res, Tt$S, Tt$strides)
                Mst <- Matrix::sparseMatrix(i = a + 1, j = b + 1, x = g$m, dims = c(Ts$Lt, Tt$Lt))
                blk <- as.matrix(crossprod(Ts$C, as.matrix(Mst %*% Tt$C)))
                xtx[Ts$cols, Tt$cols] <- xtx[Ts$cols, Tt$cols] + blk
                xtx[Tt$cols, Ts$cols] <- xtx[Tt$cols, Ts$cols] + t(blk)
            }
        }
    }
    na <- dm_sql_na_cols(D)
    if (gram) xtx[na, ] <- xtx[, na] <- NA_real_
    sumw <- sc$total(id$sw)
    if (K) {
        xty <- matrix(0, p, K, dimnames = list(D$colnames, NULL))
        for (t in live) for (k in seq_len(K)) {
            Tt <- st[[t]]
            g <- sc$get(hid[t, k])
            xty[Tt$cols, k] <- xty[Tt$cols, k] +
                crossprod(Tt$C, dm_sql_fold(dm_sql_level(g$res, Tt$S, Tt$strides), g$m, Tt$Lt))
        }
        xty[na, ] <- NA_real_
        miss <- vapply(id$nay, sc$total, 0) > 0
        xty[, miss] <- NA_real_
        if (is.null(wy)) {
            yty <- matrix(vapply(yid, sc$total, 0), K, K)
            yty[miss, ] <- yty[, miss] <- NA_real_
        }
    }
    ## a missing weight: every sum is missing, as in the dense computation
    if (hasw && sc$total(id$naw) > 0) {
        if (gram) xtx[] <- NA_real_
        if (K) xty[] <- NA_real_
        if (!is.null(yty)) yty[] <- NA_real_
        sumw <- NA_real_
    }
    list(xtx = xtx, xty = xty, yty = yty, n = sc$total(id$n), sumw = sumw)
}

## Column sums (or means) of the design, within the values of `group` if
## given, with the weights `w` (SQL expressions over the source `s`): one
## scan (dm_sql_scan()).
dm_sql_colsums <- function(D, group, w, mean) {
    st <- D$sql_terms
    hasw <- !is.null(w)
    grouped <- !is.null(group)
    bad <- dm_sql_na_terms(D)
    pbad <- dm_sql_na_parts(D)
    extra <- c(`__dm_g` = group, `__dm_w` = w,
               if (grouped && any(bad)) stats::setNames(
                   sprintf("CASE WHEN COALESCE(%s, FALSE) THEN 0 ELSE 1 END",
                           vapply(which(bad), function(k) dm_sql_term_ok(D, k), "")), paste0("__dm_bad", which(bad))))
    sc <- dm_sql_scan(D, dm_sql_base(D, extra, asis = "__dm_g"),
                      vapply(D$vars, function(v) max(1L, length(v$levels)), 1L),
                      by = if (grouped) "__dm_g" else character())
    idn <- sc$add(integer(), "1")
    idw <- if (hasw) sc$add(integer(), dm_sql_missing("__dm_w"))
    idb <- if (grouped) vapply(which(bad), function(k) sc$add(integer(), sprintf("__dm_bad%d", k)), 1L)
    ## (within groups, a term with missing values is summed where it has none)
    use <- which(vapply(st, function(t) length(t$cols) > 0L, NA) & (grouped | !pbad))
    idt <- integer(length(st))
    for (t in use) idt[t] <- sc$add(st[[t]]$S, dm_sql_prod(st[[t]]$M, hasw))
    sc$run()
    naw <- hasw && sc$total(idw) > 0
    if (!grouped) {
        out <- stats::setNames(numeric(D$p), D$colnames)
        for (t in use) {
            Tt <- st[[t]]
            g <- sc$get(idt[t])
            out[Tt$cols] <- out[Tt$cols] + crossprod(Tt$C, dm_sql_fold(dm_sql_level(g$res, Tt$S, Tt$strides), g$m, Tt$Lt))
        }
        out[dm_sql_na_cols(D)] <- NA_real_
        if (naw) out[] <- NA_real_
        return(if (mean) out / D$n else out)
    }
    ug <- unique(sc$get(idn)$res$`__dm_g`)
    if (anyNA(ug)) warning("missing values for 'group'", call. = FALSE)
    ug <- sort(ug, na.last = TRUE, method = "quick")
    G <- length(ug)
    ## sums of an item within the groups
    by_group <- function(id) {
        g <- sc$get(id)
        dm_sql_fold(match(g$res$`__dm_g`, ug) - 1, g$m, G)
    }
    out <- matrix(0, G, D$p, dimnames = list(as.character(ug), D$colnames))
    for (t in use) {
        Tt <- st[[t]]
        g <- sc$get(idt[t])
        a <- dm_sql_level(g$res, Tt$S, Tt$strides)
        ## (rows with a missing factor: their groups are set to NA below)
        ok <- !is.na(a)
        A <- Matrix::sparseMatrix(i = match(g$res$`__dm_g`, ug)[ok], j = a[ok] + 1, x = g$m[ok], dims = c(G, Tt$Lt))
        out[, Tt$cols] <- out[, Tt$cols] + as.matrix(A %*% Tt$C)
    }
    for (j in seq_along(idb)) out[by_group(idb[j]) > 0, which(D$assign == D$terms[[which(bad)[j]]]$index)] <- NA_real_
    if (naw) out[] <- NA_real_
    if (mean) out <- out / by_group(idn)
    out
}

dm_sql_predict <- function(D, B, name, outcols, keep, overwrite) {
    con <- dm_con(D)
    st <- D$sql_terms
    K <- ncol(B)
    qk <- dm_qi(D, keep)
    base <- dm_sql_base(D, stats::setNames(paste0("s.", qk), qk), asis = qk)
    joins <- character()
    pieces <- rep(list(character()), K)
    for (k in seq_along(st)) {
        Tt <- st[[k]]
        if (!length(Tt$cols)) next
        V <- Tt$C %*% B[Tt$cols, , drop = FALSE]
        N <- if (length(Tt$M)) paste0(paste0("b.", sort(Tt$M), collapse = " * "), " * ") else ""
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

## Translate a logical R expression (for `subset`) into SQL over `s`:
## comparisons, &, |, !, %in% with constants, is.na(), logical columns and
## constants; numeric operands as in dm_sql_num(). R's three-valued logic
## with NA is SQL's with NULL, and rows where the condition is NULL are
## dropped, as subset() drops rows where it is NA. Returns list(sql = ) or
## list(reason = ).
dm_sql_bool <- function(expr, types, con) {
    fail <- function(msg) {
        stop(structure(class = c("dm_untranslatable", "error", "condition"),
                       list(message = msg, call = NULL)))
    }
    col <- function(nm) sprintf("s.%s", DBI::dbQuoteIdentifier(con, nm))
    ## an operand of a comparison: a string, a text column, or a number
    operand <- function(e) {
        if (is.character(e) && length(e) == 1L && !is.na(e)) return(as.character(DBI::dbQuoteString(con, e)))
        if (is.name(e) && isTRUE(types[as.character(e)] %in% c("character", "factor")))
            return(sprintf("CAST(%s AS VARCHAR)", col(as.character(e))))
        tr <- dm_sql_num(e, types, con)
        if (is.null(tr$sql)) fail(tr$reason)
        tr$sql
    }
    tr <- function(e) {
        if (is.logical(e) && length(e) == 1L && !is.na(e)) return(if (e) "TRUE" else "FALSE")
        if (is.name(e)) {
            nm <- as.character(e)
            if (isTRUE(types[nm] == "logical")) return(col(nm))
            fail(sprintf("`%s` is not a logical column of the remote data", nm))
        }
        if (!is.call(e) || !is.name(e[[1L]])) fail("an unsupported expression")
        fn <- as.character(e[[1L]])
        args <- as.list(e)[-1L]
        if (fn == "(" && length(args) == 1L) return(sprintf("(%s)", tr(args[[1L]])))
        if (fn == "!" && length(args) == 1L) return(sprintf("(NOT %s)", tr(args[[1L]])))
        if (fn %in% c("&", "&&") && length(args) == 2L) return(sprintf("(%s AND %s)", tr(args[[1L]]), tr(args[[2L]])))
        if (fn %in% c("|", "||") && length(args) == 2L) return(sprintf("(%s OR %s)", tr(args[[1L]]), tr(args[[2L]])))
        ops <- c("==" = "=", "!=" = "<>", "<" = "<", ">" = ">", "<=" = "<=", ">=" = ">=")
        if (fn %in% names(ops) && length(args) == 2L)
            return(sprintf("(%s %s %s)", operand(args[[1L]]), ops[[fn]], operand(args[[2L]])))
        if (fn == "%in%" && length(args) == 2L) {
            v <- args[[2L]]
            if (is.call(v) && identical(v[[1L]], as.name("c"))) v <- unlist(as.list(v)[-1L])
            if (!(is.numeric(v) || is.character(v)) || anyNA(v) || !length(v))
                fail("%in% needs a vector of constants")
            vals <- if (is.character(v)) DBI::dbQuoteString(con, v) else dm_sql_double(v)
            return(sprintf("COALESCE(%s IN (%s), FALSE)", operand(args[[1L]]), paste(vals, collapse = ", ")))
        }
        if (fn == "is.na" && length(args) == 1L) {
            a <- args[[1L]]
            if (is.name(a) && isTRUE(types[as.character(a)] %in% c("character", "factor", "logical", "integer")))
                return(sprintf("(%s IS NULL)", col(as.character(a))))
            x <- operand(a)
            return(sprintf("(%s IS NULL OR isnan(%s))", x, x))
        }
        fail(sprintf("%s() has no SQL translation", fn))
    }
    tryCatch(list(sql = tr(expr)), dm_untranslatable = function(e) list(reason = conditionMessage(e)))
}

## A `weights`, `offset` or `subset` argument for remote data, as SQL over
## the columns of the source `s`: NULL; the name of a column (a string, or
## a bare name that is a column); a DBI::SQL() expression; or an R
## expression with a SQL translation (dm_sql_num(), or dm_sql_bool() for
## `subset`). A bare name that is not a column is evaluated in `env`, and
## its value resolved in turn (a string, DBI::SQL()). Returns NULL or
## list(sql, inputs).
dm_resolve_arg <- function(x, env, types, con, what, bool = FALSE) {
    if (is.null(x)) return(NULL)
    wrap <- function(sql) if (bool) sprintf("(%s)", sql) else sprintf("CAST((%s) AS DOUBLE)", sql)
    if (inherits(x, "SQL")) return(list(sql = wrap(as.character(x)), inputs = character()))
    if (is.character(x)) {
        if (length(x) != 1L || !x %in% names(types))
            stop(sprintf("`%s`: no column %s in the remote data.", what, paste0("`", x, "`", collapse = ", ")),
                 call. = FALSE)
        if (bool && types[[x]] != "logical")
            stop(sprintf("`%s`: column `%s` is not logical.", what, x), call. = FALSE)
        if (!bool && !types[[x]] %in% c("double", "integer", "logical"))
            stop(sprintf("`%s`: column `%s` is not numeric.", what, x), call. = FALSE)
        q <- sprintf("s.%s", DBI::dbQuoteIdentifier(con, x))
        return(list(sql = if (bool) q else sprintf("CAST(%s AS DOUBLE)", q), inputs = x))
    }
    if (is.name(x) && !as.character(x) %in% names(types)) {
        val <- tryCatch(eval(x, env), error = function(e) NULL)
        if (inherits(val, "SQL") || is.character(val)) return(dm_resolve_arg(val, env, types, con, what, bool))
        stop(sprintf("`%s` must name a column of the remote data, or be a DBI::SQL() or R expression over its columns.",
                     what), call. = FALSE)
    }
    if (is.call(x) && deparse2(x[[1L]]) %in% c("SQL", "DBI::SQL"))
        return(dm_resolve_arg(eval(x, env), env, types, con, what, bool))
    if (!is.language(x))
        stop(sprintf("`%s` must name a column of the remote data, or be a DBI::SQL() or R expression over its columns (values in R cannot be used).",
                     what), call. = FALSE)
    tr <- if (bool) dm_sql_bool(x, types, con) else dm_sql_num(x, types, con)
    if (is.null(tr$sql))
        stop(sprintf("`%s` has no SQL translation: %s.", what, tr$reason), call. = FALSE)
    list(sql = tr$sql, inputs = intersect(all.vars(x), names(types)))
}

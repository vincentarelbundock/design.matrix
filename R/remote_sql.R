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
## and M_st is a GROUP BY over the factors of S_s and S_t (their columns as
## stored; levels become codes in R, on the cells). X'Wy and grouped column
## sums are the same with one term. Predictions read each term's C_t beta_t
## at the row's level (xtx_sql_predict()).
##
## C_t comes from the set-up of the prototype design of a few rows
## (xtx_sql_layout()): the contrasts and columns that set-up chose, which
## are exactly those of model.matrix().

## Largest C_t (levels x columns) the SQL engine builds.
xtx_sql_max_cells <- 2^23

## A numeric constant, exactly (the string is parsed with correct rounding).
xtx_sql_double <- function(x) {
    x <- as.double(x)
    s <- ifelse(is.nan(x), "NaN", ifelse(x == Inf, "Infinity", ifelse(x == -Inf, "-Infinity",
                sprintf("%.17g", x))))
    sprintf("CAST('%s' AS DOUBLE)", s)
}

## Translate a numeric R expression into SQL over the columns of `s`.
## Returns list(sql = ) or list(reason = ).
xtx_sql_num <- function(expr, types, con) {
    fail <- function(msg) {
        stop(structure(class = c("xtx_untranslatable", "error", "condition"),
                       list(message = msg, call = NULL)))
    }
    unary <- c(exp = "exp", sqrt = "sqrt", abs = "abs")
    ## a constant base of log(): a number, its negation, Inf or NaN (where no
    ## column has that name, as model.frame() finds them), NA
    const <- function(e) {
        if (is.call(e) && identical(e[[1L]], as.name("-")) && length(e) == 2L) {
            v <- const(e[[2L]])
            return(if (is.null(v)) NULL else -v)
        }
        if (is.name(e) && as.character(e) %in% c("Inf", "NaN") && is.na(types[as.character(e)]))
            return(get(as.character(e), baseenv()))
        if ((is.numeric(e) || is.logical(e)) && length(e) == 1L && is.null(attributes(e))) return(as.double(e))
        NULL
    }
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
            return(xtx_sql_double(e))
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
            ## as R, also at 0, below, and for any base (xtx_sql_log())
            if (nargs == 1L) return(xtx_sql_rlog(tr(x)))
            b <- const(args[[ib]])
            if (is.null(b)) fail("log() with a base that is not a number")
            return(xtx_sql_log(tr(x), b))
        }
        if (fn %in% c("log2", "log10", "log1p") && nargs == 1L && is.null(an)) {
            v <- tr(args[[1L]])
            return(switch(fn, log2 = xtx_sql_log(v, 2), log10 = xtx_sql_log(v, 10), log1p = xtx_sql_log1p(v)))
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
    tryCatch(list(sql = tr(expr)), xtx_untranslatable = function(e) list(reason = conditionMessage(e)))
}

## Term structures (S, M, L, strides, cols, C) of a remote design. A term
## with a matrix variable (a variable computed in R with several columns,
## such as poly(x, 2)) is the sum of one part per column of the matrix (per
## combination of columns, with several such variables): each part has the
## columns of the term, its own C and one scalar numeric per variable, so
## the parts of a design may share columns, and their products add up.
## `term`: the index of the part's term in D$terms; `M`: the names of its
## numeric columns in the base query (xtx_sql_base()).
##
## The C of each part comes from the set-up of the prepared prototype
## design (`setup`: xtx_sql_layout()), the one source of the columns of the
## design. Only when there is no such set-up is it read off a design built
## on a grid of levels (xtx_sql_layout_grid()).
xtx_sql_prepare <- function(D, setup = NULL) {
    vi <- match(xtx_terms_vars(D$plan$terms)$names, names(D$vars))
    kinds <- vapply(D$vars, function(v) v$kind, "")
    nk <- vapply(D$vars, function(v) as.integer(v[["k"]] %||% 1L), 1L)
    out <- list()
    for (k in seq_along(D$terms)) {
        t <- D$terms[[k]]
        S <- t$vars[kinds[t$vars] %in% c("factor", "logical")]
        Mv <- t$vars[kinds[t$vars] == "numeric"]
        L <- vapply(S, function(j) length(D$vars[[j]]$levels), 1L)
        parts <- if (length(Mv)) as.matrix(expand.grid(lapply(nk[Mv], seq_len), KEEP.OUT.ATTRS = FALSE))
                 else matrix(0L, 1L, 0L)
        for (r in seq_len(nrow(parts))) {
            u <- parts[r, ]
            out[[length(out) + 1L]] <- list(term = k, S = S, L = L, Lt = as.integer(prod(L)),
                                            M = if (length(Mv)) xtx_sql_vnames(Mv, nk[Mv], u) else character(),
                                            strides = cumprod(c(1, L))[seq_along(L)],
                                            cols = which(D$assign == t$index), Mv = Mv, u = u)
        }
    }
    C <- xtx_sql_layout(D, setup, out, vi, nk)
    ## (which of the two built the layout: "setup" or "grid")
    D$sql_layout <- if (is.null(C)) "grid" else "setup"
    if (is.null(C)) C <- xtx_sql_layout_grid(D, out, vi, nk)
    D$sql_terms <- Map(function(part, Ct) {
        part$C <- Ct
        part$Mv <- part$u <- NULL
        part
    }, out, C)
    D
}

## The matrices C of the parts (xtx_sql_prepare()) from the set-up of the
## prototype design (the list of C_xtx_setup: see src/modelmatrix.c), which
## fixed its columns: `factors` (variables x terms: 0, not in the term; 1,
## coded by contrasts, `contr1`; 2, by dummy variables, `contr2`), `nlevs`
## (0 for a numeric variable), `columns` (those of a numeric variable) and
## `count` (the columns of each term; 0 for a term with none). Variable i of
## the set-up is variable i of the terms of the plan, D$vars[[vi[i]]].
##
## Within a term, the column of a combination of one column per variable is
## sum_i c_i * stride_i (c_i 0-based; stride_i the product of the columns
## of the variables of the term before i, factors and numerics alike), so
## the columns of a part (the columns `u` of its numeric variables) need
## not be consecutive. The rows of C are the combined levels, the first
## factor varying fastest: C is the Kronecker product of the factors'
## contrast matrices, written into the columns of the part.
##
## Nothing here reads the values of the prototype's variables: the layout
## does not depend on whether the in-memory kernels would accept them.
## NULL when there is no set-up (a prototype built by stats::model.matrix()
## because this version of R is not recognised: xtx_prepare()), or when it
## does not describe the terms as the design does (not expected).
xtx_sql_layout <- function(D, setup, parts, vi, nk) {
    if (is.null(setup)) return(NULL)
    nVar <- setup$nVar
    nterms <- setup$nterms
    fac <- matrix(as.integer(setup$factors), nVar, nterms)
    out <- vector("list", length(parts))
    for (r in seq_along(parts)) {
        part <- parts[[r]]
        t <- D$terms[[part$term]]
        if (t$index == 0L) {
            if (!identical(setup$intercept, 1L)) return(NULL)
            out[[r]] <- matrix(1, 1L, 1L)
            next
        }
        if (t$index > nterms) return(NULL)
        count <- setup$count[t$index]
        if (length(part$cols) != count) return(NULL)
        if (count == 0L) {
            out[[r]] <- matrix(0, part$Lt, 0L)
            next
        }
        f <- fac[, t$index]
        ii <- which(f != 0L)
        ## the variables of the term, in the order of the set-up
        if (length(ii) > length(vi) || !identical(sort(vi[ii]), sort(as.integer(t$vars)))) return(NULL)
        isf <- setup$nlevs[ii] > 0L
        cm <- lapply(seq_along(ii), function(a) {
            if (!isf[a]) return(NULL)
            m <- (if (f[ii[a]] == 1L) setup$contr1 else setup$contr2)[[ii[a]]]
            if (!is.matrix(m)) m <- as.matrix(m)
            if (!is.double(m)) storage.mode(m) <- "double"
            if (is.null(dimnames(m))) m else unname(m)
        })
        nc <- vapply(seq_along(ii), function(a) if (isf[a]) ncol(cm[[a]]) else setup$columns[ii[a]], 1L)
        if (!identical(as.integer(vi[ii[isf]]), as.integer(part$S)) ||
            !identical(as.integer(vi[ii[!isf]]), as.integer(part$Mv)) ||
            !identical(unname(vapply(cm[isf], nrow, 1L)), unname(as.integer(part$L))) ||
            !identical(unname(as.integer(nc[!isf])), unname(as.integer(nk[part$Mv]))) ||
            prod(as.numeric(nc)) != count) return(NULL)
        stride <- cumprod(c(1, nc))[seq_along(nc)]
        ## the Kronecker product of the contrasts, and the 0-based column
        ## of each of its columns among those of the term
        K <- NULL
        col <- sum((as.numeric(part$u) - 1) * stride[!isf])
        for (a in which(isf)) {
            K <- if (is.null(K)) cm[[a]] else kronecker(cm[[a]], K)
            col <- as.vector(outer(col, (seq_len(nc[a]) - 1) * stride[a], `+`))
        }
        if (is.null(K)) K <- matrix(1, 1L, 1L)
        if (identical(as.numeric(col), as.numeric(seq_len(count) - 1L))) {
            Ct <- K
        } else {
            Ct <- matrix(0, part$Lt, count)
            Ct[, col + 1] <- K
        }
        out[[r]] <- Ct
    }
    out
}

## FALLBACK, used only when the prototype design has no set-up
## (xtx_sql_layout() returned NULL): the matrices C of the parts read off
## the in-memory design of prototype rows (every combined level of the
## factors of the term, numeric variables set to 1), built by applying the
## plan once per part. Also the reference of xtx_sql_layout() in the tests.
xtx_sql_layout_grid <- function(D, parts, vi, nk) {
    plan <- D$plan
    lapply(parts, function(part) {
        t <- D$terms[[part$term]]
        if (t$index == 0L) return(matrix(1, 1L, 1L))
        S <- part$S
        Mv <- part$Mv
        Lt <- part$Lt
        grid <- expand.grid(lapply(part$L, seq_len), KEEP.OUT.ATTRS = FALSE)
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
                m[, part$u[match(j, Mv)]] <- 1
            }
            m
        })
        nd <- structure(list(), row.names = c(NA_integer_, -Lt), class = "data.frame")
        X <- as.matrix(xtx_apply_general(plan, nd, stats::na.pass, vals = vals))
        unname(X[, attr(X, "assign") == t$index, drop = FALSE])
    })
}

## Names, in the base query, of the numeric columns of variables j (with
## k columns each): column u of a matrix variable is __xtx_v<j>_<u>.
xtx_sql_vnames <- function(j, k, u) {
    n <- max(length(j), length(u))
    j <- rep_len(j, n)
    k <- rep_len(k, n)
    ifelse(k == 1L, paste0("__xtx_v", j), paste0("__xtx_v", j, "_", rep_len(u, n)))
}

## The query computing, for each row used, the factor keys (__xtx_c<j>) and
## numeric variables (__xtx_v<j>), plus `extra`: named SQL expressions over
## the source `s` (the responses, weights and groups of a product), each
## computed as a column of that name (doubles, except those in `asis`).
## The key of a factor or logical variable is its column as stored (no
## join to a table of codes: GROUP BY over a raw column is cheapest), and
## xtx_sql_decode() turns the keys of a small aggregated result into codes;
## xtx_sql_code() is the code as SQL, where a row needs it.
xtx_sql_base <- function(D, extra = character(), asis = character()) {
    items <- if (length(extra)) sprintf(ifelse(names(extra) %in% asis, "(%s) AS %s", "CAST((%s) AS DOUBLE) AS %s"),
                                        extra, names(extra)) else character()
    for (j in seq_along(D$vars)) {
        v <- D$vars[[j]]
        if (v$kind %in% c("factor", "logical")) {
            items <- c(items, sprintf("%s AS __xtx_c%d", xtx_sql_key(D, v), j))
        } else if (v$kind == "numeric") {
            items <- c(items, sprintf("%s AS %s", v$sql, xtx_sql_vnames(j, length(v$sql), seq_along(v$sql))))
        }
    }
    if (!length(items)) items <- "1 AS __xtx_one"
    sprintf("SELECT %s FROM (%s) AS s WHERE %s", paste(items, collapse = ", "), D$src$select, D$filter)
}

## The key of a factor or logical variable, as SQL over the source `s`: the
## column itself, when R reads its values as the strings that match the
## levels (text, ENUM) or as logicals; else its text, which is what the
## levels match.
xtx_sql_key <- function(D, v) {
    q <- sprintf("s.%s", xtx_qi(D, v$col))
    ty <- D$src$types[[v$col]]
    if (v$kind == "logical") {
        if (identical(ty, "logical")) return(q)
        return(sprintf("CASE WHEN %s IS NULL THEN NULL WHEN %s THEN 2 ELSE 1 END", q, q))
    }
    if (isTRUE(ty %in% c("character", "factor"))) q else sprintf("CAST(%s AS VARCHAR)", q)
}

## The keys __xtx_c<j> of a query result (xtx_sql_base()) as codes: the
## position of the level (NA: no level), 1 (FALSE) or 2 (TRUE).
xtx_sql_decode <- function(D, res) {
    for (nm in grep("^__xtx_c[0-9]+$", names(res), value = TRUE)) {
        v <- D$vars[[as.integer(substring(nm, 8L))]]
        x <- res[[nm]]
        res[[nm]] <- if (v$kind == "factor") match(as.character(x), v$levels)
                     else if (is.logical(x)) as.integer(x) + 1L
                     else as.integer(x)
    }
    res
}

## Largest table of codes or values written inline in SQL (else joined).
xtx_sql_max_inline <- 65536

## The code of factor or logical variable j as SQL over its key `ref`
## (xtx_sql_key()) in constant time: the logical as 1 or 2, an ENUM by its
## position in the type. NULL where no such expression exists (text): a
## join on the text then gives the values (a CASE over the levels, or a
## search of a list, costs more per row than the join beyond a few levels).
xtx_sql_code <- function(D, j, ref) {
    v <- D$vars[[j]]
    ty <- D$src$types[[v$col]]
    lv <- v$levels
    list_of <- function(pos) sprintf("[%s]", paste(ifelse(is.na(pos), "NULL", pos), collapse = ", "))
    if (v$kind == "logical") return(if (identical(ty, "logical")) sprintf("(CAST(%s AS INTEGER) + 1)", ref) else ref)
    if (identical(ty, "factor")) {
        dict <- levels(D$src$proto[[v$col]])
        if (is.null(dict) || length(dict) > xtx_sql_max_inline) return(NULL)
        if (identical(dict, lv)) return(sprintf("(enum_code(%s) + 1)", ref))
        return(sprintf("list_extract(%s, enum_code(%s) + 1)", list_of(match(dict, lv)), ref))
    }
    NULL
}

## Product of numeric columns (and weights), as SQL over the base query.
xtx_sql_prod <- function(M, w = FALSE, y = NULL) {
    f <- c(sort(M), y, if (w) "__xtx_w")
    if (!length(f)) "CAST(1 AS DOUBLE)" else paste(f, collapse = " * ")
}

## Combined level (0-based) of the factors S from their code columns.
xtx_sql_level <- function(res, S, strides) {
    out <- numeric(nrow(res))
    for (k in seq_along(S)) out <- out + (res[[paste0("__xtx_c", S[k])]] - 1) * strides[k]
    out
}

## One scan of the rows used for many sums at once. Each item is a SQL
## expression over the base query `b`, summed within the combined levels of
## a set U of factors: add(U, e) returns its id, run() executes, and
## get(id) gives list(res, m): the cells (their code columns __xtx_c<j>, at
## least those of U) and the sums.
##
## The items are grouped by the factors of their own terms, so that what
## comes back is bounded by the cells of the terms involved and not by
## those of all the factors of the design together: the scan groups by the
## largest sets U (a set within another is summed within the larger one,
## and summed again by the caller: xtx_sql_fold(), or a sparse matrix, which
## adds duplicates), all in one statement and one read of the data, with
## GROUP BY GROUPING SETS; the rows of each set are told apart by
## GROUPING(). One set, of all the factors, when that has no more cells
## than those sets together (and always with a single largest set: a plain
## GROUP BY). When even that would return too many values, each item is
## grouped by its own U, one query per set.
##
## `by`: columns of the base query that split every sum (the groups of
## colSums(), or the numeric variables of compressed products:
## xtx_sql_cells()); their number is unknown, so one statement is then
## used only for few cells. `limit`: for the trial of compressed products
## (one set of cells): when the scan has more cells than that, or as many
## as rows (xtx_sql_limit()), it returns none; they are counted in the
## database, and nothing of a trial that fails comes back.
## add(u, expr, agg = TRUE): `expr` is summed, or with agg = FALSE it is
## an aggregate expression itself.
xtx_sql_scan <- function(D, base, levels, by = character(), limit = NULL) {
    cap <- if (length(by)) 256 else xtx_sql_max_cells
    U <- list()
    e <- character()
    out <- NULL
    cells <- function(u) prod(as.numeric(levels[u]))
    within <- function(a, b) all(a %in% b)
    kcol <- function(u) if (length(u)) paste0("__xtx_c", u) else character()
    ## one statement: the items `i` summed within the sets of factors `sets`
    statement <- function(i, sets) {
        ge <- unique(e[i])
        keys <- sort(unique(unlist(sets)))
        kc <- c(by, kcol(keys))
        mask <- NULL
        if (length(sets) > 1L) {
            ## GROUPING(): one bit per key, set where the key is not in the
            ## grouping set, the first key in the highest bit
            mask <- vapply(sets, function(u) sum(2^(length(keys) - which(!keys %in% u))), 0)
            sel <- c(kc, sprintf("CAST(GROUPING(%s) AS INTEGER) AS __xtx_gs", paste(kcol(keys), collapse = ", ")))
            group <- sprintf(" GROUP BY GROUPING SETS (%s)",
                             paste(vapply(sets, function(u) sprintf("(%s)", paste(c(by, kcol(u)), collapse = ", ")), ""),
                                   collapse = ", "))
        } else {
            sel <- kc
            group <- if (length(kc)) " GROUP BY ALL" else ""
        }
        sql <- sprintf("SELECT %s FROM b%s", paste(c(sel, sprintf("%s AS __xtx_m%d", ge, seq_along(ge))), collapse = ", "),
                       group)
        if (length(kc) && !is.null(limit)) {
            lim <- xtx_sql_limit(D, limit)
            sql <- sprintf("WITH b AS (%s), c AS MATERIALIZED (%s LIMIT %.0f) SELECT * FROM c WHERE (SELECT COUNT(*) FROM c) <= %.0f",
                           base, sql, lim + 1, lim)
        } else {
            sql <- sprintf("WITH b AS (%s) %s", base, sql)
        }
        ## the set of each item: the smallest that has its factors
        size <- vapply(sets, cells, 0)
        owner <- vapply(U[i], function(u) {
            ok <- which(vapply(sets, function(v) within(u, v), NA))
            ok[which.min(size[ok])]
        }, 1L)
        list(items = i, col = match(e[i], ge), sql = sql, owner = owner, mask = mask)
    }
    plan <- function() {
        if (!length(e)) return(list())
        sets <- unique(U)
        top <- sets[!vapply(sets, function(a) any(vapply(sets, function(b) length(b) > length(a) && within(a, b), NA)), NA)]
        all <- sort(unique(unlist(U)))
        if (length(top) > 1L && cells(all) <= sum(vapply(top, cells, 0))) top <- list(all)
        ## (a trial with `limit` has one set; GROUPING() is read as an integer)
        if (sum(vapply(top, cells, 0)) * length(unique(e)) <= cap && length(all) <= 30L &&
            (is.null(limit) || length(top) == 1L))
            return(list(statement(seq_along(e), top)))
        key <- vapply(U, paste, "", collapse = ",")
        lapply(unique(key), function(k) {
            i <- which(key == k)
            statement(i, U[i[1L]])
        })
    }
    list(
        add = function(u, expr, agg = TRUE) {
            U[[length(U) + 1L]] <<- sort(unique(as.integer(u)))
            e <<- c(e, if (agg) sprintf("SUM(%s)", expr) else expr)
            length(e)
        },
        sql = function() vapply(plan(), function(g) g$sql, ""),
        run = function() {
            out <<- vector("list", length(e))
            for (g in plan()) {
                res <- xtx_sql_decode(D, xtx_query(D, g$sql))
                ## the rows of each set
                parts <- if (is.null(g$mask)) list(res)
                         else lapply(g$mask, function(m) res[res$`__xtx_gs` == m, , drop = FALSE])
                for (k in seq_along(g$items)) {
                    r <- parts[[g$owner[k]]]
                    m <- as.double(r[[paste0("__xtx_m", g$col[k])]])
                    m[is.na(m)] <- 0
                    out[[g$items[k]]] <<- list(res = r, m = m)
                }
            }
            invisible()
        },
        get = function(id) out[[id]],
        total = function(id) sum(out[[id]]$m))
}

## The largest number of cells a trial of compressed products may have
## (`cap`: xtx_compress_cap()): at most the cap, and fewer than the rows
## (when every row is its own cell there is nothing to compress).
xtx_sql_limit <- function(D, cap) {
    n <- D$cache$n %||% D$n
    if (n > 1) min(cap, n - 1) else cap
}

## 1 where a (double) expression of the base query is missing, to sum.
xtx_sql_missing <- function(e) sprintf("CASE WHEN %s IS NULL OR isnan(%s) THEN 1 ELSE 0 END", e, e)

## Sums of x within the 0-based levels `a` (of L levels).
xtx_sql_fold <- function(a, x, L) {
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
xtx_sql_term_ok <- function(D, k) {
    vars <- D$terms[[k]]$vars
    if (!length(vars)) return("TRUE")
    paste(sprintf("(%s)", vapply(D$vars[vars], function(v) xtx_remote_var_ok(D, v), "")), collapse = " AND ")
}

## Which terms have a missing value in a row used? (Cached.)
xtx_sql_na_terms <- function(D) {
    nt <- length(D$terms)
    if (!identical(D$na, "pass") || !nt) return(rep(FALSE, nt))
    if (!is.null(D$cache$na_terms)) return(D$cache$na_terms)
    items <- sprintf("SUM(CASE WHEN COALESCE(%s, FALSE) THEN 0 ELSE 1 END) AS __xtx_bad%d",
                     vapply(seq_len(nt), function(k) xtx_sql_term_ok(D, k), ""), seq_len(nt))
    res <- xtx_query(D, sprintf("SELECT %s FROM (%s) AS s", paste(items, collapse = ", "), D$src$select))
    out <- vapply(seq_len(nt), function(k) isTRUE(as.numeric(res[[k]]) > 0), NA)
    D$cache$na_terms <- out
    out
}

## The same for each part of D$sql_terms.
xtx_sql_na_parts <- function(D) xtx_sql_na_terms(D)[vapply(D$sql_terms, function(t) t$term, 1L)]

xtx_sql_na_cols <- function(D) {
    unique(unlist(lapply(D$sql_terms[xtx_sql_na_parts(D)], function(t) t$cols)))
}

## X'WX, X'WY, Y'WY, the number of rows and the sum of the weights, in one
## scan (xtx_sql_scan()). Columns are SQL expressions over the source `s`
## (see xtx_remote_cols()): `ys` the columns Y (a sum with a missing value is
## NA), `w` the weights, `wy` other weights for the sums with Y (X'[wy]Y;
## Y'Y is then not computed). With sql_only, the queries instead. With
## `cap`, compressed (xtx_sql_cells(), not with `wy`): NULL when there are
## more cells than that.
xtx_sql_moments <- function(D, ys = character(), w = NULL, gram = TRUE, sql_only = FALSE, wy = NULL, cap = NULL, sum_y = FALSE) {
    st <- D$sql_terms
    p <- D$p
    K <- length(ys)
    yc <- list(exprs = if (K) paste0("__xtx_y", seq_len(K)) else character())
    hasw <- !is.null(w)
    wf <- if (hasw) " * __xtx_w" else ""
    bad <- xtx_sql_na_parts(D)
    live <- which(!bad & vapply(st, function(t) length(t$cols) > 0L, NA))
    ck <- xtx_sql_cells(D, cap)
    sc <- xtx_sql_scan(D, xtx_sql_base(D, c(stats::setNames(ys, yc$exprs), `__xtx_w` = w, `__xtx_wy` = wy)),
                      vapply(D$vars, function(v) max(1L, length(v$levels)), 1L), by = ck$by, limit = cap)
    yprod <- function(M, y) if (is.null(wy)) xtx_sql_prod(M, hasw, y) else xtx_sql_prod(M, FALSE, c(y, "__xtx_wy"))
    ## (compressed, every sum is within the cells)
    u0 <- ck$S
    id <- list(n = sc$add(u0, "1"), sw = sc$add(u0, paste0("CAST(1 AS DOUBLE)", wf)),
               naw = if (hasw) sc$add(u0, xtx_sql_missing("__xtx_w")),
               nay = vapply(yc$exprs, function(e) sc$add(u0, xtx_sql_missing(e)), 1L))
    if (sum_y) sy <- sc$add(u0, yprod(character(), yc$exprs[[1L]]))
    gid <- list()
    hid <- matrix(0L, length(st), K)
    if (!is.null(cap)) {
        ## the cells' sums of the weighted responses
        cy <- vapply(yc$exprs, function(e) sc$add(u0, yprod(character(), e)), 1L)
    } else {
        if (gram) for (t in live) for (s in live[live <= t])
            gid[[length(gid) + 1L]] <- c(s, t, sc$add(c(st[[s]]$S, st[[t]]$S), xtx_sql_prod(c(st[[s]]$M, st[[t]]$M), hasw)))
        for (t in live) for (k in seq_len(K))
            hid[t, k] <- sc$add(st[[t]]$S, yprod(st[[t]]$M, yc$exprs[k]))
    }
    yid <- matrix(0L, K, K)
    if (is.null(wy)) for (k in seq_len(K)) for (l in seq_len(k))
        yid[k, l] <- yid[l, k] <- sc$add(u0, sprintf("%s * %s%s", yc$exprs[l], yc$exprs[k], wf))
    if (sql_only) return(sc$sql())
    sc$run()

    xtx <- xty <- yty <- NULL
    if (!is.null(cap)) {
        X <- xtx_sql_cell_design(D, sc$get(id$n)$res, cap)
        if (is.null(X)) return(NULL)
        if (gram) xtx <- crossprod(X, sc$get(id$sw)$m * X)
        if (K) xty <- crossprod(X, matrix(vapply(cy, function(i) sc$get(i)$m, numeric(nrow(X))), nrow(X), K))
        D$cache$n <- sc$total(id$n)
    }
    if (gram && is.null(cap)) {
        xtx <- matrix(0, p, p, dimnames = list(D$colnames, D$colnames))
        for (pr in gid) {
            Ts <- st[[pr[1L]]]
            Tt <- st[[pr[2L]]]
            g <- sc$get(pr[3L])
            a <- xtx_sql_level(g$res, Ts$S, Ts$strides)
            ## (added: the parts of a term share its columns)
            if (pr[1L] == pr[2L]) {
                xtx[Ts$cols, Ts$cols] <- xtx[Ts$cols, Ts$cols] + crossprod(Ts$C, xtx_sql_fold(a, g$m, Ts$Lt) * Ts$C)
            } else {
                b <- xtx_sql_level(g$res, Tt$S, Tt$strides)
                Mst <- Matrix::sparseMatrix(i = a + 1, j = b + 1, x = g$m, dims = c(Ts$Lt, Tt$Lt))
                blk <- as.matrix(crossprod(Ts$C, as.matrix(Mst %*% Tt$C)))
                xtx[Ts$cols, Tt$cols] <- xtx[Ts$cols, Tt$cols] + blk
                xtx[Tt$cols, Ts$cols] <- xtx[Tt$cols, Ts$cols] + t(blk)
            }
        }
    }
    na <- xtx_sql_na_cols(D)
    if (gram) {
        dimnames(xtx) <- list(D$colnames, D$colnames)
        xtx[na, ] <- xtx[, na] <- NA_real_
    }
    sumw <- sc$total(id$sw)
    if (K) {
        if (is.null(cap)) {
            xty <- matrix(0, p, K)
            for (t in live) for (k in seq_len(K)) {
                Tt <- st[[t]]
                g <- sc$get(hid[t, k])
                xty[Tt$cols, k] <- xty[Tt$cols, k] +
                    crossprod(Tt$C, xtx_sql_fold(xtx_sql_level(g$res, Tt$S, Tt$strides), g$m, Tt$Lt))
            }
        }
        dimnames(xty) <- list(D$colnames, NULL)
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
    out <- list(xtx = xtx, xty = xty, yty = yty, n = sc$total(id$n), sumw = sumw)
    if (sum_y) out$sumy <- sc$total(sy)
    out
}

## Column sums (or means) of the design, within the values of `group` if
## given, with the weights `w` (SQL expressions over the source `s`): one
## scan (xtx_sql_scan()). With `cap`, compressed (xtx_sql_cells()): NULL
## when there are more cells than that.
xtx_sql_colsums <- function(D, group, w, mean, cap = NULL, group_totals = NULL) {
    st <- D$sql_terms
    hasw <- !is.null(w)
    grouped <- !is.null(group)
    bad <- xtx_sql_na_terms(D)
    pbad <- xtx_sql_na_parts(D)
    tn <- if (length(group_totals)) paste0("__xtx_total", seq_along(group_totals))
    extra <- c(`__xtx_g` = group, `__xtx_w` = w, stats::setNames(group_totals, tn),
               if (grouped && any(bad)) stats::setNames(
                   sprintf("CASE WHEN COALESCE(%s, FALSE) THEN 0 ELSE 1 END",
                           vapply(which(bad), function(k) xtx_sql_term_ok(D, k), "")), paste0("__xtx_bad", which(bad))))
    ck <- xtx_sql_cells(D, cap)
    sc <- xtx_sql_scan(D, xtx_sql_base(D, extra, asis = "__xtx_g"),
                      vapply(D$vars, function(v) max(1L, length(v$levels)), 1L),
                      by = c(if (grouped) "__xtx_g", ck$by), limit = cap)
    u0 <- ck$S
    total_ids <- vapply(tn, function(nm) sc$add(u0, nm), 1L)
    null_ids <- vapply(tn, function(nm) sc$add(u0, sprintf("CASE WHEN %s IS NULL THEN 1 ELSE 0 END", nm)), 1L)
    idn <- sc$add(u0, "1")
    idw <- if (hasw) sc$add(u0, xtx_sql_missing("__xtx_w"))
    idb <- if (grouped) vapply(which(bad), function(k) sc$add(u0, sprintf("__xtx_bad%d", k)), 1L)
    ## (within groups, a term with missing values is summed where it has none)
    use <- which(vapply(st, function(t) length(t$cols) > 0L, NA) & (grouped | !pbad))
    idt <- integer(length(st))
    if (is.null(cap)) {
        for (t in use) idt[t] <- sc$add(st[[t]]$S, xtx_sql_prod(st[[t]]$M, hasw))
    } else {
        ## (compressed: the cells' sums of the weights)
        ids <- sc$add(u0, xtx_sql_prod(character(), hasw))
    }
    sc$run()
    if (!is.null(cap)) {
        res <- sc$get(idn)$res
        X <- xtx_sql_cell_design(D, res, cap)
        if (is.null(X)) return(NULL)
        X <- sc$get(ids)$m * X
    }
    naw <- hasw && sc$total(idw) > 0
    if (!grouped) {
        out <- stats::setNames(numeric(D$p), D$colnames)
        if (!is.null(cap)) out[] <- colSums(X)
        else for (t in use) {
            Tt <- st[[t]]
            g <- sc$get(idt[t])
            out[Tt$cols] <- out[Tt$cols] + crossprod(Tt$C, xtx_sql_fold(xtx_sql_level(g$res, Tt$S, Tt$strides), g$m, Tt$Lt))
        }
        out[xtx_sql_na_cols(D)] <- NA_real_
        if (naw) out[] <- NA_real_
        return(if (mean) out / D$n else out)
    }
    ug <- unique(sc$get(idn)$res$`__xtx_g`)
    if (anyNA(ug)) warning("missing values for 'group'", call. = FALSE)
    ug <- sort(ug, na.last = TRUE, method = "quick")
    G <- length(ug)
    ## sums of an item within the groups
    by_group <- function(id) {
        g <- sc$get(id)
        xtx_sql_fold(match(g$res$`__xtx_g`, ug) - 1, g$m, G)
    }
    out <- matrix(0, G, D$p, dimnames = list(as.character(ug), D$colnames))
    if (!is.null(cap)) {
        A <- Matrix::sparseMatrix(i = match(res$`__xtx_g`, ug), j = seq_len(nrow(X)), x = 1, dims = c(G, nrow(X)))
        out[] <- as.matrix(A %*% X)
    }
    terms <- use[idt[use] > 0L]
    ## With only scalar numeric terms, the scan groups solely by __xtx_g.
    ## Every sum comes from the same result rows, one per output group.
    direct <- length(terms) && all(vapply(st[terms], function(t)
        !length(t$S) && length(t$cols) == 1L && identical(as.vector(t$C), 1), NA))
    if (direct) {
        rows <- match(ug, sc$get(idn)$res$`__xtx_g`)
        for (t in terms) out[, st[[t]]$cols] <- sc$get(idt[t])$m[rows]
    } else for (t in terms) {
        Tt <- st[[t]]
        g <- sc$get(idt[t])
        a <- xtx_sql_level(g$res, Tt$S, Tt$strides)
        ## (rows with a missing factor: their groups are set to NA below)
        ok <- !is.na(a)
        A <- Matrix::sparseMatrix(i = match(g$res$`__xtx_g`, ug)[ok], j = a[ok] + 1, x = g$m[ok], dims = c(G, Tt$Lt))
        out[, Tt$cols] <- out[, Tt$cols] + as.matrix(A %*% Tt$C)
    }
    for (j in seq_along(idb)) out[by_group(idb[j]) > 0, which(D$assign == D$terms[[which(bad)[j]]]$index)] <- NA_real_
    if (naw) out[] <- NA_real_
    if (mean) out <- out / by_group(idn)
    if (!length(tn)) return(out)
    totals <- matrix(0, G, length(tn))
    for (j in seq_along(tn)) {
        if (direct) {
            values <- sc$get(total_ids[j])$m[rows]
            missing <- sc$get(null_ids[j])$m[rows] > 0
        } else {
            values <- by_group(total_ids[j])
            missing <- by_group(null_ids[j]) > 0
        }
        values[missing] <- NA_real_
        totals[, j] <- values
    }
    list(keys = ug, sums = out, totals = totals)
}

## X B for the rows of the source `s` of the design, as SQL: list(preds, one
## expression per column of B; joins, the joins they read, or ""). Each term
## adds N_t times the row of C_t B_t of the row's combined level: read off
## a list of values by its index (codes by xtx_sql_code()), or, when a
## factor of the term has no such code (text), from a table joined on the
## keys, one join per set of factors shared by the terms that have them.
## The keys and numeric variables are those of the base query
## (xtx_sql_base()), written in place. The caller selects the expressions
## from `(source) AS s` followed by the joins (xtx_p_mult(): the node of a
## product).
xtx_sql_predict <- function(D, B) {
    con <- xtx_con(D)
    st <- D$sql_terms
    K <- ncol(B)
    joins <- character()
    pieces <- rep(list(character()), K)
    ## the numeric columns of the base query, by name, and the factor keys
    num <- unlist(lapply(seq_along(D$vars), function(j) {
        v <- D$vars[[j]]
        if (v$kind != "numeric") return(NULL)
        stats::setNames(sprintf("(%s)", v$sql), xtx_sql_vnames(j, length(v$sql), seq_along(v$sql)))
    }))
    keys <- list()
    key <- function(j) {
        if (is.null(keys[[as.character(j)]])) keys[[as.character(j)]] <<- sprintf("(%s)", xtx_sql_key(D, D$vars[[j]]))
        keys[[as.character(j)]]
    }
    codes <- list()
    code <- function(j) {
        if (is.null(codes[[as.character(j)]]))
            codes[[as.character(j)]] <<- xtx_sql_code(D, j, key(j)) %||% NA_character_
        codes[[as.character(j)]]
    }
    ## the joined tables: per set of factors, its key columns and values
    tabs <- list()
    for (k in seq_along(st)) {
        Tt <- st[[k]]
        if (!length(Tt$cols)) next
        V <- Tt$C %*% B[Tt$cols, , drop = FALSE]
        N <- if (length(Tt$M)) paste0(paste0(num[sort(Tt$M)], collapse = " * "), " * ") else ""
        if (Tt$Lt == 1L) {
            for (j in seq_len(K)) pieces[[j]] <- c(pieces[[j]], paste0(N, xtx_sql_double(V[1L, j])))
            next
        }
        cs <- vapply(Tt$S, code, "")
        if (!anyNA(cs) && Tt$Lt <= xtx_sql_max_inline) {
            ## (0-based combined level, as in xtx_sql_level())
            lev <- paste(sprintf("(%s - 1) * %.0f", cs, Tt$strides), collapse = " + ")
            for (j in seq_len(K))
                pieces[[j]] <- c(pieces[[j]], sprintf("%slist_extract([%s], %s + 1)", N,
                                                      paste(xtx_sql_double(V[, j]), collapse = ", "), lev))
            next
        }
        id <- paste(Tt$S, collapse = ",")
        if (is.null(tabs[[id]])) tabs[[id]] <- list(S = Tt$S, L = Tt$L, alias = sprintf("xtx_p%d", length(tabs) + 1L),
                                                    vals = list())
        tb <- tabs[[id]]
        m <- length(tb$vals)
        tabs[[id]]$vals <- c(tb$vals, lapply(seq_len(K), function(j) V[, j]))
        for (j in seq_len(K)) pieces[[j]] <- c(pieces[[j]], sprintf("%s%s.v%d", N, tb$alias, m + j))
    }
    for (tb in tabs) {
        ## the key of each factor: its code where SQL computes it, else its
        ## text against the level (first factor varying fastest)
        grid <- expand.grid(lapply(tb$L, seq_len), KEEP.OUT.ATTRS = FALSE)
        kcols <- paste0("k", seq_along(tb$S))
        kv <- on <- character()
        cols <- lapply(seq_along(tb$S), function(i) {
            j <- tb$S[i]
            cj <- code(j)
            if (!is.na(cj)) {
                on <<- c(on, sprintf("%s = %s.%s", cj, tb$alias, kcols[i]))
                as.character(grid[[i]])
            } else {
                on <<- c(on, sprintf("CAST(%s AS VARCHAR) = %s.%s", key(j), tb$alias, kcols[i]))
                as.character(DBI::dbQuoteString(con, D$vars[[j]]$levels[grid[[i]]]))
            }
        })
        vals <- do.call(paste, c(cols, lapply(tb$vals, xtx_sql_double), sep = ", "))
        joins <- c(joins, sprintf("LEFT JOIN (VALUES %s) AS %s(%s) ON %s",
                                  paste0("(", vals, ")", collapse = ", "), tb$alias,
                                  paste(c(kcols, paste0("v", seq_along(tb$vals))), collapse = ", "),
                                  paste(on, collapse = " AND ")))
    }
    preds <- vapply(pieces, function(x) if (length(x)) paste(x, collapse = " + ") else "CAST(0 AS DOUBLE)", "")
    list(preds = preds, joins = if (length(joins)) paste0(" ", paste(joins, collapse = " ")) else "")
}

## Translate a logical R expression (for `subset`) into SQL over `s`:
## comparisons, &, |, !, %in% with constants, is.na(), logical columns and
## constants; numeric operands as in xtx_sql_num(). R's three-valued logic
## with NA is SQL's with NULL, and rows where the condition is NULL are
## dropped, as subset() drops rows where it is NA. Returns list(sql = ) or
## list(reason = ).
xtx_sql_bool <- function(expr, types, con) {
    fail <- function(msg) {
        stop(structure(class = c("xtx_untranslatable", "error", "condition"),
                       list(message = msg, call = NULL)))
    }
    col <- function(nm) sprintf("s.%s", DBI::dbQuoteIdentifier(con, nm))
    ## an operand of a comparison: a string, a text column, or a number
    operand <- function(e) {
        if (is.character(e) && length(e) == 1L && !is.na(e)) return(as.character(DBI::dbQuoteString(con, e)))
        if (is.name(e) && isTRUE(types[as.character(e)] %in% c("character", "factor")))
            return(sprintf("CAST(%s AS VARCHAR)", col(as.character(e))))
        tr <- xtx_sql_num(e, types, con)
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
            vals <- if (is.character(v)) DBI::dbQuoteString(con, v) else xtx_sql_double(v)
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
    tryCatch(list(sql = tr(expr)), xtx_untranslatable = function(e) list(reason = conditionMessage(e)))
}

## A `weights`, `offset` or `subset` argument for remote data, as SQL over
## the columns of the source `s`: NULL; the name of a column (a string, or
## a bare name that is a column); a DBI::SQL() expression; or an R
## expression with a SQL translation (xtx_sql_num(), or xtx_sql_bool() for
## `subset`). A bare name that is not a column is evaluated in `env`, and
## its value resolved in turn (a string, DBI::SQL()). Returns NULL or
## list(sql, inputs).
xtx_resolve_arg <- function(x, env, types, con, what, bool = FALSE) {
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
        if (inherits(val, "SQL") || is.character(val)) return(xtx_resolve_arg(val, env, types, con, what, bool))
        stop(sprintf("`%s` must name a column of the remote data, or be a DBI::SQL() or R expression over its columns.",
                     what), call. = FALSE)
    }
    if (is.call(x) && deparse2(x[[1L]]) %in% c("SQL", "DBI::SQL"))
        return(xtx_resolve_arg(eval(x, env), env, types, con, what, bool))
    if (!is.language(x))
        stop(sprintf("`%s` must name a column of the remote data, or be a DBI::SQL() or R expression over its columns (values in R cannot be used).",
                     what), call. = FALSE)
    tr <- if (bool) xtx_sql_bool(x, types, con) else xtx_sql_num(x, types, con)
    if (is.null(tr$sql))
        stop(sprintf("`%s` has no SQL translation: %s.", what, tr$reason), call. = FALSE)
    list(sql = tr$sql, inputs = intersect(all.vars(x), names(types)))
}

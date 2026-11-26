## Fixed effects on remote data: design_lm(y ~ x | fe, con, table = ...).
##
## The estimator is that of fe.R (dm_fe_core()); this file is its remote
## backend. The rows used are copied once to a temporary table, with the
## FE (and cluster) codes as dense ranks; singletons are deleted there
## iteratively (DELETE ... WHERE code IN (... HAVING COUNT(*) = 1) until no
## row goes), and the codes ranked again. The design is a remote design on
## that table (SQL or chunked engine), and the sums of the estimator are:
##
##  * X'WX, X'Wy, y'Wy: the moments of the remote design (dm_sql_moments(),
##    dm_chunked_moments());
##  * S_k = F_k'WX: colSums(w * D, group = code_k) of the remote design;
##  * counts n_k, s_k = F_k'Wy, crosstabs N_kl, total sum of squares:
##    GROUP BY queries;
##  * residuals e = y - X b - sum_k a_k[code_k]: X b is the prediction of the
##    remote design (a view or table), the FE values temporary tables joined
##    on the codes; then X'We, F_k'We, sum w e^2, the cluster scores
##    colSums(w e X, group = cl) and the crosstabs of w e over (cluster,
##    FE_k) are again products of the remote design and GROUP BY queries.
##
## Numeric variables are centered when that is exact (dm_center(), as for
## remote fits without FE). The response must have a SQL translation.

dm_fe_remote <- function(pf, source, w, contrasts, tol, clv, fixef.rm, dots) {
    dm_need_duckdb()
    nms <- names(dots)
    if (is.null(nms)) nms <- rep("", length(dots))
    keys <- c("table", "sql", "engine", "chunk_size")
    if (any(!nms %in% keys))
        stop("Unknown arguments: ", paste(nms[!nms %in% keys], collapse = ", "), call. = FALSE)
    if (length(attr(stats::terms(pf$main), "offset")))
        stop("Offsets are not supported with fixed effects.", call. = FALSE)
    ## the rows used by the main formula (and the weights)
    D0 <- suppressMessages(dm_remote(pf$main, source, contrasts, NULL, dots, weights = w, drop_unused = TRUE))
    if (is.null(D0$response)) stop("The formula has no response.", call. = FALSE)
    if (D0$vars[[D0$response]]$kind != "numeric")
        stop("The response of a remote fit must be one numeric variable with a SQL translation.", call. = FALSE)
    fe_vars <- unique(unlist(pf$fe, use.names = FALSE))
    miss <- setdiff(c(fe_vars, clv), names(D0$src$types))
    if (length(miss)) stop(sprintf("No column `%s` in the remote data.", miss[1L]), call. = FALSE)
    nall <- D0$n
    tab <- dm_fe_remote_table(D0, pf, fe_vars, w, clv, fixef.rm)
    con <- dm_con(D0)
    rest <- dots[!nms %in% c("table", "sql")]
    D <- dm_remote(pf$main, con, contrasts, NULL, c(list(table = tab$name), rest),
                   weights = w, drop_unused = TRUE)
    ## the temporary table lives on D0's connection (a private one for Parquet)
    D$src$holder$parent <- D0$src$holder
    nfe <- length(pf$fe)
    fecols <- paste0("__dm_fe", seq_len(nfe))
    clcol <- NULL
    Gc <- NULL
    if (!is.null(clv)) {
        m <- match(clv, names(pf$fe))
        clcol <- if (is.na(m)) "__dm_cl" else fecols[m]
    }
    gmax <- function(col) as.integer(dm_query(D, sprintf("SELECT MAX(%s) AS g FROM %s", dm_qi(D, col),
                                                         dm_qi(D, tab$name)))$g)
    G <- vapply(fecols, gmax, 0L, USE.NAMES = FALSE)
    if (!is.null(clcol)) {
        Gc <- gmax(clcol)
        if (Gc < 2L) stop("Clustered standard errors need at least two clusters.", call. = FALSE)
    }
    lv <- lapply(seq_len(nfe), function(k) {
        vars <- pf$fe[[k]]
        r <- dm_query(D, sprintf("SELECT DISTINCT %s, %s FROM %s ORDER BY 1", dm_qi(D, fecols[k]),
                                 paste(dm_qi(D, vars), collapse = ", "), dm_qi(D, tab$name)))
        out <- r[vars]
        for (v in vars) if (is.factor(out[[v]])) out[[v]] <- as.character(out[[v]])
        out
    })
    ctr <- if (isTRUE(getOption("design.matrix.fe_center", TRUE)) &&
               isTRUE(getOption("design.matrix.remote_center", TRUE))) dm_center(D) else NULL
    D2 <- NULL
    if (!is.null(ctr)) {
        D2 <- D
        D2$cache <- new.env(parent = emptyenv())
        D2$id <- dm_tmp_name("dm")
        if (D$engine == "sql") {
            for (nm in names(ctr$shift))
                D2$vars[[nm]]$sql <- sprintf("(%s - %s)", D$vars[[nm]]$sql, dm_sql_double(ctr$shift[[nm]]))
        } else {
            D2$shift <- ctr$shift
        }
    }
    est <- dm_fe_estimate(function(X) dm_fe_remote_backend(X, w, fecols, G, clcol, Gc),
                          D, ctr, D2, tol, !is.null(clv))
    fit <- dm_fe_finish(est, D$colnames, pf, lv, clv, tab$nsingle, NULL, nall)
    if (!is.null(w)) {
        fit$weighted <- TRUE
        fit$sumlogw <- as.double(dm_query(D, sprintf("SELECT SUM(ln(CAST(%s AS DOUBLE))) AS s FROM %s",
                                                      dm_qi(D, w), dm_qi(D, tab$name)))$s)
    }
    fit$terms <- D$terms_full
    fit$contrasts <- D$plan$contrasts
    fit$xlevels <- D$plan$xlevels
    fit$spec <- D$plan
    fit$remote <- TRUE
    fit$remote_design <- D
    fit$remote_weights <- w
    fit$fe_cols <- fecols
    fit$method <- if (D$engine == "sql") "SQL aggregates" else "aggregates computed in batches"
    fit
}

## The temporary table of the rows used: the columns the model reads, the
## FE codes __dm_fe<k> (dense ranks of the FE variables) and, for a cluster
## variable that is not an FE, its codes __dm_cl; singletons deleted.
dm_fe_remote_table <- function(D0, pf, fe_vars, w, clv, fixef.rm) {
    q <- function(x) dm_qi(D0, x)
    cols <- unique(c(dm_remote_inputs(D0), fe_vars, clv, w))
    types <- D0$src$types
    present <- vapply(unique(c(fe_vars, clv)), function(col) {
        if (types[[col]] == "double") sprintf("s.%s IS NOT NULL AND NOT isnan(s.%s)", q(col), q(col))
        else sprintf("s.%s IS NOT NULL", q(col))
    }, "")
    where <- c(sprintf("(%s)", D0$filter), present,
               if (!is.null(w)) sprintf("CAST(s.%s AS DOUBLE) <> 0", q(w)))
    nfe <- length(pf$fe)
    fecols <- paste0("__dm_fe", seq_len(nfe))
    ranks <- vapply(seq_len(nfe), function(k) {
        sprintf("DENSE_RANK() OVER (ORDER BY %s) AS %s", paste0("s.", q(pf$fe[[k]]), collapse = ", "), q(fecols[k]))
    }, "")
    clrank <- if (!is.null(clv) && is.na(match(clv, names(pf$fe))))
        sprintf("DENSE_RANK() OVER (ORDER BY s.%s) AS %s", q(clv), q("__dm_cl"))
    name <- dm_tmp_name("dm_fe")
    dm_execute(D0, sprintf("CREATE TEMPORARY TABLE %s AS SELECT %s FROM (%s) AS s WHERE %s",
                           q(name), paste(c(paste0("s.", q(cols)), ranks, clrank), collapse = ", "),
                           D0$src$select, paste(where, collapse = " AND ")))
    nsingle <- 0
    if (fixef.rm == "singleton") {
        conds <- vapply(fecols, function(col) {
            sprintf("%s IN (SELECT %s FROM %s GROUP BY %s HAVING COUNT(*) = 1)", q(col), q(col), q(name), q(col))
        }, "")
        repeat {
            k <- dm_execute(D0, sprintf("DELETE FROM %s WHERE %s", q(name), paste(conds, collapse = " OR ")))
            if (!k) break
            nsingle <- nsingle + k
        }
        if (nsingle) {
            rr <- c(sprintf("DENSE_RANK() OVER (ORDER BY %s) AS %s", q(fecols), q(fecols)),
                    if (!is.null(clrank)) sprintf("DENSE_RANK() OVER (ORDER BY %s) AS %s", q("__dm_cl"), q("__dm_cl")))
            name2 <- dm_tmp_name("dm_fe")
            dm_execute(D0, sprintf("CREATE TEMPORARY TABLE %s AS SELECT * REPLACE (%s) FROM %s",
                                   q(name2), paste(rr, collapse = ", "), q(name)))
            dm_execute(D0, sprintf("DROP TABLE %s", q(name)))
            name <- name2
        }
    }
    n <- as.numeric(dm_query(D0, sprintf("SELECT COUNT(*) AS n FROM %s", q(name)))$n)
    if (n == 0) stop("No observation left after removing singletons.", call. = FALSE)
    list(name = name, nsingle = as.integer(nsingle))
}

## The remote design over the same rows, with `select` as its source (a
## query giving the original source's columns, plus others).
dm_fe_with_select <- function(D, select) {
    D$src$select <- select
    D$src$proto <- DBI::dbGetQuery(dm_con(D), sprintf("SELECT * FROM (%s) AS s LIMIT 0", select))
    D$src$types <- vapply(D$src$proto, dm_col_type, "")
    D$src$base_table <- NULL
    D
}

## The remote design with the columns `.fitted` (X b plus the FE values),
## `.resid` (y - `.fitted`) and `__dm_we` (w times `.resid`). The FE values
## are written to temporary tables (code, v) and joined on the codes.
dm_fe_resid_design <- function(D, beta, fe, fecols, w) {
    taken <- intersect(c(".fitted", ".resid"), names(D$src$types))
    if (length(taken)) stop(sprintf("The remote data already has a column `%s`.", taken[1L]), call. = FALSE)
    con <- dm_con(D)
    q <- function(x) dm_qi(D, x)
    view <- dm_remote_predict(D, beta, column = "__dm_xb")
    joins <- character()
    terms <- character()
    for (k in seq_along(fe)) {
        tn <- dm_tmp_name("dm_fev")
        DBI::dbWriteTable(con, tn, data.frame(code = seq_along(fe[[k]]), v = unname(as.double(fe[[k]]))),
                          temporary = TRUE)
        joins <- c(joins, sprintf("LEFT JOIN %s AS f%d ON s.%s = f%d.code", q(tn), k, q(fecols[k]), k))
        terms <- c(terms, sprintf("f%d.v", k))
    }
    fit <- paste(c("s.\"__dm_xb\"", terms), collapse = " + ")
    ysql <- D$vars[[D$response]]$sql
    e <- sprintf("(%s - (%s))", ysql, fit)
    we <- if (is.null(w)) e else sprintf("%s * CAST(s.%s AS DOUBLE)", e, q(w))
    select <- sprintf("SELECT s.*, %s AS \".fitted\", %s AS \".resid\", %s AS \"__dm_we\" FROM (SELECT * FROM %s) AS s %s",
                      fit, e, we, q(unclass(view)[1L]), paste(joins, collapse = " "))
    out <- dm_fe_with_select(D, select)
    out$src$view <- view
    out
}

dm_fe_remote_backend <- function(D, w, fecols, G, clcol, Gc) {
    q <- function(x) dm_qi(D, x)
    yy <- list(cols = D$response, response = TRUE)
    sql <- D$engine == "sql"
    ysql <- D$vars[[D$response]]$sql
    wsql <- if (is.null(w)) "CAST(1 AS DOUBLE)" else sprintf("CAST(s.%s AS DOUBLE)", q(w))
    from <- function(Dx) sprintf("(%s) AS s WHERE %s", Dx$src$select, Dx$filter)
    m <- if (sql) dm_sql_moments(D, yy, w) else dm_chunked_moments(D, yy, w)
    if (!all(is.finite(m$xtx)) || !all(is.finite(m$xty))) stop("NA/NaN/Inf in 'x'", call. = FALSE)
    tot <- dm_query(D, sprintf("SELECT SUM(%s) AS sw, SUM(%s * %s) AS swy FROM %s", wsql, wsql, ysql, from(D)))
    ybar <- as.double(tot$swy) / as.double(tot$sw)
    tss <- as.double(dm_query(D, sprintf("SELECT SUM(%s * pow(%s - %s, 2)) AS t FROM %s", wsql, ysql,
                                         dm_sql_double(ybar), from(D)))$t)
    mom <- list(xtx = unname(m$xtx), xty = unname(m$xty[, 1L]), yty = m$yty[1L, 1L], tss = tss)
    ## grouped sums by the codes of `col` (1..Gk, all present)
    gsum <- function(Dx, col, Gk, exprs) {
        r <- dm_query(Dx, sprintf("SELECT %s AS g, %s FROM %s GROUP BY 1", q(col),
                                  paste(sprintf("SUM(%s) AS v%d", exprs, seq_along(exprs)), collapse = ", "),
                                  from(Dx)))
        out <- matrix(0, Gk, length(exprs))
        for (j in seq_along(exprs)) out[r$g, j] <- as.double(r[[paste0("v", j)]])
        out
    }
    colsums <- function(Dx, col, Gk, wcol) {
        S <- if (sql) dm_sql_colsums(Dx, col, wcol, FALSE) else dm_chunked_colsums(Dx, col, wcol, FALSE)
        out <- matrix(0, Gk, ncol(S))
        out[as.integer(rownames(S)), ] <- unname(S)
        out
    }
    gs <- lapply(seq_along(fecols), function(k) {
        ns <- gsum(D, fecols[k], G[k], c(wsql, sprintf("%s * %s", wsql, ysql)))
        list(nk = ns[, 1L], S = colsums(D, fecols[k], G[k], w), s = ns[, 2L])
    })
    crosstab <- function(Dx, a, b, Ga, Gb, wexpr) {
        r <- dm_query(Dx, sprintf("SELECT %s AS a, %s AS b, SUM(%s) AS v FROM %s GROUP BY 1, 2",
                                  q(a), q(b), wexpr, from(Dx)))
        Matrix::sparseMatrix(i = as.integer(r$b), j = as.integer(r$a), x = as.double(r$v),
                             dims = c(as.integer(Gb), as.integer(Ga)), repr = "C")
    }
    state <- new.env(parent = emptyenv())
    be <- list(n = m$n, p = D$p, icol = which(D$colnames == "(Intercept)"), G = G, moments = mom, gs = gs)
    be$crosstab <- function(a, b) crosstab(D, fecols[a], fecols[b], G[a], G[b], wsql)
    be$resid <- function(beta, fe) {
        De <- dm_fe_resid_design(D, beta, fe, fecols, w)
        state$De <- De
        ry <- list(cols = ".resid", response = FALSE)
        xte <- if (sql) dm_sql_xty(De, ry, w) else dm_chunked_moments(De, ry, w, gram = FALSE)$xty
        Be <- lapply(seq_along(fecols), function(k) gsum(De, fecols[k], G[k], "s.\"__dm_we\"")[, 1L])
        rss <- as.double(dm_query(De, sprintf("SELECT SUM(s.\"__dm_we\" * s.\".resid\") AS r FROM %s",
                                              from(De)))$r)
        list(xte = unname(xte[, 1L]), Be = Be, rss = rss)
    }
    be$cluster <- function() {
        De <- state$De
        U <- colsums(De, clcol, Gc, "__dm_we")
        nested <- vapply(seq_along(fecols), function(k) {
            if (identical(fecols[k], clcol)) return(TRUE)
            r <- dm_query(D, sprintf("SELECT COUNT(*) AS n FROM (SELECT %s FROM %s GROUP BY 1 HAVING COUNT(DISTINCT %s) > 1) AS t",
                                     q(fecols[k]), from(D), q(clcol)))
            as.numeric(r$n) == 0
        }, NA)
        list(G = Gc, U = U, nested = nested,
             E = function(k) crosstab(De, clcol, fecols[k], Gc, G[k], "s.\"__dm_we\""))
    }
    be
}

## residuals() of a remote FE fit: the remote design with the columns
## `.fitted` (including the FE) and `.resid`.
dm_fe_remote_residuals <- function(object) {
    D <- object$remote_design
    b <- object$coefficients
    b[is.na(b)] <- 0
    beta <- numeric(D$p)
    beta[match(names(b), D$colnames)] <- b
    dm_fe_resid_design(D, beta, object$fixef, object$fe_cols, object$remote_weights)
}

## predict(fit, newdata = <remote>): X b plus the FE values of the rows'
## levels (NULL for levels not in the fit), in a temporary view whose name
## is returned (as dm_remote_predict() does).
dm_fe_remote_predict <- function(object, newdata, dots) {
    if (is.null(object$spec)) stop("This design cannot be applied to new data.", call. = FALSE)
    Dn <- do.call(design_matrix, c(list(object$spec, newdata, na.action = stats::na.pass), dots))
    con <- dm_con(Dn)
    q <- function(x) dm_qi(Dn, x)
    b <- object$coefficients
    b[is.na(b)] <- 0
    beta <- numeric(Dn$p)
    beta[match(names(b), Dn$colnames)] <- b
    if (".fitted" %in% names(Dn$src$types))
        stop("The remote data already has a column `.fitted`.", call. = FALSE)
    view <- dm_remote_predict(Dn, beta, column = "__dm_xb")
    joins <- character()
    terms <- character()
    for (k in seq_along(object$fixef)) {
        vars <- object$fe_vars[[k]]
        miss <- setdiff(vars, names(Dn$src$types))
        if (length(miss)) stop(sprintf("`newdata` has no fixed-effect variable `%s`.", miss[1L]), call. = FALSE)
        lv <- object$fe_levels[[k]]
        lcols <- paste0("__dm_l", seq_along(vars))
        tab <- stats::setNames(lv, lcols)
        tab$v <- unname(as.double(object$fixef[[k]]))
        tn <- dm_tmp_name("dm_fel")
        DBI::dbWriteTable(con, tn, tab, temporary = TRUE)
        on <- vapply(seq_along(vars), function(j) {
            if (is.character(lv[[j]])) sprintf("CAST(s.%s AS VARCHAR) = f%d.%s", q(vars[j]), k, q(lcols[j]))
            else sprintf("s.%s = f%d.%s", q(vars[j]), k, q(lcols[j]))
        }, "")
        joins <- c(joins, sprintf("LEFT JOIN %s AS f%d ON %s", q(tn), k, paste(on, collapse = " AND ")))
        terms <- c(terms, sprintf("f%d.v", k))
    }
    name <- dm_tmp_name("dm_pred")
    dm_execute(Dn, sprintf("CREATE TEMPORARY VIEW %s AS SELECT s.* EXCLUDE (\"__dm_xb\"), %s AS \".fitted\" FROM (SELECT * FROM %s) AS s %s",
                           q(name), paste(c("s.\"__dm_xb\"", terms), collapse = " + "),
                           q(unclass(view)[1L]), paste(joins, collapse = " ")))
    structure(name, class = "dm_remote_table", connection = con, holder = Dn$src$holder)
}

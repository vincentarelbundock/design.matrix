## The primitive layer of the fitting core (Kata design.matrix#vavc).
##
## Fitting code (linear models, and on this layer GLMs and fixed effects) is
## written once against the operations below, which take a *design* and
## *columns*, and dispatch to two backends:
##
##   design                     column
##   ------------------------   -----------------------------------------
##   in memory: design_matrix   numeric vector (or matrix: several columns)
##   remote:    dm_remote       dm_rcols: SQL expressions over the rows of a
##                              frame (see below)
##
## Nothing here is exported. All functions are internal; their names start
## with dmp_ (operations on designs), dmc_ (operations on columns) or dm_.
##
## ---------------------------------------------------------------------------
## Columns
## ---------------------------------------------------------------------------
##
## A column holds one value per row of a design (in the order of its rows in
## memory; remotely, rows have no order). In memory it is a double vector of
## length nrow(x), or an n x k matrix for k columns. Remotely it is
##
##   dm_rcols(frame, sql)
##
## with `frame` a dm_remote (its source `s` holds what the expressions read,
## its filter selects its rows) and `sql` a character vector, one SQL
## expression over `s` per column (k = length(sql)). Columns are made by
##
##   dmp_column(x, sql)    a column of design x from SQL over its source
##                         (in memory: `sql` is returned, i.e. pass a vector)
##   dmp_mult(x, B)        products with the design (below)
##   arithmetic on columns
##
## Operations on columns (same R code for both backends):
##   +  -  *  /  ^         (numbers recycle; columns must have the same k or 1)
##   ==  !=  <  >  <=  >=  &  |  !
##   exp log sqrt abs log1p expm1 log2 log10 cos sin tan floor ceiling sign
##                         (Math group; log has one argument)
##   dmc_pmin dmc_pmax dmc_ifelse dmc_rowsums dmc_cbind
##   dmc_pnorm dmc_qnorm dmc_dnorm dmc_pcauchy
##                         (remotely the DuckDB macros of dm_sql_macros(), as
##                         R computes them: pnorm/dnorm/pcauchy to rounding,
##                         qnorm to one ulp; dmc_mexpm1m(u) = -expm1(-u))
##
## NA semantics: SQL NULL is R's NA, IEEE NaN is NaN; arithmetic propagates
## both as R does (NULL op x is NULL). Sums over rows (dm_colsum) return NA
## when a value is NULL, as sum() does, and NaN when a value is NaN.
## Comparisons with NULL are NULL (R: NA).
##
## Columns from successive products live on successive frames: each product
## builds a new frame (a view, or a table with the chunked engine) on top of
## the frame of the columns passed as `along`, and keeps all their columns.
## Combining columns takes the later frame (they must come from the same
## design). In memory `along` is ignored.
##
## Getting values into R: dmp_collect(col) returns the values (in memory the
## vector itself; remotely a SELECT of all rows: only for small data).
##
## ---------------------------------------------------------------------------
## Operations on designs
## ---------------------------------------------------------------------------
##
##   dmp_scale(x, s)          the design with rows scaled by column s (s * x);
##                            scales multiply. Same rows and columns as x.
##   dmp_crossprod(x)         X'X, p x p matrix
##   dmp_crossprod(x, y)      y a design (x itself, possibly scaled: X'SY);
##                            or a column (k columns): X'c, p x k matrix
##   dmp_colsums(x, group)    column sums of x (possibly scaled): vector of
##                            length p; with group (a column, values used as
##                            group labels; NA a group of its own, with
##                            rowsum()'s warning): G x p matrix, groups sorted
##                            as rowsum() sorts them, row names the labels
##   dmp_groupsum(col, group) sums of a column (k columns) within groups:
##                            G x k matrix, as rowsum(col, group)
##   dmp_crosstab(a, b, w)    for two grouping columns, the sums of w (counts
##                            without w) for each pair of groups: a
##                            Matrix::dgCMatrix, rows the sorted groups of a,
##                            columns those of b
##   dmp_mult(x, B, along)    X B (B a p-vector or p x k matrix) as a column
##   dmp_quad(x, V, along)    rowSums((X V) * X), V p x p symmetric positive
##                            semi-definite, as a column
##   dm_colsum(col)           sums over rows: a number per column
##   dmp_colnorms(col)        Euclidean norms of the columns
##   dm_nrows(x)              number of rows of a design or column
##   dmp_patterns(x, z, w, compress)
##                            compressed least squares: when rows with the
##                            same design inputs are worth grouping, a list
##                            with x (the design of the patterns, in memory),
##                            z (the weighted mean of the column z per
##                            pattern), w (the sum of the weights per pattern,
##                            NULL when all 1), within (the weighted sum of
##                            squares of z around the pattern means) and rows
##                            (the number of rows); otherwise NULL. compress:
##                            "auto", TRUE or FALSE. A weighted least squares
##                            fit of z on x with weights w gives the fit on
##                            every row; its residual sum of squares plus
##                            `within` is that of the fit on every row.
##   dmp_output(x, cols)      row-wise results of an algorithm: in memory a
##                            data frame of the columns (one row per row of x);
##                            remotely a temporary view with the columns of
##                            the remote data and the named columns `cols`,
##                            a string of class dm_remote_table. One code path
##                            computes the columns; only the container
##                            differs.
##
## Remote methods: crossprod, colSums, products and patterns run in the
## database (or in batches with the chunked engine); dmp_groupsum and
## dmp_crosstab are single GROUP BY queries.
##
## ---------------------------------------------------------------------------
## User-facing operations of a remote design (operator goal: every matrix
## operation of an in-memory design works on a remote one, through this
## layer; status as of this contract)
## ---------------------------------------------------------------------------
##
##   operation                    remote meaning                    status
##   dim, dimnames                rows, columns                     done
##   as.matrix(D, force)          collected; refused above          done
##                                remote_max_rows unless force
##   D %*% b, D %*% B             table: data columns + predictions done
##   crossprod(D), (D, s * D),    p x p / p x k matrices            done
##   (D, "col" / SQL() / column)
##   crossprod(D, numeric w)      no row alignment remotely: error; done (error)
##                                use a column or s * D
##   colSums/colMeans(s * D,      vectors / G x p matrices          done
##     group = col)
##   rowSums/rowMeans(D)          table (column .rowsum / .rowmean) to do
##   quadratic_diagonal(D, V)     table (column .quad), dmp_quad    to do
##   qr(D)                        Cholesky of X'X with lm()'s rule, to do (WIP)
##                                suspects judged from products
##   qr.coef, qr.R, solve         numbers (as in memory)            to do
##   qr.fitted, qr.resid, qr.qy,  tables (row-wise results)         to do
##     qr.Q
##   qr.qty                       as in memory (p values)           to do
##   as(D, "CsparseMatrix")       collected, same guard as          to do
##                                as.matrix
##   D[i, ]                       i a logical remote column or a    to do
##                                DBI::SQL() predicate: a remote
##                                design over those rows; numeric
##                                indices refused (no row order)
##   s * D, D / s, -D             row scaling (number, column       done
##                                name, SQL(), column)
##   D1 + D2, D1 - D2             proposal: refused (combinations of to decide
##                                two remote designs)
##   finite_rows(D) (internal)    a logical column                  to do

## ---- remote columns ----

dm_rcols <- function(frame, sql) {
    ## (parenthesized, so that expressions compose with any operator)
    sql <- as.character(sql)
    sql <- sprintf("(%s)", sql)
    structure(list(frame = frame, sql = sql), class = "dm_rcols")
}

dm_is_rcols <- function(x) inherits(x, "dm_rcols")

#' @export
length.dm_rcols <- function(x) length(unclass(x)$sql)

dm_frame_depth <- function(D) D$depth %||% 0L

## Of two frames of the same design, the later one (which keeps the columns
## of the earlier one).
dm_frame_later <- function(a, b) {
    if (is.null(a)) return(b)
    if (is.null(b)) return(a)
    if (!identical(a$id, b$id))
        stop("Columns of different remote designs cannot be combined.", call. = FALSE)
    if (dm_frame_depth(b) > dm_frame_depth(a)) b else a
}

dm_rcols_frame <- function(...) {
    f <- NULL
    for (v in list(...)) if (dm_is_rcols(v)) f <- dm_frame_later(f, unclass(v)$frame)
    f
}

## SQL of an operand (a remote column, or a number)
dm_rcols_operand <- function(x) {
    if (dm_is_rcols(x)) return(unclass(x)$sql)
    if (is.logical(x) && length(x) == 1L && !is.na(x)) return(if (x) "TRUE" else "FALSE")
    if (is.numeric(x) && length(x) == 1L) return(dm_sql_double(x))
    stop("Remote columns combine only with numbers and other remote columns.", call. = FALSE)
}

#' @export
Ops.dm_rcols <- function(e1, e2) {
    op <- .Generic
    if (nargs() == 1L) {
        a <- dm_rcols_operand(e1)
        sql <- switch(op, "-" = sprintf("(-%s)", a), "+" = a, "!" = sprintf("(NOT %s)", a),
                      stop(sprintf("`%s` is not supported for remote columns.", op), call. = FALSE))
        return(dm_rcols(unclass(e1)$frame, sql))
    }
    frame <- dm_rcols_frame(e1, e2)
    a <- dm_rcols_operand(e1)
    b <- dm_rcols_operand(e2)
    k <- max(length(a), length(b))
    if (!(length(a) %in% c(1L, k)) || !(length(b) %in% c(1L, k)))
        stop("Remote columns with different numbers of columns.", call. = FALSE)
    a <- rep_len(a, k)
    b <- rep_len(b, k)
    sql <- switch(op,
        "+" = , "-" = , "*" = , "/" = sprintf("(%s %s %s)", a, op, b),
        "^" = sprintf("pow(%s, %s)", a, b),
        "==" = sprintf("(%s = %s)", a, b), "!=" = sprintf("(%s <> %s)", a, b),
        "<" = , ">" = , "<=" = , ">=" = sprintf("(%s %s %s)", a, op, b),
        "&" = sprintf("(%s AND %s)", a, b), "|" = sprintf("(%s OR %s)", a, b),
        stop(sprintf("`%s` is not supported for remote columns.", op), call. = FALSE))
    dm_rcols(frame, sql)
}

#' @export
Math.dm_rcols <- function(x, ...) {
    a <- unclass(x)$sql
    sql <- switch(.Generic,
        exp = sprintf("exp(%s)", a), log = sprintf("ln(%s)", a), sqrt = sprintf("sqrt(%s)", a),
        abs = sprintf("abs(%s)", a), log1p = sprintf("ln(1 + %s)", a), expm1 = sprintf("(exp(%s) - 1)", a),
        log2 = sprintf("log2(%s)", a), log10 = sprintf("log10(%s)", a),
        cos = sprintf("cos(%s)", a), sin = sprintf("sin(%s)", a), tan = sprintf("tan(%s)", a),
        floor = sprintf("floor(%s)", a), ceiling = sprintf("ceil(%s)", a), sign = sprintf("sign(%s)", a),
        stop(sprintf("%s() is not supported for remote columns.", .Generic), call. = FALSE))
    dm_rcols(unclass(x)$frame, sql)
}

dmc_pmax <- function(x, y) {
    f <- dm_rcols_frame(x, y)
    if (is.null(f)) return(pmax(x, y))
    dm_rcols(f, sprintf("greatest(%s, %s)", dm_rcols_operand(x), dm_rcols_operand(y)))
}

dmc_pmin <- function(x, y) {
    f <- dm_rcols_frame(x, y)
    if (is.null(f)) return(pmin(x, y))
    dm_rcols(f, sprintf("least(%s, %s)", dm_rcols_operand(x), dm_rcols_operand(y)))
}

dmc_ifelse <- function(test, yes, no) {
    f <- dm_rcols_frame(test, yes, no)
    if (is.null(f)) return(ifelse(test, yes, no))
    dm_rcols(f, sprintf("(CASE WHEN %s THEN %s ELSE %s END)", dm_rcols_operand(test),
                        dm_rcols_operand(yes), dm_rcols_operand(no)))
}

## Columns side by side (cbind).
dmc_cbind <- function(...) {
    args <- list(...)
    f <- do.call(dm_rcols_frame, args)
    if (is.null(f)) return(do.call(cbind, args))
    dm_rcols(f, unlist(lapply(args, dm_rcols_operand)))
}

dmc_rowsums <- function(x) {
    if (!dm_is_rcols(x)) return(if (is.matrix(x)) rowSums(x) else x)
    dm_rcols(unclass(x)$frame, sprintf("(%s)", paste(unclass(x)$sql, collapse = " + ")))
}

## distribution functions: stats in memory, DuckDB macros remotely
dmc_macro <- function(x, fun, macro) {
    if (!dm_is_rcols(x)) return(fun(x))
    dm_sql_macros(unclass(x)$frame)
    dm_rcols(unclass(x)$frame, sprintf("%s(%s)", macro, unclass(x)$sql))
}
dmc_pnorm <- function(x) dmc_macro(x, stats::pnorm, "dm_pnorm")
dmc_qnorm <- function(x) dmc_macro(x, stats::qnorm, "dm_qnorm")
dmc_dnorm <- function(x) dmc_macro(x, stats::dnorm, "dm_dnorm")
dmc_pcauchy <- function(x) dmc_macro(x, stats::pcauchy, "dm_pcauchy")
dmc_mexpm1m <- function(x) dmc_macro(x, function(u) -expm1(-u), "dm_mexpm1m")

## The values of a column in R (remotely: all its rows; for small data).
dmp_collect <- function(x) {
    if (!dm_is_rcols(x)) return(x)
    fr <- unclass(x)$frame
    sql <- unclass(x)$sql
    res <- dm_query(fr, sprintf("SELECT %s FROM (%s) AS s WHERE %s",
                                paste(sprintf("CAST(%s AS DOUBLE) AS c%d", sql, seq_along(sql)), collapse = ", "),
                                fr$src$select, fr$filter))
    if (length(sql) == 1L) res[[1L]] else as.matrix(res)
}

## ---- sums over rows ----

dm_colsum <- function(x) {
    if (!dm_is_rcols(x)) {
        if (is.matrix(x)) return(colSums(x))
        return(sum(x))
    }
    fr <- unclass(x)$frame
    sql <- unclass(x)$sql
    k <- length(sql)
    items <- c(sprintf("SUM(CAST(%s AS DOUBLE)) AS s%d", sql, seq_len(k)),
               sprintf("SUM(CASE WHEN (%s) IS NULL THEN 1 ELSE 0 END) AS n%d", sql, seq_len(k)))
    res <- dm_query(fr, sprintf("SELECT %s FROM (%s) AS s WHERE %s", paste(items, collapse = ", "),
                                fr$src$select, fr$filter))
    vapply(seq_len(k), function(j) {
        v <- as.double(res[[paste0("s", j)]])
        if (as.numeric(res[[paste0("n", j)]]) > 0) return(NA_real_)
        if (is.na(v)) 0 else v
    }, 0)
}

dm_nrows <- function(x) {
    if (inherits(x, "dm_remote")) return(x$cache$n %||% x$n)
    if (dm_is_rcols(x)) return(dm_nrows(unclass(x)$frame))
    NROW(x)
}

dmp_colnorms <- function(x) {
    if (!dm_is_rcols(x)) return(.Call(C_dm_col_norms, if (is.matrix(x)) x else matrix(x, ncol = 1L)))
    sqrt(dm_colsum(x * x))
}

## ---- designs ----

dmp_column <- function(x, sql) {
    if (is.null(sql)) return(NULL)
    if (inherits(x, "dm_remote")) return(dm_rcols(x, sql))
    sql
}

## The remote design over the frame of columns (same plan, rows and scale).
dm_rebase <- function(x, frame) {
    if (is.null(frame) || identical(x$src$select, frame$src$select)) return(x)
    if (!identical(x$id, frame$id)) stop("Columns of another remote design.", call. = FALSE)
    x$src <- frame$src
    x$depth <- dm_frame_depth(frame)
    x
}

dm_along_frame <- function(x, along) {
    f <- do.call(dm_rcols_frame, as.list(along))
    if (is.null(f)) x else dm_rebase(x, f)
}

dmp_scale <- function(x, s) {
    if (is.null(s)) return(x)
    if (!dm_is_rcols(s)) return(s * x)
    x <- dm_rebase(x, unclass(s)$frame)
    x$scale <- dm_remote_scale_sql(x, DBI::SQL(unclass(s)$sql), "*")
    x
}

## A remote design with the columns `col` added to its source (named nm).
dm_with_columns <- function(x, col, nm) {
    fr <- unclass(col)$frame
    if (dm_frame_depth(x) > dm_frame_depth(fr)) fr <- x
    x <- dm_rebase(x, fr)
    sc <- x$scale
    x <- dm_remote_extend(dm_remote_unscaled(x), stats::setNames(unclass(col)$sql, nm), "column")
    x$scale <- sc
    x
}

dmp_crossprod <- function(x, y = NULL) {
    if (!inherits(x, "dm_remote")) return(if (is.null(y)) crossprod(x) else crossprod(x, y))
    if (is.null(y)) return(crossprod(x))
    if (inherits(y, "dm_remote")) {
        if (dm_frame_depth(y) > dm_frame_depth(x)) x <- dm_rebase(x, y)
        if (dm_frame_depth(x) > dm_frame_depth(y)) y <- dm_rebase(y, x)
        return(crossprod(x, y))
    }
    if (!dm_is_rcols(y)) stop("`y` must be a design or a column.", call. = FALSE)
    nm <- paste0("__dm_y", seq_along(y))
    out <- crossprod(dm_with_columns(x, y, nm), nm)
    dimnames(out) <- list(colnames(x), NULL)
    out
}

dmp_colsums <- function(x, group = NULL) {
    if (!inherits(x, "dm_remote")) {
        if (is.null(group)) return(colSums(x))
        return(colSums(x, group = group))
    }
    if (is.null(group)) return(colSums(x))
    if (!dm_is_rcols(group)) stop("`group` of a remote design must be a remote column.", call. = FALSE)
    colSums(dm_with_columns(x, group, "__dm_g"), group = "__dm_g")
}

## SQL GROUP BY over the frame of columns: list(groups (a data frame of the
## group values, sorted as rowsum() sorts them), sums (matrix)).
dm_group_query <- function(groups, values) {
    fr <- do.call(dm_rcols_frame, c(groups, values))
    gs <- vapply(groups, function(g) unclass(g)$sql, "")
    vs <- unlist(lapply(values, function(v) unclass(v)$sql))
    items <- c(sprintf("%s AS g%d", gs, seq_along(gs)),
               if (length(vs)) sprintf("SUM(CAST(%s AS DOUBLE)) AS v%d", vs, seq_along(vs)) else "COUNT(*) AS v1")
    res <- dm_query(fr, sprintf("SELECT %s FROM (%s) AS s WHERE %s GROUP BY ALL",
                                paste(items, collapse = ", "), fr$src$select, fr$filter))
    res
}

dmp_groupsum <- function(col, group) {
    if (!dm_is_rcols(col) && !dm_is_rcols(group)) {
        out <- rowsum(if (is.matrix(col)) col else matrix(col, ncol = 1L), group, reorder = TRUE)
        return(out)
    }
    res <- dm_group_query(list(group), list(col))
    ug <- sort(unique(res$g1), na.last = TRUE, method = "quick")
    if (anyNA(ug)) warning("missing values for 'group'", call. = FALSE)
    i <- match(res$g1, ug)
    k <- length(col)
    out <- matrix(0, length(ug), k, dimnames = list(as.character(ug), NULL))
    for (j in seq_len(k)) out[i, j] <- as.double(res[[paste0("v", j)]])
    out
}

dmp_crosstab <- function(a, b, w = NULL) {
    if (!dm_is_rcols(a) && !dm_is_rcols(b) && !dm_is_rcols(w)) {
        ua <- sort(unique(a), na.last = TRUE, method = "quick")
        ub <- sort(unique(b), na.last = TRUE, method = "quick")
        return(Matrix::sparseMatrix(i = match(a, ua), j = match(b, ub),
                                    x = if (is.null(w)) rep(1, length(a)) else as.double(w),
                                    dims = c(length(ua), length(ub)),
                                    dimnames = list(as.character(ua), as.character(ub))))
    }
    res <- dm_group_query(list(a, b), if (is.null(w)) list() else list(w))
    ua <- sort(unique(res$g1), na.last = TRUE, method = "quick")
    ub <- sort(unique(res$g2), na.last = TRUE, method = "quick")
    Matrix::sparseMatrix(i = match(res$g1, ua), j = match(res$g2, ub), x = as.double(res$v1),
                         dims = c(length(ua), length(ub)),
                         dimnames = list(as.character(ua), as.character(ub)))
}

dmp_mult <- function(x, B, along = list()) {
    if (!inherits(x, "dm_remote")) return(x %*% B)
    x <- dm_along_frame(x, along)
    B <- if (is.matrix(B)) B else matrix(B, ncol = 1L)
    nm <- paste0("__dm_p", substr(dm_tmp_name("p"), 3L, 10L), "_", seq_len(ncol(B)))
    colnames(B) <- nm
    view <- dm_remote_predict(x, B, column = nm[1L])
    fr <- dm_remote_unscaled(x)
    fr$src$select <- sprintf("SELECT * FROM %s", dm_qi(x, unclass(view)[1L]))
    fr$src$base_table <- NULL
    fr$src$view <- c(fr$src$view, list(view))
    proto <- DBI::dbGetQuery(dm_con(x), sprintf("SELECT * FROM (%s) AS s LIMIT 0", fr$src$select))
    fr$src$proto <- proto
    fr$src$types <- vapply(proto, dm_col_type, "")
    fr$depth <- dm_frame_depth(x) + 1L
    dm_rcols(fr, sprintf("s.%s", dm_qi(x, nm)))
}

dmp_quad <- function(x, V, along = list()) {
    if (!inherits(x, "dm_remote")) return(quadratic_diagonal(x, V))
    zero <- dm_rcols(dm_along_frame(x, along), "CAST(0 AS DOUBLE)")
    if (ncol(V) == 0L) return(zero)
    e <- eigen(V, symmetric = TRUE)
    keep <- e$values > 0
    if (!any(keep)) return(zero)
    L <- e$vectors[, keep, drop = FALSE] %*% diag(sqrt(e$values[keep]), sum(keep))
    P <- dmp_mult(x, L, along)
    dmc_rowsums(P * P)
}

## Compressed least squares (see the header); the in-memory method groups
## identical rows (lm.R), the remote one is one GROUP BY (remote_compress.R).
dmp_patterns <- function(x, z, w = NULL, compress = "auto") {
    if (inherits(x, "dm_remote")) {
        mode <- if (isTRUE(compress)) TRUE else if (isFALSE(compress)) FALSE else dm_compress_mode()
        pt <- dm_remote_patterns(x, unclass(z)$sql, if (!is.null(w)) unclass(w)$sql, within = TRUE, mode = mode)
        if (is.null(pt)) return(NULL)
        pd <- dm_pattern_design(x, pt$P[dm_compress_keys(x)])
        k <- pd$keep
        sw <- pt$sw[k]
        pos <- which(sw > 0)
        ws <- sw[pos]
        if (is.null(w) && all(ws == 1)) ws <- NULL
        return(list(x = pd$design[pos, ], z = pt$sy[k, 1L][pos] / sw[pos], w = ws,
                    within = sum(pt$within[k]), rows = sum(pt$n[k])))
    }
    g <- dmfit_groups(x, compress)
    if (is.null(g)) return(NULL)
    ss <- dmfit_suffstats(z, w, g$map, length(g$first))
    ## sums that overflow (w * z beyond the range of doubles): no compression
    if (!all(is.finite(ss$w), is.finite(ss$mean), is.finite(ss$within))) return(NULL)
    keep <- ss$w > 0
    ws <- ss$w[keep]
    if (is.null(w) && all(ws == 1)) ws <- NULL
    list(x = x[g$first[keep], ], z = ss$mean[keep], w = ws, within = ss$within, rows = nrow(x))
}

## Row-wise results: a data frame in memory, a remote table remotely.
dmp_output <- function(x, cols) {
    if (!inherits(x, "dm_remote")) return(as.data.frame(lapply(cols, function(v) as.vector(v))))
    fr <- dm_rebase(dm_remote_unscaled(x), do.call(dm_rcols_frame, unname(cols)))
    srccols <- names(x$src$types)
    items <- vapply(names(cols), function(nm) {
        v <- cols[[nm]]
        sprintf("%s AS %s", if (dm_is_rcols(v)) unclass(v)$sql else dm_rcols_operand(v), dm_qi(fr, nm))
    }, "")
    sel <- sprintf("SELECT %s FROM (%s) AS s WHERE %s",
                   paste(c(paste0("s.", dm_qi(fr, srccols)), items), collapse = ", "), fr$src$select, fr$filter)
    name <- dm_tmp_name("dm_out")
    dm_execute(fr, sprintf("CREATE TEMPORARY VIEW %s AS %s", dm_qi(fr, name), sel))
    structure(name, class = "dm_remote_table", connection = dm_con(fr), holder = fr$src$holder)
}

## DuckDB macros for the distribution functions of remote columns, as R
## computes them (from the GLM work, kata/402b-glm): dm_pnorm, dm_qnorm
## (one ulp), dm_dnorm, dm_pcauchy and dm_mexpm1m(u) = -expm1(-u). Created
## once per connection.
dm_sql_macros <- function(D) {
    if (isTRUE(D$src$holder$dm_macros)) return(invisible())
    horner <- function(x, co) {
        s <- sprintf("%.17g", co[1L])
        for (k in co[-1L]) s <- sprintf("(%s * %s + %.17g)", s, x, k)
        s
    }
    a <- c(2.2352520354606839287, 161.02823106855587881, 1067.6894854603709582, 18154.981253343561249,
           0.065682337918207449113)
    b <- c(47.20258190468824187, 976.09855173777669322, 10260.932208618978205, 45507.789335026729956)
    cc <- c(0.39894151208813466764, 8.8831497943883759412, 93.506656132177855979, 597.27027639480026226,
            2494.5375852903726711, 6848.1904505362823326, 11602.651437647350124, 9842.7148383839780218,
            1.0765576773720192317e-8)
    d <- c(22.266688044328115691, 235.38790178262499861, 1519.377599407554805, 6485.558298266760755,
           18615.571640885098091, 34900.952721145977266, 38912.003286093271411, 19685.429676859990727)
    p <- c(0.21589853405795699, 0.1274011611602473639, 0.022235277870649807, 0.001421619193227893466,
           2.9112874951168792e-5, 0.02307344176494017303)
    q <- c(1.28426009614491121, 0.468238212480865118, 0.0659881378689285515, 0.00378239633202758244,
           7.29751555083966205e-5)
    num1 <- sprintf("(x * %s)", horner("x*x", c(a[5L], a[1:4])))
    den1 <- horner("x*x", c(1, b[1:4]))
    num2 <- horner("abs(x)", c(cc[9L], cc[1:8]))
    den2 <- horner("abs(x)", c(1, d[1:8]))
    ix <- "(1.0/(x*x))"
    num3 <- sprintf("(%s * %s)", ix, horner(ix, c(p[6L], p[1:5])))
    den3 <- horner(ix, c(1, q[1:5]))
    xs <- "(trunc(abs(x)*16)/16)"
    tail <- function(temp) sprintf("(exp(-%s*%s*0.5) * exp(-((abs(x) - %s)*(abs(x) + %s))*0.5) * %s)", xs, xs, xs, xs, temp)
    t2 <- tail(sprintf("(%s / %s)", num2, den2))
    t3 <- tail(sprintf("((0.398942280401432677939946059934 - %s / %s) / abs(x))", num3, den3))
    pn <- sprintf(paste("CASE WHEN abs(x) <= 0.67448975 THEN 0.5 + %s / %s",
                        "WHEN abs(x) <= 5.656854249492380195206754896838 THEN (CASE WHEN x > 0 THEN 1 - %s ELSE %s END)",
                        "WHEN abs(x) < 50 THEN (CASE WHEN x > 0 THEN 1 - %s ELSE %s END)",
                        "ELSE (CASE WHEN x > 0 THEN 1.0 ELSE 0.0 END) END"), num1, den1, t2, t2, t3, t3)
    qa <- c(2509.0809287301226727, 33430.575583588128105, 67265.770927008700853, 45921.953931549871457,
            13731.693765509461125, 1971.5909503065514427, 133.14166789178437745, 3.387132872796366608)
    qb <- c(5226.495278852545925, 28729.085735721942674, 39307.89580009271061, 21213.794301586595867,
            5394.1960214247511077, 687.1870074920579083, 42.313330701600911252, 1)
    qc <- c(7.7454501427834140764e-4, .0227238449892691845833, .24178072517745061177, 1.27045825245236838258,
            3.64784832476320460504, 5.7694972214606914055, 4.6303378461565452959, 1.42343711074968357734)
    qd <- c(1.05075007164441684324e-9, 5.475938084995344946e-4, .0151986665636164571966, .14810397642748007459,
            .68976733498510000455, 1.6763848301838038494, 2.05319162663775882187, 1)
    qe <- c(2.01033439929228813265e-7, 2.71155556874348757815e-5, .0012426609473880784386, .026532189526576123093,
            .29656057182850489123, 1.7848265399172913358, 5.4637849111641143699, 6.6579046435011037772)
    qf <- c(2.04426310338993978564e-15, 1.4215117583164458887e-7, 1.8463183175100546818e-5, 7.868691311456132591e-4,
            .0148753612908506148525, .13692988092273580531, .59983220655588793769, 1)
    r1 <- "(0.180625 - (p - 0.5)*(p - 0.5))"
    rr <- "sqrt(-ln(least(p, 1 - p)))"
    r2 <- sprintf("(%s - 1.6)", rr)
    r3 <- sprintf("(%s - 5)", rr)
    qn <- sprintf(paste("CASE WHEN abs(p - 0.5) <= 0.425 THEN (p - 0.5) * %s / %s",
                        "ELSE (CASE WHEN p < 0.5 THEN -1 ELSE 1 END) * (CASE WHEN %s <= 5 THEN %s / %s ELSE %s / %s END) END"),
                  horner(r1, qa), horner(r1, qb), rr, horner(r2, qc), horner(r2, qd), horner(r3, qe), horner(r3, qf))
    ## dnorm as dnorm.c: for |x| >= 5, exp() of a split square
    x1 <- "(round(abs(x) * 65536) / 65536)"
    dn <- sprintf(paste("CASE WHEN abs(x) < 5 THEN 0.398942280401432677939946059934 * exp(-0.5 * x * x)",
                        "ELSE 0.398942280401432677939946059934 * exp(-0.5 * %s * %s) * exp((-0.5 * (abs(x) - %s) - %s) * (abs(x) - %s)) END"),
                  x1, x1, x1, x1, x1)
    pc <- "CASE WHEN abs(x) > 1 THEN (CASE WHEN x > 0 THEN 1 - atan(1 / x) / pi() ELSE -atan(1 / x) / pi() END) ELSE 0.5 + atan(x) / pi() END"
    ## -expm1(-u) for u >= 0, Kahan's trick: (1 - w) u / -ln(w), w = exp(-u)
    em <- "CASE WHEN u > 40 THEN 1.0 WHEN exp(-u) = 1 THEN u ELSE (1 - exp(-u)) * u / (-ln(exp(-u))) END"
    for (m in list(c("dm_pnorm", "x", pn), c("dm_qnorm", "p", qn), c("dm_dnorm", "x", dn),
                   c("dm_pcauchy", "x", pc), c("dm_mexpm1m", "u", em)))
        dm_execute(D, sprintf("CREATE OR REPLACE TEMP MACRO %s(%s) AS (%s)", m[1L], m[2L], m[3L]))
    D$src$holder$dm_macros <- TRUE
    invisible()
}

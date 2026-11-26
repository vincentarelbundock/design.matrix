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
##   dmc_pmin dmc_pmax dmc_ifelse dmc_rowsums dmc_cbind dmc_col(x, j) dmc_ncol
##   dmc_const(like, value) dmc_lgamma dmc_round (half to even)
##   dmc_is_finite (NaN-safe) dmc_notnan (not NA or NaN)
##   dmc_on(col, x) (same SQL, rows of x)
##   dmc_count(cond, negate) (1 where cond is TRUE, to sum)
##   dmc_freeze(col, table) (computed once: remotely, columns of a new
##                         frame; with table = TRUE, a temporary table of the
##                         rows used with every column of the frame)
##   dmc_pnorm dmc_qnorm dmc_dnorm dmc_pcauchy dmc_qcauchy
##   dmc_logit_inv dmc_logit_mu_eta dmc_binomial_dev(y, mu, w)
##                         (in memory the C code of binomial())
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
## builds a new frame (a view) on top of
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
##                            G x k matrix, as rowsum(col, group) (NA for a
##                            group with a missing value)
##   dmp_y(x, y)              the response of the qr() methods as a column
##   dmp_min_nonzero(x)       a lower bound on the non-zero entries (0 if none)
##   dmc_codes(codes, G)      a grouping as integer codes 1..G: grouped sums
##                            then return all G groups in code order, without
##                            labels, from one-pass C kernels in memory
##   dmp_colsums_groups(x, groups, extra)
##                            several coded groupings at once: per grouping,
##                            list(S = G x p, n = G sums of the scale, extra =
##                            G x k sums of scale * extra); one pass of the
##                            structured kernel in memory
##   dmp_crosstab(a, b, w)    for two grouping columns, the sums of w (counts
##                            without w) for each pair of groups: a
##                            Matrix::dgCMatrix, rows the sorted groups of a,
##                            columns those of b (labels as.character() of
##                            the values, as in memory); with dmc_codes(), the
##                            Ga x Gb matrix in code order
##   dmp_mult(x, B, along)    X B (B a p-vector or p x k matrix) as a column
##                            (in memory a vector for a vector B)
##   dmp_quad(x, V, along)    rowSums((X V) * X), V p x p symmetric positive
##                            semi-definite, as a column
##   dm_colsum(col)           sums over rows: a number per column
##   dmc_sums(...)            sums over rows of several columns (one query)
##   dmc_mean(col)            the mean of a column
##   dmp_colnorms(col)        Euclidean norms of the columns
##   dm_nrows(x)              number of rows of a design or column
##   dm_engine(x)             "sql" for a remote design, ""
##                            in memory (for the description of a fit)
##   dmp_patterns(x, z, w, compress, keys = NULL)
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
##                            `within` is that of the fit on every row. With
##                            `keys` (a column, k columns), patterns are also
##                            split by their values, returned as `keys`
##                            (one row per pattern).
##   dmp_materialize(x)       the design over a copy of its rows: remotely a
##                            temporary table of the rows used with every
##                            column of its frame (later queries scan it
##                            instead of recomputing the frame's views and
##                            joins); a new design (re-make columns with
##                            dmc_on()); in memory x
##   dmp_subset(x, keep)      the design of the rows where the logical column
##                            keep is TRUE (remotely also a DBI::SQL()
##                            predicate); a new remote design
##   dmp_resolve(x, exprs, env, logical, data)
##                            columns for the arguments of a fit (response,
##                            weights, offset, subset): column names, SQL(),
##                            R expressions (in SQL when translatable, else
##                            evaluated in R into a temporary table); returns
##                            list(x, cols, evaluated)
##   dmp_output(x, cols, memory, na, design, attrs)
##                            row-wise results of an algorithm: in memory
##                            memory(d), d the data frame of the named
##                            columns `cols` (one row per row of x; default
##                            d); remotely a temporary view with the columns
##                            of the remote data and `cols`, a string of
##                            class dm_remote_table (padded by napredict(na),
##                            with attributes attrs), or with design = TRUE
##                            the remote design with `cols` added to its
##                            data. One code path computes the columns; only
##                            the container differs.
##
## Remote methods: crossprod, colSums, products and patterns run in the
## database; dmp_groupsum and
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
##   rowSums/rowMeans(D)          table (column .rowsum / .rowmean) done
##   quadratic_diagonal(D, V)     table (column .quad)              done
##   qr(D)                        dmqr_chol on X'X with lm()'s      done
##                                rule, suspects from products
##   qr.coef, qr.qty, qr.R, solve numbers (as in memory)            done
##   qr.fitted, qr.resid, qr.qy,  tables (.fitted, .resid, .qy,     done
##     qr.Q                       .Q1...)
##   as(D, "CsparseMatrix")       collected, same guard as          done
##                                as.matrix (option to raise)
##   D[i, ]                       i a logical remote column or a    done
##                                DBI::SQL() predicate; numeric or
##                                logical vectors refused
##   s * D, D / s, -D             row scaling (number, column       done
##                                name, SQL(), column)
##   D1 + D2, D1 - D2             refused (operator decision)       done (error)
##   finite_rows(D) (internal)    a logical column                  done
##
## Tests: test_primitives.R (the layer), test_parity_ops.R (every operation
## above on the same data, in memory and remote: a table, Parquet).

## ---- remote columns ----

## (`notnull`: columns known never to be NULL, whose sums need no count of
## missing values)
dm_rcols <- function(frame, sql, notnull = FALSE) {
    ## (parenthesized, so that expressions compose with any operator)
    sql <- as.character(sql)
    sql <- sprintf("(%s)", sql)
    structure(list(frame = frame, sql = sql, notnull = rep_len(notnull, length(sql))), class = "dm_rcols")
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
    if (length(x) == 1L && is.na(x)) return("CAST(NULL AS DOUBLE)")
    if (is.logical(x) && length(x) == 1L) return(if (x) "TRUE" else "FALSE")
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

## (in memory a numeric vector: `no` with `yes` where test is TRUE and NA
## where it is NA, as ifelse() but without its copies and attributes)
dmc_ifelse <- function(test, yes, no) {
    f <- dm_rcols_frame(test, yes, no)
    if (is.null(f)) {
        n <- length(test)
        out <- rep_len(as.double(no), n)
        i <- which(test)
        out[i] <- if (length(yes) == 1L) yes else yes[i]
        if (anyNA(test)) out[is.na(test)] <- NA
        return(out)
    }
    dm_rcols(f, sprintf("(CASE WHEN %s THEN %s ELSE %s END)", dm_rcols_operand(test),
                        dm_rcols_operand(yes), dm_rcols_operand(no)))
}

## Columns side by side (cbind).
dmc_cbind <- function(...) {
    args <- list(...)
    f <- do.call(dm_rcols_frame, args)
    if (is.null(f)) return(do.call(cbind, args))
    dm_rcols(f, unlist(lapply(args, dm_rcols_operand)),
             notnull = unlist(lapply(args, function(a) if (dm_is_rcols(a)) unclass(a)$notnull else !is.na(a))))
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
## The logit link of binomial(): inverse, derivative, and the binomial
## deviance residuals (in memory the C code of R's family object, in one
## pass; remotely the same arithmetic in SQL)
dmc_logit_inv <- function(e) {
    if (!dm_is_rcols(e)) return(stats::binomial()$linkinv(e))
    eps <- .Machine$double.eps
    dmc_ifelse(e < -30, eps / (1 + eps), dmc_ifelse(e > 30, (1 / eps) / (1 + 1 / eps), exp(e) / (1 + exp(e))))
}
dmc_logit_mu_eta <- function(e) {
    if (!dm_is_rcols(e)) return(stats::binomial()$mu.eta(e))
    dmc_ifelse(e > 30 | e < -30, .Machine$double.eps, exp(e) / ((1 + exp(e)) * (1 + exp(e))))
}
dmc_binomial_dev <- function(y, mu, w) {
    if (is.null(dm_rcols_frame(y, mu, w))) return(stats::binomial()$dev.resids(y, mu, w))
    xlogx <- function(a, b) dmc_ifelse(a != 0, a * log(a / b), 0)
    2 * w * (xlogx(y, mu) + xlogx(1 - y, 1 - mu))
}

## qcauchy (remotely through tan(): to rounding, where R uses tanpi())
dmc_qcauchy <- function(x) {
    if (!dm_is_rcols(x)) return(stats::qcauchy(x))
    dmc_ifelse(x <= 0.5, -1 / tan(pi * x), 1 / tan(pi * (1 - x)))
}

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
    nn <- unclass(x)$notnull %||% rep(FALSE, k)
    items <- c(sprintf("SUM(CAST(%s AS DOUBLE)) AS s%d", sql, seq_len(k)),
               sprintf("SUM(CASE WHEN (%s) IS NULL THEN 1 ELSE 0 END) AS n%d", sql, seq_len(k))[!nn])
    res <- dm_query(fr, sprintf("SELECT %s FROM (%s) AS s WHERE %s", paste(items, collapse = ", "),
                                fr$src$select, fr$filter))
    vapply(seq_len(k), function(j) {
        v <- as.double(res[[paste0("s", j)]])
        if (!nn[j] && as.numeric(res[[paste0("n", j)]]) > 0) return(NA_real_)
        if (is.na(v)) 0 else v
    }, 0)
}

## Sums over rows of several columns (vectors or numbers in memory; one
## query remotely): a number per argument.
dmc_sums <- function(...) {
    a <- list(...)
    if (any(vapply(a, dm_is_rcols, NA))) return(dm_colsum(do.call(dmc_cbind, a)))
    vapply(a, function(v) sum(v), 0)
}

## The mean of a column (in memory as mean() computes it).
dmc_mean <- function(x) if (dm_is_rcols(x)) dm_colsum(x) / dm_nrows(x) else mean(x)

## "sql" for a remote design, "" in memory.
dm_engine <- function(x) if (inherits(x, "dm_remote")) "sql" else ""

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
    if (is.null(y)) return(crossprod(x))
    if (inherits(x, "dm_remote") && !inherits(y, "dm_remote") && !dm_is_rcols(y))
        stop("`y` must be a design or a column.", call. = FALSE)
    crossprod(x, y)
}

dmp_colsums <- function(x, group = NULL) {
    if (is.null(group)) colSums(x) else colSums(x, group = group)
}

## SQL GROUP BY over the frame of columns: list(groups (a data frame of the
## group values, sorted as rowsum() sorts them), sums (matrix)).
dm_group_query <- function(groups, values, nulls = FALSE) {
    fr <- do.call(dm_rcols_frame, c(groups, values))
    gs <- vapply(groups, function(g) unclass(g)$sql, "")
    vs <- unlist(lapply(values, function(v) unclass(v)$sql))
    items <- c(sprintf("%s AS g%d", gs, seq_along(gs)),
               if (length(vs)) sprintf("SUM(CAST(%s AS DOUBLE)) AS v%d", vs, seq_along(vs)) else "COUNT(*) AS v1",
               if (nulls && length(vs)) sprintf("SUM(CASE WHEN (%s) IS NULL THEN 1 ELSE 0 END) AS n%d", vs, seq_along(vs)))
    res <- dm_query(fr, sprintf("SELECT %s FROM (%s) AS s WHERE %s GROUP BY ALL",
                                paste(items, collapse = ", "), fr$src$select, fr$filter))
    res
}

dmp_groupsum <- function(col, group) {
    if (dm_is_codes(group)) return(dmp_groupsum_codes(col, group))
    if (!dm_is_rcols(col) && !dm_is_rcols(group)) {
        return(rowsum(if (is.matrix(col)) col else matrix(col, ncol = 1L), group, reorder = TRUE))
    }
    res <- dm_group_query(list(group), list(col), nulls = TRUE)
    ug <- sort(unique(res$g1), na.last = TRUE, method = "quick")
    if (anyNA(ug)) warning("missing values for 'group'", call. = FALSE)
    i <- match(res$g1, ug)
    k <- length(col)
    out <- matrix(0, length(ug), k, dimnames = list(as.character(ug), NULL))
    for (j in seq_len(k)) {
        v <- as.double(res[[paste0("v", j)]])
        ## a group with a missing value sums to NA, as rowsum() gives
        v[as.numeric(res[[paste0("n", j)]]) > 0] <- NA_real_
        out[i, j] <- v
    }
    out
}

dmp_crosstab <- function(a, b, w = NULL) {
    if (dm_is_codes(a) && dm_is_codes(b)) return(dmp_crosstab_codes(a, b, w))
    if (!dm_is_rcols(a) && !dm_is_rcols(b) && !dm_is_rcols(w)) {
        ua <- sort(unique(a), na.last = TRUE, method = "quick")
        ub <- sort(unique(b), na.last = TRUE, method = "quick")
        return(Matrix::sparseMatrix(i = match(a, ua), j = match(b, ub),
                                    x = if (is.null(w)) rep(1, length(a)) else as.double(w),
                                    dims = c(length(ua), length(ub)),
                                    dimnames = list(as.character(ua), as.character(ub))))
    }
    res <- dm_group_query(list(a, b), if (is.null(w)) list() else list(w))
    ## the group values as in memory: sorted in R, labels as.character()
    ua <- sort(unique(res$g1), na.last = TRUE, method = "quick")
    ub <- sort(unique(res$g2), na.last = TRUE, method = "quick")
    Matrix::sparseMatrix(i = match(res$g1, ua), j = match(res$g2, ub), x = as.double(res$v1),
                         dims = c(length(ua), length(ub)),
                         dimnames = list(as.character(ua), as.character(ub)))
}

dmp_mult <- function(x, B, along = list()) {
    if (!inherits(x, "dm_remote")) return(if (is.matrix(B)) x %*% B else drop(x %*% B))
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
    if (!inherits(x, "dm_remote")) {
        J <- as_jac(x)
        storage.mode(V) <- "double"
        if (!jac_kernel_ok(J)) {
            X <- as.matrix(J)
            return(unname(rowSums((X %*% V) * X)))
        }
        g <- dm_groups(J, cost = dm_row_cost(J)^2)
        q <- do.call(.Call, c(list(C_dm_jac_quad), jac_args(J),
                              list(V, if (is.null(g)) NULL else g$first, design_threads())))
        return(if (is.null(g)) q else q[g$map])
    }
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
dmp_patterns <- function(x, z, w = NULL, compress = "auto", keys = NULL) {
    if (inherits(x, "dm_remote")) {
        mode <- if (isTRUE(compress)) TRUE else if (isFALSE(compress)) FALSE else dm_compress_mode()
        ## (extra key columns: more keys of the GROUP BY)
        nm <- if (!is.null(keys)) paste0("__dm_g", seq_along(keys))
        if (!is.null(keys)) x <- dm_remote_cols(x, keys, "keys")$D
        pt <- dm_remote_patterns(x, unclass(z)$sql, if (!is.null(w)) unclass(w)$sql, within = TRUE,
                                 mode = mode, group = if (!is.null(keys)) unclass(keys)$sql)
        if (is.null(pt)) return(NULL)
        pd <- dm_pattern_design(x, pt$P[dm_compress_keys(x)])
        k <- pd$keep
        sw <- pt$sw[k]
        pos <- which(sw > 0)
        ws <- sw[pos]
        if (is.null(w) && all(ws == 1)) ws <- NULL
        return(list(x = pd$design[pos, ], z = pt$sy[k, 1L][pos] / sw[pos], w = ws,
                    within = sum(pt$within[k]), rows = sum(pt$n[k]),
                    keys = if (!is.null(keys)) as.matrix(pt$P[k[pos], nm, drop = FALSE])))
    }
    g <- dmfit_groups(x, compress)
    if (is.null(g)) return(NULL)
    if (!is.null(keys)) {
        ## patterns within the values of the key columns
        K <- if (is.matrix(keys)) keys else matrix(keys, ncol = 1L)
        key <- do.call(paste, c(list(g$map), lapply(seq_len(ncol(K)), function(j) sprintf("%a", K[, j])), sep = "\r"))
        map <- match(key, unique(key))
        g <- list(map = map, first = match(seq_len(max(map)), map))
        if (length(g$first) == nrow(x) || (!isTRUE(compress) && length(g$first) > 0.5 * nrow(x))) return(NULL)
    }
    ss <- dmfit_suffstats(z, w, g$map, length(g$first))
    ## sums that overflow (w * z beyond the range of doubles): no compression
    if (!all(is.finite(ss$w), is.finite(ss$mean), is.finite(ss$within))) return(NULL)
    keep <- ss$w > 0
    ws <- ss$w[keep]
    if (is.null(w) && all(ws == 1)) ws <- NULL
    list(x = x[g$first[keep], ], z = ss$mean[keep], w = ws, within = ss$within, rows = nrow(x),
         keys = if (!is.null(keys)) (if (is.matrix(keys)) keys else matrix(keys, ncol = 1L))[g$first[keep], , drop = FALSE])
}

## Row-wise results: in memory memory(d), d the data frame of the columns
## (default: d); remotely a table (a temporary view with the columns of the
## remote data and `cols`), padded by napredict(na) for a fit with
## na.exclude, with attributes `attrs`; with design = TRUE, the remote
## design with the columns added to its data (so that they can scale its
## rows).
dmp_output <- function(x, cols, memory = identity, na = NULL, design = FALSE, attrs = list()) {
    if (!inherits(x, "dm_remote")) return(memory(as.data.frame(lapply(cols, function(v) as.vector(v)))))
    taken <- intersect(names(cols), names(x$src$types))
    if (length(taken)) stop(sprintf("The remote data already has a column `%s`.", taken[1L]), call. = FALSE)
    if (design) return(dm_with_columns(x, do.call(dmc_cbind, unname(cols)), names(cols)))
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
    out <- stats::napredict(na, structure(name, class = "dm_remote_table", connection = dm_con(fr),
                                          holder = fr$src$holder))
    for (a in names(attrs)) attr(out, a) <- attrs[[a]]
    out
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


## ---- coded groups, multi-group sums (the kernels of the fixed effects) ----

## A grouping given as integer codes 1..G (no missing values): `codes` an
## integer vector in memory, a remote column of integer codes remotely.
## Grouped sums with coded groups return all G groups in code order
## (zeros for empty groups), without labels, and use one-pass C kernels in
## memory (src/fe.c).
dmc_codes <- function(codes, G) {
    if (!dm_is_rcols(codes)) codes <- as.integer(codes)
    structure(list(codes = codes, G = as.integer(G)), class = "dm_codes")
}

dm_is_codes <- function(x) inherits(x, "dm_codes")

## sums of columns within coded groups: G x k
dmp_groupsum_codes <- function(col, group) {
    G <- group$G
    if (!dm_is_rcols(col) && !dm_is_rcols(group$codes)) {
        W <- if (is.matrix(col)) col else matrix(col, ncol = 1L)
        storage.mode(W) <- "double"
        return(.Call(C_dm_fe_gsum, group$codes, G, W))
    }
    res <- dm_group_query(list(group$codes), list(col), nulls = TRUE)
    k <- length(col)
    out <- matrix(0, G, k)
    i <- as.integer(res$g1)
    for (j in seq_len(k)) {
        v <- as.double(res[[paste0("v", j)]])
        v[as.numeric(res[[paste0("n", j)]]) > 0] <- NA_real_
        out[i, j] <- v
    }
    out
}

## weighted counts by pair of coded groups: Ga x Gb dgCMatrix
dmp_crosstab_codes <- function(a, b, w = NULL) {
    if (!dm_is_rcols(a$codes) && !dm_is_rcols(b$codes) && !dm_is_rcols(w)) {
        ## (the kernel returns the transpose of the crosstab of its first
        ## two arguments, and is fastest when the first is sorted)
        ct <- .Call(C_dm_fe_crosstab, b$codes, a$codes, if (is.null(w)) NULL else as.double(w), b$G, a$G)
        return(methods::new("dgCMatrix", i = ct[[2L]], p = ct[[1L]], x = ct[[3L]], Dim = c(a$G, b$G)))
    }
    res <- dm_group_query(list(a$codes, b$codes), if (is.null(w)) list() else list(w))
    Matrix::sparseMatrix(i = as.integer(res$g1), j = as.integer(res$g2), x = as.double(res$v1),
                         dims = c(a$G, b$G), repr = "C")
}

## Column sums of a design within several coded groups at once:
## dmp_colsums_groups(x, groups, extra = NULL), x possibly scaled (s * x),
## groups a list of dmc_codes(), extra a column (k columns) or NULL. Returns
## one list(S = G x p sums of the scaled design, n = G sums of the scale (the
## counts without one), extra = G x k sums of scale * extra) per grouping.
## In memory: one pass of the structured kernel for all groupings when it
## applies (dm_fe_struct_gsum), else the row kernels; remotely, one GROUP BY
## per grouping.
dmp_colsums_groups <- function(x, groups, extra = NULL) {
    if (inherits(x, "dm_remote")) {
        s <- x$scale
        one <- dm_rcols(dm_remote_unscaled(x), if (is.null(s)) "CAST(1 AS DOUBLE)" else s)
        return(lapply(groups, function(g) {
            gc <- g$codes
            S <- matrix(0, g$G, ncol(x))
            cs <- dmp_colsums(x, group = gc)
            S[as.integer(rownames(cs)), ] <- unname(cs)
            ns <- dmp_groupsum_codes(if (is.null(extra)) one else dmc_cbind(one, one * extra), g)
            list(S = S, n = ns[, 1L], extra = if (!is.null(extra)) ns[, -1L, drop = FALSE])
        }))
    }
    ## in memory: the plain design and its row scale
    w <- NULL
    D <- x
    if (dm_is_comb(x)) {
        if (!is.null(x$design2)) stop("A combination of two designs has no grouped sums here.", call. = FALSE)
        D <- x$design
        w <- if (is.null(x$scale)) NULL else rep_len(as.double(x$scale), nrow(D))
    }
    G <- vapply(groups, function(g) g$G, 1L)
    codes <- lapply(groups, function(g) g$codes)
    Y <- if (is.null(extra)) NULL else { E <- if (is.matrix(extra)) extra else matrix(extra, ncol = 1L); storage.mode(E) <- "double"; E }
    if (dm_colsums_struct_ok(D, G)) {
        M <- .Call(C_dm_fe_struct_gsum, D$setup, D$vars, w, codes, as.integer(G), Y, TRUE)
        info <- dm_struct_info(D)
        icpt <- isTRUE(D$setup$intercept > 0L)
        p <- ncol(D)
        return(lapply(seq_along(groups), function(k) {
            m <- M[[k]]
            S <- matrix(0, G[k], p)
            q <- 1L
            for (t in seq_along(info)) {
                Tt <- info[[t]]
                if (t == 1L && icpt) {
                    S[, Tt$cols[1L, 1L]] <- m[1L, ]
                    next
                }
                L <- nrow(Tt$C)
                ns <- nrow(Tt$cols)
                for (u in seq_len(ns)) {
                    rows <- q + (u - 1L) * L + seq_len(L)
                    S[, Tt$cols[u, ]] <- crossprod(m[rows, , drop = FALSE], Tt$C)
                }
                q <- q + L * ns
            }
            list(S = S, n = m[1L, ], extra = if (is.null(Y)) NULL else t(m[q + seq_len(ncol(Y)), , drop = FALSE]))
        }))
    }
    one <- if (is.null(w)) rep(1, nrow(D)) else w
    lapply(seq_along(groups), function(k) {
        J <- as_jac(D)
        W <- matrix(one, ncol = 1L)
        S <- if (jac_kernel_ok(J)) t(dm_crossprod_kernel(J, W, codes[[k]], G[k]))
             else .Call(C_dm_fe_gsum, codes[[k]], G[k], one * as.matrix(D))
        ns <- .Call(C_dm_fe_gsum, codes[[k]], G[k], cbind(one, if (!is.null(Y)) one * Y))
        list(S = S, n = ns[, 1L], extra = if (!is.null(Y)) ns[, -1L, drop = FALSE])
    })
}

dm_colsums_struct_ok <- function(D, G) {
    if (!isTRUE(getOption("design.matrix.structured", TRUE))) return(FALSE)
    if (!dm_struct_ok(D)) return(FALSE)
    info <- dm_struct_info(D)
    Q <- 2 + sum(vapply(info, function(Tt) nrow(Tt$C) * nrow(Tt$cols), 0))
    all(Q * as.double(G) <= 2^25)
}

## ---- more column functions (requested by the GLM work) ----

## the number of columns of a column
dmc_ncol <- function(x) if (dm_is_rcols(x)) length(x) else NCOL(x)

## one column of a multi-column column
dmc_col <- function(x, j) {
    if (!dm_is_rcols(x)) return(if (is.matrix(x)) x[, j] else x)
    dm_rcols(unclass(x)$frame, unclass(x)$sql[j])
}

## a constant column over the rows of a column's frame (or a design)
dmc_const <- function(like, value) {
    if (inherits(like, "dm_remote")) return(dm_rcols(like, dm_sql_double(value)))
    if (!dm_is_rcols(like)) return(rep(as.double(value), NROW(like)))
    dm_rcols(unclass(like)$frame, dm_sql_double(value))
}

dmc_lgamma <- function(x) {
    if (!dm_is_rcols(x)) return(lgamma(x))
    dm_rcols(unclass(x)$frame, sprintf("lgamma(%s)", unclass(x)$sql))
}

## round half to even, as R's round(x)
dmc_round <- function(x) {
    if (!dm_is_rcols(x)) return(round(x))
    dm_rcols(unclass(x)$frame, sprintf("round_even(%s, 0)", unclass(x)$sql))
}

## is.finite, NaN-safe (in DuckDB, NaN = NaN and NaN > 0 are true)
dmc_is_finite <- function(x) {
    if (!dm_is_rcols(x)) return(is.finite(x))
    dm_rcols(unclass(x)$frame, sprintf("COALESCE(isfinite(CAST(%s AS DOUBLE)), FALSE)", unclass(x)$sql))
}

## TRUE where x is not missing (NA or NaN; infinite values count as present).
dmc_notnan <- function(x) dmc_is_finite(x) | (abs(x) == Inf)

## The rows where a logical column is TRUE (with negate, where it is not
## TRUE: FALSE or missing), as a column to sum: 1 or 0. In memory a logical
## vector, which sums as 0/1.
dmc_count <- function(cond, negate = FALSE) {
    if (!dm_is_rcols(cond)) {
        if (anyNA(cond)) cond <- !is.na(cond) & cond
        return(if (negate) !cond else cond)
    }
    sql <- unclass(cond)$sql
    dm_rcols(unclass(cond)$frame, sprintf(if (negate) "CASE WHEN %s THEN 0 ELSE 1 END" else "CASE WHEN %s THEN 1 ELSE 0 END", sql),
             notnull = TRUE)
}

## A column computed once: remotely, its expressions become columns of a
## new frame (a projection over the frame of the column, one level deeper,
## as products are), so that later expressions refer to them instead of
## repeating them (DuckDB expands macros such as dm_qnorm() inline, and
## nested expressions would grow without bound). Columns frozen separately
## from one frame live on sibling frames, which cannot be combined: freeze
## them together (dmc_cbind()). In memory, the column itself. With table =
## TRUE, the new frame is a temporary table of the rows used (computed once
## for all later queries, which then scan it: for columns that many queries
## read, such as the working weights and responses of a GLM with fixed
## effects); its rows keep the filter of the design.
dmc_freeze <- function(col, table = FALSE) {
    if (!dm_is_rcols(col)) return(col)
    fr <- unclass(col)$frame
    sql <- unclass(col)$sql
    nm <- paste0("__dm_f", substr(dm_tmp_name("f"), 3L, 10L), "_", seq_along(sql))
    x <- dm_remote_extend(dm_remote_unscaled(fr), stats::setNames(sql, nm), "column")
    if (table) x <- dm_remote_to_table(x)
    x$depth <- dm_frame_depth(fr) + 1L
    dm_rcols(x, sprintf("s.%s", dm_qi(x, nm)))
}

## The frame of a remote design as a temporary table of its rows used (with
## every column of the frame; the design's filter stays valid on them).
dm_remote_to_table <- function(x) {
    tn <- dm_qi(x, dm_tmp_name("dm_mat"))
    dm_execute(x, sprintf("CREATE TEMPORARY TABLE %s AS SELECT * FROM (%s) AS s WHERE %s", tn, x$src$select, x$filter))
    x$src$select <- sprintf("SELECT * FROM %s", tn)
    x$src$base_table <- tn
    x
}

dmp_materialize <- function(x) {
    if (!inherits(x, "dm_remote")) return(x)
    sc <- x$scale
    x <- dm_remote_to_table(dm_remote_unscaled(x))
    x$filter <- "TRUE"
    x$depth <- 0L
    x$cache <- new.env(parent = emptyenv())
    x$id <- dm_tmp_name("dm")
    x$scale <- sc
    x
}

## ---- inputs of the methods of a design ----

## Which rows of a design matrix are finite: rowSums(!is.finite(as.matrix(x)))
## == 0, computed without forming the dense matrix (TRUE for each row whose
## entries are all finite, FALSE for a row with any NA, NaN or infinite entry,
## including entries that overflow, such as the product of two large numeric
## variables, and, for a combination of designs, scaled sums that overflow).
finite_rows <- function(x) {
    ## a remote design: a logical column (every entry of the row finite)
    if (inherits(x, "dm_remote")) {
        P <- dmp_mult(x, diag(1, ncol(x)))
        return(dm_rcols(unclass(P)$frame,
                        paste(sprintf("COALESCE(isfinite(%s), FALSE)", unclass(P)$sql), collapse = " AND ")))
    }
    checkmate::assert_class(x, "design_matrix")
    J <- as_jac(x)
    if (!jac_kernel_ok(J)) {
        return(unname(rowSums(!is.finite(as.matrix(x))) == 0))
    }
    do.call(.Call, c(list(C_dm_jac_finite), jac_args(J), list(design_threads())))
}

## A lower bound on the absolute value of the non-zero entries of x (0 when
## no bound is available: always for a remote design).
dmp_min_nonzero <- function(x) {
    if (inherits(x, "dm_remote")) return(0)
    terms <- dm_terms(x)
    out <- Inf
    for (t in terms) {
        d <- t$design
        b <- if (kernel_ok(d)) .Call(C_dm_min_nonzero, d$setup, d$vars) else {
            X <- as.matrix(d)
            a <- abs(X[X != 0 & is.finite(X)])
            if (length(a)) min(a) else 1
        }
        s <- t$scale
        if (!is.null(s)) {
            a <- abs(s[s != 0 & is.finite(s)])
            b <- b * (if (length(a)) min(a) else 1)
        }
        out <- min(out, b)
    }
    ## a sum of two designs can cancel to small values
    if (length(terms) > 1L) 0 else out
}

## The `y` of the qr() methods as a column of design x: in memory a finite
## double matrix (from a vector or matrix with nrow(x) rows); remotely a
## remote column from column names, a DBI::SQL() expression or a remote
## column. Returns list(Y, vec (was y one column?), rows, cols (names)).
dmp_y <- function(x, y) {
    if (!inherits(x, "dm_remote")) return(dmqr_y(y, nrow(x)))
    D <- dm_remote_unscaled(x)
    if (dm_is_rcols(y)) return(list(Y = y, vec = length(y) == 1L, cols = NULL))
    if (inherits(y, "SQL")) return(list(Y = dm_rcols(D, as.character(y)), vec = TRUE, cols = NULL))
    if (is.character(y)) {
        dm_remote_column(D, y, "y", several = TRUE)
        return(list(Y = dm_rcols(D, sprintf("CAST(s.%s AS DOUBLE)", dm_qi(D, y))), vec = length(y) == 1L,
                    cols = if (length(y) > 1L) y))
    }
    stop("For a remote design, `y` must name columns of its data, or be a DBI::SQL() expression.", call. = FALSE)
}

## ---- rows of a design ----

## The design of the rows where `keep` (a logical column, NA counting as
## FALSE) is TRUE: x[keep, ] in memory; remotely, the same design with its
## filter restricted (a new design: columns of x must be re-made on it,
## for instance with dmc_on()).
dmp_subset <- function(x, keep) {
    if (!inherits(x, "dm_remote")) {
        if (dm_is_rcols(keep)) stop("A remote column cannot select rows of a design in memory.", call. = FALSE)
        return(x[!is.na(keep) & keep, ])
    }
    if (inherits(keep, "SQL")) keep <- dm_rcols(x, as.character(keep))
    if (!dm_is_rcols(keep))
        stop("Rows of a remote design are selected by a logical remote column or a DBI::SQL() predicate.",
             call. = FALSE)
    fr <- unclass(keep)$frame
    if (dm_frame_depth(x) > dm_frame_depth(fr)) fr <- x
    out <- dm_rebase(x, fr)
    out$filter <- sprintf("(%s) AND COALESCE(%s, FALSE)", out$filter, unclass(keep)$sql)
    out$n <- as.numeric(dm_query(out, sprintf("SELECT COUNT(*) AS n FROM (%s) AS s WHERE %s",
                                              out$src$select, out$filter))$n)
    out$cache <- new.env(parent = emptyenv())
    out$id <- dm_tmp_name("dm")
    out
}

## The same column (its SQL) over the rows of design x (remote; in memory,
## the column itself). Valid when its expressions read columns of the
## source of x.
dmc_on <- function(col, x) {
    if (!dm_is_rcols(col)) return(col)
    dm_rcols(x, unclass(col)$sql)
}

## ---- arguments of a fit (round 7 of 2hba) ----

## Columns for the arguments of a fit over the rows of a remote design:
## exprs, a named list of unevaluated arguments (each a column name, a
## DBI::SQL() expression, or an R expression over the columns, evaluated in
## env); logical: the names of logical arguments (a subset). R expressions
## with a SQL translation run in the database; others are evaluated in R,
## batch by batch, into a temporary table (dm_materialize()), which then is
## the source of the design. Returns list(x = the design (possibly on that
## new source), cols = named list of columns (NULL for NULL arguments),
## evaluated = names evaluated in R, with their reasons).
## In memory, the arguments are evaluated in `data` (a data frame, the rows
## of x) and returned as vectors.
dmp_resolve <- function(x, exprs, env, logical = character(), data = NULL) {
    if (!inherits(x, "dm_remote")) {
        cols <- lapply(exprs, function(e) {
            if (is.null(e)) return(NULL)
            if (is.character(e) && length(e) == 1L) return(data[[e]])
            eval(e, data, env)
        })
        return(list(x = x, cols = cols, evaluated = list()))
    }
    con <- dm_con(x)
    items <- list()
    later <- list()
    sqls <- list()
    for (nm in names(exprs)) {
        e <- exprs[[nm]]
        if (is.null(e)) next
        a <- dm_resolve_or_eval(e, env, x$src$types, con, nm, bool = nm %in% logical)
        if (!is.null(a$expr)) {
            items[[nm]] <- a
            later[[paste0("__dm_eval_", nm)]] <- a$expr
        } else {
            sqls[[nm]] <- a$sql
        }
    }
    if (length(later)) {
        x$src <- dm_materialize(x$src, later, env, logical = paste0("__dm_eval_", logical), chunk_size = x$chunk_size)
        for (nm in names(items)) {
            col <- DBI::dbQuoteIdentifier(con, paste0("__dm_eval_", nm))
            sqls[[nm]] <- if (nm %in% logical) sprintf("s.%s", col) else sprintf("CAST(s.%s AS DOUBLE)", col)
        }
    }
    cols <- lapply(names(exprs), function(nm) if (is.null(sqls[[nm]])) NULL else dm_rcols(x, sqls[[nm]]))
    names(cols) <- names(exprs)
    list(x = x, cols = cols, evaluated = items)
}

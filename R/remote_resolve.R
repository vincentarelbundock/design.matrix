## Arguments of a fit as columns of a remote design: a formula and remote
## data (see remote.R) with a response, weights, an offset and a subset,
## turned into the remote design of the rows used and columns over them (in
## SQL, or evaluated in R, batch by batch, into a temporary table):
## design_remote_frame(), the remote counterpart of design_frame(), on which the
## input stage of the fits is built (design_problem() of xtx.fit). Also
## design_resolve() (arguments as columns of an existing design) and
## is_remote_source().

## The arguments of remote data that a fit on a formula takes in `...`.
xtx_remote_keys <- c("table", "sql", "engine", "chunk_size")

## The remote design of a fit and the columns of its arguments (one code
## path for lm_xtx(), glm_xtx() and fixed effects).
##   formula: the model formula; the design is that of its right-hand side
##     (offset() terms included; their values are columns too), and its
##     response is given by the caller in `args`.
##   src: an opened source (design_source()).
##   args: named list of unevaluated arguments (NULL: absent), resolved in
##     env: column names, DBI::SQL() or R expressions over the columns.
##     Numeric, except those named in `logical` (`subset`, which restricts
##     the rows).
##   present: columns that must not be missing (FE and cluster variables, a
##     factor response).
##   na: "omit" (or "exclude": the same rows), "fail" or "pass" (lm()'s
##     errors when a row has a missing value).
##   mresp: the error for a response `y` evaluated in R that gives a
##     matrix with several columns (NULL: no check).
## R expressions with no SQL translation are evaluated in R into a
## temporary table (`engine = "sql"` refuses). The
## rows used are those of the subset where no model variable, argument or
## `present` column is missing, and factor levels are those of these rows,
## as model.frame(drop.unused.levels = TRUE) gives them in memory.
## Returns list(x = the design, cols = named list of columns (NULL for
## absent arguments), off = the offset (offset() terms plus the `offset`
## argument) or NULL, offk = the largest number of columns of an offset
## evaluated in R (a matrix: the caller reports its length, as in memory),
## wk = the number of columns of the weights evaluated in R (likewise)).
#' @rdname design_source
#' @export
design_remote_frame <- function(formula, src, args, env, dots = list(), contrasts = NULL, logical = character(),
                          present = character(), na = "omit", mresp = NULL) {
    miss <- setdiff(present, names(src$types))
    if (length(miss)) stop(sprintf("No column `%s` in the remote data.", miss[1L]), call. = FALSE)
    tt <- stats::terms(formula, data = src$proto)
    tv <- xtx_terms_vars(tt)
    rhs <- if (attr(tt, "response") > 0L) stats::delete.response(tt) else tt
    attr(rhs, "predvars") <- NULL
    rhs <- stats::formula(rhs)
    environment(rhs) <- environment(formula)
    ## offset() terms are arguments too: their SQL, or the column computed
    ## in R, replaces the variable of the design (`override`)
    offt <- list()
    for (k in attr(tt, "offset")) offt[[tv$names[[k]]]] <- tv$exprs[[k]][[2L]]
    all <- c(Filter(Negate(is.null), args), offt)
    label <- function(nm) if (nm %in% names(offt)) "offset term" else if (nm %in% c("y", "a", "b")) "response" else nm
    con <- src$holder$con
    q <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
    res <- list()
    items <- list()
    exprs <- list()
    mat <- list()
    for (nm in names(all)) {
        a <- xtx_resolve_or_eval(all[[nm]], env, src$types, con, label(nm), bool = nm %in% logical)
        if (is.null(a$expr)) {
            res[[nm]] <- a
            next
        }
        col <- paste0("__xtx_eval", length(exprs) + 1L)
        items[[label(nm)]] <- a
        exprs[[col]] <- a$expr
        if (nm %in% c("offset", "weights") || nm %in% names(offt)) mat[[col]] <- "count"
        if (nm == "y" && !is.null(mresp)) mat[[col]] <- mresp
        res[[nm]] <- list(sql = if (nm %in% logical) sprintf("s.%s", q(col)) else sprintf("CAST(s.%s AS DOUBLE)", q(col)),
                          inputs = col, logical = nm %in% logical)
    }
    if (length(exprs)) {
        dots <- xtx_untranslated_engine(items, dots)
        lg <- unlist(lapply(res, function(r) if (isTRUE(r$logical)) r$inputs))
        src <- xtx_materialize(src, exprs, env, logical = lg,
                              chunk_size = dots$chunk_size %||% getOption("xtx.chunk_size", 250000),
                              matrix = mat)
    }
    ## (the number of columns of the offsets and weights evaluated in R)
    nc <- attr(src, "ncol")
    ncols <- function(nms) max(c(1L, unlist(nc[unlist(lapply(res[intersect(nms, names(res))], `[[`, "inputs"))])))
    offk <- ncols(c("offset", names(offt)))
    wk <- ncols("weights")
    dots$table <- dots$sql <- NULL
    nonlog <- setdiff(names(res), c(logical, names(offt)))
    extra <- stats::setNames(lapply(res[nonlog], function(r) list(sql = r$sql, inputs = r$inputs, kind = "numeric")),
                             sprintf("(%s)", nonlog))
    for (col in present) extra[[paste0("(present:", col, ")")]] <- list(inputs = col, kind = "present")
    override <- lapply(res[names(offt)], function(r) list(sql = r$sql, inputs = r$inputs))
    if (offk > 1L && length(offt)) {
        ## (an offset() term gave a matrix, which the caller reports with
        ## the number of rows used: the offset() terms are then columns
        ## that only select rows, and the design that of the other terms)
        for (nm in names(offt)) extra[[sprintf("(%s)", nm)]] <- c(override[[nm]], kind = "numeric")
        override <- list()
        tr <- stats::terms(rhs, data = src$proto)
        tl <- attr(tr, "term.labels")
        rhs <- stats::reformulate(if (length(tl)) tl else "1", intercept = attr(tr, "intercept") > 0L,
                                  env = environment(rhs))
    }
    sub <- if (any(logical %in% names(res))) paste(sprintf("(%s)", vapply(res[intersect(logical, names(res))], `[[`, "", "sql")), collapse = " AND ")
    ## (no rows left: the error of lm())
    D <- tryCatch(xtx_remote(rhs, src, contrasts, NULL, dots, drop_unused = TRUE, subset = sub, extra = extra,
                            override = override),
                  xtx_no_rows = function(e) stop("0 (non-NA) cases", call. = FALSE))
    if (na %in% c("fail", "pass")) xtx_remote_fit_na(D, res, na, nonlog, names(offt))
    cols <- lapply(names(args), function(nm) if (!is.null(res[[nm]])) xtx_column(D, res[[nm]]$sql))
    names(cols) <- names(args)
    off <- c(res[["offset"]]$sql, vapply(res[names(offt)], `[[`, "", "sql"))
    list(x = D, cols = cols, off = if (length(off)) xtx_column(D, paste(sprintf("(%s)", off), collapse = " + ")),
         offk = offk, wk = wk)
}

## na.fail and na.pass: the errors of lm() when rows have missing values.
xtx_remote_fit_na <- function(D, res, na, args, offt) {
    total <- as.numeric(xtx_query(D, sprintf("SELECT COUNT(*) AS n FROM (%s) AS s", D$src$select))$n)
    if (total == D$n) return(invisible())
    if (na == "fail") stop("missing values in object", call. = FALSE)
    bad <- function(nm) {
        r <- res[[nm]]
        if (is.null(r)) return(FALSE)
        as.numeric(xtx_query(D, sprintf("SELECT COUNT(*) AS n FROM (%s) AS s WHERE (%s) IS NULL OR isnan(%s)",
                                       D$src$select, r$sql, r$sql))$n) > 0
    }
    if (any(vapply(c("y", "a", "b", "offset", offt), bad, NA))) stop("NA/NaN/Inf in 'y'", call. = FALSE)
    if (bad("weights")) stop("missing or negative weights not allowed", call. = FALSE)
    stop("NA/NaN/Inf in 'x'", call. = FALSE)
}

## R expressions over the remote data with no SQL translation (a response,
## an offset() term, the `weights`, `offset` or `subset` of a fit): they are
## evaluated in R, batch by batch, and written with all the columns of the
## source to a temporary table, whose new columns the fit then reads like
## any other. One pass over the data; the rows, levels, missing values and
## compressed sums are then exact, as for translated expressions.
##   exprs: named list of R expressions (the names of the new columns);
##   logical: names of the expressions to keep as logical (subset);
##   model: for the expressions that are model variables (named as in
##     exprs), list(inputs, factor, label), see xtx_eval_model(). Their
##     columns and kinds are returned in the attribute "model" of the
##     result: list(kind, cols) per expression.
##   matrix: for expressions that may give a matrix with one row per row
##     (named as in exprs): "count" (an offset or weights: its first column
##     is written, and its number of columns returned in the attribute
##     "ncol" of the result, for the caller's error), or the message of the
##     error for a matrix with several columns. A one-column matrix is its
##     column.
## Returns the new source (design_source() of the temporary table),
## which keeps the original one alive.
xtx_materialize <- function(src, exprs, env, logical = character(), chunk_size = 250000, model = list(),
                           matrix = list()) {
    info <- list()
    ncols <- list()
    con <- src$holder$con
    q <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
    from <- src$base_table
    raw <- NULL
    if (is.null(from)) {
        raw <- q(xtx_tmp_name("xtx_raw"))
        DBI::dbExecute(con, sprintf("CREATE TEMPORARY TABLE %s AS %s", raw, src$select))
        on.exit(DBI::dbExecute(con, sprintf("DROP TABLE IF EXISTS %s", raw)), add = TRUE)
        from <- raw
    }
    r <- DBI::dbGetQuery(con, sprintf("SELECT MIN(rowid) AS lo, MAX(rowid) AS hi FROM %s", from))
    out <- xtx_tmp_name("xtx_eval")
    created <- FALSE
    write <- function(df) {
        if (!created) {
            DBI::dbWriteTable(con, out, df, temporary = TRUE)
            created <<- TRUE
        } else {
            DBI::dbAppendTable(con, out, df)
        }
    }
    evaluate <- function(chunk) {
        for (nm in names(exprs)) {
            if (!is.null(model[[nm]])) {
                r <- xtx_eval_model(exprs[[nm]], chunk, env, model[[nm]], nm, info[[nm]])
                info[[nm]] <<- r$info
                for (cn in names(r$cols)) chunk[[cn]] <- r$cols[[cn]]
                next
            }
            v <- eval(exprs[[nm]], chunk, env)
            if (length(v) == 1L && nrow(chunk) != 1L) v <- rep(v, nrow(chunk))
            if (!is.null(matrix[[nm]]) && length(dim(v)) == 2L && nrow(v) == nrow(chunk)) {
                if (ncol(v) > 1L) {
                    if (!identical(matrix[[nm]], "count")) stop(matrix[[nm]], call. = FALSE)
                    ncols[[nm]] <<- max(ncols[[nm]] %||% 1L, ncol(v))
                }
                v <- v[, 1L]
            }
            if (length(v) != nrow(chunk) || !is.null(dim(v)))
                stop(sprintf("`%s` gives %d values for %d rows.", deparse2(exprs[[nm]]), length(v), nrow(chunk)),
                     call. = FALSE)
            if (nm %in% logical) {
                if (!is.logical(v)) stop(sprintf("`%s` is not logical.", deparse2(exprs[[nm]])), call. = FALSE)
            } else {
                if (!is.numeric(v) && !is.logical(v))
                    stop(sprintf("`%s` is not numeric.", deparse2(exprs[[nm]])), call. = FALSE)
                v <- as.double(v)
            }
            chunk[[nm]] <- v
        }
        chunk
    }
    lo <- as.numeric(r$lo)
    if (!is.na(lo)) {
        hi <- as.numeric(r$hi)
        a <- lo
        while (a <= hi) {
            chunk <- DBI::dbGetQuery(con, sprintf("SELECT * FROM %s WHERE rowid >= %.0f AND rowid < %.0f",
                                                  from, a, a + chunk_size))
            ## (evaluated before the write, so that its errors are not
            ## reported as errors of dbWriteTable())
            if (nrow(chunk)) {
                chunk <- evaluate(chunk)
                write(chunk)
            }
            a <- a + chunk_size
        }
    }
    if (!created) {
        chunk <- evaluate(src$proto)
        write(chunk)
    }
    new <- design_source(con, out, NULL)
    new$parent <- src$holder
    new$label <- paste(src$label, "(with columns computed in R)")
    attr(new, "model") <- info
    attr(new, "ncol") <- ncols
    new
}

## The value of a model variable (the R expression `expr`) on a batch of
## rows, as columns named after `nm`: one numeric column, or one per column
## of a numeric matrix (nm_1, nm_2, ...: poly(), ns()); a text column for a
## factor or character value, or when the specification holds levels for
## the variable (m$factor); a logical column. The expression sees only the
## rows where its input columns (m$inputs) are present; the others get
## missing values. m$class: the class of the variable in training (one
## named element of the plan's data classes), checked on each batch.
## `prev`: the kind and columns found on earlier batches,
## which every batch must give. Returns list(cols, info = list(kind, cols)).
xtx_eval_model <- function(expr, chunk, env, m, nm, prev = NULL) {
    n <- nrow(chunk)
    ok <- if (length(m$inputs)) stats::complete.cases(chunk[m$inputs]) else rep(TRUE, n)
    bad <- function(msg) stop(sprintf("`%s` %s.", m$label, msg), call. = FALSE)
    if (!any(ok) && !is.null(prev)) {
        v <- NULL
        info <- prev
    } else {
        v <- eval(expr, chunk[ok, , drop = FALSE], env)
        ## a variable of another type than in training: the error of
        ## predict.lm() (a character variable with levels is a factor, as
        ## model.frame() makes it)
        if (!is.null(m$class))
            stats::.checkMFClasses(m$class, stats::setNames(list(if (is.character(v) && isTRUE(m$factor)) factor(v) else v),
                                                            names(m$class)))
        k <- if (is.matrix(v)) ncol(v) else 1L
        if (!is.matrix(v) && length(v) == 1L && sum(ok) != 1L) v <- rep(v, sum(ok))
        if (NROW(v) != sum(ok)) bad(sprintf("gives %d values for %d rows", NROW(v), sum(ok)))
        kind <- if (is.matrix(v)) {
            if (!is.numeric(v) && !is.logical(v)) bad("gives a matrix that is not numeric")
            "numeric"
        } else if (isTRUE(m$factor) || is.factor(v) || is.character(v)) {
            "factor"
        } else if (is.logical(v)) {
            "logical"
        } else if (is.numeric(v)) {
            "numeric"
        } else {
            bad(sprintf("gives values of class %s, which a remote design cannot use", class(v)[1L]))
        }
        info <- list(kind = kind, cols = if (k == 1L) nm else paste0(nm, "_", seq_len(k)))
        if (!is.null(prev) && !identical(prev, info)) bad("gives values of different kinds on different rows")
    }
    fill <- function(x, na) {
        out <- rep(na, n)
        if (!is.null(v)) out[ok] <- x
        out
    }
    cols <- lapply(seq_along(info$cols), function(u) {
        switch(info$kind,
               numeric = fill(if (is.matrix(v)) as.double(v[, u]) else as.double(v), NA_real_),
               factor = fill(as.character(v), NA_character_),
               logical = fill(v, NA))
    })
    list(cols = stats::setNames(cols, info$cols), info = info)
}

## Resolve an argument (xtx_resolve_arg()), or, when its R expression has no
## SQL translation, return it for evaluation in R: list(sql, inputs) or
## list(expr, reason).
xtx_resolve_or_eval <- function(x, env, types, con, what, bool = FALSE) {
    if (is.null(x)) return(NULL)
    if (is.call(x) && !deparse2(x[[1L]]) %in% c("SQL", "DBI::SQL")) {
        tr <- if (bool) xtx_sql_bool(x, types, con) else xtx_sql_num(x, types, con)
        if (is.null(tr$sql)) return(list(expr = x, reason = tr$reason))
    }
    xtx_resolve_arg(x, env, types, con, what, bool)
}

## The expressions of a fit that SQL cannot compute, as one message:
## "weights `myfun(x)` cannot be translated to SQL (myfun() has no SQL
## translation)".
xtx_untranslated_msg <- function(items) {
    nm <- names(items)[1L]
    sprintf("%s `%s` cannot be translated to SQL (%s)", nm, deparse2(items[[1L]]$expr), items[[1L]]$reason)
}

## A fit with expressions evaluated in R: a message names the first one;
## `engine = "sql"` refuses.
xtx_untranslated_engine <- function(items, dots) {
    if (!length(items)) return(dots)
    why <- xtx_untranslated_msg(items)
    if (identical(dots$engine, "sql")) stop("The SQL engine cannot compute this fit: ", why, ".", call. = FALSE)
    message("Evaluated in R: ", why, ".")
    dots
}

#' Remote data sources and arguments as columns
#'
#' `is_remote_source()` tells whether `x` is remote data that
#' [design_matrix()] and the fits accept: a DuckDB connection (`DBI`
#' connection) or Parquet paths (a character vector of paths ending in
#' `.parquet`, or a directory of them). `design_resolve()` turns arguments of
#' a fit, given unevaluated, into columns over the rows of a design: for a
#' remote design, column names, `DBI::SQL()` expressions, or R expressions
#' over the columns (those with a SQL translation run in the database; the
#' others are evaluated in R, batch by batch, into a temporary table, which
#' becomes the source of the design returned); in memory, the arguments are
#' evaluated in `data`.
#'
#' Names may change before the first release.
#'
#' @param x For `is_remote_source()`, any value. For `design_resolve()`, a
#'   design (a `design_matrix` or a remote design).
#' @param exprs A named list of unevaluated arguments (a column name, a
#'   `DBI::SQL()` expression, an R expression), `NULL` for an absent one.
#' @param env The environment in which R expressions are evaluated.
#' @param logical The names of the arguments that are logical (such as a
#'   subset); the others are numeric.
#' @param data In memory, the data frame of the rows of `x` in which the
#'   arguments are evaluated.
#' @return `is_remote_source()`: `TRUE` or `FALSE`. `design_resolve()`: a
#'   list with `x` (the design, on a new source when an argument was
#'   evaluated in R), `cols` (a named list of columns: vectors in memory,
#'   [design_column()]s remotely; `NULL` for absent arguments) and
#'   `evaluated` (the arguments evaluated in R, with the reason).
#' @examples
#' is_remote_source(data.frame(x = 1))
#' d <- data.frame(y = rnorm(5), x = 1:5)
#' D <- design_matrix(y ~ x, d)
#' design_resolve(D, list(w = quote(x / 2)), environment(), data = d)$cols$w
#' @export
is_remote_source <- function(x) xtx_is_remote_source(x)

#' @rdname is_remote_source
#' @export
design_resolve <- function(x, exprs, env = parent.frame(), logical = character(), data = NULL) {
    checkmate::assert_multi_class(x, c("xtx_matrix", "xtx_remote"))
    checkmate::assert_list(exprs, names = "unique")
    xtx_p_resolve(x, exprs, env, logical = logical, data = data)
}

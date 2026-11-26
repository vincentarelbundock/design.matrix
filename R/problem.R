## The input stage of the fits (design_problem()): a formula and a data frame
## or remote data (see remote.R), a design matrix or a remote design, turned
## into the design of the rows used and the response, weights, offset and
## subset as columns over them (vectors in memory; in SQL remotely, or
## evaluated in R into a temporary table), with one code path for linear
## models, GLMs and fixed effects. Also design_resolve() (arguments of a
## fit as columns of a design) and design_exclude() (na.exclude remotely).

## The arguments of remote data that a fit on a formula takes in `...`.
xtx_remote_keys <- c("table", "sql", "engine", "chunk_size")

#' The input stage of a fit
#'
#' `design_problem()` turns the arguments of a fit (a formula and a data frame
#' or remote data, a design matrix, or a remote design, with the response,
#' weights, offset and subset) into the design of the rows used and the
#' response, weights and offset as columns over those rows: numeric vectors in
#' memory, [design_column()]s remotely (expressions with no SQL translation are
#' evaluated in R, in batches, into a temporary table). One code path serves
#' [design_lm()], [design_glm()] and fits with fixed effects: a fitting
#' function calls it with its matched call and runs the same code on what it
#' returns, whatever the source of the data.
#'
#' The sources:
#' * a formula and a data frame: a model frame, as [lm()] and [glm()] build
#'   it (`subset`, `weights`, `na.action` and `offset` of the call evaluated
#'   as [model.frame()] evaluates them);
#' * a formula and remote data (a DuckDB connection, with `table` or `sql` in
#'   `dots`, or Parquet paths): a remote design of the rows where no model
#'   variable, argument or `fe`/`cluster` variable is missing;
#' * a design matrix, with `y`, `weights` and `offset` vectors (taken from
#'   `frame`);
#' * a remote design, with `y`, `weights` and `offset` columns or expressions
#'   over its data.
#'
#' Messages are those of [lm()] (a numeric response, no `fe`) or [glm()]
#' (otherwise), in memory and remotely alike.
#'
#' Names may change before the first release.
#'
#' @param call The matched call of the fit ([match.call()]): `weights`,
#'   `offset`, `subset`, `na.action` and `y` are taken from it, unevaluated.
#' @param env The environment of the call (the fit's [parent.frame()]).
#' @param formula A model formula (with fixed effects, the formula of the
#'   main terms), a design matrix or a remote design.
#' @param data A data frame, remote data, or `NULL`.
#' @param contrasts As in [stats::model.matrix()].
#' @param dots The arguments of remote data given in `...` of the fit
#'   (`table`, `sql`, `engine`, `chunk_size`), as a named list; other names
#'   are an error ("Unknown arguments").
#' @param family `NULL` for a numeric response (a linear model), or a
#'   [column_family()]: the response is then taken as [glm()] takes it (a
#'   factor, a two-column matrix of successes and failures), and its columns
#'   are returned as `cols`.
#' @param fe `NULL`, or a named list of character vectors: the variables of
#'   each fixed effect (several for an interaction). Rows with a missing value
#'   in them are dropped, and their values over the rows are returned as
#'   `keys` (for [design_absorb_rows()]).
#' @param cluster `NULL`, or the name of a cluster variable (as `fe`).
#' @param frame The fit's own environment ([environment()]), from which the
#'   vectors given with a design matrix are taken.
#' @return A list of class `xtx_problem`: `x` (the design of the rows used),
#'   `y`, `w` and `off` (the response, weights and offset, or `NULL`), `cols`
#'   (with `family`: the response columns as [glm()]'s `initialize` sees
#'   them: `y`, `w`, `n` (trials), `off`, `cbind`, `a`, `b`, `weighted`,
#'   `y0`), `keys` (with `fe`), `frame` (the model frame, in memory),
#'   `na.action`, `exclude` (remotely, `na.exclude`), `rows`, `nall` (the
#'   number of rows used), `terms`, `formula`, `contrasts`, `xlevels`, `env`
#'   (the formula's environment), `intercept`, `offset_call` (the `offset`
#'   argument as given), `design` (whether the first argument was a design),
#'   `fe` and `cluster`.
#' @examples
#' d <- data.frame(y = rnorm(20), x = rnorm(20), w = runif(20))
#' fit <- function(formula, data, weights, ...) {
#'     design_problem(match.call(), parent.frame(), formula, data)
#' }
#' P <- fit(y ~ x, d, weights = w)
#' dim(P$x)
#' all.equal(P$w, d$w)
#' @export
design_problem <- function(call, env, formula, data = NULL, contrasts = NULL, dots = list(), family = NULL,
                           fe = NULL, cluster = NULL, frame = NULL) {
    P <- if (inherits(formula, "xtx_matrix")) {
        xtx_problem_design(xtx_frame_arg(call, frame, "y"), xtx_frame_arg(call, frame, "weights"),
                          xtx_frame_arg(call, frame, "offset"), formula, family)
    } else if (inherits(formula, "xtx_remote")) {
        xtx_problem_rdesign(formula, if (is.null(family)) call[["y"]] else xtx_frame_arg(call, frame, "y"),
                           call[["weights"]], call[["offset"]], env, family)
    } else if (!is.null(data) && xtx_is_remote_source(data)) {
        xtx_problem_remote(call, env, formula, data, contrasts, dots, family, fe, cluster)
    } else {
        if (length(dots)) stop("Unknown arguments: ", paste(names(dots), collapse = ", "), call. = FALSE)
        xtx_problem_memory(call, env, formula, contrasts, family, fe, cluster)
    }
    P$design <- inherits(formula, c("xtx_matrix", "xtx_remote"))
    if (!P$design) P$offset_call <- call$offset
    if (!is.null(fe)) P$keys <- xtx_fe_keys(P$x, fe, cluster, P$frame)
    if (!is.null(family)) {
        P$cols <- if (inherits(P$x, "xtx_remote")) {
            xtx_glm_cols(P$y, P$resp$a, P$resp$b, if (is.null(P$w)) xtx_c_const(P$x, 1) else P$w, P$off, family,
                       P$resp$cbind, weighted = !is.null(P$w))
        } else {
            xtx_glm_mem_cols(P$y, P$w, P$off, family)
        }
        ## (on a design no rows are dropped: missing values are errors)
        if (P$design) xtx_glm_check_design(P$cols)
    }
    P$fe <- fe
    P$cluster <- cluster
    class(P) <- "xtx_problem"
    P
}

## The value of argument `nm` of the fit whose environment is `frame` (NULL
## when it was not given): the vectors given with a design matrix, as the
## fit received them.
xtx_frame_arg <- function(call, frame, nm) {
    if (is.null(call[[nm]])) return(NULL)
    get(nm, envir = frame)
}

## A design matrix and the vectors `y`, `weights` and `offset` (as lm.fit()
## and glm.fit() take them; for a GLM, `y` as glm() takes the response: a
## factor, a two-column matrix).
xtx_problem_design <- function(y, weights, offset, x, family) {
    checkmate::assert_class(x, "xtx_matrix")
    n <- nrow(x)
    if (is.null(family)) {
        if (is.matrix(y) && ncol(y) == 1L) y <- drop(y)
        checkmate::assert_numeric(y, len = n)
        y <- as.double(y)
    } else {
        if (is.null(y)) stop("Give the response as `y`.", call. = FALSE)
        if (NROW(y) != n) stop(sprintf("The response must have %d values.", n), call. = FALSE)
    }
    if (!is.null(weights)) checkmate::assert_numeric(weights, len = n)
    if (!is.null(offset)) checkmate::assert_numeric(offset, len = n)
    tt <- x$training_terms
    list(x = x, y = y, w = weights, off = offset, exclude = FALSE, nall = n,
         intercept = if (!is.null(tt)) attr(tt, "intercept") > 0L else "(Intercept)" %in% colnames(x))
}

## A remote design and its arguments `y`, `weights` and `offset`
## (unevaluated: columns, DBI::SQL() or R expressions over its data,
## resolved in env; a GLM's `y`: the names of one column, or two
## (successes and failures); NULL: the response of the design's formula).
xtx_problem_rdesign <- function(D, y, weights, offset, env, family) {
    if (!is.null(D$scale)) stop("Fit an unscaled remote design, with `weights` naming a column.", call. = FALSE)
    if (is.null(family)) {
        if (is.null(y)) {
            if (is.null(D$response))
                stop("With a remote design, give `y`, the name of one column (or use a design with a response).",
                     call. = FALSE)
            v <- D$vars[[D$response]]
            if (v$kind != "numeric")
                stop("The response of a remote fit must be a numeric variable with a SQL translation.", call. = FALSE)
            y <- DBI::SQL(v$sql)
        }
        resp <- list(args = list(y = y), present = character(), cbind = FALSE)
    } else {
        if (!identical(D$na, "omit"))
            stop("Fit a remote design built with `na.action = na.omit` (or from a formula).", call. = FALSE)
        lhs <- if (!is.null(y)) {
            xtx_remote_column(D, y, "y", several = TRUE)
            if (length(y) == 2L) call("cbind", as.name(y[1L]), as.name(y[2L]))
            else if (length(y) == 1L) as.name(y)
            else stop("`y` must name one column, or two (successes and failures).", call. = FALSE)
        } else if (!is.null(D$response)) {
            D$vars[[D$response]]$expr
        } else {
            stop("With a remote design, give `y`, the name of the response column (or use a design with a response).",
                 call. = FALSE)
        }
        resp <- xtx_glm_remote_response(lhs, D$src$types, xtx_con(D))
        xtx_glm_check_arity(resp)
    }
    intercept <- attr(D$plan$terms, "intercept") > 0L && !isTRUE(D$plan$drop_intercept)
    r <- xtx_p_resolve(D, c(resp$args, list(weights = weights, offset = offset)), env)
    ## (a linear model names the expressions evaluated in R)
    if (is.null(family) && length(r$evaluated))
        message(xtx_untranslated_msg(r$evaluated), ": evaluated in R, in batches of rows.")
    D <- r$x
    xtx_glm_check_factor(resp, family)
    list(x = D, y = if (length(resp$present)) xtx_glm_remote_factor(D, resp$present) else r$cols$y,
         w = r$cols$weights, off = r$cols$offset, resp = list(a = r$cols$a, b = r$cols$b, cbind = resp$cbind),
         exclude = FALSE, nall = D$n, contrasts = D$plan$contrasts, xlevels = D$plan$xlevels, env = env,
         intercept = intercept, offset_call = offset)
}

## In memory: the model frame (with the FE and cluster variables, so that
## rows with a missing value in any of them are dropped), its terms (of
## the main terms with FE), response (as lm() takes it; a GLM's as in the
## frame: a vector, factor or matrix), weights, offset (offset() terms and
## the `offset` argument) and the design (with the frame's predvars).
xtx_problem_memory <- function(call, env, formula, contrasts, family, fe, cluster) {
    lm <- is.null(family) && is.null(fe)
    mf <- call[c(1L, match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(call), 0L))]
    if (!is.null(fe)) {
        aug <- xtx_fe_augment(formula, unique(c(unlist(fe, use.names = FALSE), cluster)))
        environment(aug) <- environment(formula)
        mf$formula <- aug
    }
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(xtx::design_frame)
    mf <- eval(mf, env)
    mt <- attr(mf, "terms")
    if (!is.null(fe)) mt <- xtx_fe_main_terms(mt, formula)
    if (attr(attr(mf, "terms"), "response") == 0L) stop("The formula has no response.", call. = FALSE)
    if (lm) {
        y <- xtx_fit_response(mf)
        if (is.matrix(y)) stop("Multivariate responses are not supported.", call. = FALSE)
    } else {
        y <- mf[[1L]]
        if (!is.null(fe) && (is.numeric(y) || is.logical(y))) y <- xtx_fit_response(mf)
    }
    w <- as.vector(stats::model.weights(mf))
    if (!is.null(w) && !is.numeric(w)) stop("'weights' must be a numeric vector", call. = FALSE)
    ## (weights with several columns: the same error remotely)
    if (!is.null(w) && length(w) != nrow(mf)) stop(xtx_scale_length_msg(nrow(mf)), call. = FALSE)
    off <- as.vector(stats::model.offset(mf))
    if (!is.null(off)) xtx_offset_length(length(off), NROW(y), lm)
    if (!lm) {
        if (!is.null(w)) w <- as.double(w)
        if (!is.null(off)) off <- as.double(off)
    }
    X <- design_matrix(mt, mf, contrasts.arg = contrasts)
    list(x = X, y = y, w = w, off = off, frame = mf, na.action = attr(mf, "na.action"), exclude = FALSE,
         rows = seq_len(nrow(mf)), nall = nrow(mf), terms = mt, formula = stats::formula(mt),
         contrasts = X$contrasts, xlevels = stats::.getXlevels(mt, mf), env = environment(mt),
         intercept = attr(mt, "intercept") > 0L)
}

## An offset of the wrong length (a matrix): the error of lm() for a
## linear model without FE, of glm() otherwise.
xtx_offset_length <- function(noff, n, lm) {
    if (noff != n)
        stop(gettextf(if (lm) "number of offsets is %d, should equal %d (number of observations)"
                      else "number of offsets is %d should equal %d (number of observations)",
                      noff, n), call. = FALSE)
}

## The terms of the main formula of an FE fit, from those of its model
## frame (which has the FE and cluster variables too).
xtx_fe_main_terms <- function(mt_aug, main) {
    labs <- attr(mt_aug, "term.labels")
    drop <- which(!labs %in% attr(stats::terms(main), "term.labels"))
    mt <- if (!length(drop)) {
        mt_aug
    } else if (length(drop) < length(labs)) {
        stats::drop.terms(mt_aug, drop, keep.response = TRUE)
    } else {
        ## no main term left: a formula with the response, offset() terms
        ## and the intercept
        tt <- stats::terms(main)
        nv <- length(attr(tt, "variables"))
        attr(tt, "predvars") <- attr(mt_aug, "predvars")[seq_len(nv)]
        attr(tt, "dataClasses") <- attr(mt_aug, "dataClasses")[seq_len(nv - 1L)]
        tt
    }
    environment(mt) <- environment(main)
    mt
}

## Remotely: the remote design of the right-hand side, and the response,
## weights, offset and subset (linear models without FE) resolved as
## columns (xtx_remote_fit()), on the rows where the model, FE and cluster
## variables, the response, weights and offset are present, as the model
## frame in memory. A GLM's response: xtx_glm_remote_response() (a factor
## as 0/1, the columns of cbind(a, b) as `resp`). The errors are those of
## the in-memory path, in its order (xtx_remote_types()): the types of the
## columns before their SQL translation; the lengths (arguments evaluated
## in R that give matrices) and the response of a GLM or an FE fit once
## the rows are known.
xtx_problem_remote <- function(call, env, formula, source, contrasts, dots, family, fe, cluster) {
    ## (the callers check the kind: checks.R)
    na <- if (is.null(call$na.action)) "omit" else xtx_problem_na_kind(eval(call$na.action, env))
    ## (unknown arguments first, as in memory)
    nms <- names(dots) %||% rep("", length(dots))
    if (any(!nms %in% xtx_remote_keys))
        stop("Unknown arguments: ", paste(nms[!nms %in% xtx_remote_keys], collapse = ", "), call. = FALSE)
    xtx_need_duckdb()
    src <- xtx_remote_source(source, dots$table, dots$sql)
    lm <- is.null(family) && is.null(fe)
    ## (the response, weights, offset and subset are evaluated in the
    ## formula's environment, as model.frame() does in memory)
    fenv <- environment(formula) %||% env
    tt <- stats::terms(formula, data = src$proto)
    if (attr(tt, "response") == 0L) stop("The formula has no response.", call. = FALSE)
    lhs <- attr(tt, "variables")[[attr(tt, "response") + 1L]]
    resp <- if (!is.null(family)) xtx_glm_remote_response(lhs, src$types, src$holder$con)
            else list(args = list(y = lhs), present = character(), cbind = FALSE)
    ## a linear model has one response, as in memory (a matrix response
    ## evaluated in R gives the same error, xtx_remote_fit())
    mresp <- if (is.null(family)) "Multivariate responses are not supported."
    if (!is.null(mresp) && is.call(lhs) && identical(lhs[[1L]], as.name("cbind")) && length(lhs) > 2L)
        stop(mresp, call. = FALSE)
    tv <- xtx_terms_vars(tt)
    offt <- lapply(attr(tt, "offset"), function(k) tv$exprs[[k]][[2L]])
    late <- xtx_remote_types(src$types, resp$args$y, call$weights, call$offset, offt, family, lm)
    if (!is.null(late)) {
        ## (a response column of the wrong type, reported below as in
        ## memory: until then it only selects the rows where it is present)
        resp$extra <- c(resp$extra, as.character(resp$args$y))
        resp$args$y <- NULL
    }
    args <- c(resp$args, list(weights = call$weights, offset = call$offset, subset = if (is.null(family)) call$subset))
    r <- xtx_remote_fit(formula, src, args, fenv, dots, contrasts, logical = "subset",
                       present = unique(c(unlist(fe, use.names = FALSE), cluster, resp$present, resp$extra)),
                       na = na, mresp = mresp)
    D <- r$x
    ## weights or an offset with several columns: the errors of the
    ## in-memory path, with the number of rows used
    if (r$wk > 1L) stop(xtx_scale_length_msg(D$n), call. = FALSE)
    if (r$offk > 1L) xtx_offset_length(r$offk * D$n, D$n, lm)
    xtx_glm_check_arity(resp)
    xtx_glm_check_factor(resp, family)
    if (!is.null(late)) stop(late, call. = FALSE)
    y <- if (length(resp$present)) xtx_glm_remote_factor(D, resp$present) else r$cols$y
    list(x = D, y = y, w = r$cols$weights, off = r$off, resp = list(a = r$cols$a, b = r$cols$b, cbind = resp$cbind),
         na.action = NULL, exclude = na == "exclude", rows = NULL, nall = D$n,
         terms = tt, formula = formula, contrasts = D$plan$contrasts,
         xlevels = D$plan$xlevels, env = fenv, intercept = attr(tt, "intercept") > 0L)
}

## The types of the columns that a fit on remote data takes as its
## response (`y`), `weights` and `offset` (and the expressions of its
## offset() terms, `offt`), checked before their SQL translation, with the
## errors of the in-memory path (lm(), glm(), model.offset()) in its order:
## a linear model's response (without FE), the weights, the offset. Only
## bare names of columns are checked here (an expression over a column of
## the wrong type fails when evaluated, as in memory). Returns the error of
## a response of the wrong type that the in-memory path reports later
## (with FE or for a GLM: after the lengths of the weights and offset), or
## NULL.
xtx_remote_types <- function(types, y, weights, offset, offt, family, lm) {
    ty <- function(e) if (is.name(e) && as.character(e) %in% names(types)) types[[as.character(e)]] else NA_character_
    num <- c("double", "integer")
    yt <- ty(y)
    ## (a GLM's character, factor and logical responses are not `y`:
    ## xtx_glm_remote_response())
    ymsg <- if (is.na(yt)) NULL
            else if (is.null(family)) { if (!yt %in% c(num, "logical")) "The response must be numeric." }
            else if (yt == "other") "The response must be numeric, logical or a factor."
    if (lm && !is.null(ymsg)) stop(ymsg, call. = FALSE)
    wt <- ty(weights)
    if (!is.na(wt) && !wt %in% num) stop("'weights' must be a numeric vector", call. = FALSE)
    omsg <- xtx_offset_type_msg(ty(offset), vapply(offt, ty, ""))
    if (!is.null(omsg)) stop(omsg, call. = FALSE)
    ymsg
}

## The error of model.offset() for an `offset` argument and offset() terms
## of these column types (NA: not a column, taken as numeric), or NULL: the
## argument alone must be numeric; with offset() terms, model.offset()
## adds them up (`+` fails on text, gives NA on a factor) and the sum must
## be numeric.
xtx_offset_type_msg <- function(arg, terms) {
    num <- c("double", "integer")
    if (!length(terms)) return(if (!is.na(arg) && !arg %in% num) "'offset' must be numeric")
    s <- if (is.na(arg)) "double" else arg
    for (t in terms) {
        if (is.na(t)) t <- "double"
        if ("character" %in% c(s, t)) return("non-numeric argument to binary operator")
        s <- if ("factor" %in% c(s, t)) "logical" else if ("other" %in% c(s, t)) "other"
             else if (s == "logical" && t == "logical") "integer" else "double"
    }
    if (!s %in% num) "'offset' must be numeric"
}

## The remote design of a fit and the columns of its arguments (one code
## path for design_lm(), design_glm() and fixed effects).
##   formula: the model formula; the design is that of its right-hand side
##     (offset() terms included; their values are columns too), and its
##     response is given by the caller in `args`.
##   src: an opened source (xtx_remote_source()).
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
xtx_remote_fit <- function(formula, src, args, env, dots, contrasts = NULL, logical = character(),
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

## The response of a remote GLM as arguments: list(args, present (a
## factor response), cbind, ncol (the columns of cbind()), extra (columns
## of cbind() beyond two, which only select rows)). cbind() of one column
## is that column, as a one-column matrix in memory (xtx_glm_mem_cols();
## but cbind() of a factor gives its codes there: refused); of more than
## two, an error once the rows are known (xtx_glm_check_arity()), as in
## memory.
xtx_glm_remote_response <- function(lhs, types, con) {
    cols <- function(es) vapply(Filter(function(e) is.name(e) && as.character(e) %in% names(types), es),
                                as.character, "")
    if (is.call(lhs) && identical(lhs[[1L]], as.name("cbind"))) {
        if (length(lhs) == 2L && !(is.name(lhs[[2L]]) && isTRUE(types[as.character(lhs[[2L]])] == "factor")))
            return(xtx_glm_remote_response(lhs[[2L]], types, con))
        if (length(lhs) < 3L)
            return(list(args = list(), present = character(), cbind = TRUE, ncol = length(lhs) - 1L,
                        extra = cols(as.list(lhs)[-1L])))
        return(list(args = list(a = lhs[[2L]], b = lhs[[3L]]), present = character(), cbind = TRUE,
                    ncol = length(lhs) - 1L, extra = cols(as.list(lhs)[-(1:3)])))
    }
    if (is.name(lhs) && isTRUE(types[as.character(lhs)] %in% c("character", "factor")))
        return(list(args = list(), present = as.character(lhs), cbind = FALSE))
    if (is.name(lhs) && isTRUE(types[as.character(lhs)] == "logical"))
        lhs <- DBI::SQL(sprintf("CAST(s.%s AS DOUBLE)", DBI::dbQuoteIdentifier(con, as.character(lhs))))
    list(args = list(y = lhs), present = character(), cbind = FALSE)
}

## cbind() with more than two columns (xtx_glm_remote_response()): the
## error of the in-memory path (xtx_glm_mem_cols()).
xtx_glm_check_arity <- function(resp) {
    if (isTRUE(resp$cbind) && resp$ncol != 2L)
        stop("A two-column response needs two columns (successes and failures).", call. = FALSE)
}

## A character or factor response (xtx_glm_remote_response()) needs a
## binomial family, as in memory (xtx_glm_mem_cols()).
xtx_glm_check_factor <- function(resp, k) {
    if (length(resp$present) && !is.null(k) && !k$binomial)
        stop("A factor response needs a binomial family.", call. = FALSE)
}

## A character or factor response as a column: 0 for the first level of the
## rows used, 1 for the others.
xtx_glm_remote_factor <- function(D, col) {
    q <- xtx_qi(D, col)
    vals <- xtx_query(D, sprintf("SELECT DISTINCT s.%s AS v FROM (%s) AS s WHERE %s", q, D$src$select, D$filter))$v
    proto <- D$src$proto[[col]]
    f <- if (is.factor(proto)) factor(as.character(vals), levels = levels(proto)) else factor(vals)
    first <- DBI::dbQuoteString(xtx_con(D), levels(droplevels(f))[1L])
    xtx_column(D, sprintf("CASE WHEN CAST(s.%s AS VARCHAR) = %s THEN 0 ELSE 1 END", q, first))
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
## Returns the new source (xtx_remote_source() of the temporary table),
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
    new <- xtx_remote_source(con, out, NULL)
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


## The response of a model frame, as model.response(mf, "numeric") but
## without names: on a large frame, turning compact row names into a
## character vector costs more than the fit.
xtx_fit_response <- function(mf) {
    if (attr(attr(mf, "terms"), "response") == 0L) return(NULL)
    v <- mf[[1L]]
    if (!is.numeric(v) && !is.logical(v))
        stop("The response must be numeric.", call. = FALSE)
    if (is.null(dim(v))) {
        ## (no copy of a plain double vector)
        if (!is.null(attributes(v))) attributes(v) <- NULL
        if (!is.double(v)) storage.mode(v) <- "double"
    } else if (ncol(v) == 1L) {
        v <- as.double(v)
    }
    v
}


## The formula with variables added as terms (FE and cluster variables, for
## the model frame and the rows used).
xtx_fe_augment <- function(main, vars) {
    rhs <- main[[length(main)]]
    for (v in vars) rhs <- call("+", rhs, as.name(v))
    main[[length(main)]] <- rhs
    main
}

## The keys of the FE (and of the clusters, when they are not an FE): the
## values of their variables over the rows of x (from `data` in memory).
xtx_fe_keys <- function(x, fe, cluster, data = NULL) {
    list(fe = lapply(fe, function(v) xtx_c_key(x, v, data)),
         cl = if (!is.null(cluster) && is.na(match(cluster, names(fe)))) xtx_c_key(x, cluster, data))
}


#' na.exclude for remote fits
#'
#' `design_exclude()` is the `na.action` of a remote fit with
#' `na.action = na.exclude`: results on the rows a fit used (remote tables,
#' [design_output()]) are padded with the rows of the data that it dropped for
#' missing values, as [napredict()] and [naresid()] pad vectors in memory. Its
#' [stats::napredict()] and [stats::naresid()] methods return a temporary view
#' with the rows of the table and the dropped rows (with `NULL` results);
#' other values are returned unchanged.
#'
#' Names may change before the first release.
#'
#' @param x For `design_exclude()`, a remote design: the design of the rows
#'   the fit could use, over the data (before rows were dropped for other
#'   reasons). For the methods, a result of a fit: a `xtx_remote_table`, or
#'   any other value.
#' @param omit A `xtx_remote_exclude` object.
#' @param ... Not used.
#' @return `design_exclude()`: an object of class `xtx_remote_exclude`.
#'   The methods: a `xtx_remote_table` (the name of a temporary view), or `x`.
#' @export
design_exclude <- function(x) {
    if (!inherits(x, "xtx_remote")) stop("`x` must be a remote design.", call. = FALSE)
    xtx_remote_exclude(x)
}

## (see design_exclude())
xtx_remote_exclude <- function(D) structure(list(design = D), class = "xtx_remote_exclude")

#' @rdname design_exclude
#' @exportS3Method stats::napredict
napredict.xtx_remote_exclude <- function(omit, x, ...) {
    if (!inherits(x, "xtx_remote_table")) return(x)
    D0 <- omit$design
    con <- attr(x, "connection")
    sel <- sprintf("SELECT * FROM %s UNION ALL BY NAME SELECT s.* FROM (%s) AS s WHERE NOT COALESCE(%s, FALSE)",
                   DBI::dbQuoteIdentifier(con, unclass(x)[1L]), D0$src$select, D0$filter)
    name <- xtx_tmp_name("xtx_out")
    DBI::dbExecute(con, sprintf("CREATE TEMPORARY VIEW %s AS %s", DBI::dbQuoteIdentifier(con, name), sel))
    structure(name, class = "xtx_remote_table", connection = con, holder = attr(x, "holder"))
}

#' @rdname design_exclude
#' @exportS3Method stats::naresid
naresid.xtx_remote_exclude <- function(omit, x, ...) napredict.xtx_remote_exclude(omit, x, ...)


## The response, prior weights, trials and offset as glm()'s initialize sees
## them (columns: vectors in memory, SQL remotely): y in [0, 1] and weights
## times trials for a two-column binomial response (a successes, b
## failures), y = 0 where the weight is 0. `weighted`: whether weights
## were given (w is then theirs, else 1); a two-column response is
## weighted by its trials.
xtx_glm_cols <- function(y, a, b, w, off, k, cb, weighted) {
    y0 <- y
    if (cb) {
        if (!k$binomial)
            stop("A two-column response (successes and failures) needs a binomial family.", call. = FALSE)
        nt <- a + b
        y <- xtx_c_ifelse(nt == 0, 0, a / nt)
        w <- w * nt
    } else {
        nt <- xtx_c_const(y, 1)
        if (k$binomial) y <- xtx_c_ifelse(w == 0, 0, y)
    }
    ## (trials weight the rows of a two-column response; y0 the response as
    ## given)
    list(y = y, w = w, n = nt, off = off, cbind = cb, a = a, b = b, weighted = weighted || cb, y0 = y0)
}

## The response of a model frame in memory (a vector, a factor or character
## vector, or a two-column matrix) as columns.
xtx_glm_mem_cols <- function(Y, w, off, k) {
    if (length(dim(Y)) == 1L) dim(Y) <- NULL
    if (is.matrix(Y) && ncol(Y) == 1L) Y <- Y[, 1L]
    n <- NROW(Y)
    weighted <- !is.null(w)
    w <- if (is.null(w)) rep.int(1, n) else as.double(w)
    off <- if (!is.null(off)) as.double(off)
    if (is.matrix(Y)) {
        if (ncol(Y) != 2L) stop("A two-column response needs two columns (successes and failures).", call. = FALSE)
        return(xtx_glm_cols(NULL, as.double(Y[, 1L]), as.double(Y[, 2L]), w, off, k, TRUE, weighted))
    }
    ## a character response is a factor, as remotely
    if (is.character(Y)) Y <- factor(Y)
    if (is.factor(Y)) {
        if (!k$binomial) stop("A factor response needs a binomial family.", call. = FALSE)
        Y <- Y != levels(Y)[1L]
    }
    if (!is.numeric(Y) && !is.logical(Y)) stop("The response must be numeric, logical or a factor.", call. = FALSE)
    y <- as.double(Y)
    attributes(y) <- NULL
    xtx_glm_cols(y, NULL, NULL, w, off, k, FALSE, weighted)
}

## The checks of a fit on a design (no rows dropped): no missing response
## or weight, a finite offset. A missing weight has the error of a weight
## that is not finite ("NA/NaN/Inf in 'weights'", xtx_glm_checks()), in
## memory as remotely: glm() has no error of its own for it, and lm()'s
## "missing or negative weights not allowed" would contradict the GLM's
## "negative weights not allowed" for negative weights.
xtx_glm_check_design <- function(cols) {
    y <- if (cols$cbind) cols$a + cols$b else cols$y0
    s <- xtx_c_sums(xtx_c_count(xtx_c_present(y), negate = TRUE),
                  if (cols$weighted && !cols$cbind) xtx_c_count(xtx_c_present(cols$w), negate = TRUE) else 0,
                  if (!is.null(cols$off)) xtx_c_count(xtx_c_is_finite(cols$off), negate = TRUE) else 0)
    bad <- is.na(s) | s > 0
    if (bad[1L]) stop("NA/NaN/Inf in 'y'", call. = FALSE)
    if (bad[2L]) stop("NA/NaN/Inf in 'weights'", call. = FALSE)
    if (bad[3L]) stop("NA/NaN/Inf in 'offset'", call. = FALSE)
}


## na.action of remote data as one of "omit", "exclude", "fail", "pass"
## (the fits check it first, with their messages).
xtx_problem_na_kind <- function(na.action) {
    if (is.null(na.action)) return("omit")
    for (nm in c("omit", "exclude", "fail", "pass")) {
        f <- get(paste0("na.", nm), envir = asNamespace("stats"))
        if (identical(na.action, f) || identical(na.action, paste0("na.", nm))) return(nm)
    }
    stop("`na.action` must be na.omit, na.exclude, na.fail, na.pass.", call. = FALSE)
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

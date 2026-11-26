## The input stage of the fits (design_problem()): a formula and a data frame
## or remote data (see remote.R), a design matrix or a remote design, turned
## into the design of the rows used and the response, weights, offset and
## subset as columns over them (vectors in memory; in SQL remotely, or
## evaluated in R into a temporary table), with one code path for linear
## models, GLMs and fixed effects. Also design_resolve() (arguments of a
## fit as columns of a design) and design_exclude() (na.exclude remotely).

## The arguments of remote data that a fit on a formula takes in `...`.
xtx_remote_keys <- c("table", "sql", "engine", "chunk_size")

## The input stage of the fits (design_lm(), design_glm() and fixed
## effects): the call's arguments as a design and the columns of the
## response, weights and offset, with one path per source of the design:
##  * a formula and a data frame: a model frame, as lm() and glm() build it
##    (xtx_problem_memory());
##  * a formula and remote data: a remote design, xtx_remote_fit()
##    (xtx_problem_remote());
##  * a design matrix, with `y`, `weights` and `offset` vectors
##    (xtx_problem_design());
##  * a remote design, with `y`, `weights` and `offset` columns or
##    expressions over its data (xtx_problem_rdesign()).
## The fits then run the same code on any of them.
##   call: the matched call of the fit (its arguments are taken from it,
##     unevaluated), env its parent frame; formula and data their values
##     (formula: a formula or a design; data NULL when missing); dots the
##     arguments in `...`.
##   family: NULL for a linear model; for a GLM, its xtx_glm_kind(): the
##     response as glm() takes it (a factor, a two-column response), and
##     its columns (xtx_glm_cols()) as `cols`.
##   fe, clv: with fixed effects, the parsed formula (xtx_fe_parse()) and
##     the cluster variable (or NULL): the design is that of the main terms,
##     rows with a missing FE or cluster value are dropped, and `keys` holds
##     their values (xtx_fe_keys()).
##   frame: the fit's own environment (its environment()), from which the
##     values of `y`, `weights` and `offset` given with a design matrix
##     (and a GLM's `y` with a remote design) are taken.
## The response, weights, offset and subset of a formula are evaluated in
## the data and the formula's environment (as model.frame() does),
## na.action in env (its kind is checked by the caller: checks.R).
## Messages are the same in memory and remotely: for the offset length,
## those of lm() for a linear model without FE, of glm() otherwise.
## Returns list(x (the design), y (the response), w, off (or NULL), cols
## (GLM), resp (remotely: the columns of a two-column response), keys,
## frame (the model frame in memory), na.action, exclude (remotely,
## na.exclude), rows, nall, terms, formula, contrasts, xlevels, env (the
## formula's environment), intercept (whether the design has one, for the
## fit), offset_call (the `offset` argument as given, for new data), design
## (whether the first argument was a design)); see xtx_na_result() for the
## na.action of the results.
design_problem <- function(call, env, formula, data, contrasts, dots, family = NULL, fe = NULL, clv = NULL,
                           frame = NULL) {
    P <- if (inherits(formula, "xtx_matrix")) {
        xtx_problem_design(xtx_frame_arg(call, frame, "y"), xtx_frame_arg(call, frame, "weights"),
                          xtx_frame_arg(call, frame, "offset"), formula, family)
    } else if (inherits(formula, "xtx_remote")) {
        xtx_problem_rdesign(formula, if (is.null(family)) call[["y"]] else xtx_frame_arg(call, frame, "y"),
                           call[["weights"]], call[["offset"]], env, family)
    } else if (!is.null(data) && xtx_is_remote_source(data)) {
        xtx_problem_remote(call, env, formula, data, contrasts, dots, family, fe, clv)
    } else {
        if (length(dots)) stop("Unknown arguments: ", paste(names(dots), collapse = ", "), call. = FALSE)
        xtx_problem_memory(call, env, contrasts, family, fe, clv)
    }
    P$design <- inherits(formula, c("xtx_matrix", "xtx_remote"))
    if (!P$design) P$offset_call <- call$offset
    if (!is.null(fe)) P$keys <- xtx_fe_keys(P$x, fe, clv, P$frame)
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
xtx_problem_memory <- function(call, env, contrasts, family, fe, clv) {
    lm <- is.null(family) && is.null(fe)
    mf <- call[c(1L, match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(call), 0L))]
    if (!is.null(fe)) {
        aug <- xtx_fe_augment(fe$main, unique(c(unlist(fe$fe, use.names = FALSE), clv)))
        environment(aug) <- environment(fe$main)
        mf$formula <- aug
    }
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(xtx::design_frame)
    mf <- eval(mf, env)
    mt <- attr(mf, "terms")
    if (!is.null(fe)) mt <- xtx_fe_main_terms(mt, fe$main)
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
xtx_problem_remote <- function(call, env, formula, source, contrasts, dots, family, fe, clv) {
    ## (the callers check the kind: checks.R)
    na <- if (is.null(call$na.action)) "omit" else xtx_na_kind(eval(call$na.action, env))
    ## (unknown arguments first, as in memory)
    nms <- names(dots) %||% rep("", length(dots))
    if (any(!nms %in% xtx_remote_keys))
        stop("Unknown arguments: ", paste(nms[!nms %in% xtx_remote_keys], collapse = ", "), call. = FALSE)
    xtx_need_duckdb()
    src <- xtx_remote_source(source, dots$table, dots$sql)
    if (!is.null(fe)) formula <- fe$main
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
                       present = unique(c(unlist(fe$fe, use.names = FALSE), clv, resp$present, resp$extra)),
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
xtx_fe_keys <- function(x, pf, clv, data = NULL) {
    list(fe = lapply(pf$fe, function(v) xtx_c_key(x, v, data)),
         cl = if (!is.null(clv) && is.na(match(clv, names(pf$fe)))) xtx_c_key(x, clv, data))
}


## na.exclude for a remote fit: results on the rows used, padded with the
## rows of the data that the fit dropped for missing values (NULL results),
## as napredict() pads in memory (see napredict.xtx_remote_exclude()).
xtx_remote_exclude <- function(D) structure(list(design = D), class = "xtx_remote_exclude")

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

#' @exportS3Method stats::naresid
naresid.xtx_remote_exclude <- function(omit, x, ...) napredict.xtx_remote_exclude(omit, x, ...)


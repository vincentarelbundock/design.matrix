## The input side of the fits (design_lm(), design_glm() and fixed
## effects): design_problem() turns a formula and a data frame or remote
## data (see remote.R) into the design of the rows used and the response,
## weights, offset and subset as columns over them (vectors in memory; in
## SQL remotely, or evaluated in R into a temporary table), with one code
## path for the three fits. The fits themselves are the code of lm.R,
## glm.R and fe.R on the primitives of primitives.R.

#' @name design_remote
#' @section Linear models:
#' [design_lm()] works on remote data as on a data frame:
#' `design_lm(y ~ x + f, con, table = "flights", weights = w, subset = year > 2000)`,
#' or on a remote design: `design_lm(D, y = "col", weights = "w")`.
#' `weights`, `offset`, `subset` and `y` are column names (strings or bare
#' names), `DBI::SQL()` expressions, or R expressions over the columns, as
#' are the response and `offset()` terms of the formula. R expressions with
#' a SQL translation (arithmetic as for the SQL engine; for `subset`, also
#' comparisons, `&`, `|`, `!`, `%in%` with constants and `is.na()`) run in
#' the database. Any other R expression, such as `weights = myfun(x)`, is
#' evaluated in R, batch by batch, and its values are written with the
#' columns of the remote data to a temporary table, on which the fit runs
#' (one more pass over the data, with a message naming the expression;
#' `engine = "sql"` refuses). Compression stays exact, since the computed values are columns
#' of that table. `na.action` can be [na.omit()] (default), [na.exclude()],
#' [na.fail()] or [na.pass()], with the results or errors of [lm()] (with
#' fixed effects, `na.omit` or `na.exclude`, as in memory). The
#' rows used are those with no missing model variable, response, weight or
#' offset, and factor levels are those of these rows, as in memory.
#'
#' The fit is the same algorithm as in memory, on remote primitives: with
#' few design patterns it is compressed (see *Compression*); otherwise it
#' runs [qr()][design_qr] of the remote design (`X'WX` in the database,
#' columns judged with [lm()]'s rule, suspect columns from products with
#' the data), with the factor's values computed with numeric variables
#' shifted when that is an exact transform of the design (so that `x + 1e7`
#' keeps its accuracy). Coefficients, aliased columns, covariance, residual
#' standard error, degrees of freedom, `logLik()`, `AIC()`, `BIC()` and
#' `summary()` are those of [lm()] on the same rows, up to rounding.
#'
#' Results that have one value per row are tables in the database (names of
#' temporary views, of class `dm_remote_table`, with all the columns of the
#' remote data): `fitted(fit)` (columns `.fitted` and `.resid`; with
#' `na.exclude`, the rows with missing values too, with `NULL`), and
#' `predict(fit)` or `predict(fit, newdata = con, table = ...)` (column
#' `.fitted`; with `se.fit = TRUE`, `.se.fit`, and the attributes `df` and
#' `residual.scale`; with intervals, `.lwr` and `.upr`; for a rank-deficient
#' fit, `.nonestim`, the non-estimable cases flagged by `rankdeficient` as in
#' [stats::predict.lm()], with its warnings), computed by the same code as
#' in memory, in SQL.
#' `residuals(fit)` returns the remote design of the fit with two more
#' columns, `.fitted` and `.resid`, which can scale its rows: for
#' heteroskedasticity-consistent and clustered standard errors,
#' ```
#' R <- residuals(fit)
#' meat_hc1 <- crossprod(R, DBI::SQL('".resid" * ".resid"') * R)
#' scores <- colSums(".resid" * R, group = "cluster")
#' meat_cl <- crossprod(scores)
#' ```
#' (with `bread <- vcov(fit) / sigma(fit)^2`, the HC1 covariance is
#' `n / (n - k) * bread %*% meat_hc1 %*% bread`).
#'
#' **Differences from in-memory fits.**
#' * `model.frame()` is not available for a remote fit, and `model` is
#'   ignored: the data stay in the database (`fit$data$design` is the
#'   remote design of the fit).
#' * Results with one value per row are tables in the database, not
#'   vectors: `fitted()` (`.fitted`, `.resid`), `predict()` (`.fitted`;
#'   `se.fit` as a column `.se.fit` with attributes `df` and
#'   `residual.scale`, instead of a list; `.nonestim` instead of the
#'   attribute `"non-estim"`), and `residuals()`, a remote design.
#' * These tables have no row order: rows are matched to the data through
#'   its columns, which every table keeps.
#' * `compress = "auto"` decides from other evidence than in memory (see
#'   *Compression*); the results are the same either way.
#' @section GLMs:
#' [design_glm()] accepts remote data like [design_lm()]:
#' `design_glm(y ~ x + f, binomial, con, table = "d", weights = w, offset = log(e))`,
#' with `weights` and `offset` column names, `DBI::SQL()` or R expressions
#' over the columns, `offset()` terms in the formula, and a remote design
#' with `y` naming the response column (two columns for successes and
#' failures). The fit is the same code as in memory, on the remote methods
#' of the products: with few design patterns (see *Compression*), the
#' iterations run on the patterns; otherwise each iteration computes its
#' products in the database. `predict(fit, newdata = con,
#' table = ..., type = "response")`, `fitted()` and `residuals(fit, type =
#' )` return a temporary view with a `.fitted` (and `.resid`) column next to
#' the columns of the remote data. Fixed effects (`y ~ x | fe`, see
#' [design_glm()]) work remotely as in memory: each iteration's sums over
#' rows are queries, the working weights and responses a temporary table.
NULL

## The arguments of remote data that a fit on a formula takes in `...`.
dm_remote_keys <- c("table", "sql", "engine", "chunk_size")

## The arguments in `...` of a fit on a formula, as a named list: those of
## remote data (dm_remote_keys) evaluated, the others NULL, unevaluated,
## for the check of unknown arguments (design_problem(); so that
## design_glm(subset = x > 0) reports `subset`, not that `x` is not found).
dm_dots <- function(...) {
    n <- ...length()
    nm <- ...names() %||% rep("", n)
    out <- vector("list", n)
    names(out) <- nm
    for (i in which(nm %in% dm_remote_keys)) out[i] <- list(...elt(i))
    out
}

## The input stage of the fits (design_lm(), design_glm() and fixed
## effects): the call's arguments as a design and the columns of the
## response, weights and offset, with one path per source of the design:
##  * a formula and a data frame: a model frame, as lm() and glm() build it
##    (dm_problem_memory());
##  * a formula and remote data: a remote design, dm_remote_fit()
##    (dm_problem_remote());
##  * a design matrix, with `y`, `weights` and `offset` vectors
##    (dm_problem_design());
##  * a remote design, with `y`, `weights` and `offset` columns or
##    expressions over its data (dm_problem_rdesign()).
## The fits then run the same code on any of them.
##   call: the matched call of the fit (its arguments are taken from it,
##     unevaluated), env its parent frame; formula and data their values
##     (formula: a formula or a design; data NULL when missing); dots the
##     arguments in `...`.
##   family: NULL for a linear model; for a GLM, its dmglm_kind(): the
##     response as glm() takes it (a factor, a two-column response), and
##     its columns (dmglm_cols()) as `cols`.
##   fe, clv: with fixed effects, the parsed formula (dm_fe_parse()) and
##     the cluster variable (or NULL): the design is that of the main terms,
##     rows with a missing FE or cluster value are dropped, and `keys` holds
##     their values (dm_fe_keys()).
##   frame: the fit's own environment (its environment()), from which the
##     values of `y`, `weights` and `offset` given with a design matrix
##     (and a GLM's `y` with a remote design) are taken.
## The response, weights, offset and subset of a formula are evaluated in
## the data and the formula's environment (as model.frame() does),
## na.action in env.
## Messages are the same in memory and remotely: for the offset length,
## those of lm() for a linear model without FE, of glm() otherwise.
## Returns list(x (the design), y (the response), w, off (or NULL), cols
## (GLM), resp (remotely: the columns of a two-column response), keys,
## frame (the model frame in memory), na.action, exclude (remotely,
## na.exclude), rows, nall, terms, formula, contrasts, xlevels, env (the
## formula's environment), intercept (whether the design has one, for the
## fit), offset_call (the `offset` argument as given, for new data), design
## (whether the first argument was a design)); see dm_na_result() for the
## na.action of the results.
design_problem <- function(call, env, formula, data, contrasts, dots, family = NULL, fe = NULL, clv = NULL,
                           frame = NULL) {
    P <- if (inherits(formula, "design_matrix")) {
        dm_problem_design(dm_frame_arg(call, frame, "y"), dm_frame_arg(call, frame, "weights"),
                          dm_frame_arg(call, frame, "offset"), formula, family)
    } else if (inherits(formula, "dm_remote")) {
        dm_problem_rdesign(formula, if (is.null(family)) call[["y"]] else dm_frame_arg(call, frame, "y"),
                           call[["weights"]], call[["offset"]], env, family)
    } else if (!is.null(data) && dm_is_remote_source(data)) {
        dm_problem_remote(call, env, formula, data, contrasts, dots, family, fe, clv)
    } else {
        if (length(dots)) stop("Unknown arguments: ", paste(names(dots), collapse = ", "), call. = FALSE)
        dm_problem_memory(call, env, contrasts, family, fe, clv)
    }
    P$design <- inherits(formula, c("design_matrix", "dm_remote"))
    if (!P$design) P$offset_call <- call$offset
    if (!is.null(fe)) P$keys <- dm_fe_keys(P$x, fe, clv, P$frame)
    if (!is.null(family)) {
        P$cols <- if (inherits(P$x, "dm_remote")) {
            dmglm_cols(P$y, P$resp$a, P$resp$b, if (is.null(P$w)) dmc_const(P$x, 1) else P$w, P$off, family,
                       P$resp$cbind, weighted = !is.null(P$w))
        } else {
            dmglm_mem_cols(P$y, P$w, P$off, family)
        }
        ## (on a design no rows are dropped: missing values are errors)
        if (P$design) dmglm_check_design(P$cols)
    }
    P
}

## The value of argument `nm` of the fit whose environment is `frame` (NULL
## when it was not given): the vectors given with a design matrix, as the
## fit received them.
dm_frame_arg <- function(call, frame, nm) {
    if (is.null(call[[nm]])) return(NULL)
    get(nm, envir = frame)
}

## A design matrix and the vectors `y`, `weights` and `offset` (as lm.fit()
## and glm.fit() take them; for a GLM, `y` as glm() takes the response: a
## factor, a two-column matrix).
dm_problem_design <- function(y, weights, offset, x, family) {
    checkmate::assert_class(x, "design_matrix")
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
dm_problem_rdesign <- function(D, y, weights, offset, env, family) {
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
            dm_remote_column(D, y, "y", several = TRUE)
            if (length(y) == 2L) call("cbind", as.name(y[1L]), as.name(y[2L]))
            else if (length(y) == 1L) as.name(y)
            else stop("`y` must name one column, or two (successes and failures).", call. = FALSE)
        } else if (!is.null(D$response)) {
            D$vars[[D$response]]$expr
        } else {
            stop("With a remote design, give `y`, the name of the response column (or use a design with a response).",
                 call. = FALSE)
        }
        resp <- dmglm_remote_response(lhs, D$src$types, dm_con(D))
    }
    intercept <- attr(D$plan$terms, "intercept") > 0L && !isTRUE(D$plan$drop_intercept)
    r <- dmp_resolve(D, c(resp$args, list(weights = weights, offset = offset)), env)
    ## (a linear model names the expressions evaluated in R)
    if (is.null(family) && length(r$evaluated))
        message(dm_untranslated_msg(r$evaluated), ": evaluated in R, in batches of rows.")
    D <- r$x
    dmglm_check_factor(resp, family)
    list(x = D, y = if (length(resp$present)) dmglm_remote_factor(D, resp$present) else r$cols$y,
         w = r$cols$weights, off = r$cols$offset, resp = list(a = r$cols$a, b = r$cols$b, cbind = resp$cbind),
         exclude = FALSE, nall = D$n, contrasts = D$plan$contrasts, xlevels = D$plan$xlevels, env = env,
         intercept = intercept, offset_call = offset)
}

## In memory: the model frame (with the FE and cluster variables, so that
## rows with a missing value in any of them are dropped), its terms (of
## the main terms with FE), response (as lm() takes it; a GLM's as in the
## frame: a vector, factor or matrix), weights, offset (offset() terms and
## the `offset` argument) and the design (with the frame's predvars).
dm_problem_memory <- function(call, env, contrasts, family, fe, clv) {
    lm <- is.null(family) && is.null(fe)
    mf <- call[c(1L, match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(call), 0L))]
    if (!is.null(fe)) {
        aug <- dm_fe_augment(fe$main, unique(c(unlist(fe$fe, use.names = FALSE), clv)))
        environment(aug) <- environment(fe$main)
        mf$formula <- aug
    }
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(design.matrix::design_frame)
    mf <- eval(mf, env)
    mt <- attr(mf, "terms")
    if (!is.null(fe)) mt <- dm_fe_main_terms(mt, fe$main)
    if (attr(attr(mf, "terms"), "response") == 0L) stop("The formula has no response.", call. = FALSE)
    if (lm) {
        y <- dmfit_response(mf)
        if (is.matrix(y)) stop("Multivariate responses are not supported.", call. = FALSE)
    } else {
        y <- mf[[1L]]
        if (!is.null(fe) && (is.numeric(y) || is.logical(y))) y <- dmfit_response(mf)
    }
    w <- as.vector(stats::model.weights(mf))
    if (!is.null(w) && !is.numeric(w)) stop("'weights' must be a numeric vector", call. = FALSE)
    off <- as.vector(stats::model.offset(mf))
    if (!is.null(off)) dm_offset_length(length(off), NROW(y), lm)
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
dm_offset_length <- function(noff, n, lm) {
    if (noff != n)
        stop(gettextf(if (lm) "number of offsets is %d, should equal %d (number of observations)"
                      else "number of offsets is %d should equal %d (number of observations)",
                      noff, n), call. = FALSE)
}

## The terms of the main formula of an FE fit, from those of its model
## frame (which has the FE and cluster variables too).
dm_fe_main_terms <- function(mt_aug, main) {
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
## columns (dm_remote_fit()), on the rows where the model, FE and cluster
## variables, the response, weights and offset are present, as the model
## frame in memory. A GLM's response: dmglm_remote_response() (a factor
## as 0/1, the columns of cbind(a, b) as `resp`).
dm_problem_remote <- function(call, env, formula, source, contrasts, dots, family, fe, clv) {
    ## (with FE, the callers allow na.omit and na.exclude only)
    na <- if (is.null(call$na.action)) "omit" else dm_na_kind(eval(call$na.action, env), what = "For remote data, ")
    if (!is.null(fe)) {
        nms <- names(dots) %||% rep("", length(dots))
        if (any(!nms %in% dm_remote_keys))
            stop("Unknown arguments: ", paste(nms[!nms %in% dm_remote_keys], collapse = ", "), call. = FALSE)
    }
    dm_need_duckdb()
    src <- dm_remote_source(source, dots$table, dots$sql)
    if (!is.null(fe)) formula <- fe$main
    ## (the response, weights, offset and subset are evaluated in the
    ## formula's environment, as model.frame() does in memory)
    fenv <- environment(formula) %||% env
    tt <- stats::terms(formula, data = src$proto)
    if (attr(tt, "response") == 0L) stop("The formula has no response.", call. = FALSE)
    lhs <- attr(tt, "variables")[[attr(tt, "response") + 1L]]
    resp <- if (!is.null(family)) dmglm_remote_response(lhs, src$types, src$holder$con)
            else list(args = list(y = lhs), present = character(), cbind = FALSE)
    ## a linear model has one response, as in memory (a matrix response
    ## evaluated in R gives the same error, dm_remote_fit())
    mresp <- if (is.null(family)) "Multivariate responses are not supported."
    if (!is.null(mresp) && is.call(lhs) && identical(lhs[[1L]], as.name("cbind")) && length(lhs) > 2L)
        stop(mresp, call. = FALSE)
    args <- c(resp$args, list(weights = call$weights, offset = call$offset, subset = if (is.null(family)) call$subset))
    r <- dm_remote_fit(formula, src, args, fenv, dots, contrasts, logical = "subset",
                       present = unique(c(unlist(fe$fe, use.names = FALSE), clv, resp$present)), na = na,
                       mresp = mresp)
    D <- r$x
    ## an offset with several columns: its number of values, as in memory
    if (r$offk > 1L) dm_offset_length(r$offk * D$n, D$n, is.null(family) && is.null(fe))
    dmglm_check_factor(resp, family)
    y <- if (length(resp$present)) dmglm_remote_factor(D, resp$present) else r$cols$y
    list(x = D, y = y, w = r$cols$weights, off = r$off, resp = list(a = r$cols$a, b = r$cols$b, cbind = resp$cbind),
         na.action = NULL, exclude = na == "exclude", rows = NULL, nall = D$n,
         terms = tt, formula = formula, contrasts = D$plan$contrasts,
         xlevels = D$plan$xlevels, env = fenv, intercept = attr(tt, "intercept") > 0L)
}

## The remote design of a fit and the columns of its arguments (one code
## path for design_lm(), design_glm() and fixed effects).
##   formula: the model formula; the design is that of its right-hand side
##     (offset() terms included; their values are columns too), and its
##     response is given by the caller in `args`.
##   src: an opened source (dm_remote_source()).
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
## evaluated in R (a matrix: the caller reports its length, as in memory)).
dm_remote_fit <- function(formula, src, args, env, dots, contrasts = NULL, logical = character(),
                          present = character(), na = "omit", mresp = NULL) {
    miss <- setdiff(present, names(src$types))
    if (length(miss)) stop(sprintf("No column `%s` in the remote data.", miss[1L]), call. = FALSE)
    tt <- stats::terms(formula, data = src$proto)
    tv <- dm_terms_vars(tt)
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
        a <- dm_resolve_or_eval(all[[nm]], env, src$types, con, label(nm), bool = nm %in% logical)
        if (is.null(a$expr)) {
            res[[nm]] <- a
            next
        }
        col <- paste0("__dm_eval", length(exprs) + 1L)
        items[[label(nm)]] <- a
        exprs[[col]] <- a$expr
        if (nm == "offset" || nm %in% names(offt)) mat[[col]] <- "offset"
        if (nm == "y" && !is.null(mresp)) mat[[col]] <- mresp
        res[[nm]] <- list(sql = if (nm %in% logical) sprintf("s.%s", q(col)) else sprintf("CAST(s.%s AS DOUBLE)", q(col)),
                          inputs = col, logical = nm %in% logical)
    }
    if (length(exprs)) {
        dots <- dm_untranslated_engine(items, dots)
        lg <- unlist(lapply(res, function(r) if (isTRUE(r$logical)) r$inputs))
        src <- dm_materialize(src, exprs, env, logical = lg,
                              chunk_size = dots$chunk_size %||% getOption("design.matrix.chunk_size", 250000),
                              matrix = mat)
    }
    offk <- max(c(1L, unlist(attr(src, "ncol"))))
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
    D <- tryCatch(dm_remote(rhs, src, contrasts, NULL, dots, drop_unused = TRUE, subset = sub, extra = extra,
                            override = override),
                  dm_no_rows = function(e) stop("0 (non-NA) cases", call. = FALSE))
    if (na %in% c("fail", "pass")) dm_remote_fit_na(D, res, na, nonlog, names(offt))
    cols <- lapply(names(args), function(nm) if (!is.null(res[[nm]])) dm_rcols(D, res[[nm]]$sql))
    names(cols) <- names(args)
    off <- c(res[["offset"]]$sql, vapply(res[names(offt)], `[[`, "", "sql"))
    list(x = D, cols = cols, off = if (length(off)) dm_rcols(D, paste(sprintf("(%s)", off), collapse = " + ")),
         offk = offk)
}

## na.fail and na.pass: the errors of lm() when rows have missing values.
dm_remote_fit_na <- function(D, res, na, args, offt) {
    total <- as.numeric(dm_query(D, sprintf("SELECT COUNT(*) AS n FROM (%s) AS s", D$src$select))$n)
    if (total == D$n) return(invisible())
    if (na == "fail") stop("missing values in object", call. = FALSE)
    bad <- function(nm) {
        r <- res[[nm]]
        if (is.null(r)) return(FALSE)
        as.numeric(dm_query(D, sprintf("SELECT COUNT(*) AS n FROM (%s) AS s WHERE (%s) IS NULL OR isnan(%s)",
                                       D$src$select, r$sql, r$sql))$n) > 0
    }
    if (any(vapply(c("y", "a", "b", "offset", offt), bad, NA))) stop("NA/NaN/Inf in 'y'", call. = FALSE)
    if (bad("weights")) stop("missing or negative weights not allowed", call. = FALSE)
    stop("NA/NaN/Inf in 'x'", call. = FALSE)
}

## The response of a remote GLM as arguments: list(args, present, cbind).
dmglm_remote_response <- function(lhs, types, con) {
    if (is.call(lhs) && identical(lhs[[1L]], as.name("cbind"))) {
        if (length(lhs) != 3L) stop("A two-column response needs two columns (successes and failures).", call. = FALSE)
        return(list(args = list(a = lhs[[2L]], b = lhs[[3L]]), present = character(), cbind = TRUE))
    }
    if (is.name(lhs) && isTRUE(types[as.character(lhs)] %in% c("character", "factor")))
        return(list(args = list(), present = as.character(lhs), cbind = FALSE))
    if (is.name(lhs) && isTRUE(types[as.character(lhs)] == "logical"))
        lhs <- DBI::SQL(sprintf("CAST(s.%s AS DOUBLE)", DBI::dbQuoteIdentifier(con, as.character(lhs))))
    list(args = list(y = lhs), present = character(), cbind = FALSE)
}

## A character or factor response (dmglm_remote_response()) needs a
## binomial family, as in memory (dmglm_mem_cols()).
dmglm_check_factor <- function(resp, k) {
    if (length(resp$present) && !is.null(k) && !k$binomial)
        stop("A factor response needs a binomial family.", call. = FALSE)
}

## A character or factor response as a column: 0 for the first level of the
## rows used, 1 for the others.
dmglm_remote_factor <- function(D, col) {
    q <- dm_qi(D, col)
    vals <- dm_query(D, sprintf("SELECT DISTINCT s.%s AS v FROM (%s) AS s WHERE %s", q, D$src$select, D$filter))$v
    proto <- D$src$proto[[col]]
    f <- if (is.factor(proto)) factor(as.character(vals), levels = levels(proto)) else factor(vals)
    first <- DBI::dbQuoteString(dm_con(D), levels(droplevels(f))[1L])
    dm_rcols(D, sprintf("CASE WHEN CAST(s.%s AS VARCHAR) = %s THEN 0 ELSE 1 END", q, first))
}

## The offsets of new remote data D for a fit: offset() terms of the plan
## and the `offset` argument of the fit, resolved against the columns of D
## (SQL), or NULL. Returns list(D, off): offsets with no SQL translation are
## evaluated in R into new columns of the data of D (dm_materialize()).
dm_remote_new_offset <- function(object, D) {
    con <- dm_con(D)
    tt <- D$plan$terms
    nms <- dm_terms_vars(tt)$names[attr(tt, "offset")]
    pieces <- character()
    exprs <- list()
    add <- function(e) {
        col <- paste0("__dm_offn", length(exprs) + 1L)
        exprs[[col]] <<- e
        pieces <<- c(pieces, sprintf("CAST(s.%s AS DOUBLE)", DBI::dbQuoteIdentifier(con, col)))
    }
    for (nm in nms) {
        v <- D$vars[[nm]]
        if (v$kind == "numeric") pieces <- c(pieces, v$sql) else add(v$expr[[2L]])
    }
    arg <- object$offset_call
    env <- object$env %||% environment(object$terms) %||% globalenv()
    if (!is.null(arg)) {
        a <- dm_resolve_or_eval(arg, env, D$src$types, con, "offset")
        if (is.null(a$expr)) pieces <- c(pieces, a$sql) else add(a$expr)
    }
    if (length(exprs)) D$src <- dm_materialize(D$src, exprs, env, chunk_size = D$chunk_size)
    list(D = D, off = if (length(pieces)) dm_rcols(D, paste(sprintf("(%s)", pieces), collapse = " + ")))
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
##     exprs), list(inputs, factor, label), see dm_eval_model(). Their
##     columns and kinds are returned in the attribute "model" of the
##     result: list(kind, cols) per expression.
##   matrix: for expressions that may give a matrix with one row per row
##     (named as in exprs): "offset" (its first column is written, and its
##     number of columns returned in the attribute "ncol" of the result, for
##     the caller's error), or the message of the error for a matrix with
##     several columns. A one-column matrix is its column.
## Returns the new source (dm_remote_source() of the temporary table),
## which keeps the original one alive.
dm_materialize <- function(src, exprs, env, logical = character(), chunk_size = 250000, model = list(),
                           matrix = list()) {
    info <- list()
    ncols <- list()
    con <- src$holder$con
    q <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
    from <- src$base_table
    raw <- NULL
    if (is.null(from)) {
        raw <- q(dm_tmp_name("dm_raw"))
        DBI::dbExecute(con, sprintf("CREATE TEMPORARY TABLE %s AS %s", raw, src$select))
        on.exit(DBI::dbExecute(con, sprintf("DROP TABLE IF EXISTS %s", raw)), add = TRUE)
        from <- raw
    }
    r <- DBI::dbGetQuery(con, sprintf("SELECT MIN(rowid) AS lo, MAX(rowid) AS hi FROM %s", from))
    out <- dm_tmp_name("dm_eval")
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
                r <- dm_eval_model(exprs[[nm]], chunk, env, model[[nm]], nm, info[[nm]])
                info[[nm]] <<- r$info
                for (cn in names(r$cols)) chunk[[cn]] <- r$cols[[cn]]
                next
            }
            v <- eval(exprs[[nm]], chunk, env)
            if (length(v) == 1L && nrow(chunk) != 1L) v <- rep(v, nrow(chunk))
            if (!is.null(matrix[[nm]]) && length(dim(v)) == 2L && nrow(v) == nrow(chunk)) {
                if (ncol(v) > 1L) {
                    if (!identical(matrix[[nm]], "offset")) stop(matrix[[nm]], call. = FALSE)
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
    new <- dm_remote_source(con, out, NULL)
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
dm_eval_model <- function(expr, chunk, env, m, nm, prev = NULL) {
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

## Resolve an argument (dm_resolve_arg()), or, when its R expression has no
## SQL translation, return it for evaluation in R: list(sql, inputs) or
## list(expr, reason).
dm_resolve_or_eval <- function(x, env, types, con, what, bool = FALSE) {
    if (is.null(x)) return(NULL)
    if (is.call(x) && !deparse2(x[[1L]]) %in% c("SQL", "DBI::SQL")) {
        tr <- if (bool) dm_sql_bool(x, types, con) else dm_sql_num(x, types, con)
        if (is.null(tr$sql)) return(list(expr = x, reason = tr$reason))
    }
    dm_resolve_arg(x, env, types, con, what, bool)
}

## The expressions of a fit that SQL cannot compute, as one message:
## "weights `myfun(x)` cannot be translated to SQL (myfun() has no SQL
## translation)".
dm_untranslated_msg <- function(items) {
    nm <- names(items)[1L]
    sprintf("%s `%s` cannot be translated to SQL (%s)", nm, deparse2(items[[1L]]$expr), items[[1L]]$reason)
}

## A fit with expressions evaluated in R: a message names the first one;
## `engine = "sql"` refuses.
dm_untranslated_engine <- function(items, dots) {
    if (!length(items)) return(dots)
    why <- dm_untranslated_msg(items)
    if (identical(dots$engine, "sql")) stop("The SQL engine cannot compute this fit: ", why, ".", call. = FALSE)
    message("Evaluated in R: ", why, ".")
    dots
}

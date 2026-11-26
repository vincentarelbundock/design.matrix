#' Compile a formula or a fitted model into a reusable plan
#'
#' A plan records everything needed to build the model matrix of new data
#' the way [stats::predict.lm()] does: the terms (with the parameters of
#' transformations such as [stats::poly()] or [splines::ns()] fixed by the
#' training data), the factor levels, and the contrasts.
#'
#' From a formula and training data, every piece comes from stats functions:
#' [stats::terms()], [stats::model.frame()], [stats::.getXlevels()] and
#' `contrasts<-`. From a fitted model, the pieces are those the model stores,
#' as returned by [design_parts()]; the training data is not needed.
#'
#' Apply a plan with `design_matrix(plan, newdata)`. Its dense matrix,
#' `as.matrix(design_matrix(plan, newdata))`, is identical (in the strict
#' sense described in [design_matrix()]) to
#' ```
#' model.matrix(plan$terms,
#'              model.frame(plan$terms, newdata, xlev = plan$xlevels,
#'                          na.action = na.action),
#'              contrasts.arg = plan$contrasts)
#' ```
#' where `plan$terms` is the training terms with the response removed.
#'
#' When `newdata` is a data frame whose variables have the same types,
#' levels and column names as the training data, the plan skips
#' [stats::model.frame()] and the contrast set-up entirely, which makes
#' small calls much faster. Anything else goes through the general path
#' above, so results never depend on which path was taken. A plan built from
#' a formula prepares this fast path from the training data; a plan built
#' from a fitted model prepares it from the first data it is applied to.
#'
#' @param object A model formula or `terms` object, or a fitted model with a
#'   [design_parts()] method (the default method covers most models with a
#'   [stats::terms()] method).
#' @param data The training data frame, for a formula.
#' @param contrasts.arg,xlev As in [stats::model.matrix()], for a formula.
#' @param ... Passed to [design_parts()] for a fitted model, for example
#'   `component = "precision"` for a `betareg` model.
#' @return A `design_plan` object.
#' @keywords internal
design_plan <- function(object, ...) {
    UseMethod("design_plan")
}

#' @rdname design_plan
#' @export
design_plan.formula <- function(object, data, contrasts.arg = NULL, xlev = NULL, ...) {
    checkmate::assert_data_frame(data)
    checkmate::assert_list(contrasts.arg, null.ok = TRUE)
    checkmate::assert_list(xlev, null.ok = TRUE)

    mf <- design_frame(object, data, xlev = xlev)
    tt <- attr(mf, "terms")
    xlevels <- stats::.getXlevels(tt, mf)
    prep <- dm_prepare(tt, mf, contrasts.arg, NULL, FALSE)
    contr <- if (is.null(prep)) {
        attr(stats::model.matrix(tt, mf, contrasts.arg = contrasts.arg), "contrasts")
    } else {
        prep$contrasts
    }
    dm_new_plan(tt, xlevels, contr, data)
}

#' @rdname design_plan
#' @export
design_plan.default <- function(object, ...) {
    parts <- design_parts(object, ...)
    checkmate::assert_list(parts)
    checkmate::assert_class(parts$terms, "terms", .var.name = "design_parts()$terms")
    checkmate::assert_list(parts$xlevels, null.ok = TRUE, .var.name = "design_parts()$xlevels")
    checkmate::assert_list(parts$contrasts, null.ok = TRUE, .var.name = "design_parts()$contrasts")
    checkmate::assert_flag(parts$drop_intercept, null.ok = TRUE, .var.name = "design_parts()$drop_intercept")
    ## levels of variables that are not in the terms (those of another
    ## component of the model, say) would make model.frame() warn
    xlevels <- parts$xlevels
    if (length(xlevels)) {
        vnames <- vapply(as.list(attr(parts$terms, "variables"))[-1L], deparse2, "")
        xlevels <- xlevels[names(xlevels) %in% vnames]
        if (!length(xlevels)) xlevels <- NULL
    }
    drop_intercept <- isTRUE(parts$drop_intercept)
    key <- list(parts$terms, xlevels, parts$contrasts, drop_intercept)
    plan <- dm_plan_cache_get(key)
    if (is.null(plan)) {
        plan <- dm_new_plan(parts$terms, xlevels, parts$contrasts, NULL,
                            drop_intercept = drop_intercept)
        dm_plan_cache_set(key, plan)
    }
    plan
}

## Plans of fitted models, memoized on the pieces that define them. A plan is
## a function of its terms (with their environment and predvars), levels,
## contrasts and intercept flag only, so two models with identical pieces
## share one plan, and the fast path a plan learns from its first data
## (stored by reference in `plan$lazy`) is reused by later calls, such as
## repeated `design_matrix(fit, newdata = batch)`. Keys are compared with
## identical(), which compares environments by address; the cache keeps the
## most recently used `dm_plan_cache_size` plans.
dm_plan_cache <- new.env(parent = emptyenv())
dm_plan_cache$entries <- list()
dm_plan_cache_size <- 16L

dm_plan_cache_get <- function(key) {
    entries <- dm_plan_cache$entries
    for (i in seq_along(entries)) {
        if (identical(entries[[i]]$key, key)) {
            if (i > 1L) dm_plan_cache$entries <- c(entries[i], entries[-i])
            return(entries[[i]]$plan)
        }
    }
    NULL
}

dm_plan_cache_clear <- function() {
    dm_plan_cache$entries <- list()
    invisible(NULL)
}

dm_plan_cache_set <- function(key, plan) {
    entries <- c(list(list(key = key, plan = plan)), dm_plan_cache$entries)
    if (length(entries) > dm_plan_cache_size) entries <- entries[seq_len(dm_plan_cache_size)]
    dm_plan_cache$entries <- entries
    invisible(plan)
}

## Assemble the plan on first reuse, without evaluating the training formula again.
dm_plan_for_design <- function(design) {
    if (!is.null(design$plan)) return(design$plan)
    if (!is.null(design$cache$plan)) return(design$cache$plan)
    tt <- design$training_terms
    if (is.null(tt))
        stop("Reusing a design requires a design built from a formula or a plan with training metadata.", call. = FALSE)

    frame <- attr(design$vars, "terms")
    if (is.null(frame)) frame <- design$frame_terms
    ## The supplied frame may contain extra variables or use a different order.
    vars <- as.list(attr(tt, "variables"))[-1L]
    frame_vars <- as.list(attr(frame, "variables"))[-1L]
    predvars <- attr(frame, "predvars")
    if (!is.null(predvars)) {
        idx <- match(vars, frame_vars)
        if (!anyNA(idx)) attr(tt, "predvars") <- as.call(c(list(quote(list)), as.list(predvars)[idx + 1L]))
    }
    xlevels <- stats::.getXlevels(tt, design$vars)
    ## model.matrix converts logicals to factors after building the model frame.
    classes <- attr(frame, "dataClasses")
    xlevels <- xlevels[!names(xlevels) %in% names(classes)[classes == "logical"]]
    plan <- dm_new_plan(tt, xlevels, design$contrasts, data = NULL)
    design$cache$plan <- plan
    plan
}

## `data`: the training data, from which the fast path is prepared now, or
## NULL to prepare it from the first data the plan is applied to.
## `drop_intercept`: remove the intercept column after the matrix is built,
## as model.matrix() does for coxph models.
dm_new_plan <- function(tt, xlevels, contrasts, data, drop_intercept = FALSE) {
    ttx <- stats::delete.response(tt)
    ## With the response also on the right-hand side (y ~ y + x),
    ## delete.response() leaves a term with no variables. R's C routine
    ## then reserves a column it never fills, shifting the columns after
    ## it. Refuse rather than reproduce that.
    fac <- attr(ttx, "factors")
    if (length(attr(ttx, "term.labels")) &&
        (!length(fac) || any(colSums(fac != 0) == 0))) {
        stop("The response also appears on the right-hand side of the formula; ",
             "remove it before building a plan.", call. = FALSE)
    }
    plan <- structure(list(terms = ttx, xlevels = xlevels,
                           contrasts = contrasts, fast = NULL, lazy = NULL,
                           drop_intercept = drop_intercept),
                      class = "design_plan")
    if (is.null(data)) {
        plan$lazy <- new.env(parent = emptyenv())
        plan$lazy$fast <- NULL
        plan$lazy$failures <- 0L
    } else {
        ## Skeleton for the fast path: the set-up of the training data on the
        ## plan's own terms, built through the general path.
        plan$fast <- tryCatch(dm_fast_spec(plan, data), error = function(e) NULL)
    }
    plan
}

#' @export
print.design_plan <- function(x, ...) {
    cat("<design_plan>\n")
    cat("terms:   ", deparse2(stats::formula(x$terms)), "\n")
    if (length(x$xlevels))
        cat("factors: ", paste(names(x$xlevels), collapse = ", "), "\n")
    fast <- if (!is.null(x$fast) || !is.null(x$lazy$fast)) {
        "yes"
    } else if (!is.null(x$lazy) && x$lazy$failures < dm_lazy_attempts) {
        "on first use"
    } else {
        "no"
    }
    cat("fast path:", fast, "\n")
    invisible(x)
}

## The general path: exactly the reference computation.
##
## `vals`, when given, holds the variables the fast path already evaluated
## (the value of the plan's predvars on `newdata`). They are handed to
## model.frame() through placeholder predvars, so that no transformation is
## evaluated twice. Everything else model.frame() does (variable names,
## xlev, na.action, row names) is unchanged, because the variable names
## come from the "variables" attribute, which is untouched.
dm_apply_general <- function(plan, newdata, na.action, vals = NULL) {
    tt <- plan$terms
    data <- newdata
    ## (with no variables, as in y ~ 1, there is nothing to reuse, and
    ## paste0() would still return one placeholder name)
    if (length(vals)) {
        ph <- paste0(".dm_var_", seq_along(vals))
        attr(tt, "predvars") <- as.call(c(list(as.name("list")), lapply(ph, as.name)))
        data <- vals
        names(data) <- ph
        attr(data, "row.names") <- .row_names_info(newdata, 0L)
        class(data) <- "data.frame"
    }
    mf <- design_frame(tt, data, xlev = plan$xlevels, na.action = na.action)
    d <- design_matrix(plan$terms, mf, contrasts.arg = plan$contrasts)
    ## the fast path reuses this set-up, so it drops the column too
    if (isTRUE(plan$drop_intercept)) d <- dm_drop_intercept(d)
    d
}

## Remove the intercept column of a design. The contrasts were chosen when
## the set-up was built, with the intercept, so the other columns are those
## of the model with an intercept; the fill and the kernels only read
## `intercept` to decide whether to write the first column.
dm_drop_intercept <- function(d) {
    if (identical(d$path, "dense")) {
        X <- d$dense
        a <- attr(X, "assign")
        if (length(a) && a[1L] == 0L) {
            att <- attributes(X)
            X <- X[, -1L, drop = FALSE]
            attr(X, "assign") <- a[-1L]
            if (!is.null(att$contrasts)) attr(X, "contrasts") <- att$contrasts
            d$dense <- X
        }
        return(d)
    }
    s <- d$setup
    if (identical(s$intercept, 1L)) {
        s$intercept <- 0L
        s$xnames <- s$xnames[-1L]
        s$assign <- s$assign[-1L]
        d$setup <- s
    }
    d
}

## A plan without training data prepares its fast path from the data it is
## applied to. Data that warns may be unusual, so the next data is tried
## instead; after this many data sets that do not allow a fast path, the plan
## stops trying.
dm_lazy_attempts <- 3L

dm_apply_plan <- function(plan, newdata, na.action) {
    eligible <- is.data.frame(newdata) &&
        (identical(na.action, stats::na.pass) || identical(na.action, "na.pass"))
    if (!eligible) return(dm_apply_general(plan, newdata, na.action))
    spec <- if (is.null(plan$lazy)) plan$fast else plan$lazy$fast
    if (!is.null(spec)) {
        d <- dm_apply_fast(plan, spec, newdata)
        if (inherits(d, "design_matrix")) return(d)
        if (!is.null(d)) return(dm_apply_general(plan, newdata, na.action, vals = d))
    } else if (!is.null(plan$lazy) && plan$lazy$failures < dm_lazy_attempts) {
        g <- dm_apply_general_checked(plan, newdata)
        if (!g$warned) {
            spec <- tryCatch(dm_fast_spec_from(plan, g$design, g$vals),
                             error = function(e) NULL)
            if (is.null(spec)) {
                plan$lazy$failures <- plan$lazy$failures + 1L
            } else {
                plan$lazy$fast <- spec
            }
        }
        return(g$design)
    }
    dm_apply_general(plan, newdata, na.action)
}

## The general path with na.pass, noting whether it warned (the warnings
## still reach the caller). The variables are evaluated once, here, and
## handed to the general path, so that the fast path can be described from
## the same values without evaluating any transformation again.
dm_apply_general_checked <- function(plan, data) {
    warned <- FALSE
    vals <- NULL
    d <- withCallingHandlers({
        vals <- eval(dm_plan_predvars(plan), data, environment(plan$terms))
        dm_apply_general(plan, data, stats::na.pass, vals = vals)
    }, warning = function(w) {
        warned <<- TRUE
    })
    list(design = d, warned = warned, vals = vals)
}

dm_plan_predvars <- function(plan) {
    pv <- attr(plan$terms, "predvars")
    if (is.null(pv)) pv <- attr(plan$terms, "variables")
    pv
}

## Describe each variable of the training model frame, as the fast path
## must reproduce it. Returns NULL when the fast path cannot be used.
dm_fast_spec <- function(plan, data) {
    g <- dm_apply_general_checked(plan, data)
    if (g$warned) return(NULL)
    dm_fast_spec_from(plan, g$design, g$vals)
}

## `d`: the design that the general path built, without warnings, from
## `vals`, the plan's variables evaluated on the data.
dm_fast_spec_from <- function(plan, d, vals) {
    ttx <- plan$terms
    pv <- dm_plan_predvars(plan)
    exprs <- as.list(pv)[-1L]
    if (!length(exprs)) return(NULL)
    if (!identical(d$path, "general")) return(NULL)
    if (length(vals) != length(exprs) || length(d$vars) != length(exprs)) return(NULL)
    env <- environment(ttx)

    vars <- vector("list", length(exprs))
    for (j in seq_along(exprs)) {
        raw <- vals[[j]]
        prepared <- d$vars[[j]]
        nm <- names(d$vars)[j]
        v <- list(expr = exprs[[j]], name = nm,
                  symbol = if (is.name(exprs[[j]])) as.character(exprs[[j]]) else NULL)
        if (is.factor(prepared)) {
            if (!is.null(attr(raw, "contrasts"))) return(NULL)
            if (is.character(raw)) {
                v$kind <- "character"
            } else if (is.factor(raw)) {
                if (!identical(names(attributes(raw)), c("levels", "class"))) return(NULL)
                v$kind <- if (is.ordered(raw)) "ordered" else "factor"
            } else if (is.logical(raw) && is.null(attributes(raw))) {
                v$kind <- "logical"
            } else {
                return(NULL)
            }
            ## a character variable takes its levels from the data unless
            ## the plan fixes them, so other data may encode it differently
            ## (model.matrix() always gives logicals the levels FALSE, TRUE)
            if (identical(v$kind, "character") && is.null(plan$xlevels[[nm]]))
                return(NULL)
            ## contrasts the plan does not fix follow options("contrasts")
            ## at each call, which a cached set-up would not
            if (is.null(plan$contrasts[[nm]])) return(NULL)
            v$levels <- levels(prepared)
            v$class <- class(prepared)
            v$contrasts <- attr(prepared, "contrasts")
            v$attributes <- list(levels = v$levels, class = v$class,
                                 contrasts = v$contrasts)
            if (identical(v$kind, "logical") && !identical(v$levels, c("FALSE", "TRUE")))
                return(NULL)
        } else if (is.numeric(prepared) && !is.object(prepared) ||
                   is.numeric(prepared) && (inherits(prepared, "AsIs") || is.matrix(prepared))) {
            if (is.factor(raw) || !is.numeric(raw)) return(NULL)
            v$kind <- "numeric"
            v$ncol <- NCOL(prepared)
            v$is_matrix <- is.matrix(prepared)
            v$colnames <- if (is.matrix(prepared)) colnames(prepared) else NULL
        } else {
            return(NULL)
        }
        vars[[j]] <- v
    }
    symbols <- vapply(vars, function(v) !is.null(v$symbol), NA)
    list(vars = vars, setup = d$setup, contrasts = d$contrasts,
         names = names(d$vars), env = env, predvars = pv,
         all_symbols = all(symbols))
}

## Build the list that the general path would hand to the C routine, or
## return NULL to fall back to the general path.
dm_apply_fast <- function(plan, spec, newdata) {
    n <- .row_names_info(newdata, 2L)
    if (n < 1L) return(NULL)
    specs <- spec$vars
    ## Evaluate all the variables at once, as model.frame() does, so that
    ## expressions see the same data mask. Variables that are plain columns
    ## of the data are read directly. (If a later check fails, the general
    ## path evaluates the variables again, as model.frame() would.)
    vals <- if (spec$all_symbols) NULL else eval(spec$predvars, newdata, spec$env)
    ## when a check fails: NULL (nothing evaluated yet), or the evaluated
    ## variables for the general path to reuse
    give_up <- vals
    out <- vector("list", length(specs))
    for (j in seq_along(specs)) {
        v <- specs[[j]]
        x <- if (is.null(vals)) .subset2(newdata, v$symbol) else vals[[j]]
        if (is.null(x)) return(give_up)
        if (NROW(x) != n) return(give_up)
        kind <- v$kind
        if (kind == "numeric") {
            if (!is.numeric(x) || is.factor(x)) return(give_up)
            if (v$is_matrix) {
                if (!is.matrix(x) || ncol(x) != v$ncol ||
                    !identical(colnames(x), v$colnames)) return(give_up)
            } else if (!is.null(dim(x))) {
                return(give_up)
            }
        } else if (kind == "factor" || kind == "ordered") {
            ## same levels and class, nothing else: the model frame would
            ## keep the codes as they are
            if (!identical(names(attributes(x)), c("levels", "class")) ||
                !identical(oldClass(x), v$class) ||
                !identical(attr(x, "levels"), v$levels)) return(give_up)
            attr(x, "contrasts") <- v$contrasts
        } else {
            if (!is.null(attributes(x))) return(give_up)
            if (kind == "character") {
                if (!is.character(x)) return(give_up)
                codes <- match(x, v$levels)
                ## new levels: the general path raises the error
                if (anyNA(codes) && any(is.na(codes) & !is.na(x))) return(give_up)
            } else {
                if (!is.logical(x)) return(give_up)
                codes <- as.integer(x) + 1L
            }
            attributes(codes) <- v$attributes
            x <- codes
        }
        out[[j]] <- x
    }
    names(out) <- spec$names
    attr(out, "row.names") <- .row_names_info(newdata, 0L)
    dm_new_design(plan$terms, out, spec$contrasts, setup = spec$setup, path = "fast")
}

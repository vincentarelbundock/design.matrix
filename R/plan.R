## Plans: a formula or a fitted model compiled for reuse on new data
## (internal: users get the plan of a model with design_matrix(fit, newdata =
## NULL), or reuse a design, design_matrix(D, newdata = batch); see "Reusing
## a design or a fitted model on new data" in ?design_matrix).
##
## A plan records what building the model matrix of new data needs, as
## stats::predict.lm() does it: the terms (with the parameters of
## transformations such as poly() or ns() fixed by the training data), the
## factor levels and the contrasts. From a formula and training data, every
## piece comes from stats functions (terms(), model.frame(), .getXlevels(),
## `contrasts<-`); from a fitted model, the pieces are those of
## design_parts(), and the training data is not needed.
##
## When new data is a data frame whose variables have the same types, levels
## and column names as the training data, the plan skips model.frame() and
## the contrast set-up (the fast path); anything else goes through the
## general path, so results never depend on the path taken. A plan built
## from a formula prepares the fast path from the training data; one built
## from a fitted model, from the first data it is applied to. Plans of
## fitted models are cached (below).
xtx_plan_of <- function(object, ...) {
    if (inherits(object, "formula")) xtx_plan_formula(object, ...) else xtx_plan_model(object, ...)
}

## From a formula (or a terms object) and training data.
xtx_plan_formula <- function(object, data, contrasts.arg = NULL, xlev = NULL, ...) {
    checkmate::assert_data_frame(data)
    checkmate::assert_list(contrasts.arg, null.ok = TRUE)
    checkmate::assert_list(xlev, null.ok = TRUE)

    mf <- design_frame(object, data, xlev = xlev)
    tt <- attr(mf, "terms")
    xlevels <- stats::.getXlevels(tt, mf)
    prep <- xtx_prepare(tt, mf, contrasts.arg, NULL, FALSE)
    contr <- if (is.null(prep)) {
        attr(stats::model.matrix(tt, mf, contrasts.arg = contrasts.arg), "contrasts")
    } else {
        prep$contrasts
    }
    xtx_new_plan(tt, xlevels, contr, data)
}

## From a fitted model: the pieces of design_parts(), checked; cached.
xtx_plan_model <- function(object, ...) {
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
    plan <- xtx_plan_cache_get(key)
    if (is.null(plan)) {
        plan <- xtx_new_plan(parts$terms, xlevels, parts$contrasts, NULL,
                            drop_intercept = drop_intercept)
        xtx_plan_cache_set(key, plan)
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
## most recently used `xtx_plan_cache_size` plans.
xtx_plan_cache <- new.env(parent = emptyenv())
xtx_plan_cache$entries <- list()
xtx_plan_cache_size <- 16L

xtx_plan_cache_get <- function(key) {
    entries <- xtx_plan_cache$entries
    for (i in seq_along(entries)) {
        if (identical(entries[[i]]$key, key)) {
            if (i > 1L) xtx_plan_cache$entries <- c(entries[i], entries[-i])
            return(entries[[i]]$plan)
        }
    }
    NULL
}

xtx_plan_cache_clear <- function() {
    xtx_plan_cache$entries <- list()
    invisible(NULL)
}

xtx_plan_cache_set <- function(key, plan) {
    entries <- c(list(list(key = key, plan = plan)), xtx_plan_cache$entries)
    if (length(entries) > xtx_plan_cache_size) entries <- entries[seq_len(xtx_plan_cache_size)]
    xtx_plan_cache$entries <- entries
    invisible(plan)
}

## Assemble the plan on first reuse, without evaluating the training formula again.
xtx_plan_for_design <- function(design) {
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
    plan <- xtx_new_plan(tt, xlevels, design$contrasts, data = NULL,
                        data_classes = classes)
    design$cache$plan <- plan
    plan
}

## `data`: the training data, from which the fast path is prepared now, or
## NULL to prepare it from the first data the plan is applied to.
## `drop_intercept`: remove the intercept column after the matrix is built,
## as model.matrix() does for coxph models.
## `data_classes`: the classes of the training variables (the "dataClasses"
## attribute of the training model frame's terms), which new data must
## match as in predict.lm(); NULL when unknown (no check).
xtx_new_plan <- function(tt, xlevels, contrasts, data, drop_intercept = FALSE,
                        data_classes = attr(tt, "dataClasses")) {
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
                           drop_intercept = drop_intercept,
                           data_classes = data_classes),
                      class = "xtx_plan")
    if (is.null(data)) {
        plan$lazy <- new.env(parent = emptyenv())
        plan$lazy$fast <- NULL
        plan$lazy$failures <- 0L
    } else {
        ## Skeleton for the fast path: the set-up of the training data on the
        ## plan's own terms, built through the general path.
        plan$fast <- tryCatch(xtx_fast_spec(plan, data), error = function(e) NULL)
    }
    plan
}

#' @export
print.xtx_plan <- function(x, ...) {
    cat("<design_plan>\n")
    cat("terms:   ", deparse2(stats::formula(x$terms)), "\n")
    if (length(x$xlevels))
        cat("factors: ", paste(names(x$xlevels), collapse = ", "), "\n")
    fast <- if (!is.null(x$fast) || !is.null(x$lazy$fast)) {
        "yes"
    } else if (!is.null(x$lazy) && x$lazy$failures < xtx_lazy_attempts) {
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
xtx_apply_general <- function(plan, newdata, na.action, vals = NULL) {
    tt <- plan$terms
    data <- newdata
    ## (with no variables, as in y ~ 1, there is nothing to reuse, and
    ## paste0() would still return one placeholder name)
    if (length(vals)) {
        ph <- paste0(".xtx_var_", seq_along(vals))
        attr(tt, "predvars") <- as.call(c(list(as.name("list")), lapply(ph, as.name)))
        data <- vals
        names(data) <- ph
        attr(data, "row.names") <- .row_names_info(newdata, 0L)
        class(data) <- "data.frame"
    }
    mf <- design_frame(tt, data, xlev = plan$xlevels, na.action = na.action)
    ## variables of another type than in training (a numeric variable
    ## supplied as a factor, say): the error of predict.lm(). The fast path
    ## needs no check: it accepts a variable only with the kind it had in
    ## training, and hands anything else to this path.
    if (!is.null(plan$data_classes)) stats::.checkMFClasses(plan$data_classes, mf)
    d <- design_matrix(plan$terms, mf, contrasts.arg = plan$contrasts)
    ## the fast path reuses this set-up, so it drops the column too
    if (isTRUE(plan$drop_intercept)) d <- xtx_drop_intercept(d)
    d
}

## Remove the intercept column of a design. The contrasts were chosen when
## the set-up was built, with the intercept, so the other columns are those
## of the model with an intercept; the fill and the kernels only read
## `intercept` to decide whether to write the first column.
xtx_drop_intercept <- function(d) {
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
xtx_lazy_attempts <- 3L

xtx_apply_plan <- function(plan, newdata, na.action) {
    eligible <- is.data.frame(newdata) &&
        (identical(na.action, stats::na.pass) || identical(na.action, "na.pass"))
    if (!eligible) return(xtx_apply_general(plan, newdata, na.action))
    spec <- if (is.null(plan$lazy)) plan$fast else plan$lazy$fast
    if (!is.null(spec)) {
        d <- xtx_apply_fast(plan, spec, newdata)
        if (inherits(d, "xtx_matrix")) return(d)
        if (!is.null(d)) return(xtx_apply_general(plan, newdata, na.action, vals = d))
    } else if (!is.null(plan$lazy) && plan$lazy$failures < xtx_lazy_attempts) {
        g <- xtx_apply_general_checked(plan, newdata)
        if (!g$warned) {
            spec <- tryCatch(xtx_fast_spec_from(plan, g$design, g$vals),
                             error = function(e) NULL)
            if (is.null(spec)) {
                plan$lazy$failures <- plan$lazy$failures + 1L
            } else {
                plan$lazy$fast <- spec
            }
        }
        return(g$design)
    }
    xtx_apply_general(plan, newdata, na.action)
}

## The general path with na.pass, noting whether it warned (the warnings
## still reach the caller). The variables are evaluated once, here, and
## handed to the general path, so that the fast path can be described from
## the same values without evaluating any transformation again.
xtx_apply_general_checked <- function(plan, data) {
    warned <- FALSE
    vals <- NULL
    d <- withCallingHandlers({
        vals <- eval(xtx_plan_predvars(plan), data, environment(plan$terms))
        xtx_apply_general(plan, data, stats::na.pass, vals = vals)
    }, warning = function(w) {
        warned <<- TRUE
    })
    list(design = d, warned = warned, vals = vals)
}

xtx_plan_predvars <- function(plan) {
    pv <- attr(plan$terms, "predvars")
    if (is.null(pv)) pv <- attr(plan$terms, "variables")
    pv
}

## Describe each variable of the training model frame, as the fast path
## must reproduce it. Returns NULL when the fast path cannot be used.
xtx_fast_spec <- function(plan, data) {
    g <- xtx_apply_general_checked(plan, data)
    if (g$warned) return(NULL)
    xtx_fast_spec_from(plan, g$design, g$vals)
}

## `d`: the design that the general path built, without warnings, from
## `vals`, the plan's variables evaluated on the data.
xtx_fast_spec_from <- function(plan, d, vals) {
    ttx <- plan$terms
    pv <- xtx_plan_predvars(plan)
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
xtx_apply_fast <- function(plan, spec, newdata) {
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
    xtx_new_design(plan$terms, out, spec$contrasts, setup = spec$setup, path = "fast")
}

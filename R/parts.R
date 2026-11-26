#' Terms, factor levels and contrasts of a fitted model
#'
#' `design_parts()` extracts from a fitted model the pieces that
#' [design_plan()] needs to build the model matrix of new data: the terms, the
#' levels of the factors and the contrasts. Each method reads the pieces that
#' the model's own `predict()` method uses.
#'
#' Methods are provided for models from \pkg{stats} (`lm`, `glm`, and the
#' classes that inherit from them, such as `MASS::rlm()`, `MASS::glm.nb()` and
#' `survey::svyglm()`), \pkg{MASS} (`polr`; `glmmPQL` through `lme`),
#' \pkg{nnet} (`multinom`), \pkg{nlme} (`lme`, `gls`), \pkg{lme4} (`merMod`),
#' \pkg{glmmTMB}, \pkg{ordinal} (`clm`, `clmm`), \pkg{survival} (`coxph`,
#' `survreg`, and the \pkg{survey} classes that inherit from them),
#' \pkg{survey} (`svyolr`), \pkg{betareg}, \pkg{estimatr} (`lm_robust`,
#' `iv_robust`), \pkg{ivreg} and \pkg{AER} (`ivreg`), \pkg{pscl} (`zeroinfl`,
#' `hurdle`) and \pkg{afex} (`afex_aov`).
#'
#' The default method, for other classes, takes the terms from
#' [stats::terms()]; the levels from an `xlevels` element of the model, or else
#' from [stats::.getXlevels()] on [stats::model.frame()]; and the contrasts
#' from a `contrasts` element, or else from the `"contrasts"` attribute of
#' [stats::model.matrix()].
#'
#' Some models have several linear predictors. For these, `component` selects
#' one: `"mean"` or `"precision"` for `betareg`; `"count"` or `"zero"` for
#' `zeroinfl` and `hurdle`; `"cond"`, `"zi"` or `"disp"` for `glmmTMB`. For
#' instrumental-variable models (`ivreg`, `iv_robust`), the parts are those of
#' the regressors. For `coxph`, `strata()` and `cluster()` terms are removed
#' and the intercept column is dropped, as in `model.matrix()` for these
#' models.
#'
#' @section Methods for other classes:
#' `design_parts()` is an S3 generic: a package (or a user) adds support for
#' its own model class by writing a method `design_parts.<class>()` and
#' registering it (with `S3method(design.matrix::design_parts, <class>)` in
#' the NAMESPACE of a package that only suggests \pkg{design.matrix}, or
#' `S3method(design_parts, <class>)` when it imports it). [design_plan()],
#' and therefore `design_matrix(fit, newdata = ...)`, then work for that
#' class. A method returns a list with elements:
#' * `terms`: a `terms` object, with the response or without. Its
#'   `"predvars"` attribute must keep the parameters of data-dependent
#'   transformations such as [stats::poly()], [splines::ns()] or [scale()],
#'   as the `terms` of a model frame do; its environment is where the
#'   variables and functions of the formula are looked up.
#' * `xlevels`: a named list of character vectors, the levels of the
#'   factors (as in the `xlev` argument of [stats::model.frame()]), or
#'   `NULL`. Levels of variables that are not in `terms` are ignored.
#' * `contrasts`: a named list as in the `contrasts.arg` argument of
#'   [stats::model.matrix()] (contrast names or matrices), or `NULL`.
#' * `drop_intercept` (optional): `TRUE` removes the intercept column after
#'   the matrix is built, so that the factors are still coded as in a model
#'   with an intercept (as `model.matrix()` does for `coxph`).
#'
#' The method reads these pieces from the fitted object, as its `predict()`
#' method does; it does not need the training data. Before writing one,
#' check whether the default method already works: it does for most classes
#' that store `terms`, `xlevels` and `contrasts`, or have
#' [stats::terms()], [stats::model.frame()] and [stats::model.matrix()]
#' methods.
#'
#' @param model A fitted model.
#' @param component For models with several linear predictors, which one.
#'   See Details.
#' @param ... Passed to methods.
#' @return A list with elements `terms`, `xlevels` and `contrasts`, and
#'   optionally `drop_intercept` (see the section *Methods for other
#'   classes*).
#' @seealso [design_plan()], which compiles these pieces into a reusable
#'   plan.
#' @examples
#' fit <- lm(mpg ~ poly(hp, 2) + factor(cyl), data = mtcars)
#' parts <- design_parts(fit)
#' names(parts)
#' parts$xlevels
#'
#' # a method for a new class: a fit that stores its formula and levels
#' toy <- structure(list(formula = mpg ~ wt + factor(gear),
#'                       levels = list(`factor(gear)` = c("3", "4", "5"))),
#'                  class = "toyfit")
#' design_parts.toyfit <- function(model, ...) {
#'     list(terms = terms(model$formula), xlevels = model$levels,
#'          contrasts = list(`factor(gear)` = "contr.treatment"))
#' }
#' # (in a package: S3method(design.matrix::design_parts, toyfit))
#' .S3method("design_parts", "toyfit", design_parts.toyfit)
#' p <- design_plan(toy)
#' as.matrix(design_matrix(p, head(mtcars)))
#' @export
design_parts <- function(model, ...) {
    UseMethod("design_parts")
}

## The pieces a model stores, read with `[[` (no partial matching).
dm_parts <- function(terms, xlevels, contrasts, drop_intercept = NULL) {
    if (!inherits(terms, "terms")) {
        stop("The model does not store a terms object.", call. = FALSE)
    }
    out <- list(terms = terms, xlevels = xlevels, contrasts = contrasts)
    if (!is.null(drop_intercept)) out$drop_intercept <- drop_intercept
    out
}

## Some models store the terms of a component without the "predvars" that
## keep the parameters of transformations such as poly(); their predict()
## methods build the model frame from the full terms instead. Copy the
## predvars of the full terms to the variables of the component.
dm_inherit_predvars <- function(tt, full) {
    if (!inherits(tt, "terms") || !is.null(attr(tt, "predvars")) ||
        !inherits(full, "terms") || is.null(attr(full, "predvars"))) return(tt)
    vars <- as.list(attr(tt, "variables"))[-1L]
    fvars <- vapply(as.list(attr(full, "variables"))[-1L], deparse2, "")
    fpred <- as.list(attr(full, "predvars"))[-1L]
    idx <- match(vapply(vars, deparse2, ""), fvars)
    if (anyNA(idx)) {
        stop("Could not match the variables of the model component to those ",
             "of the full model.", call. = FALSE)
    }
    attr(tt, "predvars") <- as.call(c(list(as.name("list")), fpred[idx]))
    tt
}

dm_parts_stored <- function(model) {
    dm_parts(model[["terms"]], model[["xlevels"]], model[["contrasts"]])
}

#' @rdname design_parts
#' @export
design_parts.default <- function(model, ...) {
    tt <- tryCatch(stats::terms(model), error = function(e) NULL)
    if (!inherits(tt, "terms")) {
        stop("`object` must be a formula, a terms object, or a fitted model ",
             "with a terms() method.", call. = FALSE)
    }
    ## `[[` on lists only: no partial matching, and S4 models have no `$`
    stored <- function(nm) if (is.list(model)) model[[nm]] else NULL
    xlevels <- stored("xlevels")
    if (is.null(xlevels)) {
        xlevels <- stats::.getXlevels(tt, stats::model.frame(model))
    }
    contrasts <- stored("contrasts")
    if (is.null(contrasts)) {
        contrasts <- attr(stats::model.matrix(model), "contrasts")
    }
    dm_parts(tt, xlevels, contrasts)
}

## Models that store `terms`, `xlevels` and `contrasts`, which their
## predict() methods use. (`contrasts` is absent when there are no factors.)

#' @rdname design_parts
#' @export
design_parts.lm <- function(model, ...) dm_parts_stored(model)

#' @rdname design_parts
#' @export
design_parts.glm <- function(model, ...) dm_parts_stored(model)

#' @rdname design_parts
#' @export
design_parts.polr <- function(model, ...) dm_parts_stored(model)

#' @rdname design_parts
#' @export
design_parts.multinom <- function(model, ...) dm_parts_stored(model)

#' @rdname design_parts
#' @export
design_parts.svyolr <- function(model, ...) dm_parts_stored(model)

#' @rdname design_parts
#' @export
design_parts.clm <- function(model, ...) dm_parts_stored(model)

#' @rdname design_parts
#' @export
design_parts.clmm <- function(model, ...) dm_parts_stored(model)

#' @rdname design_parts
#' @export
design_parts.survreg <- function(model, ...) dm_parts_stored(model)

#' @rdname design_parts
#' @export
design_parts.lm_robust <- function(model, ...) dm_parts_stored(model)

#' @rdname design_parts
#' @export
design_parts.iv_robust <- function(model, ...) {
    ## terms() folds the instruments into the formula; the regressors have
    ## their own terms, which predict.iv_robust() uses
    tt <- dm_inherit_predvars(model[["terms_regressors"]], model[["terms"]])
    dm_parts(tt, model[["xlevels"]], model[["contrasts"]])
}

#' @rdname design_parts
#' @export
design_parts.ivreg <- function(model, ...) {
    ## as predict.ivreg() in ivreg and AER, except that the parameters of
    ## transformations such as poly() are those of the training data: ivreg
    ## keeps them only in the terms of the stored model frame
    full <- attr(model[["model"]], "terms")
    if (is.null(attr(full, "predvars"))) full <- model[["terms"]][["full"]]
    tt <- dm_inherit_predvars(model[["terms"]][["regressors"]], full)
    dm_parts(tt, model[["levels"]], model[["contrasts"]][["regressors"]])
}

#' @rdname design_parts
#' @export
design_parts.coxph <- function(model, ...) {
    ## as survival's model.matrix.coxph(): cluster() terms, and strata()
    ## terms that do not interact with others, are removed before the matrix
    ## is built, and the intercept column is dropped after
    tt <- model[["terms"]]
    if (!inherits(tt, "terms")) {
        stop("The model does not store a terms object.", call. = FALSE)
    }
    tx <- stats::delete.response(tt)
    specials <- attr(tx, "specials")
    drop <- integer()
    if (length(specials$cluster)) {
        drop <- survival::untangle.specials(tx, "cluster")$terms
    }
    if (length(specials$strata)) {
        st <- survival::untangle.specials(tx, "strata", 1)
        fac <- attr(tx, "factors")
        ord <- attr(tx, "order")
        interacts <- any(vapply(st$vars, function(v) any(ord[fac[v, ] > 0] > 1), NA))
        if (interacts) {
            stop("coxph models whose strata() terms interact with other ",
                 "terms are not supported.", call. = FALSE)
        }
        drop <- c(drop, st$terms)
    }
    if (length(drop)) tt <- tt[-drop]
    dm_parts(tt, model[["xlevels"]], model[["contrasts"]], drop_intercept = TRUE)
}

#' @rdname design_parts
#' @export
design_parts.lme <- function(model, ...) {
    ## as nlme's predict.lme() and predict.gls(): the levels are the row
    ## names of the stored contrast matrices (nlme stores no contrasts for
    ## character variables, so new data cannot drop their levels)
    contrasts <- model[["contrasts"]]
    dm_parts(model[["terms"]], lapply(contrasts, rownames), contrasts)
}

#' @rdname design_parts
#' @export
design_parts.gls <- design_parts.lme

#' @rdname design_parts
#' @export
design_parts.merMod <- function(model, ...) {
    ## lme4 stores neither levels nor contrasts in the model object: the
    ## levels come from the model frame, the contrasts from the fixed-effects
    ## model matrix
    tt <- stats::terms(model, fixed.only = TRUE)
    dm_parts(tt, stats::.getXlevels(tt, stats::model.frame(model)),
             attr(stats::model.matrix(model), "contrasts"))
}

#' @rdname design_parts
#' @export
design_parts.glmmTMB <- function(model, component = c("cond", "zi", "disp"), ...) {
    component <- match.arg(component)
    info <- model[["modelInfo"]]
    tt <- info[["terms"]][[component]][["fixed"]]
    if (!inherits(tt, "terms")) {
        stop(sprintf("The model has no fixed effects for component \"%s\".", component),
             call. = FALSE)
    }
    ## glmmTMB stores only the contrasts the user supplied; the others are
    ## those model.matrix() chose when the model was fitted
    mf <- stats::model.frame(model)
    contrasts <- info[["contrasts"]]
    for (v in vapply(as.list(attr(tt, "variables"))[-1L], deparse2, "")) {
        x <- mf[[v]]
        if (!is.null(contrasts[[v]]) || !(is.factor(x) || is.character(x) || is.logical(x)))
            next
        contrasts[[v]] <- if (!is.null(attr(x, "contrasts"))) attr(x, "contrasts")
                          else getOption("contrasts")[[if (is.ordered(x)) 2L else 1L]]
    }
    dm_parts(tt, stats::.getXlevels(tt, mf), contrasts)
}

#' @rdname design_parts
#' @export
design_parts.betareg <- function(model, component = c("mean", "precision"), ...) {
    component <- match.arg(component)
    dm_parts(model[["terms"]][[component]], model[["levels"]][[component]],
             model[["contrasts"]][[component]])
}

#' @rdname design_parts
#' @export
design_parts.zeroinfl <- function(model, component = c("count", "zero"), ...) {
    ## as predict.zeroinfl() and predict.hurdle(): one set of levels, for
    ## all the variables of the model
    component <- match.arg(component)
    tt <- dm_inherit_predvars(model[["terms"]][[component]], model[["terms"]][["full"]])
    dm_parts(tt, model[["levels"]], model[["contrasts"]][[component]])
}

#' @rdname design_parts
#' @export
design_parts.hurdle <- design_parts.zeroinfl

#' @rdname design_parts
#' @export
design_parts.afex_aov <- function(model, ...) {
    design_parts(model[["lm"]], ...)
}

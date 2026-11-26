## Argument checks of design_lm() and design_glm(), on every path (a design
## matrix or remote design, a formula, a formula with fixed effects), in one
## order:
##
##  1. how the fit is called: with a design, no `data`, `subset`,
##     `na.action`, `contrasts` or unknown arguments; with a formula, no
##     `y` (xtx_check_call());
##  2. with a formula, `na.action` (xtx_check_na());
##  3. the options of the fit (xtx_check_fit_args()): `vcov`, `compress`,
##     what FE do not support, `tol`, `model`.
##
## The checks of the data (response, weights, offset: xtx_fit_check(),
## xtx_glm_checks()) are not here: they are sums over the rows, computed in
## the same pass (remotely, the same query) as the sums the fit needs.

## Step 1. `design`: whether the first argument is a design; `given`: which
## arguments that apply only to a formula were given (a named logical);
## `y` the response argument; `...` the unknown arguments, unevaluated.
xtx_check_call <- function(design, given, y, ...) {
    if (!design) {
        if (!is.null(y)) stop("`y` applies only when the first argument is a design matrix.", call. = FALSE)
        return(invisible())
    }
    if (any(given))
        stop("With a design matrix, give the response as `y`; ",
             xtx_and(sprintf("`%s`", names(given))), " do not apply.", call. = FALSE)
    if (...length()) stop("Unknown arguments: ", paste(...names(), collapse = ", "), call. = FALSE)
    invisible()
}

## "a, b and c"
xtx_and <- function(x) {
    if (length(x) < 2L) return(x)
    paste(paste(x[-length(x)], collapse = ", "), "and", x[length(x)])
}

## Step 2. `na.action`, given with a formula: the na.action functions that
## remote data can apply, refused otherwise in memory too (na_action.R):
## with fixed effects, na.omit or na.exclude; otherwise na.omit, na.exclude,
## na.fail and, for a linear model (`glm = FALSE`), na.pass.
xtx_check_na <- function(na.action, fe, glm = FALSE) {
    if (fe) xtx_fe_na_kind(na.action) else xtx_na_kind(na.action, c("omit", "exclude", "fail", if (!glm) "pass"))
    invisible()
}

## Step 3. `tol`: NULL for a GLM (its tolerance comes from `control`);
## `fe`: whether the formula has fixed effects; `subset`: whether `subset`
## was given and `start` its value (neither is supported with FE).
xtx_check_fit_args <- function(compress, model, tol = NULL, vcov = "iid", fe = FALSE,
                              subset = FALSE, start = NULL) {
    if (!fe && !identical(vcov, "iid"))
        stop("`vcov` other than \"iid\" requires fixed effects (y ~ x | fe).", call. = FALSE)
    xtx_check_compress(compress)
    if (fe) {
        if (isTRUE(compress)) stop("`compress = TRUE` is not supported with fixed effects.", call. = FALSE)
        if (subset) stop("`subset` is not supported with fixed effects.", call. = FALSE)
        if (!is.null(start)) stop("`start` is not supported with fixed effects.", call. = FALSE)
    }
    if (!is.null(tol)) checkmate::assert_number(tol, lower = 0, finite = TRUE)
    checkmate::assert_flag(model)
    invisible()
}

xtx_check_compress <- function(compress) {
    if (!(isTRUE(compress) || isFALSE(compress) || identical(compress, "auto")))
        stop("`compress` must be \"auto\", TRUE or FALSE.", call. = FALSE)
}

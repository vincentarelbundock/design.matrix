#' Build a design without forming the model matrix
#'
#' `design_matrix()` takes the same arguments as [stats::model.matrix()] and
#' returns a `design_matrix` object: the prepared variables plus everything
#' needed to build the matrix (column names, contrasts, term layout). The
#' dense matrix is built by [as.matrix()], and algebra such as
#' `design %*% beta` runs directly on the variables.
#'
#' @param object A model formula or `terms` object, a reusable design, or a fitted model.
#' @param data A data frame, as in [stats::model.matrix()]. For a design or plan, the
#'   new data to apply it to. A DuckDB connection (with `table` or `sql` in
#'   `...`) or Parquet paths give a remote design: see [design_remote].
#' @param newdata A new batch of observations. For a fitted model,
#'   `design_matrix(fit, newdata = batch)` uses its training encoding.
#'   Use `newdata = NULL` to prepare a specification for later batches.
#'   A formula specification learns its encoding from training data when
#'   first applied. Reuse that resulting design for subsequent batches.
#'   Omitting both data arguments for a fitted model returns its original matrix.
#' @param contrasts.arg,xlev As in [stats::model.matrix()].
#' @param na.action Only used when reusing a design or plan: how the model frame
#'   treats missing values. Defaults to [stats::na.pass()], so the output
#'   has one row per row of `data`.
#' @param ... With a fitted model and `newdata`, arguments selecting the model
#'   component, such as `component = "precision"` for a beta regression.
#'   Without `newdata`, passed to [stats::model.matrix()] for non-formula objects;
#'   otherwise ignored, as in `model.matrix.default()`. For remote data,
#'   `table`, `sql`, `engine` and `chunk_size` (see [design_remote]).
#'
#' @return `design_matrix()` returns a design matrix object or, with explicit
#'   `newdata = NULL`, a reusable specification.
#'
#' @details
#' For a formula, `as.matrix(design_matrix(f, d))` is identical to
#' `stats::model.matrix(f, d)` in the strict sense of
#' `identical(a, b, num.eq = FALSE, single.NA = FALSE, attrib.as.set = FALSE)`:
#' same values bit for bit (including signed zeros and the payloads of `NA`
#' and `NaN`), same names and attributes in the same order. The C routine is
#' a port of the one in R's stats package; its fill visits rows in chunks,
#' optionally on several threads (see [design_threads()]), and performs exactly
#' the same floating point operations as R for every cell.
#'
#' **Options.**
#' * `design.matrix.skip_na_copy` (default `TRUE`): build model frames with
#'   [design_frame()], which returns the model frame of
#'   [stats::model.frame()] without the copy of the whole frame that
#'   [stats::na.omit()] and [stats::na.exclude()] make when they remove
#'   nothing, and without recoding factors that already have the levels
#'   asked for in `xlev`. The result is identical either way; set `FALSE`
#'   to build the frame with `model.frame()` itself, for example in
#'   benchmarks against functions that call it.
#' * `design.matrix.structured` (default `TRUE`): use the column-structured
#'   kernels for products with finite designs (see [design_algebra]).
#' * `design.matrix.dedup`: see [design_algebra].
#'
#' @examples
#' d <- data.frame(y = rnorm(10), x = rnorm(10), f = gl(2, 5))
#' X <- as.matrix(design_matrix(y ~ x * f, d))
#' identical(X, model.matrix(y ~ x * f, d))
#'
#' D <- design_matrix(y ~ x * f, d)
#' D %*% c(1, 2, 3, 4)
#' D_new <- design_matrix(D, d[1:3, ])
#' as.matrix(D_new)
#' @export
design_matrix <- function(object, data = environment(object),
                          contrasts.arg = NULL, xlev = NULL,
                          na.action = stats::na.pass, ..., newdata = NULL) {
    checkmate::qassert(contrasts.arg, c("0", "l"))
    checkmate::qassert(xlev, c("0", "l"))

    ## remote data (see remote.R)
    if (inherits(object, "dm_remote")) object <- object$plan
    if (!missing(data) && dm_is_remote_source(data))
        return(dm_remote(object, data, contrasts.arg, xlev, list(...)))
    if (!missing(newdata) && dm_is_remote_source(newdata)) {
        if (!missing(data)) stop("Supply only one of `data` and `newdata`.", call. = FALSE)
        return(dm_remote(object, newdata, contrasts.arg, xlev, list(...)))
    }

    if (!missing(newdata)) {
        if (!missing(data)) stop("Supply only one of `data` and `newdata`.", call. = FALSE)
        if (inherits(object, "formula") && is.null(newdata)) {
            return(structure(list(formula = object, contrasts = contrasts.arg, xlevels = xlev),
                             class = "design_formula"))
        }
        if (!inherits(object, c("formula", "design_matrix", "design_formula", "design_plan")))
            object <- design_plan(object, ...)
        if (is.null(newdata)) {
            if (inherits(object, "design_matrix")) return(dm_plan_for_design(object))
            return(object)
        }
        data <- newdata
    }
    if (inherits(object, "design_formula")) {
        if (missing(data)) stop("Supply training data to apply a formula specification.", call. = FALSE)
        return(design_matrix(object$formula, data, contrasts.arg = object$contrasts, xlev = object$xlevels))
    }
    if (inherits(object, "design_matrix")) {
        if (missing(data)) stop("`data` is required when reusing a design.", call. = FALSE)
        object <- dm_plan_for_design(object)
    }
    if (inherits(object, "design_plan")) {
        if (missing(data)) {
            stop("`data` is required when `object` is a plan.", call. = FALSE)
        }
        design <- dm_apply_plan(object, data, na.action)
        design$plan <- object
        return(design)
    }

    ## Objects that stats would not send to model.matrix.default() (fitted
    ## models, formula subclasses with their own method) go to stats, with
    ## exactly the arguments the caller supplied.
    supplied <- list()
    if (!missing(data)) supplied["data"] <- list(data)
    if (!missing(contrasts.arg)) supplied["contrasts.arg"] <- list(contrasts.arg)
    if (!missing(xlev)) supplied["xlev"] <- list(xlev)
    if (!dm_uses_default_method(object, parent.frame())) {
        return(dm_dense_design(dm_call_stats(object, supplied, list(...), parent.frame())))
    }

    prep <- dm_prepare(object, data, contrasts.arg, xlev, missing(data))
    if (is.null(prep)) {
        ## This version of R's model.matrix.default() is not recognised:
        ## use it directly.
        return(dm_dense_design(dm_call_stats(object, supplied, list(...), parent.frame())))
    }
    design <- dm_new_design(prep$t, prep$data, prep$contrasts, rho = prep$rho)
    design$training_terms <- prep$t
    ## Subsetting a supplied model frame can discard its trained predvars.
    if (is.null(attr(prep$data, "terms")) && !missing(data))
        design$frame_terms <- attr(data, "terms")
    design
}

## Would stats::model.matrix(), called from `env`, dispatch to its default
## method? UseMethod() looks for a method in the calling environment, then
## in the S3 methods table of stats; so do we (getS3method() is ~10x slower).
dm_uses_default_method <- function(object, env) {
    if (!inherits(object, c("formula", "terms"))) return(FALSE)
    tbl <- asNamespace("stats")[[".__S3MethodsTable__."]]
    for (cl in class(object)) {
        nm <- paste0("model.matrix.", cl)
        if (!is.null(get0(nm, envir = env, mode = "function")) ||
            !is.null(get0(nm, envir = tbl, inherits = FALSE)))
            return(FALSE)
    }
    TRUE
}

## Call stats::model.matrix() with the values design_matrix() received.
## `supplied` holds only the arguments the caller gave, so omitted ones stay
## omitted, and nothing the caller wrote is evaluated a second time: the
## values are bound in a fresh environment whose parent is the caller's.
dm_call_stats <- function(object, supplied, dots, env) {
    e <- new.env(parent = env)
    e$object <- object
    args <- list(quote(object))
    for (nm in names(supplied)) {
        assign(nm, supplied[[nm]], envir = e)
        args[[nm]] <- as.name(nm)
    }
    if (length(dots)) {
        e$.dm_dots <- dots
        nms <- names(dots)
        if (is.null(nms)) nms <- rep("", length(dots))
        extra <- lapply(seq_along(dots), function(k) call("[[", quote(.dm_dots), k))
        names(extra) <- nms
        args <- c(args, extra)
    }
    eval(as.call(c(list(quote(stats::model.matrix)), args)), e)
}

deparse2 <- function(x) {
    paste(deparse(x, width.cutoff = 500L, backtick = !is.symbol(x) && is.language(x)),
          collapse = " ")
}

## `vars` is the list given to R's C routine. dm_setup() runs every check,
## contrast evaluation and warning of that routine, in the same order.
## `rho` is the frame of the prepared model.matrix.default() call, where R's
## C routine evaluates contrasts().
dm_new_design <- function(t, vars, contrasts, setup = NULL, path = "general",
                          rho = environment()) {
    if (is.null(setup)) {
        setup <- .Call(C_dm_setup, t, vars, rho)
    }
    out <- list(setup = setup, vars = vars, contrasts = contrasts,
                n = NROW(.subset2(vars, 1L)), dense = NULL, path = path,
                cache = new.env(parent = emptyenv()))
    class(out) <- "design_matrix"
    out
}

dm_dense_design <- function(X) {
    out <- list(setup = NULL, vars = NULL, contrasts = attr(X, "contrasts"),
                n = nrow(X), dense = X, path = "dense")
    out["dense"] <- list(X) # keep the element even when X is NULL
    class(out) <- "design_matrix"
    out
}

#' @export
as.matrix.design_matrix <- function(x, ...) {
    if (dm_is_comb(x)) return(dm_comb_dense(x))
    if (identical(x$path, "dense")) return(x$dense)
    ans <- .Call(C_dm_fill, x$setup, x$vars, design_threads())
    if (!is.null(x$contrasts)) attr(ans, "contrasts") <- x$contrasts
    ans
}

#' @export
dim.design_matrix <- function(x) {
    if (dm_is_comb(x)) return(dim(x$design))
    if (identical(x$path, "dense")) return(dim(x$dense))
    c(x$n, length(x$setup$xnames))
}

#' @export
dimnames.design_matrix <- function(x) {
    if (dm_is_comb(x)) return(dimnames(x$design))
    if (identical(x$path, "dense")) return(dimnames(x$dense))
    rn <- attr(x$vars, "row.names")
    list(if (is.null(rn)) NULL else as.character(rn), x$setup$xnames)
}

#' @export
print.design_matrix <- function(x, ...) {
    d <- dim(x)
    cat(sprintf("<design_matrix: %d rows x %d columns>\n", d[1L], d[2L]))
    if (dm_is_comb(x)) {
        cat(if (is.null(x$design2)) "scaled rows of one design\n" else "combination of two designs\n")
    }
    cn <- dimnames(x)[[2L]]
    if (length(cn)) {
        shown <- utils::head(cn, 8L)
        more <- if (length(cn) > 8L) sprintf(" ... (%d more)", length(cn) - 8L) else ""
        cat("columns:", paste(shown, collapse = ", "), more, "\n")
    }
    invisible(x)
}

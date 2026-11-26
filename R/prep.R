## The R-level preparation of stats:::model.matrix.default, taken from the
## installed version of R at load time.
##
## model.matrix.default() builds the model frame, turns character columns
## into factors and sets contrasts, then calls its C routine:
##
##     ans <- .External2(C_modelmatrix, t, data)
##     if (any(isF)) attr(ans, "contrasts") <- lapply(data[isF], attr, "contrasts")
##     ans
##
## This code changes between R versions (for example the handling of
## `contrasts.arg` changed after R 4.5). Rather than copy one version, we cut
## the installed function just before the C call and make it return what the
## call would have received. The result keeps stats' namespace as its
## environment, so every name in it resolves exactly as in stats.
##
## If the installed function does not end as expected, `prep` is NULL and
## design_matrix() falls back to stats::model.matrix().

the <- new.env(parent = emptyenv())

xtx_make_prep <- function(f = get("model.matrix.default", envir = asNamespace("stats"))) {
    b <- body(f)
    if (!is.call(b) || !identical(b[[1L]], as.name("{"))) return(NULL)
    stmts <- as.list(b)
    call_stmt <- quote(ans <- .External2(C_modelmatrix, t, data))
    tail_expected <- list(
        quote(if (any(isF)) attr(ans, "contrasts") <- lapply(data[isF], attr, "contrasts")),
        quote(ans))
    idx <- which(vapply(stmts, function(e) identical(e, call_stmt), NA))
    if (length(idx) != 1L) return(NULL)
    tail_actual <- stmts[-seq_len(idx)]
    if (!identical(lapply(tail_actual, strip_srcref), tail_expected)) return(NULL)
    new <- stmts[seq_len(idx - 1L)]
    new[[idx]] <- quote(return(list(t = t, data = data, isF = isF, rho = environment())))
    new <- lapply(new, xtx_swap_model_frame)
    body(f) <- as.call(new)
    ## Every name resolves as in stats, except design_frame().
    env <- new.env(parent = environment(f))
    env$design_frame <- design_frame
    environment(f) <- env
    compiler::cmpfun(f)
}

## Replace the call `model.frame(object, data, xlev = xlev)` of
## model.matrix.default() by design_frame() with the same arguments (see
## frame.R: the same model frame, without redundant copies).
xtx_swap_model_frame <- function(e) {
    if (is.call(e) && identical(strip_srcref(e), quote(model.frame(object, data, xlev = xlev))))
        return(quote(design_frame(object, data, xlev = xlev)))
    if (is.call(e)) {
        for (i in seq_along(e)[-1L]) {
            ei <- e[[i]]
            if (!missing(ei) && !is.null(ei)) e[[i]] <- xtx_swap_model_frame(ei)
        }
    }
    e
}

strip_srcref <- function(e) {
    attributes(e) <- NULL
    e
}

.onLoad <- function(libname, pkgname) {
    the$prep <- tryCatch(xtx_make_prep(), error = function(e) NULL)
    the$frame <- tryCatch(xtx_init_frame(), error = function(e) NULL)
}

## Returns list(t, data, isF, rho, contrasts), or NULL when the preparation
## code of this version of R is not recognised.
xtx_prepare <- function(object, data, contrasts.arg, xlev, data_missing) {
    prep <- the$prep
    if (is.null(prep)) return(NULL)
    out <- if (data_missing) {
        prep(object, contrasts.arg = contrasts.arg, xlev = xlev)
    } else {
        prep(object, data, contrasts.arg = contrasts.arg, xlev = xlev)
    }
    out$contrasts <- if (any(out$isF)) lapply(out$data[out$isF], attr, "contrasts") else NULL
    out
}

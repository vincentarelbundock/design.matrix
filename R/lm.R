## Linear models fitted on a design_matrix, without the dense matrix.
##
## The code is split so that a later glm() fit can reuse it:
##
##  * dmfit_frame():     model frame, terms, response, weights and offset,
##                       as lm() and glm() extract them;
##  * dmfit_check():     checks of the response, weights and offset;
##  * dmfit_groups():    rows that share a design pattern (the dedup
##                       machinery of algebra.R);
##  * dmfit_suffstats(): per-pattern sufficient statistics of the response;
##  * dmfit_ls():        a weighted least squares fit with lm()'s pivoting,
##                       on rows or on patterns;
##  * dmfit_vcov():      the covariance matrix with NA for aliased columns;
##  * dmfit_spec():      the design specification kept for predictions.

#' Linear models without the model matrix
#'
#' `design_lm()` fits a linear model by least squares, like [stats::lm()],
#' without forming the dense model matrix. Given a [design_matrix()] and a
#' response `y` instead of a formula and data, it fits that design, like
#' [stats::lm.fit()].
#'
#' **Method.** The fit uses [qr()][design_qr] on the (weighted) design
#' matrix: `R` is the Cholesky factor of `X'WX`, computed from factor codes
#' and numeric variables; columns are pivoted, and aliased columns detected,
#' as in [lm()] (same tolerance, same `NA` pattern in the coefficients); the
#' coefficients are refined from the residuals. Results agree with those of
#' `lm()` up to floating point rounding.
#'
#' **Compression.** When many rows share a design pattern (the same factor
#' levels and numeric values in every variable of the design), the fit only
#' needs, for each pattern, the sum of the weights, the weighted mean of the
#' response (minus the offset) and the weighted sum of squares around that
#' mean. The coefficients, residual sum of squares and covariance matrix
#' computed from these sufficient statistics are those of the fit on every
#' row, so a fit on millions of rows of categorical data costs about as
#' much as one on its distinct patterns. Patterns are found by the
#' deduplication machinery of the products (see [design_algebra] and the
#' option `design.matrix.dedup`). With `compress = "auto"`, the fit is
#' compressed when at most half the rows are distinct patterns, and the
#' search for patterns gives up as soon as that is exceeded; the search is
#' skipped when the option `design.matrix.dedup` is `FALSE`. `compress =
#' TRUE` always compresses (unless every row is distinct) and `FALSE` never
#' does.
#'
#' **What is kept.** The fit keeps the coefficients, the factor `R`, the
#' residual and explained sums of squares, and the design specification,
#' which [predict()] applies to `newdata` with `design_matrix(spec,
#' newdata)`. Per-row residuals and fitted values are not stored: with
#' `model = TRUE` (the default, as in [lm()]) the fit keeps the model frame
#' and its design matrix, which share memory with the frame, and
#' [residuals()], [fitted()] and `predict()` without `newdata` recompute
#' them (one product with the design). With `model = FALSE` these need
#' `newdata`.
#'
#' **Methods.** [coef()], [vcov()] (with `NA` rows and columns for aliased
#' coefficients, as `vcov(lm_fit)`), [sigma()], [df.residual()], [nobs()],
#' [deviance()] (the weighted residual sum of squares), [logLik()] (hence
#' [AIC()] and [BIC()]), [formula()], [terms()], [model.frame()],
#' [residuals()], [fitted()], [predict()], `print()` and [summary()] (a
#' coefficient table, R-squared and F statistic, as `summary.lm()`).
#'
#' **Differences with `lm()`.** Multivariate responses and `method =
#' "model.frame"` are not supported. The result is not an `lm` object: it
#' has no `qr`, `effects`, `assign` or `residuals` components, so functions
#' that expect those (such as `plot.lm()` or `anova.lm()`) do not apply.
#' In `predict()`, aliased coefficients count as zero and non-estimable
#' predictions on new data are flagged, as in `predict.lm()`. Prediction
#' intervals use the variance `sigma^2 / w` at the training rows of a
#' weighted fit and `sigma^2` for new observations, with the warnings of
#' `predict.lm()`; its `weights`, `pred.var` and `scale` arguments are not
#' supported. With `na.action = na.omit`, offsets of `newdata` are taken
#' for the rows kept (`predict.lm()` would recycle them).
#'
#' **Remote data.** `design_lm()` also fits remote data (a DuckDB table or
#' query, or Parquet files) and remote designs, with the same arguments and
#' methods; see the *Linear models* section of [design_remote] for how, and
#' for the differences from in-memory fits.
#'
#' @param formula A model formula, as in [stats::lm()], or a
#'   [design_matrix()] (then give the response as `y`).
#' @param data,subset,weights,na.action,offset,contrasts As in [stats::lm()].
#'   With a design matrix, `weights` and `offset` are numeric vectors with one
#'   value per row, and `data`, `subset`, `na.action` and `contrasts` do not
#'   apply. `data` can also be remote (a DuckDB connection with `table` or
#'   `sql`, or Parquet paths; see [design_remote]); `weights`, `offset` and
#'   `subset` are then columns, `DBI::SQL()` or R expressions with a SQL
#'   translation, as `y`, `weights` and `offset` with a remote design.
#' @param y With a design matrix as first argument, the numeric response,
#'   one value per row. With a remote design (see [design_remote]), a
#'   column, `DBI::SQL()` or R expression over its data (default: the
#'   response of its formula).
#' @param compress `"auto"`, `TRUE` or `FALSE`: fit on the sufficient
#'   statistics of rows that share a design pattern (see Details).
#' @param tol In `design_lm()`, the tolerance for
#'   detecting aliased columns, as in [stats::lm.fit()]. In `predict()`, the
#'   tolerance of the check for non-estimable predictions from a
#'   rank-deficient fit, as in [stats::predict.lm()].
#' @param model Keep the model frame and its design matrix, for
#'   [residuals()], [fitted()] and [model.frame()].
#' @param ... For remote data (a DuckDB connection or Parquet paths as
#'   `data`, or as `newdata` in `predict()`): `table`, `sql`, `engine` and
#'   `chunk_size`; see [design_remote]. Otherwise not used.
#' @return An object of class `"design_lm"`: a list with, among others,
#'   `coefficients`, `sigma`, `df.residual`, `nobs`, `rank`, `pivot`, `R`
#'   (the `p` by `p` factor of [design_qr], columns in pivoted order),
#'   `rss`, `compressed` (whether the fit was compressed) and `patterns`
#'   (the number of rows or patterns the fit was computed on), `spec` (the
#'   design specification), `call` and `terms`.
#' @examples
#' set.seed(1)
#' d <- data.frame(y = rnorm(1000), x = rnorm(1000),
#'                 f = factor(sample(letters[1:5], 1000, TRUE)),
#'                 g = factor(sample(1:3, 1000, TRUE)))
#' fit <- design_lm(y ~ x * f, d)
#' all.equal(coef(fit), coef(lm(y ~ x * f, d)))
#' all.equal(vcov(fit), vcov(lm(y ~ x * f, d)))
#' summary(fit)
#'
#' ## factors only: fitted on 15 patterns instead of 1000 rows
#' fit <- design_lm(y ~ f * g, d)
#' fit$patterns
#' all.equal(sigma(fit), sigma(lm(y ~ f * g, d)))
#'
#' predict(fit, newdata = d[1:3, ])
#' @export
design_lm <- function(formula, data, subset, weights, na.action, offset,
                      contrasts = NULL, compress = "auto", tol = 1e-07,
                      model = TRUE, y = NULL, vcov = "iid", fixef.rm = "singleton", ...) {
    if (inherits(formula, c("design_matrix", "dm_remote"))) {
        if (!missing(data) || !missing(subset) || !missing(na.action) || !is.null(contrasts))
            stop("With a design matrix, give the response as `y`; `data`, `subset`, ",
                 "`na.action` and `contrasts` do not apply.", call. = FALSE)
        if (...length()) stop("Unknown arguments: ", paste(names(list(...)), collapse = ", "), call. = FALSE)
        if (!identical(vcov, "iid"))
            stop("`vcov` other than \"iid\" requires fixed effects (y ~ x | fe).", call. = FALSE)
        if (inherits(formula, "dm_remote")) {
            ## a remote design (see remote_fit.R): `y`, `weights` and `offset`
            ## are columns, DBI::SQL() or R expressions over the columns
            if (!is.null(formula$scale))
                stop("Fit an unscaled remote design, with `weights` naming a column.", call. = FALSE)
            dmfit_check_compress(compress)
            fit <- dm_remote_lm_design(formula, substitute(y), if (missing(weights)) NULL else substitute(weights),
                                       if (missing(offset)) NULL else substitute(offset), compress, tol,
                                       parent.frame())
        } else {
            fit <- dmfit_design(formula, y, if (missing(weights)) NULL else weights,
                                if (missing(offset)) NULL else offset, compress, tol)
        }
        fit$call <- match.call()
        return(fit)
    }
    if (!is.null(y)) stop("`y` applies only when the first argument is a design matrix.", call. = FALSE)
    ## fixed effects (fe.R)
    if (dm_fe_has_bar(formula))
        return(dm_fe_lm(match.call(), parent.frame(), formula, contrasts, compress, tol, missing(tol),
                        model, vcov, fixef.rm, list(...)))
    if (!identical(vcov, "iid"))
        stop("`vcov` other than \"iid\" requires fixed effects (y ~ x | fe).", call. = FALSE)
    dmfit_check_compress(compress)
    checkmate::assert_number(tol, lower = 0, finite = TRUE)
    checkmate::assert_flag(model)
    cl <- match.call()
    ## remote data (see remote_fit.R)
    if (!missing(data) && dm_is_remote_source(data)) {
        return(dm_remote_lm(formula, data,
                            weights = if (missing(weights)) NULL else substitute(weights),
                            offset = if (missing(offset)) NULL else substitute(offset),
                            subset = if (missing(subset)) NULL else substitute(subset),
                            na.action = if (missing(na.action)) NULL else na.action,
                            contrasts, compress, tol, cl, list(...), parent.frame()))
    }
    if (...length()) stop("Unknown arguments: ", paste(names(list(...)), collapse = ", "), call. = FALSE)
    ## the model frame, as lm() builds it
    mf <- dmfit_model_frame(cl, c("subset", "weights", "na.action", "offset"), parent.frame())
    fr <- dmfit_frame(mf, contrasts)
    fit <- dmfit_lm(fr$design, fr$y, fr$weights, fr$offset, compress, tol,
                    intercept = attr(fr$terms, "intercept") > 0L)
    fit$call <- cl
    fit$terms <- fr$terms
    fit$contrasts <- fr$design$contrasts
    fit$xlevels <- stats::.getXlevels(fr$terms, mf)
    fit$na.action <- attr(mf, "na.action")
    fit$offset_call <- cl$offset
    fit$env <- environment(fr$terms)
    if (model) {
        fit$model <- mf
        fit$data <- list(design = fr$design, y = fr$y, weights = fr$weights,
                         offset = fr$offset)
    }
    fit
}

## design_lm() on a design matrix x and a response y.
dmfit_design <- function(x, y, weights = NULL, offset = NULL,
                         compress = "auto", tol = 1e-07) {
    checkmate::assert_class(x, "design_matrix")
    dmfit_check_compress(compress)
    checkmate::assert_number(tol, lower = 0, finite = TRUE)
    if (is.matrix(y) && ncol(y) == 1L) y <- drop(y)
    checkmate::assert_numeric(y, len = nrow(x))
    y <- as.double(y)
    if (!is.null(weights)) checkmate::assert_numeric(weights, len = nrow(x))
    if (!is.null(offset)) checkmate::assert_numeric(offset, len = nrow(x))
    tt <- x$training_terms
    intercept <- if (!is.null(tt)) attr(tt, "intercept") > 0L else
        "(Intercept)" %in% colnames(x)
    fit <- dmfit_lm(x, y, weights, offset, compress, tol, intercept = intercept)
    fit$data <- list(design = x, y = y, weights = weights, offset = offset)
    fit
}

## The model frame of a fit, as lm() and glm() build it: design_frame() on
## the formula (or `formula`), `data` and the arguments `args` of the
## matched call, with unused levels dropped.
dmfit_model_frame <- function(cl, args, env, formula = NULL) {
    mf <- cl[c(1L, match(c("formula", "data", args), names(cl), 0L))]
    if (!is.null(formula)) mf$formula <- formula
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(design.matrix::design_frame)
    eval(mf, env)
}

dmfit_check_compress <- function(compress) {
    if (!(isTRUE(compress) || isFALSE(compress) || identical(compress, "auto")))
        stop("`compress` must be \"auto\", TRUE or FALSE.", call. = FALSE)
}

## The parts of a model frame that lm() and glm() use: terms, response,
## weights, offset (offset() terms and the `offset` argument, summed by
## model.offset()) and the design matrix.
dmfit_frame <- function(mf, contrasts = NULL) {
    mt <- attr(mf, "terms")
    y <- dmfit_response(mf)
    if (is.matrix(y))
        stop("Multivariate responses are not supported.", call. = FALSE)
    w <- as.vector(stats::model.weights(mf))
    if (!is.null(w) && !is.numeric(w))
        stop("'weights' must be a numeric vector", call. = FALSE)
    offset <- as.vector(stats::model.offset(mf))
    if (!is.null(offset) && length(offset) != NROW(y))
        stop(gettextf("number of offsets is %d, should equal %d (number of observations)",
                      length(offset), NROW(y)), call. = FALSE)
    design <- design_matrix(mt, mf, contrasts.arg = contrasts)
    list(terms = mt, y = y, weights = w, offset = offset, design = design)
}

## The response of a model frame, as model.response(mf, "numeric") but
## without names: on a large frame, turning compact row names into a
## character vector costs more than the fit.
dmfit_response <- function(mf) {
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

## Response minus offset, and checks of the response, weights and offset,
## as lm.wfit() makes them, in one pass (one query remotely). Columns:
## numeric vectors in memory, remote columns (primitives.R).
dmfit_check <- function(n, y, weights, offset) {
    if (n == 0L) stop("0 (non-NA) cases", call. = FALSE)
    z <- if (is.null(offset)) y else y - offset
    w <- weights
    s <- dmc_sums(dmc_count(dmc_is_finite(z), negate = TRUE),
                  if (!is.null(w)) dmc_count(dmc_notnan(w) & w >= 0, negate = TRUE) else 0,
                  if (!is.null(w)) dmc_count(abs(w) == Inf) else 0,
                  if (!is.null(w)) dmc_count(w != 0) else 1)
    if (s[1L] > 0) stop("NA/NaN/Inf in 'y'", call. = FALSE)
    if (s[2L] > 0) stop("missing or negative weights not allowed", call. = FALSE)
    if (s[3L] > 0) stop("NA/NaN/Inf in 'weights'", call. = FALSE)
    if (s[4L] == 0) stop("0 (non-NA) cases", call. = FALSE)
    z
}

## Number of observations (non-zero weights) and sum of the log weights.
dmfit_wstats <- function(n, weights) {
    if (is.null(weights)) return(list(nobs = n, sumlogw = 0))
    s <- dmc_sums(dmc_count(weights != 0), dmc_ifelse(weights != 0, log(weights), 0))
    list(nobs = s[1L], sumlogw = s[2L])
}

## Groups of rows that share a design pattern, or NULL when the fit should
## not be compressed. `compress = "auto"` follows the option
## `design.matrix.dedup` (its "auto" search gives up beyond n / 2
## patterns); TRUE searches in any case. Either way, a fit on more than
## half as many patterns as rows is not worth compressing.
dmfit_groups <- function(x, compress) {
    if (isFALSE(compress) || nrow(x) < 2L) return(NULL)
    J <- as_jac(x)
    if (!jac_kernel_ok(J)) return(NULL)
    g <- if (isTRUE(compress)) dm_groups(J, force = TRUE) else dm_groups(J, cost = Inf)
    if (is.null(g)) return(NULL)
    if (length(g$first) == nrow(x)) return(NULL)
    if (!isTRUE(compress) && length(g$first) > 0.5 * nrow(x)) return(NULL)
    g
}

## Per-pattern sufficient statistics of z with weights w (NULL for 1):
## the sum of the weights, the weighted mean, and the weighted sum of
## squares within patterns (computed around the means, which avoids the
## cancellation of sum(w z^2) - sum(w) mean^2).
dmfit_suffstats <- function(z, w, map, G) {
    if (is.null(w)) {
        sw <- tabulate(map, nbins = G)
        swz <- drop(rowsum(z, map, reorder = TRUE))
    } else {
        ## doubles: integer weights could overflow in the sums
        w <- as.double(w)
        sw <- drop(rowsum(w, map, reorder = TRUE))
        swz <- drop(rowsum(w * z, map, reorder = TRUE))
    }
    m <- ifelse(sw > 0, swz / sw, 0)
    within <- if (is.null(w)) {
        sum((z - m[map])^2)
    } else {
        ## rows with zero weight do not enter the fit (their squared
        ## deviation can overflow, and 0 * Inf is NaN)
        pos <- w > 0
        sum(w[pos] * (z[pos] - m[map[pos]])^2)
    }
    list(w = as.double(sw), mean = m, within = within)
}

## A least squares fit of z (response minus offset) on x with weights w,
## compressed to patterns when dmfit_groups() finds them.
dmfit_lm <- function(x, y, weights = NULL, offset = NULL, compress = "auto",
                     tol = 1e-07, intercept = TRUE) {
    n <- nrow(x)
    z <- dmfit_check(n, y, weights, offset)
    st <- dmfit_wstats(n, weights)
    nobs <- st$nobs
    ## integers, as in lm() (remote counts are doubles)
    if (is.double(nobs) && nobs < .Machine$integer.max) nobs <- as.integer(nobs)
    sumlogw <- st$sumlogw
    g <- dmp_patterns(x, z, weights, compress)
    if (!is.null(g)) {
        xs <- g$x
        zs <- g$z
        ws <- g$w
        within <- g$within
        n <- g$rows
    } else {
        xs <- x
        zs <- z
        ws <- weights
        ## (rows with zero weight count for nothing in the sums)
        within <- 0
    }
    out <- dmfit_ls(xs, zs, ws, tol, intercept)
    out$rss <- out$rss + within
    out$df.residual <- nobs - out$rank
    ## as many observations as coefficients: the fit is exact, and what is
    ## left of the residuals is rounding (lm() gives exact zeros)
    if (out$df.residual == 0L) out$rss <- 0
    out$sigma <- sqrt(out$rss / out$df.residual)
    out$nobs <- nobs
    out$sumlogw <- sumlogw
    out$weighted <- !is.null(weights)
    out$compressed <- !is.null(g)
    out$patterns <- nrow(xs)
    out$rows <- n
    out$spec <- dmfit_spec(x)
    class(out) <- "design_lm"
    out
}

## Weighted least squares of z on the rows of x (design rows or patterns),
## with weights w (NULL for 1): coefficients (NA for aliased columns),
## the pivoted factor, the weighted residual sum of squares and the
## explained sum of squares (around the weighted mean with an intercept,
## around zero without), as summary.lm() computes it.
dmfit_ls <- function(x, z, w, tol, intercept, sums = TRUE) {
    p <- ncol(x)
    cn <- colnames(x)
    if (p == 0L) {
        rss <- dm_colsum(if (is.null(w)) z^2 else dmc_ifelse(w > 0, w * z^2, 0))
        return(list(coefficients = numeric(0),
                    rank = 0L, pivot = integer(0), R = matrix(0, 0, 0),
                    tol = tol, rss = rss, mss = 0, intercept = intercept))
    }
    ## (the primitives of primitives.R: columns are numeric vectors in
    ## memory, remote columns for a remote design)
    sw <- if (is.null(w)) NULL else sqrt(w)
    Y <- if (is.null(sw)) z else sw * z
    q <- dmqr_for(dmp_scale(x, sw), tol, Y)
    b <- dmqr_coef(q, Y)
    beta <- dmqr_full(q, b)[, 1L]
    coef <- beta
    coef[setdiff(seq_len(p), q$pivot[seq_len(q$rank)])] <- NA_real_
    names(coef) <- cn
    ## (the weighted least squares steps of a GLM need only the solution)
    if (!sums) return(list(coefficients = coef, rank = q$rank, pivot = q$pivot, R = q$R, tol = tol))
    f <- dmp_mult(x, beta, along = list(z, w))
    ww <- if (is.null(w)) 1 else w
    ## w v^2, 0 for rows with zero weight (whose squares can overflow)
    wsq <- function(v) if (is.null(w)) v^2 else dmc_ifelse(w > 0, w * v^2, 0)
    mo <- attr(b, "moments")
    rss <- NULL
    if (!is.null(mo)) {
        ## from the moments of the solution, when the factorization has them
        ## (a remote design: no pass over its products): rss = z'Wz - b.h,
        ## mss = b'Gb - (sum w f)^2 / sum w, with the intercept column of G
        ## giving sum w f; dmc_ls_sums() below when cancellation would
        ## cost accuracy
        zz <- dmc_sums(wsq(z), ww + 0 * z)
        bh <- sum(mo$b * mo$h)
        rss_m <- zz[1L] - bh
        ic <- match(which(colnames(x) == "(Intercept)"), q$pivot[seq_len(q$rank)])
        if (rss_m > 1e-6 * zz[1L] && (!intercept || !is.na(ic))) {
            rss <- rss_m
            mss <- if (intercept) bh - sum(mo$b * mo$G[, ic])^2 / zz[2L] else bh
            mss <- max(mss, 0)
        }
    }
    if (is.null(rss)) {
        ## (primitives_fused.R: in memory one pass in C)
        s <- dmc_ls_sums(z, f, w, intercept)
        rss <- s[1L]
        mss <- s[2L]
    }
    R <- q$R
    list(coefficients = coef, rank = q$rank, pivot = q$pivot, R = R, tol = tol,
         rss = rss, mss = mss, intercept = intercept)
}

## sigma^2 (R11' R11)^-1 in the original column order, NA for aliased
## columns.
dmfit_vcov <- function(object, scale) {
    p <- length(object$coefficients)
    if (p == 0L) return(matrix(numeric(0), 0L, 0L))
    nm <- names(object$coefficients)
    V <- matrix(NA_real_, p, p, dimnames = list(nm, nm))
    k <- object$rank
    if (k > 0L) {
        keep <- object$pivot[seq_len(k)]
        V[keep, keep] <- scale * chol2inv(unname(object$R[seq_len(k), seq_len(k), drop = FALSE]))
    }
    V
}

## The design specification for predictions on new data (NULL when the
## design cannot be reused, such as a dense design from a fitted model).
dmfit_spec <- function(x) {
    tryCatch(design_matrix(x, newdata = NULL), error = function(e) NULL)
}

## Coefficients with aliased columns set to zero.
dmfit_beta <- function(object) {
    b <- object$coefficients
    b[is.na(b)] <- 0
    unname(b)
}

## The kept data, or an error.
dmfit_data <- function(object, what) {
    if (is.null(object$data))
        stop(sprintf("%s needs the fit's data: refit with `model = TRUE`, or use `predict(fit, newdata)`.",
                     what), call. = FALSE)
    object$data
}

## X b for a design, without the offset.
dmfit_linpred <- function(object, D, along = list()) {
    if (ncol(D) == 0L) return(dmc_const(D, 0))
    dmp_mult(D, dmfit_beta(object), along)
}

#' @export
print.design_lm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    if (!is.null(x$call)) cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
    dmfit_print_rows(x)
    dmfit_print_coef(x, digits)
    cat("\n")
    invisible(x)
}

## The coefficient table of summary(), with rows of NA for aliased
## coefficients, as print.summary.lm() prints it.
dmfit_print_coefmat <- function(x, digits, signif.stars, ...) {
    if (length(x$aliased) == 0L) {
        cat("\nNo Coefficients\n")
        return(invisible())
    }
    nsing <- sum(x$aliased)
    coefs <- x$coefficients
    if (nsing > 0L) {
        cat("\nCoefficients: (", nsing, " not defined because of singularities)\n", sep = "")
        cn <- names(x$aliased)
        coefs <- matrix(NA_real_, length(x$aliased), 4L, dimnames = list(cn, colnames(x$coefficients)))
        coefs[!x$aliased, ] <- x$coefficients[cn[!x$aliased], ]
    } else {
        cat("\nCoefficients:\n")
    }
    stats::printCoefmat(coefs, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
}

## How the fit was computed: on patterns, or (remotely) from sums over the
## rows (`method`).
dmfit_print_rows <- function(x, prefix = "\n") {
    if (isTRUE(x$compressed)) {
        cat(sprintf("%sFitted on %s design patterns of %s rows.\n", prefix, dm_fmt_count(x$patterns), dm_fmt_count(x$rows)))
    } else if (!is.null(x$method)) {
        cat(sprintf("%sFitted from %s of %s rows.\n", prefix, x$method, dm_fmt_count(x$rows)))
    }
}

dmfit_print_coef <- function(x, digits, none = "\nNo coefficients\n") {
    if (length(x$coefficients)) {
        cat("\nCoefficients:\n")
        print.default(format(x$coefficients, digits = digits), print.gap = 2L, quote = FALSE)
    } else {
        cat(none)
    }
}

#' @export
coef.design_lm <- function(object, ...) object$coefficients

#' @export
vcov.design_lm <- function(object, ...) dmfit_vcov(object, object$sigma^2)

#' @export
sigma.design_lm <- function(object, ...) object$sigma

#' @export
nobs.design_lm <- function(object, ...) object$nobs

#' @export
deviance.design_lm <- function(object, ...) object$rss

#' @export
df.residual.design_lm <- function(object, ...) object$df.residual

#' @export
logLik.design_lm <- function(object, ...) {
    N <- object$nobs
    val <- 0.5 * (object$sumlogw - N * (log(2 * pi) + 1 - log(N) + log(object$rss)))
    structure(val, nall = N, nobs = N, df = object$rank + 1L, class = "logLik")
}

#' @export
formula.design_lm <- function(x, ...) {
    if (is.null(x$terms)) return(NULL)
    stats::formula(x$terms)
}

#' @export
terms.design_lm <- function(x, ...) {
    if (is.null(x$terms)) stop("The fit has no terms: it was not fitted from a formula.", call. = FALSE)
    x$terms
}

#' @export
model.frame.design_lm <- function(formula, ...) {
    if (is.null(formula$model))
        stop("The fit has no model frame: refit `design_lm()` with `model = TRUE`.", call. = FALSE)
    formula$model
}

#' @export
residuals.design_lm <- function(object, ...) {
    ## (y - offset) - X b on the rows of the fit: a vector in memory, the
    ## remote design with .fitted and .resid columns remotely
    rw <- dmfit_predict_rows(object, NULL, NULL, list(), "residuals()")
    xb <- dmfit_linpred(object, rw$D, along = list(rw$z))
    fit <- if (is.null(rw$off)) xb else xb + rw$off
    dmp_output(rw$D, list(.fitted = fit, .resid = rw$z - xb), design = TRUE, memory = function(d) {
        r <- d$.resid
        names(r) <- dimnames(rw$D)[[1L]]
        stats::naresid(rw$na_act, r)
    })
}

#' @export
fitted.design_lm <- function(object, ...) {
    ## (one implementation with predict(): a vector in memory; a table with
    ## .fitted and .resid remotely)
    dmfit_predict(object, resid = TRUE, what = "fitted()")
}

#' @param object A `"design_lm"` object.
#' @param newdata A data frame of new observations. Omitted, predictions
#'   are the fitted values (which need `model = TRUE`).
#' @param se.fit Return standard errors of the predictions?
#' @param interval Type of interval: `"none"`, `"confidence"` or
#'   `"prediction"`.
#' @param level Confidence level of the intervals.
#' @param rankdeficient For a rank-deficient fit, how predictions on
#'   `newdata` that are not estimable are flagged, as in
#'   [stats::predict.lm()].
#' @rdname design_lm
#' @export
predict.design_lm <- function(object, newdata, se.fit = FALSE,
                              interval = c("none", "confidence", "prediction"),
                              level = 0.95, na.action = stats::na.pass,
                              rankdeficient = c("warnif", "simple", "non-estim", "NA", "NAwarn"),
                              tol = 1e-06, ...) {
    interval <- match.arg(interval)
    rankdeficient <- match.arg(rankdeficient)
    checkmate::assert_flag(se.fit)
    checkmate::assert_number(level, lower = 0, upper = 1)
    no_data <- missing(newdata) || is.null(newdata)
    dmfit_predict(object, if (no_data) NULL else newdata, na.action, list(...), se.fit = se.fit,
                  interval = interval, level = level, rankdeficient = rankdeficient, tol = tol)
}

## The rows to predict on: the design D (in memory or remote), its offset
## (a column, or NULL), the weights of the training rows (for prediction
## intervals), the response minus offsets of the training rows (for
## residuals) and their na.action. One description for both backends.
dmfit_predict_rows <- function(object, newdata, na.action, dots, what) {
    if (is.null(newdata)) {
        dat <- dmfit_data(object, what)
        z <- if (is.null(dat$offset)) dat$y else dat$y - dat$offset
        return(list(D = dat$design, off = dat$offset, w = dat$weights, z = z, training = TRUE,
                    na_act = object$na.action))
    }
    nd <- dmfit_newdata(object, newdata, na.action, dots)
    list(D = nd$D, off = nd$off, w = NULL, z = NULL, training = FALSE, na_act = NULL)
}

## The design of new data for the predictions of a fit (design_lm(), with
## or without fixed effects, and design_glm()), with its offset (offset()
## terms and the `offset` argument of the fit, evaluated on the new rows,
## as predict.lm() evaluates them): list(D, off). The data are a data frame
## (in memory) or remote data (a DuckDB connection or Parquet paths, with
## the `table`, `sql`, `engine` and `chunk_size` of `dots`).
dmfit_newdata <- function(object, newdata, na.action = stats::na.pass, dots = list()) {
    if (is.null(object$spec)) stop("This design cannot be applied to new data.", call. = FALSE)
    if (dm_is_remote_source(newdata))
        return(dm_remote_new_offset(object, do.call(design_matrix, c(list(object$spec, newdata, na.action = na.action), dots))))
    if (length(dots)) stop("Unknown arguments: ", paste(names(dots), collapse = ", "), call. = FALSE)
    D <- design_matrix(object$spec, newdata, na.action = na.action)
    list(D = D, off = dmfit_new_offset(object, newdata, D))
}

## Predictions, standard errors, intervals and non-estimable cases of a
## linear fit, computed once on the primitives (primitives.R); the result is
## a vector, matrix or list in memory (as predict.lm()), a table remotely.
dmfit_predict <- function(object, newdata = NULL, na.action = stats::na.pass, dots = list(),
                          se.fit = FALSE, interval = "none", level = 0.95,
                          rankdeficient = "warnif", tol = 1e-06, resid = FALSE, what = "predict() without `newdata`") {
    rw <- dmfit_predict_rows(object, newdata, na.action, dots, what)
    D <- rw$D
    no_data <- is.null(newdata)
    fit <- dmfit_linpred(object, D, along = list(rw$off, rw$w, rw$z))
    xb <- fit
    if (!is.null(rw$off)) fit <- fit + rw$off
    nonest <- NULL
    ## as predict.lm(), only predictions on new data are checked
    if (!no_data && object$rank < length(object$coefficients)) {
        msg <- gettext("prediction from rank-deficient fit")
        if (rankdeficient == "simple") {
            warning(gettextf("%s; consider predict(., rankdeficient=\"NA\")", msg), call. = FALSE)
        } else {
            checkmate::assert_number(tol, lower = 0)
            nonest <- dmfit_nonestimable(object, D, tol, along = list(fit))
            k <- dm_colsum(dmc_count(nonest))
            if (startsWith(rankdeficient, "NA")) {
                fit <- dmc_ifelse(nonest, NA_real_, fit)
                if (rankdeficient == "NAwarn" && k > 0)
                    warning(gettextf("%s: NAs produced for non-estimable cases", msg), call. = FALSE)
                nonest <- NULL
            } else {
                if (rankdeficient == "warnif" && k > 0)
                    warning(gettextf("%s; attr(*, \"non-estim\") has doubtful cases", msg), call. = FALSE)
                if (rankdeficient != "non-estim" && k == 0) nonest <- NULL
            }
        }
    }
    se <- lwr <- upr <- NULL
    if (se.fit || interval != "none") {
        ## the unscaled covariance with aliased rows and columns zeroed, then
        ## the residual variance (NaN for a saturated fit, as in predict.lm())
        res_var <- object$sigma^2
        V <- dmfit_vcov(object, 1)
        V[is.na(V)] <- 0
        ip <- if (ncol(D) == 0L) dmc_const(fit, 0) else dmp_quad(D, V, along = list(fit))
        se <- sqrt(ip * res_var)
        if (interval != "none") {
            pred_var <- res_var
            if (interval == "prediction") {
                ## the warnings and weights of predict.lm()
                if (no_data) {
                    warning("predictions on current data refer to _future_ responses\n", call. = FALSE)
                    if (!is.null(rw$w)) {
                        pred_var <- res_var / rw$w
                        warning("assuming prediction variance inversely proportional to weights used for fitting\n",
                                call. = FALSE)
                    }
                } else if (object$weighted) {
                    warning("Assuming constant prediction variance even though model fit is weighted\n",
                            call. = FALSE)
                }
            }
            tq <- stats::qt((1 - level) / 2, object$df.residual)
            hw <- tq * if (interval == "confidence") se else sqrt(se^2 + pred_var)
            lwr <- fit + hw
            upr <- fit - hw
        }
    }
    cols <- list(.fitted = fit)
    if (resid) cols$.resid <- rw$z - xb
    if (se.fit) cols$.se.fit <- se
    if (!is.null(lwr)) cols$.lwr <- lwr
    if (!is.null(upr)) cols$.upr <- upr
    if (!is.null(nonest)) cols$.nonestim <- nonest
    na_act <- rw$na_act
    ## remotely a table; in memory, the values of predict.lm()
    dmp_output(D, cols, na = na_act,
               attrs = if (se.fit) list(df = object$df.residual, residual.scale = object$sigma),
               memory = function(d) {
        rn <- dimnames(D)[[1L]]
        fit <- d$.fitted
        names(fit) <- rn
        if (!is.null(nonest)) {
            ne <- which(d$.nonestim)
            names(ne) <- rn[ne]
            attr(fit, "non-estim") <- ne
        }
        if (!se.fit && interval == "none") return(stats::napredict(na_act, fit))
        res <- if (interval != "none") cbind(fit = fit, lwr = d$.lwr, upr = d$.upr) else fit
        res <- stats::napredict(na_act, res)
        if (!se.fit) return(res)
        se <- d$.se.fit
        names(se) <- rn
        list(fit = res, se.fit = stats::napredict(na_act, se), df = object$df.residual,
             residual.scale = object$sigma)
    })
}

## Rows of D whose predictions are not estimable from a rank-deficient fit:
## rows with a component, relative to their norm, in the null space of the
## design (as predict.lm() determines them).
dmfit_nonestimable <- function(object, D, tol, along = list()) {
    p <- length(object$coefficients)
    k <- object$rank
    tR <- t(object$R)
    tR[cbind(seq_len(p), seq_len(p))[seq_len(p) > k, , drop = FALSE]] <- 1
    nb <- qr.Q(base::qr(tR))[, (k + 1L):p, drop = FALSE]
    B <- matrix(0, p, p - k)
    B[object$pivot, ] <- nb
    ## (a logical vector in memory, a logical column remotely)
    Xb <- dmp_mult(D, B, along)
    Xb_norm <- sqrt(dmc_rowsums(Xb * Xb))
    X_norm <- sqrt(dmp_quad(D, diag(1, p), along = list(Xb_norm)))
    tol * X_norm <= Xb_norm
}

## The offset of new data: offset() terms of the formula and the `offset`
## argument of the call, evaluated on `newdata` as predict.lm() evaluates
## them, for the rows of `D`, the design of `newdata` (which can omit rows
## with missing values).
dmfit_new_offset <- function(object, newdata, D) {
    tt <- object$terms
    if (is.null(tt)) return(NULL)
    off <- NULL
    env <- environment(tt)
    vars <- attr(tt, "variables")
    for (i in attr(tt, "offset")) {
        o <- eval(vars[[i + 1L]], newdata, env)
        off <- if (is.null(off)) o else off + o
    }
    if (!is.null(object$offset_call)) {
        o <- eval(object$offset_call, newdata, object$env %||% env)
        off <- if (is.null(off)) o else off + o
    }
    if (is.null(off)) return(NULL)
    off <- as.vector(off)
    n <- nrow(D)
    if (length(off) != n) {
        ## rows omitted by na.action: keep the offsets of the rows of D,
        ## which carry the row names of newdata
        rn <- if (is.data.frame(newdata)) row.names(newdata) else as.character(seq_along(off))
        if (length(rn) != length(off))
            stop(gettextf("number of offsets is %d, should equal %d (number of observations)",
                          length(off), n), call. = FALSE)
        off <- off[match(dimnames(D)[[1L]], rn)]
    }
    off
}

#' @export
summary.design_lm <- function(object, ...) {
    p <- object$rank
    rdf <- object$df.residual
    rss <- object$rss
    resvar <- rss / rdf
    est <- object$coefficients
    aliased <- is.na(est)
    V <- dmfit_vcov(object, resvar)
    keep <- !aliased
    se <- sqrt(diag(V)[keep])
    tval <- est[keep] / se
    coefs <- cbind(Estimate = est[keep], "Std. Error" = se, "t value" = tval,
                   "Pr(>|t|)" = 2 * stats::pt(abs(tval), rdf, lower.tail = FALSE))
    ans <- list(call = object$call, terms = object$terms, coefficients = coefs,
                aliased = aliased, sigma = sqrt(resvar), df = c(p, rdf, length(est)),
                compressed = object$compressed, patterns = object$patterns,
                rows = object$rows, method = object$method)
    df.int <- if (object$intercept) 1L else 0L
    if (p != df.int) {
        mss <- object$mss
        n <- object$nobs
        ans$r.squared <- mss / (mss + rss)
        ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - df.int) / rdf)
        ans$fstatistic <- c(value = (mss / (p - df.int)) / resvar,
                            numdf = p - df.int, dendf = rdf)
    } else {
        ans$r.squared <- ans$adj.r.squared <- 0
    }
    ans$cov.unscaled <- dmfit_vcov(object, 1)
    class(ans) <- "summary.design_lm"
    ans
}

#' @export
print.summary.design_lm <- function(x, digits = max(3L, getOption("digits") - 3L),
                                    signif.stars = getOption("show.signif.stars"), ...) {
    if (!is.null(x$call)) cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
    dmfit_print_rows(x)
    dmfit_print_coefmat(x, digits, signif.stars, ...)
    cat("\nResidual standard error:", format(signif(x$sigma, digits)),
        "on", x$df[2L], "degrees of freedom\n")
    if (!is.null(x$fstatistic)) {
        cat("Multiple R-squared: ", formatC(x$r.squared, digits = digits))
        cat(",\tAdjusted R-squared: ", formatC(x$adj.r.squared, digits = digits),
            "\nF-statistic:", formatC(x$fstatistic[1L], digits = digits), "on",
            x$fstatistic[2L], "and", x$fstatistic[3L], "DF,  p-value:",
            format.pval(stats::pf(x$fstatistic[1L], x$fstatistic[2L], x$fstatistic[3L],
                                  lower.tail = FALSE), digits = digits))
        cat("\n")
    }
    cat("\n")
    invisible(x)
}

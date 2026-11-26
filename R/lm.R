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
#' without forming the dense model matrix. `design_lm_fit()` does the same
#' from a [design_matrix()] and a response, like [stats::lm.fit()].
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
#' predictions on new data are flagged, as in `predict.lm()`; prediction
#' intervals assume a constant variance
#' `sigma^2` for new observations, also for a weighted fit.
#'
#' @param formula,data,subset,weights,na.action,offset,contrasts As in
#'   [stats::lm()].
#' @param compress `"auto"`, `TRUE` or `FALSE`: fit on the sufficient
#'   statistics of rows that share a design pattern (see Details).
#' @param tol The tolerance for detecting aliased columns, as in
#'   [stats::lm.fit()].
#' @param model Keep the model frame and its design matrix, for
#'   [residuals()], [fitted()] and [model.frame()].
#' @param x A `design_matrix`.
#' @param y A numeric response vector with one value per row of `x`.
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
                      model = TRUE) {
    dmfit_check_compress(compress)
    checkmate::assert_number(tol, lower = 0, finite = TRUE)
    checkmate::assert_flag(model)
    cl <- match.call()
    ## the model frame, as lm() builds it
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(design.matrix::design_frame)
    mf <- eval(mf, parent.frame())
    fr <- dmfit_frame(mf, contrasts)
    fit <- dmfit_lm(fr$design, fr$y, fr$weights, fr$offset, compress, tol,
                    intercept = attr(fr$terms, "intercept") > 0L)
    fit$call <- cl
    fit$terms <- fr$terms
    fit$contrasts <- fr$design$contrasts
    fit$xlevels <- stats::.getXlevels(fr$terms, mf)
    fit$na.action <- attr(mf, "na.action")
    fit$offset_call <- cl$offset
    if (model) {
        fit$model <- mf
        fit$data <- list(design = fr$design, y = fr$y, weights = fr$weights,
                         offset = fr$offset)
    }
    fit
}

#' @rdname design_lm
#' @export
design_lm_fit <- function(x, y, weights = NULL, offset = NULL,
                          compress = "auto", tol = 1e-07) {
    checkmate::assert_class(x, "design_matrix")
    dmfit_check_compress(compress)
    checkmate::assert_number(tol, lower = 0, finite = TRUE)
    if (is.matrix(y) && ncol(y) == 1L) y <- drop(y)
    checkmate::assert_numeric(y, len = nrow(x))
    if (!is.null(weights)) checkmate::assert_numeric(weights, len = nrow(x))
    if (!is.null(offset)) checkmate::assert_numeric(offset, len = nrow(x))
    tt <- x$training_terms
    intercept <- if (!is.null(tt)) attr(tt, "intercept") > 0L else
        "(Intercept)" %in% colnames(x)
    fit <- dmfit_lm(x, y, weights, offset, compress, tol, intercept = intercept)
    fit$call <- match.call()
    fit$data <- list(design = x, y = y, weights = weights, offset = offset)
    fit
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
        attributes(v) <- NULL
        storage.mode(v) <- "double"
    } else if (ncol(v) == 1L) {
        v <- as.double(v)
    }
    v
}

## Response minus offset, and checks of the response, weights and offset.
dmfit_check <- function(n, y, weights, offset) {
    if (n == 0L) stop("0 (non-NA) cases", call. = FALSE)
    if (is.null(y) || length(y) != n)
        stop(sprintf("The response must have %d values.", n), call. = FALSE)
    z <- if (is.null(offset)) as.double(y) else as.double(y) - offset
    if (!all(is.finite(z))) stop("NA/NaN/Inf in 'y'", call. = FALSE)
    if (!is.null(weights)) {
        if (anyNA(weights) || any(weights < 0))
            stop("missing or negative weights not allowed", call. = FALSE)
        if (!all(is.finite(weights))) stop("NA/NaN/Inf in 'weights'", call. = FALSE)
        if (all(weights == 0)) stop("0 (non-NA) cases", call. = FALSE)
    }
    z
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
        sw <- drop(rowsum(w, map, reorder = TRUE))
        swz <- drop(rowsum(w * z, map, reorder = TRUE))
    }
    m <- ifelse(sw > 0, swz / sw, 0)
    dev2 <- (z - m[map])^2
    within <- if (is.null(w)) sum(dev2) else sum(w * dev2)
    list(w = as.double(sw), mean = m, within = within)
}

## A least squares fit of z (response minus offset) on x with weights w,
## compressed to patterns when dmfit_groups() finds them.
dmfit_lm <- function(x, y, weights = NULL, offset = NULL, compress = "auto",
                     tol = 1e-07, intercept = TRUE) {
    n <- nrow(x)
    z <- dmfit_check(n, y, weights, offset)
    nobs <- if (is.null(weights)) n else sum(weights != 0)
    sumlogw <- if (is.null(weights)) 0 else sum(log(weights[weights != 0]))
    g <- dmfit_groups(x, compress)
    if (!is.null(g)) {
        ss <- dmfit_suffstats(z, weights, g$map, length(g$first))
        keep <- ss$w > 0
        xs <- x[g$first[keep], ]
        zs <- ss$mean[keep]
        ws <- ss$w[keep]
        within <- ss$within
        if (is.null(weights) && all(ws == 1)) ws <- NULL
    } else {
        xs <- x
        zs <- z
        ws <- weights
        within <- 0
        if (!is.null(ws) && any(ws == 0)) {
            ok <- ws != 0
            xs <- x[ok, ]
            zs <- z[ok]
            ws <- ws[ok]
        }
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
dmfit_ls <- function(x, z, w, tol, intercept) {
    p <- ncol(x)
    cn <- colnames(x)
    if (p == 0L) {
        rss <- if (is.null(w)) sum(z^2) else sum(w * z^2)
        return(list(coefficients = numeric(0),
                    rank = 0L, pivot = integer(0), R = matrix(0, 0, 0),
                    tol = tol, rss = rss, mss = 0, intercept = intercept))
    }
    sw <- if (is.null(w)) NULL else sqrt(w)
    q <- qr(if (is.null(sw)) x else sw * x, tol = tol)
    b <- dmqr_coef(q, matrix(if (is.null(sw)) z else sw * z, ncol = 1L))
    beta <- dmqr_full(q, b)[, 1L]
    f <- drop(x %*% beta)
    r <- z - f
    ww <- if (is.null(w)) 1 else w
    rss <- sum(ww * r^2)
    mss <- if (intercept) {
        fbar <- if (is.null(w)) mean(f) else sum(w * f) / sum(w)
        sum(ww * (f - fbar)^2)
    } else {
        sum(ww * f^2)
    }
    coef <- beta
    coef[setdiff(seq_len(p), q$pivot[seq_len(q$rank)])] <- NA_real_
    names(coef) <- cn
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
dmfit_linpred <- function(object, D) {
    if (ncol(D) == 0L) return(rep(0, nrow(D)))
    drop(D %*% dmfit_beta(object))
}

#' @export
print.design_lm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    if (!is.null(x$call)) cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
    if (x$compressed)
        cat(sprintf("\nFitted on %d design patterns of %d rows.\n", x$patterns, x$rows))
    if (length(x$coefficients)) {
        cat("\nCoefficients:\n")
        print.default(format(x$coefficients, digits = digits), print.gap = 2L, quote = FALSE)
    } else {
        cat("\nNo coefficients\n")
    }
    cat("\n")
    invisible(x)
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
    dat <- dmfit_data(object, "residuals()")
    D <- dat$design
    y <- dat$y
    r <- y - dmfit_linpred(object, D)
    if (!is.null(dat$offset)) r <- r - dat$offset
    r <- as.vector(r)
    names(r) <- dimnames(D)[[1L]]
    stats::naresid(object$na.action, r)
}

#' @export
fitted.design_lm <- function(object, ...) {
    dat <- dmfit_data(object, "fitted()")
    D <- dat$design
    f <- dmfit_linpred(object, D)
    if (!is.null(dat$offset)) f <- f + dat$offset
    names(f) <- dimnames(D)[[1L]]
    stats::napredict(object$na.action, f)
}

#' @param object A `"design_lm"` object.
#' @param newdata A data frame of new observations. Omitted, predictions
#'   are the fitted values (which need `model = TRUE`).
#' @param se.fit Return standard errors of the predictions?
#' @param interval Type of interval: `"none"`, `"confidence"` or
#'   `"prediction"`.
#' @param level Confidence level of the intervals.
#' @param rankdeficient,tol For a rank-deficient fit, how predictions on
#'   `newdata` that are not estimable are flagged, as in
#'   [stats::predict.lm()].
#' @param ... Not used.
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
    if (missing(newdata) || is.null(newdata)) {
        dat <- dmfit_data(object, "predict() without `newdata`")
        D <- dat$design
        off <- dat$offset
        na_act <- object$na.action
    } else {
        if (is.null(object$spec))
            stop("This design cannot be applied to new data.", call. = FALSE)
        D <- design_matrix(object$spec, newdata, na.action = na.action)
        off <- dmfit_new_offset(object, newdata, na.action)
        na_act <- NULL
    }
    fit <- dmfit_linpred(object, D)
    if (!is.null(off)) fit <- fit + off
    names(fit) <- dimnames(D)[[1L]]
    if (is.null(na_act) && object$rank < length(object$coefficients)) {
        msg <- gettext("prediction from rank-deficient fit")
        if (rankdeficient == "simple") {
            warning(gettextf("%s; consider predict(., rankdeficient=\"NA\")", msg),
                    call. = FALSE)
        } else {
            checkmate::assert_number(tol, lower = 0)
            nonest <- dmfit_nonestimable(object, D, tol)
            if (startsWith(rankdeficient, "NA")) {
                fit[nonest] <- NA
                if (rankdeficient == "NAwarn" && length(nonest))
                    warning(gettextf("%s: NAs produced for non-estimable cases", msg),
                            call. = FALSE)
            } else {
                if (rankdeficient == "warnif" && length(nonest))
                    warning(gettextf("%s; attr(*, \"non-estim\") has doubtful cases", msg),
                            call. = FALSE)
                if (rankdeficient == "non-estim" || length(nonest))
                    attr(fit, "non-estim") <- nonest
            }
        }
    }
    if (!se.fit && interval == "none") {
        return(if (is.null(na_act)) fit else stats::napredict(na_act, fit))
    }
    V <- vcov(object)
    V[is.na(V)] <- 0
    se <- if (ncol(D) == 0L) rep(0, nrow(D)) else sqrt(quadratic_diagonal(D, V))
    names(se) <- names(fit)
    res <- fit
    if (interval != "none") {
        tq <- stats::qt((1 - level) / 2, object$df.residual)
        hw <- tq * if (interval == "confidence") se else sqrt(se^2 + object$sigma^2)
        res <- cbind(fit = fit, lwr = fit + hw, upr = fit - hw)
    }
    if (!is.null(na_act)) {
        res <- stats::napredict(na_act, res)
        se <- stats::napredict(na_act, se)
    }
    if (se.fit) {
        list(fit = res, se.fit = se, df = object$df.residual, residual.scale = object$sigma)
    } else {
        res
    }
}

## Rows of D whose predictions are not estimable from a rank-deficient fit:
## rows with a component, relative to their norm, in the null space of the
## design (as predict.lm() determines them).
dmfit_nonestimable <- function(object, D, tol) {
    p <- length(object$coefficients)
    k <- object$rank
    tR <- t(object$R)
    tR[cbind(seq_len(p), seq_len(p))[seq_len(p) > k, , drop = FALSE]] <- 1
    nb <- qr.Q(base::qr(tR))[, (k + 1L):p, drop = FALSE]
    B <- matrix(0, p, p - k)
    B[object$pivot, ] <- nb
    Xb <- D %*% B
    Xb_norm <- sqrt(rowSums(Xb^2))
    X_norm <- sqrt(quadratic_diagonal(D, diag(1, p)))
    nonest <- which(tol * X_norm <= Xb_norm)
    names(nonest) <- dimnames(D)[[1L]][nonest]
    nonest
}

## The offset of new data: offset() terms of the formula and the `offset`
## argument of the call, as predict.lm() evaluates them.
dmfit_new_offset <- function(object, newdata, na.action) {
    tt <- object$terms
    off <- NULL
    if (!is.null(tt) && !is.null(attr(tt, "offset"))) {
        mf <- stats::model.frame(stats::delete.response(tt), newdata,
                                 na.action = na.action, xlev = object$xlevels)
        off <- as.vector(stats::model.offset(mf))
    }
    if (!is.null(object$offset_call)) {
        o <- eval(object$offset_call, newdata, environment(tt))
        off <- if (is.null(off)) o else off + o
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
                rows = object$rows)
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
    if (x$compressed)
        cat(sprintf("\nFitted on %d design patterns of %d rows.\n", x$patterns, x$rows))
    if (length(x$aliased) == 0L) {
        cat("\nNo Coefficients\n")
    } else {
        nsing <- sum(x$aliased)
        if (nsing > 0L) {
            cat("\nCoefficients: (", nsing, " not defined because of singularities)\n", sep = "")
            cn <- names(x$aliased)
            coefs <- matrix(NA_real_, length(x$aliased), 4L, dimnames = list(cn, colnames(x$coefficients)))
            coefs[!x$aliased, ] <- x$coefficients
        } else {
            cat("\nCoefficients:\n")
            coefs <- x$coefficients
        }
        stats::printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                            na.print = "NA", ...)
    }
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

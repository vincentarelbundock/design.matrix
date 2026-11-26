## Generalized linear models fitted by IRLS on a design_matrix, without the
## dense matrix (remote data: remote_glm.R).
##
## The iterations are those of glm.fit(), step for step: same starting
## values (family$initialize, mustart, etastart, start), same weighted least
## squares problem at each iteration (solved with qr() on the scaled design,
## lm()'s pivoting with glm.fit()'s tolerance), same step halving, the same
## convergence rule and the same warnings. Only the products change:
## X'WX, X'Wz and X beta come from the design's kernels.
##
## Compression. Rows that share a design pattern and an offset share eta
## and mu. For a family whose deviance is a quasi-deviance (every standard
## family: dev(y, mu, w) = w (a(y) + y h(mu) + k(mu))), the working weights
## and responses of a pattern at fixed eta depend on its rows only through
## the sum of the prior weights W and the weighted mean response ybar, so
## the iterations on the patterns, with weights W and response ybar, are
## those on the rows; and the deviance of the rows is that of the patterns
## plus C = sum_i dev(y_i, ybar_g(i), w_i), which does not depend on the
## coefficients. The first iteration starts from per-row values (mustart
## depends on y), so its working weights and responses are summed by
## pattern from the rows. The AIC is computed once from the rows, and the
## dispersion of quasi families from the within-pattern sums of squares.
##
## The loop (dmglm_loop()) runs on a "backend": a list of closures
##   wls():      the weighted least squares fit at the current eta
##               (list(coefficients, rank, pivot, R)), or NULL when no
##               observation is informative;
##   eval(beta): move to beta: list(dev, valid);
##   mu_range(): range of the current fitted means;
##   state:      an environment (eta, mu, w2, last).
## dmglm_units() is the backend of rows or patterns held in memory; the
## remote backends are in remote_glm.R.

#' Generalized linear models without the model matrix
#'
#' `design_glm()` fits a generalized linear model by iteratively reweighted
#' least squares (IRLS), like [stats::glm()], without forming the dense
#' model matrix. It takes a formula and a data frame, a [design_matrix()]
#' and a response `y` (like [stats::glm.fit()]), or remote data (see
#' *Remote data* below).
#'
#' **Method.** The iterations are those of [glm.fit()]: the same starting
#' values (from `family$initialize`, `mustart`, `etastart` or `start`), the
#' same weighted least squares problem at each iteration, solved with
#' [qr()][design_qr] on the scaled design (columns pivoted, and aliased
#' columns detected, as in `glm.fit()`, with its tolerance
#' `min(1e-7, epsilon / 1000)`), the same step halving when the deviance is
#' not finite or the fitted values are invalid, the same convergence rule
#' `|dev - dev_old| / (|dev| + 0.1) < epsilon`, and the same warnings
#' (non-convergence, boundary, fitted probabilities numerically 0 or 1).
#' Results agree with those of `glm()` up to floating point rounding: the
#' number of iterations is the same, and the coefficients, covariance
#' matrix, deviances and AIC agree to many more digits than `epsilon`.
#'
#' **Compression.** Rows with the same design pattern and the same offset
#' have the same fitted value. For the families of [stats::family] (and
#' [quasi()], and `MASS::negative.binomial()`), the deviance has the form
#' `w * (a(y) + y * h(mu) + k(mu))`, so at fixed coefficients a pattern
#' contributes to each iteration only through the sum of its prior weights
#' and its weighted mean response, and to the deviance through a constant
#' that is computed once. The fit then runs on the patterns and gives the
#' results of the fit on every row, with the same iterations (the first
#' iteration, which starts from per-row values, is summed from the rows).
#' The AIC is computed from the rows, and the dispersion of the families
#' that estimate it from the within-pattern sums of squares. `compress`
#' works as in [design_lm()]; other families are never compressed (and
#' `compress = TRUE` is an error for them).
#'
#' **Remote data.** `data` can be a DuckDB connection (with `table` or
#' `sql`) or Parquet paths, and the first argument a remote design (see
#' [design_remote]); `weights` and `offset` then name columns (`offset()`
#' terms of the formula are also used when they have a SQL translation).
#' When the rows have few design patterns, one `GROUP BY` returns the sums
#' that the compressed fit needs, and the iterations run in memory on the
#' patterns. Otherwise each iteration computes `X'WX` and `X'Wz` in the
#' database (SQL engine: the linear predictor from the current coefficients,
#' the working weights and responses from the family, all in SQL), or in
#' batches of rows read into R (chunked engine). See the section *GLMs* of
#' [design_remote].
#'
#' **What is kept.** As in [design_lm()]: the coefficients, the factor `R`
#' of the last weighted least squares problem (which gives the covariance
#' matrix, as in [summary.glm()]), the deviances, AIC, degrees of freedom,
#' the number of iterations, and, with `model = TRUE`, the model frame and
#' design, from which [residuals()], [fitted()] and [predict()] are
#' recomputed.
#'
#' **Methods.** [coef()], [vcov()] (`NA` for aliased coefficients, scaled
#' by the dispersion), [summary()] (coefficient table, dispersion, null and
#' residual deviance, AIC, iterations), [predict()] (`type = "link"` or
#' `"response"`, with `se.fit`), [residuals()] (types `"deviance"`,
#' `"pearson"`, `"working"` and `"response"`), [fitted()], [deviance()],
#' [df.residual()], [nobs()], [logLik()] (hence [AIC()] and [BIC()]),
#' [family()], [formula()], [terms()], [model.frame()] and `print()`.
#'
#' **Differences with `glm()`.** The result is not a `glm` object: it has
#' no `qr`, `effects`, `residuals`, `fitted.values` or `weights`
#' components, so functions that expect those (such as `anova.glm()` or
#' `plot.lm()`) do not apply. `method`, `x`, `singular.ok = FALSE` and
#' `type = "terms"` in `predict()` are not supported. The null deviance of a
#' fit with an offset and an intercept is refitted, as `glm()` does (also
#' for a design matrix, where `glm.fit()` would not). Extra arguments do not
#' go to [glm.control()]: give `control = list(maxit = 50)`.
#'
#' @param formula A model formula, as in [stats::glm()], a [design_matrix()]
#'   (then give the response as `y`) or a remote design.
#' @param family A family object, function or name, as in [stats::glm()].
#' @param data,weights,subset,na.action,start,etastart,mustart,offset,contrasts
#'   As in [stats::glm()]. With a design matrix, `weights`, `offset`,
#'   `etastart` and `mustart` are numeric vectors with one value per row, and
#'   `data`, `subset`, `na.action` and `contrasts` do not apply. For remote
#'   data, `weights` and `offset` name columns (strings or bare names), and
#'   `subset`, `etastart` and `mustart` are not supported.
#' @param control A list of arguments of [stats::glm.control()].
#' @param compress `"auto"`, `TRUE` or `FALSE`: fit on the design patterns
#'   (see Details).
#' @param model Keep the model frame and its design matrix, for
#'   [residuals()], [fitted()] and [predict()] without `newdata`.
#' @param y With a design matrix as first argument, the response: a numeric
#'   vector, a factor, or a two-column matrix of successes and failures (as
#'   in [glm.fit()]). With a remote design, the name of the response column.
#' @param ... For remote data: `table`, `sql`, `engine` and `chunk_size`; see
#'   [design_remote]. Otherwise not used.
#' @return An object of class `"design_glm"`: a list with, among others,
#'   `coefficients`, `deviance`, `null.deviance`, `aic`, `iter`, `converged`,
#'   `boundary`, `df.residual`, `df.null`, `nobs`, `rank`, `pivot`, `R`,
#'   `family`, `compressed`, `patterns`, `rows`, `spec`, `call` and `terms`.
#' @examples
#' set.seed(1)
#' n <- 2000
#' d <- data.frame(x = rnorm(n), f = factor(sample(letters[1:4], n, TRUE)),
#'                 g = factor(sample(1:3, n, TRUE)))
#' d$y <- rbinom(n, 1, plogis(d$x + as.integer(d$f) / 4))
#' fit <- design_glm(y ~ x * f, binomial, d)
#' all.equal(coef(fit), coef(glm(y ~ x * f, binomial, d)))
#' all.equal(vcov(fit), vcov(glm(y ~ x * f, binomial, d)))
#' summary(fit)
#'
#' ## factors only: fitted on 12 patterns instead of 2000 rows
#' d$k <- rpois(n, exp(0.5 + as.integer(d$g) / 3))
#' fit <- design_glm(k ~ f * g, poisson, d)
#' fit$patterns
#' all.equal(deviance(fit), deviance(glm(k ~ f * g, poisson, d)))
#'
#' predict(fit, newdata = d[1:3, ], type = "response")
#' @export
design_glm <- function(formula, family = stats::gaussian, data, weights, na.action,
                       start = NULL, offset, control = list(),
                       contrasts = NULL, compress = "auto", model = TRUE, y = NULL, ...) {
    family <- dmglm_family(family, parent.frame())
    dmglm_check_family(family)
    control <- do.call(stats::glm.control, as.list(control))
    dmfit_check_compress(compress)
    checkmate::assert_flag(model)
    cl <- match.call()
    if (inherits(formula, c("design_matrix", "dm_remote"))) {
        if (!missing(data) || !missing(na.action) || !is.null(contrasts))
            stop("With a design matrix, give the response as `y`; `data`, ",
                 "`na.action` and `contrasts` do not apply.", call. = FALSE)
        if (...length()) stop("Unknown arguments: ", paste(names(list(...)), collapse = ", "), call. = FALSE)
        w <- if (missing(weights)) NULL else weights
        off <- if (missing(offset)) NULL else offset
        if (inherits(formula, "dm_remote")) {
            fit <- dm_remote_glm_design(formula, y, w, off, family, start, control, compress)
        } else {
            fit <- dmglm_design(formula, y, w, off, start, family, control, compress)
        }
        fit$call <- cl
        return(fit)
    }
    if (!is.null(y)) stop("`y` applies only when the first argument is a design matrix.", call. = FALSE)
    ## remote data (see remote_glm.R)
    if (!missing(data) && dm_is_remote_source(data)) {
        na_fail <- !missing(na.action) && dmglm_na(na.action) == "fail"
        env <- parent.frame()
        col_arg <- function(e, what) {
            ## a bare column name, a string, or a variable holding a string
            if (is.name(e)) {
                v <- tryCatch(eval(e, env), error = function(err) err)
                e <- if (is.null(v)) NULL else if (is.character(v) && length(v) == 1L) v else as.character(e)
            } else if (!is.null(e)) {
                e <- eval(e, env)
            }
            if (!is.null(e) && !checkmate::test_string(e))
                stop(sprintf("For remote data, `%s` must name a column.", what), call. = FALSE)
            e
        }
        w <- col_arg(if (missing(weights)) NULL else substitute(weights), "weights")
        off <- col_arg(if (missing(offset)) NULL else substitute(offset), "offset")
        return(dm_remote_glm(formula, data, family, w, off, start, contrasts, compress, control,
                             cl, list(...), na_fail = na_fail))
    }
    if (...length()) stop("Unknown arguments: ", paste(names(list(...)), collapse = ", "), call. = FALSE)
    ## the model frame, as glm() builds it
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "weights", "na.action", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    if (!missing(na.action)) dmglm_na(na.action)
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(design.matrix::design_frame)
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    Y <- dmglm_response(mf)
    w <- as.vector(stats::model.weights(mf))
    if (!is.null(w) && !is.numeric(w)) stop("'weights' must be a numeric vector", call. = FALSE)
    if (!is.null(w) && any(w < 0)) stop("negative weights not allowed", call. = FALSE)
    off <- as.vector(stats::model.offset(mf))
    if (!is.null(off) && length(off) != NROW(Y))
        stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                      length(off), NROW(Y)), call. = FALSE)
    X <- design_matrix(mt, mf, contrasts.arg = contrasts)
    fit <- dmglm_fit(X, Y, w, start, off, family, control,
                     intercept = attr(mt, "intercept") > 0L, compress = compress)
    fit$call <- cl
    fit$terms <- mt
    fit$formula <- stats::formula(mt)
    fit$contrasts <- X$contrasts
    fit$xlevels <- stats::.getXlevels(mt, mf)
    fit$na.action <- attr(mf, "na.action")
    fit$offset_call <- cl$offset
    if (model) {
        fit$model <- mf
    } else {
        fit$data <- NULL
    }
    fit
}

## na.action: na.omit, na.exclude or na.fail (those remote data supports).
dmglm_na <- function(na.action) {
    for (nm in c("omit", "exclude", "fail")) {
        f <- get(paste0("na.", nm), envir = asNamespace("stats"))
        if (identical(na.action, f) || identical(na.action, paste0("na.", nm))) return(nm)
    }
    stop("`na.action` must be na.omit, na.exclude or na.fail.", call. = FALSE)
}

## The family, as glm() resolves it.
dmglm_family <- function(family, env) {
    if (is.character(family)) family <- get(family, mode = "function", envir = env)
    if (is.function(family)) family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized", call. = FALSE)
    }
    family
}

## The response of a model frame as glm() passes it to glm.fit() (a vector,
## a factor or a two-column matrix), without names.
dmglm_response <- function(mf) {
    if (attr(attr(mf, "terms"), "response") == 0L) stop("The formula has no response.", call. = FALSE)
    dmglm_y(mf[[1L]])
}

dmglm_y <- function(v) {
    if (length(dim(v)) == 1L) dim(v) <- NULL
    if (is.matrix(v)) {
        if (!is.numeric(v) && !is.logical(v)) stop("The response must be numeric.", call. = FALSE)
        if (ncol(v) == 1L) return(as.double(v))
        storage.mode(v) <- "double"
        rownames(v) <- NULL
        return(v)
    }
    if (is.factor(v)) {
        names(v) <- NULL
        return(v)
    }
    if (!is.numeric(v) && !is.logical(v))
        stop("The response must be numeric, logical or a factor.", call. = FALSE)
    attributes(v) <- NULL
    storage.mode(v) <- "double"
    v
}

## design_glm() on a design matrix x and a response y.
dmglm_design <- function(x, y, weights, offset, start, family, control, compress) {
    checkmate::assert_class(x, "design_matrix")
    if (is.null(y)) stop("Give the response as `y`.", call. = FALSE)
    y <- dmglm_y(y)
    n <- nrow(x)
    if (NROW(y) != n) stop(sprintf("The response must have %d values.", n), call. = FALSE)
    if (anyNA(y)) stop("NA/NaN/Inf in 'y'", call. = FALSE)
    if (!is.null(weights)) {
        checkmate::assert_numeric(weights, len = n)
        if (anyNA(weights) || any(weights < 0)) stop("missing or negative weights not allowed", call. = FALSE)
    }
    if (!is.null(offset)) {
        checkmate::assert_numeric(offset, len = n)
        if (!all(is.finite(offset))) stop("NA/NaN/Inf in 'offset'", call. = FALSE)
    }
    tt <- x$training_terms
    intercept <- if (!is.null(tt)) attr(tt, "intercept") > 0L else "(Intercept)" %in% colnames(x)
    dmglm_fit(x, y, weights, start, offset, family, control, intercept, compress)
}

## family$initialize, evaluated as glm.fit() evaluates it: it can change y
## (proportions from a two-column response, zero for rows with zero
## weight), the weights (times the number of trials) and sets n and mustart.
dmglm_initialize <- function(family, y, weights, nobs, start) {
    n <- NULL
    etastart <- mustart <- NULL
    eval(family$initialize)
    if (is.null(n)) n <- rep.int(1, nobs)
    list(y = as.double(y), weights = as.double(weights), n = n, mustart = mustart)
}

## Families whose deviance is a quasi-deviance (see the top of the file),
## for which the fit on patterns is exact.
dmglm_compressible <- function(family) {
    family$family %in% c("gaussian", "binomial", "poisson", "Gamma", "inverse.gaussian",
                         "quasibinomial", "quasipoisson", "quasi") ||
        startsWith(family$family, "Negative Binomial(")
}

## Families with a fixed dispersion of 1 (summary.glm()).
dmglm_fixed_dispersion <- function(family) family$family %in% c("poisson", "binomial")

## Groups of rows with the same design pattern and offset, or NULL.
dmglm_groups <- function(x, offset, family, compress) {
    if (isFALSE(compress) || ncol(x) == 0L) return(NULL)
    if (!dmglm_compressible(family)) {
        if (isTRUE(compress))
            stop(sprintf("`compress = TRUE`: the %s family cannot be compressed exactly.", family$family),
                 call. = FALSE)
        return(NULL)
    }
    g <- dmfit_groups(x, compress)
    if (is.null(g)) return(NULL)
    if (any(offset != offset[1L])) {
        G <- length(g$first)
        key <- (match(offset, unique(offset)) - 1) * G + g$map
        first <- which(!duplicated(key))
        n <- length(key)
        if (length(first) == n || (!isTRUE(compress) && length(first) > 0.5 * n)) return(NULL)
        g <- list(map = match(key, key[first]), first = first)
    }
    g
}

dmglm_rowsum <- function(v, map, G) {
    out <- drop(rowsum(as.double(v), map, reorder = TRUE))
    if (length(out) != G) stop("internal error: groups", call. = FALSE)
    unname(out)
}

## The working weights and responses of the first iteration, from the
## per-row starting values, summed by pattern (with glm.fit()'s checks).
dmglm_first <- function(eta, mu, y, weights, offset, family, map, G) {
    good <- weights > 0
    varmu <- family$variance(mu)
    if (anyNA(varmu[good])) stop("NAs in V(mu)", call. = FALSE)
    if (any(varmu[good] == 0)) stop("0s in V(mu)", call. = FALSE)
    me <- family$mu.eta(eta)
    if (any(is.na(me[good]))) stop("NAs in d(mu)/d(eta)", call. = FALSE)
    good <- (weights > 0) & (me != 0)
    if (!any(good)) return(list(none = TRUE))
    w2 <- numeric(length(y))
    wz <- numeric(length(y))
    w2[good] <- (weights * me^2 / varmu)[good]
    wz[good] <- w2[good] * ((eta - offset) + (y - mu) / me)[good]
    W1 <- dmglm_rowsum(w2, map, G)
    Z1 <- dmglm_rowsum(wz, map, G)
    list(w2 = W1, z = ifelse(W1 > 0, Z1 / W1, 0))
}

## Weighted least squares on a design (rows or patterns): w2 are the
## squared IRLS weights (0 for rows left out), z the working response.
dmglm_wls <- function(x, tol) {
    force(x)
    function(w2, z) {
        p <- ncol(x)
        sw <- sqrt(w2)
        if (!all(is.finite(sw)) || !all(is.finite(z)))
            return(list(coefficients = rep(NaN, p), rank = 0L, pivot = seq_len(p), R = matrix(0, p, p)))
        q <- qr(sw * x, tol = tol)
        b <- dmqr_coef(q, matrix(sw * z, ncol = 1L))
        list(coefficients = dmqr_full(q, b)[, 1L], rank = q$rank, pivot = q$pivot, R = unname(q$R))
    }
}

## The same for an intercept alone (the null model).
dmglm_wls_mean <- function(w2, z) {
    s <- sum(w2)
    if (!is.finite(s) || !all(is.finite(z)))
        return(list(coefficients = NaN, rank = 0L, pivot = 1L, R = matrix(0, 1L, 1L)))
    list(coefficients = if (s > 0) sum(w2 * z) / s else 0, rank = if (s > 0) 1L else 0L,
         pivot = 1L, R = matrix(sqrt(s), 1L, 1L))
}

## The backend of units (rows, or patterns) in memory: linear predictor
## `lin(beta)` (without offset), solver `wls(w2, z)`, response y, prior
## weights W, offsets, and the constant C added to the deviance. `first`:
## the working weights and responses of the first iteration (from rows),
## or `eta`: the starting linear predictor of the units.
dmglm_units <- function(lin, wls, y, W, off, C, family, first = NULL, eta = NULL) {
    st <- new.env(parent = emptyenv())
    st$first <- first
    if (!is.null(eta)) {
        st$eta <- eta
        st$mu <- family$linkinv(eta)
    }
    valideta <- family$valideta %||% function(eta) TRUE
    validmu <- family$validmu %||% function(mu) TRUE
    list(
        state = st,
        wls = function() {
            if (!is.null(st$first)) {
                f <- st$first
                st$first <- NULL
                if (isTRUE(f$none)) return(NULL)
                st$last <- "rows"
                st$w2 <- f$w2
                return(wls(f$w2, f$z))
            }
            eta <- st$eta
            mu <- st$mu
            good <- W > 0
            varmu <- family$variance(mu)
            if (anyNA(varmu[good])) stop("NAs in V(mu)", call. = FALSE)
            if (any(varmu[good] == 0)) stop("0s in V(mu)", call. = FALSE)
            me <- family$mu.eta(eta)
            if (any(is.na(me[good]))) stop("NAs in d(mu)/d(eta)", call. = FALSE)
            good <- (W > 0) & (me != 0)
            if (!any(good)) return(NULL)
            w2 <- numeric(length(W))
            z <- numeric(length(W))
            w2[good] <- (W * me^2 / varmu)[good]
            z[good] <- ((eta - off) + (y - mu) / me)[good]
            st$last <- "units"
            st$w2 <- w2
            wls(w2, z)
        },
        eval = function(beta) {
            eta <- lin(beta) + off
            mu <- family$linkinv(eta)
            st$eta <- eta
            st$mu <- mu
            list(dev = sum(family$dev.resids(y, mu, W)) + C,
                 valid = isTRUE(valideta(eta) && validmu(mu)))
        },
        mu_range = function() range(st$mu))
}

## The IRLS loop of glm.fit(), on a backend (see the top of the file).
## `nobs`: number of rows; `devold`: deviance at the starting values;
## `coefold`: the starting coefficients, if given.
dmglm_loop <- function(B, family, control, nobs, devold, coefold = NULL) {
    conv <- boundary <- FALSE
    coef <- fit <- NULL
    dev <- devold
    maxit <- control$maxit
    iter <- 0L
    for (iter in seq_len(maxit)) {
        f <- B$wls()
        if (is.null(f)) {
            warning(gettextf("no observations informative at iteration %d", iter), call. = FALSE)
            break
        }
        fit <- f
        if (any(!is.finite(f$coefficients))) {
            warning(gettextf("non-finite coefficients at iteration %d", iter), call. = FALSE)
            break
        }
        if (nobs < f$rank)
            stop(sprintf(ngettext(nobs, "X matrix has rank %d, but only %d observation",
                                  "X matrix has rank %d, but only %d observations"), f$rank, nobs),
                 call. = FALSE)
        start <- f$coefficients
        ev <- B$eval(start)
        dev <- ev$dev
        if (control$trace) cat("Deviance = ", dev, " Iterations - ", iter, "\n", sep = "")
        boundary <- FALSE
        if (!is.finite(dev)) {
            if (is.null(coefold))
                stop("no valid set of coefficients has been found: please supply starting values", call. = FALSE)
            warning("step size truncated due to divergence", call. = FALSE)
            ii <- 1
            while (!is.finite(dev)) {
                if (ii > maxit) stop("inner loop 1; cannot correct step size", call. = FALSE)
                ii <- ii + 1
                start <- (start + coefold) / 2
                ev <- B$eval(start)
                dev <- ev$dev
            }
            boundary <- TRUE
            if (control$trace) cat("Step halved: new deviance = ", dev, "\n", sep = "")
        }
        if (!ev$valid) {
            if (is.null(coefold))
                stop("no valid set of coefficients has been found: please supply starting values", call. = FALSE)
            warning("step size truncated: out of bounds", call. = FALSE)
            ii <- 1
            while (!ev$valid) {
                if (ii > maxit) stop("inner loop 2; cannot correct step size", call. = FALSE)
                ii <- ii + 1
                start <- (start + coefold) / 2
                ev <- B$eval(start)
            }
            boundary <- TRUE
            dev <- ev$dev
            if (control$trace) cat("Step halved: new deviance = ", dev, "\n", sep = "")
        }
        if (abs(dev - devold) / (0.1 + abs(dev)) < control$epsilon) {
            conv <- TRUE
            coef <- start
            break
        }
        devold <- dev
        coef <- coefold <- start
    }
    if (is.null(coef)) stop("The fit failed: no valid coefficients were found.", call. = FALSE)
    if (!conv) warning("glm.fit: algorithm did not converge", call. = FALSE)
    if (boundary) warning("glm.fit: algorithm stopped at boundary value", call. = FALSE)
    eps <- 10 * .Machine$double.eps
    r <- B$mu_range()
    if (family$family == "binomial" && (r[2L] > 1 - eps || r[1L] < eps))
        warning("glm.fit: fitted probabilities numerically 0 or 1 occurred", call. = FALSE)
    if (family$family == "poisson" && r[1L] < eps)
        warning("glm.fit: fitted rates numerically 0 occurred", call. = FALSE)
    p <- length(coef)
    if (fit$rank < p) coef[fit$pivot[seq.int(fit$rank + 1L, p)]] <- NA
    list(coef = coef, fit = fit, iter = iter, conv = conv, boundary = boundary, dev = dev)
}

## The null deviance of a fit with an offset and an intercept, refitted as
## glm() does, on units: starting from the fitted values `mu` of the model.
dmglm_null_refit <- function(y, W, off, C, family, control, nobs, mu) {
    eta <- family$linkfun(mu)
    B <- dmglm_units(function(b) rep(b[1L], length(W)), dmglm_wls_mean, y, W, off, C, family, eta = eta)
    valideta <- family$valideta %||% function(eta) TRUE
    validmu <- family$validmu %||% function(mu) TRUE
    if (!(validmu(B$state$mu) && valideta(eta)))
        stop("cannot find valid starting values: please specify some", call. = FALSE)
    devold <- sum(family$dev.resids(y, B$state$mu, W)) + C
    res <- dmglm_loop(B, family, control, nobs, devold)
    if (!res$conv)
        warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?", call. = FALSE)
    res$dev
}

## The fit in memory: design x, response y (vector, factor or two-column
## matrix), prior weights, starting values, offset (NULL for none).
dmglm_fit <- function(x, y, weights, start, offset, family, control, intercept, compress) {
    nobs <- NROW(y)
    if (nrow(x) != nobs) stop("The design and the response have different numbers of rows.", call. = FALSE)
    if (nobs == 0L) stop("0 (non-NA) cases", call. = FALSE)
    nvars <- ncol(x)
    has_offset <- !is.null(offset)
    if (is.null(weights)) weights <- rep.int(1, nobs)
    offset <- if (has_offset) as.double(offset) else rep.int(0, nobs)
    if (!is.function(family$variance) || !is.function(family$linkinv))
        stop("'family' argument seems not to be a valid family object", call. = FALSE)
    linkinv <- family$linkinv
    dev.resids <- family$dev.resids
    valideta <- family$valideta %||% function(eta) TRUE
    validmu <- family$validmu %||% function(mu) TRUE
    ini <- dmglm_initialize(family, y, weights, nobs, start)
    y <- ini$y
    weights <- ini$weights
    n <- ini$n
    tol <- min(1e-07, control$epsilon / 1000)

    if (nvars == 0L) {
        ## the empty model: no iterations, as in glm.fit()
        eta <- offset
        if (!valideta(eta)) stop("invalid linear predictor values in empty model", call. = FALSE)
        mu <- linkinv(eta)
        if (!validmu(mu)) stop("invalid fitted means in empty model", call. = FALSE)
        dev <- sum(dev.resids(y, mu, weights))
        w2 <- weights * family$mu.eta(eta)^2 / family$variance(mu)
        res <- list(coef = numeric(0), fit = list(rank = 0L, pivot = integer(0), R = matrix(0, 0, 0)),
                    iter = 0L, conv = TRUE, boundary = TRUE, dev = dev)
        units <- list(y = y, W = weights, off = offset, C = 0, map = NULL, eta = eta, mu = mu, w2 = w2,
                      last = "units")
    } else {
        coefold <- NULL
        eta <- if (!is.null(start)) {
            if (length(start) != nvars)
                stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
                              nvars, paste(deparse(colnames(x)), collapse = ", ")), call. = FALSE)
            coefold <- start
            offset + drop(x %*% start)
        } else {
            family$linkfun(ini$mustart)
        }
        mu <- linkinv(eta)
        if (!(validmu(mu) && valideta(eta)))
            stop("cannot find valid starting values: please specify some", call. = FALSE)
        devold <- sum(dev.resids(y, mu, weights))

        g <- dmglm_groups(x, offset, family, compress)
        if (!is.null(g)) {
            G <- length(g$first)
            map <- g$map
            W <- dmglm_rowsum(weights, map, G)
            ybar <- dmglm_rowsum(weights * y, map, G) / W
            ybar[W == 0] <- 0
            C <- sum(dev.resids(y, ybar[map], weights))
            ## sums that overflow: fit on the rows
            if (!all(is.finite(W), is.finite(ybar), is.finite(C))) g <- NULL
        }
        if (!is.null(g)) {
            xs <- x[g$first, ]
            offs <- offset[g$first]
            first <- dmglm_first(eta, mu, y, weights, offset, family, map, G)
            B <- dmglm_units(function(b) drop(xs %*% b), dmglm_wls(xs, tol), ybar, W, offs, C, family,
                             first = first)
        } else {
            map <- NULL
            xs <- x
            W <- weights
            ybar <- y
            offs <- offset
            C <- 0
            B <- dmglm_units(function(b) drop(x %*% b), dmglm_wls(x, tol), y, weights, offset, 0, family,
                             eta = eta)
        }
        res <- dmglm_loop(B, family, control, nobs, devold, coefold)
        st <- B$state
        units <- list(y = ybar, W = W, off = offs, C = C, map = map, eta = st$eta, mu = st$mu,
                      w2 = st$w2, last = st$last)
        if (identical(st$last, "rows")) {
            ## only the first iteration ran: its working weights vary within
            ## patterns (for the dispersion)
            me0 <- family$mu.eta(eta)
            good0 <- weights > 0 & me0 != 0
            units$w2_rows <- ifelse(good0, weights * me0^2 / family$variance(mu), 0)
        }
    }

    ## deviances, AIC, dispersion
    W <- units$W
    rank <- res$fit$rank
    n.ok <- nobs - sum(weights == 0)
    wtdmu <- if (intercept) sum(weights * y) / sum(weights) else linkinv(units$off)
    nulldev <- sum(dev.resids(units$y, wtdmu, W)) + units$C
    if (has_offset && intercept)
        nulldev <- dmglm_null_refit(units$y, W, units$off, units$C, family, control, nobs, units$mu)
    mu_rows <- if (is.null(units$map)) units$mu else units$mu[units$map]
    aic <- family$aic(y, n, mu_rows, weights, res$dev) + 2 * rank
    disp <- dmglm_pearson(units, y, weights, family)

    out <- list(coefficients = res$coef, rank = rank, pivot = res$fit$pivot, R = res$fit$R, tol = tol,
                deviance = res$dev, null.deviance = nulldev, aic = aic, iter = res$iter,
                converged = res$conv, boundary = res$boundary,
                df.residual = n.ok - rank, df.null = n.ok - as.integer(intercept), nobs = n.ok,
                rows = nobs, pearson = disp$pearson, zero_weights = disp$zero, family = family,
                intercept = intercept, prior_weighted = any(weights != 1),
                compressed = !is.null(units$map), patterns = length(W), control = control,
                spec = dmfit_spec(x))
    names(out$coefficients) <- colnames(x)
    out$data <- list(design = x, y = y, weights = weights, offset = if (has_offset) offset)
    class(out) <- "design_glm"
    out
}

## sum(w * r^2) over rows with positive working weights, as summary.glm()
## estimates the dispersion: working weights of the last weighted least
## squares problem, working residuals at the final fitted values. On
## patterns, sum_i w_i (y_i - mu)^2 = within + W (ybar - mu)^2.
dmglm_pearson <- function(units, y, weights, family) {
    me <- family$mu.eta(units$eta)
    if (!is.null(units$w2_rows)) {
        ## working weights of the rows (only the first iteration ran)
        map <- units$map
        w2 <- units$w2_rows
        r2 <- (y - units$mu[map])^2 / me[map]^2
        pos <- w2 > 0
        return(list(pearson = sum((w2 * r2)[pos]), zero = any(w2 == 0)))
    }
    w2 <- units$w2
    W <- units$W
    pos <- w2 > 0
    within <- 0
    if (!is.null(units$map) && !dmglm_fixed_dispersion(family)) {
        map <- units$map
        within <- dmglm_rowsum(weights * (y - units$y[map])^2, map, length(W))
    }
    ratio <- ifelse(W > 0, w2 / W, 0)
    s <- ratio * (within + W * (units$y - units$mu)^2) / me^2
    list(pearson = sum(s[pos]), zero = any(weights == 0) || any(w2[W > 0] == 0))
}

## The dispersion, as summary.glm() computes it (with its warning).
dmglm_dispersion <- function(object, warn = TRUE) {
    fam <- object$family
    if (!is.null(fam$dispersion) && !is.na(fam$dispersion)) return(list(value = fam$dispersion, est = FALSE))
    if (dmglm_fixed_dispersion(fam)) return(list(value = 1, est = FALSE))
    df.r <- object$df.residual
    if (df.r > 0) {
        if (warn && isTRUE(object$zero_weights))
            warning("observations with zero weight not used for calculating dispersion", call. = FALSE)
        return(list(value = object$pearson / df.r, est = TRUE))
    }
    list(value = NaN, est = TRUE)
}

## The linear predictor of the kept data (with the offset).
dmglm_eta <- function(object) {
    dat <- dmfit_data(object, "This method")
    eta <- dmfit_linpred(object, dat$design)
    if (!is.null(dat$offset)) eta <- eta + dat$offset
    eta
}

#' @export
print.design_glm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
    dmglm_print_method(x)
    if (length(x$coefficients)) {
        cat("\nCoefficients:\n")
        print.default(format(x$coefficients, digits = digits), print.gap = 2L, quote = FALSE)
    } else {
        cat("\nNo coefficients\n\n")
    }
    cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ", x$df.residual, "Residual\n")
    if (nzchar(mess <- stats::naprint(x$na.action))) cat("  (", mess, ")\n", sep = "")
    cat("Null Deviance:\t   ", format(signif(x$null.deviance, digits)),
        "\nResidual Deviance:", format(signif(x$deviance, digits)), "\tAIC:",
        format(signif(x$aic, digits)))
    cat("\n")
    invisible(x)
}

dmglm_print_method <- function(x) {
    if (x$compressed) {
        cat(sprintf("\nFitted on %s design patterns of %s rows.\n", dm_fmt_count(x$patterns), dm_fmt_count(x$rows)))
    } else if (isTRUE(x$remote)) {
        cat(sprintf("\nFitted from %s of %s rows.\n", x$method, dm_fmt_count(x$rows)))
    }
}

#' @export
coef.design_glm <- function(object, ...) object$coefficients

#' @export
vcov.design_glm <- function(object, complete = TRUE, ...) {
    V <- dmfit_vcov(object, dmglm_dispersion(object)$value)
    if (!complete) {
        keep <- !is.na(object$coefficients)
        V <- V[keep, keep, drop = FALSE]
    }
    V
}

#' @export
nobs.design_glm <- function(object, ...) object$nobs

#' @export
deviance.design_glm <- function(object, ...) object$deviance

#' @export
df.residual.design_glm <- function(object, ...) object$df.residual

#' @export
family.design_glm <- function(object, ...) object$family

#' @export
formula.design_glm <- function(x, ...) {
    if (!is.null(x$formula)) return(x$formula)
    if (is.null(x$terms)) return(NULL)
    stats::formula(x$terms)
}

#' @export
terms.design_glm <- function(x, ...) {
    if (is.null(x$terms)) stop("The fit has no terms: it was not fitted from a formula.", call. = FALSE)
    x$terms
}

#' @export
model.frame.design_glm <- function(formula, ...) {
    if (is.null(formula$model))
        stop("The fit has no model frame: refit `design_glm()` with `model = TRUE`.", call. = FALSE)
    formula$model
}

#' @export
logLik.design_glm <- function(object, ...) {
    fam <- object$family$family
    dispersion <- object$family$dispersion
    p <- object$rank
    if (!is.null(dispersion)) {
        if (is.na(dispersion)) p <- p + 1L
    } else if (fam %in% c("gaussian", "Gamma", "inverse.gaussian")) {
        p <- p + 1L
    }
    structure(p - object$aic / 2, nobs = object$rows, df = p, class = "logLik")
}

#' @export
fitted.design_glm <- function(object, ...) {
    if (isTRUE(object$remote)) return(dm_remote_glm_fitted(object, "response"))
    mu <- object$family$linkinv(dmglm_eta(object))
    names(mu) <- dimnames(object$data$design)[[1L]]
    stats::napredict(object$na.action, mu)
}

#' @export
residuals.design_glm <- function(object, type = c("deviance", "pearson", "working", "response"), ...) {
    type <- match.arg(type)
    if (isTRUE(object$remote)) return(dm_remote_glm_residuals(object, type))
    dat <- dmfit_data(object, "residuals()")
    fam <- object$family
    eta <- dmglm_eta(object)
    mu <- fam$linkinv(eta)
    y <- dat$y
    wts <- dat$weights
    res <- switch(type,
                  deviance = if (object$df.residual > 0) {
                      d.res <- sqrt(pmax(fam$dev.resids(y, mu, wts), 0))
                      ifelse(y > mu, d.res, -d.res)
                  } else rep.int(0, length(mu)),
                  pearson = (y - mu) * sqrt(wts) / sqrt(fam$variance(mu)),
                  working = (y - mu) / fam$mu.eta(eta),
                  response = y - mu)
    names(res) <- dimnames(dat$design)[[1L]]
    stats::naresid(object$na.action, res)
}

#' @param object A `"design_glm"` object.
#' @param newdata A data frame of new observations, or remote data (a DuckDB
#'   connection with `table` or `sql`, or Parquet paths). Omitted,
#'   predictions are those of the fit's rows (which need `model = TRUE`).
#' @param type `"link"` (the linear predictor) or `"response"` (the fitted
#'   means).
#' @rdname design_glm
#' @export
predict.design_glm <- function(object, newdata = NULL, type = c("link", "response"),
                               na.action = stats::na.pass, ...) {
    type <- match.arg(type)
    if (!is.null(newdata) && dm_is_remote_source(newdata))
        return(dm_remote_glm_predict(object, newdata, type, na.action, list(...)))
    if (...length()) stop("Unknown arguments: ", paste(names(list(...)), collapse = ", "), call. = FALSE)
    if (is.null(newdata) && isTRUE(object$remote)) return(dm_remote_glm_fitted(object, type))
    lmo <- object
    class(lmo) <- "design_lm"
    p <- predict.design_lm(lmo, newdata, na.action = na.action)
    if (type == "response") object$family$linkinv(p) else p
}

#' @export
summary.design_glm <- function(object, dispersion = NULL, ...) {
    if (is.null(dispersion)) {
        d <- dmglm_dispersion(object)
        dispersion <- d$value
        est.disp <- d$est
    } else {
        est.disp <- FALSE
    }
    df.r <- object$df.residual
    aliased <- is.na(object$coefficients)
    p <- object$rank
    if (p > 0) {
        keep <- object$pivot[seq_len(p)]
        coef.p <- object$coefficients[keep]
        covmat.unscaled <- chol2inv(object$R[seq_len(p), seq_len(p), drop = FALSE])
        dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
        covmat <- dispersion * covmat.unscaled
        s.err <- sqrt(diag(covmat))
        tvalue <- coef.p / s.err
        dn <- c("Estimate", "Std. Error")
        if (!est.disp) {
            coef.table <- cbind(coef.p, s.err, tvalue, 2 * stats::pnorm(-abs(tvalue)))
            dimnames(coef.table) <- list(names(coef.p), c(dn, "z value", "Pr(>|z|)"))
        } else if (df.r > 0) {
            coef.table <- cbind(coef.p, s.err, tvalue, 2 * stats::pt(-abs(tvalue), df.r))
            dimnames(coef.table) <- list(names(coef.p), c(dn, "t value", "Pr(>|t|)"))
        } else {
            coef.table <- cbind(coef.p, NaN, NaN, NaN)
            dimnames(coef.table) <- list(names(coef.p), c(dn, "t value", "Pr(>|t|)"))
        }
    } else {
        coef.table <- matrix(NA_real_, 0L, 4L, dimnames = list(NULL, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")))
        covmat.unscaled <- covmat <- matrix(NA_real_, 0L, 0L)
    }
    ans <- c(object[c("call", "terms", "family", "deviance", "aic", "contrasts", "df.residual",
                      "null.deviance", "df.null", "iter", "na.action", "compressed", "patterns",
                      "rows", "remote", "method")],
             list(coefficients = coef.table, aliased = aliased, dispersion = dispersion,
                  df = c(object$rank, df.r, length(aliased)), cov.unscaled = covmat.unscaled,
                  cov.scaled = covmat))
    class(ans) <- "summary.design_glm"
    ans
}

#' @export
print.summary.design_glm <- function(x, digits = max(3L, getOption("digits") - 3L),
                                     signif.stars = getOption("show.signif.stars"), ...) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
    dmglm_print_method(x)
    if (length(x$aliased) == 0L) {
        cat("\nNo Coefficients\n")
    } else {
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
    cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ", format(x$dispersion),
        ")\n\n", apply(cbind(paste(format(c("Null", "Residual"), justify = "right"), "deviance:"),
                             format(unlist(x[c("null.deviance", "deviance")]), digits = max(5L, digits + 1L)),
                             " on", format(unlist(x[c("df.null", "df.residual")])), " degrees of freedom\n"),
                       1L, paste, collapse = " "), sep = "")
    if (nzchar(mess <- stats::naprint(x$na.action))) cat("  (", mess, ")\n", sep = "")
    cat("AIC: ", format(x$aic, digits = max(4L, digits + 1L)), "\n\n",
        "Number of Fisher Scoring iterations: ", x$iter, "\n\n", sep = "")
    invisible(x)
}

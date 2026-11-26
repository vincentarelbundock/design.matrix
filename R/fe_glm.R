## Generalized linear models with absorbed fixed effects:
## design_glm(y ~ x | fe1 + fe2, family, data).
##
## The IRLS iterations of glm.fit() (xtx_glm_loop(), glm.R), with the
## weighted least squares step of each iteration replaced by the estimator
## of linear models with fixed effects (xtx_fe_core(), fe.R): the working
## response z on the design, with the working weights, absorbing the FE.
## The parameters of the loop are the slopes and the FE values together
## (theta = c(b, fe_1, ..., fe_K)), so that its step halving, convergence
## rule (the deviance criterion of glm.fit()) and warnings are those of
## glm(); the linear predictor is X b + sum_k fe_k[code_k] + offset, a
## product with the design and lookups of the FE values (xtx_fe_fitted()).
## Everything runs on the primitives, so the same code fits data in memory
## and remote data.
##
## The reduced system of the FE changes with the working weights at each
## iteration; conjugate gradients start from the solution of the previous
## iteration (warm start), and the numbers of iterations are kept
## (fit$fe_iter).
##
## Before the fit, rows whose FE level the FE fit perfectly are dropped
## iteratively (fixest's fixef.rm = "perfect_fit"): singletons and, when the
## link maps the real line onto (0, Inf) or (0, 1) (fixest's rule:
## linkinv(-100) ~ 0, linkinv(100) ~ 1), levels whose response is all 0 or
## all 1, whose FE would be infinite.
##
## Covariance: the inverse of X'MWX from the last weighted least squares
## step (as glm()), times the dispersion, with fixest's small sample
## correction (n - 1) / (n - K), K the number of parameters with the FE
## (fixest's ssc() defaults); clustered by one variable as for linear
## models (scores w e (MX), corrected for FE not nested in the clusters).

## The weighted least squares steps of a GLM with FE (the `solver` of
## xtx_glm_fit()): parameters theta = c(b, fe_1, ..., fe_K), b the
## coefficients of the columns of x but the intercept (absorbed).
xtx_fe_solver <- function(x, codes, tol, opts) {
    p <- ncol(x)
    xcols <- setdiff(seq_len(p), which(colnames(x) == "(Intercept)"))
    px <- length(xcols)
    G <- vapply(codes, function(g) g$G, 1L)
    offs <- c(0L, cumsum(G))
    warm <- new.env(parent = emptyenv())
    lin <- function(theta) {
        Th <- as.matrix(theta)
        beta <- matrix(0, p, ncol(Th))
        beta[xcols, ] <- Th[seq_len(px), ]
        fe <- lapply(seq_along(G), function(k) Th[px + offs[k] + seq_len(G[k]), , drop = FALSE])
        if (!is.matrix(theta)) {
            beta <- beta[, 1L]
            fe <- lapply(fe, function(v) v[, 1L])
        }
        xtx_fe_fitted(x, beta, fe, codes)
    }
    wls <- function(w2, z) {
        ## (the residuals, for clustered standard errors, only for the last
        ## step: see we())
        est <- xtx_fe_core(x, list(z = z, w = w2, codes = codes), tol, opts, warm = warm, resid_final = FALSE)
        warm$w2 <- w2
        warm$z <- z
        b0 <- est$b
        b0[is.na(b0)] <- 0
        list(coefficients = c(b0, unlist(est$fe, use.names = FALSE)), rank = est$rank + est$fe_dof,
             aliased = c(is.na(est$b), rep(FALSE, sum(G))), est = est)
    }
    ## the weighted residuals of the last step, est its fit
    we <- function(est) {
        b0 <- est$b
        b0[is.na(b0)] <- 0
        ## (remotely f binds on the frame of z and w2, with which it combines)
        f <- lin(c(b0, unlist(est$fe, use.names = FALSE)))
        warm$w2 * (warm$z - f)
    }
    list(lin = lin, wls = wls, we = we, np = px + sum(G), warm = warm, xcols = xcols, G = G, offs = offs)
}

## Does the link give infinite FE for levels with an all-0 (all-1)
## response? fixest's rule.
xtx_fe_bounds <- function(family) {
    list(zero = isTRUE(abs(family$linkinv(-100)) < 1e-15),
         one = isTRUE(abs(family$linkinv(100) - 1) < 1e-15))
}

## Entry point from design_glm(): `call` is its matched call, `env` its
## parent frame.
xtx_fe_glm <- function(call, env, formula, family, data, start, control, contrasts, compress, model,
                      vcov, fixef.rm, dots) {
    k <- xtx_glm_kind(family)
    ## (arguments checked by design_glm(): checks.R)
    pf <- xtx_fe_parse(formula)
    clv <- xtx_fe_vcov_arg(vcov)
    fixef.rm <- xtx_fe_rm_arg(fixef.rm, glm = TRUE)
    ## the design of the rows with complete variables, and the columns of the
    ## response, weights and offset (remote_fit.R)
    base <- design_problem(call, env, pf$main, data, contrasts, dots, family = k, fe = pf, clv = clv)
    gc <- base$cols
    ## the columns of the GLM other than y and w go along with the rows kept
    flags <- gc[c("cbind", "weighted")]
    base$y <- gc$y
    base$w <- gc$w
    base$more <- Filter(Negate(is.null), gc[c("n", "off", "a", "b", "y0")])
    prob <- xtx_fe_prepare(base, pf, clv, fixef.rm, xtx_fe_bounds(family), wmsg = "negative weights not allowed")
    cols <- c(list(y = prob$cols$z, w = prob$cols$w), prob$more, flags)
    if (is.null(cols$y0)) cols["y0"] <- list(NULL)
    res <- xtx_fe_glm_estimate(prob$x, cols, prob$cols$codes, prob$cols$cl, family, control)
    fit <- res$fit
    est <- res$est
    x <- prob$x
    sh <- xtx_fe_fields(est, colnames(x), pf, prob$levels, clv)
    fit[names(sh$fields)] <- sh$fields
    fit$pivot <- fit$R <- NULL
    fit$nsingletons <- prob$nsingle
    fit$nperfect <- prob$nperfect
    fit$fixef.rm <- fixef.rm
    fit$fe_iter <- res$fe_iter
    fit[names(sh$inference)] <- sh$inference
    ## clustered: kept beside the iid covariance (vcov() chooses)
    if (!is.null(clv)) fit$vcov_cluster <- sh$vcov_cluster
    ## the rows of the fit: the design (not centered) and its columns
    fit$data <- list(design = x, y = cols$y, weights = cols$w, offset = cols$off, codes = prob$cols$codes)
    fit$solver_fit <- NULL
    fit$call <- call
    fit$terms <- base$terms
    fit$formula <- formula
    fit$contrasts <- base$contrasts
    fit$xlevels <- base$xlevels
    fit$spec <- xtx_fit_spec(base$x)
    fit$offset_call <- call$offset
    fit$env <- base$env
    fit$rows <- base$nall
    ## (remotely, na.exclude pads the results with every row not used:
    ## na_action.R)
    fit$na.action <- xtx_na_result(base, prob$rows, prob$x_rows)
    fit$method <- switch(xtx_engine(x), sql = "SQL iterations", NULL)
    ## (the data of a remote fit stay in the database: always kept)
    if (!model && is.null(fit$method)) fit$data <- NULL
    class(fit) <- c("xtx_glm_fe", "xtx_glm")
    fit
}

## The fit, on the design with numeric variables centered when that is
## exact (xtx_p_center(); coefficients, covariance and FE mapped back with
## xtx_fe_uncenter()), else on the design. Returns list(fit (of
## xtx_glm_fit()), est (slopes b, FE values fe, bread, cluster covariance,
## degrees of freedom), fe_iter).
xtx_fe_glm_estimate <- function(x, cols, codes, cl, family, control) {
    opts <- xtx_fe_opts()
    ## aliasing: fixest's criterion, with feglm()'s collin.tol = 1e-9
    tol <- 1e-9
    ctr <- if (isTRUE(getOption("xtx.fe_center", TRUE))) xtx_p_center(x) else NULL
    run <- function(xx, on) {
        cc <- cols
        for (nm in c("y", "w", "n", "off", "a", "b", "y0")) if (!is.null(cc[[nm]])) cc[[nm]] <- on(cc[[nm]])
        cd <- xtx_cols_map(codes, on)
        cl2 <- xtx_cols_map(cl, on)
        solver <- xtx_fe_solver(xx, cd, tol, opts)
        fit <- xtx_glm_fit(xx, cc, family, NULL, control, TRUE, FALSE, wmsg = "negative weights not allowed",
                         solver = solver)
        est <- fit$solver_fit$est
        ## (the scores from the absorption operator of the last step)
        if (!is.null(cl2)) est <- xtx_fe_cluster(est, cl2, solver$we(est))
        th <- fit$coefficients
        px <- length(solver$xcols)
        est$b <- th[seq_len(px)]
        est$fe <- lapply(seq_along(solver$G), function(j) th[px + solver$offs[j] + seq_len(solver$G[j])])
        est$we <- est$absorb <- NULL
        list(fit = fit, est = est, fe_iter = solver$warm$iter)
    }
    if (!is.null(ctr)) {
        xc <- ctr$x
        r <- run(xc, function(col) xtx_c_on(col, xc))
        est <- xtx_fe_uncenter(r$est, ctr$T, which(colnames(xc) == "(Intercept)"))
        if (!is.null(est)) {
            r$est <- est
            return(r)
        }
    }
    run(x, identity)
}

## The linear predictor of a GLM fit on the rows of design x: X b, plus the
## FE values of the rows' levels (keys: codes, or key columns matched to
## `levels`), plus the offset.
xtx_glm_eta <- function(object, x, off = NULL, keys = NULL, levels = NULL) {
    eta <- if (is.null(object$fixef)) xtx_p_mult(x, xtx_fit_beta(object))
           else xtx_fe_fitted(x, xtx_fe_beta(object, x), object$fixef, keys, levels)
    if (is.null(off)) eta else eta + off
}

## ---------------------------------------------------------------------------
## Methods

## fixest's small sample correction of the iid covariance: (n - 1) / (n - K).
xtx_fe_glm_adj <- function(object) {
    n <- object$nobs
    (n - 1) / (n - object$rank - object$fe_dof)
}

#' @export
vcov.xtx_glm_fe <- function(object, complete = TRUE, ...) {
    V <- if (identical(object$vcov_type, "cluster")) object$vcov_cluster
         else object$cov.unscaled * xtx_glm_dispersion(object)$value * xtx_fe_glm_adj(object)
    if (!complete) {
        keep <- !is.na(object$coefficients)
        V <- V[keep, keep, drop = FALSE]
    }
    V
}

#' @export
logLik.xtx_glm_fe <- function(object, ...) {
    ll <- NextMethod()
    df <- attr(ll, "df") + object$fe_dof
    structure(as.vector(ll) + object$fe_dof, nobs = object$nobs, df = df, class = "logLik")
}

#' @export
summary.xtx_glm_fe <- function(object, ...) {
    est <- object$coefficients
    aliased <- is.na(est)
    keep <- !aliased
    disp <- xtx_glm_dispersion(object)
    V <- vcov(object)
    se <- sqrt(diag(V)[keep])
    tval <- est[keep] / se
    cluster <- identical(object$vcov_type, "cluster")
    if (!disp$est) {
        coefs <- cbind(est[keep], se, tval, 2 * stats::pnorm(-abs(tval)))
        colnames(coefs) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    } else {
        tdf <- if (cluster) object$G_cluster - 1L else object$df.residual
        coefs <- cbind(est[keep], se, tval, 2 * stats::pt(-abs(tval), tdf))
        colnames(coefs) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    }
    ans <- c(object[c("call", "terms", "family", "deviance", "aic", "contrasts", "df.residual",
                      "null.deviance", "df.null", "iter", "na.action", "compressed", "patterns",
                      "rows", "method", "nobs", "fe_sizes", "nsingletons", "nperfect", "vcov_type",
                      "cluster", "G_cluster")],
             list(coefficients = coefs, aliased = aliased, dispersion = disp$value,
                  df = c(object$rank, object$df.residual, length(aliased)), cov.scaled = V))
    class(ans) <- c("summary.xtx_glm_fe", "summary.xtx_glm")
    ans
}

#' @export
print.summary.xtx_glm_fe <- function(x, digits = max(3L, getOption("digits") - 3L),
                                        signif.stars = getOption("show.signif.stars"), ...) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
    xtx_fe_glm_print_rows(x)
    cat("Standard errors: ", if (identical(x$vcov_type, "cluster"))
        sprintf("clustered by %s (%s clusters)", x$cluster, xtx_fmt_count(x$G_cluster)) else "iid",
        "\n", sep = "")
    xtx_fit_print_coefmat(x, digits, signif.stars, ...)
    cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ", format(x$dispersion),
        ")\n\n", apply(cbind(paste(format(c("Null", "Residual"), justify = "right"), "deviance:"),
                             format(unlist(x[c("null.deviance", "deviance")]), digits = max(5L, digits + 1L)),
                             " on", format(unlist(x[c("df.null", "df.residual")])), " degrees of freedom\n"),
                       1L, paste, collapse = " "), sep = "")
    cat("AIC: ", format(x$aic, digits = max(4L, digits + 1L)), "\n\n",
        "Number of Fisher Scoring iterations: ", x$iter, "\n\n", sep = "")
    invisible(x)
}

## Observations, rows removed and fixed effects of a fit or its summary.
xtx_fe_glm_print_rows <- function(x) {
    rm <- c(if (x$nsingletons) sprintf("%s singletons", xtx_fmt_count(x$nsingletons)),
            if (x$nperfect) sprintf("%s perfectly fit by the fixed effects", xtx_fmt_count(x$nperfect)))
    cat("\nObservations: ", xtx_fmt_count(x$nobs),
        if (length(rm)) sprintf(" (removed: %s)", paste(rm, collapse = ", ")), "\n", sep = "")
    if (!is.null(x$method)) cat("Method: ", x$method, "\n", sep = "")
    cat("Fixed effects: ", paste(sprintf("%s (%s)", names(x$fe_sizes), xtx_fmt_count(x$fe_sizes)),
                                 collapse = ", "), "\n", sep = "")
}

#' @export
print.xtx_glm_fe <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
    xtx_fe_glm_print_rows(x)
    xtx_fit_print_coef(x, digits, none = "\nNo coefficients\n\n")
    cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ", x$df.residual, "Residual\n")
    cat("Null Deviance:\t   ", format(signif(x$null.deviance, digits)),
        "\nResidual Deviance:", format(signif(x$deviance, digits)), "\tAIC:",
        format(signif(x$aic, digits)), "\n")
    invisible(x)
}

## fixef(): the FE values, one named vector per FE (normalized as for
## linear models: see fixef.xtx_lm_fe()).
#' @rawNamespace S3method(fixest::fixef, xtx_glm_fe)
#' @rawNamespace S3method(nlme::fixef, xtx_glm_fe)
fixef.xtx_glm_fe <- function(object, ...) object$fixef

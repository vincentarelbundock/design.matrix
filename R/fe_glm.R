## Generalized linear models with absorbed fixed effects:
## design_glm(y ~ x | fe1 + fe2, family, data).
##
## The IRLS iterations of glm.fit() (dmglm_loop(), glm.R), with the
## weighted least squares step of each iteration replaced by the estimator
## of linear models with fixed effects (dm_fe_core(), fe.R): the working
## response z on the design, with the working weights, absorbing the FE.
## The parameters of the loop are the slopes and the FE values together
## (theta = c(b, fe_1, ..., fe_K)), so that its step halving, convergence
## rule (the deviance criterion of glm.fit()) and warnings are those of
## glm(); the linear predictor is X b + sum_k fe_k[code_k] + offset, a
## product with the design and lookups of the FE values (dm_fe_fitted()).
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
## dmglm_fit()): parameters theta = c(b, fe_1, ..., fe_K), b the
## coefficients of the columns of x but the intercept (absorbed).
dm_fe_solver <- function(x, codes, tol, opts) {
    p <- ncol(x)
    xcols <- setdiff(seq_len(p), which(colnames(x) == "(Intercept)"))
    px <- length(xcols)
    G <- vapply(codes, function(g) g$G, 1L)
    offs <- c(0L, cumsum(G))
    warm <- new.env(parent = emptyenv())
    lin <- function(theta, along = list()) {
        Th <- as.matrix(theta)
        beta <- matrix(0, p, ncol(Th))
        beta[xcols, ] <- Th[seq_len(px), ]
        fe <- lapply(seq_along(G), function(k) Th[px + offs[k] + seq_len(G[k]), , drop = FALSE])
        if (!is.matrix(theta)) {
            beta <- beta[, 1L]
            fe <- lapply(fe, function(v) v[, 1L])
        }
        dm_fe_fitted(x, beta, fe, codes, along = along)
    }
    wls <- function(w2, z) {
        ## (the residuals, for clustered standard errors, only for the last
        ## step: see we())
        est <- dm_fe_core(x, list(z = z, w = w2, codes = codes), tol, opts, warm = warm, resid_final = FALSE)
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
        f <- lin(c(b0, unlist(est$fe, use.names = FALSE)), along = list(warm$z, warm$w2))
        warm$w2 * (warm$z - f)
    }
    list(lin = lin, wls = wls, we = we, np = px + sum(G), warm = warm, xcols = xcols, G = G, offs = offs)
}

## Does the link give infinite FE for levels with an all-0 (all-1)
## response? fixest's rule.
dm_fe_bounds <- function(family) {
    list(zero = isTRUE(abs(family$linkinv(-100)) < 1e-15),
         one = isTRUE(abs(family$linkinv(100) - 1) < 1e-15))
}

## The columns of a GLM (dmglm_cols()) from the base of an FE fit: in
## memory from the response of the model frame, remotely from its columns.
dm_fe_glm_cols <- function(base, k) {
    if (!is.null(base$resp)) {
        w <- base$w
        return(dmglm_cols(base$z, base$resp$a, base$resp$b, if (is.null(w)) dmc_const(base$x, 1) else w,
                          base$off, k, base$resp$cbind, weighted = !is.null(w)))
    }
    Y <- base$z
    if (is.matrix(Y) && ncol(Y) == 1L) Y <- Y[, 1L]
    cols <- dmglm_mem_cols(Y, base$w, base$off, k)
    cols$weighted <- !is.null(base$w) || cols$cbind
    cols
}

## Entry point from design_glm(): `call` is its matched call, `env` its
## parent frame.
dm_fe_glm <- function(call, env, formula, family, data, start, control, contrasts, compress, model,
                      vcov, fixef.rm, dots, na_kind) {
    k <- dmglm_kind(family)
    if (isTRUE(compress)) stop("`compress = TRUE` is not supported with fixed effects.", call. = FALSE)
    if (!is.null(start)) stop("`start` is not supported with fixed effects.", call. = FALSE)
    if (na_kind == "fail") stop("With fixed effects, `na.action` must be na.omit or na.exclude.", call. = FALSE)
    pf <- dm_fe_parse(formula)
    clv <- dm_fe_vcov_arg(vcov)
    fixef.rm <- dm_fe_rm_arg(fixef.rm, glm = TRUE)
    ## the design of the rows with complete variables, response, weights
    ## and offset (in memory: a model frame; remotely: a remote design)
    base <- dm_fe_base(pf, data, call, env, contrasts, dots, clv, glm = TRUE)
    gc <- dm_fe_glm_cols(base, k)
    ## the columns of the GLM other than y and w go along with the rows kept
    flags <- gc[c("cbind", "weighted")]
    base$z <- gc$y
    base$w <- gc$w
    base$more <- Filter(Negate(is.null), gc[c("n", "off", "a", "b", "y0")])
    prob <- dm_fe_prepare(base, pf, clv, fixef.rm, dm_fe_bounds(family), wmsg = "negative weights not allowed")
    cols <- c(list(y = prob$cols$z, w = prob$cols$w), prob$more, flags)
    if (is.null(cols$y0)) cols["y0"] <- list(NULL)
    res <- dm_fe_glm_estimate(prob$x, cols, prob$cols$codes, prob$cols$cl, family, control)
    fit <- res$fit
    est <- res$est
    x <- prob$x
    nm <- colnames(x)[est$xcols]
    b <- est$b
    names(b) <- nm
    fe <- est$fe
    names(fe) <- names(pf$fe)
    lv <- prob$levels
    names(lv) <- names(pf$fe)
    for (j in seq_along(fe)) names(fe[[j]]) <- dm_fe_labels(lv[[j]])
    G <- est$G
    names(G) <- names(pf$fe)
    cu <- est$bread
    dimnames(cu) <- list(nm, nm)
    fit$coefficients <- b
    fit$rank <- est$rank
    fit$cov.unscaled <- cu
    fit$pivot <- fit$R <- NULL
    fit$fixef <- fe
    fit$fe_levels <- lv
    fit$fe_vars <- pf$fe
    fit$fe_sizes <- G
    fit$fe_dof <- est$fe_dof
    fit$fe_dof_exact <- est$fe_dof_exact
    fit$nsingletons <- prob$nsingle
    fit$nperfect <- prob$nperfect
    fit$fixef.rm <- fixef.rm
    fit$fe_info <- est$info
    fit$fe_iter <- res$fe_iter
    fit$vcov_type <- "iid"
    if (!is.null(clv)) {
        Vc <- est$vcov_cluster
        dimnames(Vc) <- list(nm, nm)
        fit$vcov_cluster <- Vc
        fit$vcov_type <- "cluster"
        fit$cluster <- clv
        fit$G_cluster <- est$G_cluster
    }
    ## the rows of the fit: the design (not centered) and its columns
    fit$data <- list(design = x, y = cols$y, weights = cols$w, offset = cols$off, codes = prob$cols$codes)
    fit$solver_fit <- NULL
    fit$call <- call
    fit$terms <- base$terms
    fit$formula <- formula
    fit$contrasts <- base$contrasts
    fit$xlevels <- base$xlevels
    fit$spec <- base$spec
    fit$offset_call <- call$offset
    fit$env <- env
    fit$rows <- base$nall
    ## (remotely, na.exclude pads the results with every row not used)
    fit$na.action <- if (isTRUE(base$exclude)) dm_remote_exclude(prob$x_rows)
                     else dm_fe_na_action(base$na_act, prob$rows, base$nall)
    fit$method <- switch(dm_engine(x), sql = "SQL iterations", chunked = "iterations in batches", NULL)
    ## (the data of a remote fit stay in the database: always kept)
    if (!model && is.null(fit$method)) fit$data <- NULL
    class(fit) <- c("design_glm_fe", "design_glm")
    fit
}

## The fit, on the design with numeric variables centered when that is
## exact (dmp_center(); coefficients, covariance and FE mapped back with
## dm_fe_uncenter()), else on the design. Returns list(fit (of
## dmglm_fit()), est (slopes b, FE values fe, bread, cluster covariance,
## degrees of freedom), fe_iter).
dm_fe_glm_estimate <- function(x, cols, codes, cl, family, control) {
    opts <- dm_fe_opts()
    ## aliasing: fixest's criterion, with feglm()'s collin.tol = 1e-9
    tol <- 1e-9
    ctr <- if (isTRUE(getOption("design.matrix.fe_center", TRUE))) dmp_center(x) else NULL
    run <- function(xx, on) {
        cc <- cols
        for (nm in c("y", "w", "n", "off", "a", "b", "y0")) if (!is.null(cc[[nm]])) cc[[nm]] <- on(cc[[nm]])
        cd <- lapply(codes, function(g) dmc_codes(on(g$codes), g$G))
        cl2 <- if (!is.null(cl)) dmc_codes(on(cl$codes), cl$G)
        solver <- dm_fe_solver(xx, cd, tol, opts)
        fit <- dmglm_fit(xx, cc, family, NULL, control, TRUE, FALSE, wmsg = "negative weights not allowed",
                         solver = solver)
        est <- fit$solver_fit$est
        if (!is.null(cl2)) {
            est$we <- solver$we(est)
            est <- dm_fe_cluster(xx, est, cd, cl2)
        }
        th <- fit$coefficients
        px <- length(solver$xcols)
        est$b <- th[seq_len(px)]
        est$fe <- lapply(seq_along(solver$G), function(j) th[px + solver$offs[j] + seq_len(solver$G[j])])
        est$we <- est$A <- NULL
        list(fit = fit, est = est, fe_iter = solver$warm$iter)
    }
    if (!is.null(ctr)) {
        xc <- ctr$x
        r <- run(xc, function(col) dmc_on(col, xc))
        est <- dm_fe_uncenter(r$est, ctr$T, which(colnames(xc) == "(Intercept)"))
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
dmglm_eta <- function(object, x, off = NULL, keys = NULL, levels = NULL) {
    eta <- if (is.null(object$fixef)) dmp_mult(x, dmfit_beta(object))
           else dm_fe_fitted(x, dm_fe_beta(object, x), object$fixef, keys, levels)
    if (is.null(off)) eta else eta + off
}

## ---------------------------------------------------------------------------
## Methods

## fixest's small sample correction of the iid covariance: (n - 1) / (n - K).
dm_fe_glm_adj <- function(object) {
    n <- object$nobs
    (n - 1) / (n - object$rank - object$fe_dof)
}

#' @export
vcov.design_glm_fe <- function(object, complete = TRUE, ...) {
    V <- if (identical(object$vcov_type, "cluster")) object$vcov_cluster
         else object$cov.unscaled * dmglm_dispersion(object)$value * dm_fe_glm_adj(object)
    if (!complete) {
        keep <- !is.na(object$coefficients)
        V <- V[keep, keep, drop = FALSE]
    }
    V
}

#' @export
logLik.design_glm_fe <- function(object, ...) {
    ll <- NextMethod()
    df <- attr(ll, "df") + object$fe_dof
    structure(as.vector(ll) + object$fe_dof, nobs = object$nobs, df = df, class = "logLik")
}

#' @export
summary.design_glm_fe <- function(object, ...) {
    est <- object$coefficients
    aliased <- is.na(est)
    keep <- !aliased
    disp <- dmglm_dispersion(object)
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
    class(ans) <- c("summary.design_glm_fe", "summary.design_glm")
    ans
}

#' @export
print.summary.design_glm_fe <- function(x, digits = max(3L, getOption("digits") - 3L),
                                        signif.stars = getOption("show.signif.stars"), ...) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
    dm_fe_glm_print_rows(x)
    cat("Standard errors: ", if (identical(x$vcov_type, "cluster"))
        sprintf("clustered by %s (%s clusters)", x$cluster, dm_fmt_count(x$G_cluster)) else "iid",
        "\n", sep = "")
    dmfit_print_coefmat(x, digits, signif.stars, ...)
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
dm_fe_glm_print_rows <- function(x) {
    rm <- c(if (x$nsingletons) sprintf("%s singletons", dm_fmt_count(x$nsingletons)),
            if (x$nperfect) sprintf("%s perfectly fit by the fixed effects", dm_fmt_count(x$nperfect)))
    cat("\nObservations: ", dm_fmt_count(x$nobs),
        if (length(rm)) sprintf(" (removed: %s)", paste(rm, collapse = ", ")), "\n", sep = "")
    if (!is.null(x$method)) cat("Method: ", x$method, "\n", sep = "")
    cat("Fixed effects: ", paste(sprintf("%s (%s)", names(x$fe_sizes), dm_fmt_count(x$fe_sizes)),
                                 collapse = ", "), "\n", sep = "")
}

#' @export
print.design_glm_fe <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
    dm_fe_glm_print_rows(x)
    dmfit_print_coef(x, digits, none = "\nNo coefficients\n\n")
    cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ", x$df.residual, "Residual\n")
    cat("Null Deviance:\t   ", format(signif(x$null.deviance, digits)),
        "\nResidual Deviance:", format(signif(x$deviance, digits)), "\tAIC:",
        format(signif(x$aic, digits)), "\n")
    invisible(x)
}

## fixef(): the FE values, one named vector per FE (normalized as for
## linear models: see fixef.design_lm_fe()).
#' @rawNamespace S3method(fixest::fixef, design_glm_fe)
#' @rawNamespace S3method(nlme::fixef, design_glm_fe)
fixef.design_glm_fe <- function(object, ...) object$fixef

## Linear models with absorbed fixed effects: design_lm(y ~ x | fe1 + fe2).
##
## The fit side: formula and arguments, the fields of the fit, methods and
## predictions. The estimator (Frisch-Waugh-Lovell on aggregates, the
## reduced FE system and its solvers, the absorption operator, the removal
## of singletons) is in absorb.R, reached through its exports
## (design_absorb_rows(), design_absorb_wls(), design_center(),
## design_uncenter(), design_absorb_fitted(), absorb_control()).

## ---------------------------------------------------------------------------
## Formula and arguments

## Does the formula have a fixed-effect part (y ~ x | fe)?
xtx_fe_has_bar <- function(formula) {
    if (!inherits(formula, "formula")) return(FALSE)
    rhs <- formula[[length(formula)]]
    is.call(rhs) && identical(rhs[[1L]], as.name("|"))
}

## Split `y ~ x | a + b^c` into the main formula and the FE terms (each a
## character vector of variable names, several for an interaction a^b).
xtx_fe_parse <- function(formula) {
    rhs <- formula[[length(formula)]]
    if (!xtx_fe_has_bar(formula))
        stop("The formula must have a fixed-effect part: y ~ x | fe.", call. = FALSE)
    ## (which variables `.` stands for is not defined with FE, in memory or
    ## remotely)
    if ("." %in% all.vars(formula))
        stop("`.` is not supported in a formula with fixed effects: name the variables.", call. = FALSE)
    main <- formula
    main[[length(formula)]] <- rhs[[2L]]
    split_plus <- function(e) {
        if (is.call(e) && identical(e[[1L]], as.name("+"))) {
            c(split_plus(e[[2L]]), split_plus(e[[3L]]))
        } else {
            list(e)
        }
    }
    split_hat <- function(e) {
        if (is.name(e)) return(as.character(e))
        if (is.call(e) && identical(e[[1L]], as.name("^")))
            return(c(split_hat(e[[2L]]), split_hat(e[[3L]])))
        stop("Fixed effects must be variable names or interactions a^b.", call. = FALSE)
    }
    fe <- lapply(split_plus(rhs[[3L]]), split_hat)
    names(fe) <- vapply(fe, paste, "", collapse = "^")
    fe <- fe[!duplicated(names(fe))]
    list(main = main, fe = fe)
}

## The cluster variable of `vcov`: NULL for "iid", else a column name.
xtx_fe_vcov_arg <- function(vcov) {
    msg <- "`vcov` must be \"iid\" or a one-sided formula naming a cluster variable, such as ~firm."
    if (is.character(vcov) && length(vcov) == 1L) {
        if (identical(vcov, "iid")) return(NULL)
        stop(msg, call. = FALSE)
    }
    if (inherits(vcov, "formula") && length(vcov) == 2L) {
        if (is.name(vcov[[2L]])) return(as.character(vcov[[2L]]))
        stop("Only one-way clustering is supported: `vcov = ~cluster`, one variable.", call. = FALSE)
    }
    stop(msg, call. = FALSE)
}

## fixef.rm: for linear models "singleton" (or "perfect_fit", the same
## there) or "none"; for GLMs also "infinite_coef" and "perfect_fit"
## (default: singletons and levels with an infinite FE, as fixest).
xtx_fe_rm_arg <- function(fixef.rm, glm = FALSE) {
    ok <- if (glm) c("perfect_fit", "singleton", "infinite_coef", "none") else c("singleton", "perfect_fit", "none")
    if (!checkmate::test_choice(fixef.rm, ok))
        stop(sprintf("`fixef.rm` must be %s.", paste(sprintf("\"%s\"", ok), collapse = ", ")), call. = FALSE)
    if (!glm && fixef.rm == "perfect_fit") fixef.rm <- "singleton"
    fixef.rm
}

## ---------------------------------------------------------------------------
## Labels of fixed effects (their codes: design_absorb_rows())

## Labels of FE levels, as fixest names them (a_b for interactions).
xtx_fe_labels <- function(lv) {
    if (ncol(lv) == 1L) return(as.character(lv[[1L]]))
    do.call(paste, c(lapply(lv, as.character), sep = "_"))
}

## ---------------------------------------------------------------------------
## design_lm() with fixed effects
##
## Only the first step differs between data in memory and remote data, as
## it does without FE: building the design of the rows with complete model
## variables, with the response and weights as columns (a model frame in
## memory; a remote design and xtx_p_resolve() remotely). The rows used, the
## FE codes, the removal of singletons, the estimator, the methods and the
## predictions are one code path on the primitives.

#' @section Fixed effects:
#' A formula `y ~ x | fe1 + fe2` absorbs fixed effects, as
#' [fixest::feols()]: each FE is a variable (a factor, character, integer or
#' numeric column) or an interaction `a^b`. The fit never builds the FE
#' dummies or demeans the rows: it solves the normal equations from sums
#' over rows (the design's products, sums within FE levels and crosstabs of
#' the FE), eliminating the largest FE exactly and solving for the others
#' densely (up to 1000 levels) or by conjugate gradients: on the reduced
#' system formed, with an incomplete Cholesky preconditioner, when it is
#' sparse (low mobility, such as workers who rarely change firms), and
#' without forming it, with a Jacobi preconditioner, when it would be near
#' dense (high mobility). The same code fits data in memory and remote
#' data (DuckDB tables or queries, Parquet), with identical results up to
#' rounding.
#'
#' * Rows: rows with a missing value in the model, FE or cluster variables
#'   are dropped, as are rows with zero weight and, with `fixef.rm =
#'   "singleton"` (default), singletons (rows alone in a level of some FE),
#'   iteratively, as `feols()` does. `fit$nsingletons` counts them.
#' * Coefficients: the intercept is absorbed; a column collinear with the
#'   FE or earlier columns is `NA` (fixest's rule: its residual variance is
#'   below `tol` times its variance; `tol` defaults to `1e-10` with FE, as
#'   `collin.tol`). The residual variance comes from sums over rows; for a
#'   column within rounding of that threshold (about `2e-12` times its sum
#'   of squares) it is computed again from the data, one pass over the
#'   rows, so that the decision does not depend on rounding. One step of
#'   iterative refinement from the residuals is
#'   taken, and numeric variables far from zero are centered when that is
#'   an exact transform of the design, so that coefficients and standard
#'   errors stay accurate.
#' * Standard errors: `vcov = "iid"` (default) or `vcov = ~cl` (clustered by
#'   one variable), with fixest's default small sample corrections
#'   (`ssc(K.adj = TRUE, K.fixef = "nonnested", G.adj = TRUE)`); clustered
#'   tests in `summary()` use `t(G - 1)`. The cluster scores are corrected
#'   for FE not nested in the clusters. `fit$fe_dof` is the number of FE
#'   parameters counted (`sum(levels) - (number of FE - 1)`); `fit$fe_dof_exact`
#'   the exact rank of the FE (fixest's `K.exact = TRUE`), exact for one or
#'   two FE and for three or more when the reduced system was factored
#'   (otherwise `NA`).
#' * Methods: those of a `design_lm` fit (class `c("xtx_lm_fe",
#'   "xtx_lm")`); `summary()` reports R-squared and within R-squared.
#'   `fixef()` (the generic of fixest or nlme, when loaded) returns the FE
#'   values, also in `fit$fixef`: one named vector per FE, normalized with
#'   one level per connected set of levels at zero in every FE but the
#'   largest (differences within a connected set are those of fixest).
#'   `predict(fit, newdata)` is `X b` plus the FE values of the rows' levels
#'   (`NA` for levels not in the fit). With remote data, `fitted()` and
#'   `predict()` return tables and `residuals()` a remote design, as for
#'   remote fits without FE. These methods read the source rows when called;
#'   keep the training source available and unchanged for in-sample results.
#' * Offsets (the `offset` argument and `offset()` terms) are subtracted
#'   from the response, as in [fixest::feols()]; `fitted()` and `predict()`
#'   include them.
#' * Not supported with FE: `subset`, `compress = TRUE`, `na.action`
#'   other than `na.omit` and `na.exclude` (in memory and remotely),
#'   `.` in the formula, `se.fit` and intervals in `predict()`, clustering
#'   by more than one variable.
#'   For remote data, the response and `weights` are column names,
#'   `DBI::SQL()` or R expressions, as without FE (expressions with no SQL
#'   translation are evaluated in R, in batches).
#' * Options: `xtx.fe_solver` (`"auto"`, `"dense"`, `"pcg"`, or
#'   `"chol"` for a sparse Cholesky factor), `xtx.fe_precond`
#'   (`"auto"`, default: incomplete Cholesky when the reduced system is
#'   formed, Jacobi otherwise; `"jacobi"`; or `"ichol"`, which always forms
#'   it), `xtx.fe_form_max` (when to form the reduced system:
#'   default 25, `0` never, `Inf` always), `xtx.fe_refine` (refinement steps:
#'   `"auto"`, default, one step, after which the residuals are computed
#'   again unless the step's correction is negligible; or a count of steps,
#'   each followed by the residuals) and `xtx.fe_center` (default `TRUE`).
#' @param vcov With fixed effects: `"iid"` or a one-sided formula naming a
#'   cluster variable (`~firm`). Without fixed effects only `"iid"`.
#' @param fixef.rm With fixed effects: `"singleton"` (drop singletons
#'   iteratively, as [fixest::feols()]; `"perfect_fit"` is the same for a
#'   linear model) or `"none"`.
#' @examples
#' ## fixed effects (examples below: without)
#' set.seed(2)
#' panel <- data.frame(id = rep(1:100, each = 10), t = rep(1:10, 100), x = rnorm(1000))
#' panel$y <- panel$x + panel$id / 50 + rnorm(1000)
#' fe <- design_lm(y ~ x | id + t, panel, vcov = ~id)
#' summary(fe)
#' head(fe$fixef$id)
#'
#' @rdname design_lm
#' @name design_lm
NULL

## Entry point from design_lm(): `call` is its matched call, `env` its
## parent frame.
xtx_fe_lm <- function(call, env, formula, contrasts, compress, tol, tol_missing, model,
                     vcov, fixef.rm, dots) {
    ## (arguments checked by design_lm(): checks.R)
    ## aliasing: fixest's criterion (residual variance below tol times the
    ## variance before the FE), with fixest's default collin.tol = 1e-10
    if (tol_missing) tol <- 1e-10
    pf <- xtx_fe_parse(formula)
    clv <- xtx_fe_vcov_arg(vcov)
    fixef.rm <- xtx_fe_rm_arg(fixef.rm)
    data <- if (is.null(call$data)) NULL else eval(call$data, env)
    ## the design of the rows with complete model variables, with the
    ## response and weights as columns (remote_fit.R)
    base <- design_problem(call, env, pf$main, data, contrasts, dots, fe = pf$fe, cluster = clv)
    if (is.matrix(base$y)) stop("Multivariate responses are not supported.", call. = FALSE)
    if (!is.numeric(base$y) && !inherits(base$y, "xtx_column")) stop("The response must be numeric.", call. = FALSE)
    ## an offset: the fit is that of y - offset
    if (!is.null(base$off)) {
        base$y <- base$y - base$off
        base$more <- list(off = base$off)
    }
    prob <- xtx_fe_rows(base, fixef.rm)
    est <- xtx_fe_estimate(prob$x, prob$cols, tol, !is.null(clv))
    fit <- xtx_fe_finish(est, colnames(prob$x), pf, prob$levels, clv, prob$nsingle, base$nall)
    w <- prob$cols$w
    fit$weighted <- !is.null(w)
    fit$sumlogw <- if (is.null(w)) 0 else sum(log(w))
    fit$x <- prob$x
    fit$cols <- prob$cols
    fit$cols$off <- prob$more$off
    fit$offset_call <- call$offset
    fit$env <- base$env
    fit$terms <- base$terms
    fit$contrasts <- base$contrasts
    fit$xlevels <- base$xlevels
    fit$spec <- xtx_fit_spec(base$x)
    ## (remotely, na.exclude pads the results with every row not used:
    ## na_action.R)
    fit$na.action <- xtx_na_result(base, prob$rows, prob$x_rows)
    fit$method <- if (inherits(prob$x, "xtx_remote")) "SQL aggregates"
    ## (the data of a remote fit stay in the database: always kept)
    if (!is.null(fit$method)) {
        ## nothing to drop
    } else if (model) {
        fit$model <- prob$frame
    } else {
        fit$x <- fit$cols <- NULL
    }
    fit$call <- call
    fit
}

## The fields of a fit with FE shared by design_lm() and design_glm(), from
## a fit of design_absorb_wls() (`est`), the column names of the design (cn), the
## parsed formula (pf) and the FE levels (lv): the coefficients named by
## the columns but the intercept, their unscaled covariance, the FE values
## named by level, the FE sizes and degrees of freedom (`fields`); the
## clustered covariance with dimnames (`vcov_cluster`, NULL without
## clusters) and the fields that describe it (`inference`). Each caller
## finishes its own inference fields: a clustered LM stores the clustered
## covariance as fit$vcov, a GLM as fit$vcov_cluster. `named`: dimnames
## for a covariance of the coefficients.
xtx_fe_fields <- function(est, cn, pf, lv, clv) {
    nm <- cn[est$xcols]
    named <- function(V) {
        dimnames(V) <- list(nm, nm)
        V
    }
    b <- est$b
    names(b) <- nm
    fe <- est$fe
    names(fe) <- names(pf$fe)
    names(lv) <- names(pf$fe)
    for (k in seq_along(fe)) names(fe[[k]]) <- xtx_fe_labels(lv[[k]])
    G <- est$G
    names(G) <- names(pf$fe)
    fields <- list(coefficients = b, cov.unscaled = named(est$bread), rank = est$rank, fixef = fe,
                   fe_levels = lv, fe_vars = pf$fe, fe_sizes = G, fe_dof = est$fe_dof,
                   fe_dof_exact = est$fe_dof_exact, fe_info = est$info)
    inference <- if (is.null(clv)) list(vcov_type = "iid")
                 else list(vcov_type = "cluster", cluster = clv, G_cluster = est$G_cluster)
    list(fields = fields, named = named, inference = inference,
         vcov_cluster = if (!is.null(clv)) named(est$vcov_cluster))
}

## The design_lm object of an FE fit.
xtx_fe_finish <- function(est, cn, pf, lv, clv, nsingle, nall) {
    sh <- xtx_fe_fields(est, cn, pf, lv, clv)
    sigma2 <- est$rss / est$df.residual
    fit <- c(sh$fields, list(vcov = sh$named(est$bread * sigma2), sigma = sqrt(sigma2), rss = est$rss,
                             df.residual = est$df.residual, nobs = est$nobs, nsingletons = nsingle,
                             tss = est$tss, tss_within = est$tss_within, compressed = FALSE,
                             patterns = est$nobs, rows = nall, intercept = TRUE))
    fit[names(sh$inference)] <- sh$inference
    ## clustered: the covariance of the fit is the clustered one
    if (!is.null(clv)) fit$vcov <- sh$vcov_cluster
    class(fit) <- c("xtx_lm_fe", "xtx_lm")
    fit
}

## ---------------------------------------------------------------------------
## Methods

#' @export
vcov.xtx_lm_fe <- function(object, ...) object$vcov

#' @export
logLik.xtx_lm_fe <- function(object, ...) {
    N <- object$nobs
    val <- 0.5 * (object$sumlogw - N * (log(2 * pi) + 1 - log(N) + log(object$rss)))
    structure(val, nall = N, nobs = N, df = object$rank + object$fe_dof + 1L, class = "logLik")
}

#' @export
print.xtx_lm_fe <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    if (!is.null(x$call)) cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
    cat("\nFixed effects: ", paste(sprintf("%s (%s)", names(x$fe_sizes),
                                            xtx_fit_fmt_count(x$fe_sizes)), collapse = ", "), "\n", sep = "")
    xtx_fit_print_rows(x, prefix = "")
    xtx_fit_print_coef(x, digits)
    cat("\n")
    invisible(x)
}

#' @export
summary.xtx_lm_fe <- function(object, ...) {
    est <- object$coefficients
    aliased <- is.na(est)
    keep <- !aliased
    se <- sqrt(diag(object$vcov)[keep])
    tval <- est[keep] / se
    ## clustered: t(G - 1), as fixest
    tdf <- if (identical(object$vcov_type, "cluster")) object$G_cluster - 1L else object$df.residual
    coefs <- cbind(Estimate = est[keep], "Std. Error" = se, "t value" = tval,
                   "Pr(>|t|)" = 2 * stats::pt(abs(tval), tdf, lower.tail = FALSE))
    ans <- list(call = object$call, terms = object$terms, coefficients = coefs,
                aliased = aliased, sigma = object$sigma,
                df = c(object$rank, object$df.residual, length(est)),
                r.squared = 1 - object$rss / object$tss,
                within.r.squared = if (object$tss_within > 0) 1 - object$rss / object$tss_within else NA_real_,
                vcov_type = object$vcov_type, cluster = object$cluster, G_cluster = object$G_cluster,
                fe_sizes = object$fe_sizes, nobs = object$nobs, nsingletons = object$nsingletons,
                compressed = FALSE, method = object$method,
                rows = object$rows)
    n <- object$nobs
    ans$adj.r.squared <- 1 - (1 - ans$r.squared) * (n - 1) / object$df.residual
    class(ans) <- c("summary.xtx_lm_fe", "summary.xtx_lm")
    ans
}

#' @export
print.summary.xtx_lm_fe <- function(x, digits = max(3L, getOption("digits") - 3L),
                                       signif.stars = getOption("show.signif.stars"), ...) {
    if (!is.null(x$call)) cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
    cat("\nObservations: ", xtx_fit_fmt_count(x$nobs),
        if (x$nsingletons) sprintf(" (%s singletons removed)", xtx_fit_fmt_count(x$nsingletons)), "\n", sep = "")
    cat("Fixed effects: ", paste(sprintf("%s (%s)", names(x$fe_sizes), xtx_fit_fmt_count(x$fe_sizes)),
                                 collapse = ", "), "\n", sep = "")
    cat("Standard errors: ", if (identical(x$vcov_type, "cluster"))
        sprintf("clustered by %s (%s clusters)", x$cluster, xtx_fit_fmt_count(x$G_cluster)) else "iid",
        "\n", sep = "")
    xtx_fit_print_coefmat(x, digits, signif.stars, ...)
    cat("\nResidual standard error:", format(signif(x$sigma, digits)),
        "on", x$df[2L], "degrees of freedom\n")
    cat("Multiple R-squared: ", formatC(x$r.squared, digits = digits),
        ",\tWithin R-squared: ", formatC(x$within.r.squared, digits = digits), "\n\n", sep = "")
    invisible(x)
}

## The fitted values and residuals of the rows used, as columns.
xtx_fe_columns <- function(object) {
    x <- object$x
    if (is.null(x))
        stop("This needs the fit's data: refit with `model = TRUE`, or use `predict(fit, newdata)`.",
             call. = FALSE)
    cols <- object$cols
    f <- design_absorb_fitted(x, xtx_fe_beta(object, x), object$fixef, cols$codes)
    ## (z is the response minus the offset)
    list(fitted = if (is.null(cols$off)) f else f + cols$off, resid = cols$z - f)
}

## Coefficients for the columns of design x (0 for the intercept, absorbed,
## and for aliased columns).
xtx_fe_beta <- function(object, x) {
    b <- object$coefficients
    b[is.na(b)] <- 0
    beta <- numeric(ncol(x))
    beta[match(names(b), colnames(x))] <- b
    beta
}

#' @export
residuals.xtx_lm_fe <- function(object, ...) {
    ## remote: the remote design with the columns .fitted and .resid
    r <- xtx_fe_columns(object)
    design_output(object$x, list(.fitted = r$fitted, .resid = r$resid), design = TRUE,
               memory = function(d) stats::naresid(object$na.action, d$.resid))
}

#' @export
fitted.xtx_lm_fe <- function(object, ...) {
    ## remote: a table with .fitted and .resid
    r <- xtx_fe_columns(object)
    design_output(object$x, list(.fitted = r$fitted, .resid = r$resid),
               memory = function(d) stats::napredict(object$na.action, d$.fitted))
}

#' @export
predict.xtx_lm_fe <- function(object, newdata, se.fit = FALSE,
                                 interval = c("none", "confidence", "prediction"),
                                 level = 0.95, na.action = stats::na.pass, ...) {
    interval <- match.arg(interval)
    if (isTRUE(se.fit) || interval != "none")
        stop("With fixed effects, predict() gives the fit only (no `se.fit` or intervals).", call. = FALSE)
    if (missing(newdata) || is.null(newdata)) return(fitted(object))
    if (!identical(na.action, stats::na.pass))
        stop("With fixed effects, predict() keeps every row of `newdata` (na.action = na.pass).", call. = FALSE)
    ## the design of the new rows (as predict() without FE), their offsets
    nd <- xtx_fit_newdata(object, newdata, stats::na.pass, list(...))
    Dn <- nd$D
    keys <- xtx_fe_new_keys(object, Dn, newdata)
    f <- design_absorb_fitted(Dn, xtx_fe_beta(object, Dn), object$fixef, keys)
    if (!is.null(nd$off)) f <- f + nd$off
    ## a table remotely, a named vector in memory
    design_output(Dn, list(.fitted = f), memory = function(d) stats::setNames(d$.fitted, dimnames(Dn)[[1L]]))
}

## The FE of a fit over the rows of a new design Dn: groupings coded
## against the levels of the fit (NA for levels not in the fit).
xtx_fe_new_keys <- function(object, Dn, newdata) {
    lapply(seq_along(object$fe_vars), function(j) {
        vars <- object$fe_vars[[j]]
        ## (in memory the values of newdata; remotely its columns)
        key <- if (is.data.frame(newdata)) tryCatch(newdata[vars], error = function(e) NULL) else vars
        k <- if (!is.null(key)) tryCatch(design_groups(Dn, key, levels = object$fe_levels[[j]]),
                                         error = function(e) NULL)
        if (is.null(k))
            stop(sprintf("`newdata` has no fixed-effect variable among `%s`.", paste(vars, collapse = "`, `")),
                 call. = FALSE)
        k
    })
}

## fixef(): the FE values, a list with one named vector per FE (fixest's
## normalization differs: compare sums over FE, or differences within a
## connected set of levels). The method is registered for fixest's and
## nlme's fixef() generics when those packages are loaded, without a
## dependency on them; the values are also `fit$fixef`.
#' @rawNamespace S3method(fixest::fixef, xtx_lm_fe)
#' @rawNamespace S3method(nlme::fixef, xtx_lm_fe)
fixef.xtx_lm_fe <- function(object, ...) object$fixef


## Clustered covariance of a fit of design_absorb_wls() (`est`, which keeps
## its absorption operator), clusters `cl` (a grouping of
## design_absorb_rows()), weighted residuals `we`: the raw cluster scores
## of the operator with fixest's small sample correction; adds
## vcov_cluster and G_cluster to `est`.
xtx_fe_cluster <- function(est, cl, we) {
    px <- length(est$xcols)
    kp <- est$kp
    rank <- est$rank
    G <- est$G
    n <- est$nobs
    Gc <- cl$G
    sc <- est$absorb$scores(we, cl)
    U <- sc$U[, kp, drop = FALSE]
    nn <- !sc$nested
    K_cl <- rank + if (any(nn)) sum(G[nn]) - (sum(nn) - 1L) else 1L
    adj <- Gc / (Gc - 1) * (n - 1) / (n - K_cl)
    Vc <- matrix(NA_real_, px, px)
    if (rank) {
        bk <- est$bread[kp, kp, drop = FALSE]
        Vc[kp, kp] <- adj * bk %*% crossprod(U) %*% bk
    }
    est$vcov_cluster <- Vc
    est$G_cluster <- Gc
    est
}

## The fit: on the design with numeric variables centered when that is
## exact (design_center()), else on the design.
xtx_fe_estimate <- function(x, cols, tol, cluster) {
    opts <- absorb_control()
    run <- function(x, cols) {
        est <- design_absorb_wls(x, cols$z, cols$w, cols$codes, tol, opts)
        if (cluster) est <- xtx_fe_cluster(est, cols$cl, est$we)
        est
    }
    ## (columns carried over to the centered design)
    ctr <- if (isTRUE(getOption("xtx.fe_center", TRUE))) design_center(x, cols) else NULL
    if (!is.null(ctr)) {
        est <- design_uncenter(run(ctr$x, ctr$cols), ctr)
        if (!is.null(est)) return(est)
    }
    run(x, cols)
}

## The rows of an FE fit (design_absorb_rows()): `base` the input stage,
## with the response, weights and other columns of the fit; fixef.rm as
## parsed by xtx_fe_rm_arg(); `bounds` (GLMs): whether levels with an
## all-0 or all-1 response have an infinite FE (xtx_fe_bounds()).
xtx_fe_rows <- function(base, fixef.rm, bounds = list(), wmsg = "missing or negative weights not allowed") {
    w <- base$w
    ## (is.finite() is never NA: the sum counts the negative weights)
    if (!is.null(w) && sum(is.finite(w) & w < 0) > 0) stop(wmsg, call. = FALSE)
    design_absorb_rows(base, singletons = fixef.rm %in% c("singleton", "perfect_fit"),
                       zero = isTRUE(bounds$zero) && fixef.rm %in% c("infinite_coef", "perfect_fit"),
                       one = isTRUE(bounds$one) && fixef.rm %in% c("infinite_coef", "perfect_fit"))
}


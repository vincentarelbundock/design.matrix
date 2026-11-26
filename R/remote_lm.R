## design_lm() on remote data (see remote.R, remote_compress.R, lm.R).
##
## Two ways to fit:
##
##  * compressed: when the rows have few design patterns, one GROUP BY gives
##    each pattern's sum of weights, weighted mean response and weighted sum
##    of squares around that mean (dm_remote_patterns()); the fit is that of
##    design_lm() on the patterns (dmfit_ls(), with lm()'s pivoting), whose
##    residual sum of squares plus the within-pattern sum of squares is that
##    of the fit on every row;
##  * from moments: X'WX, X'Wy and y'Wy from the SQL or chunked engine, and a
##    pivoted Cholesky factor of X'WX with lm()'s column selection
##    (dm_chol_seq()).
##    Normal equations lose accuracy with the square of the condition number,
##    and a numeric variable far from zero (x + 1e7, say) makes X'X nearly
##    singular. So the moments are those of the design with numeric
##    variables shifted by their mean in a sample of rows, when that design
##    is an exact transform X T of the original one with T unit upper
##    triangular (each column changes by multiples of earlier columns, such
##    as the intercept for x, or f for x:f), which keeps lm()'s aliasing. The
##    coefficients and factor are then mapped back: beta = T beta~,
##    R = R~ T^-1.

#' @name design_remote
#' @section Linear models:
#' [design_lm()] works on remote data as on a data frame:
#' `design_lm(y ~ x + f, con, table = "flights", weights = w, subset = year > 2000)`,
#' or on a remote design: `design_lm(D, y = "col", weights = "w")`.
#' `weights`, `offset`, `subset` and `y` are column names (strings or bare
#' names), `DBI::SQL()` expressions, or R expressions with a SQL
#' translation (arithmetic as for the SQL engine; for `subset`, also
#' comparisons, `&`, `|`, `!`, `%in%` with constants and `is.na()`).
#' `offset()` terms in the formula are supported as well. `na.action` can
#' be [na.omit()] (default), [na.exclude()], [na.fail()] or [na.pass()],
#' with the results or errors of [lm()].
#'
#' With few design patterns, the fit is compressed (see *Compression*);
#' otherwise it is computed from `X'WX`, `X'Wy` and `y'Wy` with [lm()]'s
#' column selection, numeric variables centered when that is an exact
#' transform of the design. Coefficients, covariance, residual standard
#' error, degrees of freedom, `logLik()`, `AIC()`, `BIC()` and `summary()`
#' are those of [lm()] on the same rows, up to rounding.
#'
#' Results that have one value per row are tables in the database (names of
#' temporary views, of class `dm_remote_table`, with all the columns of the
#' remote data): `fitted(fit)` (columns `.fitted` and `.resid`; with
#' `na.exclude`, the rows with missing values too, with `NULL`), and
#' `predict(fit)` or `predict(fit, newdata = con, table = ...)` (column
#' `.fitted`; with `se.fit = TRUE`, `.se.fit`, and the attributes `df` and
#' `residual.scale`; with intervals, `.lwr` and `.upr`; for a rank-deficient
#' fit, `.nonestim`, the non-estimable cases flagged by `rankdeficient` as in
#' [stats::predict.lm()], with its warnings). Standard errors, intervals and
#' non-estimability are computed in batches of rows; the fit alone in SQL.
#' `residuals(fit)` returns the remote design of the fit with two more
#' columns, `.fitted` and `.resid`, which can scale its rows: for
#' heteroskedasticity-consistent and clustered standard errors,
#' ```
#' R <- residuals(fit)
#' meat_hc1 <- crossprod(R, DBI::SQL('".resid" * ".resid"') * R)
#' scores <- colSums(".resid" * R, group = "cluster")
#' meat_cl <- crossprod(scores)
#' ```
#' (with `bread <- vcov(fit) / sigma(fit)^2`, the HC1 covariance is
#' `n / (n - k) * bread %*% meat_hc1 %*% bread`). The response, offsets,
#' weights and subset must have a SQL translation; `model.frame()` is not
#' available for a remote fit (`model` is ignored).
NULL

## na.action of a remote fit: "omit", "exclude", "fail" or "pass".
dm_remote_lm_na <- function(na.action) {
    if (is.null(na.action)) return("omit")
    if (is.character(na.action)) na.action <- get(na.action, mode = "function")
    for (nm in c("omit", "exclude", "fail", "pass"))
        if (identical(na.action, get(paste0("na.", nm), envir = asNamespace("stats")))) return(nm)
    stop("For remote data, `na.action` must be na.omit, na.exclude, na.fail or na.pass.", call. = FALSE)
}

## SQL of the offsets of a remote design: offset() terms of its plan, and
## the `offset` argument (the "(offset)" variable). NULL without offsets.
dm_remote_offset_sql <- function(D) {
    tt <- D$plan$terms
    nms <- dm_terms_vars(tt)$names[attr(tt, "offset")]
    nms <- c(nms, intersect("(offset)", names(D$vars)))
    if (!length(nms)) return(NULL)
    sqls <- vapply(D$vars[nms], function(v) {
        if (v$kind != "numeric")
            stop(sprintf("The offset `%s` has no SQL translation.", v$name), call. = FALSE)
        v$sql
    }, "")
    paste(sprintf("(%s)", sqls), collapse = " + ")
}

## design_lm(formula, <remote data>, ...). `weights`, `offset` and `subset`
## are unevaluated arguments (see dm_resolve_arg()), resolved in `env`.
dm_remote_lm <- function(formula, source, weights, offset, subset, na.action, contrasts,
                         compress, tol, call, dots, env) {
    na <- dm_remote_lm_na(na.action)
    D <- dm_remote(formula, source, contrasts, NULL, dots, weights = weights, drop_unused = TRUE,
                   subset = subset, offset = offset, env = env)
    if (is.null(D$response)) stop("The formula has no response.", call. = FALSE)
    ## rows dropped for missing values: na.fail and na.pass refuse them
    if (na %in% c("fail", "pass")) {
        total <- as.numeric(dm_query(D, sprintf("SELECT COUNT(*) AS n FROM (%s) AS s", D$src$select))$n)
        if (total > D$n) {
            if (na == "fail") stop("missing values in object", call. = FALSE)
            ## na.pass: the error of the fit on rows with missing values
            bad <- function(v) {
                if (is.null(v) || v$kind != "numeric") return(FALSE)
                as.numeric(dm_query(D, sprintf("SELECT COUNT(*) AS n FROM (%s) AS s WHERE NOT COALESCE(%s, FALSE)",
                                               D$src$select, dm_remote_var_ok(D, v)))$n) > 0
            }
            offs <- dm_terms_vars(D$plan$terms)$names[attr(D$plan$terms, "offset")]
            if (bad(D$vars[[D$response]]) || any(vapply(D$vars[c(offs, intersect("(offset)", names(D$vars)))], bad, NA)))
                stop("NA/NaN/Inf in 'y'", call. = FALSE)
            if (bad(D$vars[["(weights)"]])) stop("missing or negative weights not allowed", call. = FALSE)
            stop("NA/NaN/Inf in 'x'", call. = FALSE)
        }
    }
    v <- D$vars[[D$response]]
    if (v$kind != "numeric")
        stop("The response of a remote fit must be a numeric variable with a SQL translation.", call. = FALSE)
    off <- dm_remote_offset_sql(D)
    z <- if (is.null(off)) v$sql else sprintf("(%s) - (%s)", v$sql, off)
    w <- D$vars[["(weights)"]]$sql
    fit <- dm_remote_lm_fit(D, z, w, compress, tol)
    fit$remote_na <- na
    fit$remote_offset <- off
    fit$remote_offset_arg <- offset
    fit$remote_env <- env
    fit$call <- call
    fit$terms <- D$terms_full
    fit$contrasts <- D$plan$contrasts
    fit$xlevels <- D$plan$xlevels
    fit
}

## design_lm(<remote design>, y = , weights = , offset = ): unevaluated
## arguments resolved in `env` against the columns of the design's data.
dm_remote_lm_design <- function(D, y, weights, offset, compress, tol, env) {
    types <- D$src$types
    con <- dm_con(D)
    if (is.null(y)) {
        if (is.null(D$response))
            stop("With a remote design, give `y`, the name of one column (or use a design with a response).",
                 call. = FALSE)
        v <- D$vars[[D$response]]
        if (v$kind != "numeric")
            stop("The response of a remote fit must be a numeric variable with a SQL translation.", call. = FALSE)
        ysql <- v$sql
    } else {
        ysql <- dm_resolve_arg(y, env, types, con, "y")$sql
    }
    o <- dm_resolve_arg(offset, env, types, con, "offset")
    w <- dm_resolve_arg(weights, env, types, con, "weights")
    z <- if (is.null(o)) ysql else sprintf("(%s) - (%s)", ysql, o$sql)
    fit <- dm_remote_lm_fit(D, z, w$sql, compress, tol)
    fit$remote_na <- "omit"
    fit$remote_offset <- o$sql
    fit
}

## The fit of a remote design on z (the response minus offsets, SQL over
## `s`) with weights w (SQL over `s`, or NULL).
dm_remote_lm_fit <- function(D, z, w, compress = "auto", tol = 1e-07) {
    dmfit_check_compress(compress)
    intercept <- attr(D$plan$terms, "intercept") > 0L

    ## lm()'s checks, in SQL
    items <- c("COUNT(*) AS n",
               sprintf("SUM(CASE WHEN (%s) IS NULL OR isnan(%s) OR isinf(%s) THEN 1 ELSE 0 END) AS bad_y", z, z, z))
    if (!is.null(w)) {
        items <- c(items,
                   sprintf("SUM(CASE WHEN (%s) IS NULL OR isnan(%s) OR %s < 0 THEN 1 ELSE 0 END) AS bad_w", w, w, w),
                   sprintf("SUM(CASE WHEN isinf(%s) THEN 1 ELSE 0 END) AS inf_w", w),
                   sprintf("SUM(CASE WHEN %s <> 0 THEN 1 ELSE 0 END) AS nz", w),
                   sprintf("SUM(CASE WHEN %s <> 0 THEN ln(%s) ELSE 0 END) AS slw", w, w))
    }
    ck <- dm_query(D, sprintf("SELECT %s FROM (%s) AS s WHERE %s", paste(items, collapse = ", "),
                              D$src$select, D$filter))
    n <- as.numeric(ck$n)
    if (n == 0) stop("0 (non-NA) cases", call. = FALSE)
    if (as.numeric(ck$bad_y) > 0) stop("NA/NaN/Inf in 'y'", call. = FALSE)
    if (!is.null(w)) {
        if (as.numeric(ck$bad_w) > 0) stop("missing or negative weights not allowed", call. = FALSE)
        if (as.numeric(ck$inf_w) > 0) stop("NA/NaN/Inf in 'weights'", call. = FALSE)
        if (as.numeric(ck$nz) == 0) stop("0 (non-NA) cases", call. = FALSE)
    }

    mode <- if (isTRUE(compress)) TRUE else if (isFALSE(compress)) FALSE else dm_compress_mode()
    pt <- dm_remote_patterns(D, z, w, within = TRUE, mode = mode)
    if (!is.null(pt)) {
        pd <- dm_pattern_design(D, pt$P[dm_compress_keys(D)])
        k <- pd$keep
        sw <- pt$sw[k]
        pos <- which(sw > 0)
        xs <- pd$design[pos, ]
        zs <- pt$sy[k, 1L][pos] / sw[pos]
        ws <- sw[pos]
        if (is.null(w) && all(ws == 1)) ws <- NULL
        out <- dmfit_ls(xs, zs, ws, tol, intercept)
        out$rss <- out$rss + sum(pt$within[k])
        rows <- sum(pt$n[k])
        nobs <- sum(pt$nz[k])
        sumlogw <- sum(pt$slw[k])
        out$patterns <- length(pos)
        out$method <- "patterns"
    } else {
        De <- dm_remote_extend(D, c(`__dm_z` = z, if (!is.null(w)) c(`__dm_w` = w)), "y")
        out <- dm_remote_lm_moments(De, list(cols = "__dm_z", response = FALSE),
                                    if (is.null(w)) NULL else "__dm_w", tol, intercept)
        rows <- out$n
        nobs <- if (is.null(w)) rows else as.numeric(ck$nz)
        sumlogw <- if (is.null(w)) 0 else as.numeric(ck$slw)
        out$patterns <- rows
    }
    ## integers, as in lm(), when they fit
    if (nobs < .Machine$integer.max) nobs <- as.integer(nobs)
    out$df.residual <- nobs - out$rank
    if (out$df.residual == 0) out$rss <- 0
    out$sigma <- sqrt(out$rss / out$df.residual)
    out$nobs <- nobs
    out$sumlogw <- sumlogw
    out$weighted <- !is.null(w)
    out$compressed <- !is.null(pt)
    out$rows <- rows
    out$spec <- D$plan
    out$remote <- TRUE
    out$remote_design <- D
    out$remote_z <- z
    out$remote_w <- w
    D$cache$n <- rows
    class(out) <- "design_lm"
    out
}

## The fit from moments, on the centered design when that is exact.
dm_remote_lm_moments <- function(D, yy, weights, tol, intercept) {
    ctr <- if (isTRUE(getOption("design.matrix.remote_center", TRUE))) dm_center(D) else NULL
    D2 <- D
    if (!is.null(ctr)) {
        if (D$engine == "sql") {
            for (nm in names(ctr$shift))
                D2$vars[[nm]]$sql <- sprintf("(%s - %s)", D$vars[[nm]]$sql, dm_sql_double(ctr$shift[[nm]]))
        } else {
            D2$shift <- ctr$shift
        }
    }
    m <- if (D$engine == "sql") dm_sql_moments(D2, yy, weights) else dm_chunked_moments(D2, yy, weights)
    if (!all(is.finite(m$xtx)) || !all(is.finite(m$xty))) stop("NA/NaN/Inf in 'x'", call. = FALSE)
    dg <- diag(m$xtx)
    if (!is.null(ctr)) {
        ## squared column norms of the original design: X = X~ T^-1
        S <- backsolve(ctr$T, diag(ncol(ctr$T)))
        dg <- colSums(S * (unname(m$xtx) %*% S))
    }
    fit <- dmfit_from_moments(m$xtx, m$xty[, 1L], m$yty[1L, 1L], m$sumw, tol, intercept, m$n,
                              icol = which(D$assign == 0L), dg = dg)
    if (!is.null(ctr)) {
        back <- dm_uncenter(fit, ctr$T)
        if (is.null(back)) {
            ## the aliasing pattern does not allow mapping back: plain moments
            m <- if (D$engine == "sql") dm_sql_moments(D, yy, weights) else dm_chunked_moments(D, yy, weights)
            fit <- dmfit_from_moments(m$xtx, m$xty[, 1L], m$yty[1L, 1L], m$sumw, tol, intercept, m$n,
                                      icol = which(D$assign == 0L))
        } else {
            fit <- back
        }
    }
    names(fit$coefficients) <- D$colnames
    colnames(fit$R) <- D$colnames[fit$pivot]
    fit$n <- m$n
    fit$method <- if (D$engine == "sql") "SQL moments" else "moments computed in batches"
    fit
}

## Pivoted Cholesky factor of G = X'WX with the column selection of lm()
## (LINPACK dqrdc2): columns are taken in order, and a column is aliased
## (moved to the end) when the norm of its residual after projection on the
## columns taken before it is less than `tol` times its norm. `dg` holds the
## squared norms of the columns of the original design: those of G, or,
## for the moments of a centered design X~ = X T, those of X, since the
## residual norms are the same for X and X~ while the column norms are not.
## Returns list(R (p x p, columns in pivoted order), rank, pivot).
dm_chol_seq <- function(G, tol, dg = diag(G)) {
    p <- ncol(G)
    keep <- integer()
    al <- integer()
    Rk <- matrix(0, 0L, 0L)
    for (j in seq_len(p)) {
        k <- length(keep)
        r <- if (k) backsolve(Rk, G[keep, j], transpose = TRUE) else numeric()
        d2 <- G[j, j] - sum(r^2)
        if (dg[j] > 0 && d2 > 0 && sqrt(d2) >= tol * sqrt(dg[j])) {
            Rn <- matrix(0, k + 1L, k + 1L)
            Rn[seq_len(k), seq_len(k)] <- Rk
            Rn[seq_len(k), k + 1L] <- r
            Rn[k + 1L, k + 1L] <- sqrt(d2)
            Rk <- Rn
            keep <- c(keep, j)
        } else {
            al <- c(al, j)
        }
    }
    k <- length(keep)
    pivot <- c(keep, al)
    R <- matrix(0, p, p)
    R[seq_len(k), seq_len(k)] <- Rk
    if (k && length(al))
        R[seq_len(k), k + seq_along(al)] <- backsolve(Rk, G[keep, al, drop = FALSE], transpose = TRUE)
    list(R = R, rank = k, pivot = pivot)
}

## A least squares fit from moments: G = X'WX, h = X'Wz, zz = z'Wz, the sum
## of the weights, the number of rows, and the intercept column (for the
## explained sum of squares). The same components as dmfit_ls().
dmfit_from_moments <- function(G, h, zz, sumw, tol, intercept, n, icol = integer(), dg = diag(G)) {
    p <- ncol(G)
    if (p == 0L) {
        return(list(coefficients = numeric(0), rank = 0L, pivot = integer(0), R = matrix(0, 0, 0),
                    tol = tol, rss = zz, mss = 0, intercept = intercept))
    }
    G <- unname(G)
    ch <- dm_chol_seq(G, tol, dg)
    k <- ch$rank
    keep <- ch$pivot[seq_len(k)]
    R11 <- ch$R[seq_len(k), seq_len(k), drop = FALSE]
    b <- dmqr_bs(R11, dmqr_bs(R11, matrix(h[keep], ncol = 1L), transpose = TRUE))
    beta <- numeric(p)
    beta[keep] <- b
    Gb <- drop(G %*% beta)
    fwf <- sum(beta * Gb)
    rss <- max(zz - 2 * sum(beta * h) + fwf, 0)
    mss <- if (intercept && length(icol)) fwf - Gb[icol[1L]]^2 / sumw else fwf
    coef <- beta
    coef[setdiff(seq_len(p), keep)] <- NA_real_
    list(coefficients = coef, rank = k, pivot = ch$pivot, R = ch$R, tol = tol,
         rss = rss, mss = max(mss, 0), intercept = intercept)
}

## Centering of the numeric variables: list(shift (named by variable), T)
## with X(shifted) = X T, T unit upper triangular, or NULL (nothing to
## center, or no such T).
dm_center <- function(D) {
    plan <- D$plan
    sample <- dm_remote_sample(D)
    if (!nrow(sample)) return(NULL)
    mf <- tryCatch(design_frame(plan$terms, sample, xlev = plan$xlevels, na.action = stats::na.omit),
                   error = function(e) NULL)
    if (is.null(mf) || !nrow(mf)) return(NULL)
    cls <- vapply(mf, stats::.MFclass, "")
    num <- names(cls)[cls == "numeric"]
    if (!length(num)) return(NULL)
    mu <- vapply(num, function(nm) mean(mf[[nm]]), 0)
    sdv <- vapply(num, function(nm) stats::sd(mf[[nm]]), 0)
    sdv[!is.finite(sdv)] <- 0
    ## only variables whose mean is large compared with their spread
    big <- is.finite(mu) & abs(mu) > sdv
    if (!any(big)) return(NULL)
    shift <- mu[big]
    p <- D$p
    if (p > 2000L) return(NULL)
    ## prototype rows: sampled rows, with random levels and random values
    m <- max(4L * p, 200L)
    idx <- sample.int(nrow(mf), m, replace = TRUE)
    proto <- mf[idx, , drop = FALSE]
    for (nm in names(proto)) {
        x <- proto[[nm]]
        if (is.factor(x)) {
            proto[[nm]] <- factor(sample(levels(x), m, TRUE), levels = levels(x))
            attr(proto[[nm]], "contrasts") <- attr(x, "contrasts")
        } else if (is.logical(x)) {
            proto[[nm]] <- sample(c(FALSE, TRUE), m, TRUE)
        } else if (nm %in% num) {
            ## (T does not depend on the values: values near zero keep the
            ## prototype design well conditioned)
            proto[[nm]] <- stats::rnorm(m)
        }
    }
    attr(proto, "terms") <- attr(mf, "terms")
    shifted <- proto
    for (nm in names(shift)) shifted[[nm]] <- proto[[nm]] - shift[[nm]]
    X0 <- tryCatch(as.matrix(dm_frame_design(D, proto)), error = function(e) NULL)
    X1 <- tryCatch(as.matrix(dm_frame_design(D, shifted)), error = function(e) NULL)
    if (is.null(X0) || is.null(X1)) return(NULL)
    q <- qr(X0)
    if (q$rank < p) return(NULL)
    Tm <- unname(qr.coef(q, X1))
    scale <- max(1, abs(X1))
    if (!all(is.finite(Tm)) || max(abs(X0 %*% Tm - X1)) > 1e-8 * scale) return(NULL)
    ts <- max(1, abs(Tm))
    if (max(abs(diag(Tm) - 1)) > 1e-8 * ts || max(abs(Tm[lower.tri(Tm)])) > 1e-8 * ts) return(NULL)
    Tm[lower.tri(Tm)] <- 0
    diag(Tm) <- 1
    Tm[abs(Tm) < 1e-10 * ts] <- 0
    list(shift = shift, T = Tm)
}

## The in-memory design of a model frame of the plan's terms.
dm_frame_design <- function(D, mf) {
    plan <- D$plan
    d <- design_matrix(plan$terms, mf, contrasts.arg = plan$contrasts)
    if (isTRUE(plan$drop_intercept)) d <- dm_drop_intercept(d)
    d
}

## Map a fit on the centered design X~ = X T back to X: coefficients
## beta = T beta~ and factor R = R~ T^-1, for the kept columns. NULL when
## the columns kept are not closed under T (a kept column of X~ uses an
## aliased column of X).
dm_uncenter <- function(fit, Tm) {
    p <- ncol(Tm)
    k <- fit$rank
    piv <- fit$pivot
    keep <- piv[seq_len(k)]
    al <- piv[seq_len(p) > k]
    if (!identical(keep, sort(keep))) return(NULL)
    if (length(al) && k && any(Tm[al, keep] != 0)) return(NULL)
    Tk <- Tm[keep, keep, drop = FALSE]
    bt <- fit$coefficients[keep]
    beta <- rep(NA_real_, p)
    beta[keep] <- Tk %*% bt
    R <- fit$R
    if (k) {
        Sk <- backsolve(Tk, diag(k))
        R11 <- R[seq_len(k), seq_len(k), drop = FALSE] %*% Sk
        R[seq_len(k), seq_len(k)] <- R11
        if (length(al)) {
            ## columns of X~ for the aliased columns, through those of X
            S <- backsolve(Tm, diag(p))
            Rt11 <- fit$R[seq_len(k), seq_len(k), drop = FALSE]
            Rt12 <- fit$R[seq_len(k), k + seq_along(al), drop = FALSE]
            cf <- S[keep, al, drop = FALSE] + backsolve(Rt11, Rt12) %*% S[al, al, drop = FALSE]
            R[seq_len(k), k + seq_along(al)] <- R11 %*% Tk %*% cf
        }
    }
    fit$coefficients <- beta
    fit$R <- R
    fit
}

## The offsets of new remote data D for a remote or in-memory fit: offset()
## terms of the plan and the `offset` argument of the fit, resolved
## against the columns of D (SQL), or NULL.
dm_remote_new_offset <- function(object, D) {
    off <- tryCatch(dm_remote_offset_sql(D), error = function(e) stop(conditionMessage(e), call. = FALSE))
    arg <- if (isTRUE(object$remote)) object$remote_offset_arg else object$offset_call
    if (!is.null(arg)) {
        env <- if (isTRUE(object$remote)) object$remote_env else environment(object$terms)
        a <- dm_resolve_arg(arg, env %||% globalenv(), D$src$types, dm_con(D), "offset")
        off <- if (is.null(off)) a$sql else sprintf("%s + (%s)", off, a$sql)
    }
    off
}

## Fitted values and residuals of a remote fit, or predictions on remote
## data, as a table: the columns of the remote data, `.fitted` (with the
## offset), and `.resid` (training rows, with `resid = TRUE`), `.se.fit`,
## `.lwr`, `.upr` and `.nonestim` when asked. With na.exclude, the training
## rows with missing values are kept, with NULL. Standard errors, intervals
## and non-estimability are computed in batches (the chunked engine's
## reader), the rest in SQL when the design's engine is SQL.
dm_remote_lm_table <- function(object, newdata = NULL, na.action = stats::na.pass, dots = list(),
                               resid = FALSE, se.fit = FALSE, interval = "none", level = 0.95,
                               rankdeficient = "warnif", tol = 1e-06) {
    beta <- dmfit_beta(object)
    training <- is.null(newdata)
    if (training) {
        if (!isTRUE(object$remote)) stop("This fit was not computed on remote data.", call. = FALSE)
        D <- object$remote_design
        off <- object$remote_offset
        z <- object$remote_z
        w <- object$remote_w
    } else {
        if (is.null(object$spec)) stop("This design cannot be applied to new data.", call. = FALSE)
        D <- do.call(design_matrix, c(list(object$spec, newdata, na.action = na.action), dots))
        off <- dm_remote_new_offset(object, D)
        z <- w <- NULL
    }
    D <- dm_remote_unscaled(D)
    srccols <- names(D$src$types)
    helpers <- c(`__dm_off` = off, `__dm_z` = z, `__dm_w` = w)
    if (length(helpers)) D <- dm_remote_extend(D, helpers, "offset")
    rankdef <- !training && object$rank < length(object$coefficients)
    if (rankdef && rankdeficient == "simple")
        warning("prediction from rank-deficient fit; consider predict(., rankdeficient=\"NA\")", call. = FALSE)
    batch <- se.fit || interval != "none" || (rankdef && rankdeficient != "simple")
    has <- function(h) h %in% names(helpers)
    if (!batch) {
        v <- dm_remote_predict(D, beta)
        items <- c(sprintf("%s AS \".fitted\"", if (has("__dm_off")) "\".fitted\" + \"__dm_off\"" else "\".fitted\""),
                   if (resid) "\"__dm_z\" - \".fitted\" AS \".resid\"")
        sel <- sprintf("SELECT %s FROM %s", paste(c(dm_qi(D, srccols), items), collapse = ", "), dm_qi(D, unclass(v)[1L]))
    } else {
        res_var <- object$sigma^2
        V <- dmfit_vcov(object, 1)
        V[is.na(V)] <- 0
        tq <- stats::qt((1 - level) / 2, object$df.residual)
        nonest_n <- 0
        if (interval == "prediction") {
            if (training) {
                warning("predictions on current data refer to _future_ responses\n", call. = FALSE)
                if (!is.null(w))
                    warning("assuming prediction variance inversely proportional to weights used for fitting\n",
                            call. = FALSE)
            } else if (object$weighted) {
                warning("Assuming constant prediction variance even though model fit is weighted\n", call. = FALSE)
            }
        }
        out <- dm_tmp_name("dm_pred")
        created <- FALSE
        con <- dm_con(D)
        dm_chunks(D, unique(c(srccols, names(helpers), dm_remote_inputs(D))), function(chunk) {
            cd <- dm_chunk_design(D, chunk)
            X <- cd$design
            k <- cd$keep
            xb <- dmfit_linpred(object, X)
            fit <- if (has("__dm_off")) xb + chunk$`__dm_off`[k] else xb
            res <- chunk[k, srccols, drop = FALSE]
            if (rankdef && rankdeficient != "simple") {
                ne <- seq_along(k) %in% dmfit_nonestimable(object, X, tol)
                nonest_n <<- nonest_n + sum(ne)
                if (startsWith(rankdeficient, "NA")) fit[ne] <- NA else res$.nonestim <- ne
            }
            res$.fitted <- fit
            if (resid) res$.resid <- chunk$`__dm_z`[k] - xb
            if (se.fit || interval != "none") {
                se <- sqrt((if (ncol(X)) quadratic_diagonal(X, V) else rep(0, length(k))) * res_var)
                if (se.fit) res$.se.fit <- se
                if (interval != "none") {
                    pv <- if (interval == "prediction" && training && !is.null(w)) res_var / chunk$`__dm_w`[k] else res_var
                    hw <- tq * if (interval == "confidence") se else sqrt(se^2 + pv)
                    res$.lwr <- fit + hw
                    res$.upr <- fit - hw
                }
            }
            attr(res, "row.names") <- c(NA_integer_, -length(k))
            if (!created) {
                DBI::dbWriteTable(con, out, res, temporary = TRUE)
                created <<- TRUE
            } else {
                DBI::dbAppendTable(con, out, res)
            }
        })
        if (!created) {
            res <- D$src$proto[srccols]
            res$.fitted <- numeric()
            DBI::dbWriteTable(con, out, res, temporary = TRUE)
        }
        if (rankdef && nonest_n > 0) {
            if (rankdeficient == "warnif")
                warning("prediction from rank-deficient fit; attr(*, \"non-estim\") has doubtful cases", call. = FALSE)
            if (rankdeficient == "NAwarn")
                warning("prediction from rank-deficient fit: NAs produced for non-estimable cases", call. = FALSE)
        }
        sel <- sprintf("SELECT * FROM %s", dm_qi(D, out))
    }
    ## na.exclude: the training rows with missing values, with NULL
    if (training && identical(object$remote_na, "exclude")) {
        D0 <- object$remote_design
        sel <- sprintf("%s UNION ALL BY NAME SELECT s.* FROM (%s) AS s WHERE NOT COALESCE(%s, FALSE)",
                       sel, D0$src$select, D0$filter)
    }
    name <- dm_tmp_name("dm_pred")
    dm_execute(D, sprintf("CREATE TEMPORARY VIEW %s AS %s", dm_qi(D, name), sel))
    out <- structure(name, class = "dm_remote_table", connection = dm_con(D), holder = D$src$holder)
    if (se.fit) {
        attr(out, "df") <- object$df.residual
        attr(out, "residual.scale") <- object$sigma
    }
    out
}

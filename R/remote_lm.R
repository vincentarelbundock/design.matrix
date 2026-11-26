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
#' [design_lm()] accepts remote data like [design_matrix()]:
#' `design_lm(y ~ x + f, con, table = "flights", weights = "w")`, with
#' `weights` the name of a column (a string, or a bare name), and a remote
#' design: `design_lm(D, y = "col", weights = "w")`. With few design
#' patterns, the fit is compressed (see *Compression*): coefficients,
#' covariance, residual standard error and degrees of freedom are those of
#' [lm()] on the same rows, up to rounding. Otherwise it is computed from
#' `X'WX`, `X'Wy` and `y'Wy` with [lm()]'s column selection; numeric
#' variables are centered when that is an exact transform of the design,
#' and the accuracy is that of the normal equations on the centered design.
#' The response must have a SQL translation; `subset` and offsets are not
#' supported. `predict(fit, newdata = con, table = ...)` returns the name of
#' a temporary view or table of predictions, as `D %*% beta`; standard
#' errors and intervals are not available there. `residuals(fit)` returns
#' the remote design with two more columns, `.fitted` and `.resid`, which
#' can scale its rows: for heteroskedasticity-consistent and clustered
#' standard errors,
#' ```
#' R <- residuals(fit)
#' meat_hc1 <- crossprod(R, DBI::SQL('".resid" * ".resid"') * R)
#' scores <- colSums(".resid" * R, group = "cluster")
#' meat_cl <- crossprod(scores)
#' ```
#' (with `bread <- vcov(fit) / sigma(fit)^2`, the HC1 covariance is
#' `n / (n - k) * bread %*% meat_hc1 %*% bread`).
NULL

## design_lm(formula, <remote>, ...): `weights` is a column name or NULL.
dm_remote_lm <- function(formula, source, weights, contrasts, compress, tol, call, dots) {
    if (length(attr(stats::terms(formula), "offset")))
        stop("Offsets are not supported for remote data.", call. = FALSE)
    D <- dm_remote(formula, source, contrasts, NULL, dots, weights = weights, drop_unused = TRUE)
    if (is.null(D$response)) stop("The formula has no response.", call. = FALSE)
    fit <- dm_remote_lm_fit(D, list(cols = D$response, response = TRUE), weights, compress, tol)
    fit$call <- call
    fit$terms <- D$terms_full
    fit$contrasts <- D$plan$contrasts
    fit$xlevels <- D$plan$xlevels
    fit
}

## The fit of a remote design, response yy (dm_remote_y()) and weights (a
## column name or NULL).
dm_remote_lm_fit <- function(D, yy, weights, compress = "auto", tol = 1e-07) {
    dmfit_check_compress(compress)
    ys <- dm_compress_ys(D, yy)
    if (is.null(ys) || length(ys) != 1L)
        stop("The response of a remote fit must be one numeric variable with a SQL translation.", call. = FALSE)
    w <- if (is.null(weights)) NULL else sprintf("CAST(s.%s AS DOUBLE)", dm_qi(D, weights))
    intercept <- attr(D$plan$terms, "intercept") > 0L

    ## lm()'s checks, in SQL
    items <- c("COUNT(*) AS n",
               sprintf("SUM(CASE WHEN (%s) IS NULL OR isnan(%s) OR isinf(%s) THEN 1 ELSE 0 END) AS bad_y", ys, ys, ys))
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
    pt <- dm_remote_patterns(D, ys, w, within = TRUE, mode = mode)
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
        out <- dm_remote_lm_moments(D, yy, weights, tol, intercept)
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
    out$remote_y <- if (yy$response) NULL else yy$cols
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

## predict(fit, newdata = <remote>, ...): the linear predictor written to the
## database. `...` holds the source arguments (table, sql, engine,
## chunk_size).
dm_remote_lm_predict <- function(object, newdata, na.action, dots) {
    if (is.null(object$spec)) stop("This design cannot be applied to new data.", call. = FALSE)
    if (!is.null(object$offset_call) || length(attr(object$terms, "offset")))
        stop("Offsets are not supported for remote data.", call. = FALSE)
    D <- do.call(design_matrix, c(list(object$spec, newdata, na.action = na.action), dots))
    if (object$rank < length(object$coefficients))
        warning("prediction from rank-deficient fit; consider predict(., rankdeficient=\"NA\")", call. = FALSE)
    dm_remote_predict(D, dmfit_beta(object))
}

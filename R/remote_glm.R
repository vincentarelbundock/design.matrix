## design_glm() on remote data (see glm.R for the iterations, remote.R for
## remote designs).
##
## The response, prior weights and offset become columns of the source,
## computed in SQL: __dm_y and __dm_w are y and the weights after
## family$initialize (proportions and weights times trials for a
## two-column binomial response, y = 0 where the weight is 0), __dm_n the
## number of trials, __dm_off the offset. Rows where one of them is missing
## are dropped, as na.omit() drops them.
##
## Three ways to fit, all running the loop of glm.R (dmglm_loop()):
##
##  * patterns (dmglm_remote_patterns()): one GROUP BY over the design's
##    input columns and the offset returns, per pattern, the sums the
##    compressed fit needs (glm.R): sum of weights, weighted sum of y, the
##    beta-free part of the deviance, the sums of the AIC, and the working
##    weights and responses of the first iteration computed per row from
##    mustart. The iterations then run in memory on the patterns;
##  * SQL engine (dmglm_sql_backend()): each iteration computes, in the
##    database, the linear predictor from the current coefficients (the
##    prediction joins of remote_sql.R), mu, the working weights and
##    responses (the family translated to SQL), then X'WX and X'Wz with
##    those weights (dm_sql_gram(), dm_sql_xty_core()) and the deviance;
##  * chunked engine (dmglm_chunk_backend()): each iteration reads the rows
##    in batches, computes the same in R, and accumulates X'WX and X'Wz.
## The last two solve the normal equations with lm()'s column selection
## (dm_chol_seq()), on the design with numeric variables centered when that
## is exact (as remote_lm.R).

#' @name design_remote
#' @section GLMs:
#' [design_glm()] accepts remote data like [design_lm()]:
#' `design_glm(y ~ x + f, binomial, con, table = "d", weights = "w", offset = "o")`,
#' with `weights` and `offset` naming columns, `offset()` terms in the
#' formula, and a remote design with `y` naming the response column (two
#' columns for successes and failures). The response, weights and offset
#' must have SQL translations. Families and links are those of
#' [stats::family] (every one has a SQL translation). With few design
#' patterns (see *Compression*), one `GROUP BY` returns the sums of each
#' pattern and the iterations run in memory, with the results of the fit
#' on every row. Otherwise each iteration computes `X'WX`, `X'Wz` and the
#' deviance in the database (SQL engine) or in batches (chunked engine),
#' from the normal equations, with aliased columns detected at the
#' tolerance `max(1e-7, min(1e-7, epsilon / 1000))` (that of [lm()]: the
#' normal equations cannot resolve `glm()`'s smaller tolerance).
#' `predict(fit, newdata = con, table = ..., type = "response")`, `fitted()`
#' and `residuals(fit, type = )` return the name of a temporary view or
#' table with a `.fitted` (and `.resid`) column next to the columns of the
#' remote data.
NULL

## ---- families -------------------------------------------------------------

dmglm_link_names <- c("identity", "log", "logit", "probit", "cloglog", "cauchit", "inverse", "sqrt", "1/mu^2")

## What design_glm() supports, in memory and remote alike: the families of
## stats and quasi(), with their links. Returns the description of the
## family used for SQL translations and remote checks.
dmglm_check_family <- function(family) {
    fam <- family$family
    link <- family$link
    ok_fam <- c("gaussian", "binomial", "poisson", "Gamma", "inverse.gaussian",
                "quasibinomial", "quasipoisson", "quasi")
    if (!isTRUE(fam %in% ok_fam) || !isTRUE(link %in% dmglm_link_names))
        stop(sprintf("design_glm() supports the families %s and the links %s (family `%s`, link `%s`).",
                     paste(ok_fam, collapse = ", "), paste(dmglm_link_names, collapse = ", "),
                     format(fam), format(link)), call. = FALSE)
    k <- switch(fam,
                gaussian = list(dev = "gaussian", init = "y", var = "1", validmu = "none"),
                binomial = , quasibinomial = list(dev = "binomial", init = "binomial", var = "mu(1-mu)",
                                                  validmu = "unit"),
                poisson = , quasipoisson = list(dev = "poisson", init = "poisson", var = "mu",
                                                validmu = "positive"),
                Gamma = list(dev = "gamma", init = "y", var = "mu^2", validmu = "positive"),
                inverse.gaussian = list(dev = "invgauss", init = "y", var = "mu^3", validmu = "none"),
                quasi = {
                    v <- family$varfun
                    switch(v,
                           constant = list(dev = "gaussian", init = "y", var = "1", validmu = "none"),
                           "mu(1-mu)" = list(dev = "binomial", init = "clamp", var = "mu(1-mu)", validmu = "unit_nf"),
                           mu = list(dev = "poisson", init = "y01", var = "mu", validmu = "positive_nf"),
                           "mu^2" = list(dev = "qgamma", init = "y01", var = "mu^2", validmu = "positive_nf"),
                           "mu^3" = list(dev = "invgauss", init = "y01", var = "mu^3", validmu = "positive_nf"),
                           stop(sprintf("design_glm() does not support quasi(variance = \"%s\").", v), call. = FALSE))
                })
    k$family <- fam
    k$link <- link
    k$binomial <- fam %in% c("binomial", "quasibinomial")
    k$aic <- if (fam %in% c("quasibinomial", "quasipoisson", "quasi")) "none" else fam
    k
}

## SQL expressions of the family, as functions of SQL expressions.
dmglm_sql_family <- function(k) {
    eps <- dm_sql_double(.Machine$double.eps)
    num <- function(x) dm_sql_double(x)
    link <- switch(k$link,
        identity = list(fun = function(m) m, inv = function(e) e, me = function(e) "CAST(1 AS DOUBLE)",
                        valid = function(e) "TRUE"),
        log = list(fun = function(m) sprintf("ln(%s)", m),
                   inv = function(e) sprintf("greatest(exp(%s), %s)", e, eps),
                   me = function(e) sprintf("greatest(exp(%s), %s)", e, eps),
                   valid = function(e) "TRUE"),
        logit = list(fun = function(m) sprintf("ln(%s / (1 - %s))", m, m),
                     inv = function(e) sprintf(paste("(CASE WHEN %s < -30 THEN %s WHEN %s > 30 THEN %s",
                                                     "ELSE exp(%s) / (1 + exp(%s)) END)"),
                                               e, num(.Machine$double.eps / (1 + .Machine$double.eps)), e,
                                               num((1 / .Machine$double.eps) / (1 + 1 / .Machine$double.eps)), e, e),
                     me = function(e) sprintf("(CASE WHEN %s > 30 OR %s < -30 THEN %s ELSE exp(%s) / ((1 + exp(%s)) * (1 + exp(%s))) END)",
                                              e, e, eps, e, e, e),
                     valid = function(e) "TRUE"),
        probit = {
            thr <- num(-stats::qnorm(.Machine$double.eps))
            list(fun = function(m) sprintf("dm_glm_qnorm(%s)", m),
                 inv = function(e) sprintf("dm_glm_pnorm(least(greatest(%s, -%s), %s))", e, thr, thr),
                 me = function(e) sprintf("greatest(dm_glm_dnorm(%s), %s)", e, eps),
                 valid = function(e) "TRUE")
        },
        cloglog = list(fun = function(m) sprintf("ln(-ln(1 - %s))", m),
                       inv = function(e) sprintf("greatest(least(dm_glm_mexpm1m(exp(%s)), %s), %s)",
                                                 e, num(1 - .Machine$double.eps), eps),
                       me = function(e) sprintf("greatest(exp(least(%s, 700)) * exp(-exp(least(%s, 700))), %s)", e, e, eps),
                       valid = function(e) "TRUE"),
        cauchit = {
            thr <- num(-stats::qcauchy(.Machine$double.eps))
            list(fun = function(m) sprintf("(CASE WHEN %s <= 0.5 THEN -1 / tan(pi() * %s) ELSE 1 / tan(pi() * (1 - %s)) END)", m, m, m),
                 inv = function(e) sprintf("dm_glm_pcauchy(least(greatest(%s, -%s), %s))", e, thr, thr),
                 me = function(e) sprintf("greatest(1 / (pi() * (1 + (%s) * (%s))), %s)", e, e, eps),
                 valid = function(e) "TRUE")
        },
        inverse = list(fun = function(m) sprintf("(1 / %s)", m), inv = function(e) sprintf("(1 / %s)", e),
                       me = function(e) sprintf("(-1 / ((%s) * (%s)))", e, e),
                       valid = function(e) sprintf("(isfinite(%s) AND %s <> 0)", e, e)),
        sqrt = list(fun = function(m) sprintf("sqrt(%s)", m), inv = function(e) sprintf("((%s) * (%s))", e, e),
                    me = function(e) sprintf("(2 * %s)", e),
                    valid = function(e) sprintf("(isfinite(%s) AND %s > 0)", e, e)),
        "1/mu^2" = list(fun = function(m) sprintf("(1 / ((%s) * (%s)))", m, m), inv = function(e) sprintf("(1 / sqrt(%s))", e),
                        me = function(e) sprintf("(-1 / (2 * pow(%s, 1.5)))", e),
                        valid = function(e) sprintf("(isfinite(%s) AND %s > 0)", e, e)))
    variance <- switch(k$var, "1" = function(m) "CAST(1 AS DOUBLE)",
                       "mu(1-mu)" = function(m) sprintf("(%s * (1 - %s))", m, m),
                       mu = function(m) m,
                       "mu^2" = function(m) sprintf("(%s * %s)", m, m),
                       "mu^3" = function(m) sprintf("(%s * %s * %s)", m, m, m))
    validmu <- switch(k$validmu, none = function(m) "TRUE",
                      unit = function(m) sprintf("(isfinite(%s) AND %s > 0 AND %s < 1)", m, m, m),
                      unit_nf = function(m) sprintf("(NOT isnan(%s) AND %s > 0 AND %s < 1)", m, m, m),
                      positive = function(m) sprintf("(isfinite(%s) AND %s > 0)", m, m),
                      positive_nf = function(m) sprintf("(NOT isnan(%s) AND %s > 0)", m, m))
    dev <- switch(k$dev,
        gaussian = function(y, m, w) sprintf("(%s * (%s - %s) * (%s - %s))", w, y, m, y, m),
        binomial = function(y, m, w) sprintf(paste("(2 * %s * ((CASE WHEN %s <> 0 THEN %s * ln(%s / %s) ELSE 0 END) +",
                                                   "(CASE WHEN (1 - %s) <> 0 THEN (1 - %s) * ln((1 - %s) / (1 - %s)) ELSE 0 END)))"),
                                             w, y, y, y, m, y, y, y, m),
        poisson = function(y, m, w) sprintf("(CASE WHEN %s > 0 THEN 2 * %s * (%s * ln(%s / %s) - (%s - %s)) ELSE 2 * %s * %s END)",
                                            y, w, y, y, m, y, m, m, w),
        gamma = function(y, m, w) sprintf("(-2 * %s * (ln(CASE WHEN %s = 0 THEN 1 ELSE %s / %s END) - (%s - %s) / %s))",
                                          w, y, y, m, y, m, m),
        qgamma = function(y, m, w) sprintf("greatest(-2 * %s * (ln((CASE WHEN %s = 0 THEN 1 ELSE %s END) / %s) - (%s - %s) / %s), 0)",
                                           w, y, y, m, y, m, m),
        invgauss = function(y, m, w) sprintf("(%s * (%s - %s) * (%s - %s) / (%s * %s * %s))", w, y, m, y, m, y, m, m))
    mustart <- switch(k$init,
        y = function(y, w, n, cb) y,
        binomial = function(y, w, n, cb) if (cb) sprintf("((%s * %s + 0.5) / (%s + 1))", n, y, n)
                                         else sprintf("((%s * %s + 0.5) / (%s + 1))", w, y, w),
        poisson = function(y, w, n, cb) sprintf("(%s + 0.1)", y),
        y01 = function(y, w, n, cb) sprintf("(%s + 0.1 * (CASE WHEN %s = 0 THEN 1 ELSE 0 END))", y, y),
        clamp = function(y, w, n, cb) sprintf("greatest(0.001, least(0.999, %s))", y))
    list(linkfun = link$fun, linkinv = link$inv, mu.eta = link$me, valideta = link$valid,
         variance = variance, validmu = validmu, dev = dev, mustart = mustart)
}

## The same starting values in R (for the chunked engine).
dmglm_mustart_r <- function(k, y, w, n, cb) {
    switch(k$init, y = y,
           binomial = if (cb) (n * y + 0.5) / (n + 1) else (w * y + 0.5) / (w + 1),
           poisson = y + 0.1, y01 = y + 0.1 * (y == 0), clamp = pmax(0.001, pmin(0.999, y)))
}

## Temporary macros: R's pnorm (Cody's algorithm, as pnorm.c), qnorm
## (Wichura's AS241, as qnorm.c), dnorm, pcauchy and -expm1(-u).
dmglm_sql_macros <- function(D) {
    if (isTRUE(D$src$holder$glm_macros)) return(invisible())
    horner <- function(x, co) {
        s <- sprintf("%.17g", co[1L])
        for (k in co[-1L]) s <- sprintf("(%s * %s + %.17g)", s, x, k)
        s
    }
    a <- c(2.2352520354606839287, 161.02823106855587881, 1067.6894854603709582, 18154.981253343561249,
           0.065682337918207449113)
    b <- c(47.20258190468824187, 976.09855173777669322, 10260.932208618978205, 45507.789335026729956)
    cc <- c(0.39894151208813466764, 8.8831497943883759412, 93.506656132177855979, 597.27027639480026226,
            2494.5375852903726711, 6848.1904505362823326, 11602.651437647350124, 9842.7148383839780218,
            1.0765576773720192317e-8)
    d <- c(22.266688044328115691, 235.38790178262499861, 1519.377599407554805, 6485.558298266760755,
           18615.571640885098091, 34900.952721145977266, 38912.003286093271411, 19685.429676859990727)
    p <- c(0.21589853405795699, 0.1274011611602473639, 0.022235277870649807, 0.001421619193227893466,
           2.9112874951168792e-5, 0.02307344176494017303)
    q <- c(1.28426009614491121, 0.468238212480865118, 0.0659881378689285515, 0.00378239633202758244,
           7.29751555083966205e-5)
    num1 <- sprintf("(x * %s)", horner("x*x", c(a[5L], a[1:4])))
    den1 <- horner("x*x", c(1, b[1:4]))
    num2 <- horner("abs(x)", c(cc[9L], cc[1:8]))
    den2 <- horner("abs(x)", c(1, d[1:8]))
    ix <- "(1.0/(x*x))"
    num3 <- sprintf("(%s * %s)", ix, horner(ix, c(p[6L], p[1:5])))
    den3 <- horner(ix, c(1, q[1:5]))
    xs <- "(trunc(abs(x)*16)/16)"
    tail <- function(temp) sprintf("(exp(-%s*%s*0.5) * exp(-((abs(x) - %s)*(abs(x) + %s))*0.5) * %s)", xs, xs, xs, xs, temp)
    t2 <- tail(sprintf("(%s / %s)", num2, den2))
    t3 <- tail(sprintf("((0.398942280401432677939946059934 - %s / %s) / abs(x))", num3, den3))
    pn <- sprintf(paste("CASE WHEN abs(x) <= 0.67448975 THEN 0.5 + %s / %s",
                        "WHEN abs(x) <= 5.656854249492380195206754896838 THEN (CASE WHEN x > 0 THEN 1 - %s ELSE %s END)",
                        "WHEN abs(x) < 50 THEN (CASE WHEN x > 0 THEN 1 - %s ELSE %s END)",
                        "ELSE (CASE WHEN x > 0 THEN 1.0 ELSE 0.0 END) END"), num1, den1, t2, t2, t3, t3)
    qa <- c(2509.0809287301226727, 33430.575583588128105, 67265.770927008700853, 45921.953931549871457,
            13731.693765509461125, 1971.5909503065514427, 133.14166789178437745, 3.387132872796366608)
    qb <- c(5226.495278852545925, 28729.085735721942674, 39307.89580009271061, 21213.794301586595867,
            5394.1960214247511077, 687.1870074920579083, 42.313330701600911252, 1)
    qc <- c(7.7454501427834140764e-4, .0227238449892691845833, .24178072517745061177, 1.27045825245236838258,
            3.64784832476320460504, 5.7694972214606914055, 4.6303378461565452959, 1.42343711074968357734)
    qd <- c(1.05075007164441684324e-9, 5.475938084995344946e-4, .0151986665636164571966, .14810397642748007459,
            .68976733498510000455, 1.6763848301838038494, 2.05319162663775882187, 1)
    qe <- c(2.01033439929228813265e-7, 2.71155556874348757815e-5, .0012426609473880784386, .026532189526576123093,
            .29656057182850489123, 1.7848265399172913358, 5.4637849111641143699, 6.6579046435011037772)
    qf <- c(2.04426310338993978564e-15, 1.4215117583164458887e-7, 1.8463183175100546818e-5, 7.868691311456132591e-4,
            .0148753612908506148525, .13692988092273580531, .59983220655588793769, 1)
    r1 <- "(0.180625 - (p - 0.5)*(p - 0.5))"
    rr <- "sqrt(-ln(least(p, 1 - p)))"
    r2 <- sprintf("(%s - 1.6)", rr)
    r3 <- sprintf("(%s - 5)", rr)
    qn <- sprintf(paste("CASE WHEN abs(p - 0.5) <= 0.425 THEN (p - 0.5) * %s / %s",
                        "ELSE (CASE WHEN p < 0.5 THEN -1 ELSE 1 END) * (CASE WHEN %s <= 5 THEN %s / %s ELSE %s / %s END) END"),
                  horner(r1, qa), horner(r1, qb), rr, horner(r2, qc), horner(r2, qd), horner(r3, qe), horner(r3, qf))
    ## dnorm as dnorm.c: for |x| >= 5, exp() of a split square
    x1 <- "(round(abs(x) * 65536) / 65536)"
    dn <- sprintf(paste("CASE WHEN abs(x) < 5 THEN 0.398942280401432677939946059934 * exp(-0.5 * x * x)",
                        "ELSE 0.398942280401432677939946059934 * exp(-0.5 * %s * %s) * exp((-0.5 * (abs(x) - %s) - %s) * (abs(x) - %s)) END"),
                  x1, x1, x1, x1, x1)
    pc <- "CASE WHEN abs(x) > 1 THEN (CASE WHEN x > 0 THEN 1 - atan(1 / x) / pi() ELSE -atan(1 / x) / pi() END) ELSE 0.5 + atan(x) / pi() END"
    ## -expm1(-u) for u >= 0, Kahan's trick: (1 - w) u / -ln(w), w = exp(-u)
    em <- "CASE WHEN u > 40 THEN 1.0 WHEN exp(-u) = 1 THEN u ELSE (1 - exp(-u)) * u / (-ln(exp(-u))) END"
    for (m in list(c("dm_glm_pnorm", "x", pn), c("dm_glm_qnorm", "p", qn), c("dm_glm_dnorm", "x", dn),
                   c("dm_glm_pcauchy", "x", pc), c("dm_glm_mexpm1m", "u", em)))
        dm_execute(D, sprintf("CREATE OR REPLACE TEMP MACRO %s(%s) AS (%s)", m[1L], m[2L], m[3L]))
    D$src$holder$glm_macros <- TRUE
    invisible()
}

## ---- entry points -----------------------------------------------------------

## design_glm(formula, family, <remote>, ...): `weights` and `offset` are
## column names or NULL.
dm_remote_glm <- function(formula, source, family, weights, offset, start, contrasts, compress, control,
                          call, dots, na_fail = FALSE) {
    dm_need_duckdb()
    k <- dmglm_check_family(family)
    nms <- names(dots) %||% rep("", length(dots))
    src <- dm_remote_source(source, dots$table, dots$sql)
    con <- src$holder$con
    tt <- stats::terms(formula, data = src$proto)
    if (attr(tt, "response") == 0L) stop("The formula has no response.", call. = FALSE)
    vars <- attr(tt, "variables")
    lhs <- vars[[attr(tt, "response") + 1L]]
    off_exprs <- lapply(attr(tt, "offset"), function(i) vars[[i + 1L]][[2L]])
    labs <- attr(tt, "term.labels")
    int <- attr(tt, "intercept") > 0L
    rhs <- if (length(labs)) stats::reformulate(labs, intercept = int) else if (int) ~1 else ~0
    environment(rhs) <- environment(formula)
    cols <- dmglm_remote_cols(con, src$types, src$proto, src$select, lhs, weights, off_exprs, offset, k)
    sel <- sprintf("SELECT s.*, %s FROM (%s) AS s WHERE %s",
                   paste(sprintf("%s AS %s", cols$exprs, DBI::dbQuoteIdentifier(con, names(cols$exprs))), collapse = ", "),
                   src$select, cols$notnull)
    D <- dm_remote(rhs, con, contrasts, NULL, c(list(sql = sel), dots[!nms %in% c("table", "sql")]),
                   drop_unused = TRUE)
    ## keep a private connection (Parquet) open as long as the design
    D$src$holder <- src$holder
    if (na_fail) {
        tot <- as.numeric(dm_query(D, sprintf("SELECT COUNT(*) AS n FROM (%s) AS s", src$select))$n)
        if (tot != D$n) stop("missing values in object", call. = FALSE)
    }
    fit <- dmglm_remote_fit(D, cols, k, family, start, control, compress, int, wmsg = "negative weights not allowed")
    fit$call <- call
    fit$terms <- tt
    fit$formula <- formula
    fit$contrasts <- D$plan$contrasts
    fit$xlevels <- D$plan$xlevels
    if (!is.null(offset)) fit$offset_call <- as.name(offset)
    fit
}

## design_glm(<remote design>, y = , weights = , offset = ).
dm_remote_glm_design <- function(D, y, weights, offset, family, start, control, compress) {
    k <- dmglm_check_family(family)
    if (!is.null(D$scale)) stop("Fit an unscaled remote design, with `weights` naming a column.", call. = FALSE)
    if (!identical(D$na, "omit"))
        stop("Fit a remote design built with `na.action = na.omit` (or from a formula).", call. = FALSE)
    if (!is.null(weights)) dm_remote_column(D, weights, "weights")
    if (!is.null(offset)) dm_remote_column(D, offset, "offset")
    lhs <- if (!is.null(y)) {
        dm_remote_column(D, y, "y", several = TRUE)
        if (length(y) == 2L) call("cbind", as.name(y[1L]), as.name(y[2L]))
        else if (length(y) == 1L) as.name(y)
        else stop("`y` must name one column, or two (successes and failures).", call. = FALSE)
    } else if (!is.null(D$response)) {
        D$vars[[D$response]]$expr
    } else {
        stop("With a remote design, give `y`, the name of the response column (or use a design with a response).",
             call. = FALSE)
    }
    cols <- dmglm_remote_cols(dm_con(D), D$src$types, D$src$proto, D$src$select, lhs, weights, list(), offset, k)
    D <- dm_remote_extend(D, cols$exprs, "response")
    D$filter <- sprintf("(%s) AND %s", D$filter, cols$notnull)
    D$cache <- new.env(parent = emptyenv())
    D$n <- as.numeric(dm_query(D, sprintf("SELECT COUNT(*) AS n FROM (%s) AS s WHERE %s", D$src$select, D$filter))$n)
    intercept <- attr(D$plan$terms, "intercept") > 0L && !isTRUE(D$plan$drop_intercept)
    dmglm_remote_fit(D, cols, k, family, start, control, compress, intercept,
                     wmsg = "missing or negative weights not allowed")
}

## SQL of the response, weights, trials and offset (over the columns of
## `s`), and the condition that none is missing.
dmglm_remote_cols <- function(con, types, proto, select, lhs, weights, off_exprs, off_col, k) {
    tr <- function(e, what) {
        r <- dm_sql_num(e, types, con)
        if (!is.null(r$reason))
            stop(sprintf("The %s `%s` has no SQL translation (%s).", what, deparse2(e), r$reason), call. = FALSE)
        r$sql
    }
    nn <- function(x) sprintf("((%s) IS NOT NULL AND NOT isnan(%s))", x, x)
    notnull <- character()
    raw <- list()
    w <- "CAST(1 AS DOUBLE)"
    if (!is.null(weights)) {
        if (!isTRUE(types[weights] %in% c("double", "integer", "logical")))
            stop("'weights' must be a numeric vector", call. = FALSE)
        w <- sprintf("CAST(s.%s AS DOUBLE)", DBI::dbQuoteIdentifier(con, weights))
        notnull <- c(notnull, nn(w))
    }
    cb <- is.call(lhs) && identical(lhs[[1L]], as.name("cbind"))
    n <- "CAST(1 AS DOUBLE)"
    if (cb) {
        if (length(lhs) != 3L || !k$binomial)
            stop("A two-column response (`cbind(successes, failures)`) needs a binomial family.", call. = FALSE)
        a <- tr(lhs[[2L]], "response")
        b <- tr(lhs[[3L]], "response")
        notnull <- c(notnull, nn(a), nn(b))
        n <- sprintf("(%s + %s)", a, b)
        y <- sprintf("(CASE WHEN %s = 0 THEN 0 ELSE %s / %s END)", n, a, n)
        w <- sprintf("(%s * %s)", w, n)
        raw <- list(a = a, b = b)
    } else if (is.name(lhs) && isTRUE(types[as.character(lhs)] %in% c("character", "factor"))) {
        ## a factor response: success is any level but the first, as in binomial()
        if (!k$binomial) stop("A factor response needs a binomial family.", call. = FALSE)
        col <- as.character(lhs)
        q <- DBI::dbQuoteIdentifier(con, col)
        vals <- DBI::dbGetQuery(con, sprintf("SELECT DISTINCT s.%s AS v FROM (%s) AS s WHERE s.%s IS NOT NULL",
                                             q, select, q))$v
        f <- if (is.factor(proto[[col]])) factor(as.character(vals), levels = levels(proto[[col]])) else factor(vals)
        first <- levels(droplevels(f))[1L]
        y <- sprintf("(CASE WHEN CAST(s.%s AS VARCHAR) = %s THEN 0.0 ELSE 1.0 END)", q, DBI::dbQuoteString(con, first))
        notnull <- c(notnull, sprintf("s.%s IS NOT NULL", q))
    } else {
        y <- tr(lhs, "response")
        notnull <- c(notnull, nn(y))
    }
    ## binomial()'s initialize: y = 0 where the weight is 0
    if (k$binomial && !cb) y <- sprintf("(CASE WHEN %s = 0 THEN 0 ELSE %s END)", w, y)
    off <- NULL
    pieces <- vapply(off_exprs, function(e) tr(e, "offset"), "")
    if (!is.null(off_col)) {
        if (!isTRUE(types[off_col] %in% c("double", "integer", "logical")))
            stop("The offset must be a numeric column.", call. = FALSE)
        pieces <- c(pieces, sprintf("CAST(s.%s AS DOUBLE)", DBI::dbQuoteIdentifier(con, off_col)))
    }
    if (length(pieces)) {
        off <- paste0("(", paste(pieces, collapse = " + "), ")")
        notnull <- c(notnull, nn(off))
    }
    exprs <- c(`__dm_y` = y, `__dm_w` = w, `__dm_n` = n)
    if (!is.null(off)) exprs <- c(exprs, `__dm_off` = off)
    if (cb) exprs <- c(exprs, `__dm_ca` = raw$a, `__dm_cb` = raw$b)
    list(exprs = exprs, notnull = if (length(notnull)) paste(notnull, collapse = " AND ") else "TRUE",
         cbind = cb, offset = !is.null(off), weighted = !is.null(weights))
}

## ---- the fit ----------------------------------------------------------------

dmglm_remote_fit <- function(D, cols, k, family, start, control, compress, intercept, wmsg) {
    dmfit_check_compress(compress)
    if (identical(k$link, "probit") || identical(k$link, "cloglog") || identical(k$link, "cauchit"))
        dmglm_sql_macros(D)
    sf <- dmglm_sql_family(k)
    tol <- min(1e-07, control$epsilon / 1000)
    p <- D$p
    if (!is.null(start) && length(start) != p)
        stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
                      p, paste(deparse(D$colnames), collapse = ", ")), call. = FALSE)
    cx <- dmglm_remote_context(D, cols, k, sf, start, wmsg)
    cx$family <- family
    cx$control <- control
    mode <- if (isTRUE(compress)) TRUE else if (isFALSE(compress)) FALSE else dm_compress_mode()
    res <- if (p > 0L) dmglm_remote_patterns(D, cx, start, control, tol, intercept, mode, family) else NULL
    if (is.null(res)) {
        res <- if (D$engine == "sql") {
            dmglm_remote_iterate(D, cx, start, control, tol, intercept, family, "sql")
        } else {
            dmglm_remote_iterate(D, cx, start, control, tol, intercept, family, "chunked")
        }
    }
    nobs <- cx$ck$nz
    if (nobs < .Machine$integer.max) nobs <- as.integer(nobs)
    rank <- res$fit$rank
    out <- list(coefficients = stats::setNames(res$coef, D$colnames), rank = rank, pivot = res$fit$pivot,
                R = res$fit$R, tol = tol, deviance = res$dev, null.deviance = res$nulldev,
                aic = res$aic + 2 * rank, iter = res$iter, converged = res$conv, boundary = res$boundary,
                df.residual = nobs - rank, df.null = nobs - as.integer(intercept), nobs = nobs,
                rows = cx$ck$n, pearson = res$pearson, zero_weights = res$zero, family = family,
                intercept = intercept, prior_weighted = cols$weighted, compressed = res$compressed,
                patterns = res$patterns, control = control, spec = D$plan, remote = TRUE,
                method = res$method, remote_design = D, remote_cols = cols)
    class(out) <- "design_glm"
    out
}

## Columns as SQL over a select of the source, the checks of the family's
## initialize, and the sums every path needs.
dmglm_remote_context <- function(D, cols, k, sf, start, wmsg) {
    Y <- "__dm_y"
    W <- "__dm_w"
    from <- sprintf("(SELECT * FROM (%s) AS s WHERE %s) AS t", D$src$select, D$filter)
    items <- c("COUNT(*) AS n", sprintf("SUM(CASE WHEN %s <> 0 THEN 1 ELSE 0 END) AS nz", W),
               sprintf("SUM(%s) AS sw", W), sprintf("SUM(%s * %s) AS swy", W, Y), sprintf("AVG(%s) AS ym", Y),
               "MAX(__dm_n) AS maxn",
               sprintf("SUM(CASE WHEN %s > 0 THEN ln(%s) ELSE 0 END) AS slw", W, W),
               sprintf("SUM(CASE WHEN %s = 0 THEN 1 ELSE 0 END) AS nzero", W),
               sprintf("SUM(CASE WHEN %s < 0 THEN 1 ELSE 0 END) AS neg", W),
               sprintf("SUM(CASE WHEN isinf(%s) THEN 1 ELSE 0 END) AS infw", W),
               sprintf("SUM(CASE WHEN %s > 0 THEN %s * ln(%s) ELSE 0 END) AS ly", Y, W, Y),
               sprintf("SUM(CASE WHEN %s <> round(%s) THEN 1 ELSE 0 END) AS nonint", Y, Y))
    chk <- list()
    add <- function(cond, msg, type = "error") {
        chk[[length(chk) + 1L]] <<- list(cond = cond, msg = msg, type = type)
    }
    fam <- k$family
    if (k$binomial && !cols$cbind) add(sprintf("%s < 0 OR %s > 1", Y, Y), "y values must be 0 <= y <= 1")
    if (fam == "binomial" && !cols$cbind)
        add(sprintf("abs(%s * %s - round(%s * %s)) > 0.001", W, Y, W, Y), "non-integer #successes in a binomial glm!", "warning")
    if (fam == "binomial" && cols$cbind)
        add("abs(__dm_ca - round(__dm_ca)) > 0.001 OR abs(__dm_cb - round(__dm_cb)) > 0.001",
            "non-integer counts in a binomial glm!", "warning")
    if (fam == "poisson") add(sprintf("%s < 0", Y), "negative values not allowed for the 'Poisson' family")
    if (fam == "quasipoisson") add(sprintf("%s < 0", Y), "negative values not allowed for the 'quasiPoisson' family")
    if (fam == "Gamma") add(sprintf("%s <= 0", Y), "non-positive values not allowed for the 'Gamma' family")
    if (fam == "inverse.gaussian") add(sprintf("%s <= 0", Y), "positive values only are allowed for the 'inverse.gaussian' family")
    if (fam == "gaussian" && is.null(start) && k$link %in% c("inverse", "log"))
        add(if (k$link == "inverse") sprintf("%s = 0", Y) else sprintf("%s <= 0", Y),
            "cannot find valid starting values: please specify some")
    items <- c(items, sprintf("SUM(CASE WHEN %s THEN 1 ELSE 0 END) AS chk%d",
                              vapply(chk, function(x) x$cond, ""), seq_along(chk)))
    r <- dm_query(D, sprintf("SELECT %s FROM %s", paste(items, collapse = ", "), from))
    ck <- lapply(r, function(x) {
        x <- as.numeric(x)
        if (is.na(x)) 0 else x
    })
    if (ck$n == 0) stop("0 (non-NA) cases", call. = FALSE)
    if (cols$weighted) {
        if (ck$neg > 0) stop(wmsg, call. = FALSE)
        if (ck$infw > 0) stop("NA/NaN/Inf in 'weights'", call. = FALSE)
    }
    for (j in seq_along(chk)) {
        if (ck[[paste0("chk", j)]] > 0) {
            if (chk[[j]]$type == "error") stop(chk[[j]]$msg, call. = FALSE) else warning(chk[[j]]$msg, call. = FALSE)
        }
    }
    if (ck$sw == 0 && k$binomial) stop("0 (non-NA) cases", call. = FALSE)
    list(ck = ck, k = k, sf = sf, cols = cols, off = if (cols$offset) "__dm_off" else "CAST(0 AS DOUBLE)",
         m = if (ck$maxn > 1) "__dm_n" else "__dm_w")
}

## The working select: from a select with __dm_eta (and the source
## columns), add mu, d mu / d eta, V(mu), and the working weight and
## response of IRLS.
dmglm_sql_work <- function(base, cx) {
    sf <- cx$sf
    l2 <- sprintf("SELECT *, %s AS __dm_mu, %s AS __dm_me FROM (%s) AS t2",
                  sf$linkinv("__dm_eta"), sf$mu.eta("__dm_eta"), base)
    l3 <- sprintf("SELECT *, %s AS __dm_var FROM (%s) AS t3", sf$variance("__dm_mu"), l2)
    good <- "__dm_w > 0 AND __dm_me <> 0"
    sprintf(paste("SELECT *, CASE WHEN %s THEN __dm_w * __dm_me * __dm_me / __dm_var ELSE 0 END AS __dm_wk,",
                  "CASE WHEN %s THEN (__dm_eta - %s) + (__dm_y - __dm_mu) / __dm_me ELSE 0 END AS __dm_z",
                  "FROM (%s) AS t4"), good, good, cx$off, l3)
}

## Aggregates over a working select: deviance, validity, range of mu and
## the checks of glm.fit() at each iteration.
dmglm_sql_eval_items <- function(cx) {
    sf <- cx$sf
    c(sprintf("SUM(%s) AS dev", sf$dev("__dm_y", "__dm_mu", "__dm_w")),
      sprintf("SUM(CASE WHEN %s AND %s THEN 0 ELSE 1 END) AS bad", sf$validmu("__dm_mu"), sf$valideta("__dm_eta")),
      "MIN(__dm_mu) AS mumin", "MAX(__dm_mu) AS mumax",
      "SUM(CASE WHEN __dm_w > 0 AND (__dm_var IS NULL OR isnan(__dm_var)) THEN 1 ELSE 0 END) AS navar",
      "SUM(CASE WHEN __dm_w > 0 AND __dm_var = 0 THEN 1 ELSE 0 END) AS zvar",
      "SUM(CASE WHEN __dm_w > 0 AND (__dm_me IS NULL OR isnan(__dm_me)) THEN 1 ELSE 0 END) AS name",
      "SUM(CASE WHEN __dm_w > 0 AND __dm_me <> 0 THEN 1 ELSE 0 END) AS ngood")
}

dmglm_num_row <- function(r) {
    lapply(r, function(x) {
        x <- as.numeric(x)
        if (!length(x)) NA_real_ else x
    })
}

## glm.fit()'s checks before a weighted least squares step.
dmglm_wls_checks <- function(e) {
    if (isTRUE(e$navar > 0)) stop("NAs in V(mu)", call. = FALSE)
    if (isTRUE(e$zvar > 0)) stop("0s in V(mu)", call. = FALSE)
    if (isTRUE(e$name > 0)) stop("NAs in d(mu)/d(eta)", call. = FALSE)
    isTRUE(e$ngood > 0)
}

## ---- patterns ---------------------------------------------------------------

## The compressed fit, or NULL when compression does not pay.
dmglm_remote_patterns <- function(D, cx, start, control, tol, intercept, mode, family) {
    if (isFALSE(mode)) return(NULL)
    k <- cx$k
    sf <- cx$sf
    keys <- unique(c(dm_compress_keys(D), if (cx$cols$offset) "__dm_off"))
    n <- D$n
    if (n < 2) return(NULL)
    if (identical(mode, "auto") && length(keys) && !dm_compress_worth(D, keys, n)) return(NULL)
    cap <- if (isTRUE(mode)) n else min(n / 2, getOption("design.matrix.remote_compress_max", 1e6))
    kq <- if (length(keys)) dm_qi(D, keys) else character()
    ym <- dm_sql_double(cx$ck$ym)
    Y <- "__dm_y"
    W <- "__dm_w"
    m <- cx$m
    inner <- sprintf("SELECT * FROM (%s) AS s WHERE %s", D$src$select, D$filter)
    row_start <- is.null(start)
    if (row_start) {
        eta0 <- sf$linkfun(sf$mustart(Y, W, "__dm_n", cx$cols$cbind))
        inner <- dmglm_sql_work(sprintf("SELECT *, %s AS __dm_eta FROM (%s) AS t1", eta0, inner), cx)
    }
    aggs <- c("COUNT(*) AS __dm_cnt", sprintf("SUM(%s) AS __dm_sw", W), sprintf("SUM(%s * %s) AS __dm_swy", W, Y),
              sprintf("SUM(%s * (%s - %s) * (%s - %s)) AS __dm_syy", W, Y, ym, Y, ym),
              switch(k$dev,
                     binomial = sprintf(paste("SUM(%s * ((CASE WHEN %s <> 0 THEN %s * ln(%s) ELSE 0 END) +",
                                              "(CASE WHEN (1 - %s) <> 0 THEN (1 - %s) * ln(1 - %s) ELSE 0 END))) AS __dm_sat"),
                                        W, Y, Y, Y, Y, Y, Y),
                     poisson = sprintf("SUM(CASE WHEN %s > 0 THEN %s * %s * ln(%s) ELSE 0 END) AS __dm_sat", Y, W, Y, Y),
                     gamma = , qgamma = sprintf("SUM(%s * ln(CASE WHEN %s = 0 THEN 1 ELSE %s END)) AS __dm_sat", W, Y, Y),
                     invgauss = sprintf("SUM(%s / %s) AS __dm_sat", W, Y),
                     "SUM(CAST(0 AS DOUBLE)) AS __dm_sat"))
    if (k$aic == "binomial") {
        kk <- sprintf("round(%s * %s)", m, Y)
        mm <- sprintf("round(%s)", m)
        co <- sprintf("(CASE WHEN %s > 0 THEN %s / %s ELSE 0 END)", m, W, m)
        aggs <- c(aggs, sprintf("SUM(%s * (lgamma(%s + 1) - lgamma(%s + 1) - lgamma(%s - %s + 1))) AS __dm_a", co, mm, kk, mm, kk),
                  sprintf("SUM(%s * %s) AS __dm_k", co, kk), sprintf("SUM(%s * (%s - %s)) AS __dm_mk", co, mm, kk))
    }
    if (k$aic == "poisson") aggs <- c(aggs, sprintf("SUM(%s * lgamma(%s + 1)) AS __dm_l", W, Y))
    if (row_start) {
        aggs <- c(aggs, "SUM(__dm_wk) AS __dm_w1", "SUM(__dm_wk * __dm_z) AS __dm_wz1",
                  ## for the dispersion when only the first iteration runs
                  sprintf("SUM(__dm_wk * %s) AS __dm_w1y", Y), sprintf("SUM(__dm_wk * %s * %s) AS __dm_w1yy", Y, Y),
                  "SUM(CASE WHEN __dm_wk = 0 THEN 1 ELSE 0 END) AS __dm_w1zero",
                  sub(" AS dev$", " AS __dm_dev0", dmglm_sql_eval_items(cx)[1L]),
                  sprintf("SUM(CASE WHEN %s AND %s THEN 0 ELSE 1 END) AS __dm_bad0", sf$validmu("__dm_mu"), sf$valideta("__dm_eta")),
                  "SUM(CASE WHEN __dm_w > 0 AND (__dm_var IS NULL OR isnan(__dm_var)) THEN 1 ELSE 0 END) AS __dm_navar",
                  "SUM(CASE WHEN __dm_w > 0 AND __dm_var = 0 THEN 1 ELSE 0 END) AS __dm_zvar",
                  "SUM(CASE WHEN __dm_w > 0 AND (__dm_me IS NULL OR isnan(__dm_me)) THEN 1 ELSE 0 END) AS __dm_name",
                  "SUM(CASE WHEN __dm_w > 0 AND __dm_me <> 0 THEN 1 ELSE 0 END) AS __dm_ngood")
    }
    sql <- sprintf("SELECT %s FROM (%s) AS b%s%s", paste(c(kq, aggs), collapse = ", "), inner,
                   if (length(keys)) " GROUP BY ALL" else "",
                   if (length(keys)) sprintf(" LIMIT %.0f", cap + 1) else "")
    res <- dm_query(D, sql)
    G <- nrow(res)
    if (G > cap || (G == n && length(keys))) return(NULL)
    num <- function(x) {
        x <- as.double(x)
        x[is.na(x)] <- 0
        x
    }
    if (row_start) {
        tot <- function(nm) sum(num(res[[nm]]))
        if (tot("__dm_bad0") > 0) stop("cannot find valid starting values: please specify some", call. = FALSE)
    }
    pd <- dm_pattern_design(D, res[dm_compress_keys(D)])
    kp <- pd$keep
    xs <- pd$design
    W <- num(res$`__dm_sw`)[kp]
    ybar <- ifelse(W > 0, num(res$`__dm_swy`)[kp] / W, 0)
    off <- if (cx$cols$offset) num(res$`__dm_off`)[kp] else numeric(length(kp))
    sat <- num(res$`__dm_sat`)[kp]
    xlx <- function(x) ifelse(x > 0, x * log(x), 0)
    Cg <- switch(k$dev,
                 gaussian = num(res$`__dm_syy`)[kp] - W * (ybar - cx$ck$ym)^2,
                 binomial = 2 * sat - 2 * W * (xlx(ybar) + xlx(1 - ybar)),
                 poisson = 2 * (sat - W * xlx(ybar)),
                 gamma = , qgamma = -2 * sat + 2 * W * log(ybar),
                 invgauss = sat - W / ybar)
    Cg[W == 0] <- 0
    C <- sum(Cg)
    within <- pmax(num(res$`__dm_syy`)[kp] - W * (ybar - cx$ck$ym)^2, 0)
    if (!all(is.finite(c(W, ybar, C)))) return(NULL)
    lin <- function(b) drop(xs %*% b)
    if (row_start) {
        e0 <- list(navar = sum(num(res$`__dm_navar`)), zvar = sum(num(res$`__dm_zvar`)),
                   name = sum(num(res$`__dm_name`)), ngood = sum(num(res$`__dm_ngood`)))
        good <- dmglm_wls_checks(e0)
        W1 <- num(res$`__dm_w1`)[kp]
        first <- if (good) list(w2 = W1, z = ifelse(W1 > 0, num(res$`__dm_wz1`)[kp] / W1, 0)) else list(none = TRUE)
        B <- dmglm_units(lin, dmglm_wls(xs, tol), ybar, W, off, C, family, first = first)
        devold <- sum(num(res$`__dm_dev0`))
        coefold <- NULL
    } else {
        eta <- off + lin(start)
        B <- dmglm_units(lin, dmglm_wls(xs, tol), ybar, W, off, C, family, eta = eta)
        valideta <- family$valideta %||% function(eta) TRUE
        validmu <- family$validmu %||% function(mu) TRUE
        if (!(validmu(B$state$mu) && valideta(eta)))
            stop("cannot find valid starting values: please specify some", call. = FALSE)
        devold <- sum(family$dev.resids(ybar, B$state$mu, W)) + C
        coefold <- start
    }
    r <- dmglm_loop(B, family, control, cx$ck$n, devold, coefold)
    st <- B$state
    mu <- st$mu
    ## null deviance
    nulldev <- if (cx$cols$offset && intercept) {
        dmglm_null_refit(ybar, W, off, C, family, control, cx$ck$n, mu)
    } else {
        wtdmu <- if (intercept) cx$ck$swy / cx$ck$sw else family$linkinv(off)
        sum(family$dev.resids(ybar, wtdmu, W)) + C
    }
    ## AIC (without 2 * rank)
    sums <- switch(k$aic,
        binomial = {
            a <- num(res$`__dm_a`)[kp]
            kk <- num(res$`__dm_k`)[kp]
            mk <- num(res$`__dm_mk`)[kp]
            sum(a + ifelse(kk > 0, kk * log(mu), 0) + ifelse(mk > 0, mk * log(1 - mu), 0))
        },
        poisson = {
            swy <- W * ybar
            sum(ifelse(swy > 0, swy * log(mu), 0) - W * mu - num(res$`__dm_l`)[kp])
        },
        Gamma = {
            disp <- r$dev / cx$ck$sw
            a <- 1 / disp
            sum(W * (-lgamma(a) - a * log(mu * disp)) + (a - 1) * sat - W * ybar / (mu * disp))
        },
        NA_real_)
    aic <- dmglm_aic(k, cx$ck, r$dev, sums)
    ## dispersion
    pear <- if (identical(st$last, "rows")) {
        ## only the first iteration ran: its working weights vary within
        ## patterns, sum_i w_i (y_i - mu)^2 from their sums
        me <- family$mu.eta(st$eta)
        W1 <- num(res$`__dm_w1`)[kp]
        s <- (num(res$`__dm_w1yy`)[kp] - 2 * mu * num(res$`__dm_w1y`)[kp] + mu^2 * W1) / me^2
        list(pearson = sum(s[W1 > 0]), zero = sum(num(res$`__dm_w1zero`)) > 0)
    } else {
        me <- family$mu.eta(st$eta)
        w2 <- st$w2
        pos <- w2 > 0
        ratio <- ifelse(W > 0, w2 / W, 0)
        s <- ratio * (within + W * (ybar - mu)^2) / me^2
        list(pearson = sum(s[pos]), zero = cx$ck$nzero > 0 || any(w2[W > 0] == 0))
    }
    c(r, list(nulldev = nulldev, aic = aic, pearson = pear$pearson, zero = pear$zero, compressed = TRUE,
              patterns = length(kp), method = "patterns"))
}

## The AIC of glm.fit() (without 2 * rank) from sums over rows: `ll` is the
## sum of weighted log densities for binomial, poisson and Gamma.
dmglm_aic <- function(k, ck, dev, ll) {
    switch(k$aic,
           gaussian = {
               nobs <- ck$n
               if (ck$nzero > 0) Inf else nobs * (log(dev / nobs * 2 * pi) + 1) + 2 - ck$slw
           },
           binomial = -2 * ll,
           poisson = if (ck$nonint > 0) Inf else -2 * ll,
           Gamma = -2 * ll + 2,
           inverse.gaussian = ck$sw * (1 + log(dev / ck$sw * 2 * pi)) + 3 * ck$ly + 2,
           NA_real_)
}

## ---- iterations in the database or in batches -------------------------------

## The fit from per-iteration moments (engine "sql" or "chunked").
dmglm_remote_iterate <- function(D, cx, start, control, tol, intercept, family, engine) {
    ctr <- if (isTRUE(getOption("design.matrix.remote_center", TRUE)) && D$p > 0L) dm_center(D) else NULL
    tol_sel <- max(tol, 1e-07)
    if (engine == "chunked") {
        ## batches read a temporary table of the rows used, once
        tmp <- dm_tmp_name("dm_glm")
        dm_execute(D, sprintf("CREATE TEMPORARY TABLE %s AS SELECT * FROM (%s) AS s WHERE %s",
                              dm_qi(D, tmp), D$src$select, D$filter))
        on.exit(dm_execute(D, sprintf("DROP TABLE IF EXISTS %s", dm_qi(D, tmp))), add = TRUE)
        D$src$select <- sprintf("SELECT * FROM %s", dm_qi(D, tmp))
        D$src$base_table <- dm_qi(D, tmp)
        B <- dmglm_chunk_backend(D, cx, ctr, tol_sel, start, family)
    } else {
        B <- dmglm_sql_backend(D, cx, ctr, tol_sel, start)
        on.exit(B$cleanup(), add = TRUE)
    }
    devold <- B$init()
    r <- dmglm_loop(B, family, control, cx$ck$n, devold, start)
    fin <- B$final(r$coef, r$dev, intercept)
    c(r, fin, list(compressed = FALSE, patterns = cx$ck$n,
                   method = if (engine == "sql") "SQL iterations" else "iterations in batches"))
}

## Solve the weighted least squares step from moments of the (centered)
## design: coefficients (0 for aliased columns), rank, pivot, R. NULL when
## centering cannot be mapped back.
dmglm_solve_moments <- function(G, h, ctr, tol_sel, p) {
    if (!all(is.finite(G)) || !all(is.finite(h)))
        return(list(coefficients = rep(NaN, p), rank = 0L, pivot = seq_len(p), R = matrix(0, p, p)))
    dg <- diag(G)
    if (!is.null(ctr)) {
        S <- backsolve(ctr$T, diag(ncol(ctr$T)))
        dg <- colSums(S * (unname(G) %*% S))
    }
    fit <- dmfit_from_moments(G, h, 0, 0, tol_sel, FALSE, 0, integer(), dg)
    if (!is.null(ctr)) {
        fit <- dm_uncenter(fit, ctr$T)
        if (is.null(fit)) return(NULL)
    }
    b <- fit$coefficients
    b[is.na(b)] <- 0
    list(coefficients = b, rank = fit$rank, pivot = fit$pivot, R = unname(fit$R))
}

## The backend of the SQL engine (see the top of the file).
dmglm_sql_backend <- function(D, cx, ctr, tol_sel, start) {
    st <- new.env(parent = emptyenv())
    st$views <- character()
    sf <- cx$sf
    p <- D$p
    rows <- sprintf("SELECT * FROM (%s) AS s WHERE %s", D$src$select, D$filter)
    D2 <- D
    if (!is.null(ctr)) {
        for (nm in names(ctr$shift))
            D2$vars[[nm]]$sql <- sprintf("(%s - %s)", D$vars[[nm]]$sql, dm_sql_double(ctr$shift[[nm]]))
    }
    types <- c(D$src$types, `__dm_wk` = "double", `__dm_z` = "double")
    on_work <- function(DD, work) {
        DD$src$select <- work
        DD$src$base_table <- NULL
        DD$src$types <- types
        DD
    }
    ## a select of the rows with __dm_xb = X B (columns of B: several
    ## coefficient vectors give __dm_xb, __dm_xb2, ...)
    xb_select <- function(B) {
        B <- as.matrix(B)
        outs <- c("__dm_xb", if (ncol(B) > 1L) paste0("__dm_xb", seq_len(ncol(B))[-1L]))
        v <- dm_sql_predict(D, B, NULL, outs, names(D$src$types), FALSE)
        st$views <- c(st$views, v)
        sprintf("SELECT * FROM %s", dm_qi(D, v))
    }
    eta_select <- function(beta) {
        sprintf("SELECT *, __dm_xb + %s AS __dm_eta FROM (%s) AS t1", cx$off, xb_select(beta))
    }
    eta0 <- sf$linkfun(sf$mustart("__dm_y", "__dm_w", "__dm_n", cx$cols$cbind))
    run_eval <- function(work) {
        e <- dmglm_num_row(dm_query(D, sprintf("SELECT %s FROM (%s) AS t", paste(dmglm_sql_eval_items(cx), collapse = ", "), work)))
        st$work <- work
        st$e <- e
        e
    }
    moments <- function(DD, work) {
        Dw <- on_work(DD, work)
        G <- dm_sql_gram(Dw, "__dm_wk")
        h <- dm_sql_xty_core(Dw, "__dm_zz", "CAST(s.\"__dm_z\" AS DOUBLE) AS __dm_zz", "__dm_wk")[, 1L]
        list(G = unname(G), h = h)
    }
    list(
        state = st,
        init = function() {
            base <- if (is.null(start)) sprintf("SELECT *, %s AS __dm_eta FROM (%s) AS t1", eta0, rows)
                    else eta_select(start)
            e <- run_eval(dmglm_sql_work(base, cx))
            if (e$bad > 0) stop("cannot find valid starting values: please specify some", call. = FALSE)
            st$beta <- start
            e$dev
        },
        wls = function() {
            if (!dmglm_wls_checks(st$e)) return(NULL)
            m <- moments(D2, st$work)
            fit <- dmglm_solve_moments(m$G, m$h, ctr, tol_sel, p)
            if (is.null(fit)) {
                m <- moments(D, st$work)
                fit <- dmglm_solve_moments(m$G, m$h, NULL, tol_sel, p)
            }
            st$wls_beta <- st$beta
            st$wls_done <- TRUE
            fit
        },
        eval = function(beta) {
            st$beta <- beta
            e <- run_eval(dmglm_sql_work(eta_select(beta), cx))
            list(dev = e$dev, valid = isTRUE(e$bad == 0))
        },
        mu_range = function() c(st$e$mumin, st$e$mumax),
        final = function(coef, dev, intercept) {
            beta <- coef
            beta[is.na(beta)] <- 0
            ## linear predictors of the last weighted least squares step
            ## (__dm_eta0) and of the fit (__dm_eta)
            old <- st$wls_beta
            if (is.null(old)) {
                base <- sprintf("SELECT *, __dm_xb + %s AS __dm_eta, %s AS __dm_eta0 FROM (%s) AS t1",
                                cx$off, eta0, xb_select(beta))
            } else {
                base <- sprintf("SELECT *, __dm_xb + %s AS __dm_eta, __dm_xb2 + %s AS __dm_eta0 FROM (%s) AS t1",
                                cx$off, cx$off, xb_select(cbind(beta, old)))
            }
            r <- dmglm_sql_final(D, cx, base, dev, intercept)
            if (cx$cols$offset && intercept) {
                r$nulldev <- dmglm_sql_null_refit(D, cx, base, rows)
            }
            r
        },
        cleanup = function() {
            for (v in rev(st$views)) try(dm_execute(D, sprintf("DROP VIEW IF EXISTS %s", dm_qi(D, v))), silent = TRUE)
        })
}

## Deviance-based summaries at the end of a SQL fit: Pearson sum for the
## dispersion, AIC, null deviance (without the refit).
dmglm_sql_final <- function(D, cx, base, dev, intercept) {
    sf <- cx$sf
    k <- cx$k
    l2 <- sprintf("SELECT *, %s AS __dm_mu, %s AS __dm_me, %s AS __dm_mu0, %s AS __dm_me0 FROM (%s) AS t2",
                  sf$linkinv("__dm_eta"), sf$mu.eta("__dm_eta"), sf$linkinv("__dm_eta0"), sf$mu.eta("__dm_eta0"), base)
    l3 <- sprintf("SELECT *, %s AS __dm_var0 FROM (%s) AS t3", sf$variance("__dm_mu0"), l2)
    good0 <- "__dm_w > 0 AND __dm_me0 <> 0"
    wtd <- if (intercept) dm_sql_double(cx$ck$swy / cx$ck$sw) else sf$linkinv(cx$off)
    items <- c(sprintf("SUM(CASE WHEN %s THEN __dm_w * __dm_me0 * __dm_me0 / __dm_var0 * ((__dm_y - __dm_mu) / __dm_me) * ((__dm_y - __dm_mu) / __dm_me) ELSE 0 END) AS pearson", good0),
               sprintf("SUM(CASE WHEN %s THEN 0 ELSE 1 END) AS zero", good0),
               sprintf("SUM(%s) AS nulldev", sf$dev("__dm_y", wtd, "__dm_w")),
               sprintf("SUM(%s) AS ll", dmglm_sql_ll(cx, dev)))
    r <- dmglm_num_row(dm_query(D, sprintf("SELECT %s FROM (%s) AS t", paste(items, collapse = ", "), l3)))
    list(pearson = r$pearson, zero = r$zero > 0, nulldev = r$nulldev, aic = dmglm_aic(k, cx$ck, dev, r$ll))
}

## Per-row weighted log density for the AIC, over columns __dm_y, __dm_mu,
## __dm_w, __dm_n.
dmglm_sql_ll <- function(cx, dev) {
    Y <- "__dm_y"
    W <- "__dm_w"
    switch(cx$k$aic,
           binomial = {
               m <- cx$m
               kk <- sprintf("round(%s * %s)", m, Y)
               mm <- sprintf("round(%s)", m)
               sprintf(paste("(CASE WHEN %s > 0 THEN %s / %s ELSE 0 END) * (lgamma(%s + 1) - lgamma(%s + 1) - lgamma(%s - %s + 1)",
                             "+ (CASE WHEN %s > 0 THEN %s * ln(__dm_mu) ELSE 0 END)",
                             "+ (CASE WHEN %s - %s > 0 THEN (%s - %s) * ln(1 - __dm_mu) ELSE 0 END))"),
                       m, W, m, mm, kk, mm, kk, kk, kk, mm, kk, mm, kk)
           },
           poisson = sprintf("%s * (CASE WHEN %s = 0 THEN -__dm_mu ELSE %s * ln(__dm_mu) - __dm_mu - lgamma(%s + 1) END)", W, Y, Y, Y),
           Gamma = {
               disp <- dev / cx$ck$sw
               a <- 1 / disp
               sprintf("%s * (%s - %s * ln(__dm_mu * %s) + %s * ln(%s) - %s / (__dm_mu * %s))", W,
                       dm_sql_double(-lgamma(a)), dm_sql_double(a), dm_sql_double(disp), dm_sql_double(a - 1), Y, Y,
                       dm_sql_double(disp))
           },
           "CAST(0 AS DOUBLE)")
}

## The null deviance with an offset and an intercept, refitted in SQL as
## glm() does: intercept only, from the fitted values of the model.
dmglm_sql_null_refit <- function(D, cx, base_fit, rows) {
    sf <- cx$sf
    st <- new.env(parent = emptyenv())
    run <- function(work) {
        e <- dmglm_num_row(dm_query(D, sprintf("SELECT %s, SUM(__dm_wk) AS swk, SUM(__dm_wk * __dm_z) AS swz FROM (%s) AS t",
                                               paste(dmglm_sql_eval_items(cx), collapse = ", "), work)))
        st$e <- e
        e
    }
    first <- dmglm_sql_work(sprintf("SELECT * REPLACE (%s AS __dm_eta) FROM (%s) AS t0",
                                    sf$linkfun(sf$linkinv("__dm_eta")), base_fit), cx)
    B <- list(
        state = st,
        wls = function() {
            if (!dmglm_wls_checks(st$e)) return(NULL)
            dmglm_wls_mean_sums(st$e$swk, st$e$swz)
        },
        eval = function(beta) {
            e <- run(dmglm_sql_work(sprintf("SELECT *, %s + %s AS __dm_eta FROM (%s) AS t1",
                                            dm_sql_double(beta[1L]), cx$off, rows), cx))
            list(dev = e$dev, valid = isTRUE(e$bad == 0))
        },
        mu_range = function() c(st$e$mumin, st$e$mumax))
    e <- run(first)
    if (e$bad > 0) stop("cannot find valid starting values: please specify some", call. = FALSE)
    r <- dmglm_loop(B, cx$family, cx$control, cx$ck$n, e$dev)
    if (!r$conv)
        warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?", call. = FALSE)
    r$dev
}

dmglm_wls_mean_sums <- function(s, sz) {
    if (!is.finite(s) || !is.finite(sz))
        return(list(coefficients = NaN, rank = 0L, pivot = 1L, R = matrix(0, 1L, 1L)))
    list(coefficients = if (s > 0) sz / s else 0, rank = if (s > 0) 1L else 0L, pivot = 1L,
         R = matrix(sqrt(s), 1L, 1L))
}

## The backend of the chunked engine: one pass over the batches per
## evaluation, which also accumulates X'WX and X'Wz at the new coefficients
## (the next step needs them unless the fit has converged).
dmglm_chunk_backend <- function(D, cx, ctr, tol_sel, start, family) {
    k <- cx$k
    p <- D$p
    st <- new.env(parent = emptyenv())
    D2 <- D
    if (!is.null(ctr)) D2$shift <- ctr$shift
    Tinv <- if (is.null(ctr)) NULL else backsolve(ctr$T, diag(p))
    cols <- unique(c(dm_remote_inputs(D), "__dm_y", "__dm_w", "__dm_n", if (cx$cols$offset) "__dm_off"))
    valideta <- family$valideta %||% function(eta) TRUE
    validmu <- family$validmu %||% function(mu) TRUE
    cb <- cx$cols$cbind
    ## call FUN(X, y, w, n, off) on each batch, X the design (DD) of its rows
    each <- function(DD, FUN) {
        dm_chunks(DD, cols, function(chunk) {
            cd <- dm_chunk_design(DD, chunk)
            kp <- cd$keep
            off <- if (cx$cols$offset) chunk[["__dm_off"]][kp] else 0
            FUN(cd$design, chunk[["__dm_y"]][kp], chunk[["__dm_w"]][kp], chunk[["__dm_n"]][kp], off)
        })
    }
    eta_of <- function(X, y, w, n, off, beta, Tm) {
        if (is.null(beta)) return(family$linkfun(dmglm_mustart_r(k, y, w, n, cb)))
        b <- if (is.null(Tm)) beta else drop(Tm %*% beta)
        drop(X %*% b) + off
    }
    ## deviance, validity, glm.fit()'s checks and the moments at beta
    pass <- function(DD, Tm, beta) {
        a <- list(dev = 0, valid = TRUE, mumin = Inf, mumax = -Inf, navar = 0, zvar = 0, name = 0, ngood = 0,
                  G = matrix(0, p, p), h = numeric(p))
        each(DD, function(X, y, w, n, off) {
            eta <- eta_of(X, y, w, n, off, beta, Tm)
            mu <- family$linkinv(eta)
            a$dev <<- a$dev + sum(family$dev.resids(y, mu, w))
            a$valid <<- a$valid && isTRUE(validmu(mu) && valideta(eta))
            a$mumin <<- min(a$mumin, mu)
            a$mumax <<- max(a$mumax, mu)
            varmu <- family$variance(mu)
            me <- family$mu.eta(eta)
            g0 <- w > 0
            a$navar <<- a$navar + sum(is.na(varmu[g0]))
            a$zvar <<- a$zvar + sum(varmu[g0] == 0, na.rm = TRUE)
            a$name <<- a$name + sum(is.na(me[g0]))
            good <- g0 & !is.na(me) & me != 0
            a$ngood <<- a$ngood + sum(good)
            w2 <- ifelse(good, w * me^2 / varmu, 0)
            z <- ifelse(good, (eta - off) + (y - mu) / me, 0)
            w2[!is.finite(w2)] <- 0
            z[!is.finite(z)] <- 0
            a$G <<- a$G + as.matrix(crossprod(X, w2 * X))
            a$h <<- a$h + drop(as.matrix(crossprod(X, w2 * z)))
        })
        a
    }
    list(
        state = st,
        init = function() {
            r <- pass(D2, Tinv, start)
            if (!r$valid) stop("cannot find valid starting values: please specify some", call. = FALSE)
            st$r <- r
            st$beta <- start
            r$dev
        },
        wls = function() {
            if (!dmglm_wls_checks(st$r)) return(NULL)
            fit <- dmglm_solve_moments(unname(st$r$G), st$r$h, ctr, tol_sel, p)
            if (is.null(fit)) {
                r <- pass(D, NULL, st$beta)
                fit <- dmglm_solve_moments(unname(r$G), r$h, NULL, tol_sel, p)
            }
            st$wls_beta <- st$beta
            fit
        },
        eval = function(beta) {
            r <- pass(D2, Tinv, beta)
            st$r <- r
            st$beta <- beta
            list(dev = r$dev, valid = r$valid)
        },
        mu_range = function() c(st$r$mumin, st$r$mumax),
        final = function(coef, dev, intercept) {
            beta <- coef
            beta[is.na(beta)] <- 0
            old <- st$wls_beta
            wtd <- if (intercept) cx$ck$swy / cx$ck$sw else NULL
            disp <- dev / cx$ck$sw
            a <- list(pearson = 0, zero = FALSE, nulldev = 0, ll = 0)
            each(D2, function(X, y, w, n, off) {
                eta <- eta_of(X, y, w, n, off, beta, Tinv)
                eta0 <- eta_of(X, y, w, n, off, old, Tinv)
                mu <- family$linkinv(eta)
                mu0 <- family$linkinv(eta0)
                me0 <- family$mu.eta(eta0)
                good0 <- w > 0 & me0 != 0
                w0 <- ifelse(good0, w * me0^2 / family$variance(mu0), 0)
                r <- (y - mu) / family$mu.eta(eta)
                a$pearson <<- a$pearson + sum((w0 * r^2)[good0])
                a$zero <<- a$zero || any(!good0)
                a$nulldev <<- a$nulldev + sum(family$dev.resids(y, if (is.null(wtd)) family$linkinv(off + 0 * y) else wtd, w))
                a$ll <<- a$ll + switch(k$aic,
                    binomial = {
                        m <- if (cx$ck$maxn > 1) n else w
                        sum(ifelse(m > 0, w / m, 0) * stats::dbinom(round(m * y), round(m), mu, log = TRUE))
                    },
                    poisson = suppressWarnings(sum(stats::dpois(y, mu, log = TRUE) * w)),
                    Gamma = sum(stats::dgamma(y, 1 / disp, scale = mu * disp, log = TRUE) * w),
                    0)
            })
            out <- list(pearson = a$pearson, zero = a$zero, nulldev = a$nulldev,
                        aic = dmglm_aic(k, cx$ck, dev, a$ll))
            if (cx$cols$offset && intercept) out$nulldev <- dmglm_chunk_null_refit(each, D2, Tinv, beta, eta_of, cx, family)
            out
        })
}

## The null deviance with an offset and an intercept, refitted in batches.
dmglm_chunk_null_refit <- function(each, D2, Tinv, beta, eta_of, cx, family) {
    st <- new.env(parent = emptyenv())
    valideta <- family$valideta %||% function(eta) TRUE
    validmu <- family$validmu %||% function(mu) TRUE
    pass <- function(etafun) {
        a <- list(dev = 0, valid = TRUE, mumin = Inf, mumax = -Inf, navar = 0, zvar = 0, name = 0, ngood = 0,
                  s = 0, sz = 0)
        each(D2, function(X, y, w, n, off) {
            eta <- etafun(X, y, w, n, off)
            mu <- family$linkinv(eta)
            a$dev <<- a$dev + sum(family$dev.resids(y, mu, w))
            a$valid <<- a$valid && isTRUE(validmu(mu) && valideta(eta))
            a$mumin <<- min(a$mumin, mu)
            a$mumax <<- max(a$mumax, mu)
            varmu <- family$variance(mu)
            me <- family$mu.eta(eta)
            g0 <- w > 0
            a$navar <<- a$navar + sum(is.na(varmu[g0]))
            a$zvar <<- a$zvar + sum(varmu[g0] == 0, na.rm = TRUE)
            a$name <<- a$name + sum(is.na(me[g0]))
            good <- g0 & !is.na(me) & me != 0
            a$ngood <<- a$ngood + sum(good)
            w2 <- ifelse(good, w * me^2 / varmu, 0)
            a$s <<- a$s + sum(w2[good])
            a$sz <<- a$sz + sum((w2 * ((eta - off) + (y - mu) / me))[good])
        })
        st$r <- a
        a
    }
    r <- pass(function(X, y, w, n, off) family$linkfun(family$linkinv(eta_of(X, y, w, n, off, beta, Tinv))))
    if (!r$valid) stop("cannot find valid starting values: please specify some", call. = FALSE)
    B <- list(
        state = st,
        wls = function() {
            if (!dmglm_wls_checks(st$r)) return(NULL)
            dmglm_wls_mean_sums(st$r$s, st$r$sz)
        },
        eval = function(b) {
            r <- pass(function(X, y, w, n, off) b[1L] + off + 0 * y)
            list(dev = r$dev, valid = r$valid)
        },
        mu_range = function() c(st$r$mumin, st$r$mumax))
    res <- dmglm_loop(B, family, cx$control, cx$ck$n, r$dev)
    if (!res$conv)
        warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?", call. = FALSE)
    res$dev
}

## ---- predictions and residuals ---------------------------------------------

## A view (or table) of a remote design with __dm_xb = X beta, and the
## final columns built from it: `items` (named SQL expressions over the
## columns, with __dm_eta = __dm_xb + offset).
dmglm_remote_view <- function(D, beta, items, off, drop = character()) {
    tmp <- dm_remote_predict(D, beta, column = "__dm_xb")
    con <- dm_con(D)
    ex <- c("__dm_xb", "__dm_eta", drop)
    sel <- sprintf("SELECT * EXCLUDE (%s), %s FROM (SELECT *, __dm_xb + %s AS __dm_eta FROM %s) AS t",
                   paste(DBI::dbQuoteIdentifier(con, intersect(ex, c(names(D$src$types), "__dm_xb", "__dm_eta"))),
                         collapse = ", "),
                   paste(sprintf("%s AS %s", items, DBI::dbQuoteIdentifier(con, names(items))), collapse = ", "),
                   off, DBI::dbQuoteIdentifier(con, unclass(tmp)[1L]))
    name <- dm_tmp_name("dm_pred")
    dm_execute(D, sprintf("CREATE TEMPORARY VIEW %s AS %s", dm_qi(D, name), sel))
    structure(name, class = "dm_remote_table", connection = con, holder = D$src$holder, base = tmp)
}

## predict(fit, newdata = <remote>, type = ).
dm_remote_glm_predict <- function(object, newdata, type, na.action, dots) {
    if (is.null(object$spec)) stop("This design cannot be applied to new data.", call. = FALSE)
    k <- dmglm_check_family(object$family)
    D <- do.call(design_matrix, c(list(object$spec, newdata, na.action = na.action), dots))
    if (object$rank < length(object$coefficients))
        warning("prediction from rank-deficient fit; consider predict(., rankdeficient=\"NA\")", call. = FALSE)
    ## offsets: offset() terms of the formula and the `offset` column
    con <- dm_con(D)
    pieces <- character()
    tt <- object$terms
    if (!is.null(tt)) {
        vars <- attr(tt, "variables")
        for (i in attr(tt, "offset")) {
            r <- dm_sql_num(vars[[i + 1L]][[2L]], D$src$types, con)
            if (!is.null(r$reason)) stop(sprintf("The offset has no SQL translation (%s).", r$reason), call. = FALSE)
            pieces <- c(pieces, r$sql)
        }
    }
    if (!is.null(object$offset_call)) {
        col <- dm_remote_column(D, as.character(object$offset_call), "offset")
        pieces <- c(pieces, sprintf("CAST(s.%s AS DOUBLE)", dm_qi(D, col)))
    }
    off <- "CAST(0 AS DOUBLE)"
    drop <- character()
    if (length(pieces)) {
        D <- dm_remote_extend(D, c(`__dm_off` = paste0("(", paste(pieces, collapse = " + "), ")")), "offset")
        off <- "__dm_off"
        drop <- "__dm_off"
    }
    dmglm_sql_macros(D)
    sf <- dmglm_sql_family(k)
    ## (greatest() skips NULL: keep the predictions of rows with missing
    ## values NULL)
    item <- if (type == "response") sprintf("CASE WHEN __dm_eta IS NULL THEN NULL ELSE %s END", sf$linkinv("__dm_eta"))
            else "__dm_eta"
    dmglm_remote_view(D, dmfit_beta(object), c(.fitted = item), off, drop)
}

## Predictions on the rows of a remote fit: a view with `.fitted`.
dm_remote_glm_fitted <- function(object, type = "response") {
    D <- object$remote_design
    dmglm_sql_macros(D)
    sf <- dmglm_sql_family(dmglm_check_family(object$family))
    off <- if (object$remote_cols$offset) "__dm_off" else "CAST(0 AS DOUBLE)"
    item <- if (type == "response") sf$linkinv("__dm_eta") else "__dm_eta"
    dmglm_remote_view(D, dmfit_beta(object), c(.fitted = item), off)
}

## Residuals of a remote fit: the remote design over the same rows with
## `.fitted` (the fitted means) and `.resid` columns, as residuals() of a
## remote design_lm() fit.
dm_remote_glm_residuals <- function(object, type) {
    D <- object$remote_design
    taken <- intersect(c(".fitted", ".resid"), names(D$src$types))
    if (length(taken)) stop(sprintf("The remote data already has a column `%s`.", taken[1L]), call. = FALSE)
    dmglm_sql_macros(D)
    sf <- dmglm_sql_family(dmglm_check_family(object$family))
    off <- if (object$remote_cols$offset) "__dm_off" else "CAST(0 AS DOUBLE)"
    mu <- sf$linkinv("__dm_eta")
    res <- switch(type,
                  response = sprintf("__dm_y - %s", mu),
                  working = sprintf("(__dm_y - %s) / %s", mu, sf$mu.eta("__dm_eta")),
                  pearson = sprintf("(__dm_y - %s) * sqrt(__dm_w) / sqrt(%s)", mu, sf$variance(mu)),
                  deviance = if (object$df.residual > 0) {
                      sprintf("(CASE WHEN __dm_y > %s THEN 1 ELSE -1 END) * sqrt(greatest(%s, 0))", mu,
                              sf$dev("__dm_y", mu, "__dm_w"))
                  } else "CAST(0 AS DOUBLE)")
    view <- dmglm_remote_view(D, dmfit_beta(object), c(.fitted = mu, .resid = res), off)
    out <- D
    out$src$select <- sprintf("SELECT * FROM %s", dm_qi(D, unclass(view)[1L]))
    out$src$base_table <- NULL
    out$src$view <- view
    proto <- DBI::dbGetQuery(dm_con(D), sprintf("SELECT * FROM (%s) AS s LIMIT 0", out$src$select))
    out$src$proto <- proto
    out$src$types <- vapply(proto, dm_col_type, "")
    out
}

## design_glm(), one implementation in memory and on remote data: each case
## of glm_cases.R fitted in memory (compressed or not) and on a DuckDB table
## (compressed or not) and on Parquet (full tier), compared
## with glm() on the same rows: coefficients (with the NA pattern of
## aliased columns), vcov, deviances, AIC, degrees of freedom, dispersion
## and the number of iterations. The other checks (warnings, methods,
## responses, start values) run in every mode as well.
source("helpers.R")
source("glm_cases.R", local = TRUE)
library(design.matrix)

d <- glm_data(n = 800)
dm <- d
dm$f <- factor(dm$f)
dm$g <- factor(dm$g)

## ---- modes ----
modes <- list(memory = list(name = "memory", data = list(dm)))
has_duckdb <- requireNamespace("duckdb", quietly = TRUE) && requireNamespace("DBI", quietly = TRUE)
if (has_duckdb) {
    con <- DBI::dbConnect(duckdb::duckdb())
    DBI::dbWriteTable(con, "d", d)
    pq <- tempfile(fileext = ".parquet")
    DBI::dbExecute(con, sprintf("COPY d TO '%s' (FORMAT PARQUET)", pq))
    ## (the remote methods of the primitives land on kata/2hba-duckdb: until
    ## then, the remote modes are skipped)
    probe <- try(suppressWarnings(design_glm(yb ~ x + g, binomial, con, table = "d", compress = FALSE)), silent = TRUE)
    if (!inherits(probe, "try-error")) {
        modes$sql <- list(name = "sql", data = list(con, table = "d", engine = "sql"))
    } else {
        message("Remote GLM modes skipped: ", conditionMessage(attr(probe, "condition")))
    }
}

## design_glm() in a mode; `...` holds the other arguments (weights and
## offset as column names, which work in every mode)
fit_in <- function(mode, fo, fam, ..., data = mode$data) {
    args <- list(...)
    if (mode$name == "memory") {
        for (a in intersect(names(args), c("weights", "offset")))
            if (is.character(args[[a]])) args[[a]] <- data[[1L]][[args[[a]]]]
        return(do.call(design_glm, c(list(fo, fam, data = data[[1L]]), args)))
    }
    do.call(design_glm, c(list(fo, fam), data, args))
}

## A row-wise result as a vector, in the order of the rows (remote tables
## are sorted by `id`).
as_vec <- function(x, col) {
    if (!inherits(x, "dm_remote_table")) return(unname(as.vector(x)))
    r <- as.data.frame(x)
    r[[col]][order(r$id)]
}

dd <- if (has_duckdb) DBI::dbReadTable(con, "d") else d
case_args <- function(cs) {
    a <- list()
    if (length(cs) >= 4 && !is.null(cs[[4]])) a$weights <- cs[[4]]
    if (length(cs) >= 5 && !is.null(cs[[5]])) a$offset <- cs[[5]]
    a
}

## ---- the cases ----
for (cs in glm_cases) {
    b <- suppressWarnings(glm_ref(cs, dd))
    for (mode in modes) for (cp in list("auto", FALSE)) {
        ## (remote fits without compression: full tier; the checks below
        ## fit remotely without compression in the default tier too)
        if (!full_tests && mode$name != "memory" && isFALSE(cp)) next
        a <- suppressWarnings(do.call(fit_in, c(list(mode, cs[[2]], cs[[3]]), case_args(cs), compress = cp)))
        glm_same(a, b, paste(cs[[1]], mode$name, "compress", cp))
    }
    if (length(modes) > 1L && "parquet" %in% remote_modes && cs[[1]] %in% c("logit x*f", "binomial cbind", "poisson offset()", "Gamma inverse")) {
        for (eng in "sql") {
            a <- suppressWarnings(do.call(fit_in, c(list(modes$sql, cs[[2]], cs[[3]]), case_args(cs),
                                                   data = list(list(pq, engine = eng)))))
            glm_same(a, b, paste(cs[[1]], "parquet", eng))
        }
    }
}

for (mode in modes) {
    info <- mode$name
    ## compression where rows share patterns (and offset values)
    fit <- fit_in(mode, yb ~ xd * f + g, binomial)
    expect_true(fit$compressed, info = info)
    expect_equal(fit$patterns, 48L, info = info)
    expect_false(fit_in(mode, yb ~ xd * f + g, binomial, compress = FALSE)$compressed, info = info)
    fo <- fit_in(mode, yp ~ g + offset(log(xd)), poisson)
    expect_true(fo$compressed, info = info)
    expect_equal(fo$patterns, 12L, info = info)
    glm_same(fo, glm(yp ~ g + offset(log(xd)), poisson, dd), paste("offset", info))

    ## warnings, as glm()
    expect_identical(glm_warnings(fit_in(mode, sep ~ x, binomial))$warnings,
                     glm_warnings(glm(sep ~ x, binomial, dd))$warnings, info = paste("separation", info))
    r2 <- glm_warnings(fit_in(mode, yb ~ x + f, binomial, control = list(maxit = 2)))
    expect_identical(r2$warnings, glm_warnings(glm(yb ~ x + f, binomial, dd, control = list(maxit = 2)))$warnings,
                     info = paste("maxit", info))
    expect_false(r2$value$converged, info = info)
    ## one iteration: the dispersion uses the working weights of the rows
    for (cp in list("auto", FALSE)) {
        a1 <- suppressWarnings(fit_in(mode, yg ~ xd + f, Gamma, compress = cp, control = list(maxit = 1)))
        glm_same(a1, suppressWarnings(glm(yg ~ xd + f, Gamma, dd, control = list(maxit = 1))),
                 paste("maxit 1", info, cp))
    }
    expect_warning(fit_in(mode, prop ~ xd, binomial, weights = "w"), "non-integer", info = info)
    expect_warning(summary(fit_in(mode, yg ~ x, Gamma, weights = "w")), "zero weight", info = info)
    expect_error(fit_in(mode, yp ~ x + g, poisson, na.action = na.fail), "missing values", info = info)

    ## methods
    b <- glm(yp ~ x + g + offset(log(e)), poisson, dd)
    a <- fit_in(mode, yp ~ x + g + offset(log(e)), poisson)
    for (ty in c("deviance", "pearson", "working", "response"))
        expect_equal(as_vec(residuals(a, ty), ".resid"), unname(residuals(b, ty)), info = paste(info, ty))
    expect_equal(as_vec(fitted(a), ".fitted"), unname(fitted(b)), info = info)
    for (ty in c("link", "response")) {
        expect_equal(as_vec(predict(a, type = ty), ".fitted"), unname(predict(b, type = ty)), info = paste(info, ty))
        expect_equal(as_vec(predict(a, newdata = dd[1:20, ], type = ty), ".fitted"),
                     unname(predict(b, dd[1:20, ], type = ty)), info = paste(info, ty, "data frame"))
        if (has_duckdb)
            expect_equal(as_vec(predict(a, newdata = con, table = "d", type = ty), ".fitted"),
                         unname(predict(b, dd, type = ty)), info = paste(info, ty, "remote newdata"))
    }
    expect_equal(logLik(a), logLik(b), info = info)
    expect_equal(BIC(a), BIC(b), info = info)
    expect_identical(family(a)$family, "poisson")
    expect_equal(summary(a)$coefficients, summary(b)$coefficients, info = info)
    expect_true(any(grepl("Fisher Scoring", capture.output(print(summary(a))))), info = info)
    expect_true(any(grepl("Residual Deviance", capture.output(print(a)))), info = info)

    ## a factor response
    expect_equal(coef(fit_in(mode, fy ~ x, binomial)), coef(glm(factor(fy) ~ x, binomial, dd)), info = info)

    ## start values
    s0 <- coef(glm(yp ~ xd + g, poisson, dd)) + 0.1
    for (cp in list("auto", FALSE)) {
        a <- fit_in(mode, yp ~ xd + g, poisson, start = s0, compress = cp)
        b <- glm(yp ~ xd + g, poisson, dd, start = s0)
        expect_equal(a$iter, b$iter, info = paste("start", info, cp))
        expect_equal(coef(a), coef(b), info = paste("start", info, cp))
    }
}

## regression: nested expressions (probit: DuckDB expands the qnorm and
## pnorm macros inline) once made queries so long that DuckDB hung; the
## linear predictor, fitted means and working weights are computed once
## (dmc_freeze()), which bounds the queries (full tier)
if (length(modes) > 1L && full_tests) {
    qlen <- integer()
    withCallingHandlers({
        op <- options(design.matrix.remote_debug = TRUE)
        fp <- design_glm(yb ~ xd + f, binomial("probit"), con, table = "d", compress = FALSE)
        options(op)
    }, message = function(m) {
        qlen <<- c(qlen, nchar(conditionMessage(m)))
        invokeRestart("muffleMessage")
    })
    expect_true(max(qlen) < 1e5)
    glm_same(fp, glm(yb ~ xd + f, binomial("probit"), dd), "probit query size")
}

## initialize's checks
expect_error(design_glm(yg - 3 ~ x, Gamma, dm), "non-positive values")
expect_error(design_glm(yp - 1 ~ x, poisson, dm), "negative values")
if (length(modes) > 1L)
    expect_error(design_glm(yp ~ x, Gamma, con, sql = "SELECT x, yp - 3 AS yp FROM d"), "non-positive values")

## na.exclude pads in memory
fx <- design_glm(yp ~ x + g, poisson, dm, na.action = na.exclude)
bx <- glm(yp ~ x + g, poisson, dm, na.action = na.exclude)
expect_equal(residuals(fx), residuals(bx))
expect_equal(fitted(fx), fitted(bx))

## a design with y: in memory a vector or a two-column matrix, remotely
## column names; a missing response is an error in both (no rows are
## dropped from a design), as for design_lm()
dc <- dd[complete.cases(dd[c("x", "g", "yb")]), ]
D <- design_matrix(~ x + g, dc)
expect_equal(coef(design_glm(D, binomial, y = dc$yb)), coef(glm(yb ~ x + g, binomial, dc)))
expect_equal(coef(design_glm(D, binomial, y = cbind(dc$s, dc$fl))), coef(glm(cbind(s, fl) ~ x + g, binomial, dc)))
if (length(modes) > 1L) {
    expect_error(design_glm(design_matrix(~ x + g, con, table = "d"), binomial, y = "yb"), "NA/NaN/Inf in 'y'")
    expect_error(design_glm(design_matrix(~ x + g, dd), binomial, y = dd$yb[!is.na(dd$x)]), "NA/NaN/Inf in 'y'")
    Dr <- design_matrix(~ x + g, con, sql = "SELECT * FROM d WHERE yb IS NOT NULL")
    expect_equal(coef(design_glm(Dr, binomial, y = "yb")), coef(glm(yb ~ x + g, binomial, dd)))
    Da <- design_matrix(~ x + g, con, table = "d")
    expect_equal(coef(design_glm(Da, binomial, y = c("s", "fl"))), coef(glm(cbind(s, fl) ~ x + g, binomial, dd)))
    expect_equal(coef(design_glm(Da, poisson, y = "yp", offset = "le")),
                 coef(glm(yp ~ x + g, poisson, dd, offset = le)))
}

## what design_glm() does not support (both modes)
fam_other <- poisson()
fam_other$family <- "other"
expect_error(design_glm(yp ~ x, fam_other, dm), "supports the families")
expect_error(design_glm(yp ~ x, quasi(variance = "mu", link = power(1/3)), dm), "supports the families")
if (has_duckdb) {
    expect_error(design_glm(yp ~ x, fam_other, con, table = "d"), "supports the families")
    DBI::dbDisconnect(con)
}

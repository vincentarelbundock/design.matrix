## Parity of design_lm() in memory and on remote data: one set of cases,
## each fitted by lm(), by design_lm() on a data frame, and by design_lm()
## on a DuckDB table and a Parquet file (full tier), compared on every
## method.
if (!requireNamespace("duckdb", quietly = TRUE) || !requireNamespace("DBI", quietly = TRUE))
    exit_file("duckdb is not installed")
source("helpers.R")
library(design.matrix)

con <- test_con()
same <- function(a, b, tol = 1e-7) isTRUE(all.equal(unname(a), unname(b), tolerance = tol))

set.seed(7)
n <- 600
d <- data.frame(
    id = seq_len(n),
    x = rnorm(n),
    xd = sample(1:4, n, TRUE),
    f = sample(c("b", "a", "c"), n, TRUE),
    g = sample(c("u", "v"), n, TRUE),
    z = runif(n, 1, 2),
    w = runif(n, 0.5, 2))
d$y <- d$x + d$xd + as.integer(factor(d$f)) + d$z + rnorm(n)
d$f2 <- d$f
d$x[5] <- NA
d$y[7] <- NA
d$f[9] <- NA
d$w[11] <- NA
d$z[13] <- NaN
DBI::dbWriteTable(con, "d", d)
d <- DBI::dbReadTable(con, "d")
d <- d[order(d$id), ]
rownames(d) <- NULL
pq <- tempfile(fileext = ".parquet")
DBI::dbExecute(con, sprintf("COPY d TO '%s' (FORMAT PARQUET)", pq))

## functions with no SQL translation: evaluated in R, batch by batch
myfun <- function(v) v^2 + 1
my_transform <- function(v) 2 * v - 1

## each case: a formula and the other arguments, as unevaluated R code
cases <- list(
    plain = list(fo = y ~ x + f),
    weighted = list(fo = y ~ x * f, args = alist(weights = w)),
    compressible = list(fo = y ~ xd + f + g),
    compressible_weighted = list(fo = y ~ xd * g + f, args = alist(weights = w)),
    subset = list(fo = y ~ x + f, args = alist(subset = x > -1 & f != "c")),
    offset_term = list(fo = y ~ x + offset(z) + f),
    offset_arg = list(fo = y ~ xd + f, args = alist(offset = log(z))),
    aliased = list(fo = y ~ x + f + f2),
    exclude = list(fo = y ~ x + f, args = alist(na.action = na.exclude)),
    everything = list(fo = y ~ log(z) + I(x^2) + f, args = alist(weights = w, offset = z, subset = g == "u")),
    r_weights = list(fo = y ~ x + f, args = alist(weights = myfun(x))),
    r_subset = list(fo = y ~ x + f, args = alist(subset = myfun(z) > 2.5)),
    r_offset = list(fo = y ~ xd + f, args = alist(offset = log(myfun(z)))),
    r_offset_term = list(fo = y ~ x + offset(log(myfun(z))) + f),
    r_response = list(fo = my_transform(y) ~ xd + f + g),
    r_everything = list(fo = my_transform(y) ~ xd * g + offset(myfun(z)),
                        args = alist(weights = myfun(xd), subset = myfun(x) < 3, na.action = na.exclude)))

fit_call <- function(fun, case, data, extra = list()) {
    as.call(c(list(fun, case$fo, data = data), case$args, extra))
}

check_methods <- function(f, l, info, tol = 1e-7, formula = TRUE) {
    expect_identical(is.na(coef(f)), is.na(coef(l)), info = paste(info, "aliased"))
    expect_true(same(coef(f), coef(l), tol), info = paste(info, "coef"))
    expect_true(same(vcov(f), vcov(l), tol), info = paste(info, "vcov"))
    expect_true(same(sigma(f), sigma(l), tol), info = paste(info, "sigma"))
    expect_identical(df.residual(f), df.residual(l), info = paste(info, "df"))
    expect_identical(as.numeric(nobs(f)), as.numeric(nobs(l)), info = paste(info, "nobs"))
    expect_true(same(deviance(f), deviance(l), tol), info = paste(info, "deviance"))
    expect_true(same(logLik(f), logLik(l), tol), info = paste(info, "logLik"))
    expect_true(same(AIC(f), AIC(l), tol) && same(BIC(f), BIC(l), tol), info = paste(info, "AIC/BIC"))
    sf <- summary(f)
    sl <- summary(l)
    expect_true(same(sf$r.squared, sl$r.squared, tol) && same(sf$adj.r.squared, sl$adj.r.squared, tol),
                info = paste(info, "R2"))
    expect_true(same(sf$fstatistic, sl$fstatistic, tol), info = paste(info, "F"))
    expect_true(same(coef(sf), coef(sl), tol), info = paste(info, "coefficient table"))
    if (formula) expect_equal(formula(f), formula(l), info = paste(info, "formula"))
}

## a remote table of predictions, in the order of d (NA for missing rows)
by_id <- function(tab, col) {
    p <- as.data.frame(tab)
    p[[col]][match(d$id, p$id)]
}

for (nm in names(cases)) {
    case <- cases[[nm]]
    l <- eval(fit_call(quote(lm), case, quote(d)))
    m <- eval(fit_call(quote(design_lm), case, quote(d)))
    check_methods(m, l, paste(nm, "memory"))
    sources <- list(
        sql = list(data = quote(con), extra = list(table = "d")),
        parquet = list(data = quote(pq), extra = list()))[remote_modes]
    for (src in names(sources)) {
        info <- paste(nm, src)
        s <- sources[[src]]
        r <- suppressMessages(eval(fit_call(quote(design_lm), case, s$data, s$extra)))
        check_methods(r, l, info)
        ## fitted values and residuals, as a table (NA padding with na.exclude)
        ft <- fitted(r)
        fl <- fitted(l)
        rl <- residuals(l)
        if (!is.null(case$args$na.action)) {
            ## padded with NA over the rows of the subset
            ids <- d$id[as.integer(names(fl))]
            p <- as.data.frame(ft)
            expect_identical(sort(p$id), sort(ids), info = paste(info, "fitted rows"))
            expect_true(same(p$.fitted[match(ids, p$id)], fl), info = paste(info, "fitted"))
            expect_true(same(p$.resid[match(ids, p$id)], rl), info = paste(info, "residuals"))
        } else {
            ids <- as.integer(names(fl))
            p <- as.data.frame(ft)
            expect_identical(sort(p$id), sort(d$id[ids]), info = paste(info, "fitted rows"))
            expect_true(same(p$.fitted[match(d$id[ids], p$id)], fl), info = paste(info, "fitted"))
            expect_true(same(p$.resid[match(d$id[ids], p$id)], rl), info = paste(info, "residuals"))
        }
        ## the residual design: X'e = 0 (up to rounding), and its rows
        R <- residuals(r)
        expect_equal(nrow(R), nobs(l) + if (!is.null(l$weights)) sum(l$weights == 0) else 0, info = paste(info, "residual design"))
        ## predictions on new data, with standard errors and intervals
        p_l <- tryCatch(suppressWarnings(predict(l, newdata = d, se.fit = TRUE, interval = "confidence", level = 0.9)),
                        error = function(e) e)
        if (inherits(p_l, "error")) {
            ## new levels: the same error
            msg <- tryCatch(predict(r, newdata = con, table = "d"), error = conditionMessage)
            expect_identical(msg, conditionMessage(p_l), info = paste(info, "predict error"))
            next
        }
        p_r <- suppressWarnings(predict(r, newdata = con, table = "d", se.fit = TRUE, interval = "confidence", level = 0.9))
        expect_true(same(by_id(p_r, ".fitted"), p_l$fit[, "fit"]), info = paste(info, "predict"))
        expect_true(same(by_id(p_r, ".se.fit"), p_l$se.fit), info = paste(info, "se.fit"))
        expect_true(same(by_id(p_r, ".lwr"), p_l$fit[, "lwr"]) && same(by_id(p_r, ".upr"), p_l$fit[, "upr"]),
                    info = paste(info, "interval"))
        expect_equal(attr(p_r, "df"), p_l$df, info = paste(info, "df of predict"))
        pp_l <- suppressWarnings(predict(l, newdata = d, interval = "prediction"))
        pp_r <- suppressWarnings(predict(r, newdata = con, table = "d", interval = "prediction"))
        expect_true(same(by_id(pp_r, ".lwr"), pp_l[, "lwr"]), info = paste(info, "prediction interval"))
        ## the fit only, in SQL
        expect_true(same(by_id(predict(r, newdata = con, table = "d"), ".fitted"),
                         suppressWarnings(predict(l, newdata = d))), info = paste(info, "predict fit"))
    }
}

## warnings of predict(), as predict.lm()
lw <- lm(y ~ x + f, d, weights = w)
rw <- design_lm(y ~ x + f, con, table = "d", weights = w)
wl <- tryCatch(predict(lw, newdata = d, interval = "prediction"), warning = conditionMessage)
wr <- tryCatch(predict(rw, newdata = con, table = "d", interval = "prediction"), warning = conditionMessage)
expect_identical(wr, wl)
wl <- tryCatch(predict(lw, interval = "prediction"), warning = conditionMessage)
wr <- tryCatch(predict(rw, interval = "prediction"), warning = conditionMessage)
expect_identical(wr, wl)
## training rows with prediction intervals: variance sigma^2 / w
pl <- suppressWarnings(predict(lw, interval = "prediction"))
pr <- suppressWarnings(predict(rw, interval = "prediction"))
ids <- as.integer(rownames(pl))
pr <- as.data.frame(pr)
expect_true(same(pr$.lwr[match(d$id[ids], pr$id)], pl[, "lwr"]))

## rank-deficient fits: non-estimable predictions flagged as in predict.lm()
la <- lm(y ~ x + f + f2, d)
ra <- design_lm(y ~ x + f + f2, con, table = "d")
nd <- data.frame(id = 1:3, x = 0, f = c("a", "a", "c"), f2 = c("a", "b", "c"), y = 0)
DBI::dbWriteTable(con, "nd", nd)
for (rd in c("warnif", "NA", "non-estim", "NAwarn")) {
    wa <- character()
    a <- withCallingHandlers(predict(la, newdata = nd, rankdeficient = rd),
                             warning = function(w) { wa <<- c(wa, conditionMessage(w)); invokeRestart("muffleWarning") })
    wb <- character()
    b <- withCallingHandlers(as.data.frame(predict(ra, newdata = con, table = "nd", rankdeficient = rd)),
                             warning = function(w) { wb <<- c(wb, conditionMessage(w)); invokeRestart("muffleWarning") })
    b <- b[order(b$id), ]
    expect_true(same(b$.fitted, unclass(a)[1:3]), info = rd)
    expect_identical(wb, wa, info = paste(rd, "warnings"))
    if (rd %in% c("warnif", "non-estim")) expect_identical(which(b$.nonestim), as.integer(attr(a, "non-estim")), info = rd)
}

## na.fail and na.pass, as in memory
expect_error(design_lm(y ~ x + f, con, table = "d", na.action = na.fail), "missing values in object")
expect_identical(tryCatch(design_lm(y ~ x + f, con, table = "d", na.action = na.pass), error = conditionMessage),
                 tryCatch(design_lm(y ~ x + f, d, na.action = na.pass), error = conditionMessage))
dk <- d[!is.na(d$y), ]
DBI::dbWriteTable(con, "dk", dk)
expect_identical(tryCatch(design_lm(y ~ x + f, con, table = "dk", na.action = na.pass), error = conditionMessage),
                 tryCatch(design_lm(y ~ x + f, dk, na.action = na.pass), error = conditionMessage))

## a remote design with y, weights and offset (columns or expressions)
D <- design_matrix(~ x + f, con, table = "d", na.action = na.omit)
ok <- stats::complete.cases(d[c("x", "f", "y", "w", "z")])
l <- lm(y ~ x + f, d[ok, ], weights = w, offset = log(z))
Dk <- design_matrix(~ x + f, con, sql = "SELECT * FROM d WHERE y IS NOT NULL AND w IS NOT NULL AND NOT isnan(z)",
                    na.action = na.omit)
## (a fit on a design has no formula, in memory as remotely)
fr <- design_lm(Dk, y = y, weights = "w", offset = log(z))
check_methods(fr, l, "remote design", formula = FALSE)
Dm <- design_matrix(~ x + f, d[ok, ])
fm <- design_lm(Dm, y = d$y[ok], weights = d$w[ok], offset = log(d$z[ok]))
check_methods(fm, l, "design in memory", formula = FALSE)
expect_null(formula(fr))
expect_null(formula(fm))
expect_error(design_lm(D, y = "y"), "NA/NaN/Inf in 'y'")
## ... with R expressions
l <- lm(my_transform(y) ~ x + f, d[ok, ], weights = myfun(w), offset = log(myfun(z)))
check_methods(suppressMessages(design_lm(Dk, y = my_transform(y), weights = myfun(w), offset = log(myfun(z)))),
              l, "remote design, R expressions", formula = FALSE)

## the engine with R expressions: a message with "auto", an error with "sql"
expect_message(design_lm(y ~ x + f, con, table = "d", weights = myfun(x)),
               "Evaluated in R: weights `myfun(x)` cannot be translated", fixed = TRUE)
expect_error(design_lm(y ~ x + f, con, table = "d", subset = myfun(z) > 2, engine = "sql"),
             "SQL engine cannot compute this fit: subset `myfun(z) > 2`", fixed = TRUE)
## compression stays exact with R-evaluated weights and offsets (forced:
## "auto" does not compress a weighted lm with few numeric columns)
r <- suppressMessages(design_lm(y ~ xd + f, con, table = "d", weights = myfun(xd), offset = log(myfun(z)),
                                compress = TRUE))
expect_true(r$compressed)
check_methods(r, lm(y ~ xd + f, d, weights = myfun(xd), offset = log(myfun(z))), "compressed, R expressions")

test_con_done(con)
unlink(pq)

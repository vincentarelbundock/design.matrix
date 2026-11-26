# Fixed effects in design_lm() (y ~ x | fe) against fixest::feols():
# coefficients, iid and clustered covariance matrices, residuals, number of
# observations and predictions, with feols()'s default small sample
# corrections and singleton removal. Every case runs in memory and on
# remote data (a DuckDB table, and
# Parquet). feols() demeans iteratively; a tight fixef.tol makes its results
# exact enough for a 1e-8 comparison.

if (!requireNamespace("fixest", quietly = TRUE)) exit_file("fixest not installed")
library(design.matrix)
remote <- requireNamespace("DBI", quietly = TRUE) && requireNamespace("duckdb", quietly = TRUE)

TOL <- 1e-8

fx_fit <- function(fml, d, weights = NULL, ...) {
    if (is.null(weights)) {
        fixest::feols(fml, d, fixef.tol = 1e-11, notes = FALSE, warn = FALSE, ...)
    } else {
        fixest::feols(fml, d, weights = weights, fixef.tol = 1e-11, notes = FALSE, warn = FALSE, ...)
    }
}

plain <- function(V) {
    V <- unclass(as.matrix(V))
    attributes(V) <- list(dim = dim(V))
    V
}

if (remote) {
    con <- DBI::dbConnect(duckdb::duckdb())
    pq_dir <- tempfile("fe_pq")
    dir.create(pq_dir)
}

# Write a data set to DuckDB (table) and to a Parquet file.
put <- function(d, name) {
    if (!remote) return(NULL)
    d$.row <- seq_len(nrow(d))
    DBI::dbWriteTable(con, name, d, overwrite = TRUE)
    path <- file.path(pq_dir, paste0(name, ".parquet"))
    DBI::dbExecute(con, sprintf("COPY %s TO '%s' (FORMAT PARQUET)", name, path))
    list(table = name, parquet = path)
}

# Fits of the same model in memory and on remote data.
fe_fits <- function(fml, d, src, vcov = "iid", weights = NULL, ...) {
    fits <- list(memory = do.call(design_lm, list(fml, d, vcov = vcov, ...,
                                                  weights = if (!is.null(weights)) d[[weights]])))
    if (remote && !is.null(src)) {
        ## (do.call: `weights` names a column, as a string)
        rfit <- function(source, ...) suppressMessages(do.call(design_lm, list(fml, source, vcov = vcov,
                                                                               weights = weights, ...)))
        fits$sql <- rfit(con, table = src$table, ...)
        fits$parquet <- rfit(src$parquet, ...)
    }
    fits
}

# Compare the fits of a model with a feols() fit.
check_fe <- function(fml, d, src, cluster = NULL, weights = NULL, ..., info = "") {
    vc <- if (is.null(cluster)) "iid" else cluster
    fits <- fe_fits(fml, d, src, vc, weights, ...)
    b <- fx_fit(fml, d, weights = if (is.null(weights)) NULL else stats::as.formula(paste0("~", weights)))
    fb <- unname(predict(b, d))
    for (nm in names(fits)) {
        a <- fits[[nm]]
        inf <- paste(info, nm)
        ok <- !is.na(coef(a))
        expect_equal(names(coef(a))[ok], names(coef(b)), info = paste(inf, "names"))
        expect_equal(unname(coef(a)[ok]), unname(coef(b)), tolerance = TOL, info = paste(inf, "coef"))
        expect_equal(nobs(a), nobs(b), info = paste(inf, "nobs"))
        expect_equal(plain(vcov(a)[ok, ok]),
                     plain(if (is.null(cluster)) vcov(b, vcov = "iid") else vcov(b, cluster = cluster)),
                     tolerance = TOL, info = paste(inf, "vcov"))
        expect_equal(sigma(a), unname(sqrt(b$ssr / (nobs(b) - b$nparams))), tolerance = TOL,
                     info = paste(inf, "sigma"))
        if (nm == "memory") {
            expect_equal(as.vector(resid(a))[!is.na(resid(a))], unname(resid(b)), tolerance = TOL,
                         info = paste(inf, "resid"))
            expect_equal(unname(predict(a, d)), fb, tolerance = TOL, info = paste(inf, "predict"))
        } else {
            p <- as.data.frame(predict(a, con, table = src$table))
            p <- p[order(p$.row), ]
            expect_equal(p$.fitted, fb, tolerance = TOL, info = paste(inf, "predict"))
        }
    }
    invisible(fits)
}

set.seed(20261125)
n <- 4000
d <- data.frame(
    a = sample(600, n, TRUE),
    b = sample(40, n, TRUE),
    c = sample(7, n, TRUE),
    x1 = rnorm(n), x2 = rnorm(n),
    f = factor(sample(letters[1:4], n, TRUE)),
    h = factor(sample(c("u", "v", "w"), n, TRUE)),
    w = runif(n, 0.2, 3))
d$y <- d$x1 - 0.5 * d$x2 + d$a / 300 + d$b / 20 + as.integer(d$f) + rnorm(n)
src <- put(d, "d")

# One-way ------------------------------------------------------------------
check_fe(y ~ x1 + x2 | a, d, src, info = "one-way")
check_fe(y ~ x1 + x2 | a, d, src, cluster = ~a, info = "one-way, cluster nested")
check_fe(y ~ x1 + x2 | a, d, src, cluster = ~b, info = "one-way, cluster not nested")
check_fe(y ~ x1 + x2 | b, d, src, cluster = ~b, weights = "w", info = "one-way weighted")

# Two-way and three-way ------------------------------------------------------
check_fe(y ~ x1 + x2 | a + b, d, src, cluster = ~a, info = "two-way")
check_fe(y ~ x1 + x2 | a + b, d, src, cluster = ~c, info = "two-way, cluster not nested")
check_fe(y ~ x1 + x2 | a + b + c, d, src, cluster = ~a, info = "three-way")
check_fe(y ~ x1 + x2 | a + b + c, d, src, cluster = ~c, weights = "w", info = "three-way weighted")
check_fe(y ~ x1 + x2 | a^c + b, d, src, cluster = ~b, info = "interacted FE")

# Factor regressors and interactions --------------------------------------
check_fe(y ~ x1 * f + x2 | a + b, d, src, cluster = ~a, info = "factor interaction")
check_fe(y ~ f:h + x1 | a + c, d, src, cluster = ~c, weights = "w", info = "factor x factor")
check_fe(y ~ poly(x1, 2) + log(w) | b, d, src, cluster = ~b, info = "transformations")

# Every solver of the reduced system gives the same fit
for (s in c("dense", "chol", "pcg")) {
    op <- options(design.matrix.fe_solver = s)
    a <- design_lm(y ~ x1 + x2 | a + b + c, d, vcov = ~b)
    b <- fx_fit(y ~ x1 + x2 | a + b + c, d)
    expect_equal(a$fe_info$solver, s, info = s)
    expect_equal(unname(coef(a)), unname(coef(b)), tolerance = TOL, info = s)
    expect_equal(plain(vcov(a)), plain(vcov(b, cluster = ~b)), tolerance = TOL, info = s)
    options(op)
}
# conjugate gradients: on a reduced system that is formed, IC(0) (block
# IC(0) for 3+ FE) by default, Jacobi on request; the same fit either way
# (more in test_fe_solver.R)
op <- options(design.matrix.fe_solver = "pcg", design.matrix.fe_precond = "jacobi")
a <- design_lm(y ~ x1 + x2 | a + b, d)
expect_equal(a$fe_info$precond, "jacobi")
expect_equal(unname(coef(a)), unname(coef(fx_fit(y ~ x1 + x2 | a + b, d))), tolerance = TOL)
options(op)
for (f in list(y ~ x1 + x2 | a + b, y ~ x1 + x2 | a + b + c)) {
    op <- options(design.matrix.fe_solver = "pcg", design.matrix.fe_precond = "ichol")
    a <- design_lm(f, d, vcov = ~a)
    options(op)
    expect_equal(a$fe_info$precond, "ichol")
    b <- fx_fit(f, d)
    expect_equal(unname(coef(a)), unname(coef(b)), tolerance = TOL)
    expect_equal(plain(vcov(a)), plain(vcov(b, cluster = ~a)), tolerance = TOL)
}
op <- options(design.matrix.fe_precond = "amg")
expect_error(design_lm(y ~ x1 | a + b, d), "fe_precond")
options(op)

# Unbalanced: a few large levels, many small ones ----------------------------
du <- d
du$a <- ifelse(runif(n) < 0.5, sample(3, n, TRUE), du$a)
du$b <- ifelse(runif(n) < 0.3, 1L, du$b)
check_fe(y ~ x1 + x2 | a + b, du, put(du, "du"), cluster = ~a, info = "unbalanced")

# Disconnected components of the FE graph -----------------------------------
dd <- d
half <- dd$a > 300
dd$b <- ifelse(half, dd$b + 100L, dd$b) # a in 1:300 only meets b in 1:40, a > 300 b > 100
sdd <- put(dd, "dd")
fits <- check_fe(y ~ x1 + x2 | a + b, dd, sdd, cluster = ~a, info = "disconnected")
expect_equal(fits$memory$fe_info$ncomp, 2L)
check_fe(y ~ x1 + x2 | a + b + c, dd, sdd, cluster = ~b, info = "disconnected three-way")

# Collinear regressors --------------------------------------------------------
dc <- d
dc$x3 <- dc$x1 + 2 * dc$x2                    # collinear with x1, x2
dc$xa <- (dc$a %% 13) / 7                     # constant within a: absorbed
dc$xb <- dc$b / 3 - (dc$a %% 5)              # sum of an a and a b effect
sdc <- put(dc, "dc")
fits <- check_fe(y ~ x1 + x2 + x3 | a + b, dc, sdc, cluster = ~a, info = "collinear x")
for (a in fits) expect_true(is.na(coef(a)[["x3"]]))
fits <- check_fe(y ~ x1 + xa + x2 | a + b, dc, sdc, cluster = ~a, info = "collinear with FE")
for (a in fits) expect_true(is.na(coef(a)[["xa"]]))
fits <- check_fe(y ~ xb + x1 | a + b, dc, sdc, cluster = ~b, info = "collinear with two FE")
for (a in fits) expect_true(is.na(coef(a)[["xb"]]))

# Singletons ----------------------------------------------------------------
ds <- d
ds$a[1:50] <- 10000L + seq_len(50)            # 50 singleton levels
ds$b[51] <- 999L; ds$a[51] <- 20000L          # a chain: 51 is alone in b ...
ds$a[52] <- 20000L                            # ... then 52 is alone in a
sds <- put(ds, "ds")
fits <- check_fe(y ~ x1 + x2 | a + b, ds, sds, cluster = ~b, info = "singletons")
nb <- fx_fit(y ~ x1 + x2 | a + b, ds)$nobs
for (a in fits) {
    expect_true(a$nsingletons >= 52L)
    expect_equal(a$nsingletons, nrow(ds) - nb)
}
# kept singletons: feols(fixef.rm = "none")
fits <- fe_fits(y ~ x1 + x2 | a + b, ds, sds, ~b, fixef.rm = "none")
b <- fixest::feols(y ~ x1 + x2 | a + b, ds, fixef.rm = "none", fixef.tol = 1e-11, notes = FALSE)
for (nm in names(fits)) {
    a <- fits[[nm]]
    expect_equal(nobs(a), nrow(ds), info = nm)
    expect_equal(unname(coef(a)), unname(coef(b)), tolerance = TOL, info = nm)
    expect_equal(plain(vcov(a)), plain(vcov(b, cluster = ~b)), tolerance = TOL, info = nm)
}

# Missing values and zero weights ----------------------------------------
dm <- d
dm$x1[c(3, 9)] <- NA
dm$b[17] <- NA
dm$w[c(5, 6)] <- 0
check_fe(y ~ x1 + x2 | a + b, dm, put(dm, "dm"), cluster = ~a, weights = "w", info = "NA and zero weights")
# na.exclude pads residuals and fitted values with NA for every row not used
a <- design_lm(y ~ x1 + x2 | a + b, dm, weights = w, na.action = na.exclude)
expect_equal(length(resid(a)), nrow(dm))
expect_true(all(is.na(resid(a)[c(3, 5, 6, 9, 17)])))
expect_equal(sum(!is.na(resid(a))), nobs(a))
expect_equal(unname(resid(a) + fitted(a)), ifelse(is.na(resid(a)), NA, dm$y))

# Precision: regressors with a large mean relative to their spread ---------
# One-pass aggregates cancel like a one-pass variance; centering (exact
# here: the FE absorb the shift) and one step of refinement keep the
# coefficients and standard errors accurate.
dp <- d
dp$x1 <- 1e6 + dp$x1
dp$x2 <- -3e5 + 0.5 * dp$x2
dp$y <- dp$y + dp$x1
sdp <- put(dp, "dp")
check_fe(y ~ x1 + x2 | a + b, dp, sdp, cluster = ~a, info = "mean/sd 1e6")
# with x1:f, feols() itself loses accuracy (its demeaned x1:f columns
# cancel): the reference is feols() on the centered data, mapped back
# (f_l = f_l~ - 1e6 x1:f_l)
fits <- fe_fits(y ~ x1 * f + x2 | a + b, dp, sdp)
dpc <- transform(dp, x1 = x1 - 1e6, x2 = x2 + 3e5)
bc <- fx_fit(y ~ x1 * f + x2 | a + b, dpc)
Tm <- diag(length(coef(bc)))
dimnames(Tm) <- list(names(coef(bc)), names(coef(bc)))
for (l in c("fb", "fc", "fd")) Tm[l, paste0("x1:", l)] <- -1e6
for (nm in names(fits)) {
    a <- fits[[nm]]
    expect_equal(unname(coef(a)), unname(drop(Tm %*% coef(bc))), tolerance = TOL, info = nm)
    expect_equal(plain(vcov(a)), plain(Tm %*% vcov(bc, vcov = "iid") %*% t(Tm)), tolerance = TOL, info = nm)
    expect_equal(sigma(a), sqrt(bc$ssr / (nobs(bc) - bc$nparams)), tolerance = TOL, info = nm)
}
# without centering, the covariance loses accuracy (refinement only fixes
# the coefficients)
op <- options(design.matrix.fe_center = FALSE)
a <- design_lm(y ~ x1 + x2 | a + b, dp)
b <- fx_fit(y ~ x1 + x2 | a + b, dp)
expect_equal(unname(coef(a)), unname(coef(b)), tolerance = 1e-6)
expect_true(max(abs(sqrt(diag(vcov(a))) / fixest::se(b, vcov = "iid") - 1)) > 1e-8)
options(op)

# Exact FE degrees of freedom (fixest's K.exact = TRUE) --------------------
a <- design_lm(y ~ x1 + x2 | a + b, dd)
expect_equal(a$fe_dof_exact, sum(a$fe_sizes) - 2L) # two components
b <- fx_fit(y ~ x1 + x2 | a + b, dd)
expect_equal(plain(vcov(a) * df.residual(a) / (nobs(a) - a$rank - a$fe_dof_exact)),
             plain(vcov(b, vcov = "iid", ssc = fixest::ssc(K.exact = TRUE))), tolerance = TOL)
a <- design_lm(y ~ x1 + x2 | a + b + c, d)
expect_equal(a$fe_dof_exact, a$fe_dof)

# Methods --------------------------------------------------------------------
fits <- fe_fits(y ~ x1 + x2 | a + b, d, src, ~a)
b <- fx_fit(y ~ x1 + x2 | a + b, d, cluster = ~a)
for (nm in names(fits)) {
    a <- fits[[nm]]
    expect_inherits(a, "design_lm_fe", info = nm)
    expect_inherits(a, "design_lm", info = nm)
    s <- summary(a)
    ct <- fixest::coeftable(b)
    expect_equal(unname(s$coefficients), unname(ct[, 1:4]), tolerance = 1e-7, info = nm)
    expect_equal(s$r.squared, unname(fixest::r2(b, "r2")), tolerance = 1e-8, info = nm)
    expect_equal(s$within.r.squared, unname(fixest::r2(b, "wr2")), tolerance = 1e-7, info = nm)
    expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)), tolerance = 1e-8, info = nm)
    expect_equal(df.residual(a), nobs(a) - a$rank - a$fe_dof, info = nm)
    expect_stdout(print(a))
    expect_stdout(print(s))
    # fixef(): the fixest generic dispatches; differences within an FE are
    # identified (the normalizations differ)
    fe <- fixest::fixef(a)
    expect_identical(fe, a$fixef)
    fb <- fixest::fixef(b)
    lv <- intersect(names(fe$a), names(fb$a))
    expect_equal(unname(fe$a[lv] - fe$a[lv[1L]]), unname(fb$a[lv] - fb$a[lv[1L]]), tolerance = 1e-7,
                 info = paste(nm, "fixef a"))
    lv <- intersect(names(fe$b), names(fb$b))
    expect_equal(unname(fe$b[lv] - fe$b[lv[1L]]), unname(fb$b[lv] - fb$b[lv[1L]]), tolerance = 1e-7,
                 info = paste(nm, "fixef b"))
}
# the same FE values in memory and remote
if (remote) for (nm in setdiff(names(fits), "memory"))
    expect_equal(fits[[nm]]$fixef, fits$memory$fixef, tolerance = 1e-8, info = nm)

# predict() on new data: X b plus the FE values; unseen levels give NA
nd <- d[1:6, ]
nd$a[2] <- 99999L
nd$b[5] <- 77777L
a <- fits$memory
pa <- predict(a, nd)
pb <- predict(b, nd)
expect_equal(unname(pa), unname(pb), tolerance = TOL)
expect_true(all(is.na(pa[c(2, 5)])))
if (remote) {
    nd$.row <- seq_len(nrow(nd))
    DBI::dbWriteTable(con, "nd", nd, overwrite = TRUE)
    for (nm in names(fits)) {
        p <- as.data.frame(predict(fits[[nm]], con, table = "nd"))
        p <- p[order(p$.row), ]
        expect_equal(p$.fitted, unname(pb), tolerance = TOL, info = nm)
    }
    # a remote fit predicts data in memory
    expect_equal(unname(predict(fits$sql, nd)), unname(pb), tolerance = TOL)
    # residuals() of a remote fit: the remote design with .fitted and .resid
    R <- residuals(fits$sql)
    rr <- DBI::dbGetQuery(con, sprintf("SELECT SUM(\".resid\" * \".resid\") AS rss FROM (%s) AS s", R$src$select))
    expect_equal(rr$rss, b$ssr, tolerance = 1e-8)
    # fitted() of a remote fit: a table with .fitted and .resid
    ft <- as.data.frame(fitted(fits$sql))
    expect_equal(sum(ft$.resid^2), b$ssr, tolerance = 1e-8)
    expect_equal(sort(ft$.fitted), sort(unname(fitted(b))), tolerance = TOL)
}
# interacted FE: levels named a_c, as fixest names them
a <- design_lm(y ~ x1 | a^c, d)
b <- fx_fit(y ~ x1 | a^c, d)
expect_true(all(names(fixest::fixef(b)[[1L]]) %in% names(a$fixef[[1L]])))
nd <- d[1:5, ]
nd$c[1] <- 99L
expect_equal(unname(predict(a, nd)), unname(predict(b, nd)), tolerance = TOL)
# a factor FE, and a character FE
dq <- d
dq$a <- factor(dq$a)
dq$b <- paste0("b", dq$b)
a <- design_lm(y ~ x1 + x2 | a + b, dq, vcov = ~b)
b <- fx_fit(y ~ x1 + x2 | a + b, dq)
expect_equal(unname(coef(a)), unname(coef(b)), tolerance = TOL)
expect_equal(plain(vcov(a)), plain(vcov(b, cluster = ~b)), tolerance = TOL)
expect_equal(unname(predict(a, dq[1:9, ])), unname(predict(b, dq[1:9, ])), tolerance = TOL)
if (remote) {
    sq <- put(dq, "dq")
    a2 <- design_lm(y ~ x1 + x2 | a + b, con, table = "dq", vcov = ~b)
    expect_equal(coef(a2), coef(a), tolerance = TOL)
    p <- as.data.frame(predict(a2, con, table = "dq"))
    expect_equal(p$.fitted[order(p$.row)], unname(predict(b, dq)), tolerance = TOL)
}

# Remote: response and weights as R expressions (translated to SQL, or
# evaluated in R when they have no translation), as design_lm() without FE
if (remote) {
    myw <- function(v) sqrt(v) + 0.5        # no SQL translation of a closure
    b <- fx_fit(log(y + 10) ~ x1 + x2 | a + b, transform(d, ww = sqrt(w) + 0.5), weights = ~ww)
    a1 <- design_lm(log(y + 10) ~ x1 + x2 | a + b, con, table = "d", weights = sqrt(w) + 0.5, vcov = ~a)
    a2 <- suppressMessages(design_lm(log(y + 10) ~ x1 + x2 | a + b, con, table = "d", weights = myw(w), vcov = ~a))
    a3 <- design_lm(log(y + 10) ~ x1 + x2 | a + b, d, weights = myw(w), vcov = ~a)
    for (a in list(a1, a2, a3)) {
        expect_equal(unname(coef(a)), unname(coef(b)), tolerance = TOL)
        expect_equal(plain(vcov(a)), plain(vcov(b, cluster = ~a)), tolerance = TOL)
    }
    mylog <- function(v) log(v + 10)
    a4 <- suppressMessages(design_lm(mylog(y) ~ x1 + x2 | a + b, con, table = "d", weights = "w"))
    b4 <- fx_fit(log(y + 10) ~ x1 + x2 | a + b, d, weights = ~w)
    expect_equal(unname(coef(a4)), unname(coef(b4)), tolerance = TOL)
}

# Offsets: the `offset` argument and offset() terms, as feols(offset = ) ------
bo <- fixest::feols(y ~ x1 | a + b, d, offset = ~ log(w) + 0.5 * x2, fixef.tol = 1e-11, notes = FALSE, warn = FALSE)
fo <- list(memory = design_lm(y ~ x1 + offset(0.5 * x2) | a + b, d, offset = log(w), vcov = ~c))
if (remote) {
    fo$sql <- design_lm(y ~ x1 + offset(0.5 * x2) | a + b, con, table = "d", offset = log(w), vcov = ~c)
    fo$sqlexpr <- design_lm(y ~ x1 + offset(0.5 * x2) | a + b, con, table = "d", offset = DBI::SQL("ln(w)"),
                            vcov = ~c)
    fo$parquet <- design_lm(y ~ x1 + offset(0.5 * x2) | a + b, src$parquet, offset = log(w), vcov = ~c)
}
for (nm in names(fo)) {
    a <- fo[[nm]]
    expect_equal(unname(coef(a)), unname(coef(bo)), tolerance = TOL, info = paste("offset", nm))
    expect_equal(plain(vcov(a)), plain(vcov(bo, cluster = ~c)), tolerance = TOL, info = paste("offset", nm))
    expect_equal(sigma(a), unname(sqrt(bo$ssr / (nobs(bo) - bo$nparams))), tolerance = TOL,
                 info = paste("offset", nm))
    if (nm == "memory") {
        expect_equal(unname(fitted(a)), unname(fitted(bo)), tolerance = TOL)
        expect_equal(unname(resid(a)), unname(resid(bo)), tolerance = TOL)
        expect_equal(unname(predict(a, d)), unname(predict(bo, d)), tolerance = TOL)
    } else {
        p <- as.data.frame(predict(a, con, table = "d"))
        p <- p[order(p$.row), ]
        expect_equal(p$.fitted, unname(predict(bo, d)), tolerance = TOL, info = paste("offset predict", nm))
    }
}

# Errors -----------------------------------------------------------------------
expect_error(design_lm(y ~ x1 | a, d, vcov = ~a + b), "one-way")
expect_error(design_lm(y ~ x1 | a, d, vcov = "hc1"), "vcov")
expect_error(design_lm(y ~ x1, d, vcov = ~a), "fixed effects")
expect_error(design_lm(y ~ x1 | a, d, fixef.rm = "perfect"), "fixef.rm")
expect_error(design_lm(y ~ x1 | a, d, subset = x2 > 0), "subset")
expect_error(design_lm(y ~ x1 | log(a), d), "variable names")
expect_error(predict(design_lm(y ~ x1 | a, d), d[1:3, ], se.fit = TRUE), "fit only")
expect_error(predict(design_lm(y ~ x1 | a, d), d[1:3, c("x1", "b")]), "fixed-effect variable")
if (remote) {
    expect_error(design_lm(y ~ x1 | zz, con, table = "d"), "No column")
    expect_error(design_lm(y ~ x1 | a, con, table = "d", subset = x1 > 0), "subset")
    DBI::dbDisconnect(con)
}

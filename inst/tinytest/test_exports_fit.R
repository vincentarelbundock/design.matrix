## The API that the fitting code uses (Kata design.matrix#qzqj):
## design_problem(), design_resolve(), is_remote_source(), design_wls(),
## design_patterns(fit =, check =), design_output(), design_exclude(),
## column_family(), col_irls_working(), col_const(), cbind() and
## rowSums() of columns, col_compute(table = FALSE), col_sums() with one
## na.rm per argument, design_absorb(), design_absorb_wls(),
## design_absorb_rows(), design_absorb_fitted(), absorb_control(),
## design_center(), design_uncenter(), design_groups(levels =). In memory
## and on DuckDB, against base R, lm(), glm() families and dense dummies.
## Only exported functions.
source("helpers.R")
library(xtx)

same <- function(a, b, tol = 1e-8) isTRUE(all.equal(as.vector(a), as.vector(b), tolerance = tol))

set.seed(71)
n <- 300
d <- data.frame(id = seq_len(n), y = rnorm(n), x = rnorm(n), x2 = rnorm(n) + 5, w = runif(n, 0.5, 2),
                f = factor(sample(letters[1:4], n, TRUE)), a = sample(1:15, n, TRUE),
                b = sample(1:6, n, TRUE), stringsAsFactors = FALSE)
d$yb <- rbinom(n, 1, plogis(0.4 * d$x))
d$yp <- rpois(n, exp(0.3 * d$x))
## a singleton level of a
d$a[1] <- 99L

## ---- is_remote_source(), design_resolve() ----
expect_false(is_remote_source(d))
expect_false(is_remote_source(1))
M <- design_matrix(y ~ x + f, d)
expect_equal(design_resolve(M, list(w = quote(2 * w), z = NULL), environment(), data = d)$cols$w, 2 * d$w)

## ---- design_wls() against lm() ----
fit <- design_wls(M, d$y, d$w)
ref <- lm(y ~ x + f, d, weights = w)
expect_true(same(fit$coefficients, coef(ref)))
expect_equal(names(fit$coefficients), names(coef(ref)))
expect_true(same(fit$rss, sum(d$w * residuals(ref)^2)))
expect_identical(fit$rank, ref$rank)
## aliased columns: NA as in lm()
d$x3 <- 2 * d$x
Ma <- design_matrix(y ~ x + x3 + f, d)
fa <- design_wls(Ma, d$y)
expect_identical(is.na(fa$coefficients), is.na(coef(lm(y ~ x + x3 + f, d))))
expect_null(design_wls(M, d$y, sums = FALSE)$rss)

## ---- design_patterns(fit =, check = FALSE) ----
Mp <- design_matrix(y ~ f, d)
p1 <- design_patterns(Mp, d$y, compress = TRUE)
p2 <- design_patterns(Mp, d$y, compress = TRUE, check = FALSE)
p3 <- design_patterns(Mp, d$y, compress = TRUE, fit = "glm")
expect_equal(p1, p2)
expect_equal(p1, p3)
expect_error(design_patterns(Mp, d$y, fit = "x"))

## ---- columns: col_const(), cbind(), rowSums(), col_compute(table = FALSE) ----
expect_identical(col_const(M, 2), rep(2, n))
expect_identical(col_const(d$y, 0), rep(0, n))
expect_error(col_const(M, 1:2))
expect_identical(col_compute(d$y, table = FALSE), d$y)
expect_error(col_sums(d$y, na.rm = c(TRUE, FALSE)), "one value per argument")
## one na.rm per argument
yn <- d$y
yn[3] <- NA
expect_equal(unname(col_sums(yn, yn, na.rm = c(FALSE, TRUE))), c(NA, sum(yn, na.rm = TRUE)))
expect_equal(unname(col_sums(a = yn, b = d$y, na.rm = c(TRUE, FALSE))), c(sum(yn, na.rm = TRUE), sum(d$y)))

## ---- design_output() ----
f <- drop(M %*% coef(ref))
out <- design_output(M, list(.fitted = f, .resid = d$y - f))
expect_true(is.data.frame(out))
expect_equal(out$.fitted, f)
expect_identical(design_output(M, list(.fitted = f), memory = function(dd) dd$.fitted), f)

## ---- column_family() against the family objects of stats ----
fams <- list(binomial(), binomial("probit"), binomial("cloglog"), binomial("cauchit"), poisson(), poisson("sqrt"),
             gaussian(), gaussian("log"), Gamma(), Gamma("log"), inverse.gaussian(), quasipoisson(),
             quasi(variance = "mu^2", link = "log"))
eta <- seq(-1.5, 1.5, length.out = 25)
for (fm in fams) {
    info <- paste(fm$family, fm$link)
    cf <- column_family(fm)
    expect_inherits(cf, "xtx_column_family", info = info)
    ## (the inverse link of the inverse links needs eta > 0)
    e <- if (fm$link %in% c("inverse", "1/mu^2", "sqrt", "identity") && fm$family != "gaussian") abs(eta) + 0.5 else eta
    mu <- fm$linkinv(e)
    expect_true(same(cf$linkinv(e), mu, 1e-14), info = paste(info, "linkinv"))
    expect_true(same(cf$mu.eta(e), fm$mu.eta(e), 1e-14), info = paste(info, "mu.eta"))
    expect_true(same(cf$mu_eta(e, mu), fm$mu.eta(e), 1e-14), info = paste(info, "mu_eta"))
    expect_true(same(cf$variance(mu), fm$variance(mu), 1e-14), info = paste(info, "variance"))
    expect_true(same(cf$linkfun(mu), fm$linkfun(mu), 1e-12), info = paste(info, "linkfun"))
    yy <- if (fm$family %in% c("binomial", "quasibinomial")) rep(c(0, 1), length.out = length(mu))
          else if (fm$family %in% c("poisson", "quasipoisson")) rep(0:3, length.out = length(mu)) else mu * 1.3
    expect_true(same(cf$dev(yy, mu, 1), fm$dev.resids(yy, mu, rep(1, length(mu))), 1e-12), info = paste(info, "dev"))
}
expect_error(column_family(structure(list(family = "tweedie", link = "log"), class = "family")), "supports the families")
expect_identical(column_family("poisson")$family, "poisson")
expect_null(column_family(gaussian())$boundary)

## ---- col_irls_working() against glm.fit()'s working weights ----
cf <- column_family(poisson())
e <- log(d$yp + 0.5)
mu <- cf$linkinv(e)
v <- cf$variance(mu)
me <- cf$mu_eta(e, mu)
wk <- col_irls_working(rep(1, n), e, mu, d$yp, 0, v, me)
expect_true(same(wk[[1L]], me^2 / v, 1e-14))
expect_true(same(wk[[2L]], e + (d$yp - mu) / me, 1e-14))
expect_equal(wk[[3L]][4L], n)

## ---- design_problem() ----
fitfun <- function(formula, data, weights, offset, ...) {
    design_problem(match.call(), parent.frame(), formula, data, dots = list(...))
}
P <- fitfun(y ~ x + f, d, weights = w, offset = x2)
expect_inherits(P, "xtx_problem")
expect_equal(P$w, d$w)
expect_equal(P$off, d$x2)
expect_equal(P$nall, n)
expect_true(P$intercept)
expect_error(fitfun(y ~ x, d, foo = 1), "Unknown arguments")
gfun <- function(formula, data, weights, ...) {
    design_problem(match.call(), parent.frame(), formula, data, family = column_family(binomial()), dots = list(...))
}
G <- gfun(cbind(yp, 3) ~ x, d)
expect_true(G$cols$cbind)
expect_equal(G$cols$n, d$yp + 3)

## ---- the FE estimator against dense dummies ----
P <- design_problem(quote(f(formula = y ~ x + x2, data = d)), environment(), y ~ x + x2, d, fe = list(a = "a", b = "b"))
expect_identical(names(P$keys$fe), c("a", "b"))
rw <- design_absorb_rows(P)
expect_equal(rw$nsingle, 1L)
Ga <- c(length(unique(d$a[-1])), length(unique(d$b[-1])))
expect_identical(vapply(rw$cols$codes, function(g) g$G, 1L), Ga)
x <- rw$x
est <- design_absorb_wls(x, rw$cols$z, NULL, rw$cols$codes)
expect_inherits(est, "xtx_absorb_fit")
d1 <- d[-1, ]
ref <- lm(y ~ x + x2 + factor(a) + factor(b), d1)
expect_true(same(est$b, coef(ref)[c("x", "x2")]))
expect_true(same(est$rss, sum(residuals(ref)^2)))
expect_equal(est$df.residual, ref$df.residual)
## fitted values: X b plus the FE values of the rows' levels
beta <- c(0, est$b)
fv <- design_absorb_fitted(x, beta, est$fe, rw$cols$codes)
expect_true(same(fv, fitted(ref)))
## the operator: raw cluster scores (no small sample correction)
A <- design_absorb(x, rw$cols$z, NULL, rw$cols$codes)
expect_inherits(A, "xtx_absorb")
expect_identical(unname(A$G), Ga)
## centering: the same fit, mapped back
ctr <- design_center(x, rw$cols)
expect_false(is.null(ctr))
ec <- design_uncenter(design_absorb_wls(ctr$x, ctr$cols$z, NULL, ctr$cols$codes), ctr)
expect_true(same(ec$b, est$b))
expect_true(same(ec$bread, est$bread, 1e-7))
expect_true(same(absorb_control(solver = "pcg")$solver, "pcg"))
expect_error(absorb_control(solver = "x"), "fe_solver")

## ---- design_groups(levels =): new rows coded against training levels ----
g <- design_groups(M, d["a"])
nd <- d[1:5, ]
nd$a <- c(2L, 3L, 1000L, 2L, NA)
Mn <- design_matrix(y ~ x + f, nd)
gn <- design_groups(Mn, nd["a"], levels = g$levels)
expect_identical(gn$codes, match(nd$a, g$levels$a))
expect_equal(group_lookup(seq_len(g$G) * 10, gn), match(nd$a, g$levels$a) * 10)
expect_error(group_sums(1, gn), "no level")
expect_true(same(design_absorb_fitted(Mn, c(0, 0, 0, 0, 0), list(seq_len(g$G)), list(gn)), match(nd$a, g$levels$a)))

## ---- remote: the same on DuckDB ----
if (requireNamespace("duckdb", quietly = TRUE)) {
    con <- test_con()
    DBI::dbWriteTable(con, "d", d)
    expect_true(is_remote_source(con))
    R <- design_matrix(y ~ x + f, con, table = "d")
    yr <- design_column(R, "y")
    wr <- design_column(R, "w")
    ## design_wls()
    fr <- design_wls(R, yr, wr)
    fm <- design_wls(M, d$y, d$w)
    expect_true(same(fr$coefficients, fm$coefficients))
    expect_true(same(fr$rss, fm$rss))
    ## columns
    expect_equal(sum(col_const(R, 2)), 2 * n)
    cb <- cbind(yr, 1, wr)
    expect_equal(dim(cb), c(n, 3L))
    expect_true(same(sum(rowSums(cb)), sum(d$y + 1 + d$w)))
    pc <- col_compute(yr * 2, table = FALSE)
    expect_true(same(sum(pc), 2 * sum(d$y)))
    expect_error(col_compute(yr, name = "t", table = FALSE), "only with")
    ynr <- col_ifelse(design_column(R, "id") == 3, NA, yr)
    expect_equal(unname(col_sums(ynr, ynr, yr, na.rm = c(FALSE, TRUE, FALSE))), c(NA, sum(yn, na.rm = TRUE), sum(d$y)))
    ## NaN without NA gives NaN, as in R; a count with NA rows (never NULL
    ## after !is.na()) counts them as FALSE
    expect_true(is.nan(col_sums(col_ifelse(design_column(R, "id") == 3, NaN, yr))))
    cond <- ynr > 0
    expect_equal(unname(col_sums(cond & !is.na(cond))), sum(yn > 0, na.rm = TRUE))
    expect_true(is.na(col_sums(cond)))
    ## design_patterns(check = FALSE)
    Rp <- design_matrix(y ~ f, con, table = "d")
    q1 <- design_patterns(Rp, "y", compress = TRUE)
    q2 <- design_patterns(Rp, design_column(Rp, "y"), compress = TRUE, check = FALSE)
    expect_true(same(sort(q1$z), sort(q2$z)))
    expect_true(same(sort(q1$z), sort(p1$z)))
    ## design_output(): a remote table
    fo <- R %*% coef(lm(y ~ x + f, d))
    tab <- design_output(R, list(.fitted = fo))
    expect_inherits(tab, "xtx_remote_table")
    got <- DBI::dbGetQuery(con, sprintf("SELECT id, \".fitted\" AS v FROM %s ORDER BY id", tab))
    expect_true(same(got$v, fitted(lm(y ~ x + f, d))))
    ## design_exclude(): pads a table with the rows dropped
    DBI::dbWriteTable(con, "dn", transform(d, x = ifelse(id <= 5, NA, x)))
    Rn <- design_matrix(y ~ x, con, table = "dn")
    tn <- design_output(Rn, list(.fitted = Rn %*% c(0, 1)), na = design_exclude(Rn))
    expect_equal(DBI::dbGetQuery(con, sprintf("SELECT COUNT(*) AS n FROM %s", tn))$n, n)
    expect_error(design_exclude(M), "remote design")
    ## design_resolve(): SQL, and an expression evaluated in R
    myf <- function(v) v * 3
    rr <- design_resolve(R, list(a = quote(2 * w), b = quote(myf(w))), environment())
    expect_true(same(sum(rr$cols$a), 2 * sum(d$w)))
    expect_true(same(sum(rr$cols$b), 3 * sum(d$w)))
    expect_identical(names(rr$evaluated), "b")
    ## column_family() on remote columns
    for (fm in list(binomial("probit"), binomial(), poisson(), Gamma("log"))) {
        info <- paste("remote", fm$family, fm$link)
        cf <- column_family(fm)
        er <- 0.2 * design_column(R, "x")
        e <- 0.2 * d$x
        expect_true(same(sum(cf$linkinv(er)), sum(fm$linkinv(e)), 1e-12), info = paste(info, "linkinv"))
        expect_true(same(sum(cf$mu.eta(er)), sum(fm$mu.eta(e)), 1e-12), info = paste(info, "mu.eta"))
        yv <- if (fm$family == "binomial") "yb" else "yp"
        yy <- d[[yv]] + as.numeric(fm$family == "Gamma")
        yr2 <- design_column(R, yv) + as.numeric(fm$family == "Gamma")
        expect_true(same(sum(cf$dev(yr2, cf$linkinv(er), 1)), sum(fm$dev.resids(yy, fm$linkinv(e), 1)), 1e-10),
                    info = paste(info, "dev"))
        v <- cf$valid(er, cf$linkinv(er))
        if (!is.null(v)) expect_equal(unname(col_sums(!v)), 0, info = paste(info, "valid"))
    }
    ## design_problem() on remote data: the same columns
    Pr <- fitfun(y ~ x + f, con, weights = w, offset = x2, table = "d")
    expect_inherits(Pr$x, "xtx_remote")
    expect_true(same(sum(Pr$w), sum(d$w)))
    expect_true(same(sum(Pr$off), sum(d$x2)))
    expect_equal(Pr$nall, n)
    ## the FE estimator on remote data
    Pf <- design_problem(quote(f(formula = y ~ x + x2, data = con, table = "d")), environment(), y ~ x + x2, con,
                         dots = list(table = "d"), fe = list(a = "a", b = "b"))
    rwr <- design_absorb_rows(Pf)
    expect_equal(rwr$nsingle, 1L)
    er <- design_absorb_wls(rwr$x, rwr$cols$z, NULL, rwr$cols$codes)
    expect_true(same(er$b, est$b))
    expect_true(same(er$rss, est$rss))
    fvr <- design_absorb_fitted(rwr$x, c(0, er$b), er$fe, rwr$cols$codes)
    expect_true(same(sum(fvr), sum(fitted(ref))))
    ## design_groups(levels =) remotely
    gr <- design_groups(R, "a")
    DBI::dbWriteTable(con, "nd", nd)
    Rnd <- design_matrix(y ~ x + f, con, table = "nd")
    gnr <- design_groups(Rnd, "a", levels = gr$levels)
    lk <- as.data.frame(group_lookup(seq_len(gr$G) * 10, gnr), keep = "id")
    lk <- lk[order(lk$id), ]
    expect_equal(lk[[2L]], match(nd$a, g$levels$a)[order(nd$id)] * 10)
    test_con_done(con)
}

# The compression of a GLM is weighted only when weights are given (or a
# two-column response gives trials), in memory as remotely
kp <- column_family(poisson())
kb <- column_family(binomial())
expect_false(xtx:::xtx_glm_mem_cols(c(1, 2, 3), NULL, NULL, kp)$weighted)
expect_true(xtx:::xtx_glm_mem_cols(c(1, 2, 3), c(1, 1, 2), NULL, kp)$weighted)
expect_true(xtx:::xtx_glm_mem_cols(cbind(c(1, 2), c(0, 1)), NULL, NULL, kb)$weighted)

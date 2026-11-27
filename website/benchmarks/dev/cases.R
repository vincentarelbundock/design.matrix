# Standard cases of the development benchmark (see compare.R).
#
# Each case is a function of the number of rows `n` that simulates its data
# and returns a list:
#   fit     function(): the xtx call that is timed; returns a
#           numeric vector (coefficients, or the flattened product), named
#           when the reference has names.
#   fixest  function() (optional): the fixest call that is timed (default
#           settings).
#   ref     function() (optional): the accurate reference (fixest with tight
#           tolerances, lm/glm, or a dense computation), same shape as fit().
#   ref_name  label of the reference ("fixest", "lm", "glm", "dense").
#   remote  TRUE when fit() runs on DuckDB: queries are counted.
#   setup   function() (optional): run once before timing (options).
#
# Constructors that are called in a loop must force() their arguments (lazy
# evaluation would otherwise give every case the last loop values).
# Cases live in `cases`, keyed by name; `case_group` maps each name to a group
# (fe, glm, nofe, remote, products), which --cases selects (names work too).
# To add a case: write the function, add it to `cases` and `case_group`.
# Remote cases use `con`, a DuckDB connection with `threads` set, created by
# worker.R (bench_con()).

## Data -------------------------------------------------------------------

# Panel of the fixest / pyfixest benchmarks (Kyle Butts' fixest_benchmarks):
# individuals over 10 years, about 23 individuals per firm. "simple": each
# row's firm is drawn at random; "difficult": firms assigned round robin
# (individuals move every year, poorly separated FE, slow demeaning).
# Outcomes: y (gaussian), y_pois (Poisson), y_bin (logit).
panel_dgp <- function(n, type, k = 1, seed = 42) {
    set.seed(seed)
    nb_year <- 10
    nb_indiv <- round(n / nb_year)
    nb_firm <- round(nb_indiv / 23)
    n <- nb_indiv * nb_year
    indiv_id <- rep(seq_len(nb_indiv), each = nb_year)
    year <- rep(seq_len(nb_year), nb_indiv)
    firm_id <- if (type == "simple") sample.int(nb_firm, n, TRUE) else rep_len(seq_len(nb_firm), n)
    d <- data.frame(y = numeric(n))
    xb <- numeric(n)
    for (j in seq_len(k)) {
        x <- rnorm(n)
        d[[paste0("x", j)]] <- x
        xb <- xb + x / j
    }
    firm_fe <- rnorm(nb_firm)
    unit_fe <- rnorm(nb_indiv)
    year_fe <- rnorm(nb_year)
    fe <- firm_fe[firm_id] + unit_fe[indiv_id] + year_fe[year]
    d$y <- xb + fe + rnorm(n)
    # smaller FE for the GLMs, so that few groups are all-zero / perfect fit
    eta <- 0.5 * xb + 0.3 * fe
    d$y_pois <- rpois(n, exp(eta))
    d$y_bin <- rbinom(n, 1, plogis(eta))
    d$indiv_id <- indiv_id
    d$year <- year
    d$firm_id <- firm_id
    d
}

# Cross-section with factors (website/benchmarks/setup.R style). `g` is a
# discrete numeric (5 values), so that y ~ g * f1 + f2 has few patterns.
factor_dgp <- function(n, seed = 1) {
    set.seed(seed)
    d <- data.frame(
        x = rnorm(n),
        g = sample(1:5, n, TRUE),
        f1 = factor(sample(sprintf("a%02d", 1:50), n, TRUE)),
        f2 = factor(sample(letters[1:20], n, TRUE)),
        f3 = factor(sample(sprintf("c%03d", 1:200), n, TRUE)))
    eta <- 0.5 * d$x + 0.1 * d$g + rnorm(50, sd = 0.3)[d$f1] + rnorm(20, sd = 0.3)[d$f2]
    d$y <- eta + rnorm(n)
    d$y_bin <- rbinom(n, 1, plogis(eta))
    d
}

remote_table <- function(d, name = "d") {
    DBI::dbWriteTable(con, name, d, overwrite = TRUE)
    name
}

## Cases ------------------------------------------------------------------

fe_formula <- function(lhs, nfe) {
    fe <- c("indiv_id", "year", "firm_id")[seq_len(nfe)]
    f <- as.formula(paste(lhs, "~ x1 |", paste(fe, collapse = " + ")))
    environment(f) <- globalenv()
    f
}

case_fe_ols <- function(type, nfe, remote = FALSE) {
    force(type); force(nfe); force(remote)  # loop variables: evaluate now
    function(n) {
    d <- panel_dgp(n, type)
    f <- fe_formula("y", nfe)
    tab <- if (remote) remote_table(d)
    list(
        fit = if (remote) function() coef(xtx.fit::lm_xtx(f, con, table = tab))["x1"]
              else function() coef(xtx.fit::lm_xtx(f, d))["x1"],
        fixest = function() fixest::feols(f, d, vcov = "iid", notes = FALSE),
        ref = function() coef(fixest::feols(f, d, vcov = "iid", notes = FALSE,
                                            fixef.tol = 1e-11, fixef.iter = 1e5))["x1"],
        ref_name = "fixest", remote = remote)
    }
}

case_fe_glm <- function(family, type, nfe, remote = FALSE) {
    force(family); force(type); force(nfe); force(remote)
    function(n) {
    d <- panel_dgp(n, type)
    lhs <- if (family == "poisson") "y_pois" else "y_bin"
    fam <- if (family == "poisson") stats::poisson() else stats::binomial()
    f <- fe_formula(lhs, nfe)
    tab <- if (remote) remote_table(d)
    fx <- function(...) {
        if (family == "poisson") fixest::fepois(f, d, vcov = "iid", notes = FALSE, ...)
        else fixest::feglm(f, d, family = "logit", vcov = "iid", notes = FALSE, ...)
    }
    list(
        fit = if (remote) function() coef(xtx.fit::glm_xtx(f, fam, con, table = tab))["x1"]
              else function() coef(xtx.fit::glm_xtx(f, fam, d))["x1"],
        fixest = function() fx(),
        ref = function() coef(fx(fixef.tol = 1e-11, fixef.iter = 1e5, glm.tol = 1e-11, glm.iter = 100))["x1"],
        ref_name = "fixest", remote = remote)
    }
}

f_nofe <- y ~ x * f1 + f2
f_nofe_bin <- y_bin ~ x * f1 + f2
f_disc <- y ~ g * f1 + f2

case_lm_nofe <- function(remote = FALSE, compress = "auto", formula = f_nofe) function(n) {
    d <- factor_dgp(n)
    tab <- if (remote) remote_table(d)
    list(
        setup = function() if (remote) options(xtx.remote_compress = compress),
        fit = if (remote) function() coef(xtx.fit::lm_xtx(formula, con, table = tab))
              else function() coef(xtx.fit::lm_xtx(formula, d, compress = compress)),
        fixest = function() fixest::feols(formula, d, vcov = "iid", notes = FALSE),
        ref = function() coef(stats::lm(formula, d)),
        ref_name = "lm", remote = remote)
}

case_glm_nofe <- function(remote = FALSE) function(n) {
    d <- factor_dgp(n)
    tab <- if (remote) remote_table(d)
    list(
        fit = if (remote) function() coef(xtx.fit::glm_xtx(f_nofe_bin, stats::binomial(), con, table = tab))
              else function() coef(xtx.fit::glm_xtx(f_nofe_bin, stats::binomial(), d)),
        fixest = function() fixest::feglm(f_nofe_bin, d, family = "logit", vcov = "iid", notes = FALSE),
        ref = function() coef(stats::glm(f_nofe_bin, stats::binomial(), d,
                                         control = list(epsilon = 1e-12, maxit = 100))),
        ref_name = "glm", remote = remote)
}

# Products on a factor-interaction design: X'X and colSums(w * X, group = f3)
case_crossprod <- function(remote = FALSE) function(n) {
    d <- factor_dgp(n)
    tab <- if (remote) remote_table(d)
    list(
        fit = if (remote) function() c(crossprod(xtx::design_matrix(f_nofe, con, table = tab)))
              else function() c(crossprod(xtx::design_matrix(f_nofe, d))),
        ref = function() c(crossprod(stats::model.matrix(f_nofe, d))),
        ref_name = "dense", remote = remote)
}

case_colsums <- function(remote = FALSE) function(n) {
    d <- factor_dgp(n)
    tab <- if (remote) remote_table(d)
    list(
        fit = if (remote) function() c(colSums("x" * xtx::design_matrix(f_nofe, con, table = tab), group = "f3"))
              else function() c(colSums(d$x * xtx::design_matrix(f_nofe, d), group = d$f3)),
        ref = function() c(rowsum(d$x * stats::model.matrix(f_nofe, d), d$f3)),
        ref_name = "dense", remote = remote)
}

cases <- list()
case_group <- character()
add_case <- function(name, group, f) {
    cases[[name]] <<- f
    case_group[[name]] <<- group
}

for (type in c("simple", "difficult")) for (nfe in 2:3) {
    add_case(sprintf("fe_ols_%s_%dfe", type, nfe), "fe", case_fe_ols(type, nfe))
    for (fam in c("poisson", "logit"))
        add_case(sprintf("fe_%s_%s_%dfe", fam, type, nfe), "glm", case_fe_glm(fam, type, nfe))
}
add_case("lm_nofe", "nofe", case_lm_nofe())
add_case("glm_nofe", "nofe", case_glm_nofe())
add_case("lm_disc", "nofe", case_lm_nofe(formula = f_disc))
add_case("remote_lm_nofe", "remote", case_lm_nofe(remote = TRUE, compress = FALSE))
add_case("remote_lm_disc_compress", "remote", case_lm_nofe(remote = TRUE, compress = TRUE, formula = f_disc))
add_case("remote_lm_disc_nocompress", "remote", case_lm_nofe(remote = TRUE, compress = FALSE, formula = f_disc))
add_case("remote_glm_nofe", "remote", case_glm_nofe(remote = TRUE))
add_case("remote_fe_ols_simple_2fe", "remote", case_fe_ols("simple", 2, remote = TRUE))
add_case("remote_fe_poisson_simple_2fe", "remote", case_fe_glm("poisson", "simple", 2, remote = TRUE))
add_case("crossprod", "products", case_crossprod())
add_case("crossprod_remote", "products", case_crossprod(remote = TRUE))
add_case("colsums_group", "products", case_colsums())
add_case("colsums_group_remote", "products", case_colsums(remote = TRUE))

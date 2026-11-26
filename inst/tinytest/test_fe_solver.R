# The solver of the reduced FE system (R/fe.R: xtx_fe_system(), xtx_fe_solve(),
# xtx_fe_pcg()) on small replicas of the worker x firm panels of the fixest
# benchmarks: "simple" (a random firm per row: high mobility, the reduced
# system near dense, applied matrix-free) and "difficult" (firms assigned
# round robin: low mobility, the reduced system sparse and ill conditioned,
# formed and preconditioned by IC(0)). Each path is forced in turn
# (xtx.fe_form_max = 0: matrix-free; Inf: formed) and compared
# with feols() at a tight fixef.tol, in memory and on remote data, without
# warnings. With 132 firms (even) and 10 years, round robin makes firm
# parity determine year parity: the FE are collinear beyond the references
# of the pairwise graphs, R is singular, and the solver must drop a level.

if (!requireNamespace("fixest", quietly = TRUE)) exit_file("fixest not installed")
source("helpers.R")
library(xtx)
remote <- requireNamespace("DBI", quietly = TRUE) && requireNamespace("duckdb", quietly = TRUE)

TOL <- 1e-8

dgp <- function(n, type, seed = 42) {
    set.seed(seed)
    nb_year <- 10; nb_indiv <- round(n / nb_year); nb_firm <- round(nb_indiv / 23)
    n <- nb_indiv * nb_year
    indiv_id <- rep(seq_len(nb_indiv), each = nb_year)
    year <- rep(seq_len(nb_year), nb_indiv)
    firm_id <- if (type == "simple") sample.int(nb_firm, n, TRUE) else rep_len(seq_len(nb_firm), n)
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    y <- x1 - 0.5 * x2 + rnorm(nb_firm)[firm_id] + rnorm(nb_indiv)[indiv_id] + rnorm(nb_year)[year] + rnorm(n)
    data.frame(y, x1, x2, indiv_id, year, firm_id)
}

plain <- function(V) {
    V <- unclass(as.matrix(V))
    attributes(V) <- list(dim = dim(V))
    V
}

# The value of an expression and the warnings of the solver it gave.
warns <- function(expr) {
    ws <- character(0)
    val <- withCallingHandlers(expr, warning = function(w) {
        ws <<- c(ws, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    list(value = val, warnings = grep("onjugate|converge|singular", ws, value = TRUE))
}

if (remote) con <- test_con()

f3 <- y ~ x1 + x2 | indiv_id + year + firm_id
paths <- list(matfree = list(xtx.fe_solver = "pcg", xtx.fe_form_max = 0),
              formed = list(xtx.fe_solver = "pcg", xtx.fe_form_max = Inf),
              auto = list())
for (type in c("simple", "difficult")) {
    d <- dgp(30360, type)
    b <- fixest::feols(f3, d, fixef.tol = 1e-11, fixef.iter = 1e5, notes = FALSE, warn = FALSE)
    if (remote) DBI::dbWriteTable(con, type, d, overwrite = TRUE)
    for (nm in names(paths)) {
        op <- options(paths[[nm]])
        inf <- paste(type, nm)
        r <- warns(design_lm(f3, d, vcov = ~firm_id))
        a <- r$value
        expect_identical(r$warnings, character(0), info = inf)
        if (nm == "matfree") expect_false(a$fe_info$formed, info = inf)
        if (nm == "formed") {
            expect_true(a$fe_info$formed, info = inf)
            expect_equal(a$fe_info$precond, "ichol", info = inf)
        }
        expect_equal(unname(coef(a)), unname(coef(b)), tolerance = TOL, info = paste(inf, "coef"))
        expect_equal(plain(vcov(a)), plain(vcov(b, cluster = ~firm_id)), tolerance = TOL,
                     info = paste(inf, "vcov"))
        expect_equal(sigma(a), unname(sqrt(b$ssr / (nobs(b) - b$nparams))), tolerance = TOL,
                     info = paste(inf, "sigma"))
        if (remote && nm != "auto") {
            for (eng in "sql") {
                r <- warns(suppressMessages(design_lm(f3, con, table = type, vcov = ~firm_id,
                                                      engine = eng, chunk_size = 7000)))
                ar <- r$value
                expect_identical(r$warnings, character(0), info = paste(inf, eng))
                expect_equal(ar$fe_info$formed, a$fe_info$formed, info = paste(inf, eng))
                expect_equal(unname(coef(ar)), unname(coef(b)), tolerance = TOL, info = paste(inf, eng, "coef"))
                expect_equal(plain(vcov(ar)), plain(vcov(a)), tolerance = TOL, info = paste(inf, eng, "vcov"))
            }
        }
        options(op)
    }
}
# the default paths: the simple panel matrix-free, the difficult one formed
op <- options(xtx.fe_solver = "pcg")
expect_false(design_lm(f3, dgp(2e5, "simple"))$fe_info$formed)
a <- design_lm(f3, dgp(30360, "difficult"))
expect_true(a$fe_info$formed)
expect_equal(a$fe_info$precond, "ichol")
options(op)

# the dense solver (pivoted Cholesky) on the singular system of the
# difficult panel: one more level dropped, and the exact rank of the FE
op <- options(xtx.fe_solver = "dense")
d <- dgp(30360, "difficult")
a <- design_lm(f3, d)
options(op)
expect_equal(a$fe_info$solver, "dense")
expect_equal(a$fe_dof_exact, sum(a$fe_sizes) - 3L)  # two references and firm parity = year parity
b <- fixest::feols(f3, d, fixef.tol = 1e-11, fixef.iter = 1e5, notes = FALSE, warn = FALSE)
expect_equal(unname(coef(a)), unname(coef(b)), tolerance = TOL)

# Jacobi on request on the formed system; two FE
d <- dgp(30360, "difficult")
for (f in list(f3, y ~ x1 + x2 | indiv_id + firm_id)) {
    b <- fixest::feols(f, d, fixef.tol = 1e-11, fixef.iter = 1e5, notes = FALSE, warn = FALSE)
    op <- options(xtx.fe_solver = "pcg", xtx.fe_precond = "jacobi",
                  xtx.fe_form_max = Inf)
    r <- warns(design_lm(f, d))
    options(op)
    a <- r$value
    expect_identical(r$warnings, character(0))
    expect_equal(a$fe_info$precond, "jacobi")
    expect_equal(unname(coef(a)), unname(coef(b)), tolerance = TOL)
}

# Options
op <- options(xtx.fe_form_max = -1)
expect_error(design_lm(f3, d), "fe_form_max")
options(op)
op <- options(xtx.fe_precond = "amg")
expect_error(design_lm(f3, d), "fe_precond")
options(op)

# Warm start: xtx_fe_apply() from the solution of a nearby system gives the
# same solution in fewer iterations, on both paths
ns <- asNamespace("xtx")
d <- dgp(30360, "simple")
codes <- lapply(d[c("indiv_id", "year", "firm_id")], function(v) match(v, sort(unique(v))))
G <- vapply(codes, max, 1L)
gsum <- function(v) lapply(codes, function(cd) matrix(rowsum(v, cd, reorder = TRUE)[, 1L], ncol = 1L))
for (fm in c(0, Inf)) {
    opts <- absorb_control()
    opts$solver <- "pcg"
    opts$form_max <- fm
    for (w in list(rep(1, nrow(d)), runif(nrow(d), 0.5, 2))) {
        nk <- lapply(codes, function(cd) rowsum(w, cd, reorder = TRUE)[, 1L])
        crosstab <- function(a, b) Matrix::sparseMatrix(i = codes[[b]], j = codes[[a]], x = w,
                                                       dims = c(G[b], G[a]))
        sys <- ns$xtx_fe_system(G, nk, crosstab, opts)
        expect_equal(sys$info$formed, is.infinite(fm))
        v <- rnorm(nrow(d))
        s0 <- ns$xtx_fe_apply(sys, gsum(w * v))
        v2 <- v + 1e-3 * rnorm(nrow(d))
        cold <- ns$xtx_fe_apply(sys, gsum(w * v2))
        warm <- ns$xtx_fe_apply(sys, gsum(w * v2), Z0 = s0$Z)
        expect_true(warm$info$iter < cold$info$iter, info = paste(fm, warm$info$iter, cold$info$iter))
        expect_equal(warm$Q, cold$Q, tolerance = 1e-10)
        ## (the FE values are unique up to normalization: compare their sums)
        fsum <- function(A) Reduce(`+`, lapply(seq_along(G), function(k) A[[k]][codes[[k]], 1L]))
        expect_equal(fsum(warm$A), fsum(cold$A), tolerance = 1e-8)
        ## a starting value worse than zero is dropped
        bad <- ns$xtx_fe_apply(sys, gsum(w * v2), Z0 = 1e6 * s0$Z + 1)
        expect_equal(bad$info$iter, cold$info$iter)
        expect_equal(fsum(bad$A), fsum(cold$A), tolerance = 1e-8)
    }
}

if (remote) test_con_done(con)

# The solver of the reduced FE system (R/absorb.R: xtx_fe_system(),
# xtx_fe_apply()) on a small replica of the worker x firm panels of the
# fixest benchmarks ("simple": a random firm per row, high mobility). The
# fits through each path of the solver, compared with feols(), are tests
# of xtx.fit.
source("helpers.R")
library(xtx)

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


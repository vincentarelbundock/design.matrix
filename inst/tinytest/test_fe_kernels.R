# The C kernels of the fixed-effect fits (src/fe.c) against the R code
# they replace.
source("helpers.R")
library(design.matrix)
ns <- asNamespace("design.matrix")

# ---- FE codes (dm_fe_code()): direct indexing or hashing, rows per level
# the codes as sort(unique()) and match() give them (factors: level order)
code_ref <- function(v) {
    if (is.factor(v)) {
        used <- tabulate(as.integer(v), nlevels(v)) > 0L
        return(list(codes = cumsum(used)[as.integer(v)], G = sum(used), levels = levels(v)[used]))
    }
    if (is.logical(v)) v <- as.integer(v)
    u <- sort(unique(v), method = "radix")
    list(codes = match(v, u), G = length(u), levels = u)
}
set.seed(1)
n <- 5000
cases <- list(
    int = sample(1:50, n, TRUE),
    int_wide = sample(c(1L, 1000000000L, -5L, 7L), n, TRUE),
    dbl = sample(c(0.5, -0, 0, 3, 1e300, -2.25), n, TRUE),
    dbl_int = as.double(sample(-3:30, n, TRUE)),
    dbl_huge = sample(c(2^60, -2^60, 1), n, TRUE),
    chr = sample(c("b", "a", "Z", "B", "aa", ""), n, TRUE),
    chr_runs = rep(sprintf("id%03d", sample(200)), each = 25),
    fac = factor(sample(letters, n, TRUE), levels = rev(letters)),
    fac_unused = factor(sample(c("x", "y"), n, TRUE), levels = c("y", "q", "x")),
    lgl = sample(c(TRUE, FALSE), n, TRUE),
    const = rep(7, n),
    many = sample(1e5, n, TRUE) + 0.5)
for (nm in names(cases)) {
    v <- cases[[nm]]
    a <- ns$dm_fe_code(v)
    b <- code_ref(v)
    expect_identical(a$codes, b$codes, info = nm)
    expect_identical(a$G, b$G, info = nm)
    expect_identical(a$levels, b$levels, info = nm)
    expect_identical(a$n, tabulate(b$codes, b$G), info = nm)
}
# a string in two encodings: one level, as match() has it
lat <- rawToChar(as.raw(c(0x63, 0x61, 0x66, 0xe9)))
Encoding(lat) <- "latin1"
enc <- c(enc2utf8(lat), lat, "tea")[sample(3, 300, TRUE)]
expect_identical(ns$dm_fe_code(enc)[c("codes", "G", "levels")], code_ref(enc))
# interactions: the combined codes and their levels
d <- data.frame(a = sample(c("u", "v", "w"), n, TRUE), b = sample(c(10.5, 2, 3), n, TRUE))
ct <- ns$dm_fe_code_term(d, c("a", "b"))
key <- paste(d$a, d$b)
lv <- unique(d[order(d$a, d$b, method = "radix"), ])
expect_identical(ct$levels$a[ct$codes], d$a)
expect_identical(ct$levels$b[ct$codes], d$b)
expect_identical(nrow(unique(ct$levels)), ct$G)
expect_identical(ct$n, tabulate(ct$codes, ct$G))
expect_equal(unname(as.list(ct$levels)), unname(as.list(lv)), check.attributes = FALSE)

# ---- sums of the FE fits on several threads (dm_fe_gsum(),
# dm_fe_struct_gsum() with dense crosstabs, dm_fe_crosstabs())
old_threads <- design_threads()
set.seed(7)
n <- 40000
fd <- data.frame(id = rep(1:4000, each = 10), t = rep(1:10, 4000), firm = sample(300, n, TRUE),
                 x1 = rnorm(n), x2 = rnorm(n), w = runif(n) + 0.1)
fd$y <- fd$x1 - fd$x2 + rnorm(300)[fd$firm] + rnorm(4000)[fd$id] + rnorm(n)
fd$cnt <- rpois(n, exp(0.3 * fd$x1))
cds <- lapply(fd[c("id", "firm", "t")], function(v) ns$dmc_codes(v, max(v)))
perm <- sample(n)

# crosstabs: the sums of w by pair of levels, in row order, as the sparse
# matrix of the pairs; the same entries on any number of threads
ct_ref <- function(a, b, w) {
    M <- Matrix::sparseMatrix(i = b, j = a, x = if (is.null(w)) rep(1, length(a)) else w,
                              dims = c(max(b), max(a)), repr = "C")
    list(p = M@p, i = M@i)
}
for (wt in list(NULL, fd$w)) for (pm in list(seq_len(n), perm)) {
    a <- cds$firm$codes[pm]
    bs <- list(cds$id$codes[pm], cds$t$codes[pm])
    ww <- if (!is.null(wt)) wt[pm]
    c1 <- .Call(ns$C_dm_fe_crosstabs, a, bs, ww, 300L, c(4000L, 10L), 1L)
    c4 <- .Call(ns$C_dm_fe_crosstabs, a, bs, ww, 300L, c(4000L, 10L), 4L)
    expect_true(strict_identical(c1, c4))
    for (k in 1:2) {
        r <- ct_ref(a, bs[[k]], ww)
        expect_identical(c1[[k]][[1L]], r$p)
        expect_identical(c1[[k]][[2L]], r$i)
        M <- Matrix::sparseMatrix(i = bs[[k]], j = a, x = if (is.null(ww)) rep(1, n) else ww, repr = "C",
                                  dims = c(max(bs[[k]]), 300L))
        expect_equal(c1[[k]][[3L]], M@x, tolerance = 1e-13)
    }
}
# the dense crosstabs of the grouped-sums pass: the entries of the sparse
# kernel, bit for bit (one thread)
X <- design_matrix(~ x1 + x2, fd)
design_threads(1)
for (wt in list(NULL, fd$w)) {
    r <- ns$dm_gsums_memory(X, wt, cds, fd$y, list(c(2L, 3L), c(3L, 2L)))
    s <- .Call(ns$C_dm_fe_crosstabs, cds$firm$codes, list(cds$t$codes), wt, 300L, 10L, 1L)[[1L]]
    expect_true(strict_identical(r$ct[[1L]], ns$dm_ct_matrix(s, 10L, 300L)))
    s2 <- .Call(ns$C_dm_fe_crosstabs, cds$t$codes, list(cds$firm$codes), wt, 10L, 300L, 1L)[[1L]]
    expect_true(strict_identical(r$ct[[2L]], ns$dm_ct_matrix(s2, 300L, 10L)))
    # and the grouped sums those of the pass without crosstabs
    expect_true(strict_identical(r$gs, ns$dm_gsums_memory(X, wt, cds, fd$y)$gs))
}
# grouped sums: one thread as rowsum(); a given number of threads repeats
W <- cbind(fd$x1, fd$w)
g1 <- .Call(ns$C_dm_fe_gsum, cds$id$codes, 4000L, W, 1L)
expect_equal(g1, unname(rowsum(W, cds$id$codes)), tolerance = 1e-13)
g4 <- .Call(ns$C_dm_fe_gsum, cds$id$codes, 4000L, W, 4L)
expect_true(strict_identical(g4, .Call(ns$C_dm_fe_gsum, cds$id$codes, 4000L, W, 4L)))
expect_equal(g4, g1, tolerance = 1e-12)
for (wt in list(NULL, fd$w)) {
    design_threads(4)
    a <- ns$dm_gsums_memory(X, wt, cds, fd$y, list(c(2L, 3L)))
    b <- ns$dm_gsums_memory(X, wt, cds, fd$y, list(c(2L, 3L)))
    expect_true(strict_identical(a, b))
    design_threads(1)
    c1 <- ns$dm_gsums_memory(X, wt, cds, fd$y, list(c(2L, 3L)))
    expect_equal(a, c1, tolerance = 1e-12)
}

# FE fits: identical for a given number of threads; up to rounding across
# numbers of threads; with one thread (and the refinement step always
# taken) the values of the package before the threaded sums, stored as
# hexadecimal strings by sprintf("%a") (local runs only: BLAS): the file
# holds fe_fits() as computed at ba44828 with this file's code.
fe_fits <- function() {
    hex <- function(f) sprintf("%a", unname(c(coef(f), as.vector(vcov(f)), f$rss, unlist(f$fixef),
                                              f$deviance)))
    pd <- fd[perm, ]
    list(f3 = hex(design_lm(y ~ x1 + x2 | id + firm + t, fd)),
         f3w = hex(design_lm(y ~ x1 + x2 | id + firm + t, fd, weights = fd$w, vcov = ~firm)),
         f3perm = hex(design_lm(y ~ x1 + x2 | id + firm + t, pd, vcov = ~id)),
         f2 = hex(design_lm(y ~ x1 | firm + t, fd, vcov = ~t)),
         pois = hex(design_glm(cnt ~ x1 + x2 | id + firm, fd, family = poisson())))
}
op <- options(design.matrix.fe_refine = 1L)
design_threads(1)
r1 <- fe_fits()
expect_true(strict_identical(r1, fe_fits()))
ref_file <- "fe_threads_ref.rds"
if (at_home() && file.exists(ref_file)) {
    ref <- readRDS(ref_file)
    for (k in names(ref)) expect_identical(r1[[k]], ref[[k]], info = paste("1 thread, stored:", k))
}
if (.Call(ns$C_dm_omp_max_threads) > 1L) {
    design_threads(min(4L, .Call(ns$C_dm_omp_max_threads)))
    ra <- fe_fits()
    expect_true(strict_identical(ra, fe_fits()))
    num <- function(h) suppressWarnings(as.numeric(h))
    for (k in names(r1)) expect_equal(num(ra[[k]]), num(r1[[k]]), tolerance = 1e-9, info = k)
    expect_true(any(!mapply(identical, ra, r1)))
}
options(op)
design_threads(old_threads)

# ---- the refinement step (design.matrix.fe_refine = "auto", default): its
# correction is always applied; the residuals are computed again only when
# the correction is not negligible. Coefficients and FE values are those of
# one step with the residuals recomputed (fe_refine = 1); the residual sum
# of squares and the covariance agree to rounding.
fd$big4 <- fd$x2 + 1e4
fd$big6 <- fd$x2 + 1e6
cases <- list(simple = list(y ~ x1 + x2 | id + firm + t, "iid", TRUE),
              cluster = list(y ~ x1 + x2 | id + firm + t, ~firm, TRUE),
              uncentered4 = list(y ~ x1 + big4 | id + firm + t, "iid", FALSE),
              uncentered6 = list(y ~ x1 + big6 | id + firm + t, ~firm, FALSE))
pass <- character()
for (nm in names(cases)) {
    cs <- cases[[nm]]
    op <- options(design.matrix.fe_center = cs[[3]])
    a <- design_lm(cs[[1]], fd, vcov = cs[[2]])
    options(design.matrix.fe_refine = 1L)
    b <- design_lm(cs[[1]], fd, vcov = cs[[2]])
    options(op)
    options(design.matrix.fe_refine = NULL)
    expect_identical(coef(a), coef(b), info = nm)
    expect_identical(a$fixef, b$fixef, info = nm)
    expect_equal(a$rss, b$rss, tolerance = 1e-14, info = nm)
    expect_equal(vcov(a), vcov(b), tolerance = 1e-10, info = nm)
    expect_identical(b$fe_info$resid_pass, "done", info = nm)
    pass[nm] <- a$fe_info$resid_pass
}
# skipped on the well-conditioned panels, not where the step corrects
# much (variables far from zero, not centered: the one-pass aggregates
# lose digits, b moves by 1e-6 and 4e-3 standard errors)
expect_identical(unname(pass), c("skipped", "skipped", "done", "done"))

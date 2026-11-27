# The C kernels of the fixed-effect fits (src/fe.c) against the R code
# they replace.
source("helpers.R")
library(xtx)
ns <- asNamespace("xtx")

# ---- FE codes (xtx_fe_code()): direct indexing or hashing, rows per level
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
    a <- ns$xtx_fe_code(v)
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
expect_identical(ns$xtx_fe_code(enc)[c("codes", "G", "levels")], code_ref(enc))
# interactions: the combined codes and their levels
d <- data.frame(a = sample(c("u", "v", "w"), n, TRUE), b = sample(c(10.5, 2, 3), n, TRUE))
ct <- ns$xtx_fe_code_term(d, c("a", "b"))
key <- paste(d$a, d$b)
lv <- unique(d[order(d$a, d$b, method = "radix"), ])
expect_identical(ct$levels$a[ct$codes], d$a)
expect_identical(ct$levels$b[ct$codes], d$b)
expect_identical(nrow(unique(ct$levels)), ct$G)
expect_identical(ct$n, tabulate(ct$codes, ct$G))
expect_equal(unname(as.list(ct$levels)), unname(as.list(lv)), check.attributes = FALSE)

# ---- sums of the FE fits on several threads (xtx_fe_gsum(),
# xtx_fe_struct_gsum() with dense crosstabs, xtx_fe_crosstabs())
old_threads <- design_threads()
set.seed(7)
n <- 40000
fd <- data.frame(id = rep(1:4000, each = 10), t = rep(1:10, 4000), firm = sample(300, n, TRUE),
                 x1 = rnorm(n), x2 = rnorm(n), w = runif(n) + 0.1)
fd$y <- fd$x1 - fd$x2 + rnorm(300)[fd$firm] + rnorm(4000)[fd$id] + rnorm(n)
fd$cnt <- rpois(n, exp(0.3 * fd$x1))
cds <- lapply(fd[c("id", "firm", "t")], function(v) ns$xtx_c_codes(v, max(v)))
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
    c1 <- .Call(ns$C_xtx_fe_crosstabs, a, bs, ww, 300L, c(4000L, 10L), 1L)
    c4 <- .Call(ns$C_xtx_fe_crosstabs, a, bs, ww, 300L, c(4000L, 10L), 4L)
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
    r <- ns$xtx_gsums_memory(X, wt, cds, fd$y, list(c(2L, 3L), c(3L, 2L)))
    s <- .Call(ns$C_xtx_fe_crosstabs, cds$firm$codes, list(cds$t$codes), wt, 300L, 10L, 1L)[[1L]]
    expect_true(strict_identical(r$ct[[1L]], ns$xtx_ct_matrix(s, 10L, 300L)))
    s2 <- .Call(ns$C_xtx_fe_crosstabs, cds$t$codes, list(cds$firm$codes), wt, 10L, 300L, 1L)[[1L]]
    expect_true(strict_identical(r$ct[[2L]], ns$xtx_ct_matrix(s2, 300L, 10L)))
    # and the grouped sums those of the pass without crosstabs
    expect_true(strict_identical(r$gs, ns$xtx_gsums_memory(X, wt, cds, fd$y)$gs))
}
# grouped sums: one thread as rowsum(); a given number of threads repeats
W <- cbind(fd$x1, fd$w)
g1 <- .Call(ns$C_xtx_fe_gsum, cds$id$codes, 4000L, W, 1L)
expect_equal(g1, unname(rowsum(W, cds$id$codes)), tolerance = 1e-13)
g4 <- .Call(ns$C_xtx_fe_gsum, cds$id$codes, 4000L, W, 4L)
expect_true(strict_identical(g4, .Call(ns$C_xtx_fe_gsum, cds$id$codes, 4000L, W, 4L)))
expect_equal(g4, g1, tolerance = 1e-12)
for (wt in list(NULL, fd$w)) {
    design_threads(4)
    a <- ns$xtx_gsums_memory(X, wt, cds, fd$y, list(c(2L, 3L)))
    b <- ns$xtx_gsums_memory(X, wt, cds, fd$y, list(c(2L, 3L)))
    expect_true(strict_identical(a, b))
    design_threads(1)
    c1 <- ns$xtx_gsums_memory(X, wt, cds, fd$y, list(c(2L, 3L)))
    expect_equal(a, c1, tolerance = 1e-12)
}
design_threads(old_threads)

# Multithreaded reductions (crossprod(), colSums() and colMeans() with or
# without groups): the rows are split into fixed blocks, one per thread,
# whose sums are added in a fixed order. Results are identical for a given
# number of threads; with one thread they are those of the single-threaded
# kernels; with different numbers of threads they agree up to rounding.
source("helpers.R")
library(design.matrix)

old_threads <- design_threads()
old_opts <- options(design.matrix.dedup = FALSE)

set.seed(42)
n <- 40000
d <- data.frame(
    y = rnorm(n), x = rnorm(n), z = runif(n),
    f1 = factor(sample(sprintf("a%02d", 1:12), n, TRUE)),
    f2 = factor(sample(letters[1:5], n, TRUE)),
    g = sample(1:40, n, TRUE))
w <- runif(n)
wna <- w
wna[c(7, 20000, 39999)] <- NA
Y <- cbind(d$y, d$z)
Yna <- Y
Yna[c(11, 30001), 2] <- NA

# every reduction, on the structured and the row kernels
reductions <- function() {
    out <- list()
    for (st in c(TRUE, FALSE)) {
        options(design.matrix.structured = st)
        for (fo in list(~ x * f1 + f2, ~ f1 * f2, ~ x + z + x:z, ~ z:f1:f2 + x)) {
            D <- design_matrix(fo, d)
            k <- paste(if (st) "structured" else "row", deparse(fo))
            out[[paste(k, "X'X")]] <- crossprod(D)
            out[[paste(k, "X'WX")]] <- crossprod(w * D)
            out[[paste(k, "X'WX NA")]] <- crossprod(wna * D)
            out[[paste(k, "X'Y")]] <- crossprod(D, Y)
            out[[paste(k, "X'Y NA")]] <- crossprod(D, Yna)
            out[[paste(k, "X'Z")]] <- crossprod(D, design_matrix(~ f2 * z, d))
            out[[paste(k, "J'J")]] <- crossprod(w * D - D)
            out[[paste(k, "colSums")]] <- colSums(D)
            out[[paste(k, "colMeans")]] <- colMeans(w * D)
            out[[paste(k, "colSums group")]] <- colSums(D, group = d$g)
            out[[paste(k, "colMeans group")]] <- colMeans(w * D, group = d$g)
            out[[paste(k, "%*%")]] <- D %*% seq_len(ncol(D))
        }
    }
    options(design.matrix.structured = TRUE)
    out
}

design_threads(1)
r1 <- reductions()
r1b <- reductions()
for (k in names(r1)) expect_true(strict_identical(r1[[k]], r1b[[k]]), info = paste("1 thread, repeated:", k))

# one thread: the values of the single-threaded kernels before threading,
# stored as hexadecimal strings by sprintf("%a") (local runs only: the
# dense parts use BLAS, whose sums depend on the library). The file holds
# r1 (without %*%) as computed by the package before multithreaded
# reductions (main at 369557e), with this file's code up to `r1 <- ...`:
# lapply(r1, function(v) sprintf("%a", as.vector(v))).
ref_file <- "threads_ref.rds"
if (at_home() && file.exists(ref_file)) {
    ref <- readRDS(ref_file)
    for (k in names(ref)) {
        v <- r1[[k]]
        expect_true(identical(sprintf("%a", as.vector(v)), ref[[k]]), info = paste("1 thread, stored:", k))
    }
}

# several threads (when OpenMP is available)
if (.Call(design.matrix:::C_dm_omp_max_threads) > 1L) {
    nmax <- min(4L, .Call(design.matrix:::C_dm_omp_max_threads))
    for (nt in unique(c(2L, nmax))) {
        design_threads(nt)
        ra <- reductions()
        rb <- reductions()
        for (k in names(r1)) {
            expect_true(strict_identical(ra[[k]], rb[[k]]), info = paste(nt, "threads, repeated:", k))
            expect_equal(ra[[k]], r1[[k]], tolerance = 1e-10, info = paste(nt, "threads vs 1:", k))
        }
        # the threaded path ran: some sums differ in the last digits
        expect_true(any(!mapply(strict_identical, ra, r1)), info = paste(nt, "threads: threaded"))
    }
    # small problems run on one thread
    design_threads(nmax)
    ds <- d[1:2000, ]
    Ds <- design_matrix(~ x * f1 + f2, ds)
    s4 <- list(crossprod(Ds), colSums(Ds, group = ds$g), crossprod(Ds, ds$y))
    design_threads(1)
    s1 <- list(crossprod(Ds), colSums(Ds, group = ds$g), crossprod(Ds, ds$y))
    expect_true(strict_identical(s4, s1))
}

design_threads(old_threads)
options(old_opts)

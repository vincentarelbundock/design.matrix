## One method of one benchmark cell per process:
##   Rscript cell.R <design> <n> <iters> <mode> <outdir> <method index>
## mode "formula": every method starts from the formula and data; mode
## "prebuilt": matrices are built before timing and only the fit is timed.
## Fairness: one thread everywhere (OPENBLAS_NUM_THREADS/OMP_NUM_THREADS=1
## are set by the caller); design.matrix's na.omit-copy shortcut is off, so
## every model.frame() user pays the same preparation; each result is
## checked against lm(), and a method with a different answer is flagged.
## Memory: peak resident set size above the level before the call (all
## allocations, not just R's heap).
args <- commandArgs(TRUE)
design <- args[1]; n <- as.numeric(args[2]); iters <- as.integer(args[3])
mode <- args[4]; outdir <- args[5]; mi <- as.integer(args[6])
suppressWarnings(suppressMessages({library(design.matrix); library(fixest); library(speedglm)}))
stopifnot(Sys.getenv("OPENBLAS_NUM_THREADS") == "1", Sys.getenv("OMP_NUM_THREADS") == "1")
design_threads(1); setFixest_nthreads(1); data.table::setDTthreads(1)
options(design.matrix.dedup = FALSE, design.matrix.skip_na_copy = FALSE)
set.seed(1)
d <- data.frame(y = rnorm(n), x = rnorm(n), z = rnorm(n), v = rnorm(n),
                f1 = factor(sample(sprintf("a%02d", 1:50), n, TRUE)),
                f2 = factor(sample(letters[1:20], n, TRUE)))
d$f3 <- factor(sample(1:8, n, TRUE)); d$f3[d$f1 %in% c("a01", "a02", "a03")] <- "1"
fo <- switch(design,
  p3 = y ~ x + z, p5 = y ~ x + z + v + I(x^2), p21 = y ~ x + f2, p40 = y ~ x * f2,
  p119 = y ~ x * f1 + f2, p401 = y ~ x + f1 * f3)
y <- d$y
lmfit_coef <- function(X, y) { r <- .lm.fit(X, y); b <- rep(NA_real_, ncol(X)); k <- seq_len(r$rank); b[r$pivot[k]] <- r$coefficients[k]; b }
if (mode == "formula") {
  methods <- list(
    "lm()" = function() coef(lm(fo, d)),
    ".lm.fit()" = function() lmfit_coef(model.matrix(fo, d), y),
    "RcppEigen QR" = function() RcppEigen::fastLmPure(model.matrix(fo, d), y, method = 0L)$coefficients,
    "RcppEigen Cholesky" = function() RcppEigen::fastLmPure(model.matrix(fo, d), y, method = 2L)$coefficients,
    "RcppArmadillo" = function() drop(RcppArmadillo::fastLmPure(model.matrix(fo, d), y)$coefficients),
    "speedglm::speedlm()" = function() coef(speedlm(fo, d)),
    "MatrixModels::glm4()" = function() MatrixModels::glm4(fo, data = d, sparse = TRUE)@pred@coef,
    "fixest::feols()" = function() coef(feols(fo, d, notes = FALSE, warn = FALSE)),
    "design.matrix qr()" = function() qr.coef(qr(design_matrix(fo, d)), y))
} else {
  ## build only the matrix this method uses, before timing
  if (mi %in% c(1:5, 7)) X <- model.matrix(fo, d)
  if (mi == 6) SX <- Matrix::sparse.model.matrix(fo, d)
  if (mi == 8) D <- design_matrix(fo, d)
  methods <- list(
    ".lm.fit()" = function() lmfit_coef(X, y),
    "RcppEigen QR" = function() RcppEigen::fastLmPure(X, y, method = 0L)$coefficients,
    "RcppEigen Cholesky" = function() RcppEigen::fastLmPure(X, y, method = 2L)$coefficients,
    "RcppArmadillo" = function() drop(RcppArmadillo::fastLmPure(X, y)$coefficients),
    "speedglm::speedlm.fit()" = function() coef(speedlm.fit(y, X, intercept = FALSE)),
    "MatrixModels lm.fit.sparse" = function() MatrixModels:::lm.fit.sparse(SX, y),
    "fixest::feols.fit()" = function() coef(feols.fit(y, X, notes = FALSE, warn = FALSE)),
    "design.matrix qr()" = function() qr.coef(qr(D), y))
}
ref <- readRDS(file.path(outdir, sprintf("ref_%s_%g.rds", design, n)))
## Peak resident memory of the process while f() runs, in MB: every
## allocation counts, whether R, C or C++ made it. The kernel's peak is
## reset first (Linux /proc/self/clear_refs), after a gc().
status_kb <- function(key) {
  l <- grep(paste0("^", key, ":"), readLines("/proc/self/status"), value = TRUE)
  as.numeric(gsub("[^0-9]", "", l))
}
peak <- function(f) {
  invisible(gc()); invisible(gc())
  writeLines("5", "/proc/self/clear_refs")
  base <- status_kb("VmRSS")
  f()
  (status_kb("VmHWM") - base) / 1024
}
## Fresh process: the first call gives the peak memory (nothing freed by an
## earlier method is still resident), then the repeats are timed.
m <- names(methods)[mi]
f <- methods[[m]]
res <- NULL
mb <- tryCatch(peak(function() res <<- suppressWarnings(f())), error = function(e) e)
row <- if (inherits(mb, "error")) {
  data.frame(mode, design, n, p = length(ref), method = m, ms = NA, MB = NA,
             matches_lm = NA, note = substr(conditionMessage(mb), 1, 60))
} else {
  b <- unname(res)
  same <- length(b) == length(ref) && identical(is.na(b), unname(is.na(ref))) &&
    isTRUE(all(abs(b[!is.na(b)] - ref[!is.na(ref)]) <= 1e-8 * pmax(1, abs(ref[!is.na(ref)]))))
  t <- suppressWarnings(bench::mark(f(), iterations = iters, check = FALSE, memory = FALSE, filter_gc = FALSE))
  data.frame(mode, design, n, p = length(ref), method = m, ms = round(1000 * as.numeric(t$median), 1),
             MB = round(mb), matches_lm = same, note = "")
}
write.csv(row, file.path(outdir, sprintf("%s_%s_%g_%02d.csv", mode, design, n, mi)), row.names = FALSE)

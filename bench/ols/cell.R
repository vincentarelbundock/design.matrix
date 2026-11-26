## One method on one design, in a fresh process:
##   Rscript cell.R <family> <width> <n> <iters> <mode> <outdir> <method index>
## mode "formula": every method starts from the formula and the data;
## mode "prebuilt": the matrix is built before timing, only the fit is timed.
## Fairness: one thread everywhere (the caller sets OPENBLAS_NUM_THREADS and
## OMP_NUM_THREADS to 1 and pins the process to a core); xtx's
## na.omit-copy shortcut off; the result is checked against lm().
## Memory: peak resident memory of the first call (every allocation, R or
## C++), in a fresh process. Time: gc() before each timed call, so that
## every call starts from the same heap; 3 to `iters` calls (up to about 3
## seconds per method); the median is reported, with the minimum and maximum.
args <- commandArgs(TRUE)
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(here, "designs.R"))
family <- args[1]; width <- as.integer(args[2]); n <- as.numeric(args[3])
iters <- as.integer(args[4]); mode <- args[5]; outdir <- args[6]; mi <- as.integer(args[7])
suppressWarnings(suppressMessages({library(xtx); library(fixest); library(speedglm)}))
stopifnot(Sys.getenv("OPENBLAS_NUM_THREADS") == "1", Sys.getenv("OMP_NUM_THREADS") == "1")
design_threads(1); setFixest_nthreads(1); data.table::setDTthreads(1)
options(xtx.dedup = FALSE, xtx.skip_na_copy = FALSE)
dd <- make_data(family, width, n)
d <- dd$data; fo <- dd$formula; y <- d$y
lmfit_coef <- function(X, y) {
    r <- .lm.fit(X, y); b <- rep(NA_real_, ncol(X)); k <- seq_len(r$rank)
    b[r$pivot[k]] <- r$coefficients[k]; b
}
method_names <- if (mode == "formula") {
    c("lm()", ".lm.fit()", "RcppEigen QR", "RcppEigen Cholesky", "RcppArmadillo",
      "speedglm", "MatrixModels", "fixest", "xtx qr()")
} else {
    c(".lm.fit()", "RcppEigen QR", "RcppEigen Cholesky", "RcppArmadillo",
      "speedglm", "MatrixModels", "fixest", "xtx qr()")
}
m <- method_names[mi]
if (mode == "formula") {
    f <- switch(m,
        "lm()" = function() coef(lm(fo, d)),
        ".lm.fit()" = function() lmfit_coef(model.matrix(fo, d), y),
        "RcppEigen QR" = function() RcppEigen::fastLmPure(model.matrix(fo, d), y, method = 0L)$coefficients,
        "RcppEigen Cholesky" = function() RcppEigen::fastLmPure(model.matrix(fo, d), y, method = 2L)$coefficients,
        "RcppArmadillo" = function() drop(RcppArmadillo::fastLmPure(model.matrix(fo, d), y)$coefficients),
        "speedglm" = function() coef(speedlm(fo, d)),
        "MatrixModels" = function() MatrixModels::glm4(fo, data = d, sparse = TRUE)@pred@coef,
        "fixest" = function() coef(feols(fo, d, notes = FALSE, warn = FALSE)),
        "xtx qr()" = function() qr.coef(qr(design_matrix(fo, d)), y))
} else {
    ## build only the matrix this method uses, before timing
    if (m == "MatrixModels") SX <- Matrix::sparse.model.matrix(fo, d)
    else if (m == "xtx qr()") D <- design_matrix(fo, d)
    else X <- model.matrix(fo, d)
    f <- switch(m,
        ".lm.fit()" = function() lmfit_coef(X, y),
        "RcppEigen QR" = function() RcppEigen::fastLmPure(X, y, method = 0L)$coefficients,
        "RcppEigen Cholesky" = function() RcppEigen::fastLmPure(X, y, method = 2L)$coefficients,
        "RcppArmadillo" = function() drop(RcppArmadillo::fastLmPure(X, y)$coefficients),
        "speedglm" = function() coef(speedlm.fit(y, X, intercept = FALSE)),
        "MatrixModels" = function() MatrixModels:::lm.fit.sparse(SX, y),
        "fixest" = function() coef(feols.fit(y, X, notes = FALSE, warn = FALSE)),
        "xtx qr()" = function() qr.coef(qr(D), y))
}
ref <- readRDS(file.path(outdir, sprintf("ref_%s_%d.rds", family, width)))
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
res <- NULL
mb <- tryCatch(peak(function() res <<- suppressWarnings(f())), error = function(e) e)
row <- data.frame(mode, family, width, n, p = length(ref), method = m, runs = NA)
if (inherits(mb, "error")) {
    row <- cbind(row, ms = NA, ms_min = NA, ms_max = NA, MB = NA, matches_lm = NA,
                 note = substr(conditionMessage(mb), 1, 60))
} else {
    b <- unname(res)
    ok <- !is.na(ref)
    same <- length(b) == length(ref) && identical(is.na(b), unname(is.na(ref))) &&
        isTRUE(all(abs(b[ok] - ref[ok]) <= 1e-8 * pmax(1, abs(ref[ok]))))
    ## at least 3 and at most `iters` timed calls, stopping once the method
    ## has used about 3 seconds (slow methods vary little relative to their
    ## time)
    tm <- numeric(0)
    while (length(tm) < iters && (length(tm) < 3 || sum(tm) < 3000)) {
        invisible(gc())
        t0 <- bench::hires_time()
        invisible(suppressWarnings(f()))
        tm <- c(tm, (bench::hires_time() - t0) * 1000)
    }
    row$runs <- length(tm)
    row <- cbind(row, ms = round(median(tm), 2), ms_min = round(min(tm), 2), ms_max = round(max(tm), 2),
                 MB = round(mb, 1), matches_lm = same, note = "")
}
write.csv(row, file.path(outdir, sprintf("%s_%s_%d_%02d.csv", mode, family, width, mi)), row.names = FALSE)

# Each expression runs once before it is timed (bench measures memory on the
# first timed run), so that a plan is prepared, as in repeated use.
warm_mark <- function(...) {
    exprs <- eval(substitute(alist(...)))
    for (e in exprs) invisible(eval(e, parent.frame()))
    bench::mark(..., iterations = iters(5), check = FALSE)
}
w <- rnorm(n)
sfo <- y ~ x + f1 * f2
marks <- list(
    "Predictions on new data" = warm_mark(
        predict(fit_big, big),
        drop(design_matrix(spec_big, big) %*% coef(fit_big))),
    "Standard errors of predictions" = warm_mark(
        predict(fit_big, big, se.fit = TRUE)$se.fit,
        sqrt(quadratic_diagonal(design_matrix(spec_big, big), vcov(fit_big)))),
    "Grouped weighted sums (cluster-robust variance)" = warm_mark(
        rowsum(w * model.matrix(fo, big), big$cl),
        colSums(w * design_matrix(fo, big), group = big$cl)),
    "Sparse model matrix (vs. Matrix::sparse.model.matrix())" = warm_mark(
        Matrix::sparse.model.matrix(sfo, big),
        as(design_matrix(sfo, big), "CsparseMatrix")),
    "Least squares coefficients (vs. lm())" = warm_mark(
        coef(lm(fo, big)),
        qr.coef(qr(design_matrix(fo, big)), big$y)))
cell <- function(m, i) paste0(format(m$median[i]), ", ", format(m$mem_alloc[i]))
data.frame(
    Task = names(marks),
    "Base R" = sapply(marks, cell, i = 1),
    "design.matrix" = sapply(marks, cell, i = 2),
    check.names = FALSE)

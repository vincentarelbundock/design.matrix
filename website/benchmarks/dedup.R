dup <- data.frame(
    f1 = factor(sample.int(50, n, TRUE)),
    f3 = factor(sample.int(5, n, TRUE)),
    f2 = factor(sample.int(20, n, TRUE)))
pd <- design_matrix(~ f1 * f3 + f2, dup, contrasts.arg = list(f1 = "contr.sum", f3 = "contr.sum"))
Vd <- diag(ncol(design_matrix(pd, dup[1:2, ])))

quad <- function(mode, D = design_matrix(pd, dup)) {
    op <- options(xtx.dedup = mode)
    on.exit(options(op))
    quadratic_diagonal(D, Vd)
}
D_cached <- design_matrix(pd, dup)
invisible(quad("auto", D_cached))
bench::mark(
    "every row"              = quad(FALSE),
    "grouped, first call"    = quad("auto"),
    "grouped, groups cached" = quad("auto", D_cached),
    iterations = iters(5), check = same)

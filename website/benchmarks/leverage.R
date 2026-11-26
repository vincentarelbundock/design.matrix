XtXi <- vcov(fit_big) / sigma(fit_big)^2
base_leverage <- function(X, XtXi) {
    rowSums((X %*% XtXi) * X)
}
fast_leverage <- function(D, XtXi) {
    quadratic_diagonal(D, XtXi)
}

bench::mark(
    hatvalues(fit_big),
    base_leverage(X, XtXi),
    fast_leverage(DX, XtXi),
    iterations = iters(5), check = close)

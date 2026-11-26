XtXi <- vcov(fit_big) / sigma(fit_big)^2
base_cluster_vcov <- function(X, e, cluster, bread) {
    S <- rowsum(e * X, cluster)
    bread %*% crossprod(S) %*% bread
}
fast_cluster_vcov <- function(D, e, cluster, bread) {
    S <- colSums(e * D, group = cluster)
    bread %*% crossprod(S) %*% bread
}

e <- residuals(fit_big)
bench::mark(
    base_cluster_vcov(X, e, big$cl, XtXi),
    fast_cluster_vcov(DX, e, big$cl, XtXi),
    iterations = 5, check = close)

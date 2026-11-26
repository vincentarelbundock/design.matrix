V <- crossprod(matrix(rnorm(ncol(X)^2), ncol(X))) / ncol(X)
base_diagonal <- function(X, V) rowSums((X %*% V) * X)

bench::mark(
    base_diagonal(model.matrix(fo, big), V),
    quadratic_diagonal(design_matrix(fo, big), V),
    base_diagonal(X, V),
    quadratic_diagonal(DX, V),
    iterations = iters(5), check = close)

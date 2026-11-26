beta <- rnorm(ncol(X))
bench::mark(
    "model.matrix, then %*%"  = model.matrix(fo, big) %*% beta,
    "design_matrix, then %*%" = design_matrix(fo, big) %*% beta,
    "prebuilt X"              = X %*% beta,
    "prebuilt design matrix"  = DX %*% beta,
    iterations = iters(10), check = close)

w <- rnorm(n)
bench::mark(
    "model.matrix, then crossprod"  = crossprod(model.matrix(fo, big), w),
    "design_matrix, then crossprod" = crossprod(design_matrix(fo, big), w),
    "prebuilt X"                    = crossprod(X, w),
    "prebuilt design matrix"        = crossprod(DX, w),
    iterations = iters(10), check = close)

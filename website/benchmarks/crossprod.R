bench::mark(
    "model.matrix, then crossprod"  = crossprod(model.matrix(fo, big)),
    "design_matrix, then crossprod" = crossprod(design_matrix(fo, big)),
    "prebuilt X"                    = crossprod(X),
    "prebuilt design matrix"        = crossprod(DX),
    iterations = iters(10), check = close)

w <- rnorm(n)
bench::mark(
    "model.matrix, then rowsum"   = rowsum(w * model.matrix(fo, big), big$f2),
    "design_matrix, then colSums" = colSums(w * design_matrix(fo, big), group = big$f2),
    "prebuilt X"                  = rowsum(w * X, big$f2),
    "prebuilt design matrix"      = colSums(w * DX, group = big$f2),
    iterations = 10, check = close)

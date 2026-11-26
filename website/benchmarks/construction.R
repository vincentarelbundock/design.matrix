bench::mark(
    "stats::model.matrix" = model.matrix(fo, big),
    "design_matrix"       = design_matrix(fo, big),
    iterations = 10, check = FALSE)

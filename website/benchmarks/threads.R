fast_model_matrix <- function(formula, data, threads) {
    design_threads(threads)
    as.matrix(design_matrix(formula, data))
}

bench::mark(
    model.matrix(fo, big),
    fast_model_matrix(fo, big, threads = 1),
    fast_model_matrix(fo, big, threads = 4),
    iterations = 10, check = same)

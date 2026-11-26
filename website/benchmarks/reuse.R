tt <- delete.response(terms(fit_big))
batch <- big[1:100, ]
base_model_matrix <- function(batch) {
    mf <- model.frame(tt, batch, xlev = fit_big$xlevels, na.action = na.pass)
    model.matrix(tt, mf, contrasts.arg = fit_big$contrasts)
}
fast_model_matrix <- function(batch) {
    as.matrix(design_matrix(spec_big, batch))
}

bench::mark(
    base_model_matrix(batch),
    fast_model_matrix(batch),
    iterations = iters(200), check = same)

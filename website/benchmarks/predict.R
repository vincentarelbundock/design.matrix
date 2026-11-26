fast_predict <- function(fit, newdata, spec) {
    D <- design_matrix(spec, newdata)
    drop(D %*% coef(fit))
}

bench::mark(
    "predict(fit, newdata)" = predict(fit_big, big),
    "fast_predict"          = fast_predict(fit_big, big, spec_big),
    iterations = iters(5), check = close)

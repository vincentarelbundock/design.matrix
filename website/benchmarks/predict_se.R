base_predict_se <- function(fit, newdata) {
    predict(fit, newdata, se.fit = TRUE)$se.fit
}
fast_predict_se <- function(fit, newdata, spec) {
    D <- design_matrix(spec, newdata)
    sqrt(quadratic_diagonal(D, vcov(fit)))
}

bench::mark(
    base_predict_se(fit_big, big),
    fast_predict_se(fit_big, big, spec_big),
    iterations = 5, check = close)

# Median time to predict on new batches of rows, after one-time setup.
set.seed(27)
n_new <- min(100000, n)  # n is smaller with BENCH_SIZE=quick (setup.R)
d <- data.frame(x = rnorm(n_new), z = rnorm(n_new),
                region = factor(sample(1:20, n_new, replace = TRUE)))
train <- data.frame(x = rnorm(5000), z = rnorm(5000),
                    region = factor(sample(1:20, 5000, replace = TRUE)))
train$y <- rnorm(nrow(train))

time_model <- function(formula) {
    fit <- lm(formula, train)
    tt <- delete.response(terms(fit))
    spec <- design_matrix(fit, newdata = NULL)
    beta <- coef(fit)

    base_predict <- function(batch) {
        mf <- model.frame(tt, batch, xlev = fit$xlevels)
        X <- model.matrix(tt, mf, contrasts.arg = fit$contrasts)
        drop(X %*% beta)
    }
    fast_predict <- function(batch) {
        D <- design_matrix(spec, batch)
        drop(D %*% beta)
    }

    time_batch <- function(rows) {
        batch <- head(d, rows)
        timings <- bench::mark(
            base_predict(batch), fast_predict(batch),
            iterations = iters(15), memory = FALSE, filter_gc = FALSE,
            check = function(a, b) isTRUE(all.equal(as.numeric(a), as.numeric(b))))
        data.frame(rows = rows, ms = 1000 * as.numeric(timings$median),
                   method = c("model.matrix()", "design_matrix()"))
    }
    batches <- c(1000, 10000, 25000, 50000, 100000)
    do.call(rbind, lapply(batches[batches <= n_new], time_batch))
}

rbind(
    transform(time_model(y ~ x + z), model = "Two numeric predictors"),
    transform(time_model(y ~ x * region), model = "Numeric × 20-level factor"))

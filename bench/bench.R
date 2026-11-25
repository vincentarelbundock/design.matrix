# Benchmarks against stats. Run with: Rscript bench/bench.R
library(design.matrix)
tm <- function(e, k = 5) { e <- substitute(e); env <- parent.frame(); min(replicate(k, system.time(eval(e, env))[[3]])) }
us <- function(e, k = 2000) { e <- substitute(e); env <- parent.frame(); f <- function() eval(e, env); invisible(f()); 1e6 * system.time(for (i in seq_len(k)) f())[[3]] / k }
set.seed(1)
N <- 1e6
D <- data.frame(y = rnorm(N), x1 = rnorm(N), z = runif(N),
                f1 = factor(sample(sprintf("a%02d", 1:50), N, TRUE)),
                f2 = factor(sample(letters[1:20], N, TRUE)))
res <- list()
row <- function(task, stats, ours) res[[length(res) + 1L]] <<- data.frame(task, stats = stats, design.matrix = ours, speedup = stats / ours)

cat("dense model matrix, n = 1e6 (seconds)\n")
for (fo in list(y ~ x1 * f1 + f2, y ~ x1 * f1 + poly(z, 3))) {
    s <- tm(model.matrix(fo, D), 3)
    for (th in c(1L, 4L, 16L)) { design_threads(th); row(sprintf("dense design %s, %d thr", deparse(fo), th), s, tm(as.matrix(design_matrix(fo, D)), 3)) }
}
design_threads(1L)
d5 <- D[1:1e5, ]
s <- tm(model.matrix(~ f1 * f2, d5), 3)
for (th in c(1L, 16L)) { design_threads(th); row(sprintf("dense design ~f1*f2 (n=1e5, p=1000), %d thr", th), s, tm(as.matrix(design_matrix(~ f1 * f2, d5)), 3)) }
design_threads(1L)

cat("algebra, n = 1e6, y ~ x1 * f1 + f2 (seconds; stats column includes building X)\n")
fo <- y ~ x1 * f1 + f2
X <- model.matrix(fo, D); p <- ncol(X); b <- rnorm(p); w <- rnorm(N)
V <- crossprod(matrix(rnorm(p * p), p)) / 1e4
Dd <- design_matrix(fo, D)
tX <- tm(model.matrix(fo, D), 3)
row("X %*% b", tX + tm(X %*% b), tm(design_matrix(fo, D) %*% b, 3))
row("crossprod(X, w)", tX + tm(crossprod(X, w)), tm(crossprod(design_matrix(fo, D), w), 3))
row("rowSums((X %*% V) * X)", tX + tm(rowSums((X %*% V) * X), 2), tm(quadratic_diagonal(design_matrix(fo, D), V), 3))
design_threads(16L)
row("rowSums((X %*% V) * X), 16 thr", tX + tm(rowSums((X %*% V) * X), 2), tm(quadratic_diagonal(design_matrix(fo, D), V), 3))
design_threads(1L)
P <- design_matrix(fo, D[1:1000, ])
D1 <- transform(D, x1 = x1 + 1)
s0 <- runif(N); s1 <- runif(N)
dense_J <- function() { X0 <- model.matrix(fo, D); X1 <- model.matrix(fo, D1); J <- s1 * X1 - s0 * X0; rowSums((J %*% V) * J) }
ours_J <- function() quadratic_diagonal(s1 * design_matrix(P, D1) - s0 * design_matrix(P, D), V)
row("comparison SEs: diag(J V J')", tm(dense_J(), 2), tm(ours_J(), 3))
rm(X); invisible(gc())

cat("small repeated calls, n = 100 (microseconds per call)\n")
d100 <- D[1:100, ]
fo <- y ~ x1 * f1 + f2
P <- design_matrix(fo, D[1:5000, ])
fit <- lm(fo, D[1:5000, ])
tt <- delete.response(terms(fit))
ref <- function() model.matrix(tt, model.frame(tt, d100, xlev = fit$xlevels, na.action = na.pass), contrasts.arg = fit$contrasts)
row("model.matrix(formula, data)", us(model.matrix(fo, d100)), us(as.matrix(design_matrix(fo, d100))))
row("plan on new data (predict.lm-style)", us(ref()), us(as.matrix(design_matrix(P, d100))))
b <- rnorm(ncol(ref()))
row("plan on new data, then X %*% b", us(ref() %*% b), us(design_matrix(P, d100) %*% b))
stopifnot(identical(as.matrix(design_matrix(P, d100)), ref(), num.eq = FALSE))

out <- do.call(rbind, res)
out$stats <- signif(out$stats, 3); out$design.matrix <- signif(out$design.matrix, 3); out$speedup <- round(out$speedup, 1)
print(out, row.names = FALSE)

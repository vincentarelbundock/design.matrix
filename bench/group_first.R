# Grouping identical rows before building designs: a marginaleffects-style
# average comparison with delta-method SEs. Run with: Rscript bench/group_first.R
suppressMessages({library(xtx); library(marginaleffects)})
design_row_groups <- xtx:::design_row_groups
set.seed(1); N <- 1e6
d <- data.frame(f1 = factor(sample(sprintf("r%02d", 1:20), N, TRUE)), age = sample(18:90, N, TRUE),
                sex = sample(c("m", "f"), N, TRUE), edu = factor(sample(1:5, N, TRUE)))
d$y <- rbinom(N, 1, plogis(-1 + 0.02 * (d$age - 50) + 0.3 * (d$sex == "m") + 0.1 * as.integer(d$edu)))
fo <- y ~ f1 + age * sex + edu
fit <- glm(fo, data = d, family = binomial)
b <- coef(fit); V <- vcov(fit)
P <- design_matrix(fo, d[1:5000, ])
el <- function(f) { t <- system.time(r <- f())[[3]]; list(t = t, r = r) }

workflow <- function(data, w = NULL) {        # w: weights per row (counts) or NULL
  X0 <- design_matrix(P, data); X1 <- design_matrix(P, transform(data, age = age + 1))
  e0 <- drop(X0 %*% b); e1 <- drop(X1 %*% b)
  ww <- if (is.null(w)) rep(1 / nrow(data), nrow(data)) else w / sum(w)
  est <- sum(ww * (plogis(e1) - plogis(e0)))
  J <- dlogis(e1) * X1 - dlogis(e0) * X0
  jbar <- drop(crossprod(J, ww))
  se <- sqrt(drop(t(jbar) %*% V %*% jbar))
  unit_se <- sqrt(quadratic_diagonal(J, V))
  list(est = est, se = se, unit_se = unit_se)
}
options(xtx.dedup = FALSE); A <- el(function() workflow(d))
options(xtx.dedup = "auto"); B <- el(function() workflow(d))
C <- el(function() {
  g <- design_row_groups(design_matrix(P, d))                 # group the data once
  r <- workflow(d[g$first, ], w = g$counts)               # everything on distinct rows
  r$unit_se <- r$unit_se[g$map]                           # expand unit-level results
  r$U <- length(g$first); r })
M <- el(function() avg_comparisons(fit, variables = list(age = 1), vcov = V))
cat(sprintf("distinct rows: %d of %d (%.2f%%)\n\n", C$r$U, N, 100 * C$r$U / N))
cat(sprintf("%-38s %8s\n", "strategy", "seconds"))
cat(sprintf("%-38s %8.3f\n", "every row (dedup off)", A$t))
cat(sprintf("%-38s %8.3f\n", "groups inside kernels (auto)", B$t))
cat(sprintf("%-38s %8.3f\n", "group the data first", C$t))
cat(sprintf("%-38s %8.3f   (avg only, no unit SEs)\n", "marginaleffects::avg_comparisons", M$t))
cat("\nestimate, A vs C:", format(A$r$est, digits = 15), format(C$r$est, digits = 15), " rel diff", signif(abs(A$r$est - C$r$est) / abs(A$r$est), 2), "\n")
cat("std err,  A vs C:", format(A$r$se, digits = 15), format(C$r$se, digits = 15), "\n")
cat("unit SEs identical bit for bit, A vs B:", identical(A$r$unit_se, B$r$unit_se, num.eq = FALSE),
    " A vs C:", identical(A$r$unit_se, C$r$unit_se, num.eq = FALSE), "\n")
cat("marginaleffects estimate / se:", format(M$r$estimate, digits = 10), format(M$r$std.error, digits = 6), "\n")

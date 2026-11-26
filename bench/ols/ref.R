## Reference coefficients from lm(), computed in their own process.
args <- commandArgs(TRUE); design <- args[1]; n <- as.numeric(args[2]); outdir <- args[3]
set.seed(1)
d <- data.frame(y = rnorm(n), x = rnorm(n), z = rnorm(n), v = rnorm(n),
                f1 = factor(sample(sprintf("a%02d", 1:50), n, TRUE)),
                f2 = factor(sample(letters[1:20], n, TRUE)))
d$f3 <- factor(sample(1:8, n, TRUE)); d$f3[d$f1 %in% c("a01", "a02", "a03")] <- "1"
fo <- switch(design,
  p3 = y ~ x + z, p5 = y ~ x + z + v + I(x^2), p21 = y ~ x + f2, p40 = y ~ x * f2,
  p119 = y ~ x * f1 + f2, p401 = y ~ x + f1 * f3)
saveRDS(coef(lm(fo, d)), file.path(outdir, sprintf("ref_%s_%g.rds", design, n)))

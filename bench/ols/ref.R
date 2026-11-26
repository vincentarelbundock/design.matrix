## Reference coefficients from lm(), in their own process:
##   Rscript ref.R <family> <width> <n> <outdir>
args <- commandArgs(TRUE)
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(here, "designs.R"))
family <- args[1]; width <- as.integer(args[2]); n <- as.numeric(args[3]); outdir <- args[4]
dd <- make_data(family, width, n)
saveRDS(coef(lm(dd$formula, dd$data)), file.path(outdir, sprintf("ref_%s_%d.rds", family, width)))

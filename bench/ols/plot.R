## Figures of the OLS benchmark (bench/ols/results.csv), one per family of
## designs (the rank-deficient design is a table in the README): rows are
## time and memory, columns the two comparisons, x the number of columns. Every method is named in the legend; filled points
## match lm(), hollow points give a different answer. Run from the package
## root: Rscript bench/ols/plot.R
library(ggplot2)

ols_results <- function(file = "bench/ols/results.csv") {
    res <- read.csv(file)
    res$comparison <- factor(ifelse(res$mode == "formula", "From formula and data", "Prebuilt matrix (fit only)"),
                             levels = c("From formula and data", "Prebuilt matrix (fit only)"))
    res$Method <- factor(res$method, levels = names(ols_colours))
    res
}

## categorical palette in a fixed order (validated for colour-vision
## deficiency); xtx in near-black, drawn thicker
ols_colours <- c("xtx qr()" = "#191919", "lm()" = "#257de0", ".lm.fit()" = "#f47230",
                 "RcppEigen QR" = "#00b894", "RcppEigen Cholesky" = "#f2a600",
                 "RcppArmadillo" = "#ee7fac", "speedglm" = "#008a17",
                 "MatrixModels" = "#6a4ce2", "fixest" = "#ec4647")

ols_families <- list(
    numeric = "Numeric predictors only: y ~ x1 + ... + xk",
    mixed = "Numeric predictors and a 20-level factor: y ~ x1 + ... + xk + f",
    factor = "One numeric predictor and a factor: y ~ x + f",
    interaction = "A numeric predictor interacting with a factor: y ~ x * f")

ols_plot <- function(res, family, title) {
    x <- res[res$family == family & !is.na(res$ms), ]
    long <- rbind(
        transform(x, metric = "Median time (ms)", value = ms),
        transform(x, metric = "Peak memory (MB)", value = pmax(MB, 1)))
    long$metric <- factor(long$metric, levels = c("Median time (ms)", "Peak memory (MB)"))
    long$answer <- factor(ifelse(long$matches_lm, "matches lm()", "different answer"),
                          levels = c("matches lm()", "different answer"))
    ## a shape legend only when some method gives a different answer
    all_match <- all(long$matches_lm %in% TRUE)
    ours <- long[long$Method == "xtx qr()", ]
    others <- long[long$Method != "xtx qr()", ]
    ggplot(mapping = aes(p, value, colour = Method)) +
        geom_line(data = others, linewidth = 0.6) +
        geom_point(data = others, aes(shape = answer), size = 1.8, stroke = 0.7) +
        geom_line(data = ours, linewidth = 1.3) +
        geom_point(data = ours, aes(shape = answer), size = 2.6, stroke = 1) +
        facet_grid(metric ~ comparison, scales = "free_y", switch = "y") +
        scale_x_log10(breaks = sort(unique(x$p))) +
        scale_y_log10(labels = scales::label_comma()) +
        scale_colour_manual(values = ols_colours, name = NULL, drop = FALSE) +
        scale_shape_manual(values = c("matches lm()" = 16, "different answer" = 1), name = NULL) +
        guides(colour = guide_legend(order = 1, ncol = 1, override.aes = list(shape = 16)),
               shape = if (all_match) "none" else
                   guide_legend(order = 2, ncol = 1, override.aes = list(colour = "#444444", size = 2.6, stroke = 1))) +
        labs(x = "Columns in the model matrix (log scale)", y = NULL, title = title,
             subtitle = paste0("200,000 rows; one thread everywhere; memory below 1 MB is drawn at 1 MB.",
                               if (all_match) "\nEvery method matches lm() on these designs." else "")) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "right", strip.placement = "outside", panel.grid.minor = element_blank(),
              plot.subtitle = element_text(colour = "#666666", size = 10))
}

ols_figures <- function(res = ols_results(), dir = "bench/ols") {
    specs <- ols_families
    for (fam in names(specs)) {
        ggsave(file.path(dir, sprintf("fig-%s.png", fam)), ols_plot(res, fam, specs[[fam]]),
               width = 10.5, height = 6.5, dpi = 150, bg = "white")
    }
    invisible(NULL)
}

if (sys.nframe() == 0L) ols_figures()

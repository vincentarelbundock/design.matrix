#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: [Fixed effects])
#metadata((title: "Fixed effects")) <website-metadata>

#import "/benchmarks/bench.typ": bench-code

#calepin.setup(echo: false, eval: true, results: "verbatim", warning: false, message: false)

// Building this page never runs the benchmark: the chunks read
// results/fixed_effects.csv, written by `Rscript website/benchmarks/run.R fixed_effects`.
#calepin.chunk(echo: false, results: "hide")[
```r
options(width = 200)
library(ggplot2)
res <- read.csv("results/fixed_effects.csv")
machine <- readLines("results/fixed_effects_machine.txt")
res$fe <- factor(ifelse(res$nfe == 2, "indiv + year", "indiv + year + firm"),
                 c("indiv + year", "indiv + year + firm"))
res$regressors <- factor(paste(res$k, ifelse(res$k == 1, "regressor", "regressors")),
                         paste(sort(unique(res$k)), ifelse(sort(unique(res$k)) == 1, "regressor", "regressors")))
pkgs <- c("design_lm" = "#18756d", "design_lm (DuckDB)" = "#7fbfb3",
          "feols" = "#ad481e", "felm" = "#6b5ca5", "dbreg (DuckDB)" = "#c79a2c")
res$package <- factor(res$package, names(pkgs))
time_plot <- function(type) {
    x <- res[res$dgp == type & !is.na(res$seconds), ]
    ggplot(x, aes(n, seconds, colour = package)) +
        geom_line(linewidth = 0.8) + geom_point(size = 2) +
        facet_grid(regressors ~ fe) +
        scale_x_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
        scale_y_log10(labels = scales::label_number(drop0trailing = TRUE)) +
        scale_colour_manual(values = pkgs, drop = FALSE) +
        labs(x = "Observations (log scale)", y = "Median seconds (log scale)", colour = NULL) +
        theme_minimal(base_size = 13) + theme(legend.position = "bottom")
}
provisional <- any(res$provisional)
```
]

#title()

How fast is `design_lm()` with fixed effects, compared with the other R packages that fit linear models with high-dimensional fixed effects? This page replicates the fixed-effects OLS benchmark of #link("https://github.com/py-econometrics/pyfixest")[pyfixest], itself adapted from Kyle Butts' #link("https://github.com/kylebutts/fixest_benchmarks")[fixest_benchmarks], for R packages only.

#calepin.chunk(echo: false, results: "typst")[
```r
if (provisional) cat("*Provisional results.* These numbers come from a reduced grid (", paste(format(sort(unique(res$n)), big.mark = ",", scientific = FALSE, trim = TRUE), collapse = ", "), " observations; ", paste(sort(unique(res$k)), collapse = " and "), " regressors), run while validating the scripts. The full grid, up to 10 million observations, will replace them.\n", sep = "")
```
]

= Data and models

The data are a panel of individuals observed over 10 years, with about 23 individuals per firm: for $n$ observations, $n slash 10$ individuals and $n slash 230$ firms. The regressors $x_1, dots, x_k$ are independent standard normal, with coefficients $beta_j = 1 slash j$. Individual, firm and year effects are standard normal, as is the error:
$ y = sum_(j=1)^k x_j beta_j + alpha_"indiv" + gamma_"firm" + delta_"year" + epsilon. $

The two designs differ only in how individuals are assigned to firms:

- *Simple*: each observation's firm is drawn at random.
- *Difficult*: firms are assigned round robin, so every individual changes firm every year. The firm and individual effects are then poorly separated, which slows down the alternating projections that demeaning methods use.

Each design is fitted with two sets of fixed effects, `y ~ x1 + ... + xk | indiv_id + year` and `y ~ x1 + ... + xk | indiv_id + year + firm_id`, with $k$ = 1, 5 or 10 regressors and iid standard errors, on 1,000 to 10 million observations.

The packages:

- `design_lm()`: on a data frame, and on the same data in a DuckDB table (`design_lm(f, con, table = "d")`).
- `fixest::feols()`, with `vcov = "iid"`.
- `lfe::felm()`.
- `dbreg::dbreg()` on the DuckDB table, with its automatic strategy, which here demeans the data in SQL: in a single pass for two fixed effects (the panel is balanced in individuals and years), by alternating projections for three. Its default limit of 100 alternating projections stops with an error on the difficult design with three fixed effects, so the benchmark raises it to 10,000.

`alpaca` is not included (not installed); it fits generalized linear models, and `feglm()`'s Gaussian family would go through the same demeaning as `feols()`.

= Protocol

- Everything runs on one core and one thread: `design_threads(1)`, `setFixest_nthreads(1)`, `data.table::setDTthreads(1)`, `options(lfe.threads = 1)`, DuckDB's `SET threads = 1`, and a single-threaded BLAS, with the whole run pinned to one CPU core.
- Each combination of package, design, size, fixed effects and number of regressors runs in a fresh R process: simulate the data (same seed for every package), load the data into DuckDB for the remote fits, fit once (warmup), then time three fits and keep the median.
- Memory is the peak resident memory of the process during the warmup fit, minus the resident memory just before it. It counts every allocation, in R, C, C++ or DuckDB, but not the data itself.
- Once a package takes more than 120 seconds at some size, it is not run at larger sizes: those points are missing from the plots. A fit that takes longer than 120 seconds is timed once, not three times.
- Every result is compared with the coefficients of `feols()` run with a tight convergence tolerance (`fixef.tol = 1e-11`); the table at the end gives the largest absolute difference.

#calepin.chunk(echo: false, results: "verbatim")[
```r
cat(machine, sep = "\n")
```
]

= Simple design

#calepin.chunk(
  echo: false,
  results: "render",
  fig-width: "100%",
  fig-device-width: 10,
  fig-device-height: 8,
  fig-alt-text: "Median fitting time by number of observations on the simple design, log-log scale, for each package, by number of regressors (rows) and fixed effects (columns).",
)[
```r
time_plot("simple")
```
]

= Difficult design

#calepin.chunk(
  echo: false,
  results: "render",
  fig-width: "100%",
  fig-device-width: 10,
  fig-device-height: 8,
  fig-alt-text: "Median fitting time by number of observations on the difficult design, log-log scale, for each package, by number of regressors (rows) and fixed effects (columns).",
)[
```r
time_plot("difficult")
```
]

= Memory

Extra peak memory of one fit, with the most regressors, by number of observations.

#calepin.chunk(
  echo: false,
  results: "render",
  fig-width: "100%",
  fig-device-width: 10,
  fig-device-height: 6,
  fig-alt-text: "Extra peak resident memory of one fit by number of observations, log-log scale, for each package, by design (rows) and fixed effects (columns), with the largest number of regressors.",
)[
```r
x <- res[res$k == max(res$k) & !is.na(res$mem_mb), ]
x$dgp <- factor(x$dgp, c("simple", "difficult"))
ggplot(x, aes(n, pmax(mem_mb, 0.1), colour = package)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    facet_grid(dgp ~ fe) +
    scale_x_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    scale_y_log10(labels = scales::label_number(drop0trailing = TRUE)) +
    scale_colour_manual(values = pkgs, drop = FALSE) +
    labs(x = "Observations (log scale)", y = "Extra peak memory, MB (log scale)", colour = NULL) +
    theme_minimal(base_size = 13) + theme(legend.position = "bottom")
```
]

= Accuracy and missing points

Largest absolute difference between each package's coefficients and those of `feols()` at a tight tolerance, over every fit; and the fits that were not run (more than 120 seconds at a smaller size) or did not finish.

#calepin.chunk(echo: false, results: "typst")[
```r
acc <- aggregate(max_coef_diff ~ package, res, max)
acc$max_coef_diff <- format(acc$max_coef_diff, digits = 2)
names(acc) <- c("Package", "Largest coefficient difference")
miss <- res[res$status %in% c("skipped", "killed", "error"), ]
if (nrow(miss) > 0) {
    miss <- aggregate(n ~ package + dgp + fe + k + status, miss, min)
    acc$Missing <- vapply(as.character(acc$Package), function(p) {
        m <- miss[miss$package == p, ]
        if (nrow(m) == 0) return("")
        why <- c(skipped = "not run (too slow at a smaller size)", killed = "stopped after 30 minutes",
                 error = "failed")[m$status]
        paste(sprintf("%s, %s, k = %d: %s from %s rows", m$dgp, m$fe, m$k, why,
                      format(m$n, big.mark = ",", scientific = FALSE, trim = TRUE)), collapse = "; ")
    }, "")
}
tinytable::tt(acc)
```
]

= Commentary

These comments describe the provisional run, up to 100,000 observations; they will be revised with the full grid.

- *In memory, easy cases: a tie, or `feols()` ahead.* With two fixed effects, or three on the simple design, all three in-memory packages fit 100,000 rows in a few hundredths of a second. With one regressor, `feols()` is about twice as fast as `design_lm()` (13 against 25 ms with two fixed effects); with ten regressors they are within 20% of each other, `design_lm()` slightly ahead. `felm()` is close behind.
- *In memory, the difficult design with three fixed effects: `design_lm()` ahead.* When individuals change firms every year, demeaning by alternating projections converges slowly. `design_lm()` solves for the fixed effects directly (eliminating the individual effects exactly and solving a small system for the others), so its time hardly depends on the design: at 100,000 rows with ten regressors it takes 61 ms, against 0.9 s for `feols()` and 32 s for `felm()`.
- *Memory.* The extra memory of `design_lm()` is about 25 to 30 MB at 100,000 rows, whatever the number of regressors, since it never forms the demeaned regressors; `feols()` uses less with one regressor (11 MB) and more with ten (39 MB), `felm()` 38 to 50 MB.
- *On DuckDB.* `design_lm()` on a DuckDB table has a fixed cost of 0.3 to 1 second (several SQL queries) that dominates at these sizes, so it is much slower than in memory. `dbreg()` is faster with two fixed effects, where a single demeaning pass is exact on this balanced panel; with three, it iterates alternating projections in SQL, which is slower than `design_lm()` at 100,000 rows on the simple design, and on the difficult design did not converge within 10,000 iterations at 100,000 rows (after about 15 minutes).
- *Accuracy.* `design_lm()`, in memory and on DuckDB, `felm()` and `dbreg()` agree with the tight-tolerance `feols()` coefficients to about $10^(-14)$. `feols()` at its default tolerance differs by up to $2 times 10^(-10)$, on the difficult design with three fixed effects.

= Code

The benchmark is `website/benchmarks/fixed_effects.R`, which runs each fit in a fresh R process with this script:

#bench-code("fixed_effects/cell")

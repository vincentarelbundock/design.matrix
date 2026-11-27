#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: [Fixed effects])
#metadata((title: "Fixed effects")) <website-metadata>

#import "/benchmarks/bench.typ": bench-code

#calepin.setup(echo: false, eval: true, results: "verbatim", warning: false, message: false)

// Building this page never runs the benchmark: the chunks read
// results/fixed_effects.csv and results/fixed_effects_scale.csv (if
// present), written by `Rscript website/benchmarks/run.R fixed_effects`.
// Older results files lack the reference, kill_s and ref_kill_s columns:
// the chunks fall back to the run's thresholds.
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
pkgs <- c("lm_xtx" = "#18756d", "lm_xtx (DuckDB)" = "#7fbfb3",
          "feols" = "#ad481e", "dbreg (DuckDB)" = "#c79a2c")
# (lfe::felm() was dropped from the benchmark: older result files list it)
res <- res[res$package %in% names(pkgs), ]
res$package <- factor(res$package, names(pkgs))
# A placeholder instead of an empty plot, when no fit finished
no_data <- function(what) {
    ggplot() + annotate("text", x = 0, y = 0, label = what, size = 5) + theme_void()
}
time_plot <- function(type) {
    x <- res[res$dgp == type & !is.na(res$seconds), ]
    if (nrow(x) == 0) return(no_data(paste("No timed fits on the", type, "design.")))
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
# The thresholds of the run, saved with the results (older results files
# predate that line and used the defaults, 120 s and 1800 s)
limits <- regmatches(machine, regexec("^Limits: .*after ([0-9.e+]+) s .*after ([0-9.e+]+) s", machine))
limits <- Filter(length, limits)
limit <- if (length(limits)) as.numeric(limits[[1]][2]) else 120
kill <- if (length(limits)) as.numeric(limits[[1]][3]) else 1800
duration <- function(s) {
    if (s >= 60 && s %% 60 == 0) {
        m <- s / 60
        if (m >= 60 && m %% 60 == 0) paste(m / 60, if (m == 60) "hour" else "hours")
        else paste(m, if (m == 1) "minute" else "minutes")
    } else paste(format(s), "seconds")
}
durations <- function(s) vapply(s, duration, "")
# Results saved before the cached data and separate references (no kill_s
# column) were measured with the previous protocol, and every cell was
# stopped after `kill` seconds, references included
old_protocol <- is.null(res$kill_s)
if (old_protocol) res$kill_s <- kill
if (is.null(res$ref_kill_s)) res$ref_kill_s <- kill
# Cells stopped sooner than the others (known not to converge)
short <- unique(res[res$kill_s < kill, c("package", "dgp", "kill_s")])
fmt_n <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)
# The at-scale grid, saved in its own file once it has run
scale_file <- "results/fixed_effects_scale.csv"
sc <- if (file.exists(scale_file)) read.csv(scale_file)
if (!is.null(sc)) {
    sc$fe <- factor(ifelse(sc$nfe == 2, "indiv + year", "indiv + year + firm"),
                    c("indiv + year", "indiv + year + firm"))
    sc$package <- factor(sc$package, names(pkgs))
}
scale_machine <- if (file.exists("results/fixed_effects_scale_machine.txt")) readLines("results/fixed_effects_scale_machine.txt")
```
]

#title()

How fast is `lm_xtx()` with fixed effects, compared with the other R packages that fit linear models with high-dimensional fixed effects? This page replicates the fixed-effects OLS benchmark of #link("https://github.com/py-econometrics/pyfixest")[pyfixest], itself adapted from Kyle Butts' #link("https://github.com/kylebutts/fixest_benchmarks")[fixest_benchmarks], for R packages only.

#calepin.chunk(echo: false, results: "typst")[
```r
if (provisional) cat("*Provisional results.* These numbers come from a reduced grid (", paste(fmt_n(sort(unique(res$n))), collapse = ", "), " observations; ", paste(sort(unique(res$k)), collapse = " and "), " regressors), run while validating the scripts. The full grid, up to 1 million observations, will replace them.", if (old_protocol) " They were also measured with an earlier version of the protocol below: each process simulated its own data, the DuckDB fits read an in-database table instead of the Parquet file, and the reference fit ran in the `feols()` process, after the timed fits.", "\n", sep = "")
```
]

= Data and models

The data are a panel of individuals observed over 10 years, with about 23 individuals per firm: for $n$ observations, $n slash 10$ individuals and $n slash 230$ firms. The regressors $x_1, dots, x_k$ are independent standard normal, with coefficients $beta_j = 1 slash j$. Individual, firm and year effects are standard normal, as is the error:
$ y = sum_(j=1)^k x_j beta_j + alpha_"indiv" + gamma_"firm" + delta_"year" + epsilon. $

The two designs differ only in how individuals are assigned to firms:

- *Simple*: each observation's firm is drawn at random.
- *Difficult*: firms are assigned round robin, so every individual changes firm every year. The firm and individual effects are then poorly separated, which slows down the alternating projections that demeaning methods use.

Each design is fitted with two sets of fixed effects, `y ~ x1 + ... + xk | indiv_id + year` and `y ~ x1 + ... + xk | indiv_id + year + firm_id`, with $k$ = 1, 5 or 10 regressors and iid standard errors, on 1,000 to 1 million observations. A separate grid, at scale (below), fits 10 and 50 million observations with fewer packages.

The packages:

- `lm_xtx()` (package `xtx.fit`, on the designs of `xtx`): on a data frame, and on the same data read by DuckDB from the Parquet file (`lm_xtx(f, con, table = "d")`, where `d` is a view on the file).
- `fixest::feols()`, with `vcov = "iid"`.
- `dbreg::dbreg()` on the same DuckDB view, with its automatic strategy, which here demeans the data in SQL: in a single pass for two fixed effects (the panel is balanced in individuals and years), by alternating projections for three. Its default limit of 100 alternating projections stops with an error on the difficult design with three fixed effects, so the benchmark raises it to 10,000.

`alpaca` is not included (not installed); it fits generalized linear models, and `feglm()`'s Gaussian family would go through the same demeaning as `feols()`.

= Protocol

- Everything runs on one core and one thread: `options(xtx.threads = 1)`, `setFixest_nthreads(1)`, `data.table::setDTthreads(1)`, DuckDB's `SET threads = 1`, and a single-threaded BLAS, with the whole run pinned to one CPU core.
- The data of each design and size are simulated once (same seed for every package) and saved as a Parquet file.
- Each combination of package, design, size, fixed effects and number of regressors runs in a fresh R process: read the Parquet file (into a data frame for the in-memory fits; the DuckDB fits query the file directly, through a view, so their times include reading it), fit once (warmup), then time three fits and keep the median.
- Memory is the peak resident memory of the process during the warmup fit, minus the resident memory just before it. It counts every allocation, in R, C, C++ or DuckDB, but not the data frame of the in-memory fits, read before.
#calepin.chunk(echo: false, results: "typst")[
```r
stopped_sooner <- if (nrow(short) > 0) {
    paste0(" (", paste(sprintf("%s on the %s design, which does not converge there: %s", short$package, short$dgp,
                               durations(short$kill_s)), collapse = "; "), ")")
}
cat("- Once a package's fit takes more than ", duration(limit), " at some size, or fails after that long, it is not run at larger sizes: those points are missing from the plots. A fit that takes longer than ", duration(limit), " is timed once, not three times, and a cell that runs longer than ", duration(kill), stopped_sooner, " is stopped.\n", sep = "")
cat("- Every result is compared with the coefficients of `feols()` run with a tight convergence tolerance (`fixef.tol = 1e-11`)",
    if (!old_protocol) paste0(", computed once per design, size, fixed effects and number of regressors, in a separate process outside the timed fits",
                              if (any(!is.na(res$ref_kill_s))) paste0(", and stopped after ", duration(max(res$ref_kill_s, na.rm = TRUE)))),
    "; the table at the end gives the largest absolute difference.\n", sep = "")
```
]

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
if (nrow(x) == 0) no_data("No memory measurements.") else
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

= At scale

The simple design with ten regressors and two or three fixed effects, at 10 and 50 million observations, for three fits only: `lm_xtx()` in memory, `lm_xtx()` on DuckDB reading the Parquet file, and `feols()`. The protocol is the same, with a time limit of its own, and without reference coefficients.

#calepin.chunk(echo: false, results: "typst")[
```r
if (is.null(sc)) {
    cat("_Not run yet._ These results will come with the full run of the benchmark.\n")
} else {
    if (any(sc$provisional)) cat("*Provisional results* (reduced grid: ", paste(fmt_n(sort(unique(sc$n))), collapse = ", "), " observations).\n\n", sep = "")
    m <- sc[sc$status %in% c("skipped", "killed", "error"), ]
    if (nrow(m) > 0) {
        cat("Missing points: ", paste(sprintf("%s, %s, %s rows: %s", m$package, m$fe, fmt_n(m$n),
            ifelse(m$status == "killed", paste("stopped after", durations(m$kill_s)),
                   ifelse(m$status == "skipped", "not run (too slow at a smaller size)", "failed"))), collapse = "; "),
            ".\n", sep = "")
    }
}
```
]

#calepin.chunk(
  echo: false,
  results: "render",
  fig-width: "100%",
  fig-device-width: 10,
  fig-device-height: 6,
  fig-alt-text: "Median fitting time and extra peak memory at 10 and 50 million observations, log scales, for lm_xtx in memory, lm_xtx on DuckDB and feols, with two and three fixed effects (columns).",
)[
```r
# No figure until the at-scale grid has run
if (!is.null(sc)) {
    x <- sc[!is.na(sc$seconds), ]
    long <- if (nrow(x) > 0) rbind(
        data.frame(x[c("package", "n", "fe")], what = "Median seconds", value = x$seconds),
        data.frame(x[c("package", "n", "fe")], what = "Extra peak memory, MB", value = pmax(x$mem_mb, 0.1)))
    if (!is.null(long)) long$what <- factor(long$what, c("Median seconds", "Extra peak memory, MB"))
    print(if (nrow(x) == 0) no_data("No fit finished at scale.") else
    ggplot(long, aes(n, value, colour = package)) +
        geom_line(linewidth = 0.8) + geom_point(size = 2) +
        facet_grid(what ~ fe, scales = "free_y", switch = "y") +
        scale_x_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                      breaks = sort(unique(x$n))) +
        scale_y_log10(labels = scales::label_number(drop0trailing = TRUE, big.mark = ",")) +
        scale_colour_manual(values = pkgs[levels(droplevels(x$package))]) +
        labs(x = "Observations (log scale)", y = NULL, colour = NULL) +
        theme_minimal(base_size = 13) +
        theme(legend.position = "bottom", strip.placement = "outside"))
}
```
]

#calepin.chunk(echo: false, results: "verbatim")[
```r
if (!is.null(scale_machine)) cat(scale_machine, sep = "\n")
```
]

= Accuracy and missing points

#calepin.chunk(echo: false, results: "typst")[
```r
cat("Largest absolute difference between each package's coefficients and those of `feols()` at a tight tolerance, over every fit; and the fits that were not run (more than ",
    duration(limit), " at a smaller size) or did not finish.\n", sep = "")
```
]

#calepin.chunk(echo: false, results: "typst")[
```r
# One row per package in the results, whether or not it has a comparison
present <- levels(droplevels(res$package))
acc <- data.frame(Package = present)
acc$diff <- vapply(present, function(p) {
    d <- res$max_coef_diff[res$package == p & !is.na(res$max_coef_diff)]
    if (length(d)) format(max(d), digits = 2)
    else if (old_protocol && !"feols" %in% present) "not available (feols reference not run)"
    else "not available (no reference for these fits)"
}, "")
names(acc)[2] <- "Largest coefficient difference"
# Failed, stopped or skipped reference fits (results saved since the
# reference column). Older results give it on the feols rows only, newer
# ones on every row of a design, size, fixed effects and k: one note each.
ref_miss <- if (!is.null(res$reference)) {
    unique(res[res$reference %in% c("error", "killed", "skipped"), c("dgp", "fe", "k", "n", "reference", "ref_kill_s")])
}
ref_notes <- if (!is.null(ref_miss) && nrow(ref_miss) > 0) {
    sprintf("%s, %s, k = %d, %s rows: reference fit %s", ref_miss$dgp, ref_miss$fe, ref_miss$k, fmt_n(ref_miss$n),
            ifelse(ref_miss$reference == "killed", paste("stopped after", durations(ref_miss$ref_kill_s)),
                   ifelse(ref_miss$reference == "skipped", "not run (stopped at a smaller size)", "failed")))
}
why <- c(skipped = "not run (too slow at a smaller size)", killed = "stopped", error = "failed")
miss <- res[res$status %in% names(why), ]
if (nrow(miss) > 0) miss <- aggregate(n ~ package + dgp + fe + k + status + kill_s, miss, min)
acc$Missing <- vapply(present, function(p) {
    m <- miss[miss$package == p, ]
    notes <- if (nrow(m) > 0) {
        sprintf("%s, %s, k = %d: %s from %s rows", m$dgp, m$fe, m$k,
                ifelse(m$status == "killed", paste("stopped after", durations(m$kill_s)), why[m$status]),
                fmt_n(m$n))
    }
    if (p == "feols") notes <- c(notes, ref_notes)
    if (length(notes)) paste(notes, collapse = "; ") else "none"
}, "")
tinytable::tt(acc)
```
]

#calepin.chunk(echo: false, results: "typst")[
```r
# Without feols in the results, the reference notes have no row of their own
if (!"feols" %in% present && length(ref_notes)) cat("Reference fits missing: ", paste(ref_notes, collapse = "; "), ".\n", sep = "")
```
]

= Commentary

These comments describe the provisional run, up to 100,000 observations; they will be revised with the full grid.

- *In memory, easy cases: a tie, or `feols()` ahead.* With two fixed effects, or three on the simple design, both in-memory packages fit 100,000 rows in a few hundredths of a second. With one regressor, `feols()` is about twice as fast as `lm_xtx()` (13 against 25 ms with two fixed effects); with ten regressors they are within 20% of each other, `lm_xtx()` slightly ahead.
- *In memory, the difficult design with three fixed effects: `lm_xtx()` ahead.* When individuals change firms every year, demeaning by alternating projections converges slowly. `lm_xtx()` solves for the fixed effects directly (eliminating the individual effects exactly and solving a small system for the others), so its time hardly depends on the design: at 100,000 rows with ten regressors it takes 61 ms, against 0.9 s for `feols()`.
- *Memory.* The extra memory of `lm_xtx()` is about 25 to 30 MB at 100,000 rows, whatever the number of regressors, since it never forms the demeaned regressors; `feols()` uses less with one regressor (11 MB) and more with ten (39 MB).
- *On DuckDB.* `lm_xtx()` on a DuckDB table has a fixed cost of 0.3 to 1 second (several SQL queries) that dominates at these sizes, so it is much slower than in memory. `dbreg()` is faster with two fixed effects, where a single demeaning pass is exact on this balanced panel; with three, it iterates alternating projections in SQL, which is slower than `lm_xtx()` at 100,000 rows on the simple design, and on the difficult design did not converge within 10,000 iterations at 100,000 rows (after about 15 minutes).
- *Accuracy.* `lm_xtx()`, in memory and on DuckDB, and `dbreg()` agree with the tight-tolerance `feols()` coefficients to about $10^(-14)$. `feols()` at its default tolerance differs by up to $2 times 10^(-10)$, on the difficult design with three fixed effects.

= Code

The benchmark is `website/benchmarks/fixed_effects.R`. It simulates each dataset once with `fixed_effects/simulate.R` and this data generating process:

#bench-code("fixed_effects/dgp")

then runs each fit in a fresh R process with this script (the reference fits run separately, in `fixed_effects/reference.R`):

#bench-code("fixed_effects/cell")

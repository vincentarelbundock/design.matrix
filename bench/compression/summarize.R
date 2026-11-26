## Tables of bench/compression/ANALYSIS.md from results/results-medium.csv
## and results/queries-1e+06.csv (kata design.matrix#8jp3). Base R only.
##   Rscript bench/compression/summarize.R
## Speedup = uncompressed / compressed median time (> 1: compression wins).

dir <- file.path("bench", "compression", "results")
d <- read.csv(file.path(dir, "results-medium.csv"), stringsAsFactors = FALSE)
d$speedup <- d$uncompressed_s / d$compressed_s
ops <- c("crossprod", "colSums_group", "lm", "wlm", "glm", "fe")
main <- d[d$keytype == "text" & d$threads == "default", ]

f2 <- function(x) ifelse(is.na(x), "-", formatC(x, format = "f", digits = if (all(abs(x) < 10, na.rm = TRUE)) 2 else 1))
cell <- function(x) paste(ifelse(is.na(x), "-", sprintf("%.2g", x)), collapse = " / ")

cat("## Grid coverage (expected 61 configurations x 6 ops = 366)\n")
cat("rows:", nrow(d), "\n")
cfg <- expand.grid(n = c(1e5, 1e6, 5e6), ratio_target = c(1e-4, 1e-3, 1e-2, 0.1, 0.5),
                   kind = c("factor", "discrete", "wide"), keytype = "text", threads = "default",
                   stringsAsFactors = FALSE)
side <- expand.grid(n = c(1e6, 5e6), ratio_target = c(1e-3, 0.1), kind = c("factor", "discrete"),
                    stringsAsFactors = FALSE)
cfg <- rbind(cfg, cbind(side, keytype = "enum", threads = "default"), cbind(side, keytype = "text", threads = "1"))
want <- merge(cfg, data.frame(op = ops))
key <- function(x) paste(x$n, x$ratio_target, x$kind, x$keytype, x$threads, x$op)
miss <- want[!key(want) %in% key(d), ]
cat("missing:", nrow(miss), "\n")
print(aggregate(op ~ n + kind + keytype + threads + ratio_target, miss, function(o) paste(o, collapse = ",")),
      row.names = FALSE)

cat("\n## Load average (1 min) when the cells ran, by n\n")
print(aggregate(load ~ n, d, function(x) round(c(min = min(x), median = median(x), max = max(x)), 1)))

cat("\n## Realized pattern ratio (patterns / rows) by target, kind (median over n)\n")
print(xtabs(ratio ~ ratio_target + kind, aggregate(ratio ~ ratio_target + kind, main, median)))

cat("\n## Speedup (uncompressed / compressed) by op and target ratio; cells: n = 1e5 / 1e6 / 5e6\n")
for (k in c("factor", "discrete", "wide")) {
    cat("\n### kind:", k, "\n")
    m <- sapply(sort(unique(main$ratio_target)), function(r) sapply(ops, function(o) {
        cell(sapply(c(1e5, 1e6, 5e6), function(n) {
            v <- main$speedup[main$kind == k & main$op == o & main$ratio_target == r & main$n == n]
            if (length(v)) v else NA
        }))
    }))
    colnames(m) <- paste0("r=", sort(unique(main$ratio_target)))
    print(noquote(m))
}

cat("\n## Seconds at n = 5e6, discrete (compressed / uncompressed)\n")
x <- main[main$n == 5e6 & main$kind == "discrete", ]
m <- sapply(sort(unique(x$ratio_target)), function(r) sapply(ops, function(o) {
    v <- x[x$op == o & x$ratio_target == r, ]
    sprintf("%.3f / %.3f", v$compressed_s, v$uncompressed_s)
}))
colnames(m) <- paste0("r=", sort(unique(x$ratio_target)))
print(noquote(m))

cat("\n## The auto decision: does it compress (C) or not (-), by target ratio and kind; cells n = 1e5/1e6/5e6\n")
for (g in c(FALSE, TRUE)) {
    cat(if (g) "colSums(group = g): keys + g\n" else "other ops: design keys\n")
    y <- main[if (g) main$op == "colSums_group" else main$op == "crossprod", ]
    m <- sapply(sort(unique(y$ratio_target)), function(r) sapply(c("factor", "discrete", "wide"), function(k)
        paste(sapply(c(1e5, 1e6, 5e6), function(n) {
            v <- y$auto_compresses[y$kind == k & y$ratio_target == r & y$n == n]
            if (!length(v)) "." else if (v) "C" else "-"
        }), collapse = "")))
    colnames(m) <- paste0("r=", sort(unique(y$ratio_target)))
    print(noquote(m))
}

cat("\n## Decision queries: median seconds by n and kind (and as % of the faster mode of crossprod)\n")
y <- main[main$op == "crossprod", ]
y$pct <- 100 * y$decision_s / pmin(y$compressed_s, y$uncompressed_s)
print(aggregate(cbind(decision_s, pct) ~ n + kind, y, function(v) round(median(v), 3)))

cat("\n## What auto costs vs the better mode: auto = decision + chosen mode; regret = auto / best\n")
a <- main[main$op != "fe", ]
a$auto_s <- a$decision_s + ifelse(a$auto_compresses, a$compressed_s, a$uncompressed_s)
a$best_s <- pmin(a$compressed_s, a$uncompressed_s, na.rm = TRUE)
a$regret <- a$auto_s / a$best_s
print(aggregate(regret ~ op, a, function(v) round(c(median = median(v), max = max(v)), 2)))
w <- a[order(-a$regret), c("n", "kind", "ratio_target", "op", "auto_compresses", "decision_s", "compressed_s",
                           "uncompressed_s", "regret", "load")]
cat("worst 12:\n")
print(head(w, 12), row.names = FALSE, digits = 3)

cat("\n## Side grids: speedup with ENUM keys, and with 1 thread, vs text keys/default threads (same n, ratio, kind)\n")
s <- d[d$keytype == "enum" | d$threads == "1", ]
b <- main[c("n", "ratio_target", "kind", "op", "speedup", "compressed_s", "uncompressed_s")]
s <- merge(s, b, by = c("n", "ratio_target", "kind", "op"), suffixes = c("", ".base"))
s$variant <- ifelse(s$keytype == "enum", "enum", "1 thread")
print(aggregate(cbind(speedup, speedup.base) ~ variant + op, s, function(v) round(median(v), 2)))
cat("time relative to text/default (median over n, ratio, kind):\n")
s$comp_rel <- s$compressed_s / s$compressed_s.base
s$unc_rel <- s$uncompressed_s / s$uncompressed_s.base
print(aggregate(cbind(comp_rel, unc_rel) ~ variant + op, s, function(v) round(median(v), 2)))
cat("decision flips between variant and base: ")
fl <- merge(d[d$keytype == "enum" | d$threads == "1", c("n", "ratio_target", "kind", "op", "speedup")],
            b[c("n", "ratio_target", "kind", "op", "speedup")], by = c("n", "ratio_target", "kind", "op"))
cat(sum(sign(log(fl$speedup.x)) != sign(log(fl$speedup.y)) & pmax(fl$speedup.x, 1 / fl$speedup.x, fl$speedup.y,
                                                                 1 / fl$speedup.y) > 1.2),
    "of", nrow(fl), "cells change winner where either side differs by more than 20%\n")

cat("\n## Memory: median over cells, by op and ratio (main grid, n >= 1e6)\n")
mm <- main[main$n >= 1e6, ]
print(aggregate(cbind(compressed_mem_mb, uncompressed_mem_mb, compressed_peak_rss_mb, uncompressed_peak_rss_mb) ~
                    op + ratio_target, mm, function(v) round(median(v), 1)))
cat("largest peak RSS (MB) per op:\n")
print(aggregate(cbind(compressed_peak_rss_mb, uncompressed_peak_rss_mb) ~ op, d, function(v) round(max(v), 0)))

cat("\n## Load sensitivity: within-cell speedups at high (> 6) vs low load, glm and crossprod, n >= 1e6\n")
h <- main[main$n >= 1e6 & main$op %in% c("glm", "crossprod"), ]
h$hi <- h$load > 6
print(aggregate(speedup ~ op + hi, h, function(v) round(c(n = length(v), median = median(v)), 2)))

cat("\n## Query profile at n = 1e6 (discrete, ratio ~1e-3): statements and seconds by label\n")
q <- read.csv(file.path(dir, "queries-1e+06.csv"), stringsAsFactors = FALSE)
q <- q[q$mode != "auto", ]
tot <- unique(q[c("op", "mode", "total", "iter")])
tot$statements <- as.vector(table(paste(q$op, q$mode))[paste(tot$op, tot$mode)])
tot$sql_s <- as.vector(tapply(q$seconds, paste(q$op, q$mode), sum)[paste(tot$op, tot$mode)])
print(tot, row.names = FALSE)
qa <- aggregate(seconds ~ op + mode + label, q, function(v) c(k = length(v), s = round(sum(v), 3)))
qa <- do.call(data.frame, qa)
print(qa[qa$op %in% c("glm", "wlm", "crossprod_wide") & qa$seconds.s >= 0.01, ], row.names = FALSE)

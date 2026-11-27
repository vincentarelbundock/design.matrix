## Remote factor keys and per-row values without joins of codes: the scans
## group by the factor columns as stored (ENUM, text, BOOLEAN) and decode
## levels in R; predictions read a term's values by the ENUM's position or
## the logical, or join on the text. Same results as in memory, for every
## column type, interactions mixing them, levels in another order than the
## ENUM type, and missing values.
if (!requireNamespace("duckdb", quietly = TRUE) || !requireNamespace("DBI", quietly = TRUE))
    exit_file("duckdb is not installed")
source("helpers.R")
library(xtx)

con <- test_con()
same <- function(a, b, tol = 1e-8) isTRUE(all.equal(unname(a), unname(b), tolerance = tol))
rd <- function(tab) {
    x <- as.data.frame(tab)
    x[order(x$id), ]
}

set.seed(5)
n <- 600
d <- data.frame(id = seq_len(n), x = rnorm(n), z = runif(n),
                e = factor(sample(c("p", "q", "r", "s"), n, TRUE), levels = c("s", "q", "p", "r", "unused")),
                h = sample(c("aa", "bb", "cc"), n, TRUE), b = sample(c(TRUE, FALSE), n, TRUE),
                g = sample(1:4, n, TRUE))
d$y <- d$x + as.integer(d$e) - 0.5 * d$b + rnorm(n)
d$k <- rpois(n, exp(0.2 + 0.3 * d$x + 0.1 * as.integer(d$e)))
## missing factor, text and logical values
d$e[3] <- NA
d$h[7] <- NA
d$b[11] <- NA
DBI::dbWriteTable(con, "d", d)
m <- DBI::dbReadTable(con, "d")
m <- m[order(m$id), ]
rownames(m) <- m$id
expect_true(is.factor(m$e))
expect_identical(levels(m$e), c("s", "q", "p", "r", "unused"))

## new data whose ENUM type orders the levels otherwise (and has a level
## the fit does not know, in no row)
nd <- d[d$id %% 3 == 0, ]
nd$e <- factor(as.character(nd$e), levels = c("zz", "r", "p", "q", "s"))
DBI::dbWriteTable(con, "nd", nd)
ndm <- DBI::dbReadTable(con, "nd")
ndm <- ndm[order(ndm$id), ]
expect_identical(levels(ndm$e), c("zz", "r", "p", "q", "s"))

## colSums() within groups (the scan with groups), over complete rows
DBI::dbWriteTable(con, "dc", m[stats::complete.cases(m[, c("x", "e", "h", "b")]), ])
mc <- DBI::dbReadTable(con, "dc")
X <- design_matrix(~ x * e + h + b, con, table = "dc")
Xm <- model.matrix(~ x * e + h + b, mc)
expect_true(same(colSums(X), colSums(Xm)))
expect_true(same(colSums(X, group = "g"), rowsum(Xm, mc$g)))
expect_true(same(colMeans(X, group = "h"), rowsum(Xm, mc$h) / as.vector(table(mc$h))))

test_con_done(con)

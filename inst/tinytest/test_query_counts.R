## The exact number of SQL statements each remote fit sends, counted with
## count_queries() (every DBI and duckdb_register call, including those that
## bypass dm_query()). A change in a count is a change in the remote plan:
## update the number only when the change is intended (fewer is better).
if (!requireNamespace("duckdb", quietly = TRUE) || !requireNamespace("DBI", quietly = TRUE))
    exit_file("duckdb is not installed")
library(design.matrix)
source("helpers.R")

con <- suppressMessages(DBI::dbConnect(duckdb::duckdb()))
set.seed(42)
n <- 300
d <- data.frame(x1 = rnorm(n), xd = sample(0:3, n, TRUE), f = sample(c("a", "b", "c"), n, TRUE),
                a = sample(1:8, n, TRUE), b = sample(1:5, n, TRUE), w = sample(1:3, n, TRUE))
d$y <- d$x1 + d$xd + rnorm(n)
d$yb <- as.integer(d$y > 1)
d$yp <- rpois(n, exp(0.3 * d$x1))
DBI::dbWriteTable(con, "d", d)

nq <- function(expr) count_queries(expr)$count

## LM
## (uncompressed LM and weighted LM: 1 fewer since nh9e, the fitted values
## X b were a view made at once and never read: the sums of squares come
## from the moments of the factorization; products are now made at their
## first use)
expect_equal(nq(design_lm(y ~ x1 + f, con, table = "d", compress = FALSE)), 11L, info = "lm")
expect_equal(nq(design_lm(y ~ xd + f, con, table = "d", compress = TRUE)), 9L, info = "lm compressed")
expect_equal(nq(design_lm(y ~ xd + f, con, table = "d")), 10L, info = "lm auto")
## weighted LM
expect_equal(nq(design_lm(y ~ x1 + f, con, table = "d", weights = "w", compress = FALSE)), 13L,
             info = "weighted lm")
expect_equal(nq(design_lm(y ~ xd + f, con, table = "d", weights = "w", compress = TRUE)), 11L,
             info = "weighted lm compressed")
## GLM logit, probit
expect_equal(nq(design_glm(yb ~ x1 + f, binomial, con, table = "d", compress = FALSE)), 43L,
             info = "logit")
expect_equal(nq(design_glm(yb ~ xd + f, binomial, con, table = "d", compress = TRUE)), 21L,
             info = "logit compressed")
expect_equal(nq(design_glm(yb ~ x1 + f, binomial("probit"), con, table = "d", compress = FALSE)), 48L,
             info = "probit")
expect_equal(nq(design_glm(yb ~ xd + f, binomial("probit"), con, table = "d", compress = TRUE)), 26L,
             info = "probit compressed")
## FE-LM (two FE), with clustering; FE-GLM Poisson with clustering (FE fits
## are never compressed; since etg8 their crossprod/colSums run no
## compression decision query: 3 fewer statements each)
expect_equal(nq(design_lm(y ~ x1 + xd | a + b, con, table = "d")), 31L, info = "FE-LM")
expect_equal(nq(design_lm(y ~ x1 + xd | a + b, con, table = "d", vcov = ~a)), 36L,
             info = "FE-LM cluster")
expect_equal(nq(design_glm(yp ~ x1 + xd | a + b, poisson, con, table = "d", vcov = ~a)), 173L,
             info = "FE-GLM Poisson cluster")

## the counter sees statements dm_query() does not log
op <- options(design.matrix.remote_debug = TRUE)
logged <- 0L
cnt <- withCallingHandlers(nq(design_lm(y ~ x1 + xd | a + b, con, table = "d", vcov = ~a)),
                           message = function(m) {
                               logged <<- logged + 1L
                               invokeRestart("muffleMessage")
                           })
options(op)
expect_true(cnt > logged)
## and untraces on error
expect_error(count_queries(stop("boom")), "boom")
expect_false(inherits(DBI::dbGetQuery, "functionWithTrace"))

DBI::dbDisconnect(con, shutdown = TRUE)

## Factor levels of remote fits are those of the rows used: a level present
## only in rows dropped for a missing response, weight, offset, FE or
## cluster value is dropped, as model.frame(drop.unused.levels = TRUE) drops
## it in memory (lm(), glm(), feols()). Same fit in memory and remotely, for
## design_lm(), design_glm() and fixed effects.
if (!requireNamespace("duckdb", quietly = TRUE) || !requireNamespace("DBI", quietly = TRUE))
    exit_file("duckdb is not installed")
library(design.matrix)

con <- DBI::dbConnect(duckdb::duckdb())
same <- function(a, b, tol = 1e-7) isTRUE(all.equal(unname(a), unname(b), tolerance = tol))

set.seed(11)
n <- 400
d <- data.frame(id = seq_len(n), x = rnorm(n), f = sample(c("a", "b", "c"), n, TRUE),
                g = sample(letters[1:8], n, TRUE), cl = sample(1:10, n, TRUE),
                w = runif(n, 0.5, 2), o = runif(n, -0.2, 0.2))
d$y <- d$x + as.integer(factor(d$f)) + rnorm(n)
d$k <- rpois(n, exp(0.3 + 0.2 * d$x))
## each of the levels "r", "s", "t", "u", "v" occurs in one row only, where
## the response, the weight, the offset, the FE or the cluster is missing
d$f[1:5] <- c("r", "s", "t", "u", "v")
d$y[1] <- NA
d$k[1] <- NA
d$w[2] <- NA
d$o[3] <- NA
d$g[4] <- NA
d$cl[5] <- NA
DBI::dbWriteTable(con, "d", d)
d <- DBI::dbReadTable(con, "d")
d <- d[order(d$id), ]
pq <- tempfile(fileext = ".parquet")
DBI::dbExecute(con, sprintf("COPY d TO '%s' (FORMAT PARQUET)", pq))

remotes <- list(
    sql = list(data = quote(con), extra = list(table = "d")),
    chunked = list(data = quote(con), extra = list(table = "d", engine = "chunked", chunk_size = 64)),
    parquet = list(data = quote(pq), extra = list()))
fit_on <- function(fun, args, data, extra = list()) eval(as.call(c(list(fun), args, list(data = data), extra)))

cases <- list(
    lm_response = list(fun = quote(design_lm), args = alist(y ~ x + f), ref = quote(lm)),
    lm_weights = list(fun = quote(design_lm), args = alist(y ~ x + f, weights = w), ref = quote(lm)),
    lm_offset = list(fun = quote(design_lm), args = alist(y ~ x + f, offset = o), ref = quote(lm)),
    lm_offset_term = list(fun = quote(design_lm), args = alist(y ~ x + offset(o) + f), ref = quote(lm)),
    glm_response = list(fun = quote(design_glm), args = alist(k ~ x + f, family = poisson), ref = quote(glm)),
    glm_weights = list(fun = quote(design_glm), args = alist(k ~ x + f, family = poisson, weights = w),
                       ref = quote(glm)),
    glm_offset = list(fun = quote(design_glm), args = alist(k ~ x + f, family = poisson, offset = o),
                      ref = quote(glm)),
    fe_response = list(fun = quote(design_lm), args = alist(y ~ x + f | g)),
    fe_weights = list(fun = quote(design_lm), args = alist(y ~ x + f | g, weights = w)),
    fe_cluster = list(fun = quote(design_lm), args = alist(y ~ x + f | g, vcov = ~cl)))

for (nm in names(cases)) {
    case <- cases[[nm]]
    m <- fit_on(case$fun, case$args, quote(d))
    if (!is.null(case$ref)) {
        ## in memory, as lm() and glm()
        l <- fit_on(case$ref, case$args, quote(d))
        expect_identical(names(coef(m)), names(coef(l)), info = paste(nm, "memory columns"))
        expect_true(same(coef(m), coef(l)), info = paste(nm, "memory coef"))
    }
    for (src in names(remotes)) {
        s <- remotes[[src]]
        r <- suppressMessages(fit_on(case$fun, case$args, s$data, s$extra))
        info <- paste(nm, src)
        expect_identical(names(coef(r)), names(coef(m)), info = paste(info, "columns"))
        expect_identical(r$xlevels, m$xlevels, info = paste(info, "xlevels"))
        expect_true(same(coef(r), coef(m)), info = paste(info, "coef"))
        expect_true(same(vcov(r), vcov(m)), info = paste(info, "vcov"))
        expect_equal(nobs(r), nobs(m), info = paste(info, "nobs"))
    }
}

DBI::dbDisconnect(con)

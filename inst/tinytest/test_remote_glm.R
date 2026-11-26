## design_glm() on remote data: the cases of test_glm.R (glm_cases.R), each
## fitted on a DuckDB table with both engines, compressed or not, and on
## Parquet files, and compared with glm() on the same rows.
if (!requireNamespace("duckdb", quietly = TRUE) || !requireNamespace("DBI", quietly = TRUE))
    exit_file("duckdb is not installed")
source("glm_cases.R", local = TRUE)
library(design.matrix)

con <- DBI::dbConnect(duckdb::duckdb())
d <- glm_data(n = 800)
DBI::dbWriteTable(con, "d", d)
dd <- DBI::dbReadTable(con, "d")
pq <- tempfile(fileext = ".parquet")
DBI::dbExecute(con, sprintf("COPY d TO '%s' (FORMAT PARQUET)", pq))

fit_remote <- function(cs, ...) {
    w <- if (length(cs) >= 4) cs[[4]]
    o <- if (length(cs) >= 5) cs[[5]]
    design_glm(cs[[2]], cs[[3]], ..., weights = w, offset = o)
}

for (cs in glm_cases) {
    b <- suppressWarnings(glm_ref(cs, dd))
    for (eng in c("sql", "chunked")) {
        for (cp in list("auto", FALSE)) {
            a <- suppressWarnings(fit_remote(cs, con, table = "d", engine = eng, compress = cp))
            glm_same(a, b, paste(cs[[1]], eng, "compress", cp))
        }
    }
}
# Parquet, a few cases
for (cs in glm_cases[c(1, 7, 9, 14)]) {
    b <- suppressWarnings(glm_ref(cs, dd))
    for (eng in c("sql", "chunked")) {
        a <- suppressWarnings(fit_remote(cs, pq, engine = eng))
        glm_same(a, b, paste(cs[[1]], "parquet", eng))
    }
}

# an offset with few values: compressed by pattern and offset
fo <- design_glm(yp ~ g + offset(log(xd)), poisson, con, table = "d")
expect_identical(fo$method, "patterns")
glm_same(fo, glm(yp ~ g + offset(log(xd)), poisson, dd), "offset with few values")

# which path ran
expect_identical(design_glm(yb ~ xd * f + g, binomial, con, table = "d")$method, "patterns")
expect_identical(design_glm(yb ~ x + g, binomial, con, table = "d", engine = "sql")$method, "SQL iterations")
expect_identical(design_glm(yb ~ x + g, binomial, con, table = "d", engine = "chunked")$method,
                 "iterations in batches")

# warnings, as glm()
w1 <- glm_warnings(glm(sep ~ x, binomial, dd))$warnings
for (eng in c("sql", "chunked")) {
    expect_identical(glm_warnings(design_glm(sep ~ x, binomial, con, table = "d", engine = eng))$warnings, w1,
                     info = paste("separation", eng))
}
w2 <- glm_warnings(glm(yb ~ x + f, binomial, dd, control = list(maxit = 2)))$warnings
for (eng in c("sql", "chunked")) {
    expect_identical(glm_warnings(design_glm(yb ~ x + f, binomial, con, table = "d", engine = eng,
                                             control = list(maxit = 2)))$warnings, w2,
                     info = paste("maxit", eng))
}
expect_identical(glm_warnings(design_glm(yb ~ xd + f, binomial, con, table = "d", control = list(maxit = 2)))$warnings,
                 glm_warnings(glm(yb ~ xd + f, binomial, dd, control = list(maxit = 2)))$warnings)
expect_error(design_glm(yp ~ x, Gamma, con, sql = "SELECT x, yp - 3 AS yp FROM d"), "non-positive values")
expect_warning(design_glm(prop ~ xd, binomial, con, sql = "SELECT xd, prop * 0.9 AS prop, five FROM d",
                          weights = "five"), "non-integer")
expect_error(design_glm(yp ~ x + g, poisson, con, table = "d", na.action = na.fail), "missing values")

# predictions, fitted values and residuals, as glm()
b <- glm(yp ~ x + g + offset(log(e)), poisson, dd)
for (eng in c("sql", "chunked")) {
    a <- design_glm(yp ~ x + g + offset(log(e)), poisson, con, table = "d", engine = eng)
    for (ty in c("link", "response")) {
        p <- as.data.frame(predict(a, newdata = con, table = "d", type = ty, engine = eng))
        expect_equal(p$.fitted[order(p$id)], unname(predict(b, dd, type = ty)), info = paste(eng, ty))
        p <- as.data.frame(predict(a, type = ty))
        expect_equal(p$.fitted[order(p$id)], unname(predict(b, type = ty)), info = paste(eng, ty, "training"))
        expect_equal(predict(a, newdata = dd[1:5, ], type = ty), predict(b, dd[1:5, ], type = ty))
    }
    f <- as.data.frame(fitted(a))
    expect_equal(f$.fitted[order(f$id)], unname(fitted(b)), info = eng)
    for (ty in c("deviance", "pearson", "working", "response")) {
        r <- residuals(a, type = ty)
        expect_inherits(r, "dm_remote")
        r <- as.data.frame(r$src$view)
        expect_equal(r$.resid[order(r$id)], unname(residuals(b, type = ty)), info = paste(eng, ty))
    }
}
# an in-memory fit scores remote data
am <- design_glm(yp ~ x + g + offset(log(e)), poisson, dd)
p <- as.data.frame(predict(am, newdata = con, table = "d", type = "response"))
expect_equal(p$.fitted[order(p$id)], unname(predict(b, dd, type = "response")))

# a factor (character) response, and a remote design with y
bf <- glm(factor(fy) ~ x, binomial, dd)
expect_equal(coef(design_glm(fy ~ x, binomial, con, table = "d")), coef(bf))
D <- design_matrix(~ x + g, con, table = "d")
expect_equal(coef(design_glm(D, binomial, y = "yb")), coef(glm(yb ~ x + g, binomial, dd)))
expect_equal(coef(design_glm(D, binomial, y = c("s", "fl"))), coef(glm(cbind(s, fl) ~ x + g, binomial, dd)))
expect_equal(coef(design_glm(D, poisson, y = "yp", offset = "le")),
             coef(glm(yp ~ x + g, poisson, dd, offset = le)))

# start values
s0 <- coef(glm(yp ~ xd + g, poisson, dd)) + 0.1
for (eng in c("sql", "chunked")) for (cp in list("auto", FALSE)) {
    a <- design_glm(yp ~ xd + g, poisson, con, table = "d", start = s0, engine = eng, compress = cp)
    b <- glm(yp ~ xd + g, poisson, dd, start = s0)
    expect_equal(a$iter, b$iter, info = paste("start", eng, cp))
    expect_equal(coef(a), coef(b), info = paste("start", eng, cp))
}

DBI::dbDisconnect(con)

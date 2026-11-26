## Remote factor keys and per-row values without joins of codes: the scans
## group by the factor columns as stored (ENUM, text, BOOLEAN) and decode
## levels in R; predictions read a term's values by the ENUM's position or
## the logical, or join on the text. Same results as in memory, for every
## column type, interactions mixing them, levels in another order than the
## ENUM type, and missing values.
if (!requireNamespace("duckdb", quietly = TRUE) || !requireNamespace("DBI", quietly = TRUE))
    exit_file("duckdb is not installed")
library(design.matrix)

con <- DBI::dbConnect(duckdb::duckdb())
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

forms <- list(
    enum = y ~ x * e + z,
    text = y ~ x * h + z,
    logical = y ~ x * b + e,
    integer_factor = y ~ factor(g) + x,
    enum_text = y ~ x + e:h + z,
    enum_logical_text = y ~ e * b + h:x,
    three = y ~ x:e:h + b)

for (nm in names(forms)) {
    f <- forms[[nm]]
    for (na in c("omit", "exclude")) {
        nact <- if (na == "omit") stats::na.omit else stats::na.exclude
        lab <- paste(nm, na)
        r <- suppressMessages(design_lm(f, con, table = "d", na.action = nact))
        l <- lm(f, m, na.action = nact)
        expect_identical(names(coef(r)), names(coef(l)), info = lab)
        expect_true(same(coef(r), coef(l)), info = lab)
        ## fitted values and residuals, as tables
        fr <- rd(fitted(r))
        fl <- fitted(l)
        expect_equal(nrow(fr), length(fl), info = lab)
        expect_true(same(fr$.fitted, fl[as.character(fr$id)]), info = paste(lab, "fitted"))
        expect_true(same(fr$.resid, residuals(l)[as.character(fr$id)]), info = paste(lab, "residuals"))
        ## predictions on new data with another ENUM type
        pr <- rd(suppressMessages(predict(r, con, table = "nd")))
        pl <- predict(l, ndm)
        expect_true(same(pr[[ncol(pr)]], pl), info = paste(lab, "predict"))
    }
    ## GLM: the linear predictor of each iteration
    gf <- stats::update(f, k ~ .)
    gr <- suppressMessages(design_glm(gf, poisson, con, table = "d", compress = FALSE))
    gl <- glm(gf, poisson, m)
    expect_true(same(coef(gr), coef(gl), 1e-7), info = paste(nm, "glm"))
    expect_true(same(deviance(gr), deviance(gl), 1e-7), info = paste(nm, "glm deviance"))
    pr <- rd(suppressMessages(predict(gr, con, table = "nd", type = "response")))
    expect_true(same(pr[[ncol(pr)]], predict(gl, ndm, type = "response"), 1e-7), info = paste(nm, "glm predict"))
}

## colSums() within groups (the scan with groups), over complete rows
DBI::dbWriteTable(con, "dc", m[stats::complete.cases(m[, c("x", "e", "h", "b")]), ])
mc <- DBI::dbReadTable(con, "dc")
X <- design_matrix(~ x * e + h + b, con, table = "dc")
Xm <- model.matrix(~ x * e + h + b, mc)
expect_true(same(colSums(X), colSums(Xm)))
expect_true(same(colSums(X, group = "g"), rowsum(Xm, mc$g)))
expect_true(same(colMeans(X, group = "h"), rowsum(Xm, mc$h) / as.vector(table(mc$h))))

DBI::dbDisconnect(con, shutdown = TRUE)

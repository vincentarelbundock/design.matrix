## Remote designs on DuckDB: results agree with the design in memory.
if (!requireNamespace("duckdb", quietly = TRUE) || !requireNamespace("DBI", quietly = TRUE))
    exit_file("duckdb is not installed")
source("helpers.R")
library(design.matrix)

con <- DBI::dbConnect(duckdb::duckdb())
scratch <- tempfile("dm_remote")
dir.create(scratch)

## data as the database returns it (ordered factors become plain factors,
## row names are dropped), with an id to match rows of predictions
remote_data <- function(n, na = TRUE) {
    d <- fuzz_data(n, special = FALSE, na = na)
    d$m <- NULL
    rownames(d) <- NULL
    d$id <- seq_len(n)
    d$w <- runif(n)
    d$g <- sample(c(3, 1, 2, NA), n, TRUE, prob = c(0.3, 0.3, 0.3, 0.1))
    d
}

same <- function(a, b) isTRUE(all.equal(unname(a), unname(b), tolerance = 1e-7))

## the rows of `dd` that the design in memory keeps
kept_rows <- function(fo, dd, spec = NULL) {
    mf <- if (is.null(spec)) {
        design_frame(fo, dd)
    } else {
        design_frame(spec$terms, dd, xlev = spec$xlevels, na.action = stats::na.omit)
    }
    as.integer(rownames(mf))
}

## compare a remote design with the design in memory, on every operation
check_remote <- function(Dr, M, dd, rows, info) {
    expect_equal(dim(Dr), dim(M), info = info)
    expect_identical(colnames(Dr), colnames(M), info = info)
    expect_true(same(crossprod(Dr), crossprod(M)), info = paste(info, "X'X"))
    expect_true(same(crossprod(Dr, "w"), crossprod(M, dd$w[rows])), info = paste(info, "X'w"))
    mw <- design_moments(Dr, weights = "w")
    expect_true(same(mw$xtx, crossprod(M, dd$w[rows] * M)), info = paste(info, "X'WX"))
    expect_equal(mw$sumw, sum(dd$w[rows]), info = paste(info, "sum of weights"))
    expect_true(same(colSums(Dr), colSums(M)), info = paste(info, "colSums"))
    g <- dd$g[rows]
    wg <- if (anyNA(g)) "missing values" else NA
    a <- capture(function() colSums(Dr, group = "g"))
    b <- capture(function() colSums(M, group = g))
    expect_true(same(a$value, b$value), info = paste(info, "colSums(group)"))
    expect_identical(rownames(a$value), rownames(b$value), info = paste(info, "group names"))
    expect_identical(a$warnings, b$warnings, info = paste(info, "group warning"))
    a <- suppressWarnings(colMeans(Dr, group = "g", weights = "w"))
    b <- suppressWarnings(colMeans(dd$w[rows] * M, group = g))
    expect_true(same(a, b), info = paste(info, "colMeans(group, weights)"))
    beta <- rnorm(ncol(M))
    p <- as.data.frame(Dr %*% beta)
    p <- p[order(p$id), ]
    expect_identical(p$id, dd$id[rows], info = paste(info, "prediction rows"))
    expect_true(same(p$.fitted, drop(M %*% beta)), info = paste(info, "predictions"))
}

set.seed(20261125)
dd <- remote_data(400)
DBI::dbWriteTable(con, "d", dd)
dd <- DBI::dbReadTable(con, "d")
dd <- dd[order(dd$id), ]
rownames(dd) <- NULL

## random formulas, both engines
checked <- c(sql = 0, chunked = 0)
for (i in 1:40) {
    fo <- fuzz_formula()
    ## (a remote design needs a plan, which refuses the response on the right)
    if ("m" %in% all.vars(fo) || "y" %in% all.vars(fo[[3L]])) next
    ctr <- fuzz_contrasts(dd)
    M <- capture(function() design_matrix(fo, dd, contrasts.arg = ctr))
    for (engine in c("auto", "chunked")) {
        info <- paste(deparse(fo), engine)
        Dr <- capture(function() suppressMessages(
            design_matrix(fo, con, table = "d", contrasts.arg = ctr, engine = engine, chunk_size = 37)))
        if (inherits(M$value, "captured_error")) {
            expect_inherits(Dr$value, "captured_error", info = info)
            next
        }
        if (inherits(Dr$value, "captured_error")) {
            expect_true(FALSE, info = paste(info, Dr$value$msg))
            next
        }
        Dr <- Dr$value
        checked[Dr$engine] <- checked[Dr$engine] + 1
        check_remote(Dr, M$value, dd, kept_rows(fo, dd), info)
        m <- design_moments(Dr)
        ref <- design_moments(M$value, dd$y[kept_rows(fo, dd)])
        expect_true(same(m$xty, ref$xty) && same(m$yty, ref$yty) && m$n == ref$n, info = paste(info, "moments"))
    }
}
expect_true(all(checked > 5))

## fit in R, score in the database: a specification from part of the data
for (i in 1:25) {
    fo <- fuzz_formula(response = FALSE)
    if ("m" %in% all.vars(fo)) next
    train <- dd[1:120, ]
    spec <- capture(function() design_matrix(fo, train, contrasts.arg = fuzz_contrasts(train)))
    if (inherits(spec$value, "captured_error")) next
    spec <- spec$value
    plan <- design_matrix(spec, newdata = NULL)
    M <- capture(function() design_matrix(spec, dd, na.action = stats::na.omit))
    for (engine in c("auto", "chunked")) {
        info <- paste("spec", deparse(fo), engine)
        Dr <- capture(function() suppressMessages(design_matrix(spec, con, table = "d", engine = engine, chunk_size = 50)))
        if (inherits(M$value, "captured_error")) {
            ## new levels: the same error as in memory
            expect_inherits(Dr$value, "captured_error", info = info)
            if (inherits(Dr$value, "captured_error")) expect_identical(Dr$value$msg, M$value$msg, info = info)
            next
        }
        check_remote(Dr$value, M$value, dd, kept_rows(NULL, dd, plan), info)
    }
}

## new levels, explicitly
train <- data.frame(x = rnorm(10), ch = rep(c("u", "v"), 5))
spec <- design_matrix(~ x + ch, train)
DBI::dbWriteTable(con, "newlev", data.frame(x = 1:4, ch = c("u", "w", "a", NA)))
msg <- tryCatch(design_matrix(spec, DBI::dbReadTable(con, "newlev"), na.action = stats::na.omit),
                error = conditionMessage)
expect_error(design_matrix(spec, con, table = "newlev"), msg, fixed = TRUE)
expect_error(design_matrix(spec, con, table = "newlev", engine = "chunked"), msg, fixed = TRUE)

## NaN and NA rows, transformations giving NaN
na <- data.frame(y = c(1, 2, NA, 4, 5, 6, 7, 8), x = c(1, NaN, 3, 4, -1, 6, 7, 8),
                 f = c("a", "b", "a", NA, "b", "a", "b", "a"), w = 1, g = 1, id = 1:8)
DBI::dbWriteTable(con, "na", na)
fo <- y ~ log(x) * f
M <- suppressWarnings(design_matrix(fo, na))
Dr <- design_matrix(fo, con, table = "na")
expect_identical(Dr$engine, "sql")
expect_equal(nrow(Dr), nrow(M))
expect_true(same(crossprod(Dr), crossprod(M)))
expect_true(same(design_moments(Dr)$xty, crossprod(M, na$y[c(1, 6, 7, 8)])))
expect_error(crossprod(Dr, weights = "w"), "design_moments")
mw <- design_moments(Dr, weights = "w")
expect_true(same(mw$xtx, crossprod(M)) && same(mw$xty, crossprod(M, na$y[c(1, 6, 7, 8)])))
## a column with missing values in the rows used gives NA, as in memory
DBI::dbWriteTable(con, "na2", transform(na, w = c(1, 1, 1, 1, 1, NA, 1, 1)))
Dr2 <- design_matrix(fo, con, table = "na2")
expect_true(all(is.na(crossprod(Dr2, "w"))))
expect_true(all(is.na(suppressMessages(design_matrix(fo, con, table = "na2", engine = "chunked"))
                      |> crossprod(y = "w"))))

## engine choice
expect_message(design_matrix(y ~ poly(z, 2) + f1, con, table = "d"),
               "chunked engine: term `poly(z, 2)` cannot be translated", fixed = TRUE)
expect_error(design_matrix(y ~ poly(z, 2) + f1, con, table = "d", engine = "sql"),
             "SQL engine cannot compute this design", fixed = TRUE)
expect_silent(D <- design_matrix(y ~ x1 * f1 + I(x2^2) + log(z, 10), con, table = "d"))
expect_identical(D$engine, "sql")
expect_stdout(design_explain(D, sql = TRUE), "GROUP BY")
expect_stdout(print(D), "all terms pushed down")

## sql= source, as.matrix, reuse of a remote design
D <- design_matrix(y ~ x1 + f1, con, sql = "SELECT * FROM d WHERE id <= 200")
M <- design_matrix(y ~ x1 + f1, dd[dd$id <= 200, ])
expect_true(same(crossprod(D), crossprod(M)))
expect_equal(unname(as.matrix(D)), unname(as.matrix(M)))
op <- options(design.matrix.remote_max_rows = 10)
expect_error(as.matrix(D), "force = TRUE")
expect_equal(unname(as.matrix(D, force = TRUE)), unname(as.matrix(M)))
options(op)
expect_equal(as.matrix(design_matrix(D, dd[1:5, ])), as.matrix(design_matrix(M, dd[1:5, ])))
D2 <- design_matrix(D, con, table = "d")
expect_true(same(crossprod(D2), crossprod(design_matrix(M, dd, na.action = stats::na.omit))))

## persisted predictions
beta <- seq_len(ncol(D))
tab <- design_predict(D, beta, name = "preds", keep = "id", column = "eta")
expect_inherits(tab, "dm_remote_table")
p <- DBI::dbGetQuery(con, "SELECT * FROM preds ORDER BY id")
expect_identical(names(p), c("id", "eta"))
expect_true(same(p$eta, drop(M %*% beta)))
expect_error(design_predict(D, beta, name = "preds"))
expect_silent(design_predict(D, beta, name = "preds", overwrite = TRUE))

## several prediction columns; a fitted model scored with `newdata`
fit <- lm(y ~ x1 + f1 + splines::ns(z, 2), dd)
B <- cbind(a = coef(fit), b = -coef(fit))
nd <- dd[stats::complete.cases(dd[c("x1", "f1", "z")]), ]
ref <- as.matrix(design_matrix(fit, newdata = nd)) %*% B
for (engine in c("auto", "chunked")) {
    Df <- suppressMessages(design_matrix(fit, newdata = con, table = "d", engine = engine))
    p <- as.data.frame(design_predict(Df, B, keep = "id"))
    p <- p[order(p$id), ]
    expect_identical(names(p), c("id", "a", "b"))
    expect_true(same(as.matrix(p[c("a", "b")]), ref), info = engine)
}

## Parquet files, with a private connection
pq <- file.path(scratch, "d.parquet")
DBI::dbExecute(con, sprintf("COPY (SELECT * FROM d) TO '%s' (FORMAT PARQUET)", pq))
dp <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s') ORDER BY id", pq))
for (fo in list(y ~ x1 * f2 + ch + lg, y ~ splines::ns(z, 3) + f1)) {
    M <- design_matrix(fo, dp)
    for (engine in c("auto", "chunked")) {
        Dp <- suppressMessages(design_matrix(fo, pq, engine = engine, chunk_size = 64))
        check_remote(Dp, M, dp, kept_rows(fo, dp), paste("parquet", deparse(fo), engine))
    }
}

## errors
expect_error(design_matrix(y ~ x1, con), "exactly one of")
expect_error(design_matrix(y ~ nope, con, table = "d"))
expect_error(crossprod(design_matrix(y ~ x1, con, table = "d"), "nope"), "no column")

DBI::dbDisconnect(con)
unlink(scratch, recursive = TRUE)

## The term layout of a remote design (D$sql_terms: the factors, numeric
## columns, design columns and contrast block C of each part of each term)
## comes from the set-up of the prepared prototype design, xtx_sql_layout().
## It must reproduce, entry for entry, the layout read off a design built on
## a grid of levels (xtx_sql_layout_grid(), the former and now fallback
## extractor), and the design must be that of model.matrix().
if (!requireNamespace("duckdb", quietly = TRUE) || !requireNamespace("DBI", quietly = TRUE))
    exit_file("duckdb is not installed")
library(xtx)
source("helpers.R")
ns <- asNamespace("xtx")

con <- test_con()
set.seed(11)
n <- 300
d <- data.frame(id = seq_len(n), y = rnorm(n), x = rnorm(n), z = runif(n, 1, 2),
                f = sample(letters[1:4], n, TRUE), g = sample(LETTERS[1:3], n, TRUE),
                h = sample(c("u", "v"), n, TRUE), l = sample(c(TRUE, FALSE), n, TRUE),
                stringsAsFactors = FALSE)
## an ENUM with a dictionary in another order than the sorted levels, and
## a level in no row
d$e <- factor(sample(c("s", "q", "p"), n, TRUE), levels = c("s", "q", "p", "unused"))
d$o <- factor(sample(c("lo", "mid", "hi"), n, TRUE), levels = c("lo", "mid", "hi"), ordered = TRUE)
DBI::dbWriteTable(con, "d", d)
m <- DBI::dbReadTable(con, "d")

## the layout of D from the set-up is the layout from the grid, and D is
## the design `X` (a model matrix of the rows of `m`, in their order)
n_checked <- 0L
check_layout <- function(D, X, info) {
    expect_identical(D$sql_layout, "setup", info = paste(info, ": layout from the set-up"))
    old <- ns$xtx_sql_prepare(D, NULL)
    expect_identical(old$sql_layout, "grid", info = info)
    expect_identical(length(D$sql_terms), length(old$sql_terms), info = paste(info, ": parts"))
    for (field in c("term", "S", "L", "Lt", "M", "strides", "cols", "C"))
        expect_identical(lapply(D$sql_terms, `[[`, field), lapply(old$sql_terms, `[[`, field),
                         info = paste(info, ":", field))
    ## every column of the design belongs to the parts of its term
    cols <- lapply(D$sql_terms, function(t) t$cols)
    expect_identical(sort(unique(as.integer(unlist(cols)))), seq_len(ncol(D)), info = paste(info, ": columns"))
    for (t in D$sql_terms)
        expect_identical(dim(t$C), c(t$Lt, length(t$cols)), info = paste(info, ": dim of C"))
    ## assignment, labels and values against the dense reference
    if (!is.null(X)) {
        expect_identical(as.character(colnames(D)), as.character(colnames(X)), info = paste(info, ": names"))
        expect_identical(D$assign, attr(X, "assign"), info = paste(info, ": assign"))
        expect_equal(unname(as.matrix(D, force = TRUE)[order(m$id), , drop = FALSE]),
                     unname(X[order(m$id), , drop = FALSE]), info = paste(info, ": values"))
        if (ncol(X)) {
            expect_equal(unname(crossprod(D)), unname(crossprod(X)), info = paste(info, ": crossprod"))
            expect_equal(unname(colSums(D)), unname(colSums(X)), info = paste(info, ": colSums"))
        }
    }
    n_checked <<- n_checked + 1L
}
remote <- function(f, ...) suppressMessages(design_matrix(f, con, table = "d", ...))

## formulas: intercept only, zero columns, no intercept (which changes the
## coding of the first factor), no response, logical factors, matrix
## variables (alone, several, with factors: their parts take nonconsecutive
## columns), nested and reordered interactions, the unused level of an ENUM
formulas <- list(
    y ~ 1, y ~ 0, y ~ x, ~ x + f, ~ 0 + f, y ~ 0 + f + x, y ~ 0 + x + f, y ~ 0 + f:g + x:g,
    y ~ f * g * x, y ~ l + x:l, y ~ 0 + l, y ~ 0 + x:f + l:f, y ~ f %in% g, y ~ g:f + x:g:f + f,
    y ~ x:f + z, y ~ f * g * h * l, y ~ log(z) * f + I(x^2):g:f, y ~ e * x + e:g, y ~ 0 + e + e:f,
    y ~ poly(x, 2), y ~ poly(x, 2) * f * g, y ~ poly(x, 2):poly(z, 3):f + g, y ~ 0 + poly(x, 2):f + f,
    y ~ f:poly(x, 2):g:poly(z, 2), y ~ o + x:o, y ~ x + offset(z) + f)
contrasts_list <- list(
    default = NULL,
    named = list(f = "contr.sum", g = "contr.helmert", e = "contr.SAS", o = "contr.treatment"),
    ## a selection (treatment with another base) and dense custom matrices,
    ## with and without column names
    custom = list(f = contr.treatment(4, base = 3),
                  g = matrix(c(1, 2, 3, -1, 0.5, 7), 3, 2, dimnames = list(NULL, c("p", "q"))),
                  e = matrix(c(1, 0, 0, 0, 0.25, -3, 2, 0), 4, 2)))
for (cn in names(contrasts_list)) for (f in formulas) {
    ca <- contrasts_list[[cn]]
    ca <- ca[names(ca) %in% all.vars(f)]
    if (!length(ca)) ca <- NULL
    info <- paste(deparse(f), cn)
    D <- remote(f, contrasts.arg = ca)
    X <- model.matrix(f, m, contrasts.arg = ca)
    check_layout(D, X, info)
}
expect_identical(n_checked, 3L * length(formulas))

## altered contrast options, in force when the design is built
op <- options(contrasts = c("contr.sum", "contr.poly"))
check_layout(remote(y ~ f * x + g), model.matrix(y ~ f * x + g, m), "options(contrasts)")
options(contrasts = c("contr.helmert", "contr.treatment"))
check_layout(remote(y ~ 0 + g + f:g), model.matrix(y ~ 0 + g + f:g, m), "options(contrasts) 2")
options(op)

## a plan trained in memory (ordered and logical factors, a matrix
## variable, levels given by the training data), applied to the remote rows
train <- m
train$f <- factor(train$f, levels = c("d", "b", "a", "c", "never"))
fo <- y ~ o * x + f + l:poly(z, 2)
ca <- list(f = "contr.sum")
Dm <- design_matrix(fo, train, contrasts.arg = ca)
D <- suppressMessages(design_matrix(Dm, newdata = con, table = "d"))
check_layout(D, model.matrix(fo, train, contrasts.arg = ca), "plan with ordered factor and unused level")

## a fitted model whose method drops the intercept after the coding: the
## factor keeps its contrasts (not the dummy coding of y ~ 0 + f)
plan <- ns$xtx_new_plan(stats::terms(y ~ f * x + g), list(f = sort(unique(m$f)), g = sort(unique(m$g))),
                        list(f = "contr.treatment", g = "contr.sum"), NULL, drop_intercept = TRUE)
D <- suppressMessages(design_matrix(plan, newdata = con, table = "d"))
X <- model.matrix(y ~ f * x + g, m, contrasts.arg = list(f = "contr.treatment", g = "contr.sum"))
Xd <- X[, -1L, drop = FALSE]
attr(Xd, "assign") <- attr(X, "assign")[-1L]
expect_false(any(vapply(D$terms, function(t) t$index == 0L, NA)))
check_layout(D, Xd, "intercept dropped after coding")
D0 <- remote(y ~ 0 + f * x + g, contrasts.arg = list(f = "contr.treatment", g = "contr.sum"))
expect_false(identical(colnames(D0), colnames(D)))

## a term with no columns (a gap in the assignment) cannot reach a remote
## design: the response on the right-hand side is refused, and stats
## itself fails on a matrix variable or a contrast matrix with no column
expect_error(remote(y ~ y + x), "response also appears on the right-hand side")
z0 <- function(x) matrix(0, length(x), 0L)
expect_error(suppressWarnings(remote(y ~ x + z0(z) + f)))
expect_error(suppressWarnings(remote(y ~ x + f, contrasts.arg = list(f = matrix(0, 4, 0)))))

## the layout does not depend on the values of the prototype rows: with
## non-finite and huge values, which the in-memory structured kernels
## refuse, it is the same as without
d2 <- d
d2$x[c(1, 5)] <- c(Inf, 1e300)
d2$z[2] <- NaN
DBI::dbWriteTable(con, "d2", d2)
Da <- suppressMessages(design_matrix(y ~ x * f + poly(z, 2):g, con, table = "d"))
Db <- suppressMessages(design_matrix(Da, newdata = con, table = "d2"))
expect_identical(Db$sql_layout, "setup")
expect_identical(Db$sql_terms, ns$xtx_sql_prepare(Db, NULL)$sql_terms)
## (the variables are numbered as the plan has them, the rest is the same)
for (field in c("term", "L", "Lt", "strides", "cols", "C"))
    expect_identical(lapply(Db$sql_terms, `[[`, field), lapply(Da$sql_terms, `[[`, field), info = field)

## an interaction of more variables than the structured kernels take
## (SMAXC, 32, in src/structured.h): the layout needs no such limit
k <- 33L
dk <- as.data.frame(matrix(runif(20 * k, 1, 2), 20, k, dimnames = list(NULL, paste0("v", seq_len(k)))))
dk$f <- rep(c("a", "b"), 10)
dk$y <- rnorm(20)
DBI::dbWriteTable(con, "dk", dk)
fk <- stats::as.formula(paste("y ~ f +", paste(c("f", paste0("v", seq_len(k))), collapse = ":")))
D <- suppressMessages(design_matrix(fk, con, table = "dk"))
mk <- DBI::dbReadTable(con, "dk")
expect_identical(D$sql_layout, "setup")
expect_identical(D$sql_terms, ns$xtx_sql_prepare(D, NULL)$sql_terms)
expect_equal(unname(crossprod(D)), unname(crossprod(model.matrix(fk, mk))))

## FALLBACK: when this version of R's model.matrix.default() is not
## recognised (xtx_prepare() returns NULL), the prototype is a dense matrix
## with no set-up, and the layout is read off the grid; same design
the <- get("the", ns)
prep <- the$prep
the$prep <- NULL
D <- tryCatch(remote(y ~ f * x + poly(z, 2):g, contrasts.arg = list(f = "contr.sum")),
              finally = the$prep <- prep)
expect_identical(D$sql_layout, "grid")
Dn <- remote(y ~ f * x + poly(z, 2):g, contrasts.arg = list(f = "contr.sum"))
expect_identical(Dn$sql_layout, "setup")
expect_identical(D$sql_terms, Dn$sql_terms)
expect_identical(colnames(D), colnames(Dn))
expect_equal(crossprod(D), crossprod(Dn))

test_con_done(con)

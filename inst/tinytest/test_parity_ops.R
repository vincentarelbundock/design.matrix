## Parity of the matrix operations of a design in memory and remote: every
## operation on the same data, in memory and on a DuckDB table and a
## Parquet file (full tier). Two formulas: one that SQL computes, and one with a term
## evaluated in R into columns (also in small batches).
if (!requireNamespace("duckdb", quietly = TRUE) || !requireNamespace("DBI", quietly = TRUE))
    exit_file("duckdb is not installed")
source("helpers.R")
library(xtx)

con <- test_con()
same <- function(a, b, tol = 1e-8) isTRUE(all.equal(unname(a), unname(b), tolerance = tol, check.attributes = FALSE))
set.seed(11)
n <- 400
d <- data.frame(id = seq_len(n), y = rnorm(n), x = rnorm(n), z = runif(n, 1, 2), w = runif(n, 0.5, 2),
                f = sample(c("b", "a", "c"), n, TRUE), g = sample(1:6, n, TRUE))
d$x[3] <- NA
DBI::dbWriteTable(con, "d", d)
d <- DBI::dbReadTable(con, "d")
pq <- tempfile(fileext = ".parquet")
DBI::dbExecute(con, sprintf("COPY d TO '%s' (FORMAT PARQUET)", pq))

formulas <- list(sql = y ~ x * f + log(z), r = y ~ x * f + log(z) + splines::ns(z, 2))

## row-wise results (a remote column, read with the key `id`, or a table
## with the columns of the data), in the order of the rows of M
by_id <- function(tab, col) {
    p <- if (inherits(tab, "xtx_column")) as.data.frame(tab, keep = "id") else as.data.frame(tab)
    p[[col]][match(dk$id, p$id)]
}

sources <- list(
    table = function(f) design_matrix(f, con, table = "d", engine = "auto"),
    batches = function(f) design_matrix(f, con, table = "d", chunk_size = 97),
    parquet = function(f) design_matrix(f, pq))
if (!full_tests) sources$parquet <- NULL
for (fname in names(formulas)) for (src in names(sources)) {
    fo <- formulas[[fname]]
    M <- design_matrix(fo, d)
    X <- as.matrix(M)
    rows <- as.integer(rownames(X))
    dk <- d[rows, ]
    set.seed(3)
    B <- cbind(b1 = rnorm(ncol(X)), b2 = rnorm(ncol(X)))
    V <- crossprod(matrix(rnorm(ncol(X)^2), ncol(X)))
    D <- suppressMessages(sources[[src]](fo))
    info <- paste(fname, src)
    expect_identical(is.null(D$evaluated), fname == "sql", info = info)
    ## shape and printing
    expect_equal(dim(D), dim(M), info = info)
    expect_identical(colnames(D), colnames(M), info = info)
    expect_identical(dimnames(D), list(NULL, colnames(M)), info = info)
    expect_true(any(grepl("remote design_matrix", capture.output(print(D)))), info = info)
    expect_true(any(grepl("remote design_matrix", capture.output(print("w" * D)))), info = info)
    ## as.matrix and sparse coercion (small data; the guard otherwise)
    A <- as.matrix(D)
    expect_true(same(crossprod(A), crossprod(X)) && same(colSums(A), colSums(X)), info = paste(info, "as.matrix"))
    S <- as(D, "CsparseMatrix")
    expect_inherits(S, "CsparseMatrix")
    expect_true(same(as.matrix(crossprod(S)), crossprod(X)), info = paste(info, "sparse"))
    expect_true(same(as.matrix(as(D, "dgCMatrix")), X), info = paste(info, "dgCMatrix"))
    expect_true(same(as.matrix("w" * D), dk$w * X), info = paste(info, "as.matrix scaled"))
    op <- options(xtx.remote_max_rows = 10)
    expect_error(as.matrix(D), "force = TRUE")
    expect_error(as(D, "CsparseMatrix"), "remote_max_rows")
    expect_true(same(colSums(as.matrix(D, force = TRUE)), colSums(X)))
    options(op)
    ## products
    expect_true(same(by_id(D %*% B[, 1], ".fitted"), drop(X %*% B[, 1])), info = paste(info, "%*% vector"))
    pB <- D %*% B
    expect_inherits(pB, "xtx_column")
    expect_equal(dim(pB), c(nrow(D), 2), info = paste(info, "dim of a column"))
    expect_equal(length(pB), 2 * nrow(D), info = paste(info, "length of a column"))
    expect_true(same(cbind(by_id(pB, "b1"), by_id(pB, "b2")), X %*% B), info = paste(info, "%*% matrix"))
    ## crossprod in every form
    expect_true(same(crossprod(D), crossprod(X)), info = paste(info, "crossprod"))
    expect_true(same(crossprod(D, "w" * D), crossprod(X, dk$w * X)), info = paste(info, "X'WX"))
    expect_true(same(crossprod(D, "y"), crossprod(X, dk$y)), info = paste(info, "X'y"))
    expect_true(same(crossprod(D, c("y", "w")), crossprod(X, cbind(dk$y, dk$w))), info = paste(info, "X'[y w]"))
    expect_true(same(crossprod(D, DBI::SQL("y * w")), crossprod(X, dk$y * dk$w)), info = paste(info, "X'(SQL)"))
    expect_true(same(crossprod("w" * D, "y"), crossprod(dk$w * M, dk$y)), info = paste(info, "(WX)'y"))
    expect_true(same(crossprod(2 * D), crossprod(2 * M)), info = paste(info, "(2X)'(2X)"))
    expect_true(same(crossprod(D, D), crossprod(M, M)), info = paste(info, "X'X (two arguments)"))
    expect_true(same(crossprod("w" * D, "z" * D), crossprod(dk$w * M, dk$z * M)), info = paste(info, "(WX)'(ZX)"))
    ## sums and means
    expect_true(same(colSums(D), colSums(X)), info = paste(info, "colSums"))
    expect_true(same(colMeans(D), colMeans(X)), info = paste(info, "colMeans"))
    expect_true(same(colSums("w" * D, group = "g"), colSums(dk$w * M, group = dk$g)), info = paste(info, "colSums group"))
    expect_true(same(colMeans(D, group = "f"), colMeans(M, group = dk$f)), info = paste(info, "colMeans group"))
    expect_true(same(colSums("w" * D), colSums(dk$w * M)), info = paste(info, "colSums weighted"))
    expect_true(same(colMeans("w" * D, group = "g"), colMeans(dk$w * M, group = dk$g)),
                info = paste(info, "colMeans weighted group"))
    expect_true(same(colSums(D, group = DBI::SQL("g % 2")), colSums(M, group = dk$g %% 2)),
                info = paste(info, "colSums SQL group"))
    expect_true(same(by_id(rowSums(D), ".rowsum"), rowSums(X)), info = paste(info, "rowSums"))
    expect_true(same(by_id(rowMeans(D), ".rowmean"), rowMeans(X)), info = paste(info, "rowMeans"))
    expect_true(same(by_id(rowSums("w" * D), ".rowsum"), rowSums(dk$w * X)), info = paste(info, "rowSums scaled"))
    expect_true(same(by_id(rowMeans("w" * D), ".rowmean"), rowMeans(dk$w * X)), info = paste(info, "rowMeans scaled"))
    ## quadratic forms
    expect_true(same(by_id(quadratic_diagonal(D, V), ".quad"), quadratic_diagonal(M, V)), info = paste(info, "quad"))
    expect_true(same(by_id(quadratic_diagonal("w" * D, V), ".quad"), quadratic_diagonal(dk$w * M, V)),
                info = paste(info, "quad scaled"))
    ## arithmetic
    expect_true(same(crossprod(-D), crossprod(-M)), info = paste(info, "-D"))
    expect_true(same(crossprod(D / "w"), crossprod((1 / dk$w) * M)), info = paste(info, "D / w"))
    expect_true(same(crossprod(D * 3), crossprod(M * 3)), info = paste(info, "D * 3"))
    expect_true(same(crossprod(DBI::SQL("w * z") * D), crossprod((dk$w * dk$z) * M)), info = paste(info, "SQL * D"))
    expect_true(same(by_id(("w" * D) %*% B[, 1], ".fitted"), drop((dk$w * M) %*% B[, 1])),
                info = paste(info, "scaled %*%"))
    expect_error(D + D, "not supported")
    expect_error(D - D, "not supported")
    ## rows
    Ds <- D[DBI::SQL("x > 0"), ]
    expect_true(same(crossprod(Ds), crossprod(X[dk$x > 0, ])), info = paste(info, "D[SQL, ]"))
    expect_error(D[1:3, ], "not by position")
    expect_error(D[rep(TRUE, nrow(D)), ], "not by position")
    pos <- xtx:::xtx_p_column(D, "s.x > 0")
    expect_true(same(crossprod(D[pos, ]), crossprod(X[dk$x > 0, ])), info = paste(info, "D[column, ]"))
    expect_equal(nrow(D[pos, ]), sum(dk$x > 0), info = paste(info, "nrow(D[column, ])"))
    fr <- xtx:::finite_rows(D)
    expect_equal(unname(col_sums(col_ifelse(fr, 1, 0))), sum(xtx:::finite_rows(M)),
                 info = paste(info, "finite_rows"))
    ## qr() and its methods
    q <- qr(D)
    qm <- qr(M)
    expect_identical(q$rank, qm$rank, info = paste(info, "qr rank"))
    expect_identical(q$pivot, qm$pivot, info = paste(info, "qr pivot"))
    expect_true(same(abs(qr.R(q)), abs(qr.R(qm)), 1e-7), info = paste(info, "qr.R"))
    expect_true(same(qr.coef(q, "y"), qr.coef(qm, dk$y)), info = paste(info, "qr.coef"))
    expect_true(same(qr.coef(q, c("y", "w")), qr.coef(qm, cbind(dk$y, dk$w))), info = paste(info, "qr.coef matrix"))
    expect_true(same(by_id(qr.fitted(q, "y"), ".fitted"), qr.fitted(qm, dk$y)), info = paste(info, "qr.fitted"))
    expect_true(same(by_id(qr.resid(q, "y"), ".resid"), qr.resid(qm, dk$y)), info = paste(info, "qr.resid"))
    expect_true(same(qr.qty(q, "y"), qr.qty(qm, dk$y)), info = paste(info, "qr.qty"))
    v <- rnorm(q$rank)
    expect_true(same(by_id(qr.qy(q, v), ".qy"), qr.qy(qm, v)), info = paste(info, "qr.qy"))
    QQ <- qr.Q(q)
    Qm <- qr.Q(qm)
    expect_true(same(by_id(QQ, ".Q1"), Qm[, 1]) && same(by_id(QQ, paste0(".Q", q$rank)), Qm[, q$rank]),
                info = paste(info, "qr.Q"))
    expect_true(same(solve(q, "y"), solve(qm, dk$y)), info = paste(info, "solve"))
    expect_error(solve(q), "square")
    ## weighted qr
    qw <- qr("w" * D)
    qwm <- qr(dk$w * M)
    expect_true(same(qr.coef(qw, DBI::SQL("w * y")), qr.coef(qwm, dk$w * dk$y)), info = paste(info, "weighted qr"))
    expect_true(same(by_id(qr.resid(qw, DBI::SQL("w * y")), ".resid"), qr.resid(qwm, dk$w * dk$y)),
                info = paste(info, "weighted qr.resid"))
    expect_true(same(qr.qty(q, c("y", "w")), qr.qty(qm, cbind(dk$y, dk$w))), info = paste(info, "qr.qty matrix"))
    qf <- qr.fitted(q, c("y", "w"))
    expect_true(same(cbind(by_id(qf, ".fitted1"), by_id(qf, ".fitted2")), qr.fitted(qm, cbind(dk$y, dk$w))),
                info = paste(info, "qr.fitted matrix"))
}

test_con_done(con)
unlink(pq)

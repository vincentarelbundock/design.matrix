# Shared helpers for the tests (sourced by each test file).
#
# Test tiers. By default (R CMD check, CRAN, quick iterations) the suite
# runs a fast tier: every feature in memory and on a DuckDB table (the SQL
# backend). DESIGN_MATRIX_TESTS=full runs the full tier: also on Parquet
# files, every case of the parity loops remotely (both compression modes,
# all fixed-effect cases), the random-formula sweeps at full length and the
# slow regression tests. Run it before every merge:
#   DESIGN_MATRIX_TESTS=full Rscript -e 'tinytest::run_test_dir("inst/tinytest")'
# Files shorten their loops with `full_tests` and `remote_modes`; a file
# that only belongs to the full tier starts with
#   if (!full_tests) exit_file("full tier only: set DESIGN_MATRIX_TESTS=full")
full_tests <- identical(Sys.getenv("DESIGN_MATRIX_TESTS"), "full")
# remote data sources a test runs on: a DuckDB table, plus Parquet in full
remote_modes <- if (full_tests) c("sql", "parquet") else "sql"

strict_identical <- function(a, b) {
    identical(a, b, num.eq = FALSE, single.NA = FALSE, attrib.as.set = FALSE)
}

# Run f(), collecting warnings and the error message, if any.
capture <- function(f) {
    w <- character()
    v <- withCallingHandlers(
        tryCatch(f(), error = function(e) structure(list(msg = conditionMessage(e)), class = "captured_error")),
        warning = function(cond) {
            w <<- c(w, conditionMessage(cond))
            invokeRestart("muffleWarning")
        })
    list(value = v, warnings = w)
}

same_outcome <- function(a, b) {
    if (inherits(a$value, "captured_error") || inherits(b$value, "captured_error")) {
        return(inherits(a$value, "captured_error") && inherits(b$value, "captured_error") &&
               identical(a$value$msg, b$value$msg) && identical(a$warnings, b$warnings))
    }
    strict_identical(a$value, b$value) && identical(a$warnings, b$warnings)
}

# Random data with awkward values.
fuzz_data <- function(n, special = TRUE, na = FALSE) {
    num <- function() {
        x <- rnorm(n)
        if (special && n >= 6) x[sample(n, 4)] <- sample(c(-0, Inf, -Inf, 0, 1e300, -1e-300), 4)
        if (na && n >= 3) x[sample(n, 1)] <- sample(c(NA, NaN), 1)
        x
    }
    fac <- function(k) {
        lev <- letters[seq_len(k)]
        x <- sample(lev, n, TRUE)
        ## level order that is not alphabetical, and sometimes a level
        ## that never occurs (the bugs found in Rfast and mdl)
        if (stats::runif(1) < 0.5) lev <- sample(lev)
        if (stats::runif(1) < 0.3) lev <- append(lev, "zz", after = sample(0:k, 1))
        f <- factor(x, levels = lev)
        if (na && n >= 3) f[sample(n, 1)] <- NA
        f
    }
    d <- data.frame(
        y = rnorm(n),
        x1 = num(),
        x2 = num(),
        i1 = sample(-3:3, n, TRUE),
        f1 = fac(3),
        f2 = fac(4),
        o1 = factor(sample(c("lo", "mid", "hi"), n, TRUE), levels = c("lo", "mid", "hi"), ordered = TRUE),
        ch = sample(c("u", "v", "w"), n, TRUE),
        lg = sample(c(TRUE, FALSE), n, TRUE),
        stringsAsFactors = FALSE)
    if (na && n >= 3) d$lg[sample(n, 1)] <- NA
    d$m <- I(matrix(rnorm(2 * n), n, 2, dimnames = list(NULL, c("a", "b"))))
    d$z <- runif(n, 1, 2)
    if (stats::runif(1) < 0.3) rownames(d) <- paste0("r", sample(1e6, n))
    d
}

fuzz_terms <- c("x1", "x2", "i1", "f1", "f2", "o1", "ch", "lg", "m", "poly(z, 2)",
                "I(x2^2)", "log(z)", "splines::ns(z, 3)", "y")

fuzz_formula <- function(response = TRUE) {
    k <- sample(1:3, 1)
    tm <- sample(fuzz_terms, k)
    ops <- sample(c(" + ", " * ", ":", " + ", " %in% "), max(0, k - 1), TRUE)
    rhs <- tm[1]
    if (k > 1) for (j in 2:k) rhs <- paste0(rhs, ops[j - 1], tm[j])
    if (stats::runif(1) < 0.15) rhs <- paste0("(", rhs, ")^2")
    if (stats::runif(1) < 0.25) rhs <- paste0("0 + ", rhs)
    if (stats::runif(1) < 0.05) rhs <- sample(c("1", "0 + f1", "-1 + x1"), 1)
    stats::as.formula(paste(if (response) "y" else "", "~", rhs), env = globalenv())
}

fuzz_contrasts <- function(d) {
    if (stats::runif(1) < 0.6) return(NULL)
    pool <- list(f1 = "contr.sum", f2 = sample(c("contr.helmert", "contr.sum"), 1),
                 o1 = "contr.treatment", ch = "contr.SAS", lg = "contr.sum")
    keep <- sample(names(pool), sample(1:3, 1))
    out <- pool[keep]
    if ("f2" %in% keep && stats::runif(1) < 0.5) out$f2 <- stats::contr.poly(nlevels(d$f2))
    out
}

# Remote fixture of the fixed-effects suites (test_fe.R, test_fe_glm.R):
# write `d`, with a `.row` column (its row order), to table `name` of the
# DuckDB connection `con` and to `name`.parquet in directory `dir`; return
# the two sources, list(table, parquet).
write_fe_source <- function(con, dir, d, name) {
    d$.row <- seq_len(nrow(d))
    DBI::dbWriteTable(con, name, d, overwrite = TRUE)
    path <- file.path(dir, paste0(name, ".parquet"))
    DBI::dbExecute(con, sprintf("COPY %s TO '%s' (FORMAT PARQUET)", name, path))
    list(table = name, parquet = path)
}

# Count the SQL statements an expression sends to the database: every call
# of DBI::dbGetQuery, dbExecute, dbSendQuery, dbSendStatement, dbWriteTable,
# dbAppendTable and duckdb::duckdb_register, wherever it comes from (unlike
# the option design.matrix.remote_debug, which only sees dm_query()). A call
# made while another traced call runs (a DBI method built on another
# generic) is not counted again. Returns list(value, count); the functions
# are untraced on exit, even on error.
count_queries <- function(expr) {
    state <- new.env()
    state$count <- 0L
    state$depth <- 0L
    enter <- function() {
        if (state$depth == 0L) state$count <- state$count + 1L
        state$depth <- state$depth + 1L
    }
    leave <- function() state$depth <- state$depth - 1L
    fns <- list(DBI = c("dbGetQuery", "dbExecute", "dbSendQuery", "dbSendStatement",
                        "dbWriteTable", "dbAppendTable"),
                duckdb = "duckdb_register")
    traced <- list()
    on.exit(for (t in traced) suppressMessages(untrace(t[[2]], where = asNamespace(t[[1]]))))
    for (pkg in names(fns)) for (f in fns[[pkg]]) {
        suppressMessages(trace(f, tracer = as.call(list(enter)), exit = as.call(list(leave)), print = FALSE,
                               where = asNamespace(pkg)))
        traced[[length(traced) + 1L]] <- list(pkg, f)
    }
    value <- expr
    list(value = value, count = state$count)
}

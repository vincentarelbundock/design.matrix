# Shared helpers for the tests (sourced by each test file).
#
# Test tiers. By default (R CMD check, CRAN, quick iterations) the suite
# runs a fast tier: every feature in memory and on a DuckDB table (the SQL
# backend). XTX_TESTS=full runs the full tier: also on Parquet
# files, every case of the parity loops remotely (both compression modes,
# all fixed-effect cases), the random-formula sweeps at full length and the
# slow regression tests. Run it before every merge, the files in parallel
# (4 worker processes, each limited to a share of the threads):
#   XTX_TESTS=full Rscript -e 'cl <- parallel::makeCluster(4); invisible(parallel::clusterEvalQ(cl, { library(xtx); options(xtx.threads = max(1L, parallel::detectCores() %/% 4L)) })); tinytest::run_test_dir(normalizePath("inst/tinytest"), cluster = cl)'
# or one file after another (slower):
#   XTX_TESTS=full Rscript -e 'tinytest::run_test_dir("inst/tinytest")'
# tests/tinytest.R (R CMD check) runs the files in parallel the same way;
# see there for the number of processes (XTX_TEST_NCPU).
# Files shorten their loops with `full_tests` and `remote_modes`; a file
# that only belongs to the full tier starts with
#   if (!full_tests) exit_file("full tier only: set XTX_TESTS=full")
full_tests <- identical(Sys.getenv("XTX_TESTS"), "full")
# remote data sources a test runs on: a DuckDB table, plus Parquet in full
remote_modes <- if (full_tests) c("sql", "parquet") else "sql"

# Shared DuckDB fixture. test_con() returns the one in-memory DuckDB
# connection of this R process (opened on first use, closed when the process
# exits), switched to a new, empty database (ATTACH ':memory:'; USE), so the
# tables a file writes (d, nd, ...) never meet those of another file run
# earlier in the same process. Each call starts a new database and drops the
# previous one, also when a file stopped before test_con_done(con), which
# switches back and drops it, with the temporary tables (the package's, named
# per process, in the connection's temp catalog). A
# file that needs a connection of its own (test_query_counts.R counts every
# statement) opens one with DBI::dbConnect().
if (!exists(".xtx_test_db", envir = globalenv(), inherits = FALSE)) {
    assign(".xtx_test_db", new.env(), envir = globalenv())
    reg.finalizer(get(".xtx_test_db", envir = globalenv()), function(e) {
        if (!is.null(e$con) && DBI::dbIsValid(e$con)) DBI::dbDisconnect(e$con)
    }, onexit = TRUE)
}
test_con <- function() {
    st <- get(".xtx_test_db", envir = globalenv())
    if (is.null(st$con) || !DBI::dbIsValid(st$con)) {
        st$con <- DBI::dbConnect(duckdb::duckdb())
        st$home <- DBI::dbGetQuery(st$con, "SELECT current_database() AS db")$db
        st$db <- NULL
        st$i <- 0L
    }
    test_con_done(st$con)
    st$i <- st$i + 1L
    st$db <- sprintf("xtx_test_%d", st$i)
    DBI::dbExecute(st$con, sprintf("ATTACH ':memory:' AS %s", st$db))
    DBI::dbExecute(st$con, sprintf("USE %s", st$db))
    st$con
}
test_con_done <- function(con) {
    st <- get(".xtx_test_db", envir = globalenv())
    if (is.null(st$db)) return(invisible())
    DBI::dbExecute(con, sprintf("USE %s", DBI::dbQuoteIdentifier(con, st$home)))
    DBI::dbExecute(con, sprintf("DETACH DATABASE IF EXISTS %s", st$db))
    st$db <- NULL
    ## temporary tables hold data: drop the file's (temporary views are only
    ## queries; data frames registered with duckdb_register() stay)
    tmp <- DBI::dbGetQuery(con, "SELECT table_name FROM duckdb_tables() WHERE temporary")$table_name
    for (tn in tmp)
        DBI::dbExecute(con, sprintf("DROP TABLE IF EXISTS temp.main.%s", DBI::dbQuoteIdentifier(con, tn)))
    invisible()
}

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


# What an expression sends to the database and what comes back: every call
# of DBI::dbGetQuery, dbExecute, dbSendQuery, dbSendStatement, dbWriteTable,
# dbAppendTable and duckdb::duckdb_register, wherever it comes from (the
# package's xtx_query() and xtx_execute() end in the first two; the option
# xtx.remote_debug only sees those). A call made while another traced call
# runs (a DBI method built on another generic) is not recorded again.
# Returns list(value, count = statements, log = one row per statement: `fn`,
# `sql` (the statement, or the name of the table written), `rows` and `cols`
# of the result fetched (0 for a statement that returns none) and
# `up_rows`, `up_cols` of the data sent (dbWriteTable, dbAppendTable,
# duckdb_register), fetched = sum(rows * cols), max_rows = the most rows
# of one result). The functions are untraced on exit, even on error.
record_queries <- function(expr) {
    state <- new.env()
    state$depth <- 0L
    state$log <- list()
    enter <- function() state$depth <- state$depth + 1L
    leave <- function(fn, frame, value) {
        state$depth <- state$depth - 1L
        if (state$depth > 0L) return(invisible())
        arg <- function(nm) tryCatch(get(nm, envir = frame, inherits = FALSE), error = function(e) NULL)
        sql <- arg("statement") %||% arg("name")
        sql <- if (is.character(sql) || inherits(sql, "SQL")) paste(as.character(sql), collapse = " ") else NA_character_
        up <- if (fn %in% c("dbWriteTable", "dbAppendTable")) arg("value") else if (fn == "duckdb_register") arg("df")
        got <- if (fn == "dbGetQuery" && is.data.frame(value)) dim(value) else c(0L, 0L)
        state$log[[length(state$log) + 1L]] <- data.frame(
            fn = fn, sql = sql, rows = got[1L], cols = got[2L],
            up_rows = if (is.data.frame(up)) nrow(up) else 0L, up_cols = if (is.data.frame(up)) ncol(up) else 0L,
            stringsAsFactors = FALSE)
        invisible()
    }
    `%||%` <- function(a, b) if (is.null(a)) b else a
    fns <- list(DBI = c("dbGetQuery", "dbExecute", "dbSendQuery", "dbSendStatement",
                        "dbWriteTable", "dbAppendTable"),
                duckdb = "duckdb_register")
    traced <- list()
    on.exit(for (t in traced) suppressMessages(untrace(t[[2]], where = asNamespace(t[[1]]))))
    for (pkg in names(fns)) for (f in fns[[pkg]]) {
        ## (the exit call is evaluated in the frame of the traced function:
        ## its arguments and its value)
        suppressMessages(trace(f, tracer = as.call(list(enter)),
                               exit = as.call(list(leave, f, quote(environment()), quote(returnValue()))),
                               print = FALSE, where = asNamespace(pkg)))
        traced[[length(traced) + 1L]] <- list(pkg, f)
    }
    value <- expr
    log <- if (length(state$log)) do.call(rbind, state$log)
           else data.frame(fn = character(), sql = character(), rows = integer(), cols = integer(),
                           up_rows = integer(), up_cols = integer(), stringsAsFactors = FALSE)
    list(value = value, count = nrow(log), log = log, fetched = sum(as.numeric(log$rows) * log$cols),
         max_rows = if (nrow(log)) max(log$rows) else 0L)
}

# The plan of a formula and training data, or of a fitted model (internal
# to xtx; users reach plans through design_matrix(, newdata = )).
plan_of <- function(...) get("xtx_plan_of", asNamespace("xtx"))(...)

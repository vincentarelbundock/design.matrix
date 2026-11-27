## Column arithmetic for package developers (Kata design.matrix#tapx, tier
## T1 of the pg3k decision): base generics on the public column type
## (Summary group, mean, the Math cases round/signif/trunc/gamma/lgamma
## that primitives.R does not cover) and four exported helpers,
## col_ifelse(), col_pmin(), col_pmax() and col_sums(), thin wrappers of
## xtx_c_ifelse(), xtx_c_pmin(), xtx_c_pmax() and the sums of xtx_colsum().
##
## In memory a column is a plain numeric vector, so base R applies and the
## helpers call base R (or the xtx_c_ function, which in memory computes as
## base R does). Remotely every result has R's semantics, with one rule for
## missing values (p701, 2tr6): SQL NULL is NA, IEEE NaN is NaN.
##   - sum, prod, mean: NA if a value is NA, else NaN if one is NaN (where
##     R's result for a mix of NA and NaN depends on the platform and the
##     order, the remote one is NA).
##   - min, max, range: NA if a value is NA, else NaN if one is NaN (R:
##     "any NA trumps all NaNs"); no value left gives Inf / -Inf with R's
##     warning.
##   - any, all: three-valued, a NaN counts as NA.
##   - na.rm = TRUE drops NA and NaN, as base R does.
## Remote results are doubles (logical for any/all): a logical column sums
## to a double, where base R gives an integer.
##
## Every summary is one query, whatever the number of arguments: the
## columns are bound together (xtx_c_cbind(), whose frame takes the pending
## products and lookups of all of them) and the aggregates of each SQL
## expression are read in one SELECT, as xtx_colsum() reads its sums.

## ---- aggregates of columns: one query ----

## The aggregates `what` of every SQL expression of the remote columns in
## `cols` (a list; same design), from one query: a data frame with one row
## per expression and columns na (NULL values), nan (NaN values), and as
## asked: sum, prod, min, max (over the values that are neither; 0, 1,
## NA, NA over no value), fmin, fmax (over the finite values), n (values
## that are neither), rows (all rows), ntrue, nfalse (values that are
## neither, nonzero / zero).
xtx_column_aggregates <- function(cols, what) {
    x <- do.call(xtx_c_cbind, cols)
    fr <- xtx_frame_bind(unclass(x)$frame)
    sql <- unclass(x)$sql
    k <- length(sql)
    nan <- "CAST('NaN' AS DOUBLE)"
    d <- sprintf("CAST(%s AS DOUBLE)", sql)
    v <- sprintf("NULLIF(%s, %s)", d, nan)
    agg <- list(
        na = sprintf("SUM(CASE WHEN %s IS NULL THEN 1 ELSE 0 END)", sql),
        nan = sprintf("SUM(CASE WHEN isnan(%s) THEN 1 ELSE 0 END)", d),
        sum = sprintf("COALESCE(SUM(%s), 0)", v),
        prod = sprintf("COALESCE(PRODUCT(%s), 1)", v),
        min = sprintf("MIN(%s)", v),
        max = sprintf("MAX(%s)", v),
        fmin = sprintf("MIN(CASE WHEN isfinite(%1$s) THEN %1$s END)", d),
        fmax = sprintf("MAX(CASE WHEN isfinite(%1$s) THEN %1$s END)", d),
        n = sprintf("COUNT(%s)", v),
        rows = rep("COUNT(*)", k),
        ntrue = sprintf("SUM(CASE WHEN %s <> 0 THEN 1 ELSE 0 END)", v),
        nfalse = sprintf("SUM(CASE WHEN %s = 0 THEN 1 ELSE 0 END)", v))
    what <- unique(c("na", "nan", what))
    items <- unlist(lapply(what, function(a) sprintf("%s AS %s_%d", agg[[a]], a, seq_len(k))))
    res <- xtx_query(fr, sprintf("SELECT %s FROM (%s) AS s WHERE %s", paste(items, collapse = ", "),
                                fr$src$select, fr$filter))
    out <- lapply(what, function(a) vapply(seq_len(k), function(j) as.double(res[[sprintf("%s_%d", a, j)]]), 0))
    names(out) <- what
    as.data.frame(out)
}

## The same aggregates of a value in memory (a number or vector passed
## with remote columns, as in sum(col, 1)): one row.
xtx_value_aggregates <- function(x, what) {
    if (!is.numeric(x) && !is.logical(x))
        stop("Remote columns combine only with numbers, logical values and other remote columns.", call. = FALSE)
    x <- as.double(x)
    miss <- is.na(x)
    ok <- x[!miss]
    fin <- x[is.finite(x)]
    out <- list(na = sum(miss & !is.nan(x)), nan = sum(is.nan(x)),
                sum = sum(ok), prod = prod(ok),
                min = if (length(ok)) min(ok) else NA_real_,
                max = if (length(ok)) max(ok) else NA_real_,
                fmin = if (length(fin)) min(fin) else NA_real_,
                fmax = if (length(fin)) max(fin) else NA_real_,
                n = length(ok), rows = length(x), ntrue = sum(ok != 0), nfalse = sum(ok == 0))
    as.data.frame(lapply(out[unique(c("na", "nan", what))], as.double))
}

## The aggregates of the arguments of a summary: remote columns in one
## query, other values in R; one row per expression or value.
xtx_summary_aggregates <- function(args, what) {
    remote <- vapply(args, xtx_is_column, NA)
    parts <- list()
    if (any(remote)) parts[[1L]] <- xtx_column_aggregates(args[remote], what)
    for (a in args[!remote]) parts[[length(parts) + 1L]] <- xtx_value_aggregates(a, what)
    do.call(rbind, parts)
}

## ---- Summary group: sum, prod, min, max, range, any, all ----

#' Summaries and arithmetic of columns
#'
#' Base R functions work on the columns of a design (see [design_column]):
#' in memory a column is a numeric vector, and remotely a `design_column`
#' has methods for the same functions, which compute in the database and
#' return the same values as base R on the vector. Together with
#' `col_ifelse()`, `col_pmin()`, `col_pmax()` and `col_sums()`, which
#' cover what base R cannot dispatch on, the same code runs on a data frame
#' and on a database table.
#'
#' @section Base generics on remote columns:
#' * Arithmetic and comparisons (`+`, `-`, `*`, `/`, `^`, `==`, `<`, ...,
#'   `&`, `|`, `!`), [is.na()], [is.nan()], [is.finite()], [is.infinite()]
#'   and the `Math` functions `exp()`, `log()` (`log(x, base)`: `base` a
#'   number or a column of the design, R's `log(x, 1)` and negative bases
#'   included), `sqrt()`,
#'   `abs()`, `log1p()`, `expm1()`, `log2()`, `log10()`, `cos()`, `sin()`,
#'   `tan()`, `floor()`, `ceiling()`, `sign()`, `trunc()`, `round(x,
#'   digits)` (R's algorithm: `round(2.5)` is 2, `round(0.15, 1)` is 0.1),
#'   `signif(x, digits)`, `gamma()` and `lgamma()` build a new column,
#'   lazily (no query).
#' * `sum()`, `prod()`, `min()`, `max()`, `range()` (with `finite`),
#'   `any()`, `all()` (the `Summary` group) and `mean()` return numbers (or
#'   `TRUE`/`FALSE`/`NA`), over all the values of all their arguments, as
#'   base R (a column with several columns counts all its values; numbers
#'   among the arguments count too). Each call is one query, whatever the
#'   number of arguments. `na.rm = TRUE` drops `NA` and `NaN`.
#'
#' @section Missing values:
#' SQL `NULL` is `NA` and IEEE `NaN` is `NaN`, with R's rules: arithmetic
#' propagates both, comparisons with either give `NA`, `is.na()` is `TRUE`
#' for both and `is.nan()` for `NaN` only. `sum()`, `prod()` and `mean()`
#' give `NA` when a value is `NA`, else `NaN` when one is `NaN` (where R's
#' result for a mix of `NA` and `NaN` depends on the platform, the remote
#' one is `NA`). `min()`, `max()` and `range()` give `NA` when a value is
#' `NA`, else `NaN` (as base R: an `NA` trumps a `NaN`), and with no value
#' left, `Inf` or `-Inf` with base R's warning. `any()` and `all()` are
#' three-valued, a `NaN` counting as `NA`. `col_ifelse()` gives `NA` where
#' the test is `NA`; `col_pmin()` and `col_pmax()` give `NA` or `NaN`
#' where either value is, as `pmin()` and `pmax()`.
#'
#' Remote values are doubles: a logical remote column (a comparison) sums
#' to a double, where base R gives an integer, and `as.vector()` returns it
#' as 0 and 1. Sums over many rows can differ from R's in the last bits
#' (the database adds in another order, and R in extended precision);
#' `gamma()` and `lgamma()` use the C library, R its own code (equal to
#' rounding).
#'
#' @section Helpers:
#' * `col_ifelse(test, yes, no)`: `yes` where `test` is `TRUE`, `no` where
#'   it is `FALSE`, `NA` where it is `NA`; `yes` and `no` are columns or
#'   numbers. In memory as [ifelse()], but always a double vector, without
#'   the attributes of `test`.
#' * `col_pmin(a, b)`, `col_pmax(a, b)`: the parallel minimum and maximum of
#'   two columns, or of a column and a number, as [pmin()] and [pmax()].
#' * `col_sums(...)`: the sums of several columns (each over all its
#'   values, as `sum()`), in one query remotely: a named double vector, one
#'   value per argument, named after the argument (`col_sums(n = 1 + 0 *
#'   r, rss = r^2)`), or after the variable for a bare name. `na.rm` can
#'   differ by argument: `col_sums(dev, col_ifelse(cond, 1, 0), na.rm =
#'   c(FALSE, TRUE))` counts the rows where `cond` is `TRUE` (`NA` counting as
#'   `FALSE`) in the query that sums `dev`.
#' * `col_const(like, value)`: a constant column with the rows of `like` (a
#'   column or a design): `rep(value, n)` in memory (names may change before
#'   the first release).
#' * `cbind(...)`: columns side by side, one column with several columns
#'   (numbers recycle), as `cbind()` of vectors in memory; `x[, j]` selects
#'   them again.
#' * `rowSums(x)`: the sums across the columns of a column with several
#'   columns, a column.
#' * `col_compute(x, name = NULL)` (see [design_column]): the column stored
#'   once in a table, a temporary one or, with a `name`, a permanent table
#'   of the connection, and a column that reads it, for a column that later
#'   queries read several times. In memory it returns `x`.
#'
#' Arguments can mix remote columns and numbers. Remote columns must come
#' from the same design (columns of two designs are an error); in memory
#' the vectors recycle as in base R. The `Summary` functions dispatch on
#' their first argument: put a column first (`sum(r, 1)`, not `sum(1, r)`).
#'
#' @param test A logical column (in memory a logical vector).
#' @param missing For `col_ifelse()`: `NULL` (default: a missing `test` gives a
#'   missing value, as [ifelse()]), or one number, the value where `test` is
#'   missing.
#' @param yes,no,a,b Columns (vectors in memory, `design_column`s remotely)
#'   or numbers.
#' @param ... For `col_sums()`, columns or numbers, optionally named. For
#'   `cbind()`, columns or numbers. For the methods, as in base R.
#' @param checked For `col_pmin()` and `col_pmax()`: `TRUE` (default)
#'   follows [pmin()] and [pmax()] on missing values (`NA` or `NaN` if either
#'   operand is). `FALSE` is for operands known to have no missing values:
#'   remotely plain `least()` and `greatest()`, whose SQL is shorter; in
#'   memory the same values.
#' @param like A column, or a design: the rows of the constant column.
#' @param value One number.
#' @param deparse.level Not used.
#' @param na.rm For `col_sums()`, drop `NA` and `NaN` values (as `sum()`):
#'   `TRUE` or `FALSE`, or one value per argument (names may change before
#'   the first release). For `rowSums()`, only the default.
#' @param dims For `rowSums()`, only the default.
#' @param x A remote column.
#' @param trim Only `0` is supported for a remote column.
#' @return `col_ifelse()`, `col_pmin()`, `col_pmax()`: a column (a double
#'   vector in memory, a `design_column` remotely). `col_sums()`: a double
#'   vector. The `Summary` methods and `mean()`: a number, two for
#'   `range()`, a logical value for `any()` and `all()`.
#' @examples
#' d <- data.frame(y = c(1, 4, 6, 2, 8, 3), x = c(0, 1, 1, 2, 3, 3))
#' D <- design_matrix(y ~ x, d)
#' y <- design_column(D, "y")
#' z <- col_ifelse(y > 5, NA, y)
#' sum(z)
#' sum(z, na.rm = TRUE)
#' range(z, na.rm = TRUE)
#' mean(round(y / 4, 1))
#' col_sums(n = !is.na(z), total = col_ifelse(is.na(z), 0, z))
#' col_pmax(y, 3)
#'
#' if (requireNamespace("duckdb", quietly = TRUE)) {
#'     con <- DBI::dbConnect(duckdb::duckdb())
#'     DBI::dbWriteTable(con, "d", d)
#'     Dr <- design_matrix(y ~ x, con, table = "d")
#'     yr <- design_column(Dr, "y")
#'     zr <- col_ifelse(yr > 5, NA, yr)
#'     sum(zr)
#'     sum(zr, na.rm = TRUE)
#'     range(zr, na.rm = TRUE)
#'     mean(round(yr / 4, 1))
#'     col_sums(n = !is.na(zr), total = col_ifelse(is.na(zr), 0, zr))
#'     sort(as.vector(col_pmax(yr, 3)))
#'     DBI::dbDisconnect(con)
#' }
#' @name column_arithmetic
NULL

#' @rdname column_arithmetic
#' @export
Summary.xtx_column <- function(..., na.rm = FALSE) {
    gen <- .Generic
    args <- list(...)
    finite <- FALSE
    if (gen == "range" && "finite" %in% names(args)) {
        finite <- isTRUE(args$finite)
        args$finite <- NULL
    }
    if (!is.null(names(args)) && any(names(args) != ""))
        stop(sprintf("%s() of remote columns: unknown argument `%s`.", gen, names(args)[names(args) != ""][1L]),
             call. = FALSE)
    na.rm <- isTRUE(na.rm) || finite
    what <- switch(gen, sum = "sum", prod = "prod", min = "min", max = "max",
                   range = if (finite) c("fmin", "fmax") else c("min", "max"),
                   any = , all = c("ntrue", "nfalse"))
    s <- xtx_summary_aggregates(args, what)
    na <- sum(s$na) > 0
    nan <- sum(s$nan) > 0
    miss <- !na.rm && (na || nan)
    ## (R: NA if a value is NA, else NaN if one is NaN)
    missing_value <- if (na) NA_real_ else NaN
    extreme <- function(v, f, empty) {
        if (miss) return(missing_value)
        v <- v[!is.na(v)]
        if (!length(v)) {
            warning(sprintf("no non-missing arguments to %s; returning %s", f, if (empty > 0) "Inf" else "-Inf"),
                    call. = FALSE)
            return(empty)
        }
        if (f == "min") min(v) else max(v)
    }
    switch(gen,
        ## (not na.rm: a sum of Inf and -Inf is a NaN to keep, as in R)
        sum = if (miss) missing_value else sum(s$sum),
        prod = if (miss) missing_value else prod(s$prod),
        min = extreme(s$min, "min", Inf),
        max = extreme(s$max, "max", -Inf),
        range = if (finite) c(extreme(s$fmin, "min", Inf), extreme(s$fmax, "max", -Inf))
                else c(extreme(s$min, "min", Inf), extreme(s$max, "max", -Inf)),
        any = if (sum(s$ntrue) > 0) TRUE else if (miss) NA else FALSE,
        all = if (sum(s$nfalse) > 0) FALSE else if (miss) NA else TRUE)
}

#' @rdname column_arithmetic
#' @export
mean.xtx_column <- function(x, trim = 0, na.rm = FALSE, ...) {
    if (!identical(trim, 0) && !identical(trim, 0L))
        stop("mean() of a remote column: only `trim = 0` is supported.", call. = FALSE)
    s <- xtx_column_aggregates(list(x), c("sum", "n", "rows"))
    if (!isTRUE(na.rm)) {
        if (sum(s$na) > 0) return(NA_real_)
        if (sum(s$nan) > 0) return(NaN)
    }
    ## (an empty mean is NaN, as in R: 0 / 0; Inf and -Inf give NaN)
    total <- sum(s$sum)
    total / (if (isTRUE(na.rm)) sum(s$n) else sum(s$rows))
}

## ---- Math cases beyond primitives.R: round, signif, trunc, gamma, lgamma ----

## R_pow_di(10, n) of R's nmath, as it computes it (repeated squaring), so
## that the constants of round() and signif() are those of R.
xtx_pow_di <- function(x, n) {
    xn <- 1
    if (n == 0) return(xn)
    neg <- n < 0
    n <- abs(n)
    repeat {
        if (n %% 2 == 1) xn <- xn * x
        n <- n %/% 2
        if (n > 0) x <- x * x else break
    }
    if (neg) 1 / xn else xn
}

## The DuckDB macros of round(x, digits) and signif(x, digits), created
## once per connection: R's fround() and fprec() (src/nmath), step by step
## in double arithmetic. xtx_pow10(m) is R_pow_di(10, m) from a table of the
## values R computes (pow() can differ in the last bit for m > 22).
xtx_sql_round_macros <- function(fr) {
    if (isTRUE(fr$src$holder$xtx_round_macros)) return(invisible())
    tab <- paste(xtx_sql_double(vapply(0:400, function(m) xtx_pow_di(10, m), 0)), collapse = ", ")
    pow10 <- sprintf("CASE WHEN m >= 0 THEN ([%1$s])[CAST(m AS BIGINT) + 1] ELSE 1 / ([%1$s])[1 - CAST(m AS BIGINT)] END", tab)
    ## (M_LOG10_2; a string cast, correctly rounded, where a decimal literal
    ## might not be)
    log10_2 <- "CAST('0.301029995663981195213738894724493027' AS DOUBLE)"
    ## fround(): x10 = (|x| pow10) p10, then the nearer of floor(x10) and
    ## ceiling(x10) over pow10 and p10 (the odd one up on a tie); p10 = 1
    ## when dig <= 308. No rounding when dig + log10(|x|) > 15.
    x10 <- "((abs(x) * pw) * p1)"
    xd <- sprintf("(floor(%s) / pw / p1)", x10)
    xu <- sprintf("(ceil(%s) / pw / p1)", x10)
    g <- "(CASE WHEN isfinite(x) AND x <> 0 THEN abs(x) ELSE 1 END)"
    rnd <- sprintf(paste("CASE WHEN x IS NULL OR isnan(x) OR isinf(x) OR x = 0 THEN x",
                         "WHEN %1$s * (0.5 + floor(log2(%2$s))) + dig > 15 THEN x",
                         "ELSE sign(x) * (CASE WHEN (%4$s - abs(x)) < (abs(x) - %3$s)",
                         "OR (fmod(floor(%5$s), 2) = 1 AND (%4$s - abs(x)) = (abs(x) - %3$s)) THEN %4$s ELSE %3$s END) END"),
                   log10_2, g, xd, xu, x10)
    ## fprec(): l10 = log10(|x|), e10 = dig - 1 - floor(l10)
    la <- sprintf("log10(%s)", g)
    e10 <- sprintf("CAST(dig - 1 - floor(%s) AS BIGINT)", la)
    e2 <- sprintf("(dig + CASE WHEN %s > 0 THEN 22 ELSE -22 END)", e10)
    sig <- sprintf(paste("CASE WHEN x IS NULL OR isnan(x) OR isinf(x) OR x = 0 THEN x",
                         "ELSE sign(x) * (CASE WHEN abs(%1$s) < 306 THEN (CASE",
                         "WHEN %2$s > 308 THEN (round_even((abs(x) * xtx_pow10(308)) * xtx_pow10(%2$s - 308), 0) / xtx_pow10(308)) / xtx_pow10(%2$s - 308)",
                         "WHEN %2$s > 0 THEN round_even(abs(x) * xtx_pow10(%2$s), 0) / xtx_pow10(%2$s)",
                         "ELSE round_even(abs(x) / xtx_pow10(-%2$s), 0) * xtx_pow10(-%2$s) END)",
                         "ELSE floor((abs(x) * xtx_pow10(%3$s)) * xtx_pow10(%2$s - %3$s)",
                         "+ CASE WHEN lmax - %1$s >= thr THEN CAST(0.5 AS DOUBLE) ELSE CAST(0 AS DOUBLE) END) / xtx_pow10(%3$s) / xtx_pow10(%2$s - %3$s) END) END"),
                   la, e10, e2)
    xtx_execute(fr, sprintf("CREATE OR REPLACE TEMP MACRO xtx_pow10(m) AS (%s)", pow10))
    xtx_execute(fr, sprintf("CREATE OR REPLACE TEMP MACRO xtx_round_digits(x, dig, pw, p1) AS (%s)", rnd))
    xtx_execute(fr, sprintf("CREATE OR REPLACE TEMP MACRO xtx_signif(x, dig, thr, lmax) AS (%s)", sig))
    fr$src$holder$xtx_round_macros <- TRUE
    invisible()
}

## round(x, digits) of a column's SQL `a` (constant digits), as R's fround().
xtx_sql_round <- function(fr, a, digits) {
    x <- sprintf("CAST(%s AS DOUBLE)", a)
    if (digits > 323) return(x)
    if (digits < -308) return(sprintf("(CASE WHEN %1$s IS NULL OR isnan(%1$s) OR isinf(%1$s) OR %1$s = 0 THEN %1$s ELSE CAST(0 AS DOUBLE) END)", x))
    ## (round half to even, as nearbyint())
    if (digits == 0) return(sprintf("round_even(%s, 0)", x))
    dig <- as.integer(floor(digits + 0.5))
    if (dig <= 308) {
        pw <- xtx_pow_di(10, dig)
        p1 <- 1
    } else {
        pw <- xtx_pow_di(10, 308)
        p1 <- xtx_pow_di(10, dig - 308)
    }
    xtx_sql_round_macros(fr)
    sprintf("xtx_round_digits(%s, %d, %s, %s)", x, dig, xtx_sql_double(pw), xtx_sql_double(p1))
}

## signif(x, digits) of a column's SQL `a` (constant digits), as R's fprec().
xtx_sql_signif <- function(fr, a, digits) {
    x <- sprintf("CAST(%s AS DOUBLE)", a)
    if (!is.finite(digits)) {
        if (digits > 0) return(x)
        digits <- 1
    }
    ## ((int)round(digits) in C: half away from zero)
    dig <- as.integer(sign(digits) * floor(abs(digits) + 0.5))
    if (dig > 22) return(x)
    if (dig < 1L) dig <- 1L
    xtx_sql_round_macros(fr)
    sprintf("xtx_signif(%s, %d, %s, %s)", x, dig, xtx_sql_double(xtx_pow_di(10, -dig)),
            xtx_sql_double(log10(.Machine$double.xmax)))
}

## The SQL of the Math functions that Math.xtx_column (primitives.R)
## leaves to this file, or an error. `...`: the digits of round()/signif().
xtx_math_sql <- function(x, generic, ...) {
    a <- unclass(x)$sql
    fr <- unclass(x)$frame
    digits <- function(default) {
        args <- list(...)
        d <- if (length(args)) args[[1L]] else default
        if (xtx_is_column(d) || !is.numeric(d) || length(d) != 1L || is.na(d))
            stop(sprintf("%s() of a remote column: `digits` must be one number.", generic), call. = FALSE)
        as.double(d)
    }
    d <- sprintf("CAST(%s AS DOUBLE)", a)
    switch(generic,
        trunc = sprintf("trunc(%s)", d),
        round = xtx_sql_round(fr, a, digits(0)),
        signif = xtx_sql_signif(fr, a, digits(6)),
        ## (DuckDB refuses gamma(0) and lgamma(0); R gives NaN and Inf)
        gamma = sprintf("(CASE WHEN %1$s = 0 THEN CAST('NaN' AS DOUBLE) ELSE gamma(CASE WHEN %1$s = 0 THEN CAST(1 AS DOUBLE) ELSE %1$s END) END)", d),
        lgamma = sprintf("(CASE WHEN %1$s = 0 THEN CAST('Infinity' AS DOUBLE) ELSE lgamma(CASE WHEN %1$s = 0 THEN CAST(1 AS DOUBLE) ELSE %1$s END) END)", d),
        cumsum = , cumprod = , cummax = , cummin =
            stop(sprintf("%s() is not supported for remote columns: remote rows have no order.", generic), call. = FALSE),
        stop(sprintf("%s() is not supported for remote columns.", generic), call. = FALSE))
}

## ---- exported helpers ----

#' @rdname column_arithmetic
#' @export
col_ifelse <- function(test, yes, no, missing = NULL) {
    if (!is.null(missing)) {
        checkmate::assert_number(missing, na.ok = TRUE)
        return(xtx_ifelse_missing(test, yes, no, missing))
    }
    if (is.null(xtx_column_frame(test, yes, no))) {
        if (!is.logical(test)) {
            if (!is.numeric(test)) stop("`test` must be logical.", call. = FALSE)
            test <- as.logical(test)
        }
        if (!(is.numeric(yes) || is.logical(yes)) || !(is.numeric(no) || is.logical(no)))
            stop("`yes` and `no` must be numeric.", call. = FALSE)
        ## (a vector without attributes: xtx_c_ifelse() in memory)
        dim(test) <- NULL
        return(xtx_c_ifelse(test, as.double(yes), as.double(no)))
    }
    xtx_col_check(test, "test")
    xtx_col_check(yes, "yes")
    xtx_col_check(no, "no")
    xtx_c_ifelse(test, yes, no)
}

#' @rdname column_arithmetic
#' @export
col_pmin <- function(a, b, checked = TRUE) {
    checkmate::assert_flag(checked)
    xtx_col_check(a, "a")
    xtx_col_check(b, "b")
    if (checked) xtx_c_pmin(a, b) else xtx_fast_pmin(a, b)
}

#' @rdname column_arithmetic
#' @export
col_pmax <- function(a, b, checked = TRUE) {
    checkmate::assert_flag(checked)
    xtx_col_check(a, "a")
    xtx_col_check(b, "b")
    if (checked) xtx_c_pmax(a, b) else xtx_fast_pmax(a, b)
}

#' Functions of a column that base R cannot dispatch on
#'
#' `col_math()` applies to a column, in memory or remote, one of the
#' functions that GLM families need and that have no generic in base R:
#' the standard distribution functions, and forms of `log()` and
#' `1 - exp(-x)` for values in a known range. See [column_arithmetic] for
#' the arithmetic of columns.
#'
#' Names may change before the first release.
#'
#' @param x A column: a numeric vector in memory, a [design_column()]
#'   remotely.
#' @param fun `"pnorm"`, `"qnorm"`, `"dnorm"` and `"pcauchy"` (the standard
#'   distributions, as [stats::pnorm()] and the others compute them); and
#'   `"log_positive"`, the natural logarithm of values known to be positive
#'   (remotely without the guards of `log()` at zero and below).
#' @return A column (a double vector in memory, a `design_column` remotely).
#' @examples
#' col_math(c(-1, 0, 1), "pnorm")
#' @keywords internal
#' @export
col_math <- function(x, fun = c("pnorm", "qnorm", "dnorm", "pcauchy", "log_positive")) {
    fun <- match.arg(fun)
    xtx_col_check(x, "x")
    switch(fun,
           pnorm = xtx_c_pnorm(x), qnorm = xtx_c_qnorm(x), dnorm = xtx_c_dnorm(x), pcauchy = xtx_c_pcauchy(x),
           log_positive = xtx_fast_log(x))
}

## col_ifelse() with a value where the test is missing. In memory, the
## values of col_ifelse() with `missing` on the rows of a missing test;
## remotely one CASE: WHEN test THEN yes ELSE no when `missing` is `no`
## (the test named once, the shortest form: counts of rows), otherwise a
## third branch.
xtx_ifelse_missing <- function(test, yes, no, missing) {
    if (is.null(xtx_column_frame(test, yes, no))) {
        out <- col_ifelse(test, yes, no)
        na <- is.na(test)
        if (any(na)) out[na] <- missing
        return(out)
    }
    xtx_col_check(test, "test")
    xtx_col_check(yes, "yes")
    xtx_col_check(no, "no")
    f <- xtx_column_frame(test, yes, no)
    if (!xtx_is_column(test)) {
        ## (a value in R as the test of remote columns: `missing` where it
        ## is NA)
        if (length(test) != 1L) stop("`test` must be a column, or one value.", call. = FALSE)
        if (is.na(test)) return(xtx_column(f, xtx_column_operand(missing), notnull = !is.na(missing)))
        return(col_ifelse(test, yes, no))
    }
    y <- xtx_column_operand(yes)
    n <- xtx_column_operand(no)
    ## (a numeric test is TRUE where it is not zero; NaN is missing, where
    ## DuckDB would take it as true)
    if (!isTRUE(unclass(test)$lgl))
        test <- xtx_column_logical(xtx_column(f, sprintf("(NULLIF(CAST(%s AS DOUBLE), CAST('NaN' AS DOUBLE)) <> 0)",
                                                         unclass(test)$sql)))
    same <- !xtx_is_column(no) && length(no) == 1L && identical(as.double(no), as.double(missing))
    sql <- if (same) sprintf("(CASE WHEN %s THEN %s ELSE %s END)", xtx_column_operand(test), y, n)
           else sprintf("(CASE %s WHEN TRUE THEN %s WHEN FALSE THEN %s ELSE %s END)", xtx_column_operand(test), y, n,
                        xtx_column_operand(missing))
    xtx_column(f, sql, notnull = all(xtx_operand_notnull(yes)) && all(xtx_operand_notnull(no)) && !is.na(missing))
}

#' @rdname column_arithmetic
#' @export
col_sums <- function(..., na.rm = FALSE) {
    args <- list(...)
    if (!length(args)) return(numeric())
    if (!is.logical(na.rm) || anyNA(na.rm) || !length(na.rm) %in% c(1L, length(args)))
        stop("`na.rm` must be TRUE or FALSE, or one value per argument.", call. = FALSE)
    na.rm <- rep_len(na.rm, length(args))
    ex <- as.list(substitute(list(...)))[-1L]
    nm <- names(args) %||% rep("", length(args))
    bare <- nm == "" & vapply(ex, is.symbol, NA)
    nm[bare] <- vapply(ex[bare], as.character, "")
    for (i in seq_along(args)) xtx_col_check(args[[i]], if (nm[i] != "") nm[i] else sprintf("..%d", i))
    remote <- vapply(args, xtx_is_column, NA)
    out <- numeric(length(args))
    if (any(remote)) {
        ## (one query: the sums of every SQL expression of every remote
        ## argument; a column with k columns gives k sums)
        k <- vapply(args[remote], xtx_ncols, 1L)
        s <- xtx_column_sums(args[remote], rep(na.rm[remote], k))
        idx <- rep(which(remote), k)
        out[remote] <- vapply(which(remote), function(i) {
            v <- s[idx == i]
            if (anyNA(v) && !all(is.nan(v[is.na(v)]))) return(NA_real_)
            sum(v)
        }, 0)
    }
    out[!remote] <- vapply(which(!remote), function(i) as.double(sum(args[[i]], na.rm = na.rm[i])), 0)
    if (any(nm != "")) names(out) <- nm
    out
}

## The sums of every SQL expression of the remote columns in `cols` (a list;
## same design), from one query, `narm` one flag per expression: with narm,
## the sum of the values that are neither NULL nor NaN (0 for none);
## otherwise NA when a value is NULL, else the sum (NaN when a value is NaN,
## as in R). Each expression is named once in the SUM, and once more in the
## count of NULL values unless the column is known never NULL: the
## expressions of a fit (deviances over inverse links, which DuckDB expands
## inline) are long.
xtx_column_sums <- function(cols, narm) {
    x <- do.call(xtx_c_cbind, cols)
    fr <- xtx_frame_bind(unclass(x)$frame)
    sql <- unclass(x)$sql
    k <- length(sql)
    nn <- unclass(x)$notnull %||% rep(FALSE, k)
    nan <- "CAST('NaN' AS DOUBLE)"
    items <- ifelse(narm, sprintf("COALESCE(SUM(NULLIF(CAST(%s AS DOUBLE), %s)), 0) AS s%d", sql, nan, seq_len(k)),
                    sprintf("SUM(CAST(%s AS DOUBLE)) AS s%d", sql, seq_len(k)))
    cnt <- !narm & !nn
    items <- c(items, sprintf("SUM(CASE WHEN (%s) IS NULL THEN 1 ELSE 0 END) AS n%d", sql, seq_len(k))[cnt])
    res <- xtx_query(fr, sprintf("SELECT %s FROM (%s) AS s WHERE %s", paste(items, collapse = ", "),
                                fr$src$select, fr$filter))
    vapply(seq_len(k), function(j) {
        v <- as.double(res[[paste0("s", j)]])
        if (cnt[j] && as.numeric(res[[paste0("n", j)]]) > 0) return(NA_real_)
        if (is.na(v) && !is.nan(v)) 0 else v
    }, 0)
}

#' @rdname column_arithmetic
#' @export
col_const <- function(like, value) {
    checkmate::assert_number(value, na.ok = TRUE)
    if (!xtx_is_column(like) && !inherits(like, c("xtx_matrix", "xtx_remote")) && !is.numeric(like) && !is.logical(like))
        stop("`like` must be a column or a design.", call. = FALSE)
    xtx_c_const(like, value)
}

#' @rdname column_arithmetic
#' @export
cbind.xtx_column <- function(..., deparse.level = 1) {
    args <- list(...)
    for (i in seq_along(args)) xtx_col_check(args[[i]], sprintf("..%d", i))
    xtx_c_cbind(...)
}

methods::setOldClass("xtx_column")

#' @rdname column_arithmetic
#' @export
setMethod("rowSums", "xtx_column", function(x, na.rm = FALSE, dims = 1, ...) {
    xtx_defaults_only(na.rm, dims)
    xtx_c_rowsums(x)
})

## An argument of the helpers: a remote column, or a numeric or logical
## vector (in memory), or a number (with remote columns).
xtx_col_check <- function(x, what) {
    if (xtx_is_column(x)) return(invisible())
    if (inherits(x, c("xtx_matrix", "xtx_remote")))
        stop(sprintf("`%s` is a design, not a column: use design_column(), or a product with the design.", what),
             call. = FALSE)
    if (!is.numeric(x) && !is.logical(x))
        stop(sprintf("`%s` must be a column or a number.", what), call. = FALSE)
    invisible()
}

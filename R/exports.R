#' Compress a design and a response to sufficient statistics
#'
#' `design_patterns()` groups the rows of a design that have the same row of
#' the model matrix (and, with `keys`, the same values of the key columns)
#' into patterns, and summarizes a response `z` within each pattern: the
#' sum of the weights, the weighted mean, and the weighted sum of squares
#' around that mean. A weighted least squares fit of the pattern means on the
#' pattern design, with the summed weights as weights, has the same
#' coefficients as the fit on every row; its residual sum of squares plus
#' `within` is the residual sum of squares of the fit on every row. This is
#' the compression [design_lm()] uses when a design has few distinct rows
#' (factors and discrete numeric variables).
#'
#' In memory, the patterns are found from the factor codes and numeric
#' variables of the design, without forming the dense matrix. For a remote
#' design (see [design_remote]), the database groups the rows by the factors
#' and numeric variables of the design in one scan, and only the patterns are
#' returned to R. Both give the same patterns and the same statistics, up to
#' rounding and the order of the patterns.
#'
#' @param x A design: a `design_matrix` built from a data frame, or a remote
#'   design. Scaled designs (`s * D`) are not accepted.
#' @param z The response. In memory, a numeric vector with one value per row
#'   of `x` (such as a [design_column()] of `x`, or a product with it). For a
#'   remote design, the name of a numeric column of the remote data, a
#'   `DBI::SQL()` expression over its columns, or a [design_column()] of the
#'   design (such as `design_column(x, "y") - x %*% b`). Missing or infinite
#'   values are an error.
#' @param w Optional weights, in the same form as `z`: finite and
#'   non-negative. Rows with zero weight do not enter the statistics.
#' @param keys Optional extra grouping columns: patterns are then also split
#'   by their values. In memory, a numeric vector, matrix or data frame of
#'   numeric columns with one row per row of `x`. For a remote design, the
#'   names of numeric columns of the remote data, a `DBI::SQL()`
#'   expression, or a [design_column()] of the design.
#' @param compress `"auto"` (default), `TRUE` or `FALSE`. With `"auto"`, the
#'   patterns are returned only when compression pays (few patterns
#'   relative to rows): in memory, following the option
#'   `design.matrix.dedup` (see [design_algebra]) and at most half as many
#'   patterns as rows; remotely, following the option
#'   `design.matrix.remote_compress` (see [design_remote]). The two backends
#'   can therefore decide differently with `"auto"`. With `TRUE`, the
#'   patterns are returned unless every row is a pattern of its own. With
#'   `FALSE`, `NULL` is returned.
#'
#' @return `NULL` when compression does not apply: with `compress = FALSE`;
#'   when every row is its own pattern; with `"auto"`, when compression does
#'   not pay; with fewer than two rows; in memory, when the per-pattern sums
#'   overflow; remotely, for a design with `na.action = na.pass`.
#'   Otherwise a list with elements, one row or value per pattern (patterns
#'   whose weights sum to zero are dropped) in no particular order:
#'   \describe{
#'     \item{`x`}{The design of the patterns: an in-memory `design_matrix`
#'       with the columns of `x`, one row per pattern
#'       (`as.matrix(x)` gives the matrix).}
#'     \item{`z`}{The weighted mean of `z` in each pattern.}
#'     \item{`w`}{The sum of the weights in each pattern (the number of
#'       rows, without `w`). Always a numeric vector.}
#'     \item{`within`}{One number: the weighted sum of squares of `z` around
#'       the pattern means, summed over patterns.}
#'     \item{`rows`}{The number of rows of `x` summarized.}
#'     \item{`keys`}{A numeric matrix with one row per pattern and one column
#'       per key column (named after the key columns), or `NULL` without
#'       `keys`.}
#'   }
#'
#' @examples
#' d <- data.frame(y = rnorm(100), f = gl(4, 25), x = rep(1:5, 20))
#' D <- design_matrix(y ~ f + x, d)
#' p <- design_patterns(D, d$y, compress = TRUE)
#' nrow(p$x)
#' fit <- lm.wfit(as.matrix(p$x), p$z, p$w)
#' all.equal(unname(coef(fit)), unname(coef(lm(y ~ f + x, d))))
#' # residual sum of squares of the fit on every row
#' sum(p$w * fit$residuals^2) + p$within
#'
#' if (requireNamespace("duckdb", quietly = TRUE)) {
#'     con <- DBI::dbConnect(duckdb::duckdb())
#'     DBI::dbWriteTable(con, "d", d)
#'     R <- design_matrix(y ~ f + x, con, table = "d")
#'     pr <- design_patterns(R, "y", compress = TRUE)
#'     nrow(pr$x)
#'     DBI::dbDisconnect(con)
#' }
#' @export
design_patterns <- function(x, z, w = NULL, keys = NULL, compress = "auto") {
    if (!(isTRUE(compress) || isFALSE(compress) || identical(compress, "auto")))
        stop("`compress` must be \"auto\", TRUE or FALSE.", call. = FALSE)
    out <- if (inherits(x, "dm_remote")) {
        dm_patterns_remote(x, z, w, keys, compress)
    } else if (inherits(x, "design_matrix")) {
        dm_patterns_memory(x, z, w, keys, compress)
    } else {
        stop("`x` must be a design: a `design_matrix` or a remote design.", call. = FALSE)
    }
    if (is.null(out)) return(NULL)
    ## the stable public list
    np <- nrow(out$x)
    list(x = out$x, z = as.double(out$z), w = if (is.null(out$w)) rep(1, np) else as.double(out$w),
         within = as.double(out$within), rows = as.double(out$rows), keys = out$keys)
}

dm_patterns_memory <- function(x, z, w, keys, compress) {
    if (dm_is_comb(x)) stop("`x` must not be a scaled design.", call. = FALSE)
    n <- nrow(x)
    vec <- function(v, what) {
        ## (a one-column matrix, such as D %*% b, is a vector)
        if (is.matrix(v) && ncol(v) == 1L) v <- v[, 1L]
        if (!(is.numeric(v) || is.logical(v)) || !is.null(dim(v)) || length(v) != n)
            stop(sprintf("`%s` must be a numeric vector with one value per row of `x`.", what), call. = FALSE)
        as.double(v)
    }
    z <- vec(z, "z")
    if (!all(is.finite(z))) stop("`z` has missing or infinite values.", call. = FALSE)
    if (!is.null(w)) {
        w <- vec(w, "w")
        if (!all(is.finite(w)) || any(w < 0)) stop("`w` must be finite and non-negative.", call. = FALSE)
    }
    knames <- NULL
    if (!is.null(keys)) {
        if (is.data.frame(keys)) keys <- as.matrix(keys)
        if (!is.matrix(keys)) {
            keys <- matrix(keys, ncol = 1L)
            colnames(keys) <- "keys"
        }
        if (!(is.numeric(keys) || is.logical(keys)) || nrow(keys) != n)
            stop("`keys` must be numeric, with one row per row of `x`.", call. = FALSE)
        knames <- colnames(keys) %||% paste0("keys", seq_len(ncol(keys)))
        storage.mode(keys) <- "double"
    }
    out <- dmp_patterns(x, z, w, compress, keys = keys)
    if (!is.null(out$keys)) colnames(out$keys) <- knames
    out
}

## Remote: z, w and keys as column names, DBI::SQL() expressions or
## design_column()s, resolved as the other public remote inputs (dm_remote_cols(),
## which moves the design to the frame of a remote column), then the scan of
## dm_sql_patterns() on their SQL.
dm_patterns_remote <- function(x, z, w, keys, compress) {
    if (!is.null(x$scale)) stop("`x` must not be a scaled design.", call. = FALSE)
    zc <- dm_remote_cols(x, z, "z")
    D <- zc$D
    if (length(zc$sql) != 1L) stop("`z` must be one column.", call. = FALSE)
    ws <- NULL
    if (!is.null(w)) {
        wc <- dm_remote_cols(D, w, "w")
        D <- wc$D
        ws <- wc$sql
        if (length(ws) != 1L) stop("`w` must be one column.", call. = FALSE)
    }
    ks <- NULL
    knames <- NULL
    if (!is.null(keys)) {
        kc <- dm_remote_cols(D, keys, "keys", several = TRUE)
        D <- kc$D
        ks <- kc$sql
        knames <- kc$names %||% paste0("keys", seq_along(ks))
        if (length(knames) != length(ks)) knames <- paste0("keys", seq_along(ks))
    }
    ## the checks of the vectors in memory, in one query
    dbl <- function(e) sprintf("CAST(%s AS DOUBLE)", e)
    chk <- sprintf("SUM(CASE WHEN COALESCE(isfinite(%s), FALSE) THEN 0 ELSE 1 END) AS bad_z", dbl(zc$sql))
    if (!is.null(ws))
        chk <- c(chk, sprintf("SUM(CASE WHEN COALESCE(isfinite(%s) AND %s >= 0, FALSE) THEN 0 ELSE 1 END) AS bad_w",
                              dbl(ws), dbl(ws)))
    s <- dm_query(D, sprintf("SELECT %s FROM (%s) AS s WHERE %s", paste(chk, collapse = ", "), D$src$select, D$filter))
    if (isTRUE(as.numeric(s$bad_z) > 0)) stop("`z` has missing or infinite values.", call. = FALSE)
    if (!is.null(ws) && isTRUE(as.numeric(s$bad_w) > 0))
        stop("`w` must be finite and non-negative.", call. = FALSE)
    if (isFALSE(compress)) return(NULL)
    mode <- if (isTRUE(compress)) TRUE else dm_compress_mode()
    pt <- dm_sql_patterns(D, zc$sql, ws, keys = ks, mode = mode, fit = "lm")
    if (is.null(pt)) return(NULL)
    pos <- which(pt$sw > 0)
    K <- NULL
    if (!is.null(ks)) {
        K <- as.matrix(pt$keys[pos, , drop = FALSE])
        storage.mode(K) <- "double"
        colnames(K) <- knames
    }
    list(x = dm_dense_design(pt$x[pos, , drop = FALSE]), z = pt$sy[pos] / pt$sw[pos], w = pt$sw[pos],
         within = sum(pt$within), rows = sum(pt$n), keys = K)
}

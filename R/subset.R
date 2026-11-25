#' Select rows of a design matrix
#'
#' `x[i, ]` keeps the rows `i` of a design matrix, without forming the dense
#' matrix: the factor codes and numeric variables of the selected rows are
#' copied, and every product and summary works on the result as on any
#' other design. For a design built by applying a plan to `newdata`,
#' `design_matrix(plan, newdata)[i, ]` is identical to
#' `design_matrix(plan, newdata[i, ])`, row names included, in the strict
#' sense described in [design_matrix()]. Scaled rows and combinations of two
#' designs keep their scales for the selected rows.
#'
#' Only rows can be selected. `i` can be positive integers (repeats allowed),
#' negative integers or a logical vector, as for a matrix; missing values
#' are not allowed.
#'
#' @param x A `design_matrix`.
#' @param i Rows to keep.
#' @param j Not supported: columns cannot be selected.
#' @param ... Not used.
#' @param drop Ignored: the result is always a design matrix.
#' @return A `design_matrix` with `length(i)` rows (or the number selected).
#' @examples
#' d <- data.frame(x = rnorm(6), f = gl(3, 2))
#' D <- design_matrix(~ x * f, d)
#' as.matrix(D[c(2, 5), ])
#' @export
`[.design_matrix` <- function(x, i, j, ..., drop = FALSE) {
    if (!missing(j)) {
        stop("Only rows of a design matrix can be selected, as in `x[i, ]`.", call. = FALSE)
    }
    n_index <- nargs() - as.integer(!missing(drop))
    if (n_index < 3L) {
        stop("Select rows with `x[i, ]`.", call. = FALSE)
    }
    if (missing(i)) return(x)
    n <- nrow(x)
    rows <- dm_row_index(i, n)
    ## Row names of a design are unique, so selected names repeat exactly
    ## when the indices do: one counting pass, shared by both designs of a
    ## combination.
    dup <- length(rows) > 1L && any(tabulate(rows, nbins = n) > 1L)
    dm_rows(x, rows, dup)
}

## Positive row indices from a matrix-style subscript.
dm_row_index <- function(i, n) {
    if (is.logical(i)) {
        if (anyNA(i)) stop("Row indices cannot be missing.", call. = FALSE)
        if (length(i) != n) {
            stop(sprintf("A logical row index must have %d values (one per row).", n), call. = FALSE)
        }
        return(which(i))
    }
    if (!is.numeric(i) || is.object(i)) {
        stop("Rows must be selected with integer or logical indices.", call. = FALSE)
    }
    if (anyNA(i)) stop("Row indices cannot be missing.", call. = FALSE)
    i <- as.integer(i)
    if (length(i) && any(i < 0L)) {
        if (any(i > 0L)) stop("Cannot mix positive and negative row indices.", call. = FALSE)
        return(seq_len(n)[i])
    }
    if (length(i) && min(i) == 0L) i <- i[i != 0L]
    if (length(i) && max(i) > n) stop("Row index out of bounds.", call. = FALSE)
    i
}

## Rows `rows` (positive, within range) of a design or combination; `dup`:
## whether `rows` repeats an index.
dm_rows <- function(x, rows, dup) {
    if (dm_is_comb(x)) {
        pick <- function(s) if (length(s) > 1L) s[rows] else s
        return(dm_combine(dm_rows(x$design, rows, dup), pick(x$scale),
                          if (is.null(x$design2)) NULL else dm_rows(x$design2, rows, dup),
                          pick(x$scale2)))
    }
    out <- x
    out$cache <- new.env(parent = emptyenv())
    out$n <- length(rows)
    if (identical(x$path, "dense")) {
        X <- x$dense
        att <- attributes(X)
        Y <- X[rows, , drop = FALSE]
        for (nm in setdiff(names(att), c("dim", "dimnames"))) attr(Y, nm) <- att[[nm]]
        out["dense"] <- list(Y)
        return(out)
    }
    vars <- x$vars
    att <- attributes(vars)
    new <- lapply(vars, dm_var_rows, rows = rows)
    attributes(new) <- att
    attr(new, "row.names") <- dm_row_names(att[["row.names"]], rows, dup)
    out$vars <- new
    out
}

## Rows of one prepared variable, keeping its attributes (factor levels,
## class and contrasts; the parameters of poly() and the like).
dm_var_rows <- function(v, rows) {
    att <- attributes(v)
    if (is.matrix(v)) {
        y <- unclass(v)[rows, , drop = FALSE]
        dn <- dimnames(v)
        if (!is.null(dn[[1L]])) dn[[1L]] <- dn[[1L]][rows]
        att$dim <- dim(y)
        att$dimnames <- dn
    } else {
        y <- .subset(v, rows)
        if (!is.null(att$names)) att$names <- att$names[rows]
    }
    attributes(y) <- att
    y
}

## Row names of the selected rows, as `[.data.frame` gives them.
dm_row_names <- function(rn, rows, dup) {
    new <- if (is.integer(rn) && length(rn) == 2L && is.na(rn[1L])) rows else rn[rows]
    if (dup) new <- make.unique(as.character(new))
    new
}

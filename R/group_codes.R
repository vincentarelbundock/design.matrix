## Bare group codes, for code that sums over many groupings of the same rows
## (fits with fixed effects, in xtx.fit): the codes of design_groups()
## without its checks, and codes made from keys that can have missing
## values, which the caller drops with the rows.

#' Keys and bare group codes
#'
#' The groupings of [design_groups()] are checked and complete: every row
#' has a group. Code that drops rows and codes the groups again, such as a
#' fit with fixed effects, works on *keys* and *bare codes* instead:
#'
#' * `design_key()` returns the values of variables of the data over the
#'   rows of a design, of any type: a data frame in memory, a
#'   [design_column()] with one SQL column per variable remotely.
#'   `is.na()` of a remote key works on columns of any type, one result
#'   per variable.
#' * `key_codes()` codes the combinations of values of a key that has no
#'   missing value: the bare codes, the levels and the rows per level.
#' * `design_group_sums()` computes, for several groupings at once, the
#'   sums of the columns of a design (and of its weights and other columns)
#'   within groups, and the crosstabs of pairs of groupings.
#'
#' Bare codes (class `xtx_codes`) are a list with `codes` and `G`, the
#' number of groups. In memory `codes` is an integer vector in `1..G`.
#' Remotely it is a column of the source (the key itself, or one sort key
#' for a key of several variables) with `levels`, so that the database
#' groups on stored values and no join is needed; results come back in the
#' order of the levels. [group_sums()], [group_crosstab()], [group_lookup()]
#' and [group_nested()] accept bare codes, without checking them, as does
#' `design_group_sums()`, which also takes groupings of [design_groups()].
#'
#' Names may change before the first release.
#'
#' @param x A design: a `design_matrix` or a remote design.
#' @param vars The names of the variables of the key.
#' @param data In memory, the data frame of the rows of `x` (default: the
#'   model frame the design keeps).
#' @param key A key, from `design_key()`.
#' @param groups A list of groupings over the rows of `x`: bare codes, or
#'   groupings of [design_groups()].
#' @param w Weights, a column over the rows of `x`, or `NULL`.
#' @param extra `NULL`, or columns to sum within groups with the columns of
#'   the design.
#' @param pairs `NULL`, or a list of pairs of indices of `groups` whose
#'   crosstabs are wanted: the sums and the crosstabs then come from one
#'   fused pass in memory, with `w` the weights of every sum.
#' @return `design_key()`: a key. `key_codes()`: a list with `codes` (bare codes), `levels` (a data
#'   frame, one row per group in code order, with one column per variable of
#'   the key) and `n` (the rows per group).
#'   `design_group_sums()`: without `pairs`, a
#'   list with, per grouping, `S` (a `G` by `ncol(x)` matrix of sums of
#'   the columns of `x`, scaled by `w`), `n` (the sums of the weights) and
#'   `extra`; with `pairs`, a list with `gs` (that list) and `ct` (one
#'   sparse crosstab per pair).
#' @examples
#' d <- data.frame(y = rnorm(8), x = rnorm(8), id = c(1, 1, 2, 2, 3, 3, NA, 3))
#' D <- design_matrix(y ~ x, d)
#' k <- design_key(D, "id", d)
#' ok <- !is.na(d$id)
#' s <- design_subset(D, ok, list(k = k, y = d$y))
#' cd <- key_codes(s$cols$k)
#' cd$levels
#' group_sums(s$cols$y, cd$codes)
#' @keywords internal
#' @export
design_key <- function(x, vars, data = NULL) {
    checkmate::assert_multi_class(x, c("xtx_matrix", "xtx_remote"))
    checkmate::assert_character(vars, min.len = 1L, any.missing = FALSE)
    xtx_c_key(x, vars, data)
}

#' @rdname design_key
#' @export
key_codes <- function(key) xtx_c_dense_codes(key, raw = TRUE)

## The bare codes of a grouping of design_groups() (bare codes as they
## are): for the functions that take either.
xtx_group_codes <- function(x) {
    if (xtx_is_codes(x)) return(x)
    xtx_groups_check(x, "x", na_ok = TRUE)
    xtx_groups_codes(x)
}

#' @rdname design_key
#' @export
design_group_sums <- function(x, groups, w = NULL, extra = NULL, pairs = NULL) {
    checkmate::assert_multi_class(x, c("xtx_matrix", "xtx_remote"))
    if (!all(vapply(groups, xtx_is_codes, NA))) groups <- lapply(groups, xtx_group_codes)
    if (is.null(pairs)) return(xtx_p_colsums_groups(if (is.null(w)) x else xtx_p_scale(x, w), groups, extra))
    xtx_p_fe_sums(x, w, groups, extra = extra, pairs = pairs)
}

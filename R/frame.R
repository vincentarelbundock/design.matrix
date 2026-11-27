## design_frame(): stats::model.frame() without the work that changes nothing.
##
## On a large data frame, most of the time of model.frame() goes to two steps
## that often return what they were given:
##
##  * na.omit() and na.exclude() subset the whole frame,
##    `object[!omit, , drop = FALSE]`, even when no row is omitted: every
##    column is copied, and the row names are hashed to check that they
##    are unique;
##  * with `xlev`, every factor is recoded twice, `xi[, drop = TRUE]` then
##    `factor(xi, levels = xl, exclude = NULL)`, even when it already has the
##    requested levels.
##
## As for model.matrix.default() (see prep.R), we do not copy one version of
## model.frame.default(): we patch the installed one at load time. Two
## patches, each applied only when the code it touches is the code it was
## written for:
##
##  1. the na.action handed to C_modelframe goes through xtx_na(), which
##     replaces stats::na.omit and stats::na.exclude (and nothing else) by a
##     wrapper that returns the frame without copying it when that is what
##     the real function would produce (xtx_wrap_na());
##  2. in the `xlev` loop, the two recoding calls are skipped when
##     xtx_same_factor() shows that they are no-ops. The loop itself still
##     runs, with its warnings, errors and `data[[nm]] <-` assignment.
##
## Everything else (non-standard evaluation of `subset` and `...`, the order
## in which arguments are forced, how a missing `na.action` is resolved, the
## 'newdata' warning) is R's own code, and every other name in the patched
## function resolves in the stats namespace. A self-test at load compares
## the patched function with stats::model.frame() on small frames; a patch
## that fails it is dropped. With no patch left, design_frame() is
## stats::model.frame().

#' Build a model frame without redundant copies, or the remote frame of a model
#'
#' `design_frame()` takes the arguments of [stats::model.frame()] and returns
#' the same model frame, faster on large data. [design_matrix()] uses it to
#' build its model frames. When `data` is remote data (a DuckDB connection,
#' or Parquet paths), it returns the *remote frame* of the model instead:
#' the remote design of the rows the model uses, with its response, weights
#' and offset as columns over those rows.
#'
#' @details
#' **In memory.** The result is identical to that of `model.frame()` in the
#' strict sense described in [design_matrix()]: same values, same attributes
#' in the same order, same row names; warnings and errors are the same too.
#' `design_frame()` skips two steps of `model.frame()` when they would
#' return their input unchanged:
#'
#' * [stats::na.omit()] and [stats::na.exclude()] copy every column of the
#'   frame and check its row names for duplicates, even when no row has a
#'   missing value. `design_frame()` skips the copy when the frame has
#'   nothing to omit, and makes the check with a faster one.
#' * With `xlev`, `model.frame()` recodes every factor twice, even when it
#'   already has the requested levels. `design_frame()` leaves such factors
#'   as they are.
#'
#' Anything else runs the code of `model.frame()`: other `na.action`
#' functions, columns with a class or with attributes that subsetting
#' would drop or reorder, factors whose levels differ from `xlev`, and
#' objects with their own `model.frame()` method (the method is chosen by
#' [stats::model.frame()] itself).
#'
#' `design_frame()` is built when the package is loaded, from the
#' `model.frame.default()` of the running version of R. If that function,
#' the data frame methods of `na.omit()` and `na.exclude()` or
#' `[.data.frame` are not the ones this package was written for, or if a
#' self-test at load finds a difference, the corresponding shortcut is
#' turned off; without shortcuts, `design_frame()` is `model.frame()`.
#' Option `xtx.skip_na_copy = FALSE` turns both shortcuts off.
#'
#' **Remote data.** `data` is looked at first (it is evaluated before the
#' other arguments, after `formula`). When it is a DuckDB connection or a
#' character vector of Parquet paths, the arguments are, all by name except
#' `formula` and `data`:
#'
#' * `formula`: a model formula. The design is that of its right-hand side
#'   (see [design_remote]).
#' * `weights`, `offset`, `subset`: taken unevaluated, as `model.frame()`
#'   takes them: the name of a column, a `DBI::SQL()` expression, or an R
#'   expression over the columns. Expressions with a SQL translation run in
#'   the database; the others are evaluated in R, batch by batch, into a
#'   temporary table (with a message; `engine = "sql"` refuses). R
#'   expressions are evaluated in the environment of the formula.
#' * `na.action`: [stats::na.omit()] (the default) or [stats::na.exclude()]
#'   (the same rows), [stats::na.fail()] or [stats::na.pass()] (the errors
#'   of [stats::lm()] when a row has a missing value), or their names.
#' * `contrasts`: as the `contrasts.arg` of [stats::model.matrix()].
#' * `table`, `sql`, `engine`, `chunk_size`: as [design_matrix()] takes
#'   them for remote data (see [design_remote]).
#' * other arguments: bare names of columns that are not in the model but
#'   must not be missing on the rows used (the variables of fixed effects
#'   or clusters, say), as `model.frame()` keeps the extra variables it is
#'   given.
#'
#' The rows used are those of the subset where no variable of the model,
#' response, weight, offset or extra column is missing, and the levels of
#' factors are those of these rows, as
#' `model.frame(drop.unused.levels = TRUE)` gives them in memory. The
#' response is returned whatever its type, as `model.frame()` does; a
#' column given as `weights` or `offset` must be numeric (the errors of
#' `lm()`: `'weights' must be a numeric vector`, `'offset' must be numeric`).
#'
#' @param ... Arguments of [stats::model.frame()]: `formula`, `data`,
#'   `subset`, `na.action`, `drop.unused.levels`, `xlev`, and further
#'   variables such as `weights` or `offset`. With remote data, see
#'   *Remote data*.
#' @return In memory, a model frame, as returned by [stats::model.frame()].
#'
#'   With remote data, a list of class `xtx_remote_frame`:
#'   * `x`: the remote design of the rows used.
#'   * `response`: a list with `kind` and `ncol`, and by kind:
#'     `"numeric"` (a numeric column or an expression) and `"logical"` (a
#'     logical column, as 0/1): `y`, a [design_column()] (`ncol` is more
#'     than 1 when an expression evaluated in R gave a matrix; `y` is then
#'     its first column); `"factor"` (a character or factor column): `y`,
#'     its values as a text column, which compares with a string
#'     (`y == "a"`), and `levels`, the levels on the rows used;
#'     `"matrix"` (`cbind()` of `ncol` expressions): `cols`, a list of its
#'     first two columns (`cbind()` of one expression is that expression,
#'     but for a factor column, whose codes are then the one column of
#'     `cols`, as in memory); `"other"` (a column
#'     of another type, such as a date) and `"none"` (a one-sided formula):
#'     nothing.
#'   * `weights`: a column, or `NULL`.
#'   * `offset`: the sum of the `offset()` terms of the formula and of the
#'     `offset` argument, a column, or `NULL`.
#'   * `ncol`: the numbers of columns of the `weights` and of the `offset`,
#'     more than 1 when an expression evaluated in R gave a matrix (which a
#'     caller refuses, as [stats::lm()] does).
#'   * `types`: the kind of each column of the data (`"double"`,
#'     `"integer"`, `"logical"`, `"character"`, `"factor"` or `"other"`).
#'   * `terms`, `contrasts`, `xlevels`: the terms of the formula, and the
#'     contrasts and factor levels of the design.
#'   * `na.action`: `"omit"`, `"exclude"`, `"fail"` or `"pass"`.
#' @examples
#' d <- data.frame(y = rnorm(10), x = rnorm(10), f = gl(2, 5))
#' identical(design_frame(y ~ x * f, d), model.frame(y ~ x * f, d))
#'
#' if (requireNamespace("duckdb", quietly = TRUE)) {
#'     con <- DBI::dbConnect(duckdb::duckdb())
#'     d$w <- runif(10)
#'     d$id <- rep(1:5, 2)
#'     DBI::dbWriteTable(con, "d", d)
#'     fr <- design_frame(y ~ x * f, con, table = "d", weights = 2 * w, subset = x > -1, id)
#'     fr
#'     dim(fr$x)
#'     sum(fr$weights)
#'     b <- solve(crossprod(fr$x), crossprod(fr$x, fr$response$y))
#'     DBI::dbDisconnect(con)
#' }
#' @export
design_frame <- function(...) {
    ## stats::model.frame() does the dispatch, called from an environment
    ## on top of the caller's, so that it finds every method where it
    ## would if the caller had called it. The arguments travel in `...`
    ## only: they stay the caller's promises, nothing is forced here, and
    ## substitute() and missing() give in the method what they give under
    ## model.frame().
    env <- caller <- parent.frame()
    ## Remote data: `data` (by name, or the argument after the formula) is
    ## the only argument looked at, after the formula: the two that
    ## model.frame() evaluates first, in its order (remote_resolve.R).
    if (...length() > 1L) {
        i <- xtx_frame_data_arg(...names(), ...length())
        if (i[2L] > 0L) {
            ## (an object with a model.frame() method of its own may never
            ## evaluate its second argument: it is not looked at)
            f <- if (i[1L] > 0L) ...elt(i[1L])
            if ((!is.object(f) || inherits(f, "formula")) && xtx_is_remote_source(...elt(i[2L])))
                return(xtx_frame_remote(..., .caller = caller))
        }
    }
    frame <- the$frame
    if (!is.null(frame) && isTRUE(getOption("xtx.skip_na_copy", TRUE)) &&
        xtx_frame_default_is_stats(frame$original, caller)) {
        ## S3 dispatch looks for methods in the calling environment before
        ## the registered ones: with the patched function there under the
        ## name `model.frame.default`, it is found where stats' default
        ## method would have been.
        env <- new.env(parent = caller)
        env$model.frame.default <- frame$fun
    }
    generic <- xtx_frame_generic
    environment(generic) <- env
    generic(...)
}

## Among `n` arguments named `nm` (NULL: none is named), the positions of
## `formula` and `data` as model.frame() matches them: by name, else the
## first and the next unnamed arguments. 0: not given.
xtx_frame_data_arg <- function(nm, n) {
    if (is.null(nm)) return(c(1L, 2L))
    un <- which(is.na(nm) | !nzchar(nm))
    ## (names match partially, as in a call of model.frame(): `f =` is the
    ## formula of the generic, `dat =` the data of model.frame.default(),
    ## where `d =` could also be `drop.unused.levels`)
    pre <- function(arg, min = 1L) {
        k <- which(!is.na(nm) & nchar(nm) >= min & startsWith(arg, nm))
        if (length(k)) k[1L] else 0L
    }
    f <- pre("formula")
    d <- pre("data", 2L)
    k <- 0L
    if (f == 0L) f <- if (length(un) > k) un[k <- k + 1L] else 0L
    if (d == 0L) d <- if (length(un) > k) un[k + 1L] else 0L
    c(f, d)
}

xtx_frame_generic <- function(...) stats::model.frame(...)

## Would stats::model.frame(), called from `caller`, use stats' own default
## method when it reaches the default? Not if a method was registered over
## it, or if another `model.frame.default` is visible from the caller.
xtx_frame_default_is_stats <- function(original, caller) {
    tbl <- asNamespace("stats")[[".__S3MethodsTable__."]]
    if (!identical(tbl[["model.frame.default"]], original)) return(FALSE)
    f <- get0("model.frame.default", envir = caller, mode = "function")
    is.null(f) || identical(f, original)
}

## Patch 1: the na.action -------------------------------------------------

## Evaluated in the patched model.frame.default(), in place of `na.action`
## in the call of C_modelframe. `rho` is the frame of that call: the C code
## turns a string into a symbol and calls it there, so a string names the
## first function of that name seen from `rho`.
xtx_na <- function(na.action, rho) {
    fn <- na.action
    if (is.character(fn)) {
        if (!length(fn) || is.na(fn[1L])) return(na.action)
        fn <- get0(fn[1L], envir = rho, mode = "function")
    }
    if (!is.function(fn)) return(na.action)
    if (identical(fn, stats::na.omit)) return(the$na_omit)
    if (identical(fn, stats::na.exclude)) return(the$na_exclude)
    na.action
}

## Is column `x` of a frame of `n` rows left as it is by na.omit() when no
## row is omitted? na.omit() subsets each column with an all-TRUE logical
## index, `x[i]` or `x[i, , drop = FALSE]`, and C_modelframe then copies
## back most of the attributes that this dropped. The result is `x` itself,
## attributes in the same order, for the columns listed here; any other
## column (a class with its own `[` or is.na() method, names, a "tsp" or
## any extra attribute, a factor with contrasts, an array that is not a
## matrix, a matrix without columns, a column longer than the frame as
## `subset` makes of a 3-d array) goes to the real function. A column with
## a missing value does too.
xtx_col_unchanged <- function(x, n) {
    if (!is.atomic(x) || is.null(x)) return(FALSE)
    an <- names(attributes(x))
    if (is.null(an)) return(length(x) == n && !anyNA(x))
    cl <- oldClass(x)
    d <- attr(x, "dim")
    ok <- if (is.null(d)) {
        length(x) == n &&
            if (identical(cl, "factor") || identical(cl, c("ordered", "factor"))) {
                typeof(x) == "integer" && identical(an, c("levels", "class"))
            } else {
                identical(cl, "AsIs") && identical(an, "class")
            }
    } else {
        dn <- attr(x, "dimnames")
        length(d) == 2L && d[1L] == n && d[2L] >= 1L &&
            (is.null(cl) || identical(cl, "AsIs")) &&
            identical(an, c("dim", if (!is.null(dn)) "dimnames", if (!is.null(cl)) "class")) &&
            (is.null(dn) || all(names(attributes(dn)) %in% "names"))
    }
    ok && !anyNA(x)
}

## The replacement of stats::na.omit or stats::na.exclude (`generic`), whose
## method for data frames is `method`. C_modelframe calls it with the frame:
## a list of class "data.frame" with attributes names, terms, row.names.
##
## When no column would change, what is left of
## `object[!omit, , drop = FALSE]` is what `[.data.frame` does to the data
## frame itself: its attributes are kept in order, with row.names then
## class moved last, and row names that are NA or duplicated are rewritten
## by make.unique(). (They can be: when `data` is not a data frame they are
## the names of the response.)
##  * Row names known to be unique stay as they are, and the attributes
##    are reordered here: automatic row names, and row names that a pass
##    over them shows to be without NA or duplicates (xtx_int_unique(),
##    xtx_str_unique() in C, much faster than anyDuplicated()).
##  * Otherwise (duplicates, NA, or row names for which the pass cannot
##    tell) the real method runs on a copy of the frame without columns,
##    which gives the row names and the attributes as R makes them; the
##    columns are put back. This costs the anyDuplicated() of the row
##    names, but none of the column copies.
xtx_wrap_na <- function(method, generic) {
    force(method)
    force(generic)
    stats_ns <- asNamespace("stats")
    function(object, ...) {
        ## (C_modelframe never does this: dispatch as it would)
        if (...length() || !identical(oldClass(object), "data.frame"))
            return(do.call(generic, list(object, ...), envir = stats_ns))
        n <- .row_names_info(object, 2L)
        if (n < 1L || !length(object)) return(method(object))
        cols <- unclass(object)
        for (x in cols) if (!xtx_col_unchanged(x, n)) return(method(object))
        a <- attributes(object) # row names in full, as `[.data.frame` reads them
        rn <- .row_names_info(object, 0L)
        if (is.null(attributes(rn)) &&
            if (is.integer(rn)) {
                (length(rn) == 2L && is.na(rn[1L])) || isTRUE(.Call(C_xtx_int_unique, rn))
            } else {
                is.character(rn) && isTRUE(.Call(C_xtx_str_unique, rn))
            }) {
            keep <- !(names(a) %in% c("row.names", "class"))
            a <- c(a[keep], a["row.names"], a["class"])
        } else {
            shell <- list()
            a[["names"]] <- character()
            attributes(shell) <- a
            shell <- method(shell)
            a <- attributes(shell)
            if (length(shell) || !identical(a[["names"]], character())) return(method(object))
            a[["names"]] <- names(cols)
        }
        attributes(cols) <- a
        cols
    }
}

## Patch 2: the xlev loop -------------------------------------------------

## Is `factor(xi[, drop = TRUE], levels = xl, exclude = NULL)` the factor
## `xi` itself? Yes for a plain factor whose levels are `xl`: integer codes
## that are NA or in range, no attribute but levels and class (names and
## contrasts would be dropped or moved), class exactly "factor" or
## c("ordered", "factor") (factor() returns nothing else), levels and `xl`
## plain character vectors, identical and with the same encoding marks
## (identical() takes a string marked latin1 and the same string marked
## UTF-8 for one; the result would have the marks of `xl`), without NA
## (factor() would give the NA level to missing values) or duplicates (an
## error in R).
xtx_same_factor <- function(xi, xl) {
    if (typeof(xi) != "integer" || !identical(names(attributes(xi)), c("levels", "class")))
        return(FALSE)
    cl <- oldClass(xi)
    if (!identical(cl, "factor") && !identical(cl, c("ordered", "factor"))) return(FALSE)
    lv <- attr(xi, "levels")
    if (!is.character(xl) || !is.null(attributes(xl)) ||
        !is.character(lv) || !is.null(attributes(lv))) return(FALSE)
    if (!identical(lv, xl) || !identical(Encoding(lv), Encoding(xl)) ||
        anyNA(xl) || anyDuplicated(xl)) return(FALSE)
    .Call(C_xtx_codes_ok, xi, length(xl))
}

## For drop.unused.levels: length(unique(x[!is.na(x)])) < length(levels(x))
## for a factor x, from its codes in C when they are all NA or in range.
xtx_unused_levels <- function(x) {
    r <- if (typeof(x) == "integer") .Call(C_xtx_levels_unused, x, length(levels(x))) else NA
    if (is.na(r)) length(unique(x[!is.na(x)])) < length(levels(x)) else r
}

## Building the patched function -------------------------------------------

## The statements this code was written for. A version of R where one of
## them reads differently does not get the corresponding patch.

## model.frame.default(): the call of the C routine.
xtx_frame_c_call <- quote(.External2(C_modelframe, formula, rownames, variables, varnames,
                                    extras, extranames, subset, na.action))

## model.frame.default(): the whole statement that applies `xlev`.
xtx_frame_xlev_stmt <- quote(
    if (length(xlev)) {
        for (nm in names(xlev)) if (!is.null(xl <- xlev[[nm]])) {
            xi <- data[[nm]]
            if (is.character(xi)) xi <- as.factor(xi)
            if (!is.factor(xi) || is.null(nxl <- levels(xi)))
                warning(gettextf("variable '%s' is not a factor", nm), domain = NA)
            else {
                ctr <- attr(xi, "contrasts")
                xi <- xi[, drop = TRUE]
                nxl <- levels(xi)
                if (any(m <- is.na(match(nxl, xl))))
                    stop(sprintf(ngettext(length(m), "factor %s has new level %s",
                                          "factor %s has new levels %s"),
                                 nm, paste(nxl[m], collapse = ", ")), domain = NA)
                data[[nm]] <- factor(xi, levels = xl, exclude = NULL)
                if (!identical(attr(data[[nm]], "contrasts"), ctr))
                    warning(gettextf("contrasts dropped from factor %s", nm),
                            call. = FALSE, domain = NA)
            }
        }
    } else if (drop.unused.levels) {
        for (nm in names(data)) {
            x <- data[[nm]]
            if (is.factor(x) && length(unique(x[!is.na(x)])) < length(levels(x))) {
                ctr <- attr(x, "contrasts")
                data[[nm]] <- x[, drop = TRUE]
                if (!identical(attr(data[[nm]], "contrasts"), ctr))
                    warning(gettextf("contrasts dropped from factor %s due to missing levels", nm),
                            domain = NA, call. = FALSE)
            }
        }
    })

## na.omit.data.frame(); na.exclude.data.frame() is the same with the class
## "exclude".
xtx_frame_na_body <- function(class) {
    bquote({
        n <- length(object)
        omit <- logical(nrow(object))
        vars <- seq_len(n)
        for (j in vars) {
            x <- object[[j]]
            if (!is.atomic(x)) next
            x <- is.na(x)
            d <- dim(x)
            if (is.null(d) || length(d) != 2L) omit <- omit | x
            else for (ii in 1L:d[2L]) omit <- omit | x[, ii]
        }
        xx <- object[!omit, , drop = FALSE]
        if (any(omit > 0L)) {
            temp <- setNames(seq(omit)[omit], attr(object, "row.names")[omit])
            attr(temp, "class") <- .(class)
            attr(xx, "na.action") <- temp
        }
        xx
    })
}

## `[.data.frame`: the statements that make `x[i, , drop = FALSE]` what
## xtx_wrap_na() takes it to be (attributes copied, row.names and class set
## last, columns subset one by one, row names subset).
xtx_frame_subset_stmts <- list(
    quote(x <- .Internal(copyDFattr(xx, x))),
    quote(oldClass(x) <- attr(x, "row.names") <- NULL),
    quote(x[[j]] <- if (length(dim(xj)) != 2L) xj[i] else xj[i, , drop = FALSE]),
    quote(rows <- rows[i]),
    quote(attr(x, "row.names") <- rows),
    quote(oldClass(x) <- oldClass(xx)))

## Replace every occurrence of `from` in `e` by `to`, counting them in
## `count$n`.
xtx_replace <- function(e, from, to, count) {
    if (identical(e, from)) {
        count$n <- count$n + 1L
        return(to)
    }
    if (is.call(e)) {
        for (i in seq_along(e)) {
            ei <- e[[i]]
            if (!missing(ei) && !is.null(ei)) e[[i]] <- xtx_replace(ei, from, to, count)
        }
    }
    e
}

## Number of occurrences of `what` in `e`.
xtx_count <- function(e, what) {
    count <- new.env(parent = emptyenv())
    count$n <- 0L
    xtx_replace(e, what, what, count)
    count$n
}

## Returns list(fun, original, na, xlev): the patched copy of `f`, `f`, and
## which patches it has; or NULL when none applies. `na` and `xlev` say
## which patches to try. The other arguments are the functions the na
## patch depends on (arguments so that tests can alter them).
xtx_make_frame <- function(f = get("model.frame.default", envir = asNamespace("stats")),
                          na = TRUE, xlev = TRUE,
                          omit = get("na.omit.data.frame", envir = asNamespace("stats")),
                          exclude = get("na.exclude.data.frame", envir = asNamespace("stats")),
                          subset = get("[.data.frame", envir = baseenv())) {
    if (!is.function(f)) return(NULL)
    b <- body(f)
    if (!is.call(b) || !identical(b[[1L]], as.name("{"))) return(NULL)
    stmts <- as.list(b)
    has_na <- has_xlev <- FALSE

    if (na &&
        identical(body(omit), xtx_frame_na_body("omit")) &&
        identical(body(exclude), xtx_frame_na_body("exclude")) &&
        all(vapply(xtx_frame_subset_stmts, function(s) xtx_count(body(subset), s) == 1L, NA))) {
        count <- new.env(parent = emptyenv())
        count$n <- 0L
        to <- xtx_frame_c_call
        to[[length(to)]] <- quote(.xtx_na(na.action, environment()))
        new <- lapply(stmts, xtx_replace, xtx_frame_c_call, to, count)
        if (count$n == 1L) {
            stmts <- new
            has_na <- TRUE
        }
    }

    idx <- which(vapply(stmts, function(e) identical(e, xtx_frame_xlev_stmt), NA))
    if (xlev && length(idx) == 1L) {
        c1 <- new.env(parent = emptyenv())
        c2 <- new.env(parent = emptyenv())
        c3 <- new.env(parent = emptyenv())
        c1$n <- c2$n <- c3$n <- 0L
        ## `.xtx_same` is a variable of the patched function's frame; R's
        ## code runs unchanged whenever it is FALSE
        s <- xtx_replace(stmts[[idx]], quote(xi <- xi[, drop = TRUE]),
                        quote(xi <- if (.xtx_same <- .xtx_same_factor(xi, xl)) xi else xi[, drop = TRUE]),
                        c1)
        s <- xtx_replace(s, quote(data[[nm]] <- factor(xi, levels = xl, exclude = NULL)),
                        quote(data[[nm]] <- if (.xtx_same) xi else factor(xi, levels = xl, exclude = NULL)),
                        c2)
        ## (drop.unused.levels: whether a level is unused, from one pass
        ## over the codes in C instead of unique() of a copy)
        s <- xtx_replace(s, quote(length(unique(x[!is.na(x)])) < length(levels(x))),
                        quote(.xtx_unused_levels(x)), c3)
        if (c1$n == 1L && c2$n == 1L && c3$n == 1L) {
            stmts[[idx]] <- s
            has_xlev <- TRUE
        }
    }

    if (!has_na && !has_xlev) return(NULL)
    fun <- f
    body(fun) <- as.call(stmts)
    ## Every name resolves as in stats, except the two helpers.
    env <- new.env(parent = environment(f))
    env$.xtx_na <- xtx_na
    env$.xtx_same_factor <- xtx_same_factor
    env$.xtx_unused_levels <- xtx_unused_levels
    environment(fun) <- env
    list(fun = compiler::cmpfun(fun), original = f, na = has_na, xlev = has_xlev)
}

## Does the patched function return what stats::model.frame() returns? A
## check on small frames that go through each shortcut: it covers what the
## comparison of function bodies cannot see (the C code of model.frame()
## and of subsetting).
xtx_frame_selftest <- function(fun) {
    same <- function(...) {
        identical(stats::model.frame(...), fun(...),
                  num.eq = FALSE, single.NA = FALSE, attrib.as.set = FALSE)
    }
    probes <- function() {
        d <- data.frame(y = c(1.5, 2, 3, 4, 5), x = c(1L, 2L, 4L, 3L, 5L),
                        l = c(TRUE, FALSE, TRUE, TRUE, FALSE), s = c("u", "v", "u", "w", "v"),
                        f = factor(c("a", "b", "a", "c", "b")),
                        o = factor(c("lo", "hi", "lo", "hi", "lo"), levels = c("lo", "hi"), ordered = TRUE),
                        stringsAsFactors = FALSE)
        d$m <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), 5L, 2L, dimnames = list(NULL, c("p", "q")))
        d$a <- I(c(5, 4, 3, 2, 1))
        named <- d
        attr(named, "row.names") <- c("r1", "r2", "r3", "r4", "r5")
        holes <- d
        holes$x[2L] <- NA
        holes$f[4L] <- NA
        ## not a data frame: the row names are the names of the response
        lst <- list(y = c(a = 1, a = 2, b = 3, 4, 5), x = c(1, 2, 4, 3, 5))
        env <- baseenv()
        fo <- stats::as.formula("y ~ x + l + s + f + o + m + a", env = env)
        fx <- stats::as.formula("~ x + f + o + s", env = env)
        fl <- stats::as.formula("y ~ x", env = env)
        xl <- list(f = c("a", "b", "c"), o = c("lo", "hi"), s = c("u", "v", "w"))
        xr <- list(f = c("c", "b", "a"), o = c("lo", "hi"))
        ## automatic, increasing, unsorted, short and character row names;
        ## missing values
        frames <- list(d, d[c(1L, 3L, 4L), ], d[c(4L, 1L, 3L), ], d[1:2, ], named, holes)
        for (na in list(stats::na.omit, stats::na.exclude)) {
            for (dd in frames) if (!same(fo, dd, na.action = na)) return(FALSE)
            if (!same(fx, d, na.action = na, xlev = xl)) return(FALSE)
            if (!same(fx, holes, na.action = na, xlev = xl)) return(FALSE)
            if (!same(fl, lst, na.action = na)) return(FALSE)
            if (!same(fo, d, na.action = na, subset = c(TRUE, FALSE, TRUE, TRUE, TRUE),
                      weights = c(1, 2, 1, 2, 1)))
                return(FALSE)
        }
        if (!same(fo, d, na.action = "na.omit")) return(FALSE)
        for (na in list(stats::na.pass, NULL)) {
            if (!same(fx, d, na.action = na, xlev = xl)) return(FALSE)
            if (!same(fx, d, na.action = na, xlev = xr)) return(FALSE)
            if (!same(fx, holes, na.action = na, xlev = xl)) return(FALSE)
        }
        TRUE
    }
    ## (no probe warns under stats::model.frame(): a warning is a failure)
    isTRUE(tryCatch(probes(), warning = function(w) FALSE, error = function(e) FALSE))
}

## The patched function with every patch that applies and passes the
## self-test, or NULL.
xtx_init_frame <- function() {
    omit <- get("na.omit.data.frame", envir = asNamespace("stats"))
    exclude <- get("na.exclude.data.frame", envir = asNamespace("stats"))
    the$na_omit <- xtx_wrap_na(omit, stats::na.omit)
    the$na_exclude <- xtx_wrap_na(exclude, stats::na.exclude)
    frame <- xtx_make_frame()
    if (is.null(frame) || xtx_frame_selftest(frame$fun)) return(frame)
    ## One patch does not reproduce stats on this version of R: keep the
    ## other, if it does alone.
    for (which in list(c(TRUE, FALSE), c(FALSE, TRUE))) {
        if (!frame[[c("na", "xlev")[which]]]) next
        one <- xtx_make_frame(na = which[1L], xlev = which[2L])
        if (!is.null(one) && xtx_frame_selftest(one$fun)) return(one)
    }
    NULL
}

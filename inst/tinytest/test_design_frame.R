source("helpers.R")
library(xtx)

# design_frame() must return what stats::model.frame() returns: strictly
# identical frames, same warnings, same error messages. Most cases below
# are frames on which a shortcut that looks harmless gives something else.

same_frame <- function(f, d, ...) {
    args <- list(f, d, ...)
    same_outcome(capture(function() do.call(stats::model.frame, args)),
                 capture(function() do.call(design_frame, args)))
}

fr <- xtx:::the$frame
# On the versions of R this was written for, both shortcuts are on. (On
# another version they may be off: design_frame() is then model.frame().)
if (getRversion() >= "4.5.0" && getRversion() < "4.6.0") {
    expect_true(isTRUE(fr$na))
    expect_true(isTRUE(fr$xlev))
}

old_opt <- options(na.action = "na.omit")

x <- c(1, 2, 3, 4, 6)
y <- c(1.5, 2, 3, 4, 5)
base_d <- function() {
    data.frame(y = y, x = x, f = factor(c("a", "b", "a", "c", "b")),
               g = factor(c("u", "v", "u", "v", "u")))
}
cases <- list()
add <- function(label, d, f = y ~ x, ...) {
    cases[[length(cases) + 1L]] <<- list(label = label, d = d, f = f, args = list(...))
}

# Columns: attributes that subsetting drops, moves or keeps ------------------
d <- base_d(); attr(d$x, "label") <- "my label"; add("label attribute", d)
d <- base_d(); comment(d$x) <- "c"; add("comment", d)
d <- base_d(); contrasts(d$f) <- contr.sum(3); add("factor with contrasts", d, y ~ f)
d <- base_d(); d$f <- structure(as.integer(d$f), class = "factor", levels = c("a", "b", "c"))
add("factor with attributes (class, levels)", d, y ~ f)
d <- base_d(); attr(d$f, "label") <- "L"; add("factor with a label", d, y ~ f)
d <- base_d(); d$f <- structure(d$f, class = c("factor", "ordered")); add("class c(factor, ordered)", d, y ~ f)
d <- base_d(); d$f <- structure(d$f, class = "ordered"); add("class ordered alone", d, y ~ f)
d <- base_d(); d$f <- I(d$f); add("I(factor)", d, y ~ f)
d <- base_d(); d$x <- I(structure(d$x, foo = 1)); add("AsIs with attributes (foo, class)", d)
d <- base_d(); d$x <- I(d$x); add("AsIs", d)
d <- base_d(); d$s <- I(letters[1:5]); add("AsIs character", d, y ~ s)
d <- base_d(); add("scale()", d, y ~ scale(x))
d <- base_d(); d$m <- matrix(1:10 + 0, 5, 2); add("matrix", d, y ~ m)
d <- base_d(); d$m <- matrix(1:10 + 0, 5, 2, dimnames = list(letters[1:5], c("A", "B")))
add("matrix with dimnames", d, y ~ m)
d <- base_d(); d$m <- matrix(1:10 + 0, 5, 2, dimnames = list(r = letters[1:5], c = c("A", "B")))
add("matrix with named dimnames", d, y ~ m)
d <- base_d(); m <- structure(1:10 + 0, foo = "bar"); dim(m) <- c(5, 2); d$m <- m
add("matrix with attributes (foo, dim)", d, y ~ m)
d <- base_d(); m <- matrix(1:10 + 0, 5, 2); attr(m, "foo") <- 1; dimnames(m) <- list(NULL, c("A", "B")); d$m <- m
add("matrix with attributes (dim, foo, dimnames)", d, y ~ m)
d <- base_d(); d$m <- I(matrix(1:10 + 0, 5, 2, dimnames = list(NULL, c("A", "B"))))
add("AsIs matrix", d, y ~ m)
add("matrix without columns", list(y = y, m = matrix(numeric(0), 5, 0)), y ~ m)
add("3-d array", list(y = y, a = array(1:20 + 0, c(5, 2, 2))), y ~ a)
add("3-d array and subset", list(y = y, a = array(1:20 + 0, c(5, 2, 2))), y ~ a,
    subset = c(TRUE, FALSE, TRUE, TRUE, TRUE))
add("1-d array", list(y = y, a = array(1:5 + 0, 5)), y ~ a)
add("1-d array with an attribute", list(y = y, a = structure(array(1:5 + 0, 5), foo = 1)), y ~ a)
d <- base_d(); d$r <- as.raw(1:5); add("raw", d, y ~ r)
d <- base_d(); d$z <- complex(real = 1:5, imaginary = c(1, 2, NaN, 4, 5)); add("complex with NaN", d, y ~ z)
d <- base_d(); d$x[2] <- NaN; add("NaN", d)
d <- base_d(); d$ch <- c("a", NA, "c", "d", "e"); add("character NA", d, y ~ ch)
d <- base_d(); d$dt <- as.Date("2020-01-01") + 0:4; add("Date", d, y ~ dt)
d <- base_d(); d$t <- ts(1:5); add("ts", d, y ~ t)
add("unclassed ts (tsp attribute)", list(y = y, u = unclass(ts(x))), y ~ u)
add("list column", list(y = y, l = list(1, 2, 3, 4, 5)), y ~ l)
d <- base_d()
add("poly()", d, y ~ poly(x, 2)); add("matrix response", d, cbind(y, x) ~ f)
add("~ 1", d, ~ 1); add("y ~ 1", d, y ~ 1); add("~ 0", d, ~ 0)
add("no rows", d[0, ], y ~ x + f); add("one row", d[1, ], y ~ x + f); add("two rows", d[1:2, ], y ~ x + f)
add("three rows", d[1:3, ], y ~ x + f); add("rows 4:5", d[4:5, ], y ~ x + f)
add("x (label, names)", list(y = y, x = structure(x, label = "L", names = letters[1:5])))
add("x (names, label)", list(y = y, x = structure(x, names = letters[1:5], label = "L")))
add("factor (levels, class, names)", list(y = y, x = structure(factor(letters[1:5]), names = LETTERS[1:5])))
add("AsIs (class, names)", list(y = y, x = structure(1:5, class = "AsIs", names = letters[1:5])))
xx <- x; attr(xx, "names") <- structure(letters[1:5], foo = "bar")
add("names carrying an attribute", list(y = y, x = xx))
m <- matrix(1:10 + 0, 5, 2); attr(m, "dimnames") <- structure(list(NULL, c("A", "B")), foo = 1)
add("dimnames carrying an attribute", list(y = y, m = m), y ~ m)

# Row names -------------------------------------------------------------------
d <- base_d(); attr(d, "row.names") <- c("a", "a", "b", "c", "d"); add("duplicated character row names", d)
d <- base_d(); attr(d, "row.names") <- c(1L, 1L, 2L, 3L, 4L); add("duplicated integer row names", d)
d <- base_d(); attr(d, "row.names") <- c(NA, "b", "c", "d", "e"); add("NA row name", d)
d <- base_d(); attr(d, "row.names") <- c(NA_integer_, 5L); add("row names c(NA, +n)", d)
d <- base_d(); attr(d, "row.names") <- c(5L, 4L, 3L, 2L, 1L); add("decreasing row names", d)
d <- base_d(); attr(d, "row.names") <- c(2L, 3L, 5L, 8L, 9L); add("increasing row names", d)
d <- base_d(); attr(d, "row.names") <- c(a = 1L, b = 2L, c = 3L, d = 4L, e = 7L); add("named row names", d)
d <- base_d(); rownames(d) <- c("NA", "b", "c", "d", "e"); add("row name 'NA'", d)
d <- base_d(); rownames(d) <- paste0("r", 5:1); add("character row names", d)
d <- base_d(); attr(d, "extra") <- "kept"; add("data frame with another attribute", d)
yn <- c(a = 1, a = 2, b = 3, c = 4, b = 5)
add("duplicated names on the response, list", list(y = yn, x = x))
add("duplicated names on the response, environment", list2env(list(y = yn, x = x)))
y3 <- c(1, 2, 3, 4, 5); names(y3)[2] <- "b"
add("partly named response", list(y = y3, x = x))
ym <- matrix(1:10 + 0, 5, 2, dimnames = list(c("r", "r", "s", "t", "u"), NULL))
add("duplicated row names on a matrix response", list(ym = ym, x = x), ym ~ x)

# na.action -------------------------------------------------------------------
d <- base_d(); contrasts(d$f) <- contr.sum(3)
d_na <- base_d(); d_na$x[2] <- NA; d_na$f[4] <- NA
custom_na <- function(object, ...) object[-1, , drop = FALSE]
for (na in list("na.omit", "na.exclude", na.omit, na.exclude, na.pass, na.fail, NULL, quote(na.omit),
                custom_na, "custom_na", "nonexistent", c("na.omit", "foo"), character(0), NA_character_,
                "na.action", 1)) {
    add(paste("na.action", paste(deparse(na), collapse = "")), d, y ~ f, na.action = na)
    add(paste("na.action, missing values,", paste(deparse(na), collapse = "")), d_na, y ~ x + f, na.action = na)
    add(paste("na.action, plain,", paste(deparse(na), collapse = "")), base_d(), y ~ x + f, na.action = na)
}
add("unknown variable and unknown na.action", d, y ~ zzz, na.action = "nonexistent")
d2 <- base_d(); attr(d2, "na.action") <- "na.exclude"; d2$x[2] <- NA; add("na.action attribute, string", d2)
d2 <- base_d(); attr(d2, "na.action") <- na.omit; add("na.action attribute, function", d2)
d2 <- na.omit(d_na); add("na.action attribute left by na.omit()", d2)
add("data is a matrix", cbind(y = y, x = x)); add("data is a matrix, xlev", cbind(y = y, x = x), ~ x, xlev = list(x = "a"))
d <- base_d(); d$w <- c(1, 2, NA, 4, 5)
add("subset, weights, offset", d, y ~ x + f, subset = c(TRUE, FALSE, TRUE, TRUE, TRUE), weights = d$w, offset = x)
add("drop.unused.levels", d, y ~ x + f, subset = c(TRUE, TRUE, TRUE, FALSE, TRUE), drop.unused.levels = TRUE)

# xlev ------------------------------------------------------------------------
XL <- list(f = c("a", "b", "c"))
for (na in list(na.omit, na.pass, na.exclude, na.fail, NULL)) {
    add("xlev", base_d(), ~ x + f, xlev = XL, na.action = na)
    add("xlev, two factors", base_d(), ~ x + f + g, xlev = list(f = c("a", "b", "c"), g = c("u", "v")), na.action = na)
}
d <- base_d(); d$x[2] <- NA; add("xlev, rows omitted", d, ~ x + f, xlev = XL)
d <- base_d(); d$f[4] <- NA
add("xlev, NA in the factor", d, ~ x + f, xlev = XL); add("xlev, NA in the factor, na.pass", d, ~ x + f, xlev = XL, na.action = na.pass)
d <- base_d(); class(d$f) <- c("special", "factor"); add("xlev, factor subclass", d, ~ x + f, xlev = XL, na.action = na.pass)
d <- base_d(); d$f <- I(d$f)
add("xlev, I(factor)", d, ~ x + f, xlev = XL, na.action = na.pass); add("xlev, I(factor), na.omit", d, ~ x + f, xlev = XL)
d <- base_d(); class(d$f) <- c("factor", "ordered"); add("xlev, class c(factor, ordered)", d, ~ x + f, xlev = XL, na.action = na.pass)
d <- base_d(); d$f <- factor(d$f, ordered = TRUE)
add("xlev, ordered", d, ~ x + f, xlev = XL, na.action = na.pass); add("xlev, ordered, na.omit", d, ~ x + f, xlev = XL)
d <- base_d(); d$f <- structure(c(1L, 2L, 1L, 4L, 0L), levels = c("a", "b", "c"), class = "factor")
add("xlev, codes out of range", d, ~ x + f, xlev = XL, na.action = na.pass)
d <- base_d(); attr(d$f, "levels") <- c("a", "b", "b")
add("xlev, duplicated levels", d, ~ x + f, xlev = list(f = c("a", "b", "b")), na.action = na.pass)
d <- base_d(); attr(d$f, "levels") <- c(p = "a", q = "b", r = "c")
add("xlev, named levels", d, ~ x + f, xlev = list(f = c(p = "a", q = "b", r = "c")), na.action = na.pass)
d <- base_d(); attr(d$f, "levels") <- 1:3; add("xlev, integer levels", d, ~ x + f, xlev = list(f = 1:3), na.action = na.pass)
d <- base_d(); attr(d$f, "levels") <- structure(c("a", "b", "c"), foo = 1)
add("xlev, levels carrying an attribute", d, ~ x + f, xlev = list(f = structure(c("a", "b", "c"), foo = 1)), na.action = na.pass)
add("xlev carrying an attribute", base_d(), ~ x + f, xlev = list(f = structure(c("a", "b", "c"), foo = 1)), na.action = na.pass)
d <- base_d(); levels(d$f)[3] <- NA; add("xlev, NA level", d, ~ x + f, xlev = list(f = c("a", "b", NA)), na.action = na.pass)
d <- base_d(); d$f <- structure(d$f, names = letters[1:5])
add("xlev, named factor", list(x = x, f = d$f), ~ x + f, xlev = XL, na.action = na.pass)
d <- base_d()
add("xlev, duplicated names", d, ~ x + f, xlev = list(f = c("a", "b", "c"), f = c("c", "b", "a")), na.action = na.pass)
add("xlev, absent variable", d, ~ x + f, xlev = list(zz = c("a", "b")), na.action = na.pass)
add("xlev, numeric variable", d, ~ x + f, xlev = list(x = c("a", "b")), na.action = na.pass)
add("xlev, new level", d, ~ x + f, xlev = list(f = c("a", "b")), na.action = na.pass)
add("xlev, other order", d, ~ x + f, xlev = list(f = c("c", "b", "a")), na.action = na.pass)
add("xlev, more levels", d, ~ x + f, xlev = list(f = c("a", "b", "c", "d")), na.action = na.pass)
add("xlev, character vector", d, ~ x + f, xlev = c(f = "a"), na.action = na.pass)
add("xlev, no names", d, ~ x + f, xlev = list(c("a", "b", "c")), na.action = na.pass)
add("xlev and drop.unused.levels", d, ~ x + f + g, xlev = XL, drop.unused.levels = TRUE, na.action = na.pass)
d <- base_d(); contrasts(d$f) <- contr.sum(3)
add("xlev, contrasts dropped", d, ~ x + f, xlev = XL, na.action = na.pass); add("xlev, contrasts dropped, na.omit", d, ~ x + f, xlev = XL)
d <- base_d(); d$ch <- c("p", "q", "p", "q", "r")
add("xlev, character variable", d, ~ ch, xlev = list(ch = c("p", "q", "r")), na.action = na.pass)
add("xlev, character variable, other order", d, ~ ch, xlev = list(ch = c("r", "q", "p")))
long <- ~ I(x + 1.5 + 2.5 + 3.5 + 4.5 + 5.5 + 6.5 + 7.5 + 8.5 + 9.5 + 10.5 + 11.5 + 12.5 + 13.5 + 14.5 + 15.5) + f
add("xlev, long expression", base_d(), long, xlev = XL)
add("xlev, factor(f)", base_d(), ~ factor(f), xlev = list("factor(f)" = c("a", "b", "c")), na.action = na.pass)
na_contrasts <- function(object, ...) { contrasts(object$f) <- contr.sum(3); object }
add("xlev, na.action that sets contrasts", base_d(), ~ x + f, xlev = XL, na.action = na_contrasts)
na_character <- function(object, ...) { object$f <- as.character(object$f); object }
add("xlev, na.action that makes a character", base_d(), ~ x + f, xlev = XL, na.action = na_character)
# the same levels marked UTF-8 and latin1 (made from bytes: any locale)
e_utf8 <- rawToChar(as.raw(c(0xc3, 0xa9))); Encoding(e_utf8) <- "UTF-8"
e_latin1 <- rawToChar(as.raw(0xe9)); Encoding(e_latin1) <- "latin1"
lv_utf8 <- c("a", e_utf8); lv_latin1 <- c("a", e_latin1)
expect_true(identical(lv_utf8, lv_latin1) && !identical(Encoding(lv_utf8), Encoding(lv_latin1)))
de <- data.frame(x = 1:4 + 0, f = factor(c(1, 2, 1, 2), labels = lv_utf8))
add("xlev, other encoding", de, ~ f, xlev = list(f = lv_latin1), na.action = na.pass)

for (cs in cases) {
    expect_true(do.call(same_frame, c(list(cs$f, cs$d), cs$args)), info = cs$label)
}

# identical() takes the same string in two encodings for one: compare the marks
fe <- ~ f
expect_identical(Encoding(levels(design_frame(fe, de, xlev = list(f = lv_latin1))$f)),
                 Encoding(levels(stats::model.frame(fe, de, xlev = list(f = lv_latin1))$f)))

# Calls as users write them: non-standard evaluation, promises, dispatch ------
f1 <- y ~ x + f
d <- base_d(); d$w <- 1:5 + 0
expect_true(strict_identical(design_frame(f1, d, subset = x > 1, weights = w, offset = log(w)),
                             stats::model.frame(f1, d, subset = x > 1, weights = w, offset = log(w))))
expect_true(strict_identical(design_frame(formula = f1, data = d), stats::model.frame(formula = f1, data = d)))
expect_true(strict_identical(design_frame(data = d, f1), stats::model.frame(data = d, f1)))
fx <- ~ x
expect_true(strict_identical(design_frame(fx), stats::model.frame(fx)))

# the 'newdata' warning depends on the expression given as `data`
newdata <- data.frame(q = 1:3)
expect_true(same_outcome(capture(function() stats::model.frame(fx, newdata)),
                         capture(function() design_frame(fx, newdata))))
expect_true(length(capture(function() design_frame(fx, newdata))$warnings) == 1L)

# arguments are forced when model.frame() forces them
expect_true(same_outcome(capture(function() stats::model.frame(y ~ zzz, d, na.action = stop("boom"))),
                         capture(function() design_frame(y ~ zzz, d, na.action = stop("boom")))))
dn <- d; dn$x[1] <- NA
expect_true(same_outcome(capture(function() stats::model.frame(f1, dn, na.action = na.fail, xlev = stop("boom"))),
                         capture(function() design_frame(f1, dn, na.action = na.fail, xlev = stop("boom")))))
seen <- character()
note <- function(what, value) { seen <<- c(seen, what); value }
invisible(design_frame(note("formula", f1), note("data", d), xlev = note("xlev", NULL), na.action = note("na", na.omit)))
order_ours <- seen; seen <- character()
invisible(stats::model.frame(note("formula", f1), note("data", d), xlev = note("xlev", NULL), na.action = note("na", na.omit)))
expect_identical(order_ours, seen)
expect_true(same_outcome(capture(function() stats::model.frame()), capture(function() design_frame())))

# the option na.action: a function, a string, the name of a user's function
for (opt in list(na.omit, "na.exclude", "custom_na", NULL, "nonexistent")) {
    options(na.action = opt)
    expect_true(same_frame(f1, dn), info = paste("option na.action", paste(deparse(opt), collapse = "")))
    expect_true(same_frame(f1, d), info = paste("option na.action", paste(deparse(opt), collapse = "")))
}
options(na.action = "na.omit")

# a na.omit() method for data frames defined by the user is not the one
# model.frame() calls
na.omit.data.frame <- function(object, ...) stop("not this one")
expect_true(same_frame(f1, dn))
rm(na.omit.data.frame)

# other methods of model.frame(): fitted models, methods in the caller's scope
fit <- lm(y ~ x + f, dn, na.action = na.exclude)
expect_true(strict_identical(design_frame(fit), stats::model.frame(fit)))
expect_true(strict_identical(design_frame(fit, data = dn), stats::model.frame(fit, data = dn)))
calls <- 0L
model.frame.xtx_test_formula <- function(formula, ...) {
    calls <<- calls + 1L
    if (!identical(substitute(formula), quote(ff))) stop("not the caller's expression")
    NextMethod()
}
ff <- structure(y ~ x, class = c("xtx_test_formula", "formula"))
expect_true(strict_identical(design_frame(ff, d), stats::model.frame(ff, d)))
expect_identical(calls, 2L)
local({
    model.frame.default <- function(formula, ...) "the caller's default method"
    expect_identical(design_frame(f1, d), "the caller's default method")
    expect_identical(stats::model.frame(f1, d), "the caller's default method")
})

# option: without the shortcuts, model.frame() itself
options(xtx.skip_na_copy = FALSE)
for (cs in cases[seq(1, length(cases), by = 7)]) {
    expect_true(do.call(same_frame, c(list(cs$f, cs$d), cs$args)), info = cs$label)
}
options(xtx.skip_na_copy = NULL)

# Through design_matrix() and plans ---------------------------------------------
same_matrix <- function(f, d, ...) {
    args <- list(f, d, ...)
    same_outcome(capture(function() do.call(stats::model.matrix, args)),
                 capture(function() as.matrix(do.call(design_matrix, args))))
}
# the row names of the response, made unique by `[.data.frame` (wrong before
# design_frame(): "a" "a" "b")
yy <- c(a = 1, a = 2, b = 3); xv <- c(1, 2, 4)
expect_identical(rownames(as.matrix(design_matrix(yy ~ xv))), c("a", "a.1", "b"))
expect_true(strict_identical(as.matrix(design_matrix(yy ~ xv)), stats::model.matrix(yy ~ xv)))
for (cs in cases) {
    if (length(cs$args)) next
    expect_true(same_matrix(cs$f, cs$d), info = paste("design_matrix:", cs$label))
}
set.seed(11)
tr <- fuzz_data(60)
nd <- fuzz_data(40, na = TRUE)
fo <- y ~ x1 * f1 + o1 + ch + poly(z, 2)
plan <- xtx:::design_plan(fo, tr)
for (na in list(na.pass, na.omit, na.exclude, "na.omit")) {
    for (newd in list(nd, tr, tr[sample(60), ])) {
        ref <- capture(function() stats::model.matrix(
            plan$terms, stats::model.frame(plan$terms, newd, xlev = plan$xlevels, na.action = na),
            contrasts.arg = plan$contrasts))
        expect_true(same_outcome(capture(function() as.matrix(design_matrix(plan, newd, na.action = na))), ref))
    }
}

# Which path ----------------------------------------------------------------------
# The wrapper of na.omit() calls the real method on the whole frame, on a
# frame without columns (for the row names), or not at all.
seen <- integer()
counting <- function(object, ...) {
    seen <<- c(seen, length(object))
    stats::na.omit(object, ...)
}
wrap <- xtx:::xtx_wrap_na(counting, stats::na.omit)
frame_of <- function(d, f = y ~ x + f) {
    mf <- stats::model.frame(f, d, na.action = NULL)
    attr(mf, "terms") <- NULL
    mf
}
path <- function(d, ...) {
    seen <<- integer()
    mf <- frame_of(d, ...)
    expect_true(strict_identical(wrap(mf), stats::na.omit(mf)))
    seen
}
d <- base_d()
expect_identical(path(d), integer())                           # automatic row names
expect_identical(path(d[c(1, 3, 4), ]), integer())             # integers, all different
expect_identical(path(d[c(4, 1, 3), ]), integer())
d2 <- d; attr(d2, "row.names") <- c(-5L, 2000000000L, 7L, 0L, -2000000000L)
expect_identical(path(d2), 0L)                                 # too far apart to tell quickly
attr(d2, "row.names") <- c(3L, 1L, 3L, 2L, 5L); expect_identical(path(d2), 0L)    # duplicated
attr(d2, "row.names") <- c(3L, 1L, NA, 2L, 5L); expect_identical(path(d2), 0L)    # NA
rownames(d) <- letters[1:5]; expect_identical(path(d), integer())    # character, all different
attr(d, "row.names") <- c("a", "b", "a", "c", "d"); expect_identical(path(d), 0L)    # duplicated
attr(d, "row.names") <- c("a", "b", NA, "c", "d"); expect_identical(path(d), 0L)     # NA
attr(d, "row.names") <- c("a", "b", "c", lv_utf8[2], lv_latin1[2]); expect_identical(path(d), 0L)  # equal after translation
d <- base_d(); d$x[2] <- NA; expect_identical(path(d), 3L)     # a missing value
d <- base_d(); attr(d$x, "label") <- "L"; expect_identical(path(d), 3L)
d <- base_d(); contrasts(d$f) <- contr.sum(3); expect_identical(path(d), 3L)
d <- base_d(); d$x <- as.Date("2020-01-01") + 0:4; expect_identical(path(d), 3L)
expect_identical(path(base_d()[0, ]), 3L)                      # no rows

int_unique <- function(x) .Call(xtx:::C_xtx_int_unique, x)
expect_true(int_unique(c(5L, 1L, 3L)))
expect_true(int_unique(integer()))
expect_true(int_unique(c(-3L, 0L, 3L, -.Machine$integer.max)[1:3]))
expect_false(int_unique(c(5L, 1L, 5L)))
expect_false(int_unique(c(5L, NA, 1L)))
expect_false(int_unique(c(.Machine$integer.max, -.Machine$integer.max, NA)))
expect_identical(int_unique(c(.Machine$integer.max, -.Machine$integer.max)), NA)
expect_identical(int_unique(c(1, 2)), NA)
set.seed(3)
for (k in 1:50) {
    v <- sample(-40:40, sample(1:30, 1), replace = runif(1) < .5)
    expect_identical(int_unique(v), !anyDuplicated(v))
}

str_unique <- function(x) .Call(xtx:::C_xtx_str_unique, x)
expect_true(str_unique(c("b", "a", "")))
expect_true(str_unique(character()))
expect_true(str_unique(c(e_utf8, paste0(e_utf8, e_utf8))))
expect_false(str_unique(c("b", "a", "b")))
expect_false(str_unique(c("b", NA)))
expect_false(str_unique(c("a", "NA", NA)))
expect_identical(str_unique(c("a", lv_utf8[2])), NA)          # marks differ: not decided
expect_identical(str_unique(c(lv_utf8[2], lv_latin1[2])), NA) # the same string, twice
expect_identical(str_unique(1:3), NA)
for (k in 1:50) {
    v <- sample(c(letters, LETTERS, "", "NA"), sample(1:40, 1), replace = runif(1) < .5)
    expect_identical(str_unique(v), !anyDuplicated(v))
}
v <- as.character(sample(1e5)); expect_true(str_unique(v)); v[77777] <- v[3]; expect_false(str_unique(v))

same_factor <- xtx:::xtx_same_factor
f <- factor(c("a", "b", NA, "a"), levels = c("a", "b", "c"))
expect_true(same_factor(f, c("a", "b", "c")))
expect_true(same_factor(factor(f, levels = c("a", "b", "c"), ordered = TRUE), c("a", "b", "c")))
expect_true(same_factor(f[0], c("a", "b", "c")))
expect_false(same_factor(f, c("a", "c", "b")))
expect_false(same_factor(f, c("a", "b")))
expect_false(same_factor(f, factor(c("a", "b", "c"))))
expect_false(same_factor(structure(f, class = c("special", "factor")), c("a", "b", "c")))
expect_false(same_factor(I(f), c("a", "b", "c")))
expect_false(same_factor(structure(f, names = letters[1:4]), c("a", "b", "c")))
expect_false(same_factor(C(f, contr.sum), c("a", "b", "c")))
expect_false(same_factor(structure(c(1L, 4L), levels = c("a", "b", "c"), class = "factor"), c("a", "b", "c")))
expect_false(same_factor(structure(c(1L, 0L), levels = c("a", "b", "c"), class = "factor"), c("a", "b", "c")))
expect_false(same_factor(factor(c("a", NA), levels = c("a", NA), exclude = NULL), c("a", NA)))
expect_false(same_factor(structure(1:2, levels = c("a", "a"), class = "factor"), c("a", "a")))
expect_false(same_factor(factor(c(1, 2), labels = lv_utf8), lv_latin1))
expect_true(same_factor(factor(c(1, 2), labels = lv_utf8), lv_utf8))
expect_false(same_factor(c("a", "b"), c("a", "b")))

# A version of R with other code ----------------------------------------------------
# Each patch is applied only to the statements it was written for.
make <- xtx:::xtx_make_frame
replace <- xtx:::xtx_replace
count <- function() { e <- new.env(); e$n <- 0L; e }
if (isTRUE(fr$na) && isTRUE(fr$xlev)) {
    mfd <- stats::model.frame.default
    patches <- function(...) { p <- make(...); c(na = isTRUE(p$na), xlev = isTRUE(p$xlev)) }
    expect_identical(patches(), c(na = TRUE, xlev = TRUE))
    expect_identical(patches(na = FALSE), c(na = FALSE, xlev = TRUE))
    expect_identical(patches(xlev = FALSE), c(na = TRUE, xlev = FALSE))
    expect_null(make(na = FALSE, xlev = FALSE))

    # the xlev loop reads differently
    alt <- mfd
    body(alt) <- replace(body(mfd), quote(factor(xi, levels = xl, exclude = NULL)), quote(factor(xi, levels = xl)), count())
    expect_identical(patches(alt), c(na = TRUE, xlev = FALSE))
    # the C routine takes another argument, or is called twice
    c_call <- quote(.External2(C_modelframe, formula, rownames, variables, varnames, extras, extranames, subset, na.action))
    alt <- mfd
    body(alt) <- replace(body(mfd), c_call, as.call(c(as.list(c_call), list(NULL))), count())
    expect_identical(patches(alt), c(na = FALSE, xlev = TRUE))
    alt <- mfd
    body(alt) <- as.call(c(as.list(body(mfd)), list(call("<-", quote(data), c_call))))
    expect_identical(patches(alt), c(na = FALSE, xlev = TRUE))
    # na.omit(), na.exclude() or `[.data.frame` read differently
    omit <- stats:::na.omit.data.frame
    alt <- omit
    body(alt) <- replace(body(omit), quote(xx <- object[!omit, , drop = FALSE]), quote(xx <- object[!omit, ]), count())
    expect_identical(patches(omit = alt), c(na = FALSE, xlev = TRUE))
    expect_identical(patches(exclude = alt), c(na = FALSE, xlev = TRUE))
    expect_identical(patches(exclude = omit), c(na = FALSE, xlev = TRUE))
    sub <- base::`[.data.frame`
    alt <- sub
    body(alt) <- replace(body(sub), quote(rows <- rows[i]), quote(rows <- unname(rows[i])), count())
    expect_identical(patches(subset = alt), c(na = FALSE, xlev = TRUE))
    # nothing recognised
    expect_null(make(function(formula, ...) NULL))
    expect_null(make(function(formula, ...) { formula }))
    expect_null(make(NULL))

    # a function with one patch still returns what model.frame() returns
    for (one in list(make(na = FALSE)$fun, make(xlev = FALSE)$fun)) {
        for (cs in cases[seq(1, length(cases), by = 5)]) {
            args <- c(list(cs$f, cs$d), cs$args)
            expect_true(same_outcome(capture(function() do.call(stats::model.frame, args)),
                                     capture(function() do.call(one, args))), info = cs$label)
        }
    }
}

# The self-test at load rejects a function that returns anything else.
selftest <- xtx:::xtx_frame_selftest
expect_true(selftest(stats::model.frame.default))
if (!is.null(fr)) expect_true(selftest(fr$fun))
expect_false(selftest(function(...) { mf <- stats::model.frame(...); attr(mf, "note") <- 1; mf }))
expect_false(selftest(function(...) { mf <- stats::model.frame(...); attributes(mf) <- rev(attributes(mf)); mf }))
expect_false(selftest(function(...) stats::model.frame(..., drop.unused.levels = TRUE)))
expect_false(selftest(function(...) stop("no")))
expect_false(selftest(function(...) { warning("careful"); stats::model.frame(...) }))

# Random frames ---------------------------------------------------------------------
# (a short run of the generator used to review design_frame())
set.seed(20261125)
rcol <- function(n) {
    lev <- sample(c("a", "b", "c", "d"), sample(2:4, 1))
    fac <- function() factor(sample(lev, n, TRUE), levels = if (runif(1) < .5) lev else c(lev, "zz"))
    kind <- sample(c("num", "num", "int", "lgl", "chr", "fac", "fac", "ord", "matdn", "asis", "label",
                     "contr", "named", "subclass", "date", "tsp", "arr3", "mat0", "asis_fac"), 1)
    switch(kind,
           num = rnorm(n), int = sample(-3:3, n, TRUE), lgl = sample(c(TRUE, FALSE), n, TRUE),
           chr = sample(c("u", "v", "w"), n, TRUE), fac = fac(),
           ord = factor(sample(lev, n, TRUE), levels = lev, ordered = TRUE),
           matdn = matrix(rnorm(2 * n), n, 2, dimnames = list(NULL, c("A", "B"))),
           asis = I(rnorm(n)), label = structure(rnorm(n), label = "L"),
           contr = { f <- factor(sample(lev, n, TRUE), levels = lev); contrasts(f) <- contr.sum(length(lev)); f },
           named = structure(fac(), names = sprintf("n%d", seq_len(n))),
           subclass = structure(fac(), class = c("special", "factor")),
           date = as.Date("2020-01-01") + sample(100, n, TRUE),
           tsp = if (n) unclass(ts(rnorm(n))) else numeric(),
           arr3 = array(rnorm(4 * n), c(n, 2, 2)), mat0 = matrix(numeric(0), n, 0),
           asis_fac = I(fac()))
}
rdata <- function() {
    n <- sample(c(0, 1, 2, 3, 5, 8), 1, prob = c(.05, .1, .1, .15, .3, .3))
    k <- sample(2:4, 1)
    cols <- lapply(seq_len(k), function(j) {
        x <- rcol(n)
        if (n && runif(1) < .15 && is.null(dim(x))) x[sample(n, 1)] <- NA
        x
    })
    names(cols) <- paste0("v", seq_len(k))
    form <- sample(c("df", "df", "df", "list", "env"), 1)
    if (form == "env") return(list2env(cols))
    if (form == "list") {
        if (n && is.null(attributes(cols[[1]]))) names(cols[[1]]) <- sample(c("p", "q", "r"), n, TRUE)
        return(cols)
    }
    class(cols) <- "data.frame"
    attr(cols, "row.names") <- switch(sample(5, 1), .set_row_names(n), .set_row_names(n), paste0("r", sample(1e6, n)),
                                      sort(sample(1e6, n)), sample(c("a", "b", "c"), n, TRUE))
    cols
}
nas <- list("missing", "missing", "na.omit", na.omit, na.exclude, na.pass, na.fail, NULL, custom_na)
ok <- TRUE
for (it in 1:150) {
    d <- rdata()
    nm <- if (is.environment(d)) ls(d) else names(d)
    rhs <- paste(sample(nm, sample(seq_along(nm), 1)), collapse = sample(c(" + ", ":"), 1))
    f <- as.formula(paste(if (runif(1) < .6) sample(nm, 1) else "", "~", rhs), env = globalenv())
    args <- list(f, d)
    if (runif(1) < .6) {
        mf <- tryCatch(suppressWarnings(stats::model.frame(f, d, na.action = na.pass)), error = function(e) NULL)
        xl <- list()
        for (v in names(mf)) {
            col <- mf[[v]]
            if (is.factor(col) || is.character(col)) {
                lv <- if (is.factor(col)) levels(col) else sort(unique(col))
                xl[[v]] <- if (runif(1) < .75) lv else rev(lv)
            }
        }
        if (length(xl)) args$xlev <- xl
    }
    na <- sample(nas, 1)[[1]]
    if (!identical(na, "missing")) args["na.action"] <- list(na)
    if (runif(1) < .15 && is.data.frame(d) && nrow(d) > 2) args$subset <- sample(c(TRUE, FALSE), nrow(d), TRUE)
    same <- same_outcome(capture(function() do.call(stats::model.frame, args)),
                         capture(function() do.call(design_frame, args)))
    if (!same) {
        ok <- FALSE
        break
    }
}
expect_true(ok, info = if (ok) "" else paste("random frame", it, deparse(f)))

options(old_opt)

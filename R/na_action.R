## na.action of the fits, in one place: which na.action functions a fit
## accepts (the same in memory and remotely), the na.action a fit keeps
## for its results, and how results on the rows used are padded to the rows
## of the data (napredict() and naresid(): in memory on vectors, remotely
## on tables, with the methods of xtx_remote_exclude).

## na.action of a fit, as one of `allowed` ("omit", "exclude", "fail",
## "pass"): na.omit and the others, as functions or names. Other functions
## are refused, in memory as remotely (remote data can apply only these).
xtx_na_kind <- function(na.action, allowed = c("omit", "exclude", "fail", "pass")) {
    if (is.null(na.action)) return("omit")
    for (nm in allowed) {
        f <- get(paste0("na.", nm), envir = asNamespace("stats"))
        if (identical(na.action, f) || identical(na.action, paste0("na.", nm))) return(nm)
    }
    stop(sprintf("`na.action` must be %s.", paste(paste0("na.", allowed), collapse = ", ")), call. = FALSE)
}

## na.action of a fit with fixed effects: na.omit or na.exclude, in memory
## and remotely.
xtx_fe_na_kind <- function(na.action) {
    k <- tryCatch(xtx_na_kind(na.action), error = function(e) "")
    if (!k %in% c("omit", "exclude"))
        stop("With fixed effects, `na.action` must be na.omit or na.exclude.", call. = FALSE)
    k
}

## The na.action a fit keeps for its results, from its input stage P
## (design_problem()): remotely with na.exclude, xtx_remote_exclude() of x
## (the design of the rows the fit could use, over the data); in memory,
## the model frame's na.action, with the rows of the frame that the fit did
## not use (`rows`: the rows used, with fixed effects; NULL: all of them).
## NULL for a fit on a design (no rows dropped).
xtx_na_result <- function(P, rows = NULL, x = P$x) {
    if (isTRUE(P$exclude)) return(xtx_remote_exclude(x))
    xtx_na_dropped(P$na.action, rows, P$nall)
}

## The rows dropped: those of the model frame's na.action and the rows of
## the frame not used (zero weights, singletons with fixed effects), as
## indices in the data, with the class of the na.action ("omit" by
## default).
xtx_na_dropped <- function(na_act, rows, nall) {
    if (is.null(rows) || length(rows) == nall) return(na_act)
    orig <- if (is.null(na_act)) seq_len(nall) else seq_len(nall + length(na_act))[-na_act]
    dropped <- sort(c(as.integer(na_act), orig[-rows]))
    names(dropped) <- as.character(dropped)
    class(dropped) <- if (is.null(na_act)) "omit" else class(na_act)
    dropped
}

## A result in memory, one value per row of design x: named by the rows of
## x and padded by pad(na_act) (stats::napredict, or stats::naresid for
## residuals), as lm() pads them.
xtx_na_pad <- function(na_act, v, x, pad = stats::napredict) {
    names(v) <- dimnames(x)[[1L]]
    pad(na_act, v)
}

## na.exclude for a remote fit: results on the rows used, padded with the
## rows of the data that the fit dropped for missing values (NULL results),
## as napredict() pads in memory (see napredict.xtx_remote_exclude()).
xtx_remote_exclude <- function(D) structure(list(design = D), class = "xtx_remote_exclude")

#' @exportS3Method stats::napredict
napredict.xtx_remote_exclude <- function(omit, x, ...) {
    if (!inherits(x, "xtx_remote_table")) return(x)
    D0 <- omit$design
    con <- attr(x, "connection")
    sel <- sprintf("SELECT * FROM %s UNION ALL BY NAME SELECT s.* FROM (%s) AS s WHERE NOT COALESCE(%s, FALSE)",
                   DBI::dbQuoteIdentifier(con, unclass(x)[1L]), D0$src$select, D0$filter)
    name <- xtx_tmp_name("xtx_out")
    DBI::dbExecute(con, sprintf("CREATE TEMPORARY VIEW %s AS %s", DBI::dbQuoteIdentifier(con, name), sel))
    structure(name, class = "xtx_remote_table", connection = con, holder = attr(x, "holder"))
}

#' @exportS3Method stats::naresid
naresid.xtx_remote_exclude <- function(omit, x, ...) napredict.xtx_remote_exclude(omit, x, ...)

## The primitive layer of the fitting core (Kata design.matrix#vavc).
##
## Fitting code (linear models, and on this layer GLMs and fixed effects) is
## written once against the operations below, which take a *design* and
## *columns*, and dispatch to two backends:
##
##   design                     column
##   ------------------------   -----------------------------------------
##   in memory: design_matrix   numeric vector (or matrix: several columns)
##   remote:    dm_remote       design_column: SQL expressions over the rows
##                              of a frame (see below; made by dm_column())
##
## The remote column class, design_column, is also the public column type
## (R/column.R: design_column(), as.vector(), as.matrix(), col_compute(),
## as.data.frame(), print/format, and the methods below for arithmetic). Everything else
## here is internal; names start with dmp_ (operations on designs), dmc_
## (operations on columns) or dm_.
##
## ---------------------------------------------------------------------------
## Columns
## ---------------------------------------------------------------------------
##
## A column holds one value per row of a design (in the order of its rows in
## memory; remotely, rows have no order). In memory it is a double vector of
## length nrow(x), or an n x k matrix for k columns. Remotely it is an
## object of class design_column, made by
##
##   dm_column(frame, sql)
##
## with `frame` a dm_remote (its source `s` holds what the expressions read,
## its filter selects its rows) and `sql` a character vector, one SQL
## expression over `s` per column (k = length(sql)). Columns are made by
##
##   dmp_column(x, sql)    a column of design x from SQL over its source
##                         (in memory: `sql` is returned, i.e. pass a vector)
##   dmp_mult(x, B)        products with the design (below)
##   arithmetic on columns
##
## Operations on columns (same R code for both backends):
##   +  -  *  /  ^         (numbers recycle; columns must have the same k or 1)
##   ==  !=  <  >  <=  >=  &  |  !
##                         (NA where either side is NA or NaN, as in R)
##   is.na is.nan is.finite is.infinite (R's semantics; never NA)
##   exp log sqrt abs log1p expm1 log2 log10 cos sin tan floor ceiling sign
##                         (Math group; log(x, base) with a number or a
##                         column as base, as R: dm_column_log())
##   dmc_pmin dmc_pmax dmc_ifelse dmc_rowsums dmc_cbind dmc_col(x, j) dm_ncols
##   dmc_const(like, value) dmc_lgamma dmc_round (half to even)
##   dmc_is_finite (NaN-safe) dmc_present (not NA or NaN: !is.na)
##   dmc_on(col, x) (same SQL, rows of x)
##   dmc_count(cond, negate) (1 where cond is TRUE, to sum)
##   dmc_freeze(col, table) (computed once: remotely, columns of a new
##                         frame; with table = TRUE, a temporary table of the
##                         rows used with every column of the frame)
##   dmc_pnorm dmc_qnorm dmc_dnorm dmc_pcauchy dmc_qcauchy
##   dmc_logit_inv dmc_logit_mu_eta dmc_binomial_dev(y, mu, w)
##                         (in memory the C code of binomial())
##                         (remotely the DuckDB macros of dm_sql_macros(), as
##                         R computes them: pnorm/dnorm/pcauchy to rounding,
##                         qnorm to one ulp; dmc_mexpm1m(u) = -expm1(-u))
##
## NA semantics: SQL NULL is R's NA, IEEE NaN is NaN; arithmetic propagates
## both as R does (NULL op x is NULL). Sums over rows (dm_colsum) return NA
## when a value is NULL, as sum() does, and NaN when a value is NaN.
## Comparisons with NULL are NULL (R: NA).
##
## Products (dmp_mult, dmp_quad) and lookups (dmc_lookup, dmc_dense_codes)
## are *pending* columns: a recipe (the design and B of a product; the
## lookup table, already written, and its join of a lookup) that becomes a
## new frame (a view, or a join) only at the column's first use: an
## operation that reads rows (sums, crossprod, groups, dmp_collect,
## as.vector, as.matrix, as.data.frame, output), or that makes a frame or a design from the
## column (dmc_freeze, col_compute, scale, subset, materialize). Arithmetic on
## pending columns stays pending. At first use the pending columns of the
## operands are bound on the latest settled frame among them, one on the
## other in the order of the arguments, so that they all combine (a frame
## keeps every column of the frame below it); each recipe is bound once
## (later uses read the same frame). So `D %*% b1 + D %*% b2` stacks the
## second product on the first. Combining columns takes the later frame
## (they must come from the same design, and one frame must descend from
## the other: each frame records its lineage, see dm_frame_later()). For
## columns already settled on two branches of one frame (two products each
## read before they are combined), the products and lookups of one branch
## are replayed on the other (copies of their recipes, bound again: one
## more view per product, none for a lookup; dm_frame_replay()); columns
## frozen or stored (dmc_freeze, col_compute) on each branch cannot be
## combined. dmc_bind() binds columns together that are used one by one
## but must combine later (no replay needed). In memory a product is
## computed at once.
##
## Getting values into R: dmp_collect(col) returns the values (in memory the
## vector itself; remotely a SELECT of all rows: only for small data). The
## public as.vector() and as.matrix() call it, and col_compute() calls
## dmc_freeze(table = TRUE) (with a name, a permanent table of the same
## rows and columns).
##
## ---------------------------------------------------------------------------
## Operations on designs
## ---------------------------------------------------------------------------
##
##   dmp_scale(x, s)          the design with rows scaled by column s (s * x);
##                            scales multiply. Same rows and columns as x.
##   dmp_crossprod(x)         X'X, p x p matrix
##   dmp_crossprod(x, y)      y a design (x itself, possibly scaled: X'SY);
##                            or a column (k columns): X'c, p x k matrix
##   dmp_colsums(x, group)    column sums of x (possibly scaled): vector of
##                            length p; with group (a column, values used as
##                            group labels; NA a group of its own, with
##                            rowsum()'s warning): G x p matrix, groups sorted
##                            as rowsum() sorts them, row names the labels
##   dmp_groupsum(col, group) sums of a column (k columns) within groups:
##                            G x k matrix, as rowsum(col, group) (NA for a
##                            group with a missing value)
##   dmp_y(x, y)              the response of the qr() methods as a column
##   dmp_min_nonzero(x)       a lower bound on the non-zero entries (0 if none)
##   dmc_codes(codes, G)      a grouping as integer codes 1..G: grouped sums
##                            then return all G groups in code order, without
##                            labels, from one-pass C kernels in memory
##   dmp_colsums_groups(x, groups, extra)
##                            several coded groupings at once: per grouping,
##                            list(S = G x p, n = G sums of the scale, extra =
##                            G x k sums of scale * extra); one pass of the
##                            structured kernel in memory
##   dmp_fe_sums(x, w, groups, extra, pairs)
##                            dmp_colsums_groups() of w * x and the
##                            crosstabs of pairs of coded groupings at once
##                            (fewer passes in memory)
##   dmc_groupsums(col, groups, total)
##                            sums of col within each coded grouping (a G x k
##                            matrix each, as dmp_groupsum()) and of `total`
##                            over all rows (as dm_colsum()), in one query
##                            remotely
##   dmp_crosstab(a, b, w)    for two grouping columns, the sums of w (counts
##                            without w) for each pair of groups: a
##                            Matrix::dgCMatrix, rows the sorted groups of a,
##                            columns those of b (labels as.character() of
##                            the values, as in memory); with dmc_codes(), the
##                            Ga x Gb matrix in code order
##   dmp_mult(x, B)           X B (B a p-vector or p x k matrix) as a column
##                            (in memory a vector for a vector B; remotely
##                            pending, see above)
##   dmp_quad(x, V)           rowSums((X V) * X), V p x p symmetric positive
##                            semi-definite, as a column (remotely pending)
##   dmc_bind(cols)           the pending columns of a list bound together
##                            (one frame remotely); in memory the list
##   dmp_rebase(x, cols)      the design over the latest frame of the
##                            columns of a list (bound together); in
##                            memory x
##   dm_colsum(col)           sums over rows: a number per column
##   dmc_sums(...)            sums over rows of several columns (one query)
##   dmc_mean(col)            the mean of a column
##   dmp_colnorms(col)        Euclidean norms of the columns
##   dm_nrows(x)              number of rows of a design or column
##   dm_engine(x)             "sql" for a remote design, ""
##                            in memory (for the description of a fit)
##   dmp_patterns(x, z, w, compress, keys = NULL, fit = "lm")
##                            compressed least squares: when rows with the
##                            same design inputs are worth grouping, a list
##                            with x (the design of the patterns, in memory),
##                            z (the weighted mean of the column z per
##                            pattern), w (the sum of the weights per pattern,
##                            NULL when all 1), within (the weighted sum of
##                            squares of z around the pattern means) and rows
##                            (the number of rows); otherwise NULL. compress:
##                            "auto", TRUE or FALSE. A weighted least squares
##                            fit of z on x with weights w gives the fit on
##                            every row; its residual sum of squares plus
##                            `within` is that of the fit on every row. With
##                            `keys` (a column, k columns), patterns are also
##                            split by their values, returned as `keys`
##                            (one row per pattern). fit: "lm" or "glm",
##                            for the remote rule of "auto"
##                            (dm_compress_rule(); a weighted "lm" is
##                            "wlm").
##   dmp_materialize(x, cols = NULL)
##                            the design over a copy of its rows: remotely a
##                            temporary table of the rows used with every
##                            column of its frame (later queries scan it
##                            instead of recomputing the frame's views and
##                            joins), or, given `cols` (a list of the
##                            columns to re-make on it), the columns of the
##                            data and the internal ones (lookups) that the
##                            design and `cols` read; a new design
##                            (re-make columns with dmc_on()); in memory x
##   dmp_subset(x, keep)      the design of the rows where the logical column
##                            keep is TRUE (remotely also a DBI::SQL()
##                            predicate); a new remote design
##   dmp_resolve(x, exprs, env, logical, data)
##                            columns for the arguments of a fit (response,
##                            weights, offset, subset): column names, SQL(),
##                            R expressions (in SQL when translatable, else
##                            evaluated in R into a temporary table); returns
##                            list(x, cols, evaluated)
##   dmp_output(x, cols, memory, na, design, attrs)
##                            row-wise results of an algorithm: in memory
##                            memory(d), d the data frame of the named
##                            columns `cols` (one row per row of x; default
##                            d); remotely a temporary view with the columns
##                            of the remote data and `cols`, a string of
##                            class dm_remote_table (padded by napredict(na),
##                            with attributes attrs), or with design = TRUE
##                            the remote design with `cols` added to its
##                            data. One code path computes the columns; only
##                            the container differs.
##
## Remote methods: crossprod, colSums, products and patterns run in the
## database; dmp_groupsum and
## dmp_crosstab are single GROUP BY queries.
##
## ---------------------------------------------------------------------------
## User-facing operations of a remote design (operator goal: every matrix
## operation of an in-memory design works on a remote one, through this
## layer; status as of this contract)
## ---------------------------------------------------------------------------
##
##   operation                    remote meaning                    status
##   dim, dimnames                rows, columns                     done
##   as.matrix(D, force)          collected; refused above          done
##                                remote_max_rows unless force
##   D %*% b, D %*% B             column (design_column, k = ncol(B)) done
##   crossprod(D), (D, s * D),    p x p / p x k matrices            done
##   (D, "col" / SQL() / column)
##   crossprod(D, numeric w)      no row alignment remotely: error; done (error)
##                                use a column or s * D
##   colSums/colMeans(s * D,      vectors / G x p matrices          done
##     group = col)
##   rowSums/rowMeans(D)          column (labels .rowsum/.rowmean)  done
##   quadratic_diagonal(D, V)     column (label .quad)              done
##   design_column(D, "y"/SQL())  column of the remote data         done
##   as.vector, as.matrix,        values in R / data frame (with    done
##     as.data.frame,             `keep` columns) / temporary or
##     col_compute (of a column)  permanent table
##   qr(D)                        dmqr_chol on X'X with lm()'s      done
##                                rule, suspects from products
##   qr.coef, qr.qty, qr.R, solve numbers (as in memory)            done
##   qr.fitted, qr.resid, qr.qy,  tables (.fitted, .resid, .qy,     done
##     qr.Q                       .Q1...)
##   as(D, "CsparseMatrix")       collected, same guard as          done
##                                as.matrix (option to raise)
##   D[i, ]                       i a logical remote column or a    done
##                                DBI::SQL() predicate; numeric or
##                                logical vectors refused
##   s * D, D / s, -D             row scaling (number, column       done
##                                name, SQL(), column)
##   D1 + D2, D1 - D2             refused (operator decision)       done (error)
##   finite_rows(D) (internal)    a logical column                  done
##
## Tests: test_primitives.R (the layer), test_parity_ops.R (every operation
## above on the same data, in memory and remote: a table, Parquet).

## ---- remote columns ----

## (`notnull`: columns known never to be NULL, whose sums need no count of
## missing values. An optional `label` field, set by dm_column_label() in
## the public calls, names the values for as.data.frame().)
dm_column <- function(frame, sql, notnull = FALSE) {
    ## (parenthesized, so that expressions compose with any operator)
    sql <- as.character(sql)
    sql <- sprintf("(%s)", sql)
    structure(list(frame = frame, sql = sql, notnull = rep_len(notnull, length(sql))), class = "design_column")
}

dm_is_column <- function(x) inherits(x, "design_column")

## The number of columns k of a column: of a remote column its number of
## SQL expressions, in memory NCOL(). (dim() of a remote column is c(n, k)
## and length() is n * k, as for a matrix in memory: the fitting core counts
## columns with dm_ncols(), never with length().)
dm_ncols <- function(x) if (dm_is_column(x)) length(unclass(x)$sql) else NCOL(x)

## The number of rows of a remote column: the row count cached with its
## design (every frame keeps it), else one COUNT query, cached thereafter.
dm_column_nrows <- function(x) {
    fr <- unclass(x)$frame
    n <- dm_nrows(fr)
    if (!is.null(n)) return(n)
    dm_column_valid(x)
    ## (the rows of the settled frame: binding recipes adds columns, not rows)
    fr$pend <- NULL
    n <- as.numeric(dm_query(fr, sprintf("SELECT COUNT(*) AS n FROM (%s) AS s WHERE %s",
                                         fr$src$select, fr$filter))$n)
    if (is.environment(fr$cache)) fr$cache$n <- n
    n
}

#' @export
dim.design_column <- function(x) c(dm_column_nrows(x), dm_ncols(x))

#' @export
length.design_column <- function(x) prod(dim.design_column(x))

dm_frame_depth <- function(D) D$depth %||% 0L

## The lineage of a frame: the tokens of the frames derived one from the
## other since the design was made (a design not yet derived has none), so
## that a frame descends from another when its path extends the other's.
## NULL when unknown (a derived frame made before frames had paths).
dm_frame_path <- function(D) {
    if (!is.null(D$path)) return(D$path)
    if (dm_frame_depth(D) == 0L) character() else NULL
}

## Frame `fr`, derived from frame `parent`: its path is the parent's and a
## new token. `steps` maps the tokens of the path made by binding a recipe
## (a product or a lookup) to that recipe, so that the frame can be made
## again on another branch (dm_frame_replay()); other tokens (frozen or
## stored columns, columns added to a design) cannot be replayed.
dm_frame_child <- function(fr, parent, recipe = NULL) {
    p <- dm_frame_path(parent)
    if (is.null(p)) return(fr)
    dm_tmp$frames <- (dm_tmp$frames %||% 0) + 1
    tok <- sprintf("f%x", as.integer(dm_tmp$frames))
    fr$path <- c(p, tok)
    fr$steps <- c(parent$steps, if (!is.null(recipe)) stats::setNames(list(recipe), tok))
    fr
}

## Of two frames of the same design, the later one (which keeps the columns
## of the earlier one): the descendant. For frames on two branches (neither
## descends from the other), the steps of one branch since they split are
## replayed on top of the other (dm_frame_replay(): the result is pending
## with those recipes), else that is an error. Without a known lineage,
## the deeper frame. The frames of pending columns (dm_frame_pending())
## combine their recipes: the result is pending with the recipes of both
## (in order, each once), on the later of their settled frames.
dm_frame_later <- function(a, b) {
    if (is.null(a)) return(b)
    if (is.null(b)) return(a)
    pend <- dm_pend_union(a$pend, b$pend)
    a$pend <- NULL
    b$pend <- NULL
    out <- dm_frame_later_settled(a, b)
    ## (replayed recipes first: the pending ones may read their columns)
    pend <- dm_pend_union(out$pend, pend)
    out$pend <- NULL
    if (length(pend)) dm_frame_pending(out, pend) else out
}

## (dm_frame_later() of two settled frames; the result has a `pend` field
## only when one branch is replayed on the other)
dm_frame_later_settled <- function(a, b) {
    if (!identical(a$id, b$id))
        stop("Columns of different remote designs cannot be combined (a subset or a materialized design ",
             "is another design: make the columns again on it).", call. = FALSE)
    pa <- dm_frame_path(a)
    pb <- dm_frame_path(b)
    if (!is.null(pa) && !is.null(pb) && !identical(pa, pb)) {
        prefix <- function(u, v) length(u) < length(v) && identical(u, v[seq_along(u)])
        if (prefix(pa, pb)) return(b)
        if (prefix(pb, pa)) return(a)
        ## (products and lookups bind on top of one another at their first
        ## use, and a product or lookup already bound on a sibling frame is
        ## replayed: users meet this only with a stored or frozen column on
        ## each branch, as two columns each stored by col_compute())
        out <- dm_frame_replay(a, b)
        if (!is.null(out)) return(out)
        stop("Columns built on two different branches of the same remote design cannot be combined: ",
             "on each branch a column was stored (col_compute()) or frozen after the branches split, and ",
             "remote rows have no key to align them. Combine the columns before storing them (store one ",
             "column made from both, or one product `D %*% cbind(b1, b2)` and then `[, 1]` and `[, 2]`), or ",
             "read them with as.data.frame(keep = ) and match them in R; internally, freeze them together ",
             "(dmc_cbind()) or bind them together first (dmc_bind()).",
             call. = FALSE)
    }
    if (dm_frame_depth(b) > dm_frame_depth(a)) b else a
}

## Frames a and b of one design on two branches (paths known, neither a
## prefix of the other): the frame of one of them, pending with copies of
## the recipes that made the other branch since the two split (its steps
## past their common ancestor), when all those steps are products or
## lookups; NULL when neither branch can be replayed (a stored or frozen
## column on each). The copies (dm_recipe_replay()) bind at the next use as
## any pending recipe does, on top of the frame: one view per product, no
## statement for a lookup; a copy adds the same columns as its original
## (same names), so every column of the other branch reads the same SQL
## on the result. Steps whose columns the frame already has (an earlier
## replay of the same recipe) are left out. Of two replayable branches,
## the one with fewer products is replayed (a on a tie).
dm_frame_replay <- function(a, b) {
    pa <- a$path
    pb <- b$path
    m <- min(length(pa), length(pb))
    k <- 0L
    while (k < m && identical(pa[k + 1L], pb[k + 1L])) k <- k + 1L
    side <- function(f, p) {
        toks <- p[seq_along(p) > k]
        if (!all(toks %in% names(f$steps))) return(NULL)
        unname(f$steps[toks])
    }
    sa <- side(a, pa)
    sb <- side(b, pb)
    need <- function(st, onto) {
        have <- dm_frame_keys(onto)
        Filter(function(r) !dm_recipe_key(r) %in% have, st)
    }
    views <- function(st) sum(vapply(st, function(r) identical(r$kind, "mult"), NA))
    nb <- if (!is.null(sb)) need(sb, a)
    na <- if (!is.null(sa)) need(sa, b)
    if (!is.null(sb) && (is.null(sa) || views(nb) <= views(na))) {
        out <- a
        rep <- nb
    } else if (!is.null(sa)) {
        out <- b
        rep <- na
    } else {
        return(NULL)
    }
    out$pend <- if (length(rep)) lapply(rep, dm_recipe_replay)
    out
}

## The identity of a recipe across replays: its own id, or that of the
## recipe it copies; and those of the recipes bound in a frame's lineage.
dm_recipe_key <- function(r) r$origin %||% r$id
dm_frame_keys <- function(fr) vapply(fr$steps %||% list(), dm_recipe_key, "")

## ---- pending columns: products and lookups bound at first use ----
##
## The frame of a pending column is a settled frame (the design's, or the
## key's: the rows and columns its recipe reads) with a field `pend`, the
## recipes still to bind on top of it, in order. A recipe is an environment
## (kind "mult": design x and coefficients B, with the names of the product
## columns as colnames; kind "lookup": the SQL of its join), with `deps`
## (the recipes it reads, bound before it) and `frame`, NULL until bound,
## then the frame it made (later uses read that frame: one view per
## product). Every column made from a pending column copies its frame, so
## arithmetic stays pending; the operations that read rows or make frames
## call dm_frame_bind() (through dm_rebase() for designs).

dm_recipe <- function(kind, deps = list(), ...) {
    r <- list2env(list(...), parent = emptyenv())
    r$id <- dm_recipe_id()
    r$kind <- kind
    r$deps <- deps
    r$frame <- NULL
    r
}

dm_recipe_id <- function() {
    dm_tmp$recipes <- (dm_tmp$recipes %||% 0) + 1
    sprintf("r%x", as.integer(dm_tmp$recipes))
}

## A copy of the bound recipe r, not bound, to bind again on another branch
## (dm_frame_replay()): the same design and coefficients (or lookup table
## and join), so the same columns; `origin` names the recipe it copies. Its
## dependencies are left out: they are steps of the same branch, replayed
## before it (or below the branches' common ancestor).
dm_recipe_replay <- function(r) {
    q <- list2env(as.list(r, all.names = TRUE), parent = emptyenv())
    q$id <- dm_recipe_id()
    q$origin <- dm_recipe_key(r)
    q$deps <- list()
    q$frame <- NULL
    q
}

## the recipes of a and b, those of a first, each once (a replayed copy
## counts as its original)
dm_pend_union <- function(a, b) {
    if (!length(b)) return(a)
    if (!length(a)) return(b)
    keys <- vapply(a, dm_recipe_key, "")
    c(a, Filter(function(r) !dm_recipe_key(r) %in% keys, b))
}

## the recipes with their dependencies before them, each once
dm_pend_order <- function(pend) {
    out <- list()
    keys <- character()
    add <- function(r) {
        if (dm_recipe_key(r) %in% keys) return(invisible())
        for (d in r$deps) add(d)
        out[[length(out) + 1L]] <<- r
        keys <<- c(keys, dm_recipe_key(r))
        invisible()
    }
    for (r in pend) add(r)
    out
}

## Settled frame `fr` with the recipes `pend` on top: the recipes already
## bound (by an earlier use of another column) are folded into the settled
## frame (the later of it and their frames; for frames on two branches, one
## replayed on the other, whose copies of recipes join the pending ones, or
## an error); the others stay pending. (Copies made by a fold go before the
## recipes already pending: those may sit on the lineage they remake.)
dm_frame_pending <- function(fr, pend) {
    left <- list()
    for (r in pend) {
        if (is.null(r$frame)) {
            left <- dm_pend_union(left, list(r))
            next
        }
        fr <- dm_frame_later_settled(fr, r$frame)
        if (length(fr$pend)) {
            left <- dm_pend_union(fr$pend, left)
            fr$pend <- NULL
        }
    }
    fr$pend <- if (length(left)) left
    fr
}

## The settled frame of a (possibly pending) frame: its recipes bound one
## on the other, in order, on the latest settled frame among it and the
## frames of the recipes already bound. A recipe whose columns the frame
## already has (a replayed copy of a recipe in its lineage) is not bound
## again: it reads that frame.
dm_frame_bind <- function(fr) {
    if (is.null(fr) || !length(fr$pend)) return(fr)
    pend <- dm_pend_order(fr$pend)
    fr$pend <- NULL
    fr <- dm_frame_pending(dm_remote_unscaled(fr), pend)
    left <- fr$pend
    fr$pend <- NULL
    for (r in left) {
        if (dm_recipe_key(r) %in% dm_frame_keys(fr)) {
            r$frame <- fr
            next
        }
        fr <- switch(r$kind, mult = dm_bind_mult(r, fr), lookup = dm_bind_lookup(r, fr),
                     stop("internal error: unknown recipe", call. = FALSE))
        r$frame <- fr
    }
    fr
}

## A column on its settled frame (unchanged if it is not pending).
dm_column_bind <- function(x) {
    if (!dm_is_column(x)) return(x)
    u <- unclass(x)
    if (!length(u$frame$pend)) return(x)
    u$frame <- dm_frame_bind(u$frame)
    structure(u, class = "design_column")
}

## Does the column wait for recipes not yet bound?
dm_column_unbound <- function(x) {
    dm_is_column(x) && any(vapply(unclass(x)$frame$pend %||% list(), function(r) is.null(r$frame), NA))
}

## The product recipe r bound on settled frame fr: a view of X B over the
## rows of fr, which keeps every column of fr.
dm_bind_mult <- function(r, fr) {
    x <- r$x
    if (!identical(x$id, fr$id)) stop("internal error: a product on another design", call. = FALSE)
    x$src <- fr$src
    x$depth <- dm_frame_depth(fr)
    x$path <- fr$path
    x$steps <- fr$steps
    view <- dm_remote_predict(x, r$B, column = r$nm[1L])
    out <- dm_remote_unscaled(x)
    out$src$select <- sprintf("SELECT * FROM %s", dm_qi(x, unclass(view)[1L]))
    out$src$base_table <- NULL
    out$src$view <- c(out$src$view, list(view))
    ## (the view keeps every column of the frame, and adds the products:
    ## their type is known, no query)
    proto <- dm_proto_add(out, r$nm, double(0), sprintf("SELECT * FROM (%s) AS s LIMIT 0", out$src$select))
    out$src$proto <- proto
    out$src$types <- vapply(proto, dm_col_type, "")
    out$depth <- dm_frame_depth(x) + 1L
    dm_frame_child(out, x, r)
}

## The lookup recipe r bound on settled frame fr: the lookup table joined
## to the rows of fr (keeping every column of fr).
dm_bind_lookup <- function(r, fr) {
    out <- dm_remote_unscaled(fr)
    out$src$select <- sprintf("SELECT s.*, %s FROM (%s) AS s LEFT JOIN %s AS l ON %s",
                              r$items, fr$src$select, r$table, r$on)
    out$src$base_table <- NULL
    out$src$proto <- dm_proto_add(fr, r$col, r$type, sprintf("SELECT * FROM (%s) AS s LIMIT 0", out$src$select))
    out$src$types <- vapply(out$src$proto, dm_col_type, "")
    out$depth <- dm_frame_depth(fr) + 1L
    dm_frame_child(out, fr, r)
}

## The (possibly pending) frame of the remote columns in `cols`: a column,
## a dmc_codes() grouping, NULL, or a list of those (nested); NULL when
## none is remote.
dm_cols_frame <- function(cols) {
    f <- NULL
    walk <- function(v) {
        if (dm_is_column(v)) f <<- dm_frame_later(f, unclass(v)$frame)
        else if (dm_is_codes(v)) walk(v$codes)
        else if (is.list(v) && is.null(oldClass(v))) for (e in v) walk(e)
        invisible()
    }
    walk(cols)
    f
}

## `cols` (as dm_cols_frame() takes it) with fun() applied to each column.
dm_cols_map <- function(v, fun) {
    if (dm_is_column(v)) return(fun(v))
    if (dm_is_codes(v)) {
        v$codes <- dm_cols_map(v$codes, fun)
        return(v)
    }
    if (is.list(v) && is.null(oldClass(v))) return(lapply(v, dm_cols_map, fun))
    v
}

## Columns used one by one that must still combine with each other later
## (remotely, pending columns bound separately would sit on two branches):
## the pending columns of `cols` (as dm_cols_frame() takes it) bound
## together, one on the other in order, and every remote column of `cols`
## on that one frame (which keeps the columns of the frames below it). The
## same structure is returned; in memory, `cols` itself.
dmc_bind <- function(cols) {
    fr <- dm_cols_frame(cols)
    if (is.null(fr)) return(cols)
    fr <- dm_frame_bind(fr)
    dm_cols_map(cols, function(v) {
        u <- unclass(v)
        u$frame <- fr
        structure(u, class = "design_column")
    })
}

dm_column_frame <- function(...) {
    f <- NULL
    for (v in list(...)) if (dm_is_column(v)) f <- dm_frame_later(f, unclass(v)$frame)
    f
}

## SQL of an operand (a remote column, or a number)
dm_column_operand <- function(x) {
    if (dm_is_column(x)) return(unclass(x)$sql)
    ## (NA is NULL; NaN stays NaN, as in R arithmetic)
    if (length(x) == 1L && is.na(x) && !is.nan(x)) return("CAST(NULL AS DOUBLE)")
    if (is.logical(x) && length(x) == 1L) return(if (x) "TRUE" else "FALSE")
    if (is.numeric(x) && length(x) == 1L) return(dm_sql_double(x))
    stop("Remote columns combine only with numbers and other remote columns.", call. = FALSE)
}

#' @export
Ops.design_column <- function(e1, e2) {
    op <- .Generic
    if (nargs() == 1L) {
        a <- dm_column_operand(e1)
        sql <- switch(op, "-" = sprintf("(-%s)", a), "+" = a, "!" = sprintf("(NOT %s)", a),
                      stop(sprintf("`%s` is not supported for remote columns.", op), call. = FALSE))
        return(dm_column(unclass(e1)$frame, sql))
    }
    frame <- dm_column_frame(e1, e2)
    a <- dm_column_operand(e1)
    b <- dm_column_operand(e2)
    k <- max(length(a), length(b))
    if (!(length(a) %in% c(1L, k)) || !(length(b) %in% c(1L, k)))
        stop("Remote columns with different numbers of columns.", call. = FALSE)
    a <- rep_len(a, k)
    b <- rep_len(b, k)
    ## Comparisons with NaN are NA in R, but true or false in DuckDB (which
    ## orders NaN above every number and has NaN = NaN): compare with NaN
    ## made NULL (NULLIF(v, NaN), as NaN = NaN only for NaN, keeps v once in
    ## the expression), so that NaN or NULL on either side gives NULL.
    if (op %in% c("==", "!=", "<", ">", "<=", ">=")) {
        nn <- function(v, e) {
            if (dm_is_column(e)) return(sprintf("NULLIF(%s, CAST('NaN' AS DOUBLE))", v))
            if (is.nan(e)) "CAST(NULL AS DOUBLE)" else v
        }
        a <- nn(a, e1)
        b <- nn(b, e2)
    }
    sql <- switch(op,
        "+" = , "-" = , "*" = , "/" = sprintf("(%s %s %s)", a, op, b),
        "^" = dm_sql_pow(a, b, e2),
        "==" = sprintf("(%s = %s)", a, b), "!=" = sprintf("(%s <> %s)", a, b),
        "<" = , ">" = , "<=" = , ">=" = sprintf("(%s %s %s)", a, op, b),
        "&" = sprintf("(%s AND %s)", a, b), "|" = sprintf("(%s OR %s)", a, b),
        stop(sprintf("`%s` is not supported for remote columns.", op), call. = FALSE))
    dm_column(frame, sql)
}

## is.na() and friends, with R's semantics (NA is NULL remotely): is.na
## TRUE for NA or NaN, is.nan for NaN only, is.finite FALSE for NA, NaN and
## +-Inf, is.infinite TRUE for +-Inf only. Never NA.
#' @export
is.na.design_column <- function(x) !dmc_present(x)

#' @export
is.nan.design_column <- function(x)
    dm_column(unclass(x)$frame, sprintf("COALESCE(isnan(CAST(%s AS DOUBLE)), FALSE)", unclass(x)$sql),
             notnull = TRUE)

#' @export
is.finite.design_column <- function(x) dmc_is_finite(x)

#' @export
is.infinite.design_column <- function(x)
    dm_column(unclass(x)$frame, sprintf("COALESCE(isinf(CAST(%s AS DOUBLE)), FALSE)", unclass(x)$sql),
             notnull = TRUE)

#' @export
Math.design_column <- function(x, ...) {
    if (.Generic == "log") return(dm_column_log(x, ...))
    ## (log2 and log10 as R: -Inf at 0, NaN below, where DuckDB refuses)
    if (.Generic == "log2") return(dm_column_log(x, 2))
    if (.Generic == "log10") return(dm_column_log(x, 10))
    a <- unclass(x)$sql
    sql <- switch(.Generic,
        exp = sprintf("exp(%s)", a), sqrt = sprintf("sqrt(%s)", a),
        abs = sprintf("abs(%s)", a), log1p = dm_sql_log1p(sprintf("CAST(%s AS DOUBLE)", a)), expm1 = sprintf("(exp(%s) - 1)", a),
        cos = sprintf("cos(%s)", a), sin = sprintf("sin(%s)", a), tan = sprintf("tan(%s)", a),
        floor = sprintf("floor(%s)", a), ceiling = sprintf("ceil(%s)", a),
        ## (DuckDB's sign(NaN) is 0, R's NaN)
        sign = sprintf("CASE WHEN isnan(CAST(%s AS DOUBLE)) THEN CAST('NaN' AS DOUBLE) ELSE sign(%s) END", a, a),
        ## (round, signif, trunc, gamma, lgamma and the errors: R/columns_public.R)
        dm_math_sql(x, .Generic, ...))
    dm_column(unclass(x)$frame, sql)
}

## log(x) and log(x, base) of a remote column, as R computes them (nm8t):
## R_log(v) is ln(v) for v > 0, -Inf for 0 and NaN for v < 0 (DuckDB's ln()
## refuses 0 and negative numbers: it is reached only for v > 0, NaN or
## NULL, which it propagates); log(x, base) is log10(x) where base is 10,
## log2(x) where it is 2, else R_log(x) / R_log(base) (R's logbase()), so
## that base 1 gives x / 0 (Inf, -Inf or NaN), base 0 R_log(x) / -Inf and a
## negative base NaN; NA (NULL) in x or base gives NA, else NaN gives NaN.
## `base`: a number (its R_log() computed in R) or a column of the design
## (one column, or as many as x).
dm_column_log <- function(x, base) {
    xv <- sprintf("CAST(%s AS DOUBLE)", unclass(x)$sql)
    if (missing(base)) return(dm_column(unclass(x)$frame, dm_sql_rlog(xv)))
    if (dm_is_column(base)) {
        f <- dm_column_frame(x, base)
        bv <- sprintf("CAST(%s AS DOUBLE)", unclass(base)$sql)
        k <- max(length(xv), length(bv))
        if (!(length(xv) %in% c(1L, k)) || !(length(bv) %in% c(1L, k)))
            stop("log() of remote columns with different numbers of columns.", call. = FALSE)
        xv <- rep_len(xv, k)
        bv <- rep_len(bv, k)
        sql <- sprintf("(CASE WHEN %1$s = 10 THEN %2$s WHEN %1$s = 2 THEN %3$s ELSE %4$s / %5$s END)",
                       bv, dm_sql_rlog(xv, "log10"), dm_sql_rlog(xv, "log2"), dm_sql_rlog(xv), dm_sql_rlog(bv))
        return(dm_column(f, sql))
    }
    if (!(is.numeric(base) || is.logical(base)) || length(base) != 1L || !is.null(dim(base)))
        stop("log() of a remote column: `base` must be one number or a column of the same design.", call. = FALSE)
    dm_column(unclass(x)$frame, dm_sql_log(xv, base))
}

## R_log(v) in SQL, `v` a DOUBLE operand (named three times: a short one);
## `fun` ln, log2 or log10.
dm_sql_rlog <- function(v, fun = "ln")
    sprintf("(CASE WHEN %1$s = 0 THEN CAST('-Infinity' AS DOUBLE) WHEN %1$s < 0 THEN CAST('NaN' AS DOUBLE) ELSE %2$s(%1$s) END)",
            v, fun)

## log(v, base) in SQL as R computes it (see dm_column_log()), `v` a DOUBLE
## operand, `base` one number or logical.
dm_sql_log <- function(v, base) {
    base <- as.double(base)
    if (!is.na(base) && base == 10) return(dm_sql_rlog(v, "log10"))
    if (!is.na(base) && base == 2) return(dm_sql_rlog(v, "log2"))
    ## R_log(base), in R (NA is NULL: the result is NA)
    den <- if (is.na(base) && !is.nan(base)) "CAST(NULL AS DOUBLE)"
           else dm_sql_double(if (is.nan(base)) NaN else if (base > 0) log(base) else if (base == 0) -Inf else NaN)
    sprintf("(%s / %s)", dm_sql_rlog(v), den)
}

## log1p(v) in SQL, `v` a DOUBLE operand (named several times: a short
## one). DuckDB has no log1p(), and ln(1 + v) loses the digits of v near 0:
## with u = 1 + v rounded, ln(u) * v / (u - 1) is accurate to a few ulps
## (Goldberg 1991, Theorem 4), and v itself where u rounds to 1. As R: -Inf
## at -1, NaN below (and at -Inf), Inf at Inf, NaN for NaN, NA for NA.
dm_sql_log1p <- function(v)
    sprintf(paste("(CASE WHEN %1$s = -1 THEN CAST('-Infinity' AS DOUBLE)",
                  "WHEN %1$s < -1 THEN CAST('NaN' AS DOUBLE)",
                  "WHEN 1 + %1$s = 1 OR isinf(%1$s) THEN %1$s",
                  "ELSE ln(1 + %1$s) * %1$s / ((1 + %1$s) - 1) END)"), v)

## x^y in SQL (a, b the SQL of the operands, e2 the R operand y), with R's
## R_pow() where DuckDB's pow() differs: x^0 and 1^y are 1 for any x and y
## (NA and NaN included; pow(NULL, 0) is NULL), and a negative base with
## an infinite power, or -Inf with a power that is not an integer, gives
## NaN. A power that is a nonzero integer (as y^2) needs none of this.
dm_sql_pow <- function(a, b, e2) {
    if (!dm_is_column(e2) && is.numeric(e2) && length(e2) == 1L && is.finite(e2) && e2 != 0 && e2 == round(e2))
        return(sprintf("pow(%s, %s)", a, b))
    sprintf(paste("(CASE WHEN %2$s = 0 OR %1$s = 1 THEN CAST(1 AS DOUBLE)",
                  "WHEN %1$s < 0 AND (isinf(CAST(%2$s AS DOUBLE)) OR (isinf(CAST(%1$s AS DOUBLE)) AND %2$s <> floor(%2$s)))",
                  "THEN CAST('NaN' AS DOUBLE) ELSE pow(%1$s, %2$s) END)"), a, b)
}

## pmax/pmin with R's missing values: DuckDB's greatest() and least() skip
## NULL and order NaN above every number, while in R an NA or NaN gives NA or
## NaN (that of y where both are missing, as R's pmax takes the later
## missing argument). NULLIF(v, NaN) is NULL where v is NA or NaN. With a
## finite number on one side, the column is kept where it is missing or
## beyond the number. (The fits use dm_fast_pmax/dm_fast_pmin, below.)
dmc_pmax <- function(x, y) {
    f <- dm_column_frame(x, y)
    if (is.null(f)) return(pmax(x, y))
    dm_sql_pext(f, x, y, "greatest", "<=")
}

dmc_pmin <- function(x, y) {
    f <- dm_column_frame(x, y)
    if (is.null(f)) return(pmin(x, y))
    dm_sql_pext(f, x, y, "least", ">=")
}

dm_sql_pext <- function(f, x, y, fun, cmp) {
    num <- function(v) !dm_is_column(v) && is.numeric(v) && length(v) == 1L && is.finite(v)
    a <- dm_column_operand(x)
    b <- dm_column_operand(y)
    nan <- "CAST('NaN' AS DOUBLE)"
    if (num(x) || num(y)) {
        if (num(x)) { cc <- a; v <- b } else { cc <- b; v <- a }
        return(dm_column(f, sprintf("CASE WHEN NULLIF(%1$s, %3$s) %4$s %2$s THEN %2$s ELSE %1$s END",
                                   v, cc, nan, cmp)))
    }
    dm_column(f, sprintf(paste("CASE WHEN NULLIF(%2$s, %3$s) IS NULL THEN %2$s",
                              "WHEN NULLIF(%1$s, %3$s) IS NULL THEN %1$s ELSE %4$s(%1$s, %2$s) END"),
                        a, b, nan, fun))
}

## Internal pmax/pmin for the fits (qpr8), only where the inputs are never
## missing: remotely DuckDB's greatest() and least(), which name each
## operand once. The CASE of dmc_pmax/dmc_pmin repeats its operand, and the
## links nest them around eta inside the pnorm/dnorm macros, which repeat
## their argument again (a small probit fit: 1.3 s -> 3.4 s). They give
## pmax/pmin wherever no input is NA or NaN (greatest() skips NULL; least()
## skips NULL and NaN); each use states why its inputs are not. In memory:
## pmax/pmin.
dm_fast_pmax <- function(x, y) {
    f <- dm_column_frame(x, y)
    if (is.null(f)) return(pmax(x, y))
    dm_column(f, sprintf("greatest(%s, %s)", dm_column_operand(x), dm_column_operand(y)))
}

dm_fast_pmin <- function(x, y) {
    f <- dm_column_frame(x, y)
    if (is.null(f)) return(pmin(x, y))
    dm_column(f, sprintf("least(%s, %s)", dm_column_operand(x), dm_column_operand(y)))
}

## Internal log for the fits (qpr8), only where x is positive by
## construction (or NULL): remotely DuckDB's ln(), which names x once (the
## CASE of log(), R's log of 0 and of negative numbers, names it three
## times; at 0 and below, ln() gives -Inf and NaN in DuckDB 2.0 but older
## versions refuse them). In memory: log().
dm_fast_log <- function(x) {
    if (!dm_is_column(x)) return(log(x))
    dm_column(unclass(x)$frame, sprintf("ln(CAST(%s AS DOUBLE))", unclass(x)$sql))
}

## (in memory a numeric vector: `no` with `yes` where test is TRUE and NA
## where it is NA, as ifelse() but without its copies and attributes)
dmc_ifelse <- function(test, yes, no) {
    f <- dm_column_frame(test, yes, no)
    if (is.null(f)) {
        n <- length(test)
        ## (one pass in C, which allocates only the result)
        if (is.logical(test) && n > 0L && length(yes) %in% c(1L, n) && length(no) >= 1L &&
            is.numeric(yes) && is.numeric(no)) {
            if (!is.double(yes)) yes <- as.double(yes)
            if (!is.double(no)) no <- as.double(no)
            return(.Call(C_dm_ifelse, test, yes, no))
        }
        out <- rep_len(as.double(no), n)
        i <- which(test)
        out[i] <- if (length(yes) == 1L) yes else yes[i]
        if (anyNA(test)) out[is.na(test)] <- NA
        return(out)
    }
    ## (NULL where test is NULL, as NA in memory; the simple CASE names the
    ## test once)
    dm_column(f, sprintf("(CASE %s WHEN TRUE THEN %s WHEN FALSE THEN %s END)", dm_column_operand(test),
                        dm_column_operand(yes), dm_column_operand(no)))
}

## Columns side by side (cbind).
dmc_cbind <- function(...) {
    args <- list(...)
    f <- do.call(dm_column_frame, args)
    if (is.null(f)) return(do.call(cbind, args))
    dm_column(f, unlist(lapply(args, dm_column_operand)),
             notnull = unlist(lapply(args, function(a) if (dm_is_column(a)) unclass(a)$notnull else !is.na(a))))
}

dmc_rowsums <- function(x) {
    if (!dm_is_column(x)) return(if (is.matrix(x)) rowSums(x) else x)
    dm_column(unclass(x)$frame, sprintf("(%s)", paste(unclass(x)$sql, collapse = " + ")))
}

## distribution functions: stats in memory, DuckDB macros remotely
dmc_macro <- function(x, fun, macro) {
    if (!dm_is_column(x)) return(fun(x))
    dm_sql_macros(unclass(x)$frame)
    dm_column(unclass(x)$frame, sprintf("%s(%s)", macro, unclass(x)$sql))
}
dmc_pnorm <- function(x) dmc_macro(x, stats::pnorm, "dm_pnorm")
dmc_qnorm <- function(x) dmc_macro(x, stats::qnorm, "dm_qnorm")
dmc_dnorm <- function(x) dmc_macro(x, stats::dnorm, "dm_dnorm")
dmc_pcauchy <- function(x) dmc_macro(x, stats::pcauchy, "dm_pcauchy")
dmc_mexpm1m <- function(x) dmc_macro(x, function(u) -expm1(-u), "dm_mexpm1m")
## The logit link of binomial(): inverse, derivative, and the binomial
## deviance residuals (in memory the C code of R's family object, in one
## pass; remotely the same arithmetic in SQL)
dmc_logit_inv <- function(e) {
    if (!dm_is_column(e)) return(stats::binomial()$linkinv(e))
    eps <- .Machine$double.eps
    dmc_ifelse(e < -30, eps / (1 + eps), dmc_ifelse(e > 30, (1 / eps) / (1 + 1 / eps), exp(e) / (1 + exp(e))))
}
dmc_logit_mu_eta <- function(e) {
    if (!dm_is_column(e)) return(stats::binomial()$mu.eta(e))
    dmc_ifelse(e > 30 | e < -30, .Machine$double.eps, exp(e) / ((1 + exp(e)) * (1 + exp(e))))
}
## (`positive`: mu is in (0, 1) by construction, so that the logarithms,
## reached only for y in (0, 1], take the fast ln of dm_fast_log())
dmc_binomial_dev <- function(y, mu, w, positive = FALSE) {
    if (is.null(dm_column_frame(y, mu, w))) return(stats::binomial()$dev.resids(y, mu, w))
    lg <- if (positive) dm_fast_log else log
    xlogx <- function(a, b) dmc_ifelse(a != 0, a * lg(a / b), 0)
    2 * w * (xlogx(y, mu) + xlogx(1 - y, 1 - mu))
}

## qcauchy (remotely through tan(): to rounding, where R uses tanpi())
dmc_qcauchy <- function(x) {
    if (!dm_is_column(x)) return(stats::qcauchy(x))
    dmc_ifelse(x <= 0.5, -1 / tan(pi * x), 1 / tan(pi * (1 - x)))
}

## The values of a column in R (remotely: all its rows; for small data).
dmp_collect <- function(x) {
    if (!dm_is_column(x)) return(x)
    fr <- dm_frame_bind(unclass(x)$frame)
    sql <- unclass(x)$sql
    res <- dm_query(fr, sprintf("SELECT %s FROM (%s) AS s WHERE %s",
                                paste(sprintf("CAST(%s AS DOUBLE) AS c%d", sql, seq_along(sql)), collapse = ", "),
                                fr$src$select, fr$filter))
    if (length(sql) == 1L) res[[1L]] else as.matrix(res)
}

## ---- sums over rows ----

dm_colsum <- function(x) {
    if (!dm_is_column(x)) {
        if (is.matrix(x)) return(colSums(x))
        return(sum(x))
    }
    fr <- dm_frame_bind(unclass(x)$frame)
    sql <- unclass(x)$sql
    k <- length(sql)
    nn <- unclass(x)$notnull %||% rep(FALSE, k)
    items <- c(sprintf("SUM(CAST(%s AS DOUBLE)) AS s%d", sql, seq_len(k)),
               sprintf("SUM(CASE WHEN (%s) IS NULL THEN 1 ELSE 0 END) AS n%d", sql, seq_len(k))[!nn])
    res <- dm_query(fr, sprintf("SELECT %s FROM (%s) AS s WHERE %s", paste(items, collapse = ", "),
                                fr$src$select, fr$filter))
    vapply(seq_len(k), function(j) {
        v <- as.double(res[[paste0("s", j)]])
        if (!nn[j] && as.numeric(res[[paste0("n", j)]]) > 0) return(NA_real_)
        if (is.na(v)) 0 else v
    }, 0)
}

## Sums over rows of several columns (vectors or numbers in memory; one
## query remotely): a number per argument.
dmc_sums <- function(...) {
    a <- list(...)
    if (any(vapply(a, dm_is_column, NA))) return(dm_colsum(do.call(dmc_cbind, a)))
    vapply(a, function(v) sum(v), 0)
}

## The mean of a column (in memory as mean() computes it).
dmc_mean <- function(x) if (dm_is_column(x)) dm_colsum(x) / dm_nrows(x) else mean(x)

## "sql" for a remote design, "" in memory.
dm_engine <- function(x) if (inherits(x, "dm_remote")) "sql" else ""

dm_nrows <- function(x) {
    if (inherits(x, "dm_remote")) return(x$cache$n %||% x$n)
    if (dm_is_column(x)) return(dm_nrows(unclass(x)$frame))
    NROW(x)
}

dmp_colnorms <- function(x) {
    if (!dm_is_column(x)) return(.Call(C_dm_col_norms, if (is.matrix(x)) x else matrix(x, ncol = 1L)))
    sqrt(dm_colsum(x * x))
}

## ---- designs ----

dmp_column <- function(x, sql) {
    if (is.null(sql)) return(NULL)
    if (inherits(x, "dm_remote")) return(dm_column(x, sql))
    sql
}

## The remote design over the frame of columns (same plan, rows and scale);
## a pending frame is bound first.
dm_rebase <- function(x, frame) {
    frame <- dm_frame_bind(frame)
    if (is.null(frame) || identical(x$src$select, frame$src$select)) return(x)
    if (!identical(x$id, frame$id)) stop("Columns of another remote design.", call. = FALSE)
    x$src <- frame$src
    x$depth <- dm_frame_depth(frame)
    x$path <- frame$path
    x$steps <- frame$steps
    x
}

## The remote design over the later of its frame and `frame` (the pending
## columns of `frame` bound on top of the design's frame when it is the
## later one). (A frame of another design: as dm_rebase().)
dm_rebase_later <- function(x, frame) {
    if (!is.null(frame) && identical(frame$id, x$id)) frame <- dm_frame_later(frame, x)
    dm_rebase(x, frame)
}

## The design over the latest frame of the columns of `cols` (a list, as
## dmc_bind() takes it): their pending columns bound together, on top of
## the design's frame when it is the later one. In memory x.
dmp_rebase <- function(x, cols) {
    if (!inherits(x, "dm_remote")) return(x)
    dm_rebase_later(x, dm_cols_frame(cols))
}

dmp_scale <- function(x, s) {
    if (is.null(s)) return(x)
    if (!dm_is_column(s)) return(s * x)
    x <- dm_rebase_later(x, unclass(s)$frame)
    x$scale <- dm_remote_scale_sql(x, DBI::SQL(unclass(s)$sql), "*")
    x
}

## A remote design with the columns `col` added to its source (named nm).
dm_with_columns <- function(x, col, nm) {
    x <- dm_rebase(x, dm_frame_later(unclass(col)$frame, x))
    sc <- x$scale
    x <- dm_frame_child(dm_remote_extend(dm_remote_unscaled(x), stats::setNames(unclass(col)$sql, nm), "column"), x)
    x$scale <- sc
    x
}

dmp_crossprod <- function(x, y = NULL, compress = TRUE) {
    ## A refinement step reads the residual column once. Skip the pattern
    ## estimate and aggregate its products directly.
    if (!compress && inherits(x, "dm_remote") && dm_is_column(y)) {
        yc <- dm_remote_cols(x, y, "y")
        return(dm_remote_moments(dm_remote_unscaled(yc$D), yc$sql, yc$D$scale,
                                 gram = FALSE, compress = FALSE)$xty)
    }
    if (is.null(y)) return(crossprod(x))
    if (inherits(x, "dm_remote") && !inherits(y, "dm_remote") && !dm_is_column(y))
        stop("`y` must be a design or a column.", call. = FALSE)
    crossprod(x, y)
}

dmp_colsums <- function(x, group = NULL) {
    if (is.null(group)) colSums(x) else colSums(x, group = group)
}

## SQL GROUP BY over the frame of columns: list(groups (a data frame of the
## group values, sorted as rowsum() sorts them), sums (matrix)).
dm_group_query <- function(groups, values, nulls = FALSE) {
    fr <- dm_frame_bind(do.call(dm_column_frame, c(groups, values)))
    gs <- vapply(groups, function(g) unclass(g)$sql, "")
    vs <- unlist(lapply(values, function(v) unclass(v)$sql))
    items <- c(sprintf("%s AS g%d", gs, seq_along(gs)),
               if (length(vs)) sprintf("SUM(CAST(%s AS DOUBLE)) AS v%d", vs, seq_along(vs)) else "COUNT(*) AS v1",
               if (nulls && length(vs)) sprintf("SUM(CASE WHEN (%s) IS NULL THEN 1 ELSE 0 END) AS n%d", vs, seq_along(vs)))
    res <- dm_query(fr, sprintf("SELECT %s FROM (%s) AS s WHERE %s GROUP BY ALL",
                                paste(items, collapse = ", "), fr$src$select, fr$filter))
    res
}

dmp_groupsum <- function(col, group) {
    if (dm_is_codes(group)) return(dmp_groupsum_codes(col, group))
    if (!dm_is_column(col) && !dm_is_column(group)) {
        return(rowsum(if (is.matrix(col)) col else matrix(col, ncol = 1L), group, reorder = TRUE))
    }
    res <- dm_group_query(list(group), list(col), nulls = TRUE)
    ug <- sort(unique(res$g1), na.last = TRUE, method = "quick")
    if (anyNA(ug)) warning("missing values for 'group'", call. = FALSE)
    i <- match(res$g1, ug)
    k <- dm_ncols(col)
    out <- matrix(0, length(ug), k, dimnames = list(as.character(ug), NULL))
    for (j in seq_len(k)) {
        v <- as.double(res[[paste0("v", j)]])
        ## a group with a missing value sums to NA, as rowsum() gives
        v[as.numeric(res[[paste0("n", j)]]) > 0] <- NA_real_
        out[i, j] <- v
    }
    out
}

dmp_crosstab <- function(a, b, w = NULL) {
    if (dm_is_codes(a) && dm_is_codes(b)) return(dmp_crosstab_codes(a, b, w))
    if (!dm_is_column(a) && !dm_is_column(b) && !dm_is_column(w)) {
        ua <- sort(unique(a), na.last = TRUE, method = "quick")
        ub <- sort(unique(b), na.last = TRUE, method = "quick")
        return(Matrix::sparseMatrix(i = match(a, ua), j = match(b, ub),
                                    x = if (is.null(w)) rep(1, length(a)) else as.double(w),
                                    dims = c(length(ua), length(ub)),
                                    dimnames = list(as.character(ua), as.character(ub))))
    }
    res <- dm_group_query(list(a, b), if (is.null(w)) list() else list(w))
    ## the group values as in memory: sorted in R, labels as.character()
    ua <- sort(unique(res$g1), na.last = TRUE, method = "quick")
    ub <- sort(unique(res$g2), na.last = TRUE, method = "quick")
    Matrix::sparseMatrix(i = match(res$g1, ua), j = match(res$g2, ub), x = as.double(res$v1),
                         dims = c(length(ua), length(ub)),
                         dimnames = list(as.character(ua), as.character(ub)))
}

## (remotely a pending column: the view of the product is made at its first
## use, on the frame that use needs; B is checked now)
dmp_mult <- function(x, B) {
    if (!inherits(x, "dm_remote")) return(if (is.matrix(B)) x %*% B else drop(x %*% B))
    B <- if (is.matrix(B)) B else matrix(B, ncol = 1L)
    checkmate::assert_matrix(B, mode = "numeric", nrows = x$p, any.missing = FALSE, .var.name = "beta")
    checkmate::assert_true(all(is.finite(B)), .var.name = "all(is.finite(beta))")
    storage.mode(B) <- "double"
    nm <- paste0("__dm_p", substr(dm_tmp_name("p"), 3L, 10L), "_", seq_len(ncol(B)))
    colnames(B) <- nm
    fr <- dm_remote_unscaled(x)
    fr$pend <- list(dm_recipe("mult", x = x, B = B, nm = nm))
    dm_column(fr, sprintf("s.%s", dm_qi(x, nm)))
}

dmp_quad <- function(x, V) {
    if (!inherits(x, "dm_remote")) {
        J <- as_jac(x)
        storage.mode(V) <- "double"
        if (!jac_kernel_ok(J)) {
            X <- as.matrix(J)
            return(unname(rowSums((X %*% V) * X)))
        }
        g <- dm_groups(J, cost = dm_row_cost(J)^2)
        q <- do.call(.Call, c(list(C_dm_jac_quad), jac_args(J),
                              list(V, if (is.null(g)) NULL else g$first, design_threads())))
        return(if (is.null(g)) q else q[g$map])
    }
    zero <- dm_column(x, "CAST(0 AS DOUBLE)")
    if (ncol(V) == 0L) return(zero)
    e <- eigen(V, symmetric = TRUE)
    keep <- e$values > 0
    if (!any(keep)) return(zero)
    L <- e$vectors[, keep, drop = FALSE] %*% diag(sqrt(e$values[keep]), sum(keep))
    P <- dmp_mult(x, L)
    dmc_rowsums(P * P)
}

## Compressed least squares (see the header); the in-memory method groups
## identical rows (lm.R), the remote one is the scan of the SQL engine
## grouped by the cells of the design (remote_compress.R).
dmp_patterns <- function(x, z, w = NULL, compress = "auto", keys = NULL, fit = "lm") {
    if (inherits(x, "dm_remote")) {
        mode <- if (isTRUE(compress)) TRUE else if (isFALSE(compress)) FALSE else dm_compress_mode()
        ## (on the frame of z and w: their pending columns bound)
        x <- dm_rebase_later(x, dm_column_frame(z, w))
        ## (extra key columns: more keys of the GROUP BY)
        if (!is.null(keys)) x <- dm_remote_cols(x, keys, "keys")$D
        pt <- dm_sql_patterns(x, unclass(z)$sql, if (!is.null(w)) unclass(w)$sql,
                              keys = if (!is.null(keys)) unclass(keys)$sql, mode = mode, fit = fit)
        if (is.null(pt)) return(NULL)
        pos <- which(pt$sw > 0)
        ws <- pt$sw[pos]
        if (is.null(w) && all(ws == 1)) ws <- NULL
        return(list(x = dm_dense_design(pt$x[pos, , drop = FALSE]), z = pt$sy[pos] / pt$sw[pos], w = ws,
                    within = sum(pt$within), rows = sum(pt$n),
                    keys = if (!is.null(keys)) pt$keys[pos, , drop = FALSE]))
    }
    g <- dmfit_groups(x, compress)
    if (is.null(g)) return(NULL)
    if (!is.null(keys)) {
        ## patterns within the values of the key columns
        K <- if (is.matrix(keys)) keys else matrix(keys, ncol = 1L)
        key <- do.call(paste, c(list(g$map), lapply(seq_len(ncol(K)), function(j) sprintf("%a", K[, j])), sep = "\r"))
        map <- match(key, unique(key))
        g <- list(map = map, first = match(seq_len(max(map)), map))
        if (length(g$first) == nrow(x) || (!isTRUE(compress) && length(g$first) > 0.5 * nrow(x))) return(NULL)
    }
    ss <- dmfit_suffstats(z, w, g$map, length(g$first))
    ## sums that overflow (w * z beyond the range of doubles): no compression
    if (!all(is.finite(ss$w), is.finite(ss$mean), is.finite(ss$within))) return(NULL)
    keep <- ss$w > 0
    ws <- ss$w[keep]
    if (is.null(w) && all(ws == 1)) ws <- NULL
    list(x = x[g$first[keep], ], z = ss$mean[keep], w = ws, within = ss$within, rows = nrow(x),
         keys = if (!is.null(keys)) (if (is.matrix(keys)) keys else matrix(keys, ncol = 1L))[g$first[keep], , drop = FALSE])
}

## Row-wise results: in memory memory(d), d the data frame of the columns
## (default: d); remotely a table (a temporary view with the columns of the
## remote data and `cols`), padded by napredict(na) for a fit with
## na.exclude, with attributes `attrs`; with design = TRUE, the remote
## design with the columns added to its data (so that they can scale its
## rows).
dmp_output <- function(x, cols, memory = identity, na = NULL, design = FALSE, attrs = list()) {
    if (!inherits(x, "dm_remote")) return(memory(as.data.frame(lapply(cols, function(v) as.vector(v)))))
    taken <- intersect(names(cols), names(x$src$types))
    if (length(taken)) stop(sprintf("The remote data already has a column `%s`.", taken[1L]), call. = FALSE)
    if (design) return(dm_with_columns(x, do.call(dmc_cbind, unname(cols)), names(cols)))
    fr <- dm_rebase_later(dm_remote_unscaled(x), do.call(dm_column_frame, unname(cols)))
    srccols <- names(x$src$types)
    items <- vapply(names(cols), function(nm) {
        v <- cols[[nm]]
        sprintf("%s AS %s", if (dm_is_column(v)) unclass(v)$sql else dm_column_operand(v), dm_qi(fr, nm))
    }, "")
    sel <- sprintf("SELECT %s FROM (%s) AS s WHERE %s",
                   paste(c(paste0("s.", dm_qi(fr, srccols)), items), collapse = ", "), fr$src$select, fr$filter)
    name <- dm_tmp_name("dm_out")
    dm_execute(fr, sprintf("CREATE TEMPORARY VIEW %s AS %s", dm_qi(fr, name), sel))
    out <- stats::napredict(na, structure(name, class = "dm_remote_table", connection = dm_con(fr),
                                          holder = fr$src$holder))
    for (a in names(attrs)) attr(out, a) <- attrs[[a]]
    out
}

## DuckDB macros for the distribution functions of remote columns, as R
## computes them (from the GLM work, kata/402b-glm): dm_pnorm, dm_qnorm
## (one ulp), dm_dnorm, dm_pcauchy and dm_mexpm1m(u) = -expm1(-u). Created
## once per connection.
dm_sql_macros <- function(D) {
    if (isTRUE(D$src$holder$dm_macros)) return(invisible())
    horner <- function(x, co) {
        s <- sprintf("%.17g", co[1L])
        for (k in co[-1L]) s <- sprintf("(%s * %s + %.17g)", s, x, k)
        s
    }
    a <- c(2.2352520354606839287, 161.02823106855587881, 1067.6894854603709582, 18154.981253343561249,
           0.065682337918207449113)
    b <- c(47.20258190468824187, 976.09855173777669322, 10260.932208618978205, 45507.789335026729956)
    cc <- c(0.39894151208813466764, 8.8831497943883759412, 93.506656132177855979, 597.27027639480026226,
            2494.5375852903726711, 6848.1904505362823326, 11602.651437647350124, 9842.7148383839780218,
            1.0765576773720192317e-8)
    d <- c(22.266688044328115691, 235.38790178262499861, 1519.377599407554805, 6485.558298266760755,
           18615.571640885098091, 34900.952721145977266, 38912.003286093271411, 19685.429676859990727)
    p <- c(0.21589853405795699, 0.1274011611602473639, 0.022235277870649807, 0.001421619193227893466,
           2.9112874951168792e-5, 0.02307344176494017303)
    q <- c(1.28426009614491121, 0.468238212480865118, 0.0659881378689285515, 0.00378239633202758244,
           7.29751555083966205e-5)
    num1 <- sprintf("(x * %s)", horner("x*x", c(a[5L], a[1:4])))
    den1 <- horner("x*x", c(1, b[1:4]))
    num2 <- horner("abs(x)", c(cc[9L], cc[1:8]))
    den2 <- horner("abs(x)", c(1, d[1:8]))
    ix <- "(1.0/(x*x))"
    num3 <- sprintf("(%s * %s)", ix, horner(ix, c(p[6L], p[1:5])))
    den3 <- horner(ix, c(1, q[1:5]))
    xs <- "(trunc(abs(x)*16)/16)"
    tail <- function(temp) sprintf("(exp(-%s*%s*0.5) * exp(-((abs(x) - %s)*(abs(x) + %s))*0.5) * %s)", xs, xs, xs, xs, temp)
    t2 <- tail(sprintf("(%s / %s)", num2, den2))
    t3 <- tail(sprintf("((0.398942280401432677939946059934 - %s / %s) / abs(x))", num3, den3))
    pn <- sprintf(paste("CASE WHEN abs(x) <= 0.67448975 THEN 0.5 + %s / %s",
                        "WHEN abs(x) <= 5.656854249492380195206754896838 THEN (CASE WHEN x > 0 THEN 1 - %s ELSE %s END)",
                        "WHEN abs(x) < 50 THEN (CASE WHEN x > 0 THEN 1 - %s ELSE %s END)",
                        "ELSE (CASE WHEN x > 0 THEN 1.0 ELSE 0.0 END) END"), num1, den1, t2, t2, t3, t3)
    qa <- c(2509.0809287301226727, 33430.575583588128105, 67265.770927008700853, 45921.953931549871457,
            13731.693765509461125, 1971.5909503065514427, 133.14166789178437745, 3.387132872796366608)
    qb <- c(5226.495278852545925, 28729.085735721942674, 39307.89580009271061, 21213.794301586595867,
            5394.1960214247511077, 687.1870074920579083, 42.313330701600911252, 1)
    qc <- c(7.7454501427834140764e-4, .0227238449892691845833, .24178072517745061177, 1.27045825245236838258,
            3.64784832476320460504, 5.7694972214606914055, 4.6303378461565452959, 1.42343711074968357734)
    qd <- c(1.05075007164441684324e-9, 5.475938084995344946e-4, .0151986665636164571966, .14810397642748007459,
            .68976733498510000455, 1.6763848301838038494, 2.05319162663775882187, 1)
    qe <- c(2.01033439929228813265e-7, 2.71155556874348757815e-5, .0012426609473880784386, .026532189526576123093,
            .29656057182850489123, 1.7848265399172913358, 5.4637849111641143699, 6.6579046435011037772)
    qf <- c(2.04426310338993978564e-15, 1.4215117583164458887e-7, 1.8463183175100546818e-5, 7.868691311456132591e-4,
            .0148753612908506148525, .13692988092273580531, .59983220655588793769, 1)
    r1 <- "(0.180625 - (p - 0.5)*(p - 0.5))"
    rr <- "sqrt(-ln(least(p, 1 - p)))"
    r2 <- sprintf("(%s - 1.6)", rr)
    r3 <- sprintf("(%s - 5)", rr)
    qn <- sprintf(paste("CASE WHEN abs(p - 0.5) <= 0.425 THEN (p - 0.5) * %s / %s",
                        "ELSE (CASE WHEN p < 0.5 THEN -1 ELSE 1 END) * (CASE WHEN %s <= 5 THEN %s / %s ELSE %s / %s END) END"),
                  horner(r1, qa), horner(r1, qb), rr, horner(r2, qc), horner(r2, qd), horner(r3, qe), horner(r3, qf))
    ## dnorm as dnorm.c: for |x| >= 5, exp() of a split square
    x1 <- "(round(abs(x) * 65536) / 65536)"
    dn <- sprintf(paste("CASE WHEN abs(x) < 5 THEN 0.398942280401432677939946059934 * exp(-0.5 * x * x)",
                        "ELSE 0.398942280401432677939946059934 * exp(-0.5 * %s * %s) * exp((-0.5 * (abs(x) - %s) - %s) * (abs(x) - %s)) END"),
                  x1, x1, x1, x1, x1)
    pc <- "CASE WHEN abs(x) > 1 THEN (CASE WHEN x > 0 THEN 1 - atan(1 / x) / pi() ELSE -atan(1 / x) / pi() END) ELSE 0.5 + atan(x) / pi() END"
    ## -expm1(-u) for u >= 0, Kahan's trick: (1 - w) u / -ln(w), w = exp(-u)
    em <- "CASE WHEN u > 40 THEN 1.0 WHEN exp(-u) = 1 THEN u ELSE (1 - exp(-u)) * u / (-ln(exp(-u))) END"
    for (m in list(c("dm_pnorm", "x", pn), c("dm_qnorm", "p", qn), c("dm_dnorm", "x", dn),
                   c("dm_pcauchy", "x", pc), c("dm_mexpm1m", "u", em)))
        dm_execute(D, sprintf("CREATE OR REPLACE TEMP MACRO %s(%s) AS (%s)", m[1L], m[2L], m[3L]))
    D$src$holder$dm_macros <- TRUE
    invisible()
}


## ---- coded groups, multi-group sums (the kernels of the fixed effects) ----

## In memory, codes are integers 1..G. A remote FE grouping carries the
## source key and its level table; public coded groupings use integer codes.
## Grouped sums with coded groups return all G groups in code order
## (zeros for empty groups), without labels, and use one-pass C kernels in
## memory (src/fe.c).
dmc_codes <- function(codes, G) {
    if (!dm_is_column(codes)) codes <- as.integer(codes)
    structure(list(codes = codes, G = as.integer(G)), class = "dm_codes")
}

## Fixed effects can group directly on source keys. Translate the returned
## group labels to solver indices only after DuckDB has aggregated the rows.
dm_group_index <- function(g, values) {
    if (is.null(g$levels)) as.integer(values) else match(values, g$levels[[1L]])
}

dm_is_codes <- function(x) inherits(x, "dm_codes")

## sums of columns within coded groups: G x k
dmp_groupsum_codes <- function(col, group) {
    G <- group$G
    if (!dm_is_column(col) && !dm_is_column(group$codes)) {
        ## (the kernel takes a vector as one column: no copy of a double)
        W <- col
        if (!is.double(W)) storage.mode(W) <- "double"
        return(.Call(C_dm_fe_gsum, group$codes, G, W, design_threads()))
    }
    res <- dm_group_query(list(group$codes), list(col), nulls = TRUE)
    k <- dm_ncols(col)
    out <- matrix(0, G, k)
    i <- dm_group_index(group, res$g1)
    for (j in seq_len(k)) {
        v <- as.double(res[[paste0("v", j)]])
        v[as.numeric(res[[paste0("n", j)]]) > 0] <- NA_real_
        out[i, j] <- v
    }
    out
}

## Sums of a column (k columns) within each of several coded groupings, and
## of another column `total` (or NULL) over all rows: list(groups = one G x k
## matrix per grouping, as dmp_groupsum(); total = dm_colsum(total)). In
## memory those calls; remotely one query (GROUP BY GROUPING SETS), whose
## rows the frame computes once (the frame of residuals, with its joins).
dmc_groupsums <- function(col, groups, total = NULL) {
    remote <- dm_is_column(col) && all(vapply(groups, function(g) dm_is_column(g$codes), NA))
    if (!remote || length(groups) + (!is.null(total)) < 2L)
        return(list(groups = lapply(groups, function(g) dmp_groupsum(col, g)),
                    total = if (!is.null(total)) dm_colsum(total)))
    fr <- dm_frame_bind(do.call(dm_column_frame, c(list(col, total), lapply(groups, `[[`, "codes"))))
    gs <- vapply(groups, function(g) unclass(g$codes)$sql, "")
    vs <- unclass(col)$sql
    ts <- if (!is.null(total)) unclass(total)$sql else character()
    k <- length(vs)
    sums <- function(e, p) c(sprintf("SUM(CAST(%s AS DOUBLE)) AS %s%d", e, p, seq_along(e)),
                             sprintf("SUM(CASE WHEN (%s) IS NULL THEN 1 ELSE 0 END) AS %sn%d", e, p, seq_along(e)))
    sets <- c(sprintf("(%s)", gs), if (length(ts)) "()")
    res <- dm_query(fr, sprintf("SELECT %s FROM (%s) AS s WHERE %s GROUP BY GROUPING SETS (%s)",
                                paste(c(sprintf("%s AS g%d", gs, seq_along(gs)), sums(vs, "v"), sums(ts, "t")),
                                      collapse = ", "),
                                fr$src$select, fr$filter, paste(sets, collapse = ", ")))
    ## (codes are never missing: the rows of grouping j are those with g_j)
    gk <- lapply(seq_along(gs), function(j) !is.na(res[[paste0("g", j)]]))
    out <- lapply(seq_along(groups), function(j) {
        r <- gk[[j]]
        for (l in seq_along(gs)[-j]) r <- r & !gk[[l]]
        m <- matrix(0, groups[[j]]$G, k)
        i <- dm_group_index(groups[[j]], res[[paste0("g", j)]][r])
        for (u in seq_len(k)) {
            v <- as.double(res[[paste0("v", u)]][r])
            v[as.numeric(res[[paste0("vn", u)]][r]) > 0] <- NA_real_
            m[i, u] <- v
        }
        m
    })
    tot <- NULL
    if (length(ts)) {
        r <- Reduce(`&`, lapply(gk, `!`))
        tot <- vapply(seq_along(ts), function(u) {
            if (!any(r)) return(0)
            if (as.numeric(res[[paste0("tn", u)]][r][1L]) > 0) return(NA_real_)
            v <- as.double(res[[paste0("t", u)]][r][1L])
            if (is.na(v)) 0 else v
        }, 0)
    }
    list(groups = out, total = tot)
}

## weighted counts by pair of coded groups: Ga x Gb dgCMatrix
dmp_crosstab_codes <- function(a, b, w = NULL) {
    if (!dm_is_column(a$codes) && !dm_is_column(b$codes) && !dm_is_column(w)) {
        ## (the kernel returns the transposes of the crosstabs of its first
        ## argument with the others, and is fastest when the first is
        ## sorted)
        ct <- .Call(C_dm_fe_crosstabs, b$codes, list(a$codes), if (is.null(w)) NULL else as.double(w), b$G,
                    a$G, design_threads())[[1L]]
        return(dm_ct_matrix(ct, a$G, b$G))
    }
    res <- dm_group_query(list(a$codes, b$codes), if (is.null(w)) list() else list(w))
    ## (the pairs, one per row, through the kernel of the rows in memory)
    ct <- .Call(C_dm_fe_crosstabs, dm_group_index(b, res$g2), list(dm_group_index(a, res$g1)), as.double(res$v1), b$G,
                a$G, design_threads())[[1L]]
    dm_ct_matrix(ct, a$G, b$G)
}

## A crosstab from the kernels (compressed columns list(p, i, x)) as an
## nr x nc dgCMatrix.
dm_ct_matrix <- function(ct, nr, nc) {
    methods::new("dgCMatrix", i = ct[[2L]], p = ct[[1L]], x = ct[[3L]], Dim = c(as.integer(nr), as.integer(nc)))
}

## Column sums of a design within several coded groups at once:
## dmp_colsums_groups(x, groups, extra = NULL), x possibly scaled (s * x),
## groups a list of dmc_codes(), extra a column (k columns) or NULL. Returns
## one list(S = G x p sums of the scaled design, n = G sums of the scale (the
## counts without one), extra = G x k sums of scale * extra) per grouping.
## In memory: one pass of the structured kernel for all groupings when it
## applies (dm_fe_struct_gsum), else the row kernels; remotely, one GROUP BY
## per grouping combines the design, scale, and extra column sums.
dmp_colsums_groups <- function(x, groups, extra = NULL) {
    if (inherits(x, "dm_remote")) {
        s <- x$scale
        one <- dm_column(dm_remote_unscaled(x), if (is.null(s)) "CAST(1 AS DOUBLE)" else s)
        totals <- if (is.null(extra)) one else dmc_cbind(one, one * extra)
        x <- dmp_rebase(x, list(totals, groups))
        return(lapply(groups, function(g) {
            cs <- dm_sql_colsums(x, unclass(g$codes)$sql, s, FALSE,
                                 group_totals = unclass(totals)$sql)
            i <- dm_group_index(g, cs$keys)
            S <- matrix(0, g$G, ncol(x))
            S[i, ] <- cs$sums
            ns <- matrix(0, g$G, ncol(cs$totals))
            ns[i, ] <- cs$totals
            list(S = S, n = ns[, 1L], extra = if (!is.null(extra)) ns[, -1L, drop = FALSE])
        }))
    }
    ## in memory: the plain design and its row scale
    w <- NULL
    D <- x
    if (dm_is_comb(x)) {
        if (!is.null(x$design2)) stop("A combination of two designs has no grouped sums here.", call. = FALSE)
        D <- x$design
        s <- x$scale
        w <- if (is.null(s)) NULL
             else if (is.double(s) && length(s) == nrow(D) && is.null(attributes(s))) s
             else rep_len(as.double(s), nrow(D))
    }
    dm_gsums_memory(D, w, groups, extra)$gs
}

## In memory: the grouped sums of dmp_colsums_groups() of design D scaled by
## w (NULL for 1), and the crosstabs (as dmp_crosstab(), sums of w) of the
## pairs of groupings `pairs` (a list of c(a, b), indices in `groups`: the
## G_b x G_a matrix of t(crosstab(a, b))) whose dense table is small enough,
## all in one pass of the structured kernel when it applies; list(gs, ct),
## ct with NULL for the pairs left out.
dm_gsums_memory <- function(D, w, groups, extra = NULL, pairs = list()) {
    G <- vapply(groups, function(g) g$G, 1L)
    codes <- lapply(groups, function(g) g$codes)
    ct <- vector("list", length(pairs))
    ## (the kernels take a vector as one column: no copy of a double)
    Y <- extra
    if (!is.null(Y) && !is.double(Y)) storage.mode(Y) <- "double"
    if (dm_colsums_struct_ok(D, G)) {
        ## dense crosstabs: up to 2^22 cells, and at most two per row
        dense <- vapply(pairs, function(pr) as.double(G[pr[1L]]) * G[pr[2L]] <= min(2^22, 2 * nrow(D)), NA)
        pm <- if (any(dense)) matrix(as.integer(unlist(pairs[dense])) - 1L, nrow = 2L)
        M <- .Call(C_dm_fe_struct_gsum, D$setup, D$vars, w, codes, as.integer(G), Y, TRUE, pm, design_threads())
        for (q in seq_along(which(dense))) {
            pr <- pairs[[which(dense)[q]]]
            ct[[which(dense)[q]]] <- dm_ct_matrix(M[[length(groups) + q]], G[pr[2L]], G[pr[1L]])
        }
        info <- dm_struct_info(D)
        icpt <- isTRUE(D$setup$intercept > 0L)
        p <- ncol(D)
        gs <- lapply(seq_along(groups), function(k) {
            m <- M[[k]]
            S <- matrix(0, G[k], p)
            q <- 1L
            for (t in seq_along(info)) {
                Tt <- info[[t]]
                if (t == 1L && icpt) {
                    S[, Tt$cols[1L, 1L]] <- m[1L, ]
                    next
                }
                L <- Tt$L
                ns <- nrow(Tt$cols)
                rows <- q + seq_len(L * ns)
                S[, Tt$cols] <- t(dm_struct_contract(m[rows, , drop = FALSE], Tt, G[k]))
                q <- q + L * ns
            }
            list(S = S, n = m[1L, ], extra = if (is.null(Y)) NULL else t(m[q + seq_len(NCOL(Y)), , drop = FALSE]))
        })
        return(list(gs = gs, ct = ct))
    }
    one <- if (is.null(w)) rep(1, nrow(D)) else w
    nt <- design_threads()
    gs <- lapply(seq_along(groups), function(k) {
        J <- as_jac(D)
        W <- matrix(one, ncol = 1L)
        S <- if (jac_kernel_ok(J)) t(dm_crossprod_kernel(J, W, codes[[k]], G[k]))
             else .Call(C_dm_fe_gsum, codes[[k]], G[k], one * as.matrix(D), nt)
        ns <- .Call(C_dm_fe_gsum, codes[[k]], G[k], cbind(one, if (!is.null(Y)) one * Y), nt)
        list(S = S, n = ns[, 1L], extra = if (!is.null(Y)) ns[, -1L, drop = FALSE])
    })
    list(gs = gs, ct = ct)
}

## The sums of a fit with fixed effects (fe.R): the grouped sums of
## dmp_colsums_groups(w * x, groups, extra) (x unscaled, w NULL for 1) and
## the crosstabs t(dmp_crosstab(groups[[b]], groups[[a]], w)) (G_a x G_b
## sums of w, transposed) for each pair c(a, b) of `pairs`: list(gs, ct).
## In memory one pass of the structured kernel for the grouped sums and
## the small crosstabs (dense tables), and one pass for the crosstabs of
## each grouping a with all its b (dm_fe_crosstabs(): the rows by level of
## a); remotely the queries of dmp_colsums_groups() and dmp_crosstab().
dmp_fe_sums <- function(x, w, groups, extra = NULL, pairs = list()) {
    if (inherits(x, "dm_remote") || dm_is_comb(x)) {
        gs <- dmp_colsums_groups(if (is.null(w)) x else dmp_scale(x, w), groups, extra)
        ct <- lapply(pairs, function(pr) dmp_crosstab(groups[[pr[2L]]], groups[[pr[1L]]], w))
        return(list(gs = gs, ct = ct))
    }
    if (!is.null(w)) w <- as.double(w)
    r <- dm_gsums_memory(x, w, groups, extra, pairs)
    ## the other crosstabs: one pass per grouping a, for all its b
    left <- which(vapply(r$ct, is.null, NA))
    for (a in unique(vapply(pairs[left], `[`, 1L, 1L))) {
        qs <- left[vapply(pairs[left], `[`, 1L, 1L) == a]
        bs <- vapply(pairs[qs], `[`, 1L, 2L)
        cts <- .Call(C_dm_fe_crosstabs, groups[[a]]$codes, lapply(groups[bs], `[[`, "codes"), w,
                     groups[[a]]$G, vapply(groups[bs], function(g) g$G, 1L), design_threads())
        for (j in seq_along(qs)) r$ct[[qs[j]]] <- dm_ct_matrix(cts[[j]], groups[[bs[j]]]$G, groups[[a]]$G)
    }
    r
}

dm_colsums_struct_ok <- function(D, G) {
    if (!isTRUE(getOption("design.matrix.structured", TRUE))) return(FALSE)
    if (!dm_struct_ok(D)) return(FALSE)
    info <- dm_struct_info(D)
    Q <- 2 + sum(vapply(info, function(Tt) Tt$L * nrow(Tt$cols), 0))
    all(Q * as.double(G) <= 2^25)
}

## ---- more column functions (requested by the GLM work) ----

## one column of a multi-column column
dmc_col <- function(x, j) {
    if (!dm_is_column(x)) return(if (is.matrix(x)) x[, j] else x)
    dm_column(unclass(x)$frame, unclass(x)$sql[j])
}

## a constant column over the rows of a column's frame (or a design)
dmc_const <- function(like, value) {
    if (inherits(like, "dm_remote")) return(dm_column(like, dm_sql_double(value)))
    if (!dm_is_column(like)) return(rep(as.double(value), NROW(like)))
    dm_column(unclass(like)$frame, dm_sql_double(value))
}

dmc_lgamma <- function(x) {
    if (!dm_is_column(x)) return(lgamma(x))
    dm_column(unclass(x)$frame, sprintf("lgamma(%s)", unclass(x)$sql))
}

## round half to even, as R's round(x)
dmc_round <- function(x) {
    if (!dm_is_column(x)) return(round(x))
    dm_column(unclass(x)$frame, sprintf("round_even(%s, 0)", unclass(x)$sql))
}

## is.finite, NaN-safe (in DuckDB, NaN = NaN and NaN > 0 are true)
dmc_is_finite <- function(x) {
    if (!dm_is_column(x)) return(is.finite(x))
    dm_column(unclass(x)$frame, sprintf("COALESCE(isfinite(CAST(%s AS DOUBLE)), FALSE)", unclass(x)$sql))
}

## TRUE where x is not missing, FALSE where it is NA or NaN (NULL remotely);
## infinite values count as present (= !is.na(x); was dmc_notnan, renamed
## in tapx: it never tested NaN alone. The key version is dmc_complete().) Never NA, so that !dmc_present(x)
## counts the missing values. (In memory one pass in C for a plain vector.)
dmc_present <- function(x) {
    if (!dm_is_column(x)) {
        if (is.null(attributes(x)) && (is.double(x) || is.integer(x) || is.logical(x)))
            return(.Call(C_dm_notnan, x))
        return(!is.na(x))
    }
    dm_column(unclass(x)$frame, sprintf("COALESCE(NOT isnan(CAST(%s AS DOUBLE)), FALSE)", unclass(x)$sql),
             notnull = TRUE)
}

## The rows where a logical column is TRUE (with negate, where it is not
## TRUE: FALSE or missing), as a column to sum: 1 or 0. In memory a logical
## vector, which sums as 0/1.
dmc_count <- function(cond, negate = FALSE) {
    if (!dm_is_column(cond)) {
        if (anyNA(cond)) cond <- !is.na(cond) & cond
        return(if (negate) !cond else cond)
    }
    sql <- unclass(cond)$sql
    dm_column(unclass(cond)$frame, sprintf(if (negate) "CASE WHEN %s THEN 0 ELSE 1 END" else "CASE WHEN %s THEN 1 ELSE 0 END", sql),
             notnull = TRUE)
}

## A column computed once: remotely, its expressions become columns of a
## new frame (a projection over the frame of the column, one level deeper,
## as products are), so that later expressions refer to them instead of
## repeating them (DuckDB expands macros such as dm_qnorm() inline, and
## nested expressions would grow without bound). Columns frozen separately
## from one frame live on sibling frames, which cannot be combined (a frozen
## column is not replayed on another branch): freeze them together
## (dmc_cbind()). In memory, the column itself. With table =
## TRUE, the new frame is a temporary table of the rows used (computed once
## for all later queries, which then scan it: for columns that many queries
## read, such as the working weights and responses of a GLM with fixed
## effects); its rows keep the filter of the design.
dmc_freeze <- function(col, table = FALSE) {
    if (!dm_is_column(col)) return(col)
    fr <- dm_frame_bind(unclass(col)$frame)
    sql <- unclass(col)$sql
    nm <- paste0("__dm_f", substr(dm_tmp_name("f"), 3L, 10L), "_", seq_along(sql))
    x <- dm_remote_extend(dm_remote_unscaled(fr), stats::setNames(sql, nm), "column")
    if (table) x <- dm_remote_to_table(x)
    x$depth <- dm_frame_depth(fr) + 1L
    x <- dm_frame_child(x, fr)
    dm_column(x, sprintf("s.%s", dm_qi(x, nm)))
}

## The frame of a remote design as a temporary table of its rows used (with
## every column of the frame; the design's filter stays valid on them).
## (`keep`: the columns of the frame to copy, NULL for all; `name`: a
## permanent table of that name instead, replaced if `overwrite`, for the
## public col_compute())
dm_remote_to_table <- function(x, keep = NULL, name = NULL, overwrite = FALSE) {
    tn <- dm_qi(x, if (is.null(name)) dm_tmp_name("dm_mat") else name)
    create <- if (is.null(name)) "CREATE TEMPORARY TABLE" else if (overwrite) "CREATE OR REPLACE TABLE" else "CREATE TABLE"
    sel <- if (is.null(keep)) "*" else paste(dm_qi(x, keep), collapse = ", ")
    dm_execute(x, sprintf("%s %s AS SELECT %s FROM (%s) AS s WHERE %s", create, tn, sel, x$src$select, x$filter))
    x$src$select <- sprintf("SELECT * FROM %s", tn)
    x$src$base_table <- tn
    if (!is.null(keep) && !is.null(x$src$proto)) {
        x$src$proto <- x$src$proto[keep]
        x$src$types <- x$src$types[keep]
    }
    x
}

## The columns of the frame of remote design x to copy for the design and
## the columns `cols` (a list of columns, dmc_codes() groupings or NULL):
## the columns of the data (whose names do not start with "__dm_"; results
## on the rows keep them), and the internal columns (lookups, products of
## earlier steps) that the design (the SQL and inputs of its variables, its
## scale) or `cols` read: those whose names occur in their SQL.
dm_remote_used <- function(x, cols) {
    sql <- function(v) {
        if (is.null(v)) return(character())
        if (dm_is_codes(v)) return(sql(v$codes))
        if (dm_is_column(v)) return(unclass(v)$sql)
        if (is.list(v)) return(unlist(lapply(v, sql), use.names = FALSE))
        character()
    }
    txt <- c(sql(cols), unlist(lapply(x$vars, function(v) c(v$sql, v$col, v$inputs))), x$scale %||% character())
    txt <- paste(as.character(txt), collapse = " ")
    nm <- names(x$src$proto)
    nm[!startsWith(nm, "__dm_") | vapply(nm, grepl, NA, x = txt, fixed = TRUE)]
}

dmp_materialize <- function(x, cols = NULL) {
    if (!inherits(x, "dm_remote")) return(x)
    ## (the columns' pending recipes bound: the frame copied has them)
    if (!is.null(cols)) x <- dmp_rebase(x, cols)
    sc <- x$scale
    keep <- if (!is.null(cols) && !is.null(x$src$proto)) dm_remote_used(x, cols)
    if (!length(keep)) keep <- NULL
    x <- dm_remote_to_table(dm_remote_unscaled(x), keep)
    x$filter <- "TRUE"
    x$depth <- 0L
    x$path <- character()
    x$steps <- NULL
    x$cache <- new.env(parent = emptyenv())
    x$id <- dm_tmp_name("dm")
    x$scale <- sc
    x
}

## ---- inputs of the methods of a design ----

## Which rows of a design matrix are finite: rowSums(!is.finite(as.matrix(x)))
## == 0, computed without forming the dense matrix (TRUE for each row whose
## entries are all finite, FALSE for a row with any NA, NaN or infinite entry,
## including entries that overflow, such as the product of two large numeric
## variables, and, for a combination of designs, scaled sums that overflow).
finite_rows <- function(x) {
    ## a remote design: a logical column (every entry of the row finite)
    if (inherits(x, "dm_remote")) {
        P <- dmp_mult(x, diag(1, ncol(x)))
        return(dm_column(unclass(P)$frame,
                        paste(sprintf("COALESCE(isfinite(%s), FALSE)", unclass(P)$sql), collapse = " AND ")))
    }
    checkmate::assert_class(x, "design_matrix")
    J <- as_jac(x)
    if (!jac_kernel_ok(J)) {
        return(unname(rowSums(!is.finite(as.matrix(x))) == 0))
    }
    do.call(.Call, c(list(C_dm_jac_finite), jac_args(J), list(design_threads())))
}

## A lower bound on the absolute value of the non-zero entries of x (0 when
## no bound is available: always for a remote design).
dmp_min_nonzero <- function(x) {
    if (inherits(x, "dm_remote")) return(0)
    terms <- dm_terms(x)
    out <- Inf
    for (t in terms) {
        d <- t$design
        b <- if (kernel_ok(d)) .Call(C_dm_min_nonzero, d$setup, d$vars) else {
            X <- as.matrix(d)
            a <- abs(X[X != 0 & is.finite(X)])
            if (length(a)) min(a) else 1
        }
        s <- t$scale
        if (!is.null(s)) {
            a <- abs(s[s != 0 & is.finite(s)])
            b <- b * (if (length(a)) min(a) else 1)
        }
        out <- min(out, b)
    }
    ## a sum of two designs can cancel to small values
    if (length(terms) > 1L) 0 else out
}

## The `y` of the qr() methods as a column of design x: in memory a finite
## double matrix (from a vector or matrix with nrow(x) rows); remotely a
## remote column from column names, a DBI::SQL() expression or a remote
## column. Returns list(Y, vec (was y one column?), rows, cols (names)).
dmp_y <- function(x, y) {
    if (!inherits(x, "dm_remote")) return(dmqr_y(y, nrow(x)))
    D <- dm_remote_unscaled(x)
    if (dm_is_column(y)) return(list(Y = y, vec = dm_ncols(y) == 1L, cols = NULL))
    if (inherits(y, "SQL")) return(list(Y = dm_column(D, as.character(y)), vec = TRUE, cols = NULL))
    if (is.character(y)) {
        dm_remote_column(D, y, "y", several = TRUE)
        return(list(Y = dm_column(D, sprintf("CAST(s.%s AS DOUBLE)", dm_qi(D, y))), vec = length(y) == 1L,
                    cols = if (length(y) > 1L) y))
    }
    stop("For a remote design, `y` must name columns of its data, or be a DBI::SQL() expression.", call. = FALSE)
}

## ---- rows of a design ----

## The design of the rows where `keep` (a logical column, NA counting as
## FALSE) is TRUE: x[keep, ] in memory; remotely, the same design with its
## filter restricted (a new design: columns of x must be re-made on it,
## for instance with dmc_on()).
dmp_subset <- function(x, keep) {
    if (!inherits(x, "dm_remote")) {
        if (dm_is_column(keep)) stop("A remote column cannot select rows of a design in memory.", call. = FALSE)
        return(x[if (is.logical(keep) && !anyNA(keep)) keep else !is.na(keep) & keep, ])
    }
    if (inherits(keep, "SQL")) keep <- dm_column(x, as.character(keep))
    if (!dm_is_column(keep))
        stop("Rows of a remote design are selected by a logical remote column or a DBI::SQL() predicate.",
             call. = FALSE)
    out <- dm_rebase(x, dm_frame_later(unclass(keep)$frame, x))
    out$filter <- sprintf("(%s) AND COALESCE(%s, FALSE)", out$filter, unclass(keep)$sql)
    out$n <- as.numeric(dm_query(out, sprintf("SELECT COUNT(*) AS n FROM (%s) AS s WHERE %s",
                                              out$src$select, out$filter))$n)
    out$cache <- new.env(parent = emptyenv())
    out$id <- dm_tmp_name("dm")
    out
}

## The same column (its SQL) over the rows of design x (remote; in memory,
## the column itself). Valid when its expressions read columns of the
## source of x.
dmc_on <- function(col, x) {
    if (!dm_is_column(col)) return(col)
    ## (a pending column reads columns that its recipes make: bind it, on
    ## the frame that x copies, before)
    if (dm_column_unbound(col))
        stop("internal error: dmc_on() of a pending column (bind it first: dmp_rebase(), dmc_bind())",
             call. = FALSE)
    dm_column(x, unclass(col)$sql)
}

## ---- arguments of a fit (round 7 of 2hba) ----

## Columns for the arguments of a fit over the rows of a remote design:
## exprs, a named list of unevaluated arguments (each a column name, a
## DBI::SQL() expression, or an R expression over the columns, evaluated in
## env); logical: the names of logical arguments (a subset). R expressions
## with a SQL translation run in the database; others are evaluated in R,
## batch by batch, into a temporary table (dm_materialize()), which then is
## the source of the design. Returns list(x = the design (possibly on that
## new source), cols = named list of columns (NULL for NULL arguments),
## evaluated = names evaluated in R, with their reasons).
## In memory, the arguments are evaluated in `data` (a data frame, the rows
## of x) and returned as vectors.
dmp_resolve <- function(x, exprs, env, logical = character(), data = NULL) {
    if (!inherits(x, "dm_remote")) {
        cols <- lapply(exprs, function(e) {
            if (is.null(e)) return(NULL)
            if (is.character(e) && length(e) == 1L) return(data[[e]])
            eval(e, data, env)
        })
        return(list(x = x, cols = cols, evaluated = list()))
    }
    con <- dm_con(x)
    items <- list()
    later <- list()
    sqls <- list()
    for (nm in names(exprs)) {
        e <- exprs[[nm]]
        if (is.null(e)) next
        a <- dm_resolve_or_eval(e, env, x$src$types, con, nm, bool = nm %in% logical)
        if (!is.null(a$expr)) {
            items[[nm]] <- a
            later[[paste0("__dm_eval_", nm)]] <- a$expr
        } else {
            sqls[[nm]] <- a$sql
        }
    }
    if (length(later)) {
        x$src <- dm_materialize(x$src, later, env, logical = paste0("__dm_eval_", logical), chunk_size = x$chunk_size)
        for (nm in names(items)) {
            col <- DBI::dbQuoteIdentifier(con, paste0("__dm_eval_", nm))
            sqls[[nm]] <- if (nm %in% logical) sprintf("s.%s", col) else sprintf("CAST(s.%s AS DOUBLE)", col)
        }
    }
    cols <- lapply(names(exprs), function(nm) if (is.null(sqls[[nm]])) NULL else dm_column(x, sqls[[nm]]))
    names(cols) <- names(exprs)
    list(x = x, cols = cols, evaluated = items)
}

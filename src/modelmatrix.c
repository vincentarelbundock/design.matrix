/*
 *  design.matrix: exact, faster construction of model matrices.
 *
 *  dm_setup() and dm_fill() are a port of modelmatrix() from R's
 *  src/library/stats/src/model.c:
 *
 *    Copyright (C) 1997--2025  The R Core Team
 *    Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *
 *  The port splits R's single function in two.  dm_setup() performs every
 *  check, contrast evaluation, warning and column-name computation in the
 *  same order as R.  dm_fill() allocates and fills the matrix.  The fill
 *  performs exactly the same floating point operations, cell by cell, as
 *  R does; only the order in which rows are visited changes (row chunks,
 *  optionally on several threads), so the result is bitwise identical.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */

#include <string.h>
#include <limits.h>
#include <R.h>
#include <Rinternals.h>
#include <Rversion.h>
#include "dm.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/* Use stats' own message catalogue so messages are identical to stats. */
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext("stats", String)
#else
#define _(String) (String)
#endif

static R_INLINE int isUnordered_int(SEXP s)
{
    return (TYPEOF(s) == INTSXP && inherits(s, "factor") && !inherits(s, "ordered"));
}

static R_INLINE int isOrdered_int(SEXP s)
{
    return (TYPEOF(s) == INTSXP && inherits(s, "factor") && inherits(s, "ordered"));
}

#define BUFSIZE 4096

static char *AppendString(char *buf, const char *str)
{
    while (*str)
	*buf++ = *str++;
    *buf = '\0';
    return buf;
}

static char *AppendInteger(char *buf, int i)
{
    char str[32];
    snprintf(str, sizeof(str), "%d", i);
    return AppendString(buf, str);
}

static SEXP ColumnNames(SEXP x)
{
    SEXP dn = getAttrib(x, R_DimNamesSymbol);
    if (dn == R_NilValue || length(dn) < 2)
	return R_NilValue;
    else
	return VECTOR_ELT(dn, 1);
}

static const char *type_name(SEXP x)
{
#if R_VERSION >= R_Version(4, 4, 0)
    return R_typeToChar(x);
#else
    return type2char(TYPEOF(x));
#endif
}

/* Names of the elements of the list returned by dm_setup(). */
static const char *setup_names[] = {
    "factors", "contr1", "contr2", "nlevs", "columns", "count",
    "rhs_response", "intercept", "xnames", "assign", "nVar", "nterms", "contrast_rows", ""
};

/* .Call(C_dm_setup, terms, vars, rho)
 * The body follows R's modelmatrix() up to (but excluding) the allocation
 * of the design matrix. */
SEXP dm_setup(SEXP terms, SEXP vars, SEXP rho)
{
    SEXP factors, vnames, assign, xnames, count, contr1, contr2, nlevs,
	ordered, columns, variable, var_i, expr, x;
    int fik, first, i, j, k, kk, ll, n, nc, nterms, nVar;
    int intrcept, risponse, indx, rhs_response;
    double dk, dnc;
    char buf[BUFSIZE] = "\0";
    char *bufp;
    const char *addp;

    intrcept = asLogical(getAttrib(terms, install("intercept")));
    if (intrcept == NA_INTEGER)
	intrcept = 0;

    risponse = asLogical(getAttrib(terms, install("response")));
    if (risponse == NA_INTEGER)
	risponse = 0;

    nVar = nterms = 0;
    PROTECT(factors = duplicate(getAttrib(terms, install("factors"))));
    if (length(factors) == 0) {
	nVar = 1;
	nterms = 0;
    }
    else if (isInteger(factors) && isMatrix(factors)) {
	nVar = nrows(factors);
	nterms = ncols(factors);
    }
    else error(_("invalid '%s' argument"), "terms");

    vnames = getAttrib(factors, R_DimNamesSymbol);
    if (length(factors) > 0) {
	if (length(vnames) < 1 ||
	    (nVar - intrcept > 0 && !isString(VECTOR_ELT(vnames, 0))))
	    error(_("invalid '%s' argument"), "terms");
	vnames = VECTOR_ELT(vnames, 0);
    }
    PROTECT(vnames);

    if (!isNewList(vars) || length(vars) < nVar)
	error(_("invalid model frame"));
    if (length(vars) == 0)
	error(_("do not know how many cases"));

    n = nrows(VECTOR_ELT(vars, 0));

    PROTECT(variable = allocVector(VECSXP, nVar));
    PROTECT(nlevs = allocVector(INTSXP, nVar));
    PROTECT(ordered = allocVector(LGLSXP, nVar));
    PROTECT(columns = allocVector(INTSXP, nVar));

    for (i = 0; i < nVar; i++) {
	var_i = SET_VECTOR_ELT(variable, i, VECTOR_ELT(vars, i));
	if (nrows(var_i) != n)
	    error(_("variable lengths differ (found for variable %d)"), i);
	if (isOrdered_int(var_i)) {
	    LOGICAL(ordered)[i] = 1;
	    if ((INTEGER(nlevs)[i] = nlevels(var_i)) < 1)
		error(_("variable %d has no levels"), i + 1);
	    INTEGER(columns)[i] = ncols(var_i);
	}
	else if (isUnordered_int(var_i)) {
	    LOGICAL(ordered)[i] = 0;
	    if ((INTEGER(nlevs)[i] = nlevels(var_i)) < 1)
		error(_("variable %d has no levels"), i + 1);
	    INTEGER(columns)[i] = ncols(var_i);
	}
	else if (isLogical(var_i)) {
	    LOGICAL(ordered)[i] = 0;
	    INTEGER(nlevs)[i] = 2;
	    INTEGER(columns)[i] = ncols(var_i);
	}
	else if (isNumeric(var_i)) {
	    /* R coerces here; we only need the column count, and dm_vars()
	       performs the same coercion when the values are needed. */
	    LOGICAL(ordered)[i] = 0;
	    INTEGER(nlevs)[i] = 0;
	    INTEGER(columns)[i] = ncols(var_i);
	}
	else {
	    LOGICAL(ordered)[i] = 0;
	    INTEGER(nlevs)[i] = 0;
	    INTEGER(columns)[i] = ncols(var_i);
	}
    }

    /* No intercept: the first factor found is coded by dummy variables. */
    if (!intrcept) {
	for (j = 0; j < nterms; j++) {
	    for (i = risponse; i < nVar; i++) {
		if (INTEGER(nlevs)[i] > 1
		    && INTEGER(factors)[i + j * nVar] > 0) {
		    INTEGER(factors)[i + j * nVar] = 2;
		    goto alldone;
		}
	    }
	}
    }
 alldone:
    ;

    PROTECT(contr1 = allocVector(VECSXP, nVar));
    PROTECT(contr2 = allocVector(VECSXP, nVar));

    PROTECT(expr = lang3(install("contrasts"), R_NilValue, allocVector(LGLSXP, 1)));

    for (i = 0; i < nVar; i++) {
	if (INTEGER(nlevs)[i]) {
	    k = 0;
	    for (j = 0; j < nterms; j++) {
		if (INTEGER(factors)[i + j * nVar] == 1)
		    k |= 1;
		else if (INTEGER(factors)[i + j * nVar] == 2)
		    k |= 2;
	    }
	    SETCADR(expr, VECTOR_ELT(variable, i));
	    if (k & 1) {
		LOGICAL(CADDR(expr))[0] = 1;
		SET_VECTOR_ELT(contr1, i, eval(expr, rho));
	    }
	    if (k & 2) {
		LOGICAL(CADDR(expr))[0] = 0;
		SET_VECTOR_ELT(contr2, i, eval(expr, rho));
	    }
	}
    }

    rhs_response = -1;
    if (risponse > 0)
	for (j = 0; j < nterms; j++)
	    if (INTEGER(factors)[risponse - 1 + j * nVar]) {
		for (i = 0, k = 0; i < nVar; i++)
		    k += INTEGER(factors)[i + j * nVar] > 0;
		if (k == 1) {
		    rhs_response = j;
		    break;
		}
	    }

    PROTECT(count = allocVector(INTSXP, nterms));
    if (intrcept)
	dnc = 1;
    else
	dnc = 0;
    for (j = 0; j < nterms; j++) {
	if (j == rhs_response) {
	    warning(_("the response appeared on the right-hand side and was dropped"));
	    INTEGER(count)[j] = 0;
	    continue;
	}
	dk = 1;
	for (i = 0; i < nVar; i++) {
	    if (INTEGER(factors)[i + j * nVar]) {
		if (INTEGER(nlevs)[i]) {
		    switch (INTEGER(factors)[i + j * nVar]) {
		    case 1:
			dk *= ncols(VECTOR_ELT(contr1, i));
			break;
		    case 2:
			dk *= ncols(VECTOR_ELT(contr2, i));
			break;
		    }
		}
		else dk *= INTEGER(columns)[i];
	    }
	}
	if (dk > INT_MAX) error(_("term %d would require %.0g columns"), j + 1, dk);
	INTEGER(count)[j] = (int) dk;
	dnc = dnc + dk;
    }
    if (dnc > INT_MAX) error(_("matrix would require %.0g columns"), dnc);
    nc = (int) dnc;

    PROTECT(assign = allocVector(INTSXP, nc));
    k = 0;
    if (intrcept) INTEGER(assign)[k++] = 0;
    for (j = 0; j < nterms; j++) {
	if (INTEGER(count)[j] <= 0)
	    warning(_("problem with term %d in model.matrix: no columns are assigned"),
		    j + 1);
	for (i = 0; i < INTEGER(count)[j]; i++)
	    INTEGER(assign)[k++] = j + 1;
    }

    PROTECT(xnames = allocVector(STRSXP, nc));

    k = 0;
    if (intrcept)
	SET_STRING_ELT(xnames, k++, mkChar("(Intercept)"));

    for (j = 0; j < nterms; j++) {
	if (j == rhs_response) continue;
	for (kk = 0; kk < INTEGER(count)[j]; kk++) {
	    first = 1;
	    indx = kk;
	    bufp = &buf[0];
	    buf[0] = '\0';
	    for (i = 0; i < nVar; i++) {
		ll = INTEGER(factors)[i + j * nVar];
		if (ll) {
		    var_i = VECTOR_ELT(variable, i);
		    if (!first) {
			if (strlen(buf) + 1 < BUFSIZE)
			    bufp = AppendString(bufp, ":");
			else
			    warning(_("term names will be truncated"));
		    }
		    first = 0;
		    if (isFactor(var_i) || isLogical(var_i)) {
			if (ll == 1) {
			    x = ColumnNames(VECTOR_ELT(contr1, i));
			    ll = ncols(VECTOR_ELT(contr1, i));
			}
			else {
			    x = ColumnNames(VECTOR_ELT(contr2, i));
			    ll = ncols(VECTOR_ELT(contr2, i));
			}
			addp = translateChar(STRING_ELT(vnames, i));
			if (strlen(buf) + strlen(addp) < BUFSIZE)
			    bufp = AppendString(bufp, addp);
			else
			    warning(_("term names will be truncated"));
			if (x == R_NilValue) {
			    if (strlen(buf) + 10 < BUFSIZE)
				bufp = AppendInteger(bufp, indx % ll + 1);
			    else
				warning(_("term names will be truncated"));
			} else {
			    addp = translateChar(STRING_ELT(x, indx % ll));
			    if (strlen(buf) + strlen(addp) < BUFSIZE)
				bufp = AppendString(bufp, addp);
			    else
				warning(_("term names will be truncated"));
			}
		    } else if (isComplex(var_i)) {
			error(_("complex variables are not currently allowed in model matrices"));
		    } else if (isNumeric(var_i)) {
			x = ColumnNames(var_i);
			ll = ncols(var_i);
			addp = translateChar(STRING_ELT(vnames, i));
			if (strlen(buf) + strlen(addp) < BUFSIZE)
			    bufp = AppendString(bufp, addp);
			else
			    warning(_("term names will be truncated"));
			if (ll > 1) {
			    if (x == R_NilValue) {
				if (strlen(buf) + 10 < BUFSIZE)
				    bufp = AppendInteger(bufp, indx % ll + 1);
				else
				    warning(_("term names will be truncated"));
			    } else {
				addp = translateChar(STRING_ELT(x, indx % ll));
				if (strlen(buf) + strlen(addp) < BUFSIZE)
				    bufp = AppendString(bufp, addp);
				else
				    warning(_("term names will be truncated"));
			    }
			}
		    } else
			error(_("variables of type '%s' are not allowed in model matrices"),
			      type_name(var_i));
		    indx /= ll;
		}
	    }
	    SET_STRING_ELT(xnames, k++, mkChar(buf));
	}
    }

    /* R coerces the contrast matrices to double at fill time; doing it once
       here gives the same values. */
    for (i = 0; i < nVar; i++) {
	if (VECTOR_ELT(contr1, i) != R_NilValue)
	    SET_VECTOR_ELT(contr1, i, coerceVector(VECTOR_ELT(contr1, i), REALSXP));
	if (VECTOR_ELT(contr2, i) != R_NilValue)
	    SET_VECTOR_ELT(contr2, i, coerceVector(VECTOR_ELT(contr2, i), REALSXP));
    }

    SEXP out = PROTECT(mkNamed(VECSXP, setup_names));
    SET_VECTOR_ELT(out, 0, factors);
    SET_VECTOR_ELT(out, 1, contr1);
    SET_VECTOR_ELT(out, 2, contr2);
    SET_VECTOR_ELT(out, 3, nlevs);
    SET_VECTOR_ELT(out, 4, columns);
    SET_VECTOR_ELT(out, 5, count);
    SET_VECTOR_ELT(out, 6, ScalarInteger(rhs_response));
    SET_VECTOR_ELT(out, 7, ScalarInteger(intrcept));
    SET_VECTOR_ELT(out, 8, xnames);
    SET_VECTOR_ELT(out, 9, assign);
    SET_VECTOR_ELT(out, 10, ScalarInteger(nVar));
    SET_VECTOR_ELT(out, 11, ScalarInteger(nterms));
    /* One lazy row index per variable and contrast coding, shared by plans. */
    SET_VECTOR_ELT(out, 12, allocVector(VECSXP, 2 * (R_xlen_t) nVar));
    UNPROTECT(13);
    return out;
}

/* The variables as the fill and the kernels read them: numeric (non-factor
   integer or double) variables coerced to double, exactly as R does.  A
   variable that no term uses (the response) is not read: it is left as it
   is, without a coerced copy at every call. */
SEXP dm_vars(SEXP setup, SEXP vars)
{
    int nVar = asInteger(VECTOR_ELT(setup, 10));
    const int *nlevs = INTEGER(VECTOR_ELT(setup, 3));
    const int *count = INTEGER(VECTOR_ELT(setup, 5));
    int rhs_response = asInteger(VECTOR_ELT(setup, 6)),
	nterms = asInteger(VECTOR_ELT(setup, 11));
    const int *fac = nterms > 0 ? INTEGER(VECTOR_ELT(setup, 0)) : NULL;
    SEXP variable = PROTECT(allocVector(VECSXP, nVar));
    for (int i = 0; i < nVar; i++) {
	SEXP v = VECTOR_ELT(vars, i);
	int used = 0;
	for (int k = 0; k < nterms && !used; k++)
	    if (k != rhs_response && count[k] > 0 && fac[i + k * nVar]) used = 1;
	if (used && !nlevs[i] && !isLogical(v) && isNumeric(v) && TYPEOF(v) != REALSXP)
	    v = coerceVector(v, REALSXP);
	SET_VECTOR_ELT(variable, i, v);
    }
    UNPROTECT(1);
    return variable;
}

/* Unrolling the inner loops lets GCC at -O2 match the speed of R's own
   routine (as compiled by R's build); it does not change any operation. */
#if defined(__GNUC__) && !defined(__clang__)
# define DM_UNROLL _Pragma("GCC unroll 4")
#else
# define DM_UNROLL
#endif

/* ---- Fill kernels: R's firstfactor/addfactor/firstvar/addvar, with a
   leading dimension so that they can work on a block of rows.  Each cell
   gets the same operations as in R. */

static void firstfactor(double *x, R_xlen_t ldx, int nr, int ncx,
			const double *c, int nrc, int ncc, const int *v, int adj)
{
    for (int j = 0; j < ncc; j++) {
	double *xj = &x[j * ldx];
	const double *cj = &c[j * (R_xlen_t)nrc];
	DM_UNROLL
	for (int i = 0; i < nr; i++)
	    if (v[i] == NA_INTEGER) xj[i] = NA_REAL;
	    else xj[i] = cj[v[i] - 1 + adj];
    }
}

static void addfactor(double *x, R_xlen_t ldx, int nr, int ncx,
		      const double *c, int nrc, int ncc, const int *v, int adj)
{
    for (int k = ncc - 1; k >= 0; k--) {
	for (int j = 0; j < ncx; j++) {
	    double *xj = &x[j * ldx];
	    double *yj = &x[(k * (R_xlen_t)ncx + j) * ldx];
	    const double *ck = &c[k * (R_xlen_t)nrc];
	    DM_UNROLL
	    for (int i = 0; i < nr; i++)
		if (v[i] == NA_INTEGER) yj[i] = NA_REAL;
		else yj[i] = ck[v[i] - 1 + adj] * xj[i];
	}
    }
}

static void firstvar(double *x, R_xlen_t ldx, int nr, int ncx,
		     const double *c, R_xlen_t ldc, int ncc)
{
    for (int j = 0; j < ncc; j++) {
	double *xj = &x[j * ldx];
	const double *cj = &c[j * ldc];
	DM_UNROLL
	for (int i = 0; i < nr; i++)
	    xj[i] = cj[i];
    }
}

static void addvar(double *x, R_xlen_t ldx, int nr, int ncx,
		   const double *c, R_xlen_t ldc, int ncc)
{
    for (int k = ncc - 1; k >= 0; k--) {
	for (int j = 0; j < ncx; j++) {
	    double *xj = &x[j * ldx];
	    double *yj = &x[(k * (R_xlen_t)ncx + j) * ldx];
	    const double *ck = &c[k * ldc];
	    DM_UNROLL
	    for (int i = 0; i < nr; i++)
		yj[i] = ck[i] * xj[i];
	}
    }
}

/* Plain-C view of the variables, built on the main thread so the threaded
   fill never calls the R API (REAL()/INTEGER() can allocate for ALTREP). */
typedef struct {
    const int *code;      /* factor or logical codes, or NULL */
    const double *num;    /* numeric values, or NULL */
    int adj, ncol;
    const double *c[3];   /* contrast matrices for fik = 1, 2 */
    int nrc[3], ncc[3];
} fillvar;

static void fill_rows(double *rx, R_xlen_t n, int r0, int nr, int intrcept,
		      int nVar, int nterms, int rhs_response, const int *factors,
		      const int *nlevs, const int *columns, const int *count,
		      const fillvar *fv)
{
    int jstart, jnext;
    if ((jnext = jstart = intrcept) != 0)
	for (int i = 0; i < nr; i++)
	    rx[r0 + i] = 1.0;

    for (int k = 0; k < nterms; k++) {
	if (k == rhs_response) continue;
	/* A term with no columns (e.g. an interaction with a zero-column
	   matrix) writes nothing.  R's own routine still writes the other
	   variables' columns here, past the end of the matrix. */
	if (count[k] == 0) continue;
	for (int i = 0; i < nVar; i++) {
	    if (columns[i] == 0)
		continue;
	    int fik = factors[i + k * nVar];
	    if (!fik) continue;
	    const fillvar *v = &fv[i];
	    double *base = &rx[jstart * n + r0];
	    if (jnext == jstart) {
		if (nlevs[i] > 0) {
		    firstfactor(base, n, nr, jnext - jstart, v->c[fik],
				v->nrc[fik], v->ncc[fik], v->code + r0, v->adj);
		    jnext = jnext + v->ncc[fik];
		}
		else {
		    firstvar(base, n, nr, jnext - jstart, v->num + r0, n, v->ncol);
		    jnext = jnext + v->ncol;
		}
	    }
	    else {
		if (nlevs[i] > 0) {
		    addfactor(base, n, nr, jnext - jstart, v->c[fik],
			      v->nrc[fik], v->ncc[fik], v->code + r0, v->adj);
		    jnext = jnext + (jnext - jstart) * (v->ncc[fik] - 1);
		}
		else {
		    addvar(base, n, nr, jnext - jstart, v->num + r0, n, v->ncol);
		    jnext = jnext + (jnext - jstart) * (v->ncol - 1);
		}
	    }
	}
	jstart = jnext;
    }
}

#define FILL_CHUNK 4096

/* .Call(C_dm_fill, setup, vars, nthreads)
 * `vars` is the list given to R's C routine (it carries the row names). */
SEXP dm_fill(SEXP setup, SEXP vars, SEXP nthreads_)
{
    SEXP factors = VECTOR_ELT(setup, 0), contr1 = VECTOR_ELT(setup, 1),
	contr2 = VECTOR_ELT(setup, 2);
    const int *nlevs = INTEGER(VECTOR_ELT(setup, 3)),
	*columns = INTEGER(VECTOR_ELT(setup, 4)),
	*count = INTEGER(VECTOR_ELT(setup, 5));
    int rhs_response = asInteger(VECTOR_ELT(setup, 6)),
	intrcept = asInteger(VECTOR_ELT(setup, 7)),
	nVar = asInteger(VECTOR_ELT(setup, 10)),
	nterms = asInteger(VECTOR_ELT(setup, 11)),
	nthreads = asInteger(nthreads_);
    SEXP xnames = VECTOR_ELT(setup, 8), assign = VECTOR_ELT(setup, 9);
    int nc = LENGTH(xnames);
    const int *fac = nterms > 0 ? INTEGER(factors) : NULL;

    int n = nrows(VECTOR_ELT(vars, 0));
    SEXP rnames = PROTECT(getAttrib(vars, R_RowNamesSymbol));
    SEXP variable = PROTECT(dm_vars(setup, vars));

    /* Which variables does the fill read? */
    fillvar *fv = (fillvar *) R_alloc(nVar > 0 ? nVar : 1, sizeof(fillvar));
    for (int i = 0; i < nVar; i++) {
	fillvar *v = &fv[i];
	memset(v, 0, sizeof(fillvar));
	int used = 0;
	for (int k = 0; k < nterms; k++)
	    if (k != rhs_response && fac[i + k * nVar]) used = 1;
	if (!used || columns[i] == 0) continue;
	SEXP var_i = VECTOR_ELT(variable, i);
	v->ncol = ncols(var_i);
	if (nlevs[i] > 0) {
	    v->code = INTEGER(var_i);
	    v->adj = isLogical(var_i) ? 1 : 0;
	    for (int f = 1; f <= 2; f++) {
		SEXP cm = VECTOR_ELT(f == 1 ? contr1 : contr2, i);
		if (cm == R_NilValue) continue;
		v->c[f] = REAL(cm);
		v->nrc[f] = nrows(cm);
		v->ncc[f] = ncols(cm);
	    }
	} else {
	    v->num = REAL(var_i);
	}
    }

    SEXP x = PROTECT(allocMatrix(REALSXP, n, nc));
    double *rx = REAL(x);
    /* One thread: a single chunk, as R does. Several threads: a few chunks
       per thread, each at least FILL_CHUNK rows. */
    nthreads = dm_worker_count(nthreads, n);
    int chunk = n > 0 ? n : 1;
    if (nthreads > 1) {
	int target = (int) dm_block_count(dm_block_count(n, nthreads), 4);
	chunk = target > FILL_CHUNK ? target : FILL_CHUNK;
    }
    int nchunks = (int) dm_block_count(n, chunk);
    nthreads = dm_worker_count(nthreads, nchunks);

    if (nthreads == 1) {
	/* outside any parallel region, so the compiler sees R's plain loop */
	if (n > 0)
	    fill_rows(rx, n, 0, n, intrcept, nVar, nterms, rhs_response, fac,
		      nlevs, columns, count, fv);
    } else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
	for (int ch = 0; ch < nchunks; ch++) {
	    int r0 = (int) ((R_xlen_t) ch * chunk);
	    int nr = n - r0 < chunk ? n - r0 : chunk;
	    fill_rows(rx, n, r0, nr, intrcept, nVar, nterms, rhs_response, fac,
		      nlevs, columns, count, fv);
	}
    }

    SEXP tnames = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(tnames, 0, rnames);
    SET_VECTOR_ELT(tnames, 1, xnames);
    setAttrib(x, R_DimNamesSymbol, tnames);
    setAttrib(x, install("assign"), assign);
    UNPROTECT(4);
    return x;
}

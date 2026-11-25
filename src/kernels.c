/*
 *  design.matrix: algebra on a design without forming the dense matrix.
 *
 *  A design is described by the list returned by dm_setup() plus the
 *  variables.  Row i of the dense matrix is the intercept followed, for
 *  each term, by the product of its variables' columns, the first variable
 *  varying fastest (as in R's model.matrix).  row_entries() lists the
 *  columns of row i that can be non-zero, with their values, computing each
 *  product in the same order as the dense fill.
 *
 *  Zero contrast entries are skipped only when that cannot change the
 *  result: when every value in the term is finite and the running products
 *  cannot overflow (a bound is checked row by row).  Otherwise every column
 *  of the term is listed, so that 0 * Inf gives NaN as in the dense matrix.
 *  A missing factor code makes every column of the term NA.  With `force`,
 *  every column of the row is listed.
 *
 *  A Jacobian row is s1[i] * x1_i + s2[i] * x2_i, where x1 and x2 are rows
 *  of two designs with the same columns (s2 and x2 optional).  Entries of
 *  the two designs are combined column by column, with the same arithmetic
 *  as the dense Jacobian, before any product is taken.
 *
 *  Copyright (C) 2026 Vincent Arel-Bundock.  GPL (>= 2).
 */

#include <string.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include "dm.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define MAXC 32
#define SAFE_BOUND (DBL_MAX / 16)

typedef struct {
    int isfac;
    const int *code;
    int adj;
    const double *cm;     /* contrast matrix, nrc x nc */
    int nrc;
    const double *num;    /* numeric values, n x nc */
    int nc;
    int stride;           /* column stride of this component within the term */
    const double *sp_start;        /* sparse contrast rows: entries for level r are   */
    const int *sp_col;          /* sp_col/sp_val[sp_start[r] .. sp_start[r + 1]) */
    const double *sp_val;
    const double *rowmax;       /* max |entry| of each contrast row */
    const int *rowbad;          /* contrast row has a non-finite entry */
} kcomp;

typedef struct {
    int ncomp, off, count;
    kcomp comp[MAXC];
} kterm;

typedef struct {
    R_xlen_t n;
    int p, intercept, nterm;
    kterm *term;
} kctx;

/* Build once on the first kernel use. Store ordinary R vectors so a plan
   can share them across data sets and survive serialization. */
static SEXP contrast_rows(SEXP cm)
{
    int nr = nrows(cm), nc = ncols(cm);
    const double *v = REAL(cm);
    R_xlen_t nnz = 0;
    for (R_xlen_t i = 0; i < XLENGTH(cm); i++)
        if (v[i] != 0.0) nnz++;

    SEXP out = PROTECT(allocVector(VECSXP, 5));
    /* Double offsets also index contrast matrices with more than INT_MAX
       entries. Each row still has at most nc entries. */
    double *start = REAL(SET_VECTOR_ELT(out, 0, allocVector(REALSXP, (R_xlen_t) nr + 1)));
    int *col = INTEGER(SET_VECTOR_ELT(out, 1, allocVector(INTSXP, nnz)));
    double *val = REAL(SET_VECTOR_ELT(out, 2, allocVector(REALSXP, nnz)));
    double *rowmax = REAL(SET_VECTOR_ELT(out, 3, allocVector(REALSXP, nr)));
    int *rowbad = INTEGER(SET_VECTOR_ELT(out, 4, allocVector(INTSXP, nr)));

    nnz = 0;
    for (int r = 0; r < nr; r++) {
        start[r] = (double) nnz;
        rowmax[r] = 0.0;
        rowbad[r] = 0;
        for (int j = 0; j < nc; j++) {
            double x = v[r + (R_xlen_t) j * nr];
            if (!R_FINITE(x)) rowbad[r] = 1;
            else if (fabs(x) > rowmax[r]) rowmax[r] = fabs(x);
            if (x != 0.0) {
                col[nnz] = j;
                val[nnz++] = x;
            }
        }
    }
    start[nr] = (double) nnz;
    UNPROTECT(1);
    return out;
}

/* Build the context on the main thread; no R API calls happen later.
   Leaves one object PROTECTed (the coerced variables, which the context
   points into); the caller unprotects it. */
static void build_ctx(kctx *K, SEXP setup, SEXP vars)
{
    SEXP factors = VECTOR_ELT(setup, 0), contr1 = VECTOR_ELT(setup, 1),
	contr2 = VECTOR_ELT(setup, 2);
    const int *nlevs = INTEGER(VECTOR_ELT(setup, 3)),
	*columns = INTEGER(VECTOR_ELT(setup, 4)),
	*count = INTEGER(VECTOR_ELT(setup, 5));
    int rhs_response = asInteger(VECTOR_ELT(setup, 6)),
	nVar = asInteger(VECTOR_ELT(setup, 10)),
	nterms = asInteger(VECTOR_ELT(setup, 11));
    const int *fac = nterms > 0 ? INTEGER(factors) : NULL;

    SEXP variable = PROTECT(dm_vars(setup, vars));
    SEXP cache;
    if (XLENGTH(setup) > 12) {
        cache = VECTOR_ELT(setup, 12);
    } else {
        /* Designs saved before the cache was introduced use a per-call
           index. Keep its vectors alive with the protected variable list. */
        cache = PROTECT(allocVector(VECSXP, 2 * (R_xlen_t) nVar));
        setAttrib(variable, install("contrast_rows"), cache);
        UNPROTECT(1);
    }
    K->n = nrows(VECTOR_ELT(vars, 0));
    K->intercept = asInteger(VECTOR_ELT(setup, 7));
    K->p = LENGTH(VECTOR_ELT(setup, 8));
    K->term = (kterm *) R_alloc(nterms > 0 ? nterms : 1, sizeof(kterm));
    K->nterm = 0;

    int off = K->intercept;
    for (int k = 0; k < nterms; k++) {
	if (k == rhs_response || count[k] <= 0) continue;
	kterm *T = &K->term[K->nterm++];
	T->ncomp = 0;
	T->off = off;
	T->count = count[k];
	off += count[k];
	int stride = 1;
	for (int i = 0; i < nVar; i++) {
	    int fik = fac[i + k * nVar];
	    if (!fik || columns[i] == 0) continue;
	    if (T->ncomp >= MAXC)
		error("internal error: term with more than %d variables", MAXC);
	    kcomp *C = &T->comp[T->ncomp++];
	    memset(C, 0, sizeof(kcomp));
	    SEXP var_i = VECTOR_ELT(variable, i);
	    if (nlevs[i] > 0) {
		SEXP cm = VECTOR_ELT(fik == 1 ? contr1 : contr2, i);
		C->isfac = 1;
		C->code = INTEGER(var_i);
		C->adj = isLogical(var_i) ? 1 : 0;
		C->cm = REAL(cm);
		C->nrc = nrows(cm);
		C->nc = ncols(cm);
		R_xlen_t slot = (fik == 1 ? 0 : nVar) + (R_xlen_t) i;
		SEXP rows = VECTOR_ELT(cache, slot);
		if (isNull(rows))
		    rows = SET_VECTOR_ELT(cache, slot, contrast_rows(cm));
		C->sp_start = REAL(VECTOR_ELT(rows, 0));
		C->sp_col = INTEGER(VECTOR_ELT(rows, 1));
		C->sp_val = REAL(VECTOR_ELT(rows, 2));
		C->rowmax = REAL(VECTOR_ELT(rows, 3));
		C->rowbad = INTEGER(VECTOR_ELT(rows, 4));
	    } else {
		if (TYPEOF(var_i) != REALSXP)
		    error("internal error: variable %d is not double", i + 1);
		C->num = REAL(var_i);
		C->nc = ncols(var_i);
	    }
	    C->stride = stride;
	    stride *= C->nc;
	}
    }
}

/* Entries of row i: writes columns and values, returns how many.
   Buffers must hold at least K->p entries. */
static int row_entries(const kctx *K, R_xlen_t i, int force, int *cols, double *vals)
{
    int m = 0;
    if (K->intercept) {
	cols[m] = 0;
	vals[m++] = 1.0;
    }
    for (int t = 0; t < K->nterm; t++) {
	const kterm *T = &K->term[t];
	int nc = T->ncomp;

	/* NA codes; can zero contrast entries be skipped safely? */
	int na = 0, dense = force;
	double bound = 1.0;
	for (int c = 0; c < nc; c++) {
	    const kcomp *C = &T->comp[c];
	    if (C->isfac) {
		int code = C->code[i];
		if (code == NA_INTEGER) {
		    na = 1;
		    continue;
		}
		int r = code - 1 + C->adj;
		if (C->rowbad[r]) dense = 1;
		else bound *= (C->rowmax[r] > 1.0 ? C->rowmax[r] : 1.0);
	    } else {
		double mx = 0.0;
		for (int j = 0; j < C->nc; j++) {
		    double v = C->num[i + (R_xlen_t)j * K->n];
		    if (!R_FINITE(v)) dense = 1;
		    else if (fabs(v) > mx) mx = fabs(v);
		}
		bound *= (mx > 1.0 ? mx : 1.0);
	    }
	}
	if (na) {
	    for (int j = 0; j < T->count; j++) {
		cols[m] = T->off + j;
		vals[m++] = NA_REAL;
	    }
	    continue;
	}
	/* every running product is at most `bound` in absolute value */
	if (!(bound < SAFE_BOUND)) dense = 1;

	int nch[MAXC], pos[MAXC], lev[MAXC];
	int empty = 0;
	for (int c = 0; c < nc; c++) {
	    const kcomp *C = &T->comp[c];
	    pos[c] = 0;
	    if (C->isfac) {
		lev[c] = C->code[i] - 1 + C->adj;
		nch[c] = dense ? C->nc : C->sp_start[lev[c] + 1] - C->sp_start[lev[c]];
	    } else {
		nch[c] = C->nc;
	    }
	    if (nch[c] == 0) empty = 1;
	}
	if (empty) continue;

	/* odometer, first component fastest */
	for (;;) {
	    double val = 0.0;
	    int idx = 0;
	    for (int c = 0; c < nc; c++) {
		const kcomp *C = &T->comp[c];
		int col;
		double v;
		if (C->isfac) {
		    if (dense) {
			col = pos[c];
			v = C->cm[lev[c] + (R_xlen_t)col * C->nrc];
		    } else {
			R_xlen_t e = (R_xlen_t) C->sp_start[lev[c]] + pos[c];
			col = C->sp_col[e];
			v = C->sp_val[e];
		    }
		} else {
		    col = pos[c];
		    v = C->num[i + (R_xlen_t)col * K->n];
		}
		/* same order as the fill: new value times running product */
		val = (c == 0) ? v : v * val;
		idx += col * C->stride;
	    }
	    cols[m] = T->off + idx;
	    vals[m++] = val;
	    int c = 0;
	    while (c < nc && ++pos[c] == nch[c]) {
		pos[c] = 0;
		c++;
	    }
	    if (c == nc) break;
	}
    }
    return m;
}

/* A Jacobian: one or two designs with per-row scales. */
typedef struct {
    kctx d[2];
    const double *s[2];
    int scalar[2];
    int nd;       /* number of designs = number of objects left PROTECTed */
    int p;
} jctx;

static void build_jac(jctx *J, SEXP setup1, SEXP vars1, SEXP s1,
		      SEXP setup2, SEXP vars2, SEXP s2)
{
    build_ctx(&J->d[0], setup1, vars1);
    J->scalar[0] = !isNull(s1) && XLENGTH(s1) == 1;
    J->s[0] = isNull(s1) ? NULL : REAL(s1);
    J->nd = 1;
    if (!isNull(setup2)) {
	build_ctx(&J->d[1], setup2, vars2);
	J->scalar[1] = !isNull(s2) && XLENGTH(s2) == 1;
	J->s[1] = isNull(s2) ? NULL : REAL(s2);
	J->nd = 2;
	if (J->d[1].p != J->d[0].p || J->d[1].n != J->d[0].n)
	    error("the two designs must have the same dimensions");
    }
    J->p = J->d[0].p;
    if (J->s[0] && !J->scalar[0] && XLENGTH(s1) != J->d[0].n)
	error("'scale' must have one value per row");
    if (J->nd == 2 && J->s[1] && !J->scalar[1] && XLENGTH(s2) != J->d[0].n)
	error("'scale2' must have one value per row");
}

/* Per-thread scratch space. */
typedef struct {
    int *cols;
    double *vals;
    int *where;   /* p entries, all -1 between calls */
} scratch;

static scratch *alloc_scratch(int nthreads, int p, int nd)
{
    int cap = nd * (p > 0 ? p : 1);
    scratch *S = (scratch *) R_alloc(nthreads, sizeof(scratch));
    for (int t = 0; t < nthreads; t++) {
	S[t].cols = (int *) R_alloc(cap, sizeof(int));
	S[t].vals = (double *) R_alloc(cap, sizeof(double));
	S[t].where = NULL;
	if (nd == 2) {
	    S[t].where = (int *) R_alloc(p > 0 ? p : 1, sizeof(int));
	    for (int j = 0; j < p; j++) S[t].where[j] = -1;
	}
    }
    return S;
}

/* Row i of the Jacobian, each column at most once.  The value of a column
   present in both designs is s1 * x1 + s2 * x2, as in the dense Jacobian.
   A non-finite scale lists every column of that design (Inf * 0 is NaN). */
static int jac_entries(const jctx *J, R_xlen_t i, int force, scratch *S)
{
    int *cols = S->cols;
    double *vals = S->vals;
    const double *s0 = J->s[0];
    double scale0 = s0 ? s0[J->scalar[0] ? 0 : i] : 1.0;
    int f0 = force || !R_FINITE(scale0);
    int m = row_entries(&J->d[0], i, f0, cols, vals);
    if (s0) {
	for (int u = 0; u < m; u++) vals[u] = scale0 * vals[u];
    }
    if (J->nd == 1) return m;

    const double *s1 = J->s[1];
    double scale1 = s1 ? s1[J->scalar[1] ? 0 : i] : 1.0;
    int f1 = force || !R_FINITE(scale1);
    int *where = S->where;
    for (int u = 0; u < m; u++) where[cols[u]] = u;
    int m2 = row_entries(&J->d[1], i, f1, cols + m, vals + m);
    int out = m;
    for (int u = m; u < m + m2; u++) {
	int c = cols[u];
	double v = s1 ? scale1 * vals[u] : vals[u];
	if (where[c] >= 0) {
	    vals[where[c]] = vals[where[c]] + v;
	} else {
	    cols[out] = c;
	    vals[out] = v;
	    where[c] = out++;
	}
    }
    for (int u = 0; u < out; u++) where[cols[u]] = -1;
    return out;
}

static R_INLINE int thread_id(void)
{
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

/* J %*% B, B is p x K (finite): returns n x K */
/* Rows to compute: all (rows = NULL) or the given 1-based indices. */
static R_xlen_t nrows_of(SEXP rows, R_xlen_t n)
{
    return isNull(rows) ? n : XLENGTH(rows);
}

static R_INLINE R_xlen_t row_at(const int *rows, R_xlen_t r)
{
    return rows ? (R_xlen_t) rows[r] - 1 : r;
}

static const int *rows_ptr(SEXP rows, R_xlen_t n)
{
    if (isNull(rows)) return NULL;
    const int *r = INTEGER(rows);
    for (R_xlen_t k = 0; k < XLENGTH(rows); k++)
	if (r[k] == NA_INTEGER || r[k] < 1 || r[k] > n)
	    error("internal error: row index out of range");
    return r;
}

#define DENSE_BLOCK 64

SEXP dm_jac_mult(SEXP setup1, SEXP vars1, SEXP s1, SEXP setup2, SEXP vars2,
		 SEXP s2, SEXP B, SEXP rows_, SEXP nthreads_)
{
    jctx J;
    build_jac(&J, setup1, vars1, s1, setup2, vars2, s2);
    const int *rows = rows_ptr(rows_, J.d[0].n);
    R_xlen_t n = nrows_of(rows_, J.d[0].n);
    int p = J.p, K = ncols(B);
    R_xlen_t nblocks = dm_block_count(n, DENSE_BLOCK);
    int nthreads = dm_worker_count(asInteger(nthreads_), K == 1 ? n : nblocks);
    if (nrows(B) != p) error("non-conformable arguments");
    const double *b = REAL(B);
    SEXP out = PROTECT(allocMatrix(REALSXP, n, K));
    double *o = REAL(out);
    scratch *S = alloc_scratch(nthreads, p, J.nd);

    if (K == 1) {
#ifdef _OPENMP
#pragma omp parallel for if(nthreads > 1) num_threads(nthreads) schedule(static)
#endif
	for (R_xlen_t i = 0; i < n; i++) {
	    scratch *St = &S[thread_id()];
	    int m = jac_entries(&J, row_at(rows, i), 0, St);
	    double s = 0.0;
	    for (int u = 0; u < m; u++) s += St->vals[u] * b[St->cols[u]];
	    o[i] = s;
	}
	UNPROTECT(1 + J.nd);
	return out;
    }

    /* Several columns: B transposed (K values per row of B, contiguous),
       so that the loop over the columns of B runs over contiguous memory.
       Each output value still sums its terms in the same order.  Rows are
       computed in blocks and copied out column by column. */
    double *bt = (double *) R_alloc((size_t) p * K, sizeof(double));
    for (int k = 0; k < K; k++)
	for (int c = 0; c < p; c++) bt[(size_t) c * K + k] = b[c + (R_xlen_t) k * p];
    double *buf = (double *) R_alloc((size_t) nthreads * DENSE_BLOCK * K, sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for if(nthreads > 1) num_threads(nthreads) schedule(static)
#endif
    for (R_xlen_t blk = 0; blk < nblocks; blk++) {
	int t = thread_id();
	scratch *St = &S[t];
	double *A = buf + (size_t) t * DENSE_BLOCK * K;
	R_xlen_t i0 = blk * DENSE_BLOCK;
	int nb = (int) (n - i0 < DENSE_BLOCK ? n - i0 : DENSE_BLOCK);
	for (int r = 0; r < nb; r++) {
	    double *acc = A + (size_t) r * K;
	    for (int k = 0; k < K; k++) acc[k] = 0.0;
	    int m = jac_entries(&J, row_at(rows, i0 + r), 0, St);
	    for (int u = 0; u < m; u++) {
		double v = St->vals[u];
		const double *br = bt + (size_t) St->cols[u] * K;
		for (int k = 0; k < K; k++) acc[k] += v * br[k];
	    }
	}
	for (int k = 0; k < K; k++) {
	    double *ok = o + i0 + (R_xlen_t) k * n;
	    for (int r = 0; r < nb; r++) ok[r] = A[(size_t) r * K + k];
	}
    }
    UNPROTECT(1 + J.nd);
    return out;
}

/* t(J) %*% W, W is n x K; with `group` (codes 1..G), sums within groups:
   returns p x (K * G), column k + K * g.  A row with a non-finite weight
   lists every column (0 * NA is NA in the dense product).  Single-threaded:
   the order of the sums is fixed. */
SEXP dm_jac_crossprod(SEXP setup1, SEXP vars1, SEXP s1, SEXP setup2,
		      SEXP vars2, SEXP s2, SEXP W, SEXP group, SEXP G_, SEXP rows_)
{
    jctx J;
    build_jac(&J, setup1, vars1, s1, setup2, vars2, s2);
    const int *rows = rows_ptr(rows_, J.d[0].n);
    R_xlen_t n = nrows_of(rows_, J.d[0].n);
    int p = J.p, K = ncols(W), G = asInteger(G_);
    if (nrows(W) != n) error("non-conformable arguments");
    const double *w = REAL(W);
    const int *grp = isNull(group) ? NULL : INTEGER(group);
    if (grp && XLENGTH(group) != n) error("'group' must have one value per row");
    SEXP out = PROTECT(allocMatrix(REALSXP, p, (R_xlen_t)K * G));
    double *o = REAL(out);
    memset(o, 0, sizeof(double) * (size_t)p * K * G);
    scratch *S = alloc_scratch(1, p, J.nd);
    for (R_xlen_t i = 0; i < n; i++) {
	int g = 0;
	if (grp) {
	    g = grp[i];
	    if (g == NA_INTEGER) continue;
	    g -= 1;
	    if (g < 0 || g >= G) error("internal error: group code out of range");
	}
	int force = 0;
	for (int k = 0; k < K; k++)
	    if (!R_FINITE(w[i + (R_xlen_t)k * n])) force = 1;
	int m = jac_entries(&J, row_at(rows, i), force, S);
	for (int k = 0; k < K; k++) {
	    double wik = w[i + (R_xlen_t)k * n];
	    double *ok = o + (R_xlen_t)(k + K * g) * p;
	    for (int u = 0; u < m; u++) ok[S->cols[u]] += S->vals[u] * wik;
	}
    }
    UNPROTECT(1 + J.nd);
    return out;
}

/* diag(J V t(J)), V finite.  When a column of the row's projection J V
   could overflow, the whole row is computed densely, so that 0 * Inf
   gives NaN as in the dense computation. */
SEXP dm_jac_quad(SEXP setup1, SEXP vars1, SEXP s1, SEXP setup2, SEXP vars2,
		 SEXP s2, SEXP V, SEXP rows_, SEXP nthreads_)
{
    jctx J;
    build_jac(&J, setup1, vars1, s1, setup2, vars2, s2);
    const int *rows = rows_ptr(rows_, J.d[0].n);
    R_xlen_t n = nrows_of(rows_, J.d[0].n);
    int p = J.p, nthreads = dm_worker_count(asInteger(nthreads_), n);
    if (nrows(V) != p || ncols(V) != p) error("'V' must be a %d x %d matrix", p, p);
    const double *v = REAL(V);
    double vmax = 0.0;
    for (R_xlen_t j = 0; j < (R_xlen_t)p * p; j++)
	if (fabs(v[j]) > vmax) vmax = fabs(v[j]);
    SEXP out = PROTECT(allocVector(REALSXP, n));
    double *o = REAL(out);
    scratch *S = alloc_scratch(nthreads, p, J.nd);

#ifdef _OPENMP
#pragma omp parallel for if(nthreads > 1) num_threads(nthreads) schedule(static)
#endif
    for (R_xlen_t i = 0; i < n; i++) {
	scratch *St = &S[thread_id()];
	R_xlen_t ri = row_at(rows, i);
	int m = jac_entries(&J, ri, 0, St);
	double sumabs = 0.0;
	for (int a = 0; a < m; a++) sumabs += fabs(St->vals[a]);
	if (!(sumabs * vmax < SAFE_BOUND))
	    m = jac_entries(&J, ri, 1, St);   /* every column */
	const int *cols = St->cols;
	const double *vals = St->vals;
	/* the projection J V in double, like the matrix product; the final
	   sum in long double, like rowSums() */
	long double s = 0.0;
	for (int a = 0; a < m; a++) {
	    const double *va = v + (R_xlen_t)cols[a] * p;
	    double acc = 0.0;
	    for (int b = 0; b < m; b++) acc += vals[b] * va[cols[b]];
	    s += (long double) (vals[a] * acc);
	}
	o[i] = (double) s;
    }
    UNPROTECT(1 + J.nd);
    return out;
}

/* Compressed sparse column form of J: list(p, i, x), 0-based.  Exact
   zeros (+0 and -0) are dropped, as as(<matrix>, "CsparseMatrix") does;
   NA and NaN are kept.  Values are those of the dense matrix, bit for bit:
   each product is computed in the order of the fill. */
SEXP dm_jac_sparse(SEXP setup1, SEXP vars1, SEXP s1, SEXP setup2, SEXP vars2,
		   SEXP s2)
{
    jctx J;
    build_jac(&J, setup1, vars1, s1, setup2, vars2, s2);
    R_xlen_t n = J.d[0].n;
    int p = J.p;
    scratch *S = alloc_scratch(1, p, J.nd);

    /* pass 1: nonzeros per column */
    double *cnt = (double *) R_alloc(p + 1, sizeof(double));
    for (int j = 0; j <= p; j++) cnt[j] = 0.0;
    for (R_xlen_t i = 0; i < n; i++) {
	int m = jac_entries(&J, i, 0, S);
	for (int u = 0; u < m; u++)
	    if (!(S->vals[u] == 0.0)) cnt[S->cols[u]] += 1.0;
    }
    double total = 0.0;
    for (int j = 0; j < p; j++) total += cnt[j];
    if (total > INT_MAX)
	error("the sparse matrix would have %.0f non-zero entries, more than %d",
	      total, INT_MAX);

    SEXP colp = PROTECT(allocVector(INTSXP, (R_xlen_t)p + 1));
    int *cp = INTEGER(colp);
    cp[0] = 0;
    for (int j = 0; j < p; j++) cp[j + 1] = cp[j] + (int) cnt[j];
    int nnz = cp[p];
    SEXP ri = PROTECT(allocVector(INTSXP, nnz));
    SEXP rx = PROTECT(allocVector(REALSXP, nnz));
    int *ii = INTEGER(ri);
    double *xx = REAL(rx);

    /* pass 2: rows in order, so row indices are sorted within columns */
    int *next = (int *) R_alloc(p > 0 ? p : 1, sizeof(int));
    for (int j = 0; j < p; j++) next[j] = cp[j];
    for (R_xlen_t i = 0; i < n; i++) {
	int m = jac_entries(&J, i, 0, S);
	for (int u = 0; u < m; u++) {
	    double v = S->vals[u];
	    if (v == 0.0) continue;
	    int k = next[S->cols[u]]++;
	    ii[k] = (int) i;
	    xx[k] = v;
	}
    }

    SEXP out = PROTECT(allocVector(VECSXP, 3));
    SET_VECTOR_ELT(out, 0, colp);
    SET_VECTOR_ELT(out, 1, ri);
    SET_VECTOR_ELT(out, 2, rx);
    UNPROTECT(4 + J.nd);
    return out;
}

/* The dense matrix of J: every column of every row is computed (force),
   with the arithmetic of s1 * X1 + s2 * X2 on the dense matrices, so that
   finite values, infinities and zeros are those of the dense matrix, bit
   for bit.  Rows are computed in blocks into a row-major buffer, then
   copied out column by column, so that writes to the column-major result
   are contiguous. */
SEXP dm_jac_dense(SEXP setup1, SEXP vars1, SEXP s1, SEXP setup2, SEXP vars2,
		  SEXP s2, SEXP nthreads_)
{
    jctx J;
    build_jac(&J, setup1, vars1, s1, setup2, vars2, s2);
    R_xlen_t n = J.d[0].n;
    R_xlen_t nblocks = dm_block_count(n, DENSE_BLOCK);
    int p = J.p, nthreads = dm_worker_count(asInteger(nthreads_), nblocks);
    SEXP out = PROTECT(allocMatrix(REALSXP, n, p));
    double *o = REAL(out);
    scratch *S = alloc_scratch(nthreads, p, J.nd);
    double *buf = (double *) R_alloc((size_t) nthreads * DENSE_BLOCK * (p > 0 ? p : 1),
				     sizeof(double));
    int *nan = (int *) R_alloc(nthreads, sizeof(int));
    for (int t = 0; t < nthreads; t++) nan[t] = 0;

#ifdef _OPENMP
#pragma omp parallel for if(nthreads > 1) num_threads(nthreads) schedule(static)
#endif
    for (R_xlen_t b = 0; b < nblocks; b++) {
	int t = thread_id();
	scratch *St = &S[t];
	double *B = buf + (size_t) t * DENSE_BLOCK * p;
	R_xlen_t i0 = b * DENSE_BLOCK;
	int nb = (int) (n - i0 < DENSE_BLOCK ? n - i0 : DENSE_BLOCK);
	for (int r = 0; r < nb; r++) {
	    double *row = B + (size_t) r * p;
	    for (int j = 0; j < p; j++) row[j] = 0.0;
	    int m = jac_entries(&J, i0 + r, 1, St);
	    for (int u = 0; u < m; u++) {
		double v = St->vals[u];
		if (ISNAN(v)) nan[t] = 1;
		row[St->cols[u]] = v;
	    }
	}
	for (int j = 0; j < p; j++) {
	    double *oj = o + i0 + (R_xlen_t) j * n;
	    for (int r = 0; r < nb; r++) oj[r] = B[(size_t) r * p + j];
	}
    }
    /* Which NaN payload (NA or NaN, and its sign) survives an operation on
       two NaNs depends on the order of the operands, which the compiler
       may swap.  Flag any NaN, so that the caller can use R's arithmetic. */
    int any = 0;
    for (int t = 0; t < nthreads; t++) any |= nan[t];
    if (any) setAttrib(out, install("dm_nan"), ScalarLogical(1));
    UNPROTECT(1 + J.nd);
    return out;
}

/* Mean number of entries per row over the first k rows (a cost estimate). */
SEXP dm_jac_entries_mean(SEXP setup1, SEXP vars1, SEXP s1, SEXP setup2,
			 SEXP vars2, SEXP s2, SEXP k_)
{
    jctx J;
    build_jac(&J, setup1, vars1, s1, setup2, vars2, s2);
    R_xlen_t k = (R_xlen_t) asInteger(k_);
    if (k > J.d[0].n) k = J.d[0].n;
    scratch *S = alloc_scratch(1, J.p, J.nd);
    double tot = 0.0;
    for (R_xlen_t i = 0; i < k; i++) tot += jac_entries(&J, i, 0, S);
    UNPROTECT(J.nd);
    return ScalarReal(k > 0 ? tot / (double) k : 0.0);
}

/* Is every entry of each row finite?  Entries that the kernels skip are
   structural zeros, which are finite unless an input of the row is not
   (the caller checks the inputs); a non-finite scale lists every column of
   its design, so that Inf * 0 is seen.  This catches products and sums that
   overflow.  Returns a logical vector with one value per row. */
SEXP dm_jac_finite(SEXP setup1, SEXP vars1, SEXP s1, SEXP setup2, SEXP vars2,
		   SEXP s2, SEXP nthreads_)
{
    jctx J;
    build_jac(&J, setup1, vars1, s1, setup2, vars2, s2);
    R_xlen_t n = J.d[0].n;
    int p = J.p, nthreads = dm_worker_count(asInteger(nthreads_), n);
    SEXP out = PROTECT(allocVector(LGLSXP, n));
    int *o = LOGICAL(out);
    scratch *S = alloc_scratch(nthreads, p, J.nd);

#ifdef _OPENMP
#pragma omp parallel for if(nthreads > 1) num_threads(nthreads) schedule(static)
#endif
    for (R_xlen_t i = 0; i < n; i++) {
	scratch *St = &S[thread_id()];
	int m = jac_entries(&J, i, 0, St);
	int ok = 1;
	for (int u = 0; u < m; u++) {
	    if (!R_FINITE(St->vals[u])) { ok = 0; break; }
	}
	o[i] = ok;
    }
    UNPROTECT(1 + J.nd);
    return out;
}

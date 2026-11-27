/*
 *  xtx: kernels for groupings and grouped sums (R/primitives.R,
 *  R/primitives_lookup.R, R/groups_public.R), on which fits with fixed
 *  effects are built (package xtx.fit).
 *
 *  xtx_fe_gsum():       group sums of the columns of a matrix (codes 1..G);
 *  xtx_fe_crosstabs():  weighted crosstabs of a code vector with others, as
 *                      the compressed rows of G_a x G_b sparse matrices
 *                      with sorted, unique column indices;
 *  xtx_fe_struct_gsum(): grouped sums of the design's structured moments,
 *                      the weights and other columns, for every grouping
 *                      at once, and dense crosstabs of pairs of groupings,
 *                      in one pass;
 *  xtx_fe_code():       codes of a grouping variable and rows per level.
 *
 *  Copyright (C) 2026 Vincent Arel-Bundock.  GPL (>= 2).
 */

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>
#include <R.h>
#include <Rinternals.h>
#include "xtx.h"
#include "structured.h"

/* The kernels that index by code (accumulators, ranks, markers) check
   the codes first: an integer vector with every value in 1..G (no NA).
   One branch-free pass over the codes, small next to the kernels. */
void xtx_check_codes(SEXP codes_, int G)
{
    if (TYPEOF(codes_) != INTSXP) error("internal error: codes must be integer");
    R_xlen_t n = XLENGTH(codes_);
    const int *c = INTEGER_RO(codes_);
    unsigned int g = G > 0 ? (unsigned int) G : 0u;
    int ok = 1;
    for (R_xlen_t i = 0; i < n; i++) ok &= ((unsigned int) c[i] - 1u) < g;
    if (!ok) error("internal error: codes outside 1..G");
}

/* G x k matrix of sums of the columns of W (n x k) within codes 1..G.  On
   several threads: fixed row chunks with their own sums, added in chunk
   order (xtx_reduce_threads()); one thread is the single-threaded pass. */
static void gsum_rows(const int *codes, R_xlen_t n, int G, int k, const double *w,
		      R_xlen_t i0, R_xlen_t i1, double *o)
{
    for (int j = 0; j < k; j++) {
	const double *wj = w + (R_xlen_t)j * n;
	double *oj = o + (R_xlen_t)j * G;
	for (R_xlen_t i = i0; i < i1; i++) oj[codes[i] - 1] += wj[i];
    }
}

SEXP xtx_fe_gsum(SEXP codes_, SEXP G_, SEXP W, SEXP nthreads_)
{
    R_xlen_t n = XLENGTH(codes_);
    int G = asInteger(G_), k = ncols(W);
    if (nrows(W) != n) error("non-conformable arguments");
    xtx_check_codes(codes_, G);
    const int *codes = INTEGER(codes_);
    const double *w = REAL(W);
    SEXP out = PROTECT(allocMatrix(REALSXP, G, k));
    double *o = REAL(out);
    R_xlen_t len = (R_xlen_t) G * k;
    memset(o, 0, sizeof(double) * (size_t) len);
    int T = xtx_reduce_threads(asInteger(nthreads_), n, (double) n * k, (double) len);
    if (T == 1) {
	gsum_rows(codes, n, G, k, w, 0, n, o);
    } else {
	double **acc = xtx_chunk_acc(o, T, len);
#ifdef _OPENMP
#pragma omp parallel for num_threads(T) schedule(static)
#endif
	for (int t = 0; t < T; t++) {
	    if (t > 0) xtx_zero(acc[t], len);
	    gsum_rows(codes, n, G, k, w, xtx_chunk_start(n, T, t, 1), xtx_chunk_start(n, T, t + 1, 1), acc[t]);
	}
	xtx_combine(acc, T, len);
    }
    UNPROTECT(1);
    return out;
}

/* Crosstabs of codes a (1..Ga) with each of the codes b_1, ..., b_K
   (1..Gb_k), weights w (NULL for 1): for each k, list(p, j, x), the
   compressed rows (by a) of the Ga x Gb_k matrix of the sums of w, j
   0-based and sorted within rows.  The rows are taken by level of a (in
   row order: as they come when a is sorted, as in panel data, otherwise
   after a stable counting sort); within a level, a marker array over b_k
   finds each pair's entry, so that every entry is the sum of its rows in
   row order, and the level's entries are then sorted by b_k.  On several
   threads, each takes a range of levels of a (about n / T rows): the
   entries are the same on any number of threads. */
typedef struct { int j; double x; } ct_entry;

static int ct_cmp(const void *u, const void *v)
{
    int a = ((const ct_entry *) u)->j, b = ((const ct_entry *) v)->j;
    return (a > b) - (a < b);
}

/* sort entries jj[0..m), xx[0..m) by jj (unique) */
static void ct_sort(int *jj, double *xx, R_xlen_t m, ct_entry *tmp)
{
    int sorted = 1;
    for (R_xlen_t e = 1; e < m; e++) if (jj[e] < jj[e - 1]) { sorted = 0; break; }
    if (sorted) return;
    if (m <= 32) {
	for (R_xlen_t e = 1; e < m; e++) {
	    int jv = jj[e];
	    double xv = xx[e];
	    R_xlen_t f = e;
	    while (f > 0 && jj[f - 1] > jv) { jj[f] = jj[f - 1]; xx[f] = xx[f - 1]; f--; }
	    jj[f] = jv;
	    xx[f] = xv;
	}
	return;
    }
    for (R_xlen_t e = 0; e < m; e++) { tmp[e].j = jj[e]; tmp[e].x = xx[e]; }
    qsort(tmp, (size_t) m, sizeof(ct_entry), ct_cmp);
    for (R_xlen_t e = 0; e < m; e++) { jj[e] = tmp[e].j; xx[e] = tmp[e].x; }
}

SEXP xtx_fe_crosstabs(SEXP a_, SEXP bs_, SEXP w_, SEXP Ga_, SEXP Gb_, SEXP nthreads_)
{
    R_xlen_t n = XLENGTH(a_);
    int Ga = asInteger(Ga_), K = LENGTH(bs_);
    const int *a = INTEGER(a_), *Gb = INTEGER(Gb_);
    const double *w = isNull(w_) ? NULL : REAL(w_);
    if (w && XLENGTH(w_) != n) error("internal error: one weight per row");
    if (LENGTH(Gb_) != K) error("internal error: one G per code vector");
    xtx_check_codes(a_, Ga);
    const int **b = (const int **) R_alloc(K > 0 ? K : 1, sizeof(int *));
    for (int k = 0; k < K; k++) {
	if (XLENGTH(VECTOR_ELT(bs_, k)) != n) error("internal error: one code per row");
	xtx_check_codes(VECTOR_ELT(bs_, k), Gb[k]);
	b[k] = INTEGER(VECTOR_ELT(bs_, k));
    }
    /* the rows by level of a: start[g] .. start[g + 1] in ord (or in the
       rows themselves when a is sorted) */
    R_xlen_t *start = (R_xlen_t *) R_alloc((size_t) Ga + 1, sizeof(R_xlen_t));
    memset(start, 0, sizeof(R_xlen_t) * ((size_t) Ga + 1));
    int sorted = 1;
    for (R_xlen_t i = 0; i < n; i++) {
	start[a[i]]++;
	if (i > 0 && a[i] < a[i - 1]) sorted = 0;
    }
    for (int g = 1; g <= Ga; g++) start[g] += start[g - 1];
    int *ord = NULL;
    if (!sorted) {
	R_xlen_t *pos = (R_xlen_t *) R_alloc((size_t) Ga + 1, sizeof(R_xlen_t));
	memcpy(pos, start, sizeof(R_xlen_t) * ((size_t) Ga + 1));
	ord = (int *) R_alloc((size_t) n, sizeof(int));
	for (R_xlen_t i = 0; i < n; i++) ord[pos[a[i] - 1]++] = (int) i;
    }
#define ROW(t) (sorted ? (t) : (R_xlen_t) ord[t])
    int T = xtx_worker_count(asInteger(nthreads_), n / XTX_RED_MIN_ROWS);
    if (T > Ga) T = Ga > 0 ? Ga : 1;
    /* levels of thread t: lv[t] .. lv[t + 1], about n / T rows each */
    int *lv = (int *) R_alloc((size_t) T + 1, sizeof(int));
    lv[0] = 0;
    for (int t = 1; t < T; t++) {
	R_xlen_t target = (R_xlen_t) ((double) n * t / T);
	int g = lv[t - 1];
	while (g < Ga && start[g] < target) g++;
	lv[t] = g;
    }
    lv[T] = Ga;
    /* per thread and b_k: a marker (level of a + 1) and a position */
    int Gmax = 1;
    for (int k = 0; k < K; k++) if (Gb[k] > Gmax) Gmax = Gb[k];
    int **mark = (int **) R_alloc((size_t) T, sizeof(int *));
    R_xlen_t **pos = (R_xlen_t **) R_alloc((size_t) T, sizeof(R_xlen_t *));
    for (int t = 0; t < T; t++) {
	mark[t] = (int *) R_alloc((size_t) Gmax + 1, sizeof(int));
	pos[t] = (R_xlen_t *) R_alloc((size_t) Gmax + 1, sizeof(R_xlen_t));
    }
    SEXP out = PROTECT(allocVector(VECSXP, K));
    /* the largest level, for the sorting buffers */
    R_xlen_t mmax = 0;
    for (int g = 0; g < Ga; g++) if (start[g + 1] - start[g] > mmax) mmax = start[g + 1] - start[g];
    ct_entry **tmp = (ct_entry **) R_alloc((size_t) T, sizeof(ct_entry *));
    for (int t = 0; t < T; t++)
	tmp[t] = mmax > 32 ? (ct_entry *) R_alloc((size_t) mmax, sizeof(ct_entry)) : NULL;
    for (int k = 0; k < K; k++) {
	const int *bk = b[k];
	SEXP pk = PROTECT(allocVector(INTSXP, (R_xlen_t) Ga + 1));
	int *pp = INTEGER(pk);
	pp[0] = 0;
	/* count the pairs of each level */
#ifdef _OPENMP
#pragma omp parallel for num_threads(T) schedule(static)
#endif
	for (int t = 0; t < T; t++) {
	    int *mk = mark[t];
	    for (int g = 0; g <= Gb[k]; g++) mk[g] = 0;
	    for (int g = lv[t]; g < lv[t + 1]; g++) {
		int c = 0;
		for (R_xlen_t q = start[g]; q < start[g + 1]; q++) {
		    int bi = bk[ROW(q)];
		    if (mk[bi] != g + 1) { mk[bi] = g + 1; c++; }
		}
		pp[g + 1] = c;
	    }
	}
	double nnz = 0.0;
	for (int g = 0; g < Ga; g++) nnz += pp[g + 1];
	if (nnz > INT_MAX) error("the crosstab has too many entries");
	for (int g = 0; g < Ga; g++) pp[g + 1] += pp[g];
	SEXP jk = PROTECT(allocVector(INTSXP, pp[Ga]));
	SEXP xk = PROTECT(allocVector(REALSXP, pp[Ga]));
	int *jj = INTEGER(jk);
	double *xx = REAL(xk);
#ifdef _OPENMP
#pragma omp parallel for num_threads(T) schedule(static)
#endif
	for (int t = 0; t < T; t++) {
	    int *mk = mark[t];
	    R_xlen_t *ps = pos[t];
	    for (int g = 0; g <= Gb[k]; g++) mk[g] = 0;
	    for (int g = lv[t]; g < lv[t + 1]; g++) {
		R_xlen_t e0 = pp[g], e = e0;
		for (R_xlen_t q = start[g]; q < start[g + 1]; q++) {
		    R_xlen_t i = ROW(q);
		    int bi = bk[i];
		    if (mk[bi] != g + 1) {
			mk[bi] = g + 1;
			ps[bi] = e;
			jj[e] = bi - 1;
			xx[e] = 0.0;
			e++;
		    }
		    xx[ps[bi]] += w ? w[i] : 1.0;
		}
		ct_sort(jj + e0, xx + e0, e - e0, tmp[t]);
	    }
	}
	SEXP ok = SET_VECTOR_ELT(out, k, allocVector(VECSXP, 3));
	SET_VECTOR_ELT(ok, 0, pk);
	SET_VECTOR_ELT(ok, 1, jk);
	SET_VECTOR_ELT(ok, 2, xk);
	UNPROTECT(3);
    }
#undef ROW
    UNPROTECT(1);
    return out;
}

/* Is code b nested in code a (each level of b in one level of a)? */
SEXP xtx_fe_nested(SEXP a_, SEXP b_, SEXP Gb_)
{
    R_xlen_t n = XLENGTH(a_);
    int Gb = asInteger(Gb_);
    if (XLENGTH(b_) != n) error("internal error: one code per row");
    xtx_check_codes(b_, Gb);
    const int *a = INTEGER(a_), *b = INTEGER(b_);
    int *first = (int *) R_alloc((size_t) Gb + 1, sizeof(int));
    for (int g = 0; g <= Gb; g++) first[g] = 0;
    for (R_xlen_t i = 0; i < n; i++) {
	int *f = first + b[i];
	if (!*f) *f = a[i];
	else if (*f != a[i]) return ScalarLogical(0);
    }
    return ScalarLogical(1);
}

/* values[codes], NA for codes out of 1..length(values) or NA. */
SEXP xtx_fe_lookup(SEXP values_, SEXP codes_)
{
    /* values: a vector, or a G x K matrix (then an n x K matrix) */
    int mat = isMatrix(values_);
    R_xlen_t n = XLENGTH(codes_), G = mat ? nrows(values_) : XLENGTH(values_);
    int K = mat ? ncols(values_) : 1;
    const double *v = REAL_RO(values_);
    const int *c = INTEGER_RO(codes_);
    SEXP out = PROTECT(mat ? allocMatrix(REALSXP, n, K) : allocVector(REALSXP, n));
    double *o = REAL(out);
    for (int j = 0; j < K; j++) {
	const double *vj = v + (R_xlen_t) j * G;
	double *oj = o + (R_xlen_t) j * n;
	for (R_xlen_t i = 0; i < n; i++) {
	    int k = c[i];
	    oj[i] = (k == NA_INTEGER || k < 1 || k > G) ? NA_REAL : vj[k - 1];
	}
    }
    UNPROTECT(1);
    return out;
}

/* ------------------------------------------------------------------------
   Structured grouped sums (see R/absorb.R of xtx.fit, xf_fe_moments()).

   For each fixed effect k with codes c_k (1..G_k), a Q x G_k matrix whose
   column g holds, over the rows i with c_k(i) = g:
     row 0:                 sum w_i
     rows off_t + F + L_t u: sum w_i N_t(i)[u] over rows with F_t(i) = F
                            (term t of the design, its combined level F and
                            numeric product u; see structured.c), for every
                            term but the intercept (row 0 covers it);
     rows offY + j:         sum w_i Y_ij.
   One pass over the rows for all fixed effects: the Q sums of a row for
   effect k are contiguous, so that each row touches one or two cache lines
   per effect.  `terms = FALSE` skips the design (row 0 and Y only).

   In the same pass, crosstabs of pairs of fixed effects (`pairs`: NULL or
   a 2 x P integer matrix of 0-based indices (a, b)) accumulated densely:
   for each pair the sums of w over the rows of each (a, b) level pair, as
   the compressed columns (by a) list(p, i, x) of the G_b x G_a matrix
   (entries with no rows left out; w > 0, so every entry with rows is
   positive).  Each sum is over its rows in row order, as the sparse
   kernel (xtx_fe_crosstabs()) has them.

   On several threads: fixed row chunks (whole blocks of BLOCK rows), each
   with its own sums, added in chunk order (xtx_reduce_threads()); one
   thread is the single-threaded pass.  Returns the nk matrices, then the
   P crosstabs.
   ------------------------------------------------------------------------ */

typedef struct {
    const sctx *S;
    const int *toff;
    int nk, Q, K, offY, np;
    const int *const *codes;
    const double *w, *Y;
    const int *pa, *pb;  /* the pairs: codes indices */
    const int *G;
} gsum_ctx;

static void struct_gsum_rows(const gsum_ctx *C, R_xlen_t i0_, R_xlen_t i1, double **acc, double **pacc,
			     int **F, double **N, double *wb)
{
    const sctx *S = C->S;
    R_xlen_t n = S->n;
    int Q = C->Q, K = C->K, offY = C->offY;
    for (R_xlen_t i0 = i0_; i0 < i1; i0 += BLOCK) {
	int nb = (int) (i1 - i0 < BLOCK ? i1 - i0 : BLOCK);
	for (int t = 0; t < S->nt; t++)
	    if (C->toff[t] >= 0) term_block(&S->t[t], n, i0, nb, F[t], N[t]);
	for (int r = 0; r < nb; r++) wb[r] = C->w ? C->w[i0 + r] : 1.0;
	for (int k = 0; k < C->nk; k++) {
	    const int *ck = C->codes[k] + i0;
	    double *a = acc[k];
	    for (int r = 0; r < nb; r++) {
		double wi = wb[r];
		double *base = a + (size_t) Q * (ck[r] - 1);
		base[0] += wi;
		for (int t = 0; t < S->nt; t++) {
		    int o = C->toff[t];
		    if (o < 0) continue;
		    const sterm *Tt = &S->t[t];
		    int L = Tt->L, f = F[t][r];
		    const double *Nt = N[t] + r;
		    for (int u = 0; u < Tt->ns; u++)
			base[o + f + L * u] += wi * Nt[(size_t) u * BLOCK];
		}
		for (int j = 0; j < K; j++)
		    base[offY + j] += wi * C->Y[i0 + r + (R_xlen_t) j * n];
	    }
	}
	for (int q = 0; q < C->np; q++) {
	    const int *ca = C->codes[C->pa[q]] + i0, *cb = C->codes[C->pb[q]] + i0;
	    size_t Gb = (size_t) C->G[C->pb[q]];
	    double *a = pacc[q];
	    for (int r = 0; r < nb; r++) a[(size_t) (ca[r] - 1) * Gb + (cb[r] - 1)] += wb[r];
	}
    }
}

SEXP xtx_fe_struct_gsum(SEXP setup, SEXP vars, SEXP w_, SEXP codes_, SEXP G_, SEXP Y_, SEXP terms_,
		       SEXP pairs_, SEXP nthreads_)
{
    sctx S;
    int use_terms = asLogical(terms_);
    /* (build_sctx() leaves one object protected) */
    if (use_terms) {
	build_sctx(&S, setup, vars);
    } else {
	PROTECT(R_NilValue);
	S.nt = 0;
	S.t = NULL;
	S.n = XLENGTH(VECTOR_ELT(codes_, 0));
    }
    R_xlen_t n = S.n;
    int nk = LENGTH(codes_);
    if (LENGTH(G_) != nk) error("internal error: one G per code vector");
    const int *G = INTEGER(G_);
    const double *w = isNull(w_) ? NULL : REAL(w_);
    if (w && XLENGTH(w_) != n) error("internal error: one weight per row");
    int K = isNull(Y_) ? 0 : ncols(Y_);
    const double *Y = K ? REAL(Y_) : NULL;
    if (K && nrows(Y_) != n) error("internal error: one row of Y per row");
    int np = isNull(pairs_) ? 0 : ncols(pairs_);
    const int *pr = np ? INTEGER(pairs_) : NULL;
    /* term offsets: the intercept (no component) is covered by row 0 */
    int *toff = (int *) R_alloc(S.nt > 0 ? S.nt : 1, sizeof(int));
    int Q = 1;
    double work = 1.0;
    for (int t = 0; t < S.nt; t++) {
	const sterm *Tt = &S.t[t];
	if (Tt->nf == 0 && Tt->nn == 0) { toff[t] = -1; continue; }
	if (Tt->L < 0) error("internal error: term too large for the structured kernel");
	toff[t] = Q;
	Q += Tt->L * Tt->ns;
	work += Tt->ns + Tt->nf + Tt->nn;
    }
    int offY = Q;
    Q += K;
    const int **codes = (const int **) R_alloc(nk, sizeof(int *));
    for (int k = 0; k < nk; k++) {
	SEXP ck = VECTOR_ELT(codes_, k);
	if (XLENGTH(ck) != n) error("internal error: one code per row");
	xtx_check_codes(ck, G[k]);
	codes[k] = INTEGER(ck);
    }
    int *pa = (int *) R_alloc(np > 0 ? np : 1, sizeof(int)), *pb = (int *) R_alloc(np > 0 ? np : 1, sizeof(int));
    for (int q = 0; q < np; q++) {
	pa[q] = pr[2 * q];
	pb[q] = pr[2 * q + 1];
	if (pa[q] < 0 || pa[q] >= nk || pb[q] < 0 || pb[q] >= nk) error("internal error: pairs");
    }
    /* one accumulator: the nk matrices, then the np dense crosstabs */
    R_xlen_t *off = (R_xlen_t *) R_alloc((size_t) nk + np + 1, sizeof(R_xlen_t));
    off[0] = 0;
    for (int k = 0; k < nk; k++) off[k + 1] = off[k] + (R_xlen_t) Q * G[k];
    for (int q = 0; q < np; q++) off[nk + q + 1] = off[nk + q] + (R_xlen_t) G[pa[q]] * G[pb[q]];
    R_xlen_t len = off[nk + np];
    work = (double) n * (nk * (work + K) + np);
    int T = xtx_reduce_threads(asInteger(nthreads_), n, work, (double) len);
    SEXP out = PROTECT(allocVector(VECSXP, nk + np));
    double **mats = (double **) R_alloc(nk > 0 ? nk : 1, sizeof(double *));
    for (int k = 0; k < nk; k++) {
	SEXP m = SET_VECTOR_ELT(out, k, allocMatrix(REALSXP, Q, G[k]));
	mats[k] = REAL(m);
	memset(mats[k], 0, sizeof(double) * (size_t) Q * G[k]);
    }
    /* one thread and no crosstabs: the sums go to the matrices; otherwise
       each chunk sums into a contiguous buffer (chunk 0's: buf0), and the
       buffers are added in chunk order */
    int direct = T == 1 && np == 0;
    double *buf0 = NULL;
    if (!direct) {
	buf0 = (double *) R_alloc((size_t) (len > 0 ? len : 1), sizeof(double));
	memset(buf0, 0, sizeof(double) * (size_t) len);
    }
    double ***acc = (double ***) R_alloc(T, sizeof(double **));
    double ***pacc = (double ***) R_alloc(T, sizeof(double **));
    double **bufs = (double **) R_alloc(T, sizeof(double *));
    int ***F = (int ***) R_alloc(T, sizeof(int **));
    double ***N = (double ***) R_alloc(T, sizeof(double **));
    double **wb = (double **) R_alloc(T, sizeof(double *));
    int nta = S.nt > 0 ? S.nt : 1;
    for (int t = 0; t < T; t++) {
	acc[t] = (double **) R_alloc(nk > 0 ? nk : 1, sizeof(double *));
	pacc[t] = (double **) R_alloc(np > 0 ? np : 1, sizeof(double *));
	if (direct) {
	    bufs[t] = NULL;
	    for (int k = 0; k < nk; k++) acc[t][k] = mats[k];
	} else {
	    bufs[t] = t == 0 ? buf0 : (double *) R_alloc((size_t) (len > 0 ? len : 1), sizeof(double));
	    for (int k = 0; k < nk; k++) acc[t][k] = bufs[t] + off[k];
	    for (int q = 0; q < np; q++) pacc[t][q] = bufs[t] + off[nk + q];
	}
	F[t] = (int **) R_alloc(nta, sizeof(int *));
	N[t] = (double **) R_alloc(nta, sizeof(double *));
	for (int u = 0; u < S.nt; u++) {
	    F[t][u] = (int *) R_alloc(BLOCK, sizeof(int));
	    N[t][u] = (double *) R_alloc((size_t) BLOCK * (S.t[u].ns > 0 ? S.t[u].ns : 1), sizeof(double));
	}
	wb[t] = (double *) R_alloc(BLOCK, sizeof(double));
    }
    gsum_ctx C = { &S, toff, nk, Q, K, offY, np, (const int *const *) codes, w, Y, pa, pb, G };
    if (T == 1) {
	struct_gsum_rows(&C, 0, n, acc[0], pacc[0], F[0], N[0], wb[0]);
    } else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(T) schedule(static)
#endif
	for (int t = 0; t < T; t++) {
	    if (t > 0) xtx_zero(bufs[t], len);
	    struct_gsum_rows(&C, xtx_chunk_start(n, T, t, BLOCK), xtx_chunk_start(n, T, t + 1, BLOCK),
			     acc[t], pacc[t], F[t], N[t], wb[t]);
	}
	xtx_combine(bufs, T, len);
    }
    if (!direct)
	for (int k = 0; k < nk; k++) memcpy(mats[k], buf0 + off[k], sizeof(double) * (size_t) Q * G[k]);
    /* the dense crosstabs, compressed */
    for (int q = 0; q < np; q++) {
	int Ga = G[pa[q]], Gb = G[pb[q]];
	const double *d = buf0 + off[nk + q];
	R_xlen_t nnz = 0, tot = (R_xlen_t) Ga * Gb;
	for (R_xlen_t e = 0; e < tot; e++) nnz += d[e] != 0.0;
	if (nnz > INT_MAX) error("the crosstab has too many entries");
	SEXP pk = PROTECT(allocVector(INTSXP, (R_xlen_t) Ga + 1));
	SEXP ik = PROTECT(allocVector(INTSXP, nnz));
	SEXP xk = PROTECT(allocVector(REALSXP, nnz));
	int *pp = INTEGER(pk), *ii = INTEGER(ik);
	double *xx = REAL(xk);
	R_xlen_t e = 0;
	pp[0] = 0;
	for (int g = 0; g < Ga; g++) {
	    const double *dg = d + (size_t) g * Gb;
	    for (int h = 0; h < Gb; h++)
		if (dg[h] != 0.0) { ii[e] = h; xx[e++] = dg[h]; }
	    pp[g + 1] = (int) e;
	}
	SEXP ok = SET_VECTOR_ELT(out, nk + q, allocVector(VECSXP, 3));
	SET_VECTOR_ELT(ok, 0, pk);
	SET_VECTOR_ELT(ok, 1, ik);
	SET_VECTOR_ELT(ok, 2, xk);
	UNPROTECT(3);
    }
    UNPROTECT(2);
    return out;
}

/* Codes of the values of an FE variable (no missing values), and the rows
   per level (for the removal of singletons), in one pass over the rows:

   xtx_fe_code(v, maxrange): v integer, logical or double.  When v is
   integer valued with a range below maxrange, by direct indexing: codes
   1..G in increasing order of the values; returns list(codes, G, values,
   counts).  Otherwise, and for character vectors, by hashing:
   codes in order of first appearance; returns list(codes, G, values,
   counts, FALSE), the values in order of first appearance, which the
   caller sorts and passes to xtx_fe_recode().  NULL for a vector with
   missing values or of another type.  Doubles are equal when ==
   (0 and -0), as for unique() and match(); strings by their cached
   CHARSXP (the caller checks that the values are unique as strings,
   which they are unless the same string comes in two encodings).

   xtx_fe_recode(codes, rank): rank[codes], in place on a fresh copy. */

static inline uint64_t fe_mix(uint64_t h)
{
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccdULL;
    h ^= h >> 33;
    h *= 0xc4ceb9fe1a85ec53ULL;
    h ^= h >> 33;
    return h;
}

typedef struct { uint64_t key; int g; } fe_slot; /* g: level + 1, 0 empty */

static SEXP fe_code_hash(SEXP v)
{
    R_xlen_t n = XLENGTH(v);
    int str = TYPEOF(v) == STRSXP, isi = TYPEOF(v) == INTSXP || TYPEOF(v) == LGLSXP;
    const double *x = str || isi ? NULL : REAL_RO(v);
    const int *xi = isi ? INTEGER_RO(v) : NULL;
    const SEXP *sx = str ? STRING_PTR_RO(v) : NULL;
    /* open addressing on the keys (pointers to cached strings, integers,
       or the bits of doubles with -0 as 0: equal keys, equal values),
       load at most 1/2; a run of equal keys skips the table */
    size_t cap = 1024;
    fe_slot *slot = (fe_slot *) calloc(cap, sizeof(fe_slot));
    R_xlen_t ucap = 1024;
    R_xlen_t *first = (R_xlen_t *) malloc((size_t) ucap * sizeof(R_xlen_t));
    uint64_t *ukey = (uint64_t *) malloc((size_t) ucap * sizeof(uint64_t));
    int *cnt = (int *) malloc((size_t) ucap * sizeof(int));
#define FE_FREE() do { free(slot); free(first); free(ukey); free(cnt); } while (0)
    if (!slot || !first || !ukey || !cnt) { FE_FREE(); error("cannot allocate the codes"); }
    SEXP codes = PROTECT(allocVector(INTSXP, n));
    int *cd = INTEGER(codes);
    int G = 0, bad = 0, glast = -1;
    uint64_t klast = 0;
    for (R_xlen_t i = 0; i < n; i++) {
	uint64_t key;
	if (str) {
	    if (sx[i] == NA_STRING) { bad = 1; break; }
	    key = (uint64_t) (uintptr_t) sx[i];
	} else if (isi) {
	    if (xi[i] == NA_INTEGER) { bad = 1; break; }
	    key = (uint64_t) (int64_t) xi[i];
	} else {
	    double a = x[i];
	    if (ISNAN(a)) { bad = 1; break; }
	    if (a == 0.0) a = 0.0; /* -0 */
	    memcpy(&key, &a, sizeof(key));
	}
	int g;
	if (glast >= 0 && key == klast) {
	    g = glast;
	} else {
	    size_t mask = cap - 1, h = (size_t) fe_mix(key) & mask;
	    for (;;) {
		if (slot[h].g == 0) {
		    if (G == INT_MAX - 1) { bad = 1; break; }
		    if (G >= ucap) {
			ucap *= 2;
			R_xlen_t *f2 = (R_xlen_t *) realloc(first, (size_t) ucap * sizeof(R_xlen_t));
			if (f2) first = f2;
			uint64_t *k2 = (uint64_t *) realloc(ukey, (size_t) ucap * sizeof(uint64_t));
			if (k2) ukey = k2;
			int *c2 = (int *) realloc(cnt, (size_t) ucap * sizeof(int));
			if (c2) cnt = c2;
			if (!f2 || !k2 || !c2) { FE_FREE(); error("cannot allocate the codes"); }
		    }
		    first[G] = i;
		    ukey[G] = key;
		    cnt[G] = 0;
		    g = G++;
		    slot[h].key = key;
		    slot[h].g = g + 1;
		    if ((size_t) G * 2 > cap) {
			/* grow and rehash */
			size_t nc = cap * 2;
			fe_slot *ns = (fe_slot *) calloc(nc, sizeof(fe_slot));
			if (!ns) { FE_FREE(); error("cannot allocate the codes"); }
			for (int u = 0; u < G; u++) {
			    size_t hh = (size_t) fe_mix(ukey[u]) & (nc - 1);
			    while (ns[hh].g) hh = (hh + 1) & (nc - 1);
			    ns[hh].key = ukey[u];
			    ns[hh].g = u + 1;
			}
			free(slot);
			slot = ns;
			cap = nc;
		    }
		    break;
		}
		if (slot[h].key == key) { g = slot[h].g - 1; break; }
		h = (h + 1) & mask;
	    }
	    if (bad) break;
	    glast = g;
	    klast = key;
	}
	cd[i] = g + 1;
	cnt[g]++;
    }
    free(ukey);
    ukey = NULL;
    free(slot);
    slot = NULL;
    if (bad) {
	FE_FREE();
	UNPROTECT(1);
	return R_NilValue;
    }
    SEXP vals = PROTECT(allocVector(str ? STRSXP : isi ? INTSXP : REALSXP, G));
    SEXP counts = PROTECT(allocVector(INTSXP, G));
    for (int g = 0; g < G; g++) {
	if (str) SET_STRING_ELT(vals, g, sx[first[g]]);
	else if (isi) INTEGER(vals)[g] = xi[first[g]];
	else REAL(vals)[g] = x[first[g]];
	INTEGER(counts)[g] = cnt[g];
    }
    FE_FREE();
#undef FE_FREE
    SEXP out = PROTECT(allocVector(VECSXP, 5));
    SET_VECTOR_ELT(out, 0, codes);
    SET_VECTOR_ELT(out, 1, ScalarInteger(G));
    SET_VECTOR_ELT(out, 2, vals);
    SET_VECTOR_ELT(out, 3, counts);
    SET_VECTOR_ELT(out, 4, ScalarLogical(0));
    UNPROTECT(4);
    return out;
}

SEXP xtx_fe_code(SEXP v, SEXP maxrange_)
{
    R_xlen_t n = XLENGTH(v);
    double maxrange = asReal(maxrange_);
    if (n == 0) return R_NilValue;
    if (TYPEOF(v) == STRSXP) return fe_code_hash(v);
    double lo, hi;
    int isint = TYPEOF(v) == INTSXP || TYPEOF(v) == LGLSXP;
    if (isint) {
	const int *x = INTEGER_RO(v);
	int a = x[0], b = x[0];
	for (R_xlen_t i = 0; i < n; i++) {
	    if (x[i] == NA_INTEGER) return R_NilValue;
	    if (x[i] < a) a = x[i];
	    if (x[i] > b) b = x[i];
	}
	lo = a; hi = b;
    } else if (TYPEOF(v) == REALSXP) {
	const double *x = REAL_RO(v);
	lo = x[0]; hi = x[0];
	for (R_xlen_t i = 0; i < n; i++) {
	    double a = x[i];
	    if (ISNAN(a)) return R_NilValue;
	    if (!R_FINITE(a) || fabs(a) > 4.5e15 || a != (double)(long long) a) return fe_code_hash(v);
	    if (a < lo) lo = a;
	    if (a > hi) hi = a;
	}
    } else {
	return R_NilValue;
    }
    double range = hi - lo + 1;
    if (range > maxrange || range > INT_MAX / 2) return fe_code_hash(v);
    int R = (int) range;
    int *map = (int *) R_alloc((size_t) R, sizeof(int));
    memset(map, 0, sizeof(int) * (size_t) R);
    SEXP codes = PROTECT(allocVector(INTSXP, n));
    int *cd = INTEGER(codes);
    /* map[r]: the rows with value lo + r */
    if (isint) {
	const int *x = INTEGER_RO(v);
	int a = (int) lo;
	for (R_xlen_t i = 0; i < n; i++) { cd[i] = x[i] - a; map[cd[i]]++; }
    } else {
	const double *x = REAL_RO(v);
	for (R_xlen_t i = 0; i < n; i++) { cd[i] = (int) (x[i] - lo); map[cd[i]]++; }
    }
    int G = 0;
    for (int r = 0; r < R; r++) if (map[r]) G++;
    SEXP vals = PROTECT(allocVector(TYPEOF(v) == REALSXP ? REALSXP : INTSXP, G));
    SEXP counts = PROTECT(allocVector(INTSXP, G));
    int g = 0;
    for (int r = 0; r < R; r++) {
	if (!map[r]) continue;
	INTEGER(counts)[g] = map[r];
	if (TYPEOF(v) == REALSXP) REAL(vals)[g] = lo + r;
	else INTEGER(vals)[g] = (int) lo + r;
	map[r] = ++g;
    }
    for (R_xlen_t i = 0; i < n; i++) cd[i] = map[cd[i]];
    SEXP out = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(out, 0, codes);
    SET_VECTOR_ELT(out, 1, ScalarInteger(G));
    SET_VECTOR_ELT(out, 2, vals);
    SET_VECTOR_ELT(out, 3, counts);
    UNPROTECT(4);
    return out;
}

SEXP xtx_fe_recode(SEXP codes_, SEXP rank_)
{
    R_xlen_t n = XLENGTH(codes_);
    xtx_check_codes(codes_, LENGTH(rank_));
    const int *c = INTEGER_RO(codes_), *rk = INTEGER_RO(rank_);
    SEXP out = PROTECT(allocVector(INTSXP, n));
    int *o = INTEGER(out);
    for (R_xlen_t i = 0; i < n; i++) o[i] = rk[c[i] - 1];
    UNPROTECT(1);
    return out;
}

/*
 *  xtx: kernels for linear models with fixed effects (R/fe.R).
 *
 *  xtx_fe_gsum():       group sums of the columns of a matrix (codes 1..G);
 *  xtx_fe_crosstabs():  weighted crosstabs of a code vector with others, as
 *                      the compressed rows of G_a x G_b sparse matrices
 *                      with sorted, unique column indices;
 *  xtx_fe_components(): connected components of the bipartite graph whose
 *                      edges are the (a, b) pairs of a crosstab;
 *  xtx_fe_struct_gsum(): grouped sums of the design's structured moments,
 *                      the weights and other columns, for every FE at once,
 *                      and dense crosstabs of pairs of FE, in one pass;
 *  xtx_fe_code():       codes of an FE variable and rows per level;
 *  xtx_fe_schur():      the reduced system, formed (lower triangle);
 *  xtx_fe_ichol(), xtx_fe_pcg(): IC(0) and preconditioned conjugate
 *                      gradients for the reduced system, formed or
 *                      applied matrix-free.
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

static int uf_find(int *parent, int x)
{
    while (parent[x] != x) {
	parent[x] = parent[parent[x]];
	x = parent[x];
    }
    return x;
}

/* Components of the bipartite graph with edges (i[e], j[e]) (0-based, i in
   0..Ga-1, j in 0..Gb-1): the 1-based component of each b level, numbered
   in order of first appearance among the b levels. */
SEXP xtx_fe_components(SEXP i_, SEXP j_, SEXP Ga_, SEXP Gb_)
{
    R_xlen_t ne = XLENGTH(i_);
    int Ga = asInteger(Ga_), Gb = asInteger(Gb_);
    const int *ii = INTEGER(i_), *jj = INTEGER(j_);
    int N = Ga + Gb;
    int *parent = (int *) R_alloc((size_t)N, sizeof(int));
    for (int v = 0; v < N; v++) parent[v] = v;
    for (R_xlen_t e = 0; e < ne; e++) {
	int ra = uf_find(parent, ii[e]), rb = uf_find(parent, Ga + jj[e]);
	if (ra != rb) parent[ra] = rb;
    }
    int *label = (int *) R_alloc((size_t)N, sizeof(int));
    for (int v = 0; v < N; v++) label[v] = 0;
    SEXP out = PROTECT(allocVector(INTSXP, Gb));
    int *o = INTEGER(out), nc = 0;
    for (int v = 0; v < Gb; v++) {
	int r = uf_find(parent, Ga + v);
	if (!label[r]) label[r] = ++nc;
	o[v] = label[r];
    }
    UNPROTECT(1);
    return out;
}

/* ------------------------------------------------------------------------
   Structured grouped sums (see R/fe.R, xtx_fe_moments()).

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

/* ------------------------------------------------------------------------
   Preconditioned conjugate gradients for the reduced fixed-effect system.

   xtx_fe_ichol(): incomplete Cholesky factor with no fill, IC(0), of a
   symmetric matrix given by its lower triangle in compressed columns
   (p, i, x; row indices sorted within columns, diagonal first); with
   omega > 0, the modified factor MIC(0) (relaxed by omega): fill that the
   pattern drops is subtracted from the two diagonals, which keeps the row
   sums of L L' those of the matrix.  The
   reduced systems are (block diagonals of) weighted graph Laplacians with
   reference levels removed, which are M-matrices: IC(0) exists for them.
   Returns the factor's values (same pattern), or NULL on a breakdown
   (a non-positive pivot).

   xtx_fe_schur(): the reduced system R = A - N diag(d) N' as its lower
   triangle in compressed columns (rows sorted, diagonal first), A
   symmetric given by its lower triangle (the same layout), N a sparse
   matrix in compressed columns (rows sorted) whose rows map to the rows of
   R by `map` (0-based, increasing, -1 for a row dropped).  Gustavson's product
   by output column: column i of the lower triangle of N diag(d) N' sums
   over the columns j of N with a nonzero in row i the entries of column j
   at rows k >= i, which costs sum_j c_j (c_j + 1) / 2 for c_j the nonzeros
   of column j.  Gives up (returns NULL) as soon as the lower triangle has
   more than cap_nnz nonzeros or the sum of the squared column counts (the
   cost of IC(0) of R) exceeds cap_work: R is then applied without being
   formed (xtx_fe_pcg() with N).

   xtx_fe_pcg(): conjugate gradients for R X = B, one recursion per column
   of B (all columns advance together and share each product with R), R
   symmetric positive (semi)definite: R = A (lower triangle in compressed
   columns) or, when N is given, R = A - N diag(d) N' applied without
   forming it (two passes over N per product; map: the 0-based index in R
   of each row of N, or -1 for a row dropped).  Preconditioned by L L'
   (IC(0) factor of R or of a block diagonal of it) or by the diagonal of R
   when L is NULL.  X0: a starting value (n x m) or NULL (zero; also for
   the columns where X0 is worse than zero).  Stops
   when every column's normwise backward error ||B - R X|| / (||R|| ||X||
   + ||B||) is below tol, with max(diag R) for ||R||: the relative
   residual ||B - R X|| / ||B|| cannot go below rounding times ||R|| ||X||
   / ||B||, large when R is ill conditioned (on such systems the recursion
   then drifts and stalls).  Returns list(X, iterations (the largest over
   columns), max backward error, and NULL or, when the recursion of a
   column broke down (p'Rp <= 0: R singular, the FE collinear beyond the
   references removed), its direction p, close to a null vector of R).
   nthreads: on T > 1 threads (pcg_threads(): one for small systems) the
   products, sums and updates of the recursion are threaded (pcg_op,
   below: identical results for a given T); T = 1 runs the
   single-threaded code.
   ------------------------------------------------------------------------ */

SEXP xtx_fe_ichol(SEXP p_, SEXP i_, SEXP x_, SEXP n_, SEXP omega_)
{
    int n = asInteger(n_);
    double omega = asReal(omega_);
    const int *cp = INTEGER(p_), *ri = INTEGER(i_);
    SEXP out = PROTECT(duplicate(x_));
    double *L = REAL(out);
    int *pos = (int *) R_alloc((size_t) n + 1, sizeof(int));
    for (int r = 0; r < n; r++) pos[r] = -1;
    for (int k = 0; k < n; k++) {
	int d = cp[k];
	if (d >= cp[k + 1] || ri[d] != k) error("internal error: missing diagonal");
	double piv = L[d];
	if (!(piv > 0.0)) { UNPROTECT(1); return R_NilValue; }
	piv = sqrt(piv);
	L[d] = piv;
	for (int q = d + 1; q < cp[k + 1]; q++) L[q] /= piv;
	/* update the columns j > k that column k touches */
	for (int q = d + 1; q < cp[k + 1]; q++) {
	    int j = ri[q];
	    double ljk = L[q];
	    for (int t = cp[j]; t < cp[j + 1]; t++) pos[ri[t]] = t;
	    for (int s = q; s < cp[k + 1]; s++) {
		int t = pos[ri[s]];
		if (t >= 0) L[t] -= L[s] * ljk;
		else if (omega != 0.0) {
		    /* modified IC: the dropped fill goes to both diagonals */
		    double f = omega * L[s] * ljk;
		    L[cp[ri[s]]] -= f;
		    L[cp[j]] -= f;
		}
	    }
	    for (int t = cp[j]; t < cp[j + 1]; t++) pos[ri[t]] = -1;
	}
    }
    UNPROTECT(1);
    return out;
}

/* Solve L L' z = r in place for the nc interleaved columns of z
   (z[k * nc + c]); L lower, compressed columns, diagonal first. */
static void ic_solve(int n, int nc, const int *cp, const int *ri, const double *L, double *z)
{
    for (int k = 0; k < n; k++) {
	int d = cp[k];
	double *zk = z + (size_t) k * nc;
	for (int c = 0; c < nc; c++) zk[c] /= L[d];
	for (int q = d + 1; q < cp[k + 1]; q++) {
	    double *zr = z + (size_t) ri[q] * nc, l = L[q];
	    for (int c = 0; c < nc; c++) zr[c] -= l * zk[c];
	}
    }
    for (int k = n - 1; k >= 0; k--) {
	int d = cp[k];
	double *zk = z + (size_t) k * nc;
	for (int q = d + 1; q < cp[k + 1]; q++) {
	    const double *zr = z + (size_t) ri[q] * nc;
	    double l = L[q];
	    for (int c = 0; c < nc; c++) zk[c] -= l * zr[c];
	}
	for (int c = 0; c < nc; c++) zk[c] /= L[d];
    }
}

static int int_cmp(const void *u, const void *v)
{
    int a = *(const int *) u, b = *(const int *) v;
    return (a > b) - (a < b);
}

SEXP xtx_fe_schur(SEXP ap_, SEXP ai_, SEXP ax_, SEXP np_, SEXP ni_, SEXP nx_, SEXP map_,
		 SEXP d_, SEXP cap_nnz_, SEXP cap_work_)
{
    int m = length(ap_) - 1, G = length(np_) - 1;
    const int *ap = INTEGER(ap_), *ai = INTEGER(ai_), *np = INTEGER(np_), *ni = INTEGER(ni_);
    const int *map = INTEGER(map_);
    const double *ax = REAL(ax_), *nx = REAL(nx_), *d = REAL(d_);
    double cap_nnz = asReal(cap_nnz_), cap_work = asReal(cap_work_);
    int nnzN = np[G];
    /* the rows of N kept: for each entry, its column and its position in N
       (malloc: no pressure on R's garbage collector) */
    int *rp = (int *) R_alloc((size_t) m + 1, sizeof(int));
    int *fill = (int *) R_alloc((size_t) m + 1, sizeof(int));
    int *rj = (int *) malloc(((size_t) nnzN + 1) * sizeof(int));
    int *rq = (int *) malloc(((size_t) nnzN + 1) * sizeof(int));
    if (!rj || !rq) { free(rj); free(rq); error("cannot allocate the reduced system"); }
    for (int k = 0; k <= m; k++) rp[k] = 0;
    for (int q = 0; q < nnzN; q++) if (map[ni[q]] >= 0) rp[map[ni[q]] + 1]++;
    for (int k = 0; k < m; k++) rp[k + 1] += rp[k];
    memcpy(fill, rp, sizeof(int) * (size_t) m);
    for (int j = 0; j < G; j++)
	for (int q = np[j]; q < np[j + 1]; q++) {
	    int k = map[ni[q]];
	    if (k < 0) continue;
	    int t = fill[k]++;
	    rj[t] = j;
	    rq[t] = q;
	}
    int *mark = (int *) R_alloc((size_t) m + 1, sizeof(int));
    int *list = (int *) R_alloc((size_t) m + 1, sizeof(int));
    double *acc = (double *) R_alloc((size_t) m + 1, sizeof(double));
    for (int k = 0; k < m; k++) mark[k] = -1;
    double work = 0.0;
    int ok = 1;
    /* a probe first: the counts of 64 evenly spaced columns estimate the
       nonzeros and the work; give up at once when either is four times
       its cap (the exact checks below stay) */
    int ns = m < 64 ? m : 64;
    if (R_FINITE(cap_work) || R_FINITE(cap_nnz)) {
	double ew = 0.0, en = 0.0;
	for (int t = 0; t < ns; t++) {
	    int i = (int) (((double) t + 0.5) * m / ns), tag = -2 - t, cnt = 0;
	    for (int q = ap[i]; q < ap[i + 1]; q++)
		if (mark[ai[q]] != tag) { mark[ai[q]] = tag; cnt++; }
	    if (mark[i] != tag) { mark[i] = tag; cnt++; }
	    for (int e = rp[i]; e < rp[i + 1]; e++) {
		int j = rj[e];
		for (int q = rq[e]; q < np[j + 1]; q++) {
		    int k = map[ni[q]];
		    if (k >= 0 && mark[k] != tag) { mark[k] = tag; cnt++; }
		}
	    }
	    ew += (double) cnt * cnt;
	    en += cnt;
	}
	ew *= (double) m / ns;
	en *= (double) m / ns;
	if (ew > 4.0 * cap_work || en > 4.0 * cap_nnz) ok = 0;
	for (int k = 0; k < m; k++) mark[k] = -1;
    }
    size_t cap = (size_t) ap[m] + (size_t) m + 16, nnz = 0;
    int *oi = ok ? (int *) malloc(cap * sizeof(int)) : NULL;
    double *ox = ok ? (double *) malloc(cap * sizeof(double)) : NULL;
    if (ok && (!oi || !ox)) {
	free(oi); free(ox); free(rj); free(rq);
	error("cannot allocate the reduced system");
    }
    SEXP p_ = PROTECT(allocVector(INTSXP, m + 1));
    int *op = INTEGER(p_);
    op[0] = 0;
    for (int i = 0; i < m && ok; i++) {
	int cnt = 0;
	for (int q = ap[i]; q < ap[i + 1]; q++) {
	    int k = ai[q];
	    if (mark[k] != i) { mark[k] = i; acc[k] = 0.0; list[cnt++] = k; }
	    acc[k] += ax[q];
	}
	if (mark[i] != i) { mark[i] = i; acc[i] = 0.0; list[cnt++] = i; }
	for (int e = rp[i]; e < rp[i + 1]; e++) {
	    int j = rj[e], q0 = rq[e];
	    double wj = nx[q0] * d[j];
	    for (int q = q0; q < np[j + 1]; q++) {
		int k = map[ni[q]];
		if (k < 0) continue;
		if (mark[k] != i) { mark[k] = i; acc[k] = 0.0; list[cnt++] = k; }
		acc[k] -= wj * nx[q];
	    }
	}
	work += (double) cnt * cnt;
	if ((double) (nnz + cnt) > cap_nnz || work > cap_work) { ok = 0; break; }
	/* sorted rows: by a scan of the marks when the column is dense */
	if ((double) cnt * 16.0 > (double) (m - i)) {
	    int t = 0;
	    for (int k = i; k < m; k++) if (mark[k] == i) list[t++] = k;
	} else {
	    qsort(list, (size_t) cnt, sizeof(int), int_cmp);
	}
	if (nnz + cnt > cap) {
	    while (nnz + cnt > cap) cap *= 2;
	    int *oi2 = (int *) realloc(oi, cap * sizeof(int));
	    if (oi2) oi = oi2;
	    double *ox2 = (double *) realloc(ox, cap * sizeof(double));
	    if (ox2) ox = ox2;
	    if (!oi2 || !ox2) {
		free(oi); free(ox); free(rj); free(rq);
		error("cannot allocate the reduced system");
	    }
	}
	for (int t = 0; t < cnt; t++) {
	    oi[nnz] = list[t];
	    ox[nnz++] = acc[list[t]];
	}
	op[i + 1] = (int) nnz;
    }
    free(rj);
    free(rq);
    if (!ok) {
	free(oi);
	free(ox);
	UNPROTECT(1);
	return R_NilValue;
    }
    SEXP i_ = PROTECT(allocVector(INTSXP, nnz)), x_ = PROTECT(allocVector(REALSXP, nnz));
    memcpy(INTEGER(i_), oi, sizeof(int) * nnz);
    memcpy(REAL(x_), ox, sizeof(double) * nnz);
    free(oi);
    free(ox);
    SEXP out = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(out, 0, p_);
    SET_VECTOR_ELT(out, 1, i_);
    SET_VECTOR_ELT(out, 2, x_);
    SET_VECTOR_ELT(out, 3, ScalarReal(work));
    UNPROTECT(4);
    return out;
}

/* y = R x for nc interleaved columns: R = A (lower triangle) - N diag(d) N',
   the rows of N mapped by mp to 0..n-1, or to n (a row dropped: x has a
   zero row n, and row n of y is ignored).
   (Instances for one to four columns, where nc is a constant the compiler
   unrolls, and a generic one.) */
#define FE_OP_BODY(NC)							\
    memset(y, 0, sizeof(double) * ((size_t) n + 1) * (NC));		\
    for (int j = 0; j < n; j++) {					\
	const double *xj = x + (size_t) j * (NC);			\
	double *yj = y + (size_t) j * (NC);				\
	for (int q = ap[j]; q < ap[j + 1]; q++) {			\
	    int k = ai[q];						\
	    double a = ax[q];						\
	    double *yk = y + (size_t) k * (NC);				\
	    for (int c = 0; c < (NC); c++) yk[c] += a * xj[c];		\
	    if (k != j) {						\
		const double *xk = x + (size_t) k * (NC);		\
		for (int c = 0; c < (NC); c++) yj[c] += a * xk[c];	\
	    }								\
	}								\
    }									\
    if (np == NULL) return;						\
    for (int j = 0; j < G; j++) {					\
	double sj[4];							\
	double *sv = (NC) <= 4 ? sj : s;				\
	for (int c = 0; c < (NC); c++) sv[c] = 0.0;			\
	int q1 = np[j + 1];						\
	for (int q = np[j]; q < q1; q++) {				\
	    const double *xk = x + (size_t) mp[ni[q]] * (NC);		\
	    double v = nx[q];						\
	    for (int c = 0; c < (NC); c++) sv[c] += v * xk[c];		\
	}								\
	for (int c = 0; c < (NC); c++) sv[c] *= d[j];			\
	for (int q = np[j]; q < q1; q++) {				\
	    double *yk = y + (size_t) mp[ni[q]] * (NC);			\
	    double v = nx[q];						\
	    for (int c = 0; c < (NC); c++) yk[c] -= v * sv[c];		\
	}								\
    }

#define FE_OP_ARGS int n, const int *restrict ap, const int *restrict ai, \
	const double *restrict ax, int G, const int *restrict np,	\
	const int *restrict ni, const int *restrict mp,			\
	const double *restrict nx,						\
	const double *restrict d, const double *restrict x,		\
	double *restrict y, double *restrict s
static void fe_op1(FE_OP_ARGS) { FE_OP_BODY(1) }
static void fe_op2(FE_OP_ARGS) { FE_OP_BODY(2) }
static void fe_op3(FE_OP_ARGS) { FE_OP_BODY(3) }
static void fe_op4(FE_OP_ARGS) { FE_OP_BODY(4) }
static void fe_opn(int nc, FE_OP_ARGS) { FE_OP_BODY(nc) }

static void fe_op(int n, int nc, const int *ap, const int *ai, const double *ax,
		  int G, const int *np, const int *ni, const int *mp, const double *nx,
		  const double *d, const double *x, double *y, double *s)
{
    switch (nc) {
    case 1: fe_op1(n, ap, ai, ax, G, np, ni, mp, nx, d, x, y, s); break;
    case 2: fe_op2(n, ap, ai, ax, G, np, ni, mp, nx, d, x, y, s); break;
    case 3: fe_op3(n, ap, ai, ax, G, np, ni, mp, nx, d, x, y, s); break;
    case 4: fe_op4(n, ap, ai, ax, G, np, ni, mp, nx, d, x, y, s); break;
    default: fe_opn(nc, n, ap, ai, ax, G, np, ni, mp, nx, d, x, y, s);
    }
}

/* ------------------------------------------------------------------------
   xtx_fe_pcg() on T > 1 threads.

   Products y = R x: A as its full symmetric pattern in compressed rows
   (fp, fi, fx: row k lists (j, a) for each entry a at (k, j) or (j, k) of
   the triangle given), N by the rows of R in compressed rows (tp, tj, tx;
   the rows dropped left out).  First s = diag(d) N' x by columns of N,
   then each row of y as the sum over its row of A minus the sum over its
   row of N times s.  Every entry of s and of y is computed by one thread,
   in an order of operations fixed by the matrices alone: the products are
   the same for any T > 1 (the partition into T parts, balanced by
   nonzeros, decides only which thread computes what).

   Sums over the n rows (p'Rp, ||r||, ||x||, r'z): the rows split into T
   chunks by xtx_chunk_start() (boundaries from n and T only), one partial
   sum per chunk and column, the partial sums added in chunk order: the
   same for a given T.

   Updates of vectors are elementwise: the same for any T.

   The IC(0) solves stay on one thread (triangular solves).  No R API
   inside the parallel regions; everything is allocated before.
   ------------------------------------------------------------------------ */
typedef struct {
    int n, nc, G, T;
    const int *fp, *fi;
    const double *fx;
    const int *tp, *tj;		/* NULL without N */
    const double *tx;
    const int *np, *ni, *mp;
    const double *nx, *d;
    double *s;			/* G x nc, interleaved */
    const int *rb, *gb;		/* T + 1 bounds: rows of R, columns of N */
} pcg_op;

/* bounds b[0..T] (nondecreasing, b[0] = 0, b[T] = len) of T parts of
   `len` units balanced by the cost of unit u, 1 + cnt1[u + 1] - cnt1[u] +
   cnt2[u + 1] - cnt2[u] (cnt1, cnt2: cumulative counts, or NULL) */
static int *pcg_bounds(int len, int T, const int *cnt1, const int *cnt2)
{
    int *b = (int *) R_alloc((size_t) T + 1, sizeof(int));
    double total = (double) len;
    if (cnt1) total += (double) cnt1[len] - cnt1[0];
    if (cnt2) total += (double) cnt2[len] - cnt2[0];
    double cum = 0.0;
    int t = 1;
    b[0] = 0;
    for (int u = 0; u < len && t < T; u++) {
	while (t < T && cum >= total * t / T) b[t++] = u;
	cum += 1.0;
	if (cnt1) cum += (double) cnt1[u + 1] - cnt1[u];
	if (cnt2) cum += (double) cnt2[u + 1] - cnt2[u];
    }
    while (t < T) b[t++] = len;
    b[T] = len;
    return b;
}

static void pcg_op_build(pcg_op *P, int n, int nc, int T, const int *ap, const int *ai,
			 const double *ax, int G, const int *np, const int *ni, const int *mp,
			 const double *nx, const double *d)
{
    P->n = n;
    P->nc = nc;
    P->T = T;
    P->G = np ? G : 0;
    P->np = np;
    P->ni = ni;
    P->mp = mp;
    P->nx = nx;
    P->d = d;
    /* full A by rows */
    int *fp = (int *) R_alloc((size_t) n + 1, sizeof(int));
    int *fill = (int *) R_alloc((size_t) n + 1, sizeof(int));
    for (int k = 0; k <= n; k++) fp[k] = 0;
    for (int j = 0; j < n; j++)
	for (int q = ap[j]; q < ap[j + 1]; q++) {
	    int k = ai[q];
	    fp[k + 1]++;
	    if (k != j) fp[j + 1]++;
	}
    for (int k = 0; k < n; k++) fp[k + 1] += fp[k];
    int *fi = (int *) R_alloc((size_t) fp[n] + 1, sizeof(int));
    double *fx = (double *) R_alloc((size_t) fp[n] + 1, sizeof(double));
    memcpy(fill, fp, sizeof(int) * (size_t) n);
    for (int j = 0; j < n; j++)
	for (int q = ap[j]; q < ap[j + 1]; q++) {
	    int k = ai[q], t = fill[k]++;
	    fi[t] = j;
	    fx[t] = ax[q];
	    if (k != j) {
		t = fill[j]++;
		fi[t] = k;
		fx[t] = ax[q];
	    }
	}
    P->fp = fp;
    P->fi = fi;
    P->fx = fx;
    /* N by rows of R (the rows mapped to n, dropped, left out) */
    int *tp = NULL;
    P->tp = NULL;
    P->tj = NULL;
    P->tx = NULL;
    P->s = NULL;
    P->gb = NULL;
    if (P->G > 0) {
	tp = (int *) R_alloc((size_t) n + 1, sizeof(int));
	for (int k = 0; k <= n; k++) tp[k] = 0;
	int nnzN = np[G];
	for (int q = 0; q < nnzN; q++) {
	    int k = mp[ni[q]];
	    if (k < n) tp[k + 1]++;
	}
	for (int k = 0; k < n; k++) tp[k + 1] += tp[k];
	int *tj = (int *) R_alloc((size_t) tp[n] + 1, sizeof(int));
	double *tx = (double *) R_alloc((size_t) tp[n] + 1, sizeof(double));
	memcpy(fill, tp, sizeof(int) * (size_t) n);
	for (int j = 0; j < G; j++)
	    for (int q = np[j]; q < np[j + 1]; q++) {
		int k = mp[ni[q]];
		if (k >= n) continue;
		int t = fill[k]++;
		tj[t] = j;
		tx[t] = nx[q];
	    }
	P->tp = tp;
	P->tj = tj;
	P->tx = tx;
	P->s = (double *) R_alloc((size_t) G * nc + 1, sizeof(double));
	P->gb = pcg_bounds(G, T, np, NULL);
    }
    P->rb = pcg_bounds(n, T, fp, tp);
}

/* s[j] = d[j] (N' x)[j] for the columns j0..j1-1 of N (the rows mapped to
   n read the zero row n of x), the order of fe_op() */
static inline void pcg_gather(const pcg_op *P, int NC, const double *restrict x, int j0, int j1)
{
    const int *np = P->np, *ni = P->ni, *mp = P->mp;
    const double *nx = P->nx, *d = P->d;
    double *restrict s = P->s;
    for (int j = j0; j < j1; j++) {
	double *sj = s + (size_t) j * NC;
	for (int c = 0; c < NC; c++) sj[c] = 0.0;
	for (int q = np[j]; q < np[j + 1]; q++) {
	    const double *xk = x + (size_t) mp[ni[q]] * NC;
	    double v = nx[q];
	    for (int c = 0; c < NC; c++) sj[c] += v * xk[c];
	}
	for (int c = 0; c < NC; c++) sj[c] *= d[j];
    }
}

/* rows k0..k1-1 of y = A x - N s */
static inline void pcg_rows(const pcg_op *P, int NC, const double *restrict x, double *restrict y,
			    int k0, int k1)
{
    const int *fp = P->fp, *fi = P->fi, *tp = P->tp, *tj = P->tj;
    const double *fx = P->fx, *tx = P->tx;
    const double *s = P->s;
    for (int k = k0; k < k1; k++) {
	double *yk = y + (size_t) k * NC;
	for (int c = 0; c < NC; c++) yk[c] = 0.0;
	for (int t = fp[k]; t < fp[k + 1]; t++) {
	    const double *xj = x + (size_t) fi[t] * NC;
	    double a = fx[t];
	    for (int c = 0; c < NC; c++) yk[c] += a * xj[c];
	}
	if (tp)
	    for (int t = tp[k]; t < tp[k + 1]; t++) {
		const double *sj = s + (size_t) tj[t] * NC;
		double v = tx[t];
		for (int c = 0; c < NC; c++) yk[c] -= v * sj[c];
	    }
    }
}

/* (instances where the column count is a constant the compiler unrolls) */
static void pcg_gather_nc(const pcg_op *P, const double *x, int j0, int j1)
{
    switch (P->nc) {
    case 1: pcg_gather(P, 1, x, j0, j1); break;
    case 2: pcg_gather(P, 2, x, j0, j1); break;
    case 3: pcg_gather(P, 3, x, j0, j1); break;
    case 4: pcg_gather(P, 4, x, j0, j1); break;
    default: pcg_gather(P, P->nc, x, j0, j1);
    }
}

static void pcg_rows_nc(const pcg_op *P, const double *x, double *y, int k0, int k1)
{
    switch (P->nc) {
    case 1: pcg_rows(P, 1, x, y, k0, k1); break;
    case 2: pcg_rows(P, 2, x, y, k0, k1); break;
    case 3: pcg_rows(P, 3, x, y, k0, k1); break;
    case 4: pcg_rows(P, 4, x, y, k0, k1); break;
    default: pcg_rows(P, P->nc, x, y, k0, k1);
    }
}

/* y = R x on T threads (x: n + 1 rows, row n zero; row n of y set to
   zero).  The parts are loop iterations, not thread ids: the result is
   the same when OpenMP gives fewer threads than T. */
static void pcg_apply(const pcg_op *P, const double *x, double *y)
{
    int T = P->T, hasN = P->G > 0;
#ifdef _OPENMP
#pragma omp parallel num_threads(T)
#endif
    {
	if (hasN) {
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
	    for (int t = 0; t < T; t++) pcg_gather_nc(P, x, P->gb[t], P->gb[t + 1]);
	    /* (the implicit barrier: s complete before the rows) */
	}
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
	for (int t = 0; t < T; t++) pcg_rows_nc(P, x, y, P->rb[t], P->rb[t + 1]);
    }
    for (int c = 0; c < P->nc; c++) y[(size_t) P->n * P->nc + c] = 0.0;
}

/* out[c] = sum_k a[k, c] b[k, c] (interleaved, n rows) by T chunks, the
   partial sums (at acc + t * stride) added in chunk order */
static void pcg_dot(int n, int nc, int T, const double *a, const double *b, double *acc,
		    size_t stride, double *out)
{
#ifdef _OPENMP
#pragma omp parallel for num_threads(T) schedule(static)
#endif
    for (int t = 0; t < T; t++) {
	double *at = acc + (size_t) t * stride;
	for (int c = 0; c < nc; c++) at[c] = 0.0;
	R_xlen_t k0 = xtx_chunk_start(n, T, t, 1), k1 = xtx_chunk_start(n, T, t + 1, 1);
	for (R_xlen_t k = k0; k < k1; k++) {
	    const double *ak = a + (size_t) k * nc, *bk = b + (size_t) k * nc;
	    for (int c = 0; c < nc; c++) at[c] += ak[c] * bk[c];
	}
    }
    for (int c = 0; c < nc; c++) {
	double v = acc[c];
	for (int t = 1; t < T; t++) v += acc[(size_t) t * stride + c];
	out[c] = v;
    }
}

/* Threads for the products and sums of xtx_fe_pcg(): one for small systems
   (fewer than 2 x 4096 rows, or fewer than XTX_RED_MIN_WORK multiply-adds
   per thread and iteration).  requested < 0: -requested threads whatever
   the size (for the tests). */
static int pcg_threads(int requested, int n, int nc, double nnzA, double nnzN)
{
    if (requested == NA_INTEGER || n < 2) return 1;
    /* the full pattern of A must fit int offsets */
    if (2.0 * nnzA >= (double) INT_MAX) return 1;
    if (requested < 0) {
	int T = -requested;
	return T < n ? T : n;
    }
    double work = (double) nc * (2.0 * nnzA + 2.0 * nnzN + 8.0 * n);
    int T = xtx_worker_count(requested, n / 4096);
    while (T > 1 && work < XTX_RED_MIN_WORK * T) T--;
    return T;
}

SEXP xtx_fe_pcg(SEXP ap_, SEXP ai_, SEXP ax_, SEXP np_, SEXP ni_, SEXP nx_, SEXP map_, SEXP d_,
	       SEXP lp_, SEXP li_, SEXP lx_, SEXP B_, SEXP X0_, SEXP bmin_, SEXP tol_, SEXP maxit_,
	       SEXP nthreads_)
{
    int n = nrows(B_), nc = ncols(B_);
    const int *ap = INTEGER(ap_), *ai = INTEGER(ai_);
    const double *ax = REAL(ax_);
    int mf = !isNull(np_);
    int G = mf ? length(np_) - 1 : 0;
    const int *np = mf ? INTEGER(np_) : NULL, *ni = mf ? INTEGER(ni_) : NULL;
    const double *nx = mf ? REAL(nx_) : NULL, *d = mf ? REAL(d_) : NULL;
    int use_ic = !isNull(lx_);
    const int *lp = use_ic ? INTEGER(lp_) : NULL, *li = use_ic ? INTEGER(li_) : NULL;
    const double *lx = use_ic ? REAL(lx_) : NULL;
    double tol = asReal(tol_), bmin = asReal(bmin_);
    int maxit = asInteger(maxit_);
    const double *B = REAL(B_);
    SEXP X_ = PROTECT(allocMatrix(REALSXP, n, nc)), nullv_ = R_NilValue;
    int nprot = 1;
    double *X = REAL(X_);
    /* the rows of N: their index among the n levels, or n (dropped) */
    int *mp = NULL;
    if (mf) {
	int nr = length(map_);
	mp = (int *) R_alloc((size_t) nr + 1, sizeof(int));
	for (int k = 0; k < nr; k++) mp[k] = INTEGER(map_)[k] < 0 ? n : INTEGER(map_)[k];
    }
    /* vectors of n + 1 rows: row n (zero) stands for the rows dropped */
    size_t nn = (size_t) n * nc, len = ((size_t) n + 1) * nc;
    double *x = (double *) R_alloc(len, sizeof(double));
    double *r = (double *) R_alloc(len, sizeof(double));
    double *z = (double *) R_alloc(len, sizeof(double));
    double *pv = (double *) R_alloc(len, sizeof(double));
    double *Ap = (double *) R_alloc(len, sizeof(double));
    memset(x, 0, sizeof(double) * len);
    memset(r, 0, sizeof(double) * len);
    memset(z, 0, sizeof(double) * len);
    memset(pv, 0, sizeof(double) * len);
    double *s = (double *) R_alloc((size_t) nc + 1, sizeof(double));
    double *bn = (double *) R_alloc((size_t) nc + 1, sizeof(double));
    double *rz = (double *) R_alloc((size_t) nc + 1, sizeof(double));
    double *res = (double *) R_alloc((size_t) nc + 1, sizeof(double));
    double *alpha = (double *) R_alloc((size_t) nc + 1, sizeof(double));
    double *beta = (double *) R_alloc((size_t) nc + 1, sizeof(double));
    int *done = (int *) R_alloc((size_t) nc + 1, sizeof(int));
    int *its = (int *) R_alloc((size_t) nc + 1, sizeof(int));
    /* the diagonal of R, for the Jacobi preconditioner */
    double *dinv = (double *) R_alloc((size_t) n + 1, sizeof(double));
    double dmax = 0.0;
    for (int j = 0; j < n; j++) {
	dinv[j] = 0.0;
	for (int q = ap[j]; q < ap[j + 1]; q++) if (ai[q] == j) dinv[j] += ax[q];
    }
    if (mf)
	for (int j = 0; j < G; j++)
	    for (int q = np[j]; q < np[j + 1]; q++) dinv[mp[ni[q]]] -= nx[q] * nx[q] * d[j];
    for (int j = 0; j < n; j++) {
	if (dinv[j] > dmax) dmax = dinv[j];
	dinv[j] = dinv[j] > 0.0 ? 1.0 / dinv[j] : 1.0;
    }
    /* interleaved copies: x[k * nc + c] */
    for (int c = 0; c < nc; c++) {
	const double *b = B + (size_t) n * c;
	double t = 0.0;
	for (int k = 0; k < n; k++) t += b[k] * b[k];
	bn[c] = sqrt(t);
	for (int k = 0; k < n; k++) {
	    x[(size_t) k * nc + c] = isNull(X0_) ? 0.0 : REAL(X0_)[(size_t) n * c + k];
	    r[(size_t) k * nc + c] = b[k];
	}
	its[c] = 0;
    }
    /* T > 1: the products, sums and updates of the recursion on T threads
       (pcg_op); T = 1: the single-threaded code */
    int T = pcg_threads(asInteger(nthreads_), n, nc, (double) ap[n], mf ? (double) np[G] : 0.0);
    pcg_op P;
    memset(&P, 0, sizeof(P));
    double *acc = NULL;
    size_t stride = 0;
    if (T > 1) {
	pcg_op_build(&P, n, nc, T, ap, ai, ax, G, np, ni, mp, nx, d);
	/* per chunk: 3 nc partial sums, padded to 64 bytes (no false sharing) */
	stride = ((3 * (size_t) nc + 7) / 8) * 8;
	acc = (double *) R_alloc(stride * T, sizeof(double));
    }
    if (!isNull(X0_)) {
	if (T > 1) pcg_apply(&P, x, Ap);
	else fe_op(n, nc, ap, ai, ax, G, np, ni, mp, nx, d, x, Ap, s);
	for (size_t t = 0; t < nn; t++) r[t] -= Ap[t];
	/* a starting value worse than zero (residual above ||b||) is dropped */
	for (int c = 0; c < nc; c++) {
	    double t = 0.0;
	    for (int k = 0; k < n; k++) t += r[(size_t) k * nc + c] * r[(size_t) k * nc + c];
	    if (sqrt(t) > bn[c])
		for (int k = 0; k < n; k++) {
		    x[(size_t) k * nc + c] = 0.0;
		    r[(size_t) k * nc + c] = B[(size_t) n * c + k];
		}
	}
    }
    /* normwise backward error of column c: ||r|| / (max(diag R) ||x|| +
       max(||b||, bmin)), max(diag R) standing for ||R|| (a lower bound);
       bmin: a floor for right-hand sides at rounding level (such as those
       of a refinement step), relative to which tol is meaningless */
    #define BACKERR(c) do { \
	double rr_ = 0.0, xx_ = 0.0; \
	for (int k = 0; k < n; k++) { \
	    double rk_ = r[(size_t) k * nc + (c)], xk_ = x[(size_t) k * nc + (c)]; \
	    rr_ += rk_ * rk_; xx_ += xk_ * xk_; \
	} \
	res[c] = sqrt(rr_) / (dmax * sqrt(xx_) + (bn[c] > bmin ? bn[c] : bmin)); \
    } while (0)
    for (int c = 0; c < nc; c++) {
	if (bn[c] == 0.0) {
	    /* zero right-hand side: zero solution */
	    for (int k = 0; k < n; k++) x[(size_t) k * nc + c] = 0.0;
	    res[c] = 0.0;
	} else {
	    BACKERR(c);
	}
	done[c] = res[c] <= tol;
    }
    #define PRECOND() do { \
	if (use_ic) { memcpy(z, r, sizeof(double) * nn); ic_solve(n, nc, lp, li, lx, z); } \
	else for (int k = 0; k < n; k++) for (int c = 0; c < nc; c++) \
	    z[(size_t) k * nc + c] = dinv[k] * r[(size_t) k * nc + c]; \
    } while (0)
    PRECOND();
    memcpy(pv, z, sizeof(double) * nn);
    for (int c = 0; c < nc; c++) {
	rz[c] = 0.0;
	for (int k = 0; k < n; k++) rz[c] += r[(size_t) k * nc + c] * z[(size_t) k * nc + c];
    }
    int it = 0;
    if (T > 1) {
	double *pap = (double *) R_alloc((size_t) nc + 1, sizeof(double));
	double *rr = (double *) R_alloc((size_t) nc + 1, sizeof(double));
	double *xx = (double *) R_alloc((size_t) nc + 1, sizeof(double));
	double *rzn = (double *) R_alloc((size_t) nc + 1, sizeof(double));
	for (;;) {
	    int active = 0;
	    for (int c = 0; c < nc; c++) active += !done[c];
	    if (!active || it >= maxit) break;
	    it++;
	    pcg_apply(&P, pv, Ap);
	    pcg_dot(n, nc, T, pv, Ap, acc, stride, pap);
	    for (int c = 0; c < nc; c++) {
		alpha[c] = 0.0;
		if (done[c]) continue;
		if (!(pap[c] > 0.0)) {
		    /* p'Rp <= 0: as below */
		    done[c] = 1;
		    if (isNull(nullv_)) {
			nullv_ = PROTECT(allocVector(REALSXP, n));
			nprot++;
			for (int k = 0; k < n; k++) REAL(nullv_)[k] = pv[(size_t) k * nc + c];
		    }
		    continue;
		}
		alpha[c] = rz[c] / pap[c];
		its[c] = it;
	    }
	    /* x += alpha p, r -= alpha Rp with the sums of r^2 and x^2 (and,
	       Jacobi, z = diag(R)^-1 r with the sums of r z), by chunks */
#ifdef _OPENMP
#pragma omp parallel for num_threads(T) schedule(static)
#endif
	    for (int t = 0; t < T; t++) {
		double *a_rr = acc + (size_t) t * stride, *a_xx = a_rr + nc, *a_rz = a_xx + nc;
		for (int c = 0; c < 3 * nc; c++) a_rr[c] = 0.0;
		R_xlen_t k0 = xtx_chunk_start(n, T, t, 1), k1 = xtx_chunk_start(n, T, t + 1, 1);
		for (R_xlen_t k = k0; k < k1; k++) {
		    double dk = dinv[k];
		    for (int c = 0; c < nc; c++) {
			size_t i = (size_t) k * nc + c;
			x[i] += alpha[c] * pv[i];
			r[i] -= alpha[c] * Ap[i];
			a_rr[c] += r[i] * r[i];
			a_xx[c] += x[i] * x[i];
			if (!use_ic) {
			    z[i] = dk * r[i];
			    a_rz[c] += r[i] * z[i];
			}
		    }
		}
	    }
	    for (int c = 0; c < nc; c++) {
		double vr = acc[c], vx = acc[nc + c], vz = acc[2 * nc + c];
		for (int t = 1; t < T; t++) {
		    const double *at = acc + (size_t) t * stride;
		    vr += at[c];
		    vx += at[nc + c];
		    vz += at[2 * nc + c];
		}
		rr[c] = vr;
		xx[c] = vx;
		rzn[c] = vz;
	    }
	    for (int c = 0; c < nc; c++) {
		if (done[c]) continue;
		res[c] = sqrt(rr[c]) / (dmax * sqrt(xx[c]) + (bn[c] > bmin ? bn[c] : bmin));
		if (res[c] <= tol) done[c] = 1;
	    }
	    if (use_ic) {
		memcpy(z, r, sizeof(double) * nn);
		ic_solve(n, nc, lp, li, lx, z);
		pcg_dot(n, nc, T, r, z, acc, stride, rzn);
	    }
	    for (int c = 0; c < nc; c++) {
		beta[c] = 0.0;
		if (done[c]) continue;
		beta[c] = rzn[c] / rz[c];
		rz[c] = rzn[c];
	    }
#ifdef _OPENMP
#pragma omp parallel for num_threads(T) schedule(static)
#endif
	    for (int k = 0; k < n; k++)
		for (int c = 0; c < nc; c++) {
		    size_t i = (size_t) k * nc + c;
		    pv[i] = done[c] ? 0.0 : z[i] + beta[c] * pv[i];
		}
	    if ((it & 63) == 0) R_CheckUserInterrupt();
	}
    }
    else for (;;) {
	int active = 0;
	for (int c = 0; c < nc; c++) active += !done[c];
	if (!active || it >= maxit) break;
	it++;
	fe_op(n, nc, ap, ai, ax, G, np, ni, mp, nx, d, pv, Ap, s);
	for (int c = 0; c < nc; c++) {
	    alpha[c] = 0.0;
	    if (done[c]) continue;
	    double pap = 0.0;
	    for (int k = 0; k < n; k++) pap += pv[(size_t) k * nc + c] * Ap[(size_t) k * nc + c];
	    if (!(pap > 0.0)) {
		/* p'Rp <= 0: R is singular (rounding makes it slightly
		   indefinite) and p is close to a null vector */
		done[c] = 1;
		if (isNull(nullv_)) {
		    nullv_ = PROTECT(allocVector(REALSXP, n));
		    nprot++;
		    for (int k = 0; k < n; k++) REAL(nullv_)[k] = pv[(size_t) k * nc + c];
		}
		continue;
	    }
	    alpha[c] = rz[c] / pap;
	    its[c] = it;
	}
	for (int k = 0; k < n; k++)
	    for (int c = 0; c < nc; c++) {
		size_t t = (size_t) k * nc + c;
		x[t] += alpha[c] * pv[t];
		r[t] -= alpha[c] * Ap[t];
	    }
	for (int c = 0; c < nc; c++) {
	    if (done[c]) continue;
	    BACKERR(c);
	    if (res[c] <= tol) done[c] = 1;
	}
	PRECOND();
	for (int c = 0; c < nc; c++) {
	    beta[c] = 0.0;
	    if (done[c]) continue;
	    double t = 0.0;
	    for (int k = 0; k < n; k++) t += r[(size_t) k * nc + c] * z[(size_t) k * nc + c];
	    beta[c] = t / rz[c];
	    rz[c] = t;
	}
	for (int k = 0; k < n; k++)
	    for (int c = 0; c < nc; c++) {
		size_t t = (size_t) k * nc + c;
		pv[t] = done[c] ? 0.0 : z[t] + beta[c] * pv[t];
	    }
	if ((it & 63) == 0) R_CheckUserInterrupt();
    }
    #undef PRECOND
    #undef BACKERR
    int itmax = 0;
    double resmax = 0.0;
    for (int c = 0; c < nc; c++) {
	for (int k = 0; k < n; k++) X[(size_t) n * c + k] = x[(size_t) k * nc + c];
	if (its[c] > itmax) itmax = its[c];
	if (res[c] > resmax) resmax = res[c];
    }
    SEXP out = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(out, 0, X_);
    SET_VECTOR_ELT(out, 1, ScalarInteger(itmax));
    SET_VECTOR_ELT(out, 2, ScalarReal(resmax));
    SET_VECTOR_ELT(out, 3, nullv_);
    UNPROTECT(nprot + 1);
    return out;
}


/*
 *  design.matrix: kernels for linear models with fixed effects (R/fe.R).
 *
 *  dm_fe_gsum():       group sums of the columns of a matrix (codes 1..G);
 *  dm_fe_crosstab():   weighted crosstab of two code vectors, as the
 *                      compressed rows of a G_a x G_b sparse matrix with
 *                      sorted, unique column indices;
 *  dm_fe_components(): connected components of the bipartite graph whose
 *                      edges are the (a, b) pairs of a crosstab;
 *  dm_fe_struct_gsum(): grouped sums of the design's structured moments,
 *                      the weights and other columns, for every FE at once;
 *  dm_fe_code_int():   integer codes of an integer-valued vector;
 *  dm_fe_schur():      the reduced system, formed (lower triangle);
 *  dm_fe_ichol(), dm_fe_pcg(): IC(0) and preconditioned conjugate
 *                      gradients for the reduced system, formed or
 *                      applied matrix-free.
 *
 *  Copyright (C) 2026 Vincent Arel-Bundock.  GPL (>= 2).
 */

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <R.h>
#include <Rinternals.h>
#include "dm.h"
#include "structured.h"

/* G x k matrix of sums of the columns of W (n x k) within codes 1..G. */
SEXP dm_fe_gsum(SEXP codes_, SEXP G_, SEXP W)
{
    R_xlen_t n = XLENGTH(codes_);
    int G = asInteger(G_), k = ncols(W);
    if (nrows(W) != n) error("non-conformable arguments");
    const int *codes = INTEGER(codes_);
    const double *w = REAL(W);
    SEXP out = PROTECT(allocMatrix(REALSXP, G, k));
    double *o = REAL(out);
    memset(o, 0, sizeof(double) * (size_t)G * k);
    for (int j = 0; j < k; j++) {
	const double *wj = w + (R_xlen_t)j * n;
	double *oj = o + (R_xlen_t)j * G;
	for (R_xlen_t i = 0; i < n; i++) oj[codes[i] - 1] += wj[i];
    }
    UNPROTECT(1);
    return out;
}

/* Crosstab of codes a (1..Ga) and b (1..Gb) with weights w (NULL for 1):
   list(p, j, x), the compressed rows (by a) of the Ga x Gb matrix, j
   0-based and sorted within rows.  The rows are taken in order of a (as
   they come when a is sorted, as in panel data, otherwise after a counting
   sort); within a level of a, a marker array over b finds each pair's
   entry, and the level's entries are then sorted by b. */
typedef struct { int j; double x; } ct_entry;

static int ct_cmp(const void *u, const void *v)
{
    int a = ((const ct_entry *) u)->j, b = ((const ct_entry *) v)->j;
    return (a > b) - (a < b);
}

SEXP dm_fe_crosstab(SEXP a_, SEXP b_, SEXP w_, SEXP Ga_, SEXP Gb_)
{
    R_xlen_t n = XLENGTH(a_);
    int Ga = asInteger(Ga_), Gb = asInteger(Gb_);
    const int *a = INTEGER(a_), *b = INTEGER(b_);
    const double *w = isNull(w_) ? NULL : REAL(w_);
    int sorted = 1;
    for (R_xlen_t i = 1; i < n; i++) if (a[i] < a[i - 1]) { sorted = 0; break; }
    int *ord = NULL;
    if (!sorted) {
	R_xlen_t *cnt = (R_xlen_t *) R_alloc((size_t) Ga + 1, sizeof(R_xlen_t));
	ord = (int *) R_alloc((size_t) n, sizeof(int));
	memset(cnt, 0, sizeof(R_xlen_t) * ((size_t) Ga + 1));
	for (R_xlen_t i = 0; i < n; i++) cnt[a[i]]++;
	for (int g = 1; g <= Ga; g++) cnt[g] += cnt[g - 1];
	for (R_xlen_t i = 0; i < n; i++) ord[cnt[a[i] - 1]++] = (int) i;
    }
#define ROW(t) (sorted ? (t) : (R_xlen_t) ord[t])
    int *mark = (int *) R_alloc((size_t) Gb + 1, sizeof(int));
    for (int g = 0; g <= Gb; g++) mark[g] = 0;
    /* count the pairs */
    R_xlen_t nnz = 0;
    for (R_xlen_t t = 0; t < n; t++) {
	R_xlen_t i = ROW(t);
	if (mark[b[i]] != a[i]) { mark[b[i]] = a[i]; nnz++; }
    }
    SEXP p = PROTECT(allocVector(INTSXP, (R_xlen_t) Ga + 1));
    SEXP j = PROTECT(allocVector(INTSXP, nnz));
    SEXP x = PROTECT(allocVector(REALSXP, nnz));
    int *pp = INTEGER(p), *jj = INTEGER(j);
    double *xx = REAL(x);
    memset(pp, 0, sizeof(int) * ((size_t) Ga + 1));
    for (int g = 0; g <= Gb; g++) mark[g] = 0;
    R_xlen_t *pos = (R_xlen_t *) R_alloc((size_t) Gb + 1, sizeof(R_xlen_t));
    R_xlen_t q = 0, start = 0;
    ct_entry *tmp = NULL;
    size_t tmpn = 0;
    for (R_xlen_t t = 0; t <= n; t++) {
	R_xlen_t i = t < n ? ROW(t) : 0;
	int ai = t < n ? a[i] : -1;
	if (t > 0 && (t == n || ai != a[ROW(t - 1)])) {
	    /* close the level of the previous row: sort its entries by b */
	    R_xlen_t m = q - start;
	    int sortedb = 1;
	    for (R_xlen_t e = start + 1; e < q; e++) if (jj[e] < jj[e - 1]) { sortedb = 0; break; }
	    if (!sortedb) {
		if ((size_t) m > tmpn) {
		    tmpn = (size_t) m;
		    tmp = (ct_entry *) R_alloc(tmpn, sizeof(ct_entry));
		}
		for (R_xlen_t e = 0; e < m; e++) { tmp[e].j = jj[start + e]; tmp[e].x = xx[start + e]; }
		qsort(tmp, (size_t) m, sizeof(ct_entry), ct_cmp);
		for (R_xlen_t e = 0; e < m; e++) { jj[start + e] = tmp[e].j; xx[start + e] = tmp[e].x; }
	    }
	    pp[a[ROW(t - 1)]] = (int) m;
	    start = q;
	}
	if (t == n) break;
	int bi = b[i];
	if (mark[bi] != ai) {
	    mark[bi] = ai;
	    pos[bi] = q;
	    jj[q] = bi - 1;
	    xx[q] = 0.0;
	    q++;
	}
	xx[pos[bi]] += w ? w[i] : 1.0;
    }
#undef ROW
    for (int g = 1; g <= Ga; g++) pp[g] += pp[g - 1];
    SEXP out = PROTECT(allocVector(VECSXP, 3));
    SET_VECTOR_ELT(out, 0, p);
    SET_VECTOR_ELT(out, 1, j);
    SET_VECTOR_ELT(out, 2, x);
    UNPROTECT(4);
    return out;
}

/* Is code b nested in code a (each level of b in one level of a)? */
SEXP dm_fe_nested(SEXP a_, SEXP b_, SEXP Gb_)
{
    R_xlen_t n = XLENGTH(a_);
    int Gb = asInteger(Gb_);
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
SEXP dm_fe_lookup(SEXP values_, SEXP codes_)
{
    R_xlen_t n = XLENGTH(codes_), G = XLENGTH(values_);
    const double *v = REAL(values_);
    const int *c = INTEGER(codes_);
    SEXP out = PROTECT(allocVector(REALSXP, n));
    double *o = REAL(out);
    for (R_xlen_t i = 0; i < n; i++) {
	int k = c[i];
	o[i] = (k == NA_INTEGER || k < 1 || k > G) ? NA_REAL : v[k - 1];
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
SEXP dm_fe_components(SEXP i_, SEXP j_, SEXP Ga_, SEXP Gb_)
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
   Structured grouped sums (see R/fe.R, dm_fe_moments()).

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
   ------------------------------------------------------------------------ */


SEXP dm_fe_struct_gsum(SEXP setup, SEXP vars, SEXP w_, SEXP codes_, SEXP G_, SEXP Y_, SEXP terms_)
{
    sctx S;
    int use_terms = asLogical(terms_);
    R_xlen_t n;
    if (use_terms) {
	build_sctx(&S, setup, vars);
	n = S.n;
    } else {
	PROTECT(R_NilValue);
	S.nt = 0;
	S.t = NULL;
	n = XLENGTH(VECTOR_ELT(codes_, 0));
    }
    int nk = LENGTH(codes_);
    const int *G = INTEGER(G_);
    const double *w = isNull(w_) ? NULL : REAL(w_);
    if (w && XLENGTH(w_) != n) error("internal error: one weight per row");
    int K = isNull(Y_) ? 0 : ncols(Y_);
    const double *Y = K ? REAL(Y_) : NULL;
    if (K && nrows(Y_) != n) error("internal error: one row of Y per row");
    /* term offsets: the intercept (no component) is covered by row 0 */
    int *toff = (int *) R_alloc(S.nt > 0 ? S.nt : 1, sizeof(int));
    int Q = 1;
    for (int t = 0; t < S.nt; t++) {
	const sterm *T = &S.t[t];
	if (T->nf == 0 && T->nn == 0) { toff[t] = -1; continue; }
	if (T->L < 0) error("internal error: term too large for the structured kernel");
	toff[t] = Q;
	Q += T->L * T->ns;
    }
    int offY = Q;
    Q += K;
    const int **codes = (const int **) R_alloc(nk, sizeof(int *));
    double **acc = (double **) R_alloc(nk, sizeof(double *));
    SEXP out = PROTECT(allocVector(VECSXP, nk));
    for (int k = 0; k < nk; k++) {
	SEXP ck = VECTOR_ELT(codes_, k);
	if (XLENGTH(ck) != n) error("internal error: one code per row");
	codes[k] = INTEGER(ck);
	SEXP m = SET_VECTOR_ELT(out, k, allocMatrix(REALSXP, Q, G[k]));
	acc[k] = REAL(m);
	memset(acc[k], 0, sizeof(double) * (size_t) Q * G[k]);
    }
    int **F = (int **) R_alloc(S.nt > 0 ? S.nt : 1, sizeof(int *));
    double **N = (double **) R_alloc(S.nt > 0 ? S.nt : 1, sizeof(double *));
    for (int t = 0; t < S.nt; t++) {
	F[t] = (int *) R_alloc(BLOCK, sizeof(int));
	N[t] = (double *) R_alloc((size_t) BLOCK * (S.t[t].ns > 0 ? S.t[t].ns : 1), sizeof(double));
    }
    double *wb = (double *) R_alloc(BLOCK, sizeof(double));
    for (R_xlen_t i0 = 0; i0 < n; i0 += BLOCK) {
	int nb = (int) (n - i0 < BLOCK ? n - i0 : BLOCK);
	for (int t = 0; t < S.nt; t++)
	    if (toff[t] >= 0) term_block(&S.t[t], n, i0, nb, F[t], N[t]);
	for (int r = 0; r < nb; r++) wb[r] = w ? w[i0 + r] : 1.0;
	for (int k = 0; k < nk; k++) {
	    const int *ck = codes[k] + i0;
	    double *a = acc[k];
	    for (int r = 0; r < nb; r++) {
		double wi = wb[r];
		double *base = a + (size_t) Q * (ck[r] - 1);
		base[0] += wi;
		for (int t = 0; t < S.nt; t++) {
		    int o = toff[t];
		    if (o < 0) continue;
		    const sterm *T = &S.t[t];
		    int L = T->L, f = F[t][r];
		    const double *Nt = N[t] + r;
		    for (int u = 0; u < T->ns; u++)
			base[o + f + L * u] += wi * Nt[(size_t) u * BLOCK];
		}
		for (int j = 0; j < K; j++)
		    base[offY + j] += wi * Y[i0 + r + (R_xlen_t) j * n];
	    }
	}
    }
    UNPROTECT(2);
    return out;
}

/* Integer codes of an integer-valued vector (integer or double, no missing
   values) with a range below `maxrange`: codes 1..G in increasing order of
   the values, and the values.  NULL when the vector is not of that kind. */
SEXP dm_fe_code_int(SEXP v, SEXP maxrange_)
{
    R_xlen_t n = XLENGTH(v);
    double maxrange = asReal(maxrange_);
    if (n == 0) return R_NilValue;
    double lo, hi;
    if (TYPEOF(v) == INTSXP) {
	const int *x = INTEGER(v);
	int a = x[0], b = x[0];
	for (R_xlen_t i = 0; i < n; i++) {
	    if (x[i] == NA_INTEGER) return R_NilValue;
	    if (x[i] < a) a = x[i];
	    if (x[i] > b) b = x[i];
	}
	lo = a; hi = b;
    } else if (TYPEOF(v) == REALSXP) {
	const double *x = REAL(v);
	lo = x[0]; hi = x[0];
	for (R_xlen_t i = 0; i < n; i++) {
	    double a = x[i];
	    if (!R_FINITE(a) || a != (double)(long long) a) return R_NilValue;
	    if (a < lo) lo = a;
	    if (a > hi) hi = a;
	}
    } else {
	return R_NilValue;
    }
    double range = hi - lo + 1;
    if (range > maxrange || range > INT_MAX / 2) return R_NilValue;
    int R = (int) range;
    int *map = (int *) R_alloc((size_t) R, sizeof(int));
    memset(map, 0, sizeof(int) * (size_t) R);
    SEXP codes = PROTECT(allocVector(INTSXP, n));
    int *cd = INTEGER(codes);
    if (TYPEOF(v) == INTSXP) {
	const int *x = INTEGER(v);
	int a = (int) lo;
	for (R_xlen_t i = 0; i < n; i++) { cd[i] = x[i] - a; map[cd[i]] = 1; }
    } else {
	const double *x = REAL(v);
	for (R_xlen_t i = 0; i < n; i++) { cd[i] = (int) (x[i] - lo); map[cd[i]] = 1; }
    }
    int G = 0;
    for (int r = 0; r < R; r++) if (map[r]) map[r] = ++G;
    SEXP vals = PROTECT(allocVector(TYPEOF(v) == INTSXP ? INTSXP : REALSXP, G));
    for (int r = 0; r < R; r++) {
	if (!map[r]) continue;
	if (TYPEOF(v) == INTSXP) INTEGER(vals)[map[r] - 1] = (int) lo + r;
	else REAL(vals)[map[r] - 1] = lo + r;
    }
    for (R_xlen_t i = 0; i < n; i++) cd[i] = map[cd[i]];
    SEXP out = PROTECT(allocVector(VECSXP, 3));
    SET_VECTOR_ELT(out, 0, codes);
    SET_VECTOR_ELT(out, 1, ScalarInteger(G));
    SET_VECTOR_ELT(out, 2, vals);
    UNPROTECT(3);
    return out;
}

/* ------------------------------------------------------------------------
   Preconditioned conjugate gradients for the reduced fixed-effect system.

   dm_fe_ichol(): incomplete Cholesky factor with no fill, IC(0), of a
   symmetric matrix given by its lower triangle in compressed columns
   (p, i, x; row indices sorted within columns, diagonal first); with
   omega > 0, the modified factor MIC(0) (relaxed by omega): fill that the
   pattern drops is subtracted from the two diagonals, which keeps the row
   sums of L L' those of the matrix.  The
   reduced systems are (block diagonals of) weighted graph Laplacians with
   reference levels removed, which are M-matrices: IC(0) exists for them.
   Returns the factor's values (same pattern), or NULL on a breakdown
   (a non-positive pivot).

   dm_fe_schur(): the reduced system R = A - N diag(d) N' as its lower
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
   formed (dm_fe_pcg() with N).

   dm_fe_pcg(): conjugate gradients for R X = B, one recursion per column
   of B (all columns advance together and share each product with R), R
   symmetric positive (semi)definite: R = A (lower triangle in compressed
   columns) or, when N is given, R = A - N diag(d) N' applied without
   forming it (two passes over N per product; map: the 0-based index in R
   of each row of N, or -1 for a row dropped).  Preconditioned by L L'
   (IC(0) factor of R or of a block diagonal of it) or by the diagonal of R
   when L is NULL.  X0: a starting value (n x m) or NULL (zero).  Stops
   when every column's normwise backward error ||B - R X|| / (||R|| ||X||
   + ||B||) is below tol, with max(diag R) for ||R||: the relative
   residual ||B - R X|| / ||B|| cannot go below rounding times ||R|| ||X||
   / ||B||, large when R is ill conditioned (on such systems the recursion
   then drifts and stalls).  Returns list(X, iterations (the largest over
   columns), max backward error, and NULL or, when the recursion of a
   column broke down (p'Rp <= 0: R singular, the FE collinear beyond the
   references removed), its direction p, close to a null vector of R).
   ------------------------------------------------------------------------ */

SEXP dm_fe_ichol(SEXP p_, SEXP i_, SEXP x_, SEXP n_, SEXP omega_)
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

SEXP dm_fe_schur(SEXP ap_, SEXP ai_, SEXP ax_, SEXP np_, SEXP ni_, SEXP nx_, SEXP map_,
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

SEXP dm_fe_pcg(SEXP ap_, SEXP ai_, SEXP ax_, SEXP np_, SEXP ni_, SEXP nx_, SEXP map_, SEXP d_,
	       SEXP lp_, SEXP li_, SEXP lx_, SEXP B_, SEXP X0_, SEXP tol_, SEXP maxit_)
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
    double tol = asReal(tol_);
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
    if (!isNull(X0_)) {
	fe_op(n, nc, ap, ai, ax, G, np, ni, mp, nx, d, x, Ap, s);
	for (size_t t = 0; t < nn; t++) r[t] -= Ap[t];
    }
    /* normwise backward error of column c: ||r|| / (max(diag R) ||x|| +
       ||b||), max(diag R) standing for ||R|| (a lower bound) */
    #define BACKERR(c) do { \
	double rr_ = 0.0, xx_ = 0.0; \
	for (int k = 0; k < n; k++) { \
	    double rk_ = r[(size_t) k * nc + (c)], xk_ = x[(size_t) k * nc + (c)]; \
	    rr_ += rk_ * rk_; xx_ += xk_ * xk_; \
	} \
	res[c] = sqrt(rr_) / (dmax * sqrt(xx_) + bn[c]); \
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
    for (;;) {
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


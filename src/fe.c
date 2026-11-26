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
 *  dm_fe_ichol(), dm_fe_pcg(): IC(0) and preconditioned conjugate
 *                      gradients for the reduced system.
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

   dm_fe_pcg(): conjugate gradients for A X = B, A symmetric positive
   (semi)definite given in full compressed columns, one recursion per
   column of B, preconditioned by L L' (IC(0) factor) or by the diagonal of
   A when L is NULL.  Stops when every column's relative residual is below
   tol.  Returns list(X, iterations, max relative residual).
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

/* Solve L L' z = r in place (L lower, compressed columns, diagonal first). */
static void ic_solve(int n, const int *cp, const int *ri, const double *L, double *z)
{
    for (int k = 0; k < n; k++) {
	int d = cp[k];
	z[k] /= L[d];
	double zk = z[k];
	for (int q = d + 1; q < cp[k + 1]; q++) z[ri[q]] -= L[q] * zk;
    }
    for (int k = n - 1; k >= 0; k--) {
	int d = cp[k];
	double s = z[k];
	for (int q = d + 1; q < cp[k + 1]; q++) s -= L[q] * z[ri[q]];
	z[k] = s / L[d];
    }
}

SEXP dm_fe_pcg(SEXP ap_, SEXP ai_, SEXP ax_, SEXP lp_, SEXP li_, SEXP lx_, SEXP B_,
	       SEXP tol_, SEXP maxit_)
{
    int n = nrows(B_), m = ncols(B_);
    const int *ap = INTEGER(ap_), *ai = INTEGER(ai_);
    const double *ax = REAL(ax_);
    int use_ic = !isNull(lx_);
    const int *lp = use_ic ? INTEGER(lp_) : NULL, *li = use_ic ? INTEGER(li_) : NULL;
    const double *lx = use_ic ? REAL(lx_) : NULL;
    double tol = asReal(tol_);
    int maxit = asInteger(maxit_);
    const double *B = REAL(B_);
    SEXP X_ = PROTECT(allocMatrix(REALSXP, n, m));
    double *X = REAL(X_);
    memset(X, 0, sizeof(double) * (size_t) n * m);
    double *dinv = (double *) R_alloc((size_t) n + 1, sizeof(double));
    if (!use_ic) {
	for (int j = 0; j < n; j++) {
	    double d = 0.0;
	    for (int q = ap[j]; q < ap[j + 1]; q++) if (ai[q] == j) d += ax[q];
	    dinv[j] = d > 0.0 ? 1.0 / d : 1.0;
	}
    }
    double *r = (double *) R_alloc((size_t) n + 1, sizeof(double));
    double *z = (double *) R_alloc((size_t) n + 1, sizeof(double));
    double *pv = (double *) R_alloc((size_t) n + 1, sizeof(double));
    double *Ap = (double *) R_alloc((size_t) n + 1, sizeof(double));
    int itmax = 0;
    double resmax = 0.0;
    for (int c = 0; c < m; c++) {
	const double *b = B + (size_t) n * c;
	double *x = X + (size_t) n * c;
	double bn = 0.0;
	for (int j = 0; j < n; j++) { r[j] = b[j]; bn += b[j] * b[j]; }
	bn = sqrt(bn);
	if (bn == 0.0) continue;
	if (use_ic) { memcpy(z, r, sizeof(double) * n); ic_solve(n, lp, li, lx, z); }
	else for (int j = 0; j < n; j++) z[j] = dinv[j] * r[j];
	memcpy(pv, z, sizeof(double) * n);
	double rz = 0.0;
	for (int j = 0; j < n; j++) rz += r[j] * z[j];
	int it = 0;
	double res = 1.0;
	while (it < maxit) {
	    it++;
	    for (int j = 0; j < n; j++) Ap[j] = 0.0;
	    for (int j = 0; j < n; j++) {
		double pj = pv[j];
		for (int q = ap[j]; q < ap[j + 1]; q++) Ap[ai[q]] += ax[q] * pj;
	    }
	    double pap = 0.0;
	    for (int j = 0; j < n; j++) pap += pv[j] * Ap[j];
	    if (!(pap > 0.0)) break;
	    double alpha = rz / pap, rn = 0.0;
	    for (int j = 0; j < n; j++) {
		x[j] += alpha * pv[j];
		r[j] -= alpha * Ap[j];
		rn += r[j] * r[j];
	    }
	    res = sqrt(rn) / bn;
	    if (res <= tol) break;
	    if (use_ic) { memcpy(z, r, sizeof(double) * n); ic_solve(n, lp, li, lx, z); }
	    else for (int j = 0; j < n; j++) z[j] = dinv[j] * r[j];
	    double rz2 = 0.0;
	    for (int j = 0; j < n; j++) rz2 += r[j] * z[j];
	    double beta = rz2 / rz;
	    rz = rz2;
	    for (int j = 0; j < n; j++) pv[j] = z[j] + beta * pv[j];
	    if ((it & 63) == 0) R_CheckUserInterrupt();
	}
	if (it > itmax) itmax = it;
	if (res > resmax) resmax = res;
    }
    SEXP out = PROTECT(allocVector(VECSXP, 3));
    SET_VECTOR_ELT(out, 0, X_);
    SET_VECTOR_ELT(out, 1, ScalarInteger(itmax));
    SET_VECTOR_ELT(out, 2, ScalarReal(resmax));
    UNPROTECT(2);
    return out;
}

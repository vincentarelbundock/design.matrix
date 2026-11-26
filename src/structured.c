/*
 *  design.matrix: column-structured algebra on a plain design.
 *
 *  Row i of term t of the dense matrix is the Kronecker product of its
 *  components' rows.  Split the components into factors and numerics:
 *  the factors give one combined level F_t(i) (0 .. L_t - 1, the first
 *  factor varying fastest) and the numerics a vector N_t(i) of ns_t
 *  products (the first numeric varying fastest).  Then
 *
 *      X_t[i, colmap_t[u, k]] = N_t(i)[u] * C_t[F_t(i), k]
 *
 *  where C_t is the Kronecker product of the factors' contrast matrices.
 *  Products with the design reduce to sums over rows of N's within
 *  combined levels (moments), followed by small dense products with the
 *  C_t, done in R.  The moments are computed here, in blocks of rows so
 *  that each block's levels and products stay in cache.
 *
 *  Only for designs whose inputs are all finite, with no missing factor
 *  codes, and whose products cannot overflow (checked by dm_struct_ok):
 *  then the results agree with the dense computation up to rounding.
 *
 *  Copyright (C) 2026 Vincent Arel-Bundock.  GPL (>= 2).
 */

#include <string.h>
#include <float.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include "dm.h"

#define SMAXC 32
#define BLOCK 1024
#define SAFE (DBL_MAX / 16)

typedef struct {
    int L, ns, nf, nn, kc, count, off;
    const int *code[SMAXC];
    int adj[SMAXC], lstride[SMAXC], nrc[SMAXC], fnc[SMAXC], fstride[SMAXC];
    const double *cm[SMAXC];
    const double *num[SMAXC];
    int nnc[SMAXC], nstride[SMAXC];
} sterm;

typedef struct {
    R_xlen_t n;
    int p, nt;
    sterm *t;
} sctx;

/* Terms of a plain design (the intercept is a term with no component).
   Leaves one object PROTECTed (the coerced variables). */
static void build_sctx(sctx *S, SEXP setup, SEXP vars)
{
    SEXP factors = VECTOR_ELT(setup, 0), contr1 = VECTOR_ELT(setup, 1),
	contr2 = VECTOR_ELT(setup, 2);
    const int *nlevs = INTEGER(VECTOR_ELT(setup, 3)),
	*columns = INTEGER(VECTOR_ELT(setup, 4)),
	*count = INTEGER(VECTOR_ELT(setup, 5));
    int rhs_response = asInteger(VECTOR_ELT(setup, 6)),
	intercept = asInteger(VECTOR_ELT(setup, 7)),
	nVar = asInteger(VECTOR_ELT(setup, 10)),
	nterms = asInteger(VECTOR_ELT(setup, 11));
    const int *fac = nterms > 0 ? INTEGER(factors) : NULL;
    SEXP variable = PROTECT(dm_vars(setup, vars));
    S->n = nrows(VECTOR_ELT(vars, 0));
    S->p = LENGTH(VECTOR_ELT(setup, 8));
    S->t = (sterm *) R_alloc(nterms + 1, sizeof(sterm));
    S->nt = 0;
    int off = 0;
    if (intercept) {
	sterm *T = &S->t[S->nt++];
	memset(T, 0, sizeof(sterm));
	T->L = 1; T->ns = 1; T->kc = 1; T->count = 1; T->off = 0;
	off = 1;
    }
    for (int k = 0; k < nterms; k++) {
	if (k == rhs_response || count[k] <= 0) continue;
	sterm *T = &S->t[S->nt++];
	memset(T, 0, sizeof(sterm));
	T->off = off;
	T->count = count[k];
	off += count[k];
	int stride = 1, L = 1, ns = 1, kc = 1;
	for (int i = 0; i < nVar; i++) {
	    int fik = fac[i + k * nVar];
	    if (!fik || columns[i] == 0) continue;
	    if (T->nf + T->nn >= SMAXC) error("internal error: term with too many variables");
	    SEXP var_i = VECTOR_ELT(variable, i);
	    if (nlevs[i] > 0) {
		SEXP cm = VECTOR_ELT(fik == 1 ? contr1 : contr2, i);
		int c = T->nf++;
		T->code[c] = INTEGER(var_i);
		T->adj[c] = isLogical(var_i) ? 1 : 0;
		T->cm[c] = REAL(cm);
		T->nrc[c] = nrows(cm);
		T->fnc[c] = ncols(cm);
		T->lstride[c] = L;
		T->fstride[c] = stride;
		L *= nrows(cm);
		kc *= ncols(cm);
		stride *= ncols(cm);
	    } else {
		int c = T->nn++;
		T->num[c] = REAL(var_i);
		T->nnc[c] = ncols(var_i);
		T->nstride[c] = stride;
		ns *= ncols(var_i);
		stride *= ncols(var_i);
	    }
	}
	T->L = L; T->ns = ns; T->kc = kc;
    }
}

/* Combined levels and numeric products of rows i0 .. i0 + nb - 1. */
static void term_block(const sterm *T, R_xlen_t n, R_xlen_t i0, int nb, int *F, double *N)
{
    for (int r = 0; r < nb; r++) F[r] = 0;
    for (int c = 0; c < T->nf; c++) {
	const int *code = T->code[c] + i0;
	int a = T->adj[c] - 1, ls = T->lstride[c];
	for (int r = 0; r < nb; r++) F[r] += (code[r] + a) * ls;
    }
    int ns = T->ns;
    if (T->nn == 0) {
	for (int r = 0; r < nb; r++) N[r] = 1.0;
	return;
    }
    /* N[r * ns + u], built component by component, first fastest */
    int len = 1;
    for (int c = 0; c < T->nn; c++) {
	int nc = T->nnc[c];
	const double *v = T->num[c] + i0;
	if (c == 0) {
	    for (int j = 0; j < nc; j++)
		for (int r = 0; r < nb; r++) N[r * ns + j] = v[r + (R_xlen_t) j * n];
	} else {
	    for (int r = 0; r < nb; r++) {
		double *Nr = N + (size_t) r * ns;
		for (int j = nc - 1; j >= 0; j--) {
		    double x = v[r + (R_xlen_t) j * n];
		    for (int u = 0; u < len; u++) Nr[u + len * j] = Nr[u] * x;
		}
	    }
	}
	len *= nc;
    }
}

typedef struct {
    int *F;
    double *N;
} tblock;

static tblock *alloc_blocks(const sctx *S)
{
    tblock *B = (tblock *) R_alloc(S->nt > 0 ? S->nt : 1, sizeof(tblock));
    for (int t = 0; t < S->nt; t++) {
	B[t].F = (int *) R_alloc(BLOCK, sizeof(int));
	B[t].N = (double *) R_alloc((size_t) BLOCK * S->t[t].ns, sizeof(double));
    }
    return B;
}

/* Can the structured kernels handle this design?  All factor codes
   present, all numerics finite, products bounded, moments of a size that
   fits (`cap` doubles in all). */
SEXP dm_struct_ok(SEXP setup, SEXP vars, SEXP cap_)
{
    sctx S;
    build_sctx(&S, setup, vars);
    double cap = asReal(cap_), total = 0.0;
    int ok = 1;
    for (int t = 0; t < S.nt && ok; t++) {
	const sterm *T = &S.t[t];
	double bound = 1.0;
	for (int c = 0; c < T->nf && ok; c++) {
	    const int *code = T->code[c];
	    int lo = 1 - T->adj[c], hi = T->nrc[c] - T->adj[c];
	    for (R_xlen_t i = 0; i < S.n; i++)
		if (code[i] == NA_INTEGER || code[i] < lo || code[i] > hi) { ok = 0; break; }
	    double mx = 0.0;
	    for (R_xlen_t e = 0; e < (R_xlen_t) T->nrc[c] * T->fnc[c]; e++) {
		double x = T->cm[c][e];
		if (!R_FINITE(x)) { ok = 0; break; }
		if (fabs(x) > mx) mx = fabs(x);
	    }
	    bound *= mx > 1.0 ? mx : 1.0;
	}
	for (int c = 0; c < T->nn && ok; c++) {
	    double mx = 0.0;
	    const double *v = T->num[c];
	    for (R_xlen_t e = 0; e < S.n * T->nnc[c]; e++) {
		if (!R_FINITE(v[e])) { ok = 0; break; }
		if (fabs(v[e]) > mx) mx = fabs(v[e]);
	    }
	    bound *= mx > 1.0 ? mx : 1.0;
	}
	/* squares of entries appear in X'X */
	if (!(bound < sqrt(SAFE))) ok = 0;
	for (int s = 0; s <= t; s++) {
	    const sterm *U = &S.t[s];
	    total += (s == t) ? (double) T->L * T->ns * T->ns
		: (double) T->L * U->L * T->ns * U->ns;
	}
	total += (double) T->L * T->kc;
    }
    if (total > cap) ok = 0;
    UNPROTECT(1);
    return ScalarLogical(ok);
}

/* Term descriptions: list per term of list(off, L, C, colmap), where C is
   the L x kc combined contrast matrix and colmap the ns x kc matrix of
   0-based column indices within the term. */
SEXP dm_struct_info(SEXP setup, SEXP vars)
{
    sctx S;
    build_sctx(&S, setup, vars);
    SEXP out = PROTECT(allocVector(VECSXP, S.nt));
    for (int t = 0; t < S.nt; t++) {
	const sterm *T = &S.t[t];
	SEXP el = SET_VECTOR_ELT(out, t, allocVector(VECSXP, 3));
	SET_VECTOR_ELT(el, 0, ScalarInteger(T->off));
	SEXP C = SET_VECTOR_ELT(el, 1, allocMatrix(REALSXP, T->L, T->kc));
	SEXP M = SET_VECTOR_ELT(el, 2, allocMatrix(INTSXP, T->ns, T->kc));
	double *c = REAL(C);
	int *m = INTEGER(M);
	for (int k = 0; k < T->kc; k++) {
	    /* factor columns of combined column k (first factor fastest) */
	    int col[SMAXC], rem = k, fpart = 0;
	    for (int f = 0; f < T->nf; f++) {
		col[f] = rem % T->fnc[f];
		rem /= T->fnc[f];
		fpart += col[f] * T->fstride[f];
	    }
	    for (int a = 0; a < T->L; a++) {
		double v = 1.0;
		int r2 = a;
		for (int f = 0; f < T->nf; f++) {
		    int lev = r2 % T->nrc[f];
		    r2 /= T->nrc[f];
		    v *= T->cm[f][lev + (R_xlen_t) col[f] * T->nrc[f]];
		}
		c[a + (R_xlen_t) k * T->L] = v;
	    }
	    for (int u = 0; u < T->ns; u++) {
		int r3 = u, npart = 0;
		for (int q = 0; q < T->nn; q++) {
		    npart += (r3 % T->nnc[q]) * T->nstride[q];
		    r3 /= T->nnc[q];
		}
		m[u + k * T->ns] = npart + fpart;
	    }
	}
    }
    UNPROTECT(2);
    return out;
}

/* Moments for X'WX: list over pairs (s <= t, in order t = 0.., s = 0..t)
   of arrays: s < t: [Ls, Lt, ns, nt] sums of w N_s[u] N_t[v] by level
   pair; s == t: [Ls, ns, ns]. */
SEXP dm_struct_gram(SEXP setup, SEXP vars, SEXP w_)
{
    sctx S;
    build_sctx(&S, setup, vars);
    const double *w = isNull(w_) ? NULL : REAL(w_);
    if (w && XLENGTH(w_) != S.n) error("internal error: one weight per row");
    int np = S.nt * (S.nt + 1) / 2;
    SEXP out = PROTECT(allocVector(VECSXP, np));
    double **acc = (double **) R_alloc(np > 0 ? np : 1, sizeof(double *));
    int q = 0;
    for (int t = 0; t < S.nt; t++)
	for (int s = 0; s <= t; s++, q++) {
	    const sterm *A = &S.t[s], *B = &S.t[t];
	    R_xlen_t len = (s == t) ? (R_xlen_t) A->L * A->ns * A->ns
		: (R_xlen_t) A->L * B->L * A->ns * B->ns;
	    SEXP m = SET_VECTOR_ELT(out, q, allocVector(REALSXP, len));
	    acc[q] = REAL(m);
	    memset(acc[q], 0, sizeof(double) * (size_t) len);
	}
    tblock *Bk = alloc_blocks(&S);
    double *wN = (double *) R_alloc(BLOCK, sizeof(double));
    for (R_xlen_t i0 = 0; i0 < S.n; i0 += BLOCK) {
	int nb = (int) (S.n - i0 < BLOCK ? S.n - i0 : BLOCK);
	for (int t = 0; t < S.nt; t++) term_block(&S.t[t], S.n, i0, nb, Bk[t].F, Bk[t].N);
	q = 0;
	for (int t = 0; t < S.nt; t++) {
	    const sterm *B = &S.t[t];
	    const int *Ft = Bk[t].F;
	    const double *Nt = Bk[t].N;
	    int nt = B->ns;
	    for (int s = 0; s <= t; s++, q++) {
		const sterm *A = &S.t[s];
		const int *Fs = Bk[s].F;
		const double *Ns = Bk[s].N;
		int ns = A->ns, Ls = A->L;
		double *M = acc[q];
		if (s == t) {
		    for (int u = 0; u < ns; u++) {
			if (w) for (int r = 0; r < nb; r++) wN[r] = w[i0 + r] * Ns[r * ns + u];
			else for (int r = 0; r < nb; r++) wN[r] = Ns[r * ns + u];
			for (int v = 0; v < ns; v++) {
			    double *Muv = M + (R_xlen_t) Ls * (u + ns * v);
			    if (Ls == 1) {
				double sum = 0.0;
				for (int r = 0; r < nb; r++) sum += wN[r] * Ns[r * ns + v];
				Muv[0] += sum;
			    } else {
				for (int r = 0; r < nb; r++) Muv[Fs[r]] += wN[r] * Ns[r * ns + v];
			    }
			}
		    }
		} else {
		    int Lt = B->L;
		    R_xlen_t LL = (R_xlen_t) Ls * Lt;
		    for (int u = 0; u < ns; u++) {
			if (w) for (int r = 0; r < nb; r++) wN[r] = w[i0 + r] * Ns[r * ns + u];
			else for (int r = 0; r < nb; r++) wN[r] = Ns[r * ns + u];
			for (int v = 0; v < nt; v++) {
			    double *Muv = M + LL * (u + ns * v);
			    if (LL == 1) {
				double sum = 0.0;
				for (int r = 0; r < nb; r++) sum += wN[r] * Nt[r * nt + v];
				Muv[0] += sum;
			    } else if (Lt == 1) {
				for (int r = 0; r < nb; r++) Muv[Fs[r]] += wN[r] * Nt[r * nt + v];
			    } else if (Ls == 1) {
				for (int r = 0; r < nb; r++) Muv[Ft[r]] += wN[r] * Nt[r * nt + v];
			    } else {
				for (int r = 0; r < nb; r++)
				    Muv[Fs[r] + (R_xlen_t) Ls * Ft[r]] += wN[r] * Nt[r * nt + v];
			    }
			}
		    }
		}
	    }
	}
    }
    UNPROTECT(2);
    return out;
}

/* Moments for X'Y (Y n x K): list per term of [L, ns, K] sums of
   N[u] Y[, k] by level. */
SEXP dm_struct_xty(SEXP setup, SEXP vars, SEXP Y_)
{
    sctx S;
    build_sctx(&S, setup, vars);
    if (nrows(Y_) != S.n) error("non-conformable arguments");
    int K = ncols(Y_);
    const double *Y = REAL(Y_);
    SEXP out = PROTECT(allocVector(VECSXP, S.nt));
    double **acc = (double **) R_alloc(S.nt > 0 ? S.nt : 1, sizeof(double *));
    for (int t = 0; t < S.nt; t++) {
	R_xlen_t len = (R_xlen_t) S.t[t].L * S.t[t].ns * K;
	acc[t] = REAL(SET_VECTOR_ELT(out, t, allocVector(REALSXP, len)));
	memset(acc[t], 0, sizeof(double) * (size_t) len);
    }
    tblock *Bk = alloc_blocks(&S);
    for (R_xlen_t i0 = 0; i0 < S.n; i0 += BLOCK) {
	int nb = (int) (S.n - i0 < BLOCK ? S.n - i0 : BLOCK);
	for (int t = 0; t < S.nt; t++) {
	    const sterm *T = &S.t[t];
	    term_block(T, S.n, i0, nb, Bk[t].F, Bk[t].N);
	    const int *F = Bk[t].F;
	    const double *N = Bk[t].N;
	    int ns = T->ns, L = T->L;
	    for (int k = 0; k < K; k++) {
		const double *y = Y + i0 + (R_xlen_t) k * S.n;
		for (int u = 0; u < ns; u++) {
		    double *M = acc[t] + (R_xlen_t) L * (u + (R_xlen_t) ns * k);
		    if (L == 1) {
			double sum = 0.0;
			for (int r = 0; r < nb; r++) sum += N[r * ns + u] * y[r];
			M[0] += sum;
		    } else {
			for (int r = 0; r < nb; r++) M[F[r]] += N[r * ns + u] * y[r];
		    }
		}
	    }
	}
    }
    UNPROTECT(2);
    return out;
}

/* X B from per-term tables P_t [L, ns, K] (P_t = C_t B_t, built in R):
   row i, column k is sum over t, u of N_t(i)[u] P_t[F_t(i), u, k]. */
SEXP dm_struct_mult(SEXP setup, SEXP vars, SEXP P, SEXP K_)
{
    sctx S;
    build_sctx(&S, setup, vars);
    int K = asInteger(K_);
    if (LENGTH(P) != S.nt) error("internal error: one table per term");
    SEXP out = PROTECT(allocMatrix(REALSXP, S.n, K));
    double *o = REAL(out);
    memset(o, 0, sizeof(double) * (size_t) S.n * K);
    tblock *Bk = alloc_blocks(&S);
    for (R_xlen_t i0 = 0; i0 < S.n; i0 += BLOCK) {
	int nb = (int) (S.n - i0 < BLOCK ? S.n - i0 : BLOCK);
	for (int t = 0; t < S.nt; t++) {
	    const sterm *T = &S.t[t];
	    term_block(T, S.n, i0, nb, Bk[t].F, Bk[t].N);
	    const int *F = Bk[t].F;
	    const double *N = Bk[t].N;
	    const double *p = REAL(VECTOR_ELT(P, t));
	    int ns = T->ns, L = T->L;
	    for (int k = 0; k < K; k++) {
		double *ok = o + i0 + (R_xlen_t) k * S.n;
		for (int u = 0; u < ns; u++) {
		    const double *pu = p + (R_xlen_t) L * (u + (R_xlen_t) ns * k);
		    if (L == 1) {
			double c = pu[0];
			for (int r = 0; r < nb; r++) ok[r] += N[r * ns + u] * c;
		    } else {
			for (int r = 0; r < nb; r++) ok[r] += N[r * ns + u] * pu[F[r]];
		    }
		}
	    }
	}
    }
    UNPROTECT(2);
    return out;
}

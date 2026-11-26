/*
 *  design.matrix: helpers for the fixed-effect prototype (R/fe.R).
 *
 *  dm_fe_gsum():       group sums of the columns of a matrix (codes 1..G);
 *  dm_fe_crosstab():   weighted crosstab of two code vectors, as the
 *                      compressed rows of a G_a x G_b sparse matrix with
 *                      sorted, unique column indices;
 *  dm_fe_components(): connected components of the bipartite graph whose
 *                      edges are the (a, b) pairs of a crosstab.
 *
 *  Copyright (C) 2026 Vincent Arel-Bundock.  GPL (>= 2).
 */

#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include "dm.h"

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
   0-based and sorted within rows.  Two stable counting sorts (by b, then by
   a) order the pairs; equal pairs are then adjacent. */
SEXP dm_fe_crosstab(SEXP a_, SEXP b_, SEXP w_, SEXP Ga_, SEXP Gb_)
{
    R_xlen_t n = XLENGTH(a_);
    int Ga = asInteger(Ga_), Gb = asInteger(Gb_);
    const int *a = INTEGER(a_), *b = INTEGER(b_);
    const double *w = isNull(w_) ? NULL : REAL(w_);
    R_xlen_t *cnt = (R_xlen_t *) R_alloc((size_t)(Gb > Ga ? Gb : Ga) + 1, sizeof(R_xlen_t));
    int *ord1 = (int *) R_alloc((size_t)n, sizeof(int));
    int *ord2 = (int *) R_alloc((size_t)n, sizeof(int));
    /* sort by b */
    memset(cnt, 0, sizeof(R_xlen_t) * ((size_t)Gb + 1));
    for (R_xlen_t i = 0; i < n; i++) cnt[b[i]]++;
    for (int g = 1; g <= Gb; g++) cnt[g] += cnt[g - 1];
    for (R_xlen_t i = 0; i < n; i++) ord1[cnt[b[i] - 1]++] = (int) i;
    /* stable sort by a */
    memset(cnt, 0, sizeof(R_xlen_t) * ((size_t)Ga + 1));
    for (R_xlen_t i = 0; i < n; i++) cnt[a[i]]++;
    for (int g = 1; g <= Ga; g++) cnt[g] += cnt[g - 1];
    for (R_xlen_t t = 0; t < n; t++) {
	int i = ord1[t];
	ord2[cnt[a[i] - 1]++] = i;
    }
    /* count unique pairs */
    R_xlen_t nnz = 0;
    for (R_xlen_t t = 0; t < n; t++) {
	int i = ord2[t];
	if (t == 0 || a[i] != a[ord2[t - 1]] || b[i] != b[ord2[t - 1]]) nnz++;
    }
    SEXP p = PROTECT(allocVector(INTSXP, (R_xlen_t)Ga + 1));
    SEXP j = PROTECT(allocVector(INTSXP, nnz));
    SEXP x = PROTECT(allocVector(REALSXP, nnz));
    int *pp = INTEGER(p), *jj = INTEGER(j);
    double *xx = REAL(x);
    memset(pp, 0, sizeof(int) * ((size_t)Ga + 1));
    R_xlen_t q = -1;
    for (R_xlen_t t = 0; t < n; t++) {
	int i = ord2[t];
	if (t == 0 || a[i] != a[ord2[t - 1]] || b[i] != b[ord2[t - 1]]) {
	    q++;
	    jj[q] = b[i] - 1;
	    xx[q] = 0.0;
	    pp[a[i]]++;
	}
	xx[q] += w ? w[i] : 1.0;
    }
    for (int g = 1; g <= Ga; g++) pp[g] += pp[g - 1];
    SEXP out = PROTECT(allocVector(VECSXP, 3));
    SET_VECTOR_ELT(out, 0, p);
    SET_VECTOR_ELT(out, 1, j);
    SET_VECTOR_ELT(out, 2, x);
    UNPROTECT(4);
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

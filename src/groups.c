/*
 *  design.matrix: groups of identical rows.
 *
 *  Two rows of a design (or Jacobian) give the same row-wise results when
 *  every input the row's entries depend on has the same bits: the codes of
 *  the factors and the values of the numeric variables used by the terms,
 *  the per-row scales, and an optional group code.  dm_row_groups() hashes
 *  those bits, confirms every match by comparing the bits exactly, and
 *  returns, for each row, the 1-based id of its pattern (in order of first
 *  appearance) and the 1-based index of each pattern's first row.
 *
 *  Copyright (C) 2026 Vincent Arel-Bundock.  GPL (>= 2).
 */

#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <R.h>
#include <Rinternals.h>
#include "dm.h"

typedef struct {
    int nint, ndbl;
    const int **ic;       /* integer key columns */
    const double **dc;    /* double key columns */
} keycols;

static void add_int(keycols *K, const int *x) { K->ic[K->nint++] = x; }
static void add_dbl(keycols *K, const double *x) { K->dc[K->ndbl++] = x; }

/* Key columns of one design: variables used by a term with columns. */
static void design_keys(keycols *K, SEXP setup, SEXP variable, R_xlen_t n)
{
    SEXP factors = VECTOR_ELT(setup, DM_SETUP_FACTORS);
    const int *nlevs = INTEGER(VECTOR_ELT(setup, DM_SETUP_NLEVS)),
	*columns = INTEGER(VECTOR_ELT(setup, DM_SETUP_COLUMNS)),
	*count = INTEGER(VECTOR_ELT(setup, DM_SETUP_COUNT));
    int rhs_response = asInteger(VECTOR_ELT(setup, DM_SETUP_RHS_RESPONSE)),
	nVar = asInteger(VECTOR_ELT(setup, DM_SETUP_NVAR)),
	nterms = asInteger(VECTOR_ELT(setup, DM_SETUP_NTERMS));
    const int *fac = nterms > 0 ? INTEGER(factors) : NULL;
    for (int i = 0; i < nVar; i++) {
	int used = 0;
	for (int k = 0; k < nterms; k++)
	    if (k != rhs_response && count[k] > 0 && fac[i + k * nVar]) used = 1;
	if (!used || columns[i] == 0) continue;
	SEXP v = VECTOR_ELT(variable, i);
	if (nlevs[i] > 0) {
	    add_int(K, INTEGER(v));
	} else {
	    const double *x = REAL(v);
	    int nc = ncols(v);
	    for (int j = 0; j < nc; j++) add_dbl(K, x + (R_xlen_t)j * n);
	}
    }
}

static int design_nkeys(SEXP setup, SEXP variable)
{
    const int *columns = INTEGER(VECTOR_ELT(setup, DM_SETUP_COLUMNS));
    int nVar = asInteger(VECTOR_ELT(setup, DM_SETUP_NVAR)), total = 0;
    for (int i = 0; i < nVar; i++)
	total += columns[i] > 0 ? ncols(VECTOR_ELT(variable, i)) : 0;
    return total;
}

static R_INLINE uint64_t mix(uint64_t h, uint64_t x)
{
    /* splitmix64 finaliser on the combination */
    h ^= x + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    h ^= h >> 30; h *= 0xbf58476d1ce4e5b9ULL;
    h ^= h >> 27; h *= 0x94d049bb133111ebULL;
    h ^= h >> 31;
    return h;
}

static R_INLINE uint64_t dbl_bits(double x)
{
    uint64_t u;
    memcpy(&u, &x, sizeof(u));
    return u;
}

static R_INLINE uint64_t row_hash(const keycols *K, R_xlen_t i)
{
    uint64_t h = 0x1234567ULL;
    for (int j = 0; j < K->nint; j++) h = mix(h, (uint64_t)(uint32_t) K->ic[j][i]);
    for (int j = 0; j < K->ndbl; j++) h = mix(h, dbl_bits(K->dc[j][i]));
    return h;
}

static R_INLINE int same_row(const keycols *K, R_xlen_t a, R_xlen_t b)
{
    for (int j = 0; j < K->nint; j++)
	if (K->ic[j][a] != K->ic[j][b]) return 0;
    for (int j = 0; j < K->ndbl; j++)
	if (dbl_bits(K->dc[j][a]) != dbl_bits(K->dc[j][b])) return 0;
    return 1;
}

/* .Call(C_dm_row_groups, setup1, vars1, s1, setup2, vars2, s2, group, max_ratio)
   Returns NULL as soon as more than max_ratio * n distinct rows are seen. */
SEXP dm_row_groups(SEXP setup1, SEXP vars1, SEXP s1, SEXP setup2, SEXP vars2,
		   SEXP s2, SEXP group, SEXP max_ratio_)
{
    SEXP var1 = PROTECT(dm_vars(setup1, vars1));
    SEXP var2 = PROTECT(isNull(setup2) ? R_NilValue : dm_vars(setup2, vars2));
    R_xlen_t n = nrows(VECTOR_ELT(vars1, 0));
    if (n > INT_MAX) error("too many rows");
    double max_ratio = asReal(max_ratio_);

    int cap = design_nkeys(setup1, var1) + 4 +
	(isNull(setup2) ? 0 : design_nkeys(setup2, var2));
    keycols K;
    K.nint = K.ndbl = 0;
    K.ic = (const int **) R_alloc(cap, sizeof(int *));
    K.dc = (const double **) R_alloc(cap, sizeof(double *));
    design_keys(&K, setup1, var1, n);
    if (!isNull(setup2)) design_keys(&K, setup2, var2, n);
    if (!isNull(s1) && XLENGTH(s1) > 1) add_dbl(&K, REAL(s1));
    if (!isNull(s2) && XLENGTH(s2) > 1) add_dbl(&K, REAL(s2));
    if (!isNull(group)) add_int(&K, INTEGER(group));

    /* Each slot holds a representative row + 1. Its group id is already
       stored in map. Grow only when a new group would exceed half load. */
    R_xlen_t size = 128, mask = size - 1;
    PROTECT_INDEX table_index;
    SEXP table;
    PROTECT_WITH_INDEX(table = allocVector(INTSXP, size), &table_index);
    int *slot = INTEGER(table);
    memset(slot, 0, sizeof(int) * size);

    SEXP map = PROTECT(allocVector(INTSXP, n));
    int *mp = INTEGER(map);
    int U = 0;
    double limit = max_ratio * (double) n;

    for (R_xlen_t i = 0; i < n; i++) {
	uint64_t h = row_hash(&K, i);
	R_xlen_t s = (R_xlen_t)(h & mask);
	while (slot[s] && !same_row(&K, slot[s] - 1, i))
	    s = (s + 1) & mask;
	if (slot[s]) {
	    mp[i] = mp[slot[s] - 1];
	    continue;
	}
	if (++U > limit) {
	    UNPROTECT(4);
	    return R_NilValue;
	}
	if (U > size / 2) {
	    R_xlen_t old_size = size;
	    size *= 2;
	    for (int step = 1; step < 3 && size < 2 * n; step++) size *= 2;
	    mask = size - 1;
	    SEXP larger = PROTECT(allocVector(INTSXP, size));
	    int *next = INTEGER(larger);
	    memset(next, 0, sizeof(int) * size);
	    for (R_xlen_t j = 0; j < old_size; j++) {
		int row = slot[j];
		if (!row) continue;
		R_xlen_t k = (R_xlen_t)(row_hash(&K, row - 1) & mask);
		while (next[k]) k = (k + 1) & mask;
		next[k] = row;
	    }
	    REPROTECT(table = larger, table_index);
	    UNPROTECT(1);
	    slot = next;
	    s = (R_xlen_t)(h & mask);
	    while (slot[s]) s = (s + 1) & mask;
	}
	slot[s] = (int) i + 1;
	mp[i] = U;
    }

    SEXP fst = PROTECT(allocVector(INTSXP, U));
    int *first = INTEGER(fst);
    /* New ids occur in increasing order. Stop once the last one is found. */
    int seen = 0;
    for (R_xlen_t i = 0; seen < U; i++) {
	if (mp[i] > seen) first[seen++] = (int) i + 1;
    }
    SEXP out = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(out, 0, map);
    SET_VECTOR_ELT(out, 1, fst);
    UNPROTECT(6);
    return out;
}

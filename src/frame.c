#include <string.h>
#include <stdint.h>
#include <R.h>
#include <Rinternals.h>
#include "dm.h"

/* Are all the codes of an integer vector NA or in 1..k? (No allocation.) */
SEXP dm_codes_ok(SEXP x, SEXP k_)
{
    if (TYPEOF(x) != INTSXP) return ScalarLogical(0);
    const int *v = INTEGER_RO(x);
    int k = asInteger(k_);
    R_xlen_t n = XLENGTH(x);
    for (R_xlen_t i = 0; i < n; i++) {
	int c = v[i];
	if (c != NA_INTEGER && (c < 1 || c > k)) return ScalarLogical(0);
    }
    return ScalarLogical(1);
}

/* Are the values of an integer vector all different, none of them NA?
   TRUE or FALSE, decided with a bit set over [min, max]; NA when that
   range is too wide for one (more than 64 values per element). */
SEXP dm_int_unique(SEXP x)
{
    if (TYPEOF(x) != INTSXP) return ScalarLogical(NA_LOGICAL);
    const int *v = INTEGER_RO(x);
    R_xlen_t n = XLENGTH(x);
    if (n == 0) return ScalarLogical(1);
    int lo = v[0], hi = v[0];
    for (R_xlen_t i = 0; i < n; i++) {
	int c = v[i];
	if (c == NA_INTEGER) return ScalarLogical(0);
	if (c < lo) lo = c;
	if (c > hi) hi = c;
    }
    double range = (double) hi - (double) lo + 1.0;
    if (range > 64.0 * (double) n + 1024.0) return ScalarLogical(NA_LOGICAL);
    size_t nbytes = (size_t) (range / 8.0) + 1;
    unsigned char *seen = (unsigned char *) R_alloc(nbytes, sizeof(unsigned char));
    memset(seen, 0, nbytes);
    for (R_xlen_t i = 0; i < n; i++) {
	size_t k = (size_t) ((double) v[i] - (double) lo);
	unsigned char bit = (unsigned char) (1u << (k & 7));
	if (seen[k >> 3] & bit) return ScalarLogical(0);
	seen[k >> 3] |= bit;
    }
    return ScalarLogical(1);
}

/* Are the strings of a character vector all different, none of them NA?
   R keeps one copy of each string, so that two strings with the same
   encoding mark are equal only if they are the same object: a set of
   pointers decides. Strings with different marks can be equal after
   translation: NA (not decided here) when the marks are not all the same. */
SEXP dm_str_unique(SEXP x)
{
    if (TYPEOF(x) != STRSXP) return ScalarLogical(NA_LOGICAL);
    R_xlen_t n = XLENGTH(x);
    if (n == 0) return ScalarLogical(1);
    cetype_t ce = getCharCE(STRING_ELT(x, 0));
    size_t size = 16;
    while (size < 2 * (size_t) n) size *= 2;
    size_t mask = size - 1;
    SEXP *seen = (SEXP *) R_alloc(size, sizeof(SEXP));
    memset(seen, 0, size * sizeof(SEXP));
    for (R_xlen_t i = 0; i < n; i++) {
	SEXP s = STRING_ELT(x, i);
	if (s == NA_STRING) return ScalarLogical(0);
	if (getCharCE(s) != ce) return ScalarLogical(NA_LOGICAL);
	size_t h = (size_t) (((uintptr_t) s >> 4) * (uintptr_t) 0x9E3779B97F4A7C15ULL) & mask;
	while (seen[h] != NULL) {
	    if (seen[h] == s) return ScalarLogical(0);
	    h = (h + 1) & mask;
	}
	seen[h] = s;
    }
    return ScalarLogical(1);
}

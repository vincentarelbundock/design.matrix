/* Fused one-pass kernels of the in-memory fitting code (Kata
   design.matrix#bg3b): each replaces several R vector operations, and the
   n-length temporaries between them, by one pass that allocates only its
   result.  The arithmetic is that of the R expressions they replace, in
   the same order (one operation per statement, so that no compiler
   contracts two of them into one rounding), and sums accumulate in long
   double as sum() does: results are those of the R code. */

#include <R.h>
#include <Rinternals.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <Rmath.h>
#include "xtx.h"

/* R's sum() of doubles: long double accumulation, +-Inf beyond DBL_MAX */
static double xtx_ld_value(long double s)
{
    if (s > DBL_MAX) return R_PosInf;
    if (s < -DBL_MAX) return R_NegInf;
    return (double) s;
}

/* xtx_c_ifelse() in memory: `no` (recycled) with `yes` (recycled) where test
   is TRUE and NA where test is NA.  test logical; yes, no doubles. */
SEXP xtx_ifelse(SEXP test_, SEXP yes_, SEXP no_)
{
    R_xlen_t n = XLENGTH(test_), ny = XLENGTH(yes_), nn = XLENGTH(no_);
    const int *t = LOGICAL_RO(test_);
    const double *y = REAL_RO(yes_), *no = REAL_RO(no_);
    SEXP out = PROTECT(allocVector(REALSXP, n));
    double *o = REAL(out);
    if (n > 0 && (ny == 0 || nn == 0)) error("internal error: empty `yes` or `no`");
    for (R_xlen_t i = 0; i < n; i++) {
	int ti = t[i];
	if (ti == NA_LOGICAL) o[i] = NA_REAL;
	else if (ti) o[i] = ny == 1 ? y[0] : y[i % ny];
	else o[i] = nn == 1 ? no[0] : no[i % nn];
    }
    UNPROTECT(1);
    return out;
}

/* xtx_c_notnan() in memory: TRUE where x is present, FALSE where it is NA
   or NaN (as !is.na(x) gives). */
SEXP xtx_notnan(SEXP x_)
{
    R_xlen_t n = XLENGTH(x_);
    SEXP out = PROTECT(allocVector(LGLSXP, n));
    int *o = LOGICAL(out);
    if (TYPEOF(x_) == REALSXP) {
	const double *x = REAL_RO(x_);
	for (R_xlen_t i = 0; i < n; i++) o[i] = !ISNAN(x[i]);
    } else if (TYPEOF(x_) == INTSXP || TYPEOF(x_) == LGLSXP) {
	const int *x = TYPEOF(x_) == INTSXP ? INTEGER_RO(x_) : LOGICAL_RO(x_);
	for (R_xlen_t i = 0; i < n; i++) o[i] = x[i] != NA_INTEGER;
    } else {
	error("internal error: xtx_notnan() of a non-numeric vector");
    }
    UNPROTECT(1);
    return out;
}

/* The moments of a response in memory (xtx_c_moments()), as its R code sums them:
   c(sum(wz z), sum(w), sum(wz), sum(w (z - ybar)^2)) with wz = w z (z
   without weights; then sum(w) is n) and ybar = sum(wz) / sum(w). */
SEXP xtx_fe_moments(SEXP z_, SEXP w_)
{
    R_xlen_t n = XLENGTH(z_);
    const double *z = REAL_RO(z_);
    const double *w = isNull(w_) ? NULL : REAL_RO(w_);
    if (w && XLENGTH(w_) != n) error("internal error: one weight per row");
    long double s1 = 0.0, s2 = 0.0, s3 = 0.0, s4 = 0.0;
    for (R_xlen_t i = 0; i < n; i++) {
	double wz = w ? w[i] * z[i] : z[i];
	double q = wz * z[i];
	s1 += q;
	if (w) s2 += w[i];
	s3 += wz;
    }
    double sw = w ? xtx_ld_value(s2) : (double) n;
    double ybar = xtx_ld_value(s3) / sw;
    for (R_xlen_t i = 0; i < n; i++) {
	double d = z[i] - ybar;
	double d2 = d * d;
	double t = w ? w[i] * d2 : d2;
	s4 += t;
    }
    SEXP out = PROTECT(allocVector(REALSXP, 4));
    REAL(out)[0] = xtx_ld_value(s1);
    REAL(out)[1] = w ? xtx_ld_value(s2) : 0.0;
    REAL(out)[2] = xtx_ld_value(s3);
    REAL(out)[3] = xtx_ld_value(s4);
    UNPROTECT(1);
    return out;
}

/* The variables of a design with ALTREP numeric vectors replaced by plain
   copies (same values and attributes).  stats sets the contrasts of a
   factor on a column shared with the model frame, which makes it a
   wrapper; the kernels take writable pointers to the variables, and R
   then copies a wrapper's data at every call.  One copy here instead. */
SEXP xtx_unwrap_vars(SEXP vars)
{
    if (TYPEOF(vars) != VECSXP) return vars;
    R_xlen_t nv = XLENGTH(vars);
    SEXP out = vars;
    int nprot = 0;
    for (R_xlen_t j = 0; j < nv; j++) {
	SEXP v = VECTOR_ELT(vars, j);
	int t = TYPEOF(v);
	if (!ALTREP(v) || (t != INTSXP && t != REALSXP && t != LGLSXP)) continue;
	if (out == vars) {
	    out = PROTECT(shallow_duplicate(vars));
	    nprot++;
	}
	R_xlen_t n = XLENGTH(v);
	SEXP c = PROTECT(allocVector(t, n));
	size_t sz = t == REALSXP ? sizeof(double) : sizeof(int);
	void *dst = t == REALSXP ? (void *) REAL(c) : (t == INTSXP ? (void *) INTEGER(c) : (void *) LOGICAL(c));
	if (n > 0) memcpy(dst, DATAPTR_RO(v), sz * (size_t) n);
	DUPLICATE_ATTRIB(c, v);
	SET_VECTOR_ELT(out, j, c);
	UNPROTECT(1);
    }
    UNPROTECT(nprot);
    return out;
}

/* Does factor x (integer codes, L levels) leave a level unused, as
   length(unique(x[!is.na(x)])) < length(levels(x)) tells?  NA when a code
   is out of range (then the R code decides). */
SEXP xtx_levels_unused(SEXP x_, SEXP L_)
{
    R_xlen_t n = XLENGTH(x_);
    int L = asInteger(L_);
    const int *x = INTEGER_RO(x_);
    if (L == NA_INTEGER || L < 0) return ScalarLogical(NA_LOGICAL);
    char *seen = (char *) R_alloc((size_t) L + 1, 1);
    memset(seen, 0, (size_t) L + 1);
    int distinct = 0;
    for (R_xlen_t i = 0; i < n; i++) {
	int c = x[i];
	if (c == NA_INTEGER) continue;
	if (c < 1 || c > L) return ScalarLogical(NA_LOGICAL);
	if (!seen[c]) { seen[c] = 1; distinct++; }
    }
    return ScalarLogical(distinct < L);
}

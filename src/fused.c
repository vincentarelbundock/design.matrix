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
#include "dm.h"

/* R's sum() of doubles: long double accumulation, +-Inf beyond DBL_MAX */
static double dm_ld_value(long double s)
{
    if (s > DBL_MAX) return R_PosInf;
    if (s < -DBL_MAX) return R_NegInf;
    return (double) s;
}

/* dmc_ifelse() in memory: `no` (recycled) with `yes` (recycled) where test
   is TRUE and NA where test is NA.  test logical; yes, no doubles. */
SEXP dm_ifelse(SEXP test_, SEXP yes_, SEXP no_)
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

/* dmc_notnan() in memory: TRUE where x is present, NA where it is NA or
   NaN (as is.finite(x) | abs(x) == Inf gives). */
SEXP dm_notnan(SEXP x_)
{
    R_xlen_t n = XLENGTH(x_);
    SEXP out = PROTECT(allocVector(LGLSXP, n));
    int *o = LOGICAL(out);
    if (TYPEOF(x_) == REALSXP) {
	const double *x = REAL_RO(x_);
	for (R_xlen_t i = 0; i < n; i++) o[i] = ISNAN(x[i]) ? NA_LOGICAL : 1;
    } else if (TYPEOF(x_) == INTSXP || TYPEOF(x_) == LGLSXP) {
	const int *x = TYPEOF(x_) == INTSXP ? INTEGER_RO(x_) : LOGICAL_RO(x_);
	for (R_xlen_t i = 0; i < n; i++) o[i] = x[i] == NA_INTEGER ? NA_LOGICAL : 1;
    } else {
	error("internal error: dm_notnan() of a non-numeric vector");
    }
    UNPROTECT(1);
    return out;
}

/* f + fe_1[codes_1] + ... + fe_K[codes_K] (in that order), NA for a code
   outside 1..G_k: dm_fe_fitted() in memory, with one column. */
SEXP dm_fe_fitted_add(SEXP f_, SEXP fe_, SEXP codes_)
{
    R_xlen_t n = XLENGTH(f_);
    int K = LENGTH(fe_);
    if (LENGTH(codes_) != K) error("internal error: one code vector per FE");
    SEXP out = PROTECT(allocVector(REALSXP, n));
    double *o = REAL(out);
    memcpy(o, REAL_RO(f_), sizeof(double) * (size_t) n);
    for (int k = 0; k < K; k++) {
	SEXP v_ = VECTOR_ELT(fe_, k), c_ = VECTOR_ELT(codes_, k);
	if (XLENGTH(c_) != n) error("internal error: one code per row");
	const double *v = REAL_RO(v_);
	const int *c = INTEGER_RO(c_);
	R_xlen_t G = XLENGTH(v_);
	for (R_xlen_t i = 0; i < n; i++) {
	    int ci = c[i];
	    double a = (ci == NA_INTEGER || ci < 1 || ci > G) ? NA_REAL : v[ci - 1];
	    o[i] = o[i] + a;
	}
    }
    UNPROTECT(1);
    return out;
}

/* The residuals of dm_fe_core() in memory: f = Xb + fe_1[codes_1] + ...
   + fe_K[codes_K] (as dm_fe_fitted(), in that order), e = z - f, we = w e
   (e without weights), and in the same pass the sums of we within the
   levels of each FE (as dm_fe_gsum()) and sum(we e) (as sum()).
   list(we (an n x 1 matrix), Be (one G_k x 1 matrix per FE), rss). */
SEXP dm_fe_resid(SEXP z_, SEXP xb_, SEXP w_, SEXP fe_, SEXP codes_, SEXP G_)
{
    R_xlen_t n = XLENGTH(z_);
    int K = LENGTH(codes_);
    if (LENGTH(fe_) != K) error("internal error: one FE vector per grouping");
    const double *z = REAL_RO(z_), *xb = REAL_RO(xb_);
    const double *w = isNull(w_) ? NULL : REAL_RO(w_);
    if (XLENGTH(xb_) != n || (w && XLENGTH(w_) != n)) error("internal error: one value per row");
    const int *G = INTEGER_RO(G_);
    SEXP out = PROTECT(allocVector(VECSXP, 3));
    SEXP we_ = SET_VECTOR_ELT(out, 0, allocMatrix(REALSXP, n, 1));
    double *we = REAL(we_);
    SEXP Be = SET_VECTOR_ELT(out, 1, allocVector(VECSXP, K));
    const int **codes = (const int **) R_alloc(K > 0 ? K : 1, sizeof(int *));
    const double **fe = (const double **) R_alloc(K > 0 ? K : 1, sizeof(double *));
    R_xlen_t *nfe = (R_xlen_t *) R_alloc(K > 0 ? K : 1, sizeof(R_xlen_t));
    double **acc = (double **) R_alloc(K > 0 ? K : 1, sizeof(double *));
    for (int k = 0; k < K; k++) {
	SEXP ck = VECTOR_ELT(codes_, k);
	if (XLENGTH(ck) != n) error("internal error: one code per row");
	codes[k] = INTEGER_RO(ck);
	fe[k] = REAL_RO(VECTOR_ELT(fe_, k));
	nfe[k] = XLENGTH(VECTOR_ELT(fe_, k));
	if (nfe[k] < G[k]) error("internal error: one FE value per level");
	SEXP m = SET_VECTOR_ELT(Be, k, allocMatrix(REALSXP, G[k], 1));
	acc[k] = REAL(m);
	memset(acc[k], 0, sizeof(double) * (size_t) G[k]);
    }
    long double rss = 0.0;
    for (R_xlen_t i = 0; i < n; i++) {
	double f = xb[i];
	for (int k = 0; k < K; k++) {
	    int c = codes[k][i];
	    double a = (c == NA_INTEGER || c < 1 || c > nfe[k]) ? NA_REAL : fe[k][c - 1];
	    f = f + a;
	}
	double e = z[i] - f;
	double wei = w ? w[i] * e : e;
	double sq = wei * e;
	we[i] = wei;
	rss += sq;
	for (int k = 0; k < K; k++) acc[k][codes[k][i] - 1] += wei;
    }
    SET_VECTOR_ELT(out, 2, ScalarReal(dm_ld_value(rss)));
    UNPROTECT(1);
    return out;
}

/* The moments of dm_fe_core() in memory, as its R code sums them:
   c(sum(wz z), sum(w), sum(wz), sum(w (z - ybar)^2)) with wz = w z (z
   without weights; then sum(w) is n) and ybar = sum(wz) / sum(w). */
SEXP dm_fe_moments(SEXP z_, SEXP w_)
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
    double sw = w ? dm_ld_value(s2) : (double) n;
    double ybar = dm_ld_value(s3) / sw;
    for (R_xlen_t i = 0; i < n; i++) {
	double d = z[i] - ybar;
	double d2 = d * d;
	double t = w ? w[i] * d2 : d2;
	s4 += t;
    }
    SEXP out = PROTECT(allocVector(REALSXP, 4));
    REAL(out)[0] = dm_ld_value(s1);
    REAL(out)[1] = w ? dm_ld_value(s2) : 0.0;
    REAL(out)[2] = dm_ld_value(s3);
    REAL(out)[3] = dm_ld_value(s4);
    UNPROTECT(1);
    return out;
}

/* The working weights and responses of IRLS (dmglm_working() in memory),
   from the prior weights W, eta, mu, the response y, the offset off
   (length 1 or n) and the family's variance v and d mu / d eta me at mu:
   good = W > 0 & me != 0 (NA as R's & gives it), w2 = W me me / v and
   z = (eta - off) + (y - mu) / me where good (0 where not, NA where good
   is NA), and the counts of the checks of dmglm_working(), as its R code
   gives them in memory: 0, W > 0 and v == 0, 0, good, and rows where w2
   or z is not finite.  list(w2, z, counts). */
SEXP dm_glm_working(SEXP W_, SEXP eta_, SEXP mu_, SEXP y_, SEXP off_, SEXP v_, SEXP me_)
{
    R_xlen_t n = XLENGTH(eta_);
    SEXP args[6] = {W_, mu_, y_, off_, v_, me_};
    for (int a = 0; a < 6; a++) {
	R_xlen_t na = XLENGTH(args[a]);
	if (na != n && na != 1) error("internal error: lengths of the working columns");
    }
    const double *W = REAL_RO(W_), *eta = REAL_RO(eta_), *mu = REAL_RO(mu_), *y = REAL_RO(y_),
	*off = REAL_RO(off_), *v = REAL_RO(v_), *me = REAL_RO(me_);
    int sW = XLENGTH(W_) == n, smu = XLENGTH(mu_) == n, sy = XLENGTH(y_) == n,
	soff = XLENGTH(off_) == n, sv = XLENGTH(v_) == n, sme = XLENGTH(me_) == n;
    SEXP out = PROTECT(allocVector(VECSXP, 3));
    double *w2 = REAL(SET_VECTOR_ELT(out, 0, allocVector(REALSXP, n)));
    double *zz = REAL(SET_VECTOR_ELT(out, 1, allocVector(REALSXP, n)));
    double *cnt = REAL(SET_VECTOR_ELT(out, 2, allocVector(REALSXP, 5)));
    double c1 = 0, c3 = 0, c4 = 0;
    for (R_xlen_t i = 0; i < n; i++) {
	double Wi = W[sW ? i : 0], mei = me[sme ? i : 0], vi = v[sv ? i : 0];
	/* W > 0 and me != 0 as logicals: 1, 0 or NA */
	int wpos = ISNAN(Wi) ? NA_LOGICAL : (Wi > 0);
	int menz = ISNAN(mei) ? NA_LOGICAL : (mei != 0);
	int good = (wpos == 0 || menz == 0) ? 0 : ((wpos == NA_LOGICAL || menz == NA_LOGICAL) ? NA_LOGICAL : 1);
	/* (as in R: (W > 0) & !dmc_notnan(v) is never TRUE, dmc_notnan()
	   being TRUE or NA, so that the first and third counts stay 0) */
	if (wpos == 1 && !ISNAN(vi) && vi == 0) c1++;
	double a, b;
	if (good == NA_LOGICAL) {
	    a = NA_REAL;
	    b = NA_REAL;
	} else if (good) {
	    c3++;
	    double t1 = Wi * mei;
	    double t2 = t1 * mei;
	    a = t2 / vi;
	    double d1 = eta[i] - off[soff ? i : 0];
	    double d2 = y[sy ? i : 0] - mu[smu ? i : 0];
	    double d3 = d2 / mei;
	    b = d1 + d3;
	} else {
	    a = 0;
	    b = 0;
	}
	w2[i] = a;
	zz[i] = b;
	if (!R_FINITE(a) || !R_FINITE(b)) c4++;
    }
    cnt[0] = 0; cnt[1] = c1; cnt[2] = 0; cnt[3] = c3; cnt[4] = c4;
    UNPROTECT(1);
    return out;
}

/* The Poisson deviance residuals of dmglm_cf() in memory:
   y > 0 ? 2 w (y log(y / m) - (y - m)) : 2 m w, NA where y > 0 is NA
   (y, m, w of length 1 or n). */
SEXP dm_poisson_dev(SEXP y_, SEXP m_, SEXP w_)
{
    R_xlen_t ny = XLENGTH(y_), nm = XLENGTH(m_), nw = XLENGTH(w_);
    R_xlen_t n = ny > nm ? ny : nm;
    if (nw > n) n = nw;
    if ((ny != n && ny != 1) || (nm != n && nm != 1) || (nw != n && nw != 1))
	error("internal error: lengths of the deviance columns");
    const double *y = REAL_RO(y_), *m = REAL_RO(m_), *w = REAL_RO(w_);
    SEXP out = PROTECT(allocVector(REALSXP, n));
    double *o = REAL(out);
    for (R_xlen_t i = 0; i < n; i++) {
	double yi = y[ny == 1 ? 0 : i], mi = m[nm == 1 ? 0 : i], wi = w[nw == 1 ? 0 : i];
	if (ISNAN(yi)) {
	    o[i] = NA_REAL;
	} else if (yi > 0) {
	    double r = yi / mi;
	    /* R's log() */
	    double l = r > 0 ? log(r) : (r == 0 ? R_NegInf : R_NaN);
	    double a = yi * l;
	    double b = yi - mi;
	    double c = a - b;
	    double tw = 2 * wi;
	    o[i] = tw * c;
	} else {
	    double tm = 2 * mi;
	    o[i] = tm * wi;
	}
    }
    UNPROTECT(1);
    return out;
}

/* The variables of a design with ALTREP numeric vectors replaced by plain
   copies (same values and attributes).  stats sets the contrasts of a
   factor on a column shared with the model frame, which makes it a
   wrapper; the kernels take writable pointers to the variables, and R
   then copies a wrapper's data at every call.  One copy here instead. */
SEXP dm_unwrap_vars(SEXP vars)
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

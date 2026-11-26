#ifndef DESIGN_MATRIX_H
#define DESIGN_MATRIX_H
#include <Rinternals.h>

/* Ceiling division for nonnegative row counts and positive block sizes. */
static inline R_xlen_t dm_block_count(R_xlen_t rows, R_xlen_t size)
{
    return rows == 0 ? 0 : 1 + (rows - 1) / size;
}

static inline int dm_worker_count(int requested, R_xlen_t units)
{
#ifdef _OPENMP
    if (requested > 1 && units > 1)
        return requested < units ? requested : (int) units;
#endif
    return 1;
}

/* Reductions (sums over rows) on several threads: the rows are split into
   T chunks whose boundaries depend only on the number of rows and T, each
   chunk is summed into its own accumulator, and the accumulators are added
   in chunk order afterwards.  Results are identical for a given T (which
   depends only on the requested threads and the problem's dimensions);
   with T = 1 the code path is the single-threaded one.

   dm_reduce_threads(): T for a reduction over `rows` rows costing `work`
   multiply-adds in all, with accumulators of `acc` doubles.  One thread for
   small problems (fewer than DM_RED_MIN_ROWS rows or DM_RED_MIN_WORK
   multiply-adds per thread), and fewer threads when the extra accumulators
   would take more than DM_RED_MAX_BYTES, or adding them up would cost more
   than a quarter of the work. */
#define DM_RED_MIN_ROWS 8192
#define DM_RED_MIN_WORK 2.5e5
#define DM_RED_MAX_BYTES 268435456.0

static inline int dm_reduce_threads(int requested, R_xlen_t rows, double work, double acc)
{
    int T = dm_worker_count(requested, rows / DM_RED_MIN_ROWS);
    while (T > 1) {
	double extra = (double) (T - 1) * acc;
	if (work >= DM_RED_MIN_WORK * T && extra * sizeof(double) <= DM_RED_MAX_BYTES &&
	    4.0 * extra <= work)
	    break;
	T--;
    }
    return T;
}

/* First row of chunk t of T (t = T gives `rows`): chunks are whole blocks
   of `unit` rows, as even as possible. */
static inline R_xlen_t dm_chunk_start(R_xlen_t rows, int T, int t, R_xlen_t unit)
{
    R_xlen_t nb = dm_block_count(rows, unit);
    R_xlen_t i = (nb * t / T) * unit;
    return i < rows ? i : rows;
}

/* acc[0][j] += acc[1][j] + ... + acc[T - 1][j], in that order, j < len. */
static inline void dm_combine(double **acc, int T, R_xlen_t len)
{
    if (T <= 1) return;
#ifdef _OPENMP
#pragma omp parallel for num_threads(T) schedule(static) if(len >= 65536)
#endif
    for (R_xlen_t j = 0; j < len; j++) {
	double s = acc[0][j];
	for (int t = 1; t < T; t++) s += acc[t][j];
	acc[0][j] = s;
    }
}

/* Accumulators for chunks 1 .. T - 1 (chunk 0 uses `first`): acc[t]
   points to `len` doubles, which the chunk zeroes (dm_zero()) on its
   thread. */
static inline double **dm_chunk_acc(double *first, int T, R_xlen_t len)
{
    double **acc = (double **) R_alloc(T, sizeof(double *));
    acc[0] = first;
    for (int t = 1; t < T; t++)
	acc[t] = (double *) R_alloc(len > 0 ? len : 1, sizeof(double));
    return acc;
}

static inline void dm_zero(double *x, R_xlen_t len)
{
    for (R_xlen_t j = 0; j < len; j++) x[j] = 0.0;
}

SEXP dm_setup(SEXP terms, SEXP vars, SEXP rho);
SEXP dm_vars(SEXP setup, SEXP vars);
SEXP dm_fill(SEXP setup, SEXP vars, SEXP nthreads);
SEXP dm_jac_mult(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_jac_crossprod(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_jac_quad(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_jac_sparse(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_jac_dense(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_row_groups(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_jac_entries_mean(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_jac_finite(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_jac_gram(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_chol_pivot(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_struct_ok(SEXP, SEXP, SEXP);
SEXP dm_struct_info(SEXP, SEXP);
SEXP dm_struct_gram(SEXP, SEXP, SEXP, SEXP);
SEXP dm_struct_xty(SEXP, SEXP, SEXP, SEXP);
SEXP dm_struct_mult(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_kron_mode(SEXP, SEXP, SEXP);
SEXP dm_min_nonzero(SEXP, SEXP);
SEXP dm_all_finite(SEXP);
SEXP dm_col_norms(SEXP);
SEXP dm_in_range(SEXP, SEXP, SEXP);
SEXP dm_codes_ok(SEXP, SEXP);
SEXP dm_int_unique(SEXP);
SEXP dm_str_unique(SEXP);
SEXP dm_fe_gsum(SEXP, SEXP, SEXP);
SEXP dm_fe_crosstab(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_fe_components(SEXP, SEXP, SEXP, SEXP);
SEXP dm_fe_struct_gsum(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_fe_code(SEXP, SEXP);
SEXP dm_fe_recode(SEXP, SEXP);
SEXP dm_fe_nested(SEXP, SEXP, SEXP);
SEXP dm_fe_lookup(SEXP, SEXP);
/* fused.c */
SEXP dm_ifelse(SEXP, SEXP, SEXP);
SEXP dm_notnan(SEXP);
SEXP dm_fe_fitted_add(SEXP, SEXP, SEXP);
SEXP dm_fe_resid(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_fe_moments(SEXP, SEXP);
SEXP dm_glm_working(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_poisson_dev(SEXP, SEXP, SEXP);
SEXP dm_unwrap_vars(SEXP);
SEXP dm_levels_unused(SEXP, SEXP);
SEXP dm_ls_sums(SEXP, SEXP, SEXP, SEXP);
SEXP dm_binomial_ll(SEXP, SEXP, SEXP, SEXP);
SEXP dm_fe_ichol(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_fe_pcg(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_fe_schur(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_omp_max_threads(void);
#endif

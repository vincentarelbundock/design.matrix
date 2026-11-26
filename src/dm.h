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

SEXP dm_setup(SEXP terms, SEXP vars, SEXP rho);
SEXP dm_vars(SEXP setup, SEXP vars);
SEXP dm_fill(SEXP setup, SEXP vars, SEXP nthreads);
SEXP dm_jac_mult(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_jac_crossprod(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_jac_quad(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_jac_sparse(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_jac_dense(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_row_groups(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_jac_entries_mean(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_jac_finite(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_jac_gram(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_chol_pivot(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_struct_ok(SEXP, SEXP, SEXP);
SEXP dm_struct_info(SEXP, SEXP);
SEXP dm_struct_gram(SEXP, SEXP, SEXP);
SEXP dm_struct_xty(SEXP, SEXP, SEXP);
SEXP dm_struct_mult(SEXP, SEXP, SEXP, SEXP);
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
SEXP dm_fe_code_int(SEXP, SEXP);
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
SEXP dm_fe_ichol(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_fe_pcg(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_fe_schur(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dm_omp_max_threads(void);
#endif

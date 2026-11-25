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
SEXP dm_omp_max_threads(void);
#endif

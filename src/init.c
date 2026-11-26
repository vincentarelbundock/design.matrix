#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "dm.h"
#ifdef _OPENMP
#include <omp.h>
#endif

SEXP dm_omp_max_threads(void)
{
#ifdef _OPENMP
    return ScalarInteger(omp_get_max_threads());
#else
    return ScalarInteger(1);
#endif
}

static const R_CallMethodDef CallEntries[] = {
    {"dm_setup", (DL_FUNC) &dm_setup, 3},
    {"dm_fill", (DL_FUNC) &dm_fill, 3},
    {"dm_jac_mult", (DL_FUNC) &dm_jac_mult, 9},
    {"dm_jac_crossprod", (DL_FUNC) &dm_jac_crossprod, 10},
    {"dm_jac_quad", (DL_FUNC) &dm_jac_quad, 9},
    {"dm_jac_sparse", (DL_FUNC) &dm_jac_sparse, 6},
    {"dm_jac_dense", (DL_FUNC) &dm_jac_dense, 7},
    {"dm_row_groups", (DL_FUNC) &dm_row_groups, 8},
    {"dm_jac_entries_mean", (DL_FUNC) &dm_jac_entries_mean, 7},
    {"dm_jac_finite", (DL_FUNC) &dm_jac_finite, 7},
    {"dm_jac_gram", (DL_FUNC) &dm_jac_gram, 14},
    {"dm_chol_pivot", (DL_FUNC) &dm_chol_pivot, 5},
    {"dm_struct_ok", (DL_FUNC) &dm_struct_ok, 3},
    {"dm_struct_info", (DL_FUNC) &dm_struct_info, 2},
    {"dm_struct_gram", (DL_FUNC) &dm_struct_gram, 3},
    {"dm_struct_xty", (DL_FUNC) &dm_struct_xty, 3},
    {"dm_struct_mult", (DL_FUNC) &dm_struct_mult, 4},
    {"dm_omp_max_threads", (DL_FUNC) &dm_omp_max_threads, 0},
    {NULL, NULL, 0}
};

void R_init_design_matrix(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

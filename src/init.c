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
    {"dm_jac_crossprod", (DL_FUNC) &dm_jac_crossprod, 11},
    {"dm_jac_quad", (DL_FUNC) &dm_jac_quad, 9},
    {"dm_jac_sparse", (DL_FUNC) &dm_jac_sparse, 6},
    {"dm_jac_dense", (DL_FUNC) &dm_jac_dense, 7},
    {"dm_row_groups", (DL_FUNC) &dm_row_groups, 8},
    {"dm_jac_entries_mean", (DL_FUNC) &dm_jac_entries_mean, 7},
    {"dm_jac_finite", (DL_FUNC) &dm_jac_finite, 7},
    {"dm_jac_gram", (DL_FUNC) &dm_jac_gram, 15},
    {"dm_chol_pivot", (DL_FUNC) &dm_chol_pivot, 5},
    {"dm_struct_ok", (DL_FUNC) &dm_struct_ok, 3},
    {"dm_struct_info", (DL_FUNC) &dm_struct_info, 2},
    {"dm_struct_gram", (DL_FUNC) &dm_struct_gram, 4},
    {"dm_struct_xty", (DL_FUNC) &dm_struct_xty, 4},
    {"dm_struct_mult", (DL_FUNC) &dm_struct_mult, 5},
    {"dm_kron_mode", (DL_FUNC) &dm_kron_mode, 3},
    {"dm_min_nonzero", (DL_FUNC) &dm_min_nonzero, 2},
    {"dm_all_finite", (DL_FUNC) &dm_all_finite, 1},
    {"dm_col_norms", (DL_FUNC) &dm_col_norms, 1},
    {"dm_in_range", (DL_FUNC) &dm_in_range, 3},
    {"dm_codes_ok", (DL_FUNC) &dm_codes_ok, 2},
    {"dm_int_unique", (DL_FUNC) &dm_int_unique, 1},
    {"dm_str_unique", (DL_FUNC) &dm_str_unique, 1},
    {"dm_fe_gsum", (DL_FUNC) &dm_fe_gsum, 3},
    {"dm_fe_crosstab", (DL_FUNC) &dm_fe_crosstab, 5},
    {"dm_fe_components", (DL_FUNC) &dm_fe_components, 4},
    {"dm_fe_struct_gsum", (DL_FUNC) &dm_fe_struct_gsum, 7},
    {"dm_fe_code_int", (DL_FUNC) &dm_fe_code_int, 2},
    {"dm_fe_nested", (DL_FUNC) &dm_fe_nested, 3},
    {"dm_fe_lookup", (DL_FUNC) &dm_fe_lookup, 2},
    {"dm_fe_ichol", (DL_FUNC) &dm_fe_ichol, 5},
    {"dm_fe_pcg", (DL_FUNC) &dm_fe_pcg, 16},
    {"dm_fe_schur", (DL_FUNC) &dm_fe_schur, 10},
    {"dm_omp_max_threads", (DL_FUNC) &dm_omp_max_threads, 0},
    {NULL, NULL, 0}
};

void R_init_design_matrix(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

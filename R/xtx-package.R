## Names after the rename design.matrix -> xtx (Kata design.matrix#c8wd).
## Public functions kept their names (design_matrix(), design_column(),
## col_*(), group_*(), ...). Everything else moved to the xtx prefix:
##
##   classes   design_matrix -> xtx_matrix, dm_remote -> xtx_remote,
##             dm_remote_table -> xtx_remote_table, design_column ->
##             xtx_column, design_groups -> xtx_groups, design_plan ->
##             xtx_plan, design_qr -> xtx_qr, design_formula -> xtx_formula,
##             design_lm, design_glm, design_lm_fe, design_glm_fe -> xtx_lm,
##             xtx_glm, xtx_lm_fe, xtx_glm_fe (and summary.*); every other
##             dm_* class (dm_remote_exclude, dm_fe_op, ...) -> xtx_*
##   internal  dm_ -> xtx_, dmp_ -> xtx_p_ (primitives), dmc_ -> xtx_c_
##             (remote columns), dmfit_ -> xtx_fit_, dmglm_ -> xtx_glm_,
##             dmqr_ -> xtx_qr_; two names that would collide keep an infix:
##             dmglm_working() -> xtx_glm_working_memory() (the C entry
##             dm_glm_working is xtx_glm_working), and the row-dedup
##             dm_groups() -> xtx_dedup_groups() (xtx_groups is the class)
##   C         dm_* entry points -> xtx_*, DM_* macros -> XTX_*, dm.h -> xtx.h
##   options   design.matrix.* -> xtx.*; environment DESIGN_MATRIX_* -> XTX_*
##   DuckDB    temporary names __dm_* -> __xtx_*, dm_* macros/tables -> xtx_*

#' xtx: Exact and Fast Model Matrices
#'
#' Builds model matrices that are identical, bit for bit, to those of
#' [stats::model.matrix()], faster, and offers algebra on the design
#' without forming the dense matrix.
#'
#' @useDynLib xtx, .registration = TRUE, .fixes = "C_"
#' @import stats
#' @keywords internal
"_PACKAGE"

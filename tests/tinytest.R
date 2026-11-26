# Default (fast) tier: every feature in memory and on a DuckDB table. Set
# DESIGN_MATRIX_TESTS=full for the full tier (Parquet, every parity loop,
# the slow regression tests); see inst/tinytest/helpers.R.
if (requireNamespace("tinytest", quietly = TRUE)) {
  tinytest::test_package("design.matrix")
}

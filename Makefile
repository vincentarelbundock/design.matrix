.PHONY: install reference website serve bench bench-ols bench-quick

install: ## Install the package
	R CMD INSTALL .

reference: ## Convert man/*.Rd to website/reference/*.typ
	./website/reference.sh

# Benchmarks never run when the site builds: the pages read saved results.
# BENCH_SIZE=full (default) saves the published results; BENCH_SIZE=quick
# runs small grids and saves them in $BENCH_OUT/website and $BENCH_OUT/ols
# (default BENCH_OUT: $TMPDIR/design.matrix-bench-quick), never over the
# published ones. The fixed-effects data cache is outside the repository too
# (FE_BENCH_CACHE, see website/benchmarks/fixed_effects.R).
BENCH_SIZE ?= full

bench: install ## Run website benchmarks (all, or B="crossprod qr"; BENCH_SIZE=quick|full)
	BENCH_SIZE=$(BENCH_SIZE) Rscript website/benchmarks/run.R $(B)

bench-ols: install ## Run the least squares comparison (fresh process per method; BENCH_SIZE=quick|full)
	BENCH_SIZE=$(BENCH_SIZE) ./bench/ols/run.sh

bench-quick: install ## Quick check of every published benchmark; results outside the repository
	BENCH_SIZE=quick Rscript website/benchmarks/run.R $(B)
	BENCH_SIZE=quick ./bench/ols/run.sh

website: install reference ## Build the site into _site/
	calepin compile website _site

serve: ## Serve the built site
	calepin serve _site --open

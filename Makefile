.PHONY: install reference website serve bench bench-ols

install: ## Install the package
	R CMD INSTALL .

reference: ## Convert man/*.Rd to website/reference/*.typ
	./website/reference.sh

# Benchmarks never run when the site builds: the pages read saved results.
bench: install ## Run website benchmarks (all, or B="crossprod qr")
	Rscript website/benchmarks/run.R $(B)

bench-ols: install ## Run the least squares comparison (fresh process per method)
	./bench/ols/run.sh

website: install reference ## Build the site into _site/
	calepin compile website _site

serve: ## Serve the built site
	calepin serve _site --open

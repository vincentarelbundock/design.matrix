#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: [Architecture])
#metadata((title: "Architecture")) <website-metadata>

#calepin.setup(echo: true, eval: false)

#title()

This page describes how the package computes, for readers who want to know why the results are what they are, where the time goes, or where to look in the source. File names refer to the #link("https://github.com/vincentarelbundock/xtx")[repository].

= The idea

A model matrix is highly redundant. A factor with 1,000 levels fills 999 columns, but each row holds one integer code. An interaction of a factor and a numeric variable fills as many columns again, but each row holds one code and one number. The package keeps that compact description (factor codes, numeric variables, contrast matrices and the layout of the columns) and computes with it.

Three choices follow from this, and the rest of the page describes them:

+ *A design is a description of the matrix.* The dense matrix is built only by `as.matrix()`, and then bit for bit as `model.matrix()` builds it.
+ *Estimation is a sequence of sums over rows.* $X'X$, $X'y$, grouped column sums and crosstabs are computed from codes and numbers. Linear models, GLMs and fixed effects are fitted from these sums: the rows are never transformed (no demeaning, no dense $Q$).
+ *The fitting code is written once.* It calls a small set of primitives, each with a method for data in memory (C) and for remote data (SQL), through exported functions only. There is no second implementation for out-of-core data.

= Layers

#table(
  columns: 3,
  table.header([Layer], [Role], [Source]),
  [Frame and plan], [The model frame of R, and the fixed transformations, levels and contrasts that new data must reuse], [`R/frame.R`, `R/prep.R`, `R/plan.R`],
  [Design], [Checks, contrasts and column names as in R; the dense fill], [`src/modelmatrix.c`, `R/design_matrix.R`],
  [Kernels], [Products with the design from codes and numbers], [`src/kernels.c`, `src/structured.c`, `R/algebra.R`],
  [Primitives], [Operations on designs and columns, in memory and remote; the public column type], [`R/primitives.R`, `R/primitives_lookup.R`, `R/column.R`, `R/columns_public.R`],
  [Least squares], [`qr()` of a design: pivoting, the solve and its refinement, exact centering], [`R/qr.R`, `R/center.R`],
  [Fitting layers (package `xtx.fit`)], [Weighted least squares, the input stage, GLM families on columns, the fixed-effects estimator and its solver], [`xtx.fit`: `R/wls.R`, `R/problem.R`, `R/family.R`, `R/kernels.R`, `R/absorb.R`, `src/fe_solver.c`, `src/kernels.c`],
  [Fits (package `xtx.fit`)], [`lm_xtx()`, `glm_xtx()`: IRLS, checks, covariances, methods], [`xtx.fit`: `R/lm.R`, `R/glm.R`, `R/fe.R`, `R/fe_glm.R`, `R/remote_fit.R`, `R/checks.R`, `R/na_action.R`],
  [Remote], [Designs over DuckDB tables, queries and Parquet files], [`R/remote.R`, `R/remote_sql.R`, `R/remote_compress.R`],
)

Each layer uses only the ones above it in the table. Everything that is about fitting a model lives in a package of its own, `xtx.fit`, which reaches `xtx` only through its exported functions; `xtx` itself stops at the design, its products and its decompositions. The fitting code does not know which backend it runs on.

= A design is a description

== R's own preparation

`model.matrix()` and `model.frame()` do a good deal of work in R before they reach C: non-standard evaluation, character columns turned into factors, contrasts, `xlev`, `na.action`. This code changes between R versions. The package does not copy it. At load time it takes the installed `stats:::model.matrix.default()`, cuts it just before its call to C, and makes it return what that call would have received (`R/prep.R`). The same is done for `model.frame.default()`, with two patches that skip copies which change nothing (`R/frame.R`). Each patch is applied only when the installed code is the code it was written for, and a self-test at load drops a patch that fails; the fallback is the function of stats.

== Setup and fill

R's C routine `modelmatrix()` is ported and split in two (`src/modelmatrix.c`):

- `xtx_setup()` performs every check, evaluates the contrasts, emits the warnings and computes the column names, in the order of R;
- `xtx_fill()` allocates the matrix and fills it.

A `design_matrix` is the result of the setup plus the variables. The fill performs the same floating point operations as R, cell by cell. Only the order in which rows are visited changes (chunks of rows, on several threads when asked), so `as.matrix(D)` is identical to `model.matrix()`, including `NA`, `NaN`, `Inf` and the sign of zero.

== Plans

A plan keeps what a fitted design fixes for later data: the terms, the data-dependent transformations (`poly()`, `ns()`, `scale()`) with their training parameters, the factor levels and the contrasts (`R/plan.R`). `design_matrix(plan, batch)` applies it to new rows without repeating the checks. Remote designs carry a plan too: it fixes their terms, levels and contrasts, so that their columns are those of `model.matrix()`. `design_matrix(fit, newdata = NULL)` returns the plan of a fitted model, and `design_parts()` is an exported generic: a package supports its own model class with a `design_parts()` method that returns the terms, levels and contrasts its `predict()` method uses (`R/parts.R`).

= Products from the structure

== Terms

Take term $t$ of the design, with factors and numeric variables. On row $i$, let $F_t (i)$ be the combined level of its factors, $N_t (i)$ the product of its numeric variables, and $C_t$ the Kronecker product of the contrast matrices of its factors, with one row per combined level. The block of columns of the term is

$ X_t [i, dot] = N_t (i) thin C_t [F_t (i), dot]. $

A product with the design then separates into a sum over rows, which involves only codes and numbers, and a small dense product with the $C_t$. For two terms $s$ and $t$,

$ (X' W X)_(s t) = C_s ' M_(s t) C_t, quad M_(s t) [a, b] = sum_(i : F_s (i) = a, thin F_t (i) = b) w_i N_s (i) N_t (i). $

$M_(s t)$ is a weighted crosstab. In memory it is one pass over the rows in C (`src/structured.c`); in a database it is a `GROUP BY` over the factors of the two terms (`R/remote_sql.R`). $X' W y$, column sums and grouped column sums have the same form with one term.

== Two sets of kernels

- The *structured kernels* (`src/structured.c`) use the factorization above. They are used when every input is finite and no product can overflow; the results then agree with the dense computation up to rounding.
- The *row kernels* (`src/kernels.c`) list, for each row, the columns that can be non-zero, with their values computed in the order of the dense fill. They handle every other case (so that `0 * Inf` is `NaN` where the dense matrix has it), Jacobians of two designs, and results with one value per row ($X beta$, $op("diag")(X V X')$).

Both can run on several threads (`design_threads()`; as in fixest, the default is half of the available threads, and one during `R CMD check`). Results with one value per row do not depend on the number of threads. Sums over rows split the rows into fixed blocks, one per thread, added in a fixed order: they are identical for a given number of threads, and can differ in the last digits between numbers of threads.

== Patterns

Rows with the same factor codes and the same numeric values give the same design row. When a design has far fewer such patterns than rows, the kernels compute once per pattern and expand, or sum the weights within patterns first (`R/algebra.R`, option `xtx.dedup`). The fits use the same idea as compression: a weighted least squares fit on the patterns, with the sum of the weights and the weighted mean response of each, is the fit on the rows, and the within-pattern sum of squares completes the residual sum of squares (`R/lm.R` of `xtx.fit`). For a GLM, the iterations on the patterns are those on the rows for every supported family (`R/glm.R` of `xtx.fit`). `design_patterns()` exports this compression, in memory and in a database, as one list: the pattern design, the mean response and summed weights of each pattern, and the within-pattern sum of squares (`R/exports.R`).

= Least squares

`qr()` of a design (`R/qr.R`) decomposes $X P = Q R$ without the dense $X$ or the $n times p$ matrix $Q$:

- $R$ is the Cholesky factor of $X'X$, from `crossprod()`;
- columns are pivoted with the rule of `lm()` (LINPACK's limited pivoting). From $X'X$ a residual norm is accurate to about the square root of machine precision, so the columns near the tolerance are judged from the data, with one product per block of columns;
- coefficients solve the normal equations and are refined from the residuals $y - X b$, each step costing one $X b$ and one $X' r$;
- numeric variables far from zero are centered when that is an exact linear transform $X T$ of the design, which keeps $X'X$ accurate: the columns kept come from $X'X$, the values of $R = tilde(R) T^(-1)$ from the centered design. Remote designs are always centered; in memory, only when $X'X$ is badly conditioned ($kappa^2 > 10^6$, columns scaled to unit norm), so that other fits take no extra pass.

The cost is a few passes over the rows and memory of order $p^2$. The price is that of the normal equations: on badly conditioned designs the coefficients agree with a Householder QR to fewer digits, and when $X'X$ cannot represent the scale of the columns `qr()` falls back to the dense decomposition.

= One fitting code, two backends

== The primitive contract

`R/primitives.R` defines the operations that fitting code may use. They take a *design* and *columns* (one value per row):

#table(
  columns: 3,
  table.header([], [In memory], [Remote]),
  [Design], [`design_matrix`], [`xtx_remote`: a source, a plan, and how to compute each variable in SQL],
  [Column], [numeric vector or matrix], [`design_column`: SQL expressions over the rows of a frame],
)

The operations are products with the design ($X'X$, $X'c$, $X B$, $op("rowSums")((X V) circle.small X)$), sums and grouped sums of columns, crosstabs of two groupings, lookups of a value per level, subsets, and arithmetic on columns (`+`, `*`, `exp()`, `log()`, comparisons, the functions the GLM families need). Arithmetic on a remote column builds SQL; nothing runs until a sum is asked for. Products and lookups are lazy too: a product (or a lookup of values per level) is a recipe, which becomes a view (or a join) over the rows only at its first use, on the latest frame among the columns of that use, the products of one use stacked one on the other. Each view keeps the columns of the view below it, so `D %*% b1 + D %*% b2` is one chain of two views, and the fitting code never says on which frame a product goes. A product already bound on a sibling frame (read on its own before the combination) is replayed on the frame it joins (see below). Missing values follow R: SQL `NULL` is `NA`, and sums return `NA` where `sum()` does.

`lm.R`, `glm.R`, `fe.R` and `fe_glm.R` contain no branch on the backend. A fit on a DuckDB table is the same sequence of operations as in memory, which is why its coefficients, aliased columns, covariance and degrees of freedom are those of `lm()` or `glm()` on the same rows.

== One public column type

The remote column of the fitting core is also the column users get (`R/column.R`). `D %*% b`, `rowSums(D)`, `rowMeans(D)` and `quadratic_diagonal(D, V)` of a remote design return a `design_column` where a design in memory returns a vector, and `design_column(D, "y")` or `design_column(D, DBI::SQL(...))` makes one from the data, so that public calls chain: `exp(D %*% b) - design_column(D, "y")` is a column like the ones the fits build. The base generics `as.vector()` and `as.matrix()` read the values into R (a `SELECT` of every row), `as.data.frame()` reads them with columns of the data (`keep`, for instance a key to match rows), and `col_compute()` stores the column in a temporary table (the frame of `xtx_c_freeze(table = TRUE)`), or with `name =` in a permanent table of the connection that holds the same rows and columns (the columns of the data, internal columns and the values), so that the column returned stays aligned with the rows of its design and the table can be used outside the package. In memory a column is a vector, and base R applies. As for a matrix in memory, `dim()` of a remote column is $(n, k)$, with $n$ the row count kept with the design (no query), and `length()` is $n k$; the fitting core counts columns with `xtx_ncols()`, never with `length()`. A column is valid while its connection is open; the views and temporary tables under it are temporaries of that session. Frames, their depth and lineage stay internal: users see columns, which combine when they come from the same design, including separate products such as `(D %*% b) / sqrt(quadratic_diagonal(D, V))`. Two products each read before they are combined sit on two branches of views, whose rows the database cannot align: each frame records the recipes that made it (a product's design and coefficients, a lookup's table and join), so the product (or lookup) of one branch is *replayed* on top of the other, one more view for a product and no statement for a lookup, under the same column names, and the columns of both branches then read one frame. What cannot be replayed is a stored or frozen column (`col_compute()`, `xtx_c_freeze()`): two such columns stored separately on two branches remain an error, as do columns of two designs (a subset or a materialized design is another design). Results of fits with one value per row (`fitted()`, `predict()`, `qr.fitted()`) remain tables that keep the columns of the data.

Base generics complete the type (`R/columns_public.R`): `sum()`, `prod()`, `min()`, `max()`, `range()`, `any()`, `all()` and `mean()` of remote columns are one query each, whatever the number of arguments, with R's rule for missing values (`NA` if a value is `NA`, else `NaN` if one is `NaN`; `na.rm` drops both); `round(x, digits)` and `signif()` follow R's algorithms step by step in SQL. What base R cannot dispatch on is exported as `col_ifelse()`, `col_pmin()`, `col_pmax()` and `col_sums()` (several sums, one query), thin wrappers of the primitives. `e * D` scales the rows of a remote design by a column `e` (R's `chooseOpsMethod()` picks the design's method), so that cluster scores are `colSums(e * D, group = cl)` on both backends.

== Groupings

In memory, fixed effects use a dense integer code per observation. DuckDB fits aggregate on the source keys, return group totals and pair counts, and match those summary rows to solver indices in R. Interacted fixed effects use a single sort key built from the tuple. Refinement joins the effects directly on these keys. The exported grouping functions (`R/groups_public.R`) expose integer codes on both backends. `design_groups(D, key)` returns the codes (an integer vector in memory, a `design_column` remotely), $G$ and the levels (a data frame, one row per code, sorted as the levels of fixed effects: radix order, strings in the C locale); `group_sums(col, g)` is a $G times k$ matrix of sums (one `GROUP BY` remotely, one pass of a C kernel in memory), `group_crosstab(g1, g2, w)` a sparse $G_1 times G_2$ matrix of weighted counts, `group_lookup(values, g)` the column `values[codes]`, and `group_nested(g1, g2)` whether each level of `g1` lies within one level of `g2`. In memory the codes come from one pass in C (`xtx_fe_code()`); remotely from a `SELECT DISTINCT` of the key, sorted in R and written to a temporary table, which a join turns into codes. Each use of a remote grouping makes a fresh lookup of that table (a new recipe, no statement), so the join binds on top of the other columns of the use, whatever the grouping was used with before: a grouping combines with any column and any other grouping of its design. A missing key is an error on both backends (the fits drop those rows first; a grouping has no rows of its own), and codes are checked to lie in $1..G$ before the kernels index by them. The within transformation of a one-way fixed effect, `y - group_lookup(group_sums(y, g) / group_sums(1, g), g)`, is then the same code on a data frame and on a DuckDB table.

== The SQL engine

For a remote design, the products of the previous section are `GROUP BY` aggregates over factor columns as stored. Levels become codes in R, on the cells that come back, and the $C_t$ come from the set-up of the in-memory design of a few prototype rows (the contrasts it chose for each factor of each term, and the column of each combination), so that columns, contrasts and their order are those of `model.matrix()`. The sums of a term, or of a pair of terms for $X' W X$, are grouped by the factors of that term or pair, in one statement (`GROUP BY GROUPING SETS`): what comes back is bounded by the cells of the terms involved, not by the cells of all the factors of the design together, and never grows with the rows. R applies the contrasts to those sums.

A variable that SQL cannot compute (an R function with no translation) is evaluated in R once, batch by batch, into a column of a temporary table, and is then a column like any other. Results of fits with one value per row (fitted values, residuals, predictions) are views in the database, which keep the columns of the data; products with the design are columns.

== Compression in the database

When the data have few design patterns, the scan of the SQL engine also groups by the numeric variables of the design, and returns each cell with its number of rows and the sums of the weights and responses; the design of the cells is built in memory from the engine's own $C_t$ and the fit runs there (`R/remote_compress.R`). An estimate of the number of patterns, from a sample and from `approx_count_distinct()`, decides whether to try. The cells of a trial are counted in the database: when there are more than the cap, or as many as rows, the statement returns none of them and the ordinary scan runs.

== The API of the fitting code

The fitting code is the package `xtx.fit`, which exports `lm_xtx()` and `glm_xtx()` only: the fits (`lm.R`, `glm.R`, `fe.R`, `fe_glm.R`, `remote_fit.R`, `checks.R`, `na_action.R`) and the layers they are built on, which are internal to `xtx.fit` (they are documented in its sources, and its tests reach them through the namespace):

#table(
  columns: 2,
  table.header([Internal function of `xtx.fit`], [Role]),
  [`design_problem()`], [The input stage: a formula and a data frame or remote data, a design matrix or a remote design, as the design of the rows used and the response, weights and offset as columns (class `xtx_problem`)],
  [`design_exclude()`], [`na.exclude` for remote fits (`napredict()` and `naresid()` methods)],
  [`design_wls()`], [Weighted least squares with `lm()`'s pivoting, in memory or remote],
  [`column_family()`, `col_irls_working()`], [GLM families as functions of columns (links, variances, deviances, starting values, log densities); the working weights and responses of IRLS],
  [`design_absorb_rows()`], [The rows of a fit with fixed effects and the codes of the FE (missing values, singletons, levels perfectly fit)],
  [`design_absorb_wls()`, `design_absorb()`, `absorb_control()`], [The estimator with absorbed fixed effects: coefficients, bread, FE values, degrees of freedom; the absorption operator, whose `scores()` are the raw cluster scores (the fit applies fixest's small sample corrections); solver options],
  [`design_absorb_fitted()`, `design_uncenter()`], [$X b$ plus the FE values of the rows' levels; the map of a fit on a centered design back to the design],
)

`xtx.fit` reaches `xtx` only through its documented, exported functions (no internal function is exported for it). The functions that exist only for fitting code are exported and documented but kept out of the help index (`@keywords internal`): `design_resolve()`, `design_output()`, `design_center()`, `design_moments()`, `pivoted_chol()`, `design_subset()`, `design_key()`, `key_codes()`, `design_group_sums()` and `col_math()`. Besides base generics on designs and columns (`%*%`, arithmetic, `sum()`, `is.na()`, `is.finite()`, `round()`, `log()`, `cbind()`, `x[, j]`, `rowSums()`, `==` of a text column and a string), it uses:

#table(
  columns: 2,
  table.header([Function of `xtx`], [Role]),
  [`design_matrix()`, `design_frame()`], [The design and the model frame. With remote data (a DuckDB connection or Parquet paths as `data`), `design_frame()` is the single entry point of a model: it opens the data, takes `weights`, `offset` and `subset` unevaluated (SQL, or evaluated in R into a table), the variables of fixed effects and clusters as extra columns that must be present, and returns a remote frame (class `xtx_remote_frame`): the remote design of the rows used, the response with its kind (numeric, logical, a factor with its levels, the columns of `cbind()`, or a type no fit takes), the weights and offset as columns, and the column types. What a fit makes of the response (a factor needs a binomial family, its first level is the reference, the trials of `cbind()`) is decided in `xtx.fit`, which branches on the class of the result],
  [`design_resolve(response =)`], [Arguments of a fit as columns of an existing design (SQL, or evaluated in R into a table), and a response resolved as in a remote frame],
  [`design_patterns(check = FALSE, fit =)`], [Compression: rows that share a design pattern and their sufficient statistics],
  [`design_output()`], [Results with one value per row: vectors in memory, tables remotely],
  [`col_sums()`, `col_ifelse()`, `col_pmax()`, `col_const()`, `col_compute(table = FALSE)`], [Sums of several columns in one query (one `na.rm` per argument), conditionals, constants, columns computed once in a projection],
  [`design_center()`, `design_groups(levels =)`], [Exact centering of numeric variables; groupings of new rows coded against the levels of a fit],
  [`design_moments()`, `qr(D, y =)`, `pivoted_chol()`], [The weighted moments of a design and a response; the factorization with $X'y$ from the same scan; the pivoted Cholesky factor with `lm()`'s rule and the recheck of suspect columns],
  [`design_subset()`], [The rows of a design where a column is true, with the columns that go with the rows; with `keep = NULL`, no row is dropped and the design is moved to the frame of the columns computed from it],
  [`design_key()`, `is.na()` of a key, `key_codes()`], [Keys that can have missing values (of any type: `is.na()` of a stored column that is not double is `IS NULL`), and bare group codes: remotely the stored values themselves, so that the database groups without a join],
  [`group_sums()`, `group_crosstab()`, `group_lookup()`, `group_nested()`, `design_group_sums()`], [Sums, crosstabs, lookups and nesting of groupings or bare codes; the grouped sums of a design for several groupings (bare codes or groupings of `design_groups()`) in one pass],
  [`col_math()`, `col_pmin()`, `col_pmax(checked = FALSE)`], [The functions of columns that base R cannot dispatch on: distributions, and unchecked forms for columns known to have no missing values],
)

The estimators, the IRLS loop, the checks and messages of `lm()` and `glm()`, the covariance matrices, the small sample corrections and the methods of the fits are all in `xtx.fit`.

= Generalized linear models

`glm_xtx()` runs the iterations of `glm.fit()` step for step: the same starting values, the same weighted least squares problem at each iteration (`qr()` of the design scaled by the working weights), the same step halving, convergence rule and warnings (`R/glm.R` of `xtx.fit`). The linear predictor, means, working weights and working responses are columns, and the deviance is a sum of a column, so the loop runs unchanged on remote data.

= Fixed effects

`lm_xtx(y ~ x | f1 + f2)` absorbs the fixed effects without demeaning (`R/absorb.R` and `src/fe_solver.c` of `xtx.fit`, on the grouped sums and crosstabs of `xtx`, `src/fe.c`; the fit in `R/fe.R`). Write the model $y = X b + F a + e$, with $F$ the dummies of the $K$ fixed effects and $W$ the weights. With

$ S = F' W X, quad s = F' W y, quad K = F' W F, $

the Frisch-Waugh-Lovell moments are

$ X' M X = X' W X - S' K^(-) S, quad X' M y = X' W y - S' K^(-) s. $

Every input is a sum over rows: $X' W X$ and $X' W y$ are products with the design, $S$ and $s$ are grouped sums, and $K$ has the level counts on its diagonal blocks and the crosstabs of pairs of fixed effects off the diagonal. The remaining work is the solve with $K$, whose size is the number of levels and not the number of rows:

- the largest fixed effect is eliminated analytically, since its block of $K$ is diagonal;
- the Schur complement over the other fixed effects is a sparse, symmetric, singular system whose diagonal blocks are weighted graph Laplacians. One level per connected component is set to zero;
- it is solved by a pivoted dense Cholesky factorization when small (which also gives the rank), and otherwise by conjugate gradients: on the formed matrix with an incomplete Cholesky preconditioner when mobility between levels is low and the matrix is sparse, or without forming it, with a diagonal preconditioner, when it would be near dense. A probe of a few columns chooses between the two.

The same solve gives the values of the fixed effects. The coefficients are then refined from residuals computed row by row, which corrects the rounding of the one-pass sums. Cluster-robust scores are also sums: grouped column sums of $w e X$, minus a correction from the crosstab of clusters and fixed effects.

For a GLM with fixed effects, the weighted least squares step of each IRLS iteration is this estimator, with the working weights (`R/fe_glm.R`). The reduced system changes a little at each iteration, so conjugate gradients start from the previous solution.

Remotely, the database returns the grouped sums and crosstabs and R solves the reduced system. The fit scales with the number of levels and of non-empty crosstab cells, which must fit in memory, and not with the number of rows.

= Compared with other packages

*fixest* removes fixed effects by alternating projections: it demeans the response and every regressor, sweeping over the rows until convergence, then runs least squares on the demeaned data. The iterations here are on a system over the levels of the fixed effects, after a fixed number of passes over the rows. fixest works on the data directly, where this package goes through the normal equations and recovers accuracy by centering and refinement. fixest also does far more (several cluster variables, varying slopes, instrumental variables, and many other estimators); here, fits with fixed effects take `"iid"` or one cluster variable.

*dbreg* fits regressions in a database with a choice of strategies: compression by `GROUP BY`, moments ($X'X$ and $X'y$ in SQL, without fixed effects), a within transformation, with alternating projections in SQL when double demeaning is not exact (unbalanced or weighted two-way panels, three or more fixed effects), and Mundlak regressors. Here compression is an optimization inside one algorithm, fixed effects are solved in R from aggregates in every case, and the same code covers GLMs. The design itself is also translated: factors, contrasts, interactions and aliased columns follow `model.matrix()` and `lm()` on remote data as in memory.

= What this costs

- Estimates agree with `lm()`, `glm()` and fixest up to rounding, not bit for bit. Only the dense matrix is exact.
- Least squares goes through $X'X$. Centering, the check of suspect columns against the data, and refinement cover the usual cases, with a dense fallback at the extremes of scale.
- Memory is of order $p^2$ for the moments, plus the levels and crosstabs of the fixed effects. A design with very many columns outside the fixed effects is not the intended use.
- A feature that cannot be computed from sums on both backends is left out.

#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: [Architecture])
#metadata((title: "Architecture")) <website-metadata>

#calepin.setup(echo: true, eval: false)

#title()

This page describes how the package computes, for readers who want to know why the results are what they are, where the time goes, or where to look in the source. File names refer to the #link("https://github.com/vincentarelbundock/design.matrix")[repository].

= The idea

A model matrix is highly redundant. A factor with 1,000 levels fills 999 columns, but each row holds one integer code. An interaction of a factor and a numeric variable fills as many columns again, but each row holds one code and one number. The package keeps that compact description (factor codes, numeric variables, contrast matrices and the layout of the columns) and computes with it.

Three choices follow from this, and the rest of the page describes them:

+ *A design is a description of the matrix.* The dense matrix is built only by `as.matrix()`, and then bit for bit as `model.matrix()` builds it.
+ *Estimation is a sequence of sums over rows.* $X'X$, $X'y$, grouped column sums and crosstabs are computed from codes and numbers. Linear models, GLMs and fixed effects are fitted from these sums: the rows are never transformed (no demeaning, no dense $Q$).
+ *The fitting code is written once.* It calls a small set of internal primitives, each with a method for data in memory (C) and for remote data (SQL). There is no second implementation for out-of-core data.

= Layers

#table(
  columns: 3,
  table.header([Layer], [Role], [Source]),
  [Frame and plan], [The model frame of R, and the fixed transformations, levels and contrasts that new data must reuse], [`R/frame.R`, `R/prep.R`, `R/plan.R`],
  [Design], [Checks, contrasts and column names as in R; the dense fill], [`src/modelmatrix.c`, `R/design_matrix.R`],
  [Kernels], [Products with the design from codes and numbers], [`src/kernels.c`, `src/structured.c`, `R/algebra.R`],
  [Primitives], [Operations on designs and columns, in memory and remote], [`R/primitives.R`, `R/primitives_lookup.R`],
  [Fits], [Least squares, IRLS, fixed effects], [`R/qr.R`, `R/lm.R`, `R/glm.R`, `R/fe.R`, `R/fe_glm.R`, `src/fe.c`, `src/fused.c`],
  [Remote], [Designs over DuckDB tables, queries and Parquet files], [`R/remote.R`, `R/remote_sql.R`, `R/remote_fit.R`, `R/remote_compress.R`],
)

Each layer uses only the ones above it in the table. The fits do not know which backend they run on.

= A design is a description

== R's own preparation

`model.matrix()` and `model.frame()` do a good deal of work in R before they reach C: non-standard evaluation, character columns turned into factors, contrasts, `xlev`, `na.action`. This code changes between R versions. The package does not copy it. At load time it takes the installed `stats:::model.matrix.default()`, cuts it just before its call to C, and makes it return what that call would have received (`R/prep.R`). The same is done for `model.frame.default()`, with two patches that skip copies which change nothing (`R/frame.R`). Each patch is applied only when the installed code is the code it was written for, and a self-test at load drops a patch that fails; the fallback is the function of stats.

== Setup and fill

R's C routine `modelmatrix()` is ported and split in two (`src/modelmatrix.c`):

- `dm_setup()` performs every check, evaluates the contrasts, emits the warnings and computes the column names, in the order of R;
- `dm_fill()` allocates the matrix and fills it.

A `design_matrix` is the result of the setup plus the variables. The fill performs the same floating point operations as R, cell by cell. Only the order in which rows are visited changes (chunks of rows, on several threads when asked), so `as.matrix(D)` is identical to `model.matrix()`, including `NA`, `NaN`, `Inf` and the sign of zero.

== Plans

A `design_plan` keeps what a fitted design fixes for later data: the terms, the data-dependent transformations (`poly()`, `ns()`, `scale()`) with their training parameters, the factor levels and the contrasts (`R/plan.R`). `design_matrix(plan, batch)` applies it to new rows without repeating the checks. Remote designs carry a plan too: it fixes their terms, levels and contrasts, so that their columns are those of `model.matrix()`.

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

Both can run on several threads (`design_threads()`; the default is one). Results with one value per row do not depend on the number of threads. Sums over rows split the rows into fixed blocks, one per thread, added in a fixed order: they are identical for a given number of threads, and can differ in the last digits between numbers of threads.

== Patterns

Rows with the same factor codes and the same numeric values give the same design row. When a design has far fewer such patterns than rows, the kernels compute once per pattern and expand, or sum the weights within patterns first (`R/algebra.R`, option `design.matrix.dedup`). The fits use the same idea as compression: a weighted least squares fit on the patterns, with the sum of the weights and the weighted mean response of each, is the fit on the rows, and the within-pattern sum of squares completes the residual sum of squares (`R/lm.R`). For a GLM, the iterations on the patterns are those on the rows for every supported family (`R/glm.R`).

= Least squares

`qr()` of a design (`R/qr.R`) decomposes $X P = Q R$ without the dense $X$ or the $n times p$ matrix $Q$:

- $R$ is the Cholesky factor of $X'X$, from `crossprod()`;
- columns are pivoted with the rule of `lm()` (LINPACK's limited pivoting). From $X'X$ a residual norm is accurate to about the square root of machine precision, so the columns near the tolerance are judged from the data, with one product per block of columns;
- coefficients solve the normal equations and are refined from the residuals $y - X b$, each step costing one $X b$ and one $X' r$;
- numeric variables far from zero are centered first when that is an exact linear transform of the design, which keeps $X'X$ accurate.

The cost is a few passes over the rows and memory of order $p^2$. The price is that of the normal equations: on badly conditioned designs the coefficients agree with a Householder QR to fewer digits, and when $X'X$ cannot represent the scale of the columns `qr()` falls back to the dense decomposition.

= One fitting code, two backends

== The primitive contract

`R/primitives.R` defines the operations that fitting code may use. They take a *design* and *columns* (one value per row):

#table(
  columns: 3,
  table.header([], [In memory], [Remote]),
  [Design], [`design_matrix`], [`dm_remote`: a source, a plan, and how to compute each variable in SQL],
  [Column], [numeric vector or matrix], [SQL expressions over the rows of a frame],
)

The operations are products with the design ($X'X$, $X'c$, $X B$, $op("rowSums")((X V) circle.small X)$), sums and grouped sums of columns, crosstabs of two groupings, lookups of a value per level, subsets, and arithmetic on columns (`+`, `*`, `exp()`, `log()`, comparisons, the functions the GLM families need). Arithmetic on a remote column builds SQL; nothing runs until a sum is asked for. Missing values follow R: SQL `NULL` is `NA`, and sums return `NA` where `sum()` does.

`lm.R`, `glm.R`, `fe.R` and `fe_glm.R` contain no branch on the backend. A fit on a DuckDB table is the same sequence of operations as in memory, which is why its coefficients, aliased columns, covariance and degrees of freedom are those of `lm()` or `glm()` on the same rows.

== The SQL engine

For a remote design, the products of the previous section are `GROUP BY` aggregates over factor columns as stored. Levels become codes in R, on the cells that come back, and the $C_t$ are read off the in-memory design of a few prototype rows (every combined level, numeric variables set to 1), so that columns, contrasts and their order are those of `model.matrix()`.

A variable that SQL cannot compute (an R function with no translation) is evaluated in R once, batch by batch, into a column of a temporary table, and is then a column like any other. Results with one value per row (fitted values, residuals, predictions) are views in the database, which keep the columns of the data.

== Compression in the database

When the data have few design patterns, the scan of the SQL engine also groups by the numeric variables of the design, and returns each cell with its number of rows and the sums of the weights and responses; the design of the cells is built in memory from the engine's own $C_t$ and the fit runs there (`R/remote_compress.R`). An estimate of the number of patterns, from a sample and from `approx_count_distinct()`, decides whether to try.

= Generalized linear models

`design_glm()` runs the iterations of `glm.fit()` step for step: the same starting values, the same weighted least squares problem at each iteration (`qr()` of the design scaled by the working weights), the same step halving, convergence rule and warnings (`R/glm.R`). The linear predictor, means, working weights and working responses are columns, and the deviance is a sum of a column, so the loop runs unchanged on remote data.

= Fixed effects

`design_lm(y ~ x | f1 + f2)` absorbs the fixed effects without demeaning (`R/fe.R`, `src/fe.c`). Write the model $y = X b + F a + e$, with $F$ the dummies of the $K$ fixed effects and $W$ the weights. With

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

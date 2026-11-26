# The input stage of the formula fits (design_lm(), design_glm(), fixed
# effects; in memory and on a DuckDB table): every error and warning
# message, exactly (fixed = TRUE), and the error reported when a call has
# two problems (validation order). These pin current behaviour, including
# differences between paths that are not decided yet:
# * offsets: lm() says "is %d, should equal", glm() and FE "is %d should equal";
# * weights: LM "missing or negative weights not allowed", GLM "negative
#   weights not allowed";
# * a two-column response, a subset leaving no rows: other errors remotely
#   than in memory;
# * na.action = na.pass with FE: fits in memory, refused remotely.
if (!requireNamespace("duckdb", quietly = TRUE) || !requireNamespace("DBI", quietly = TRUE))
    exit_file("duckdb is not installed")
library(design.matrix)

con <- DBI::dbConnect(duckdb::duckdb())
set.seed(1)
n <- 40
d <- data.frame(y = rnorm(n), x = rnorm(n), z = runif(n, 1, 2),
                f = factor(sample(c("a", "b", "c"), n, TRUE)), ch = sample(c("u", "v"), n, TRUE),
                id = rep(1:8, each = 5), cl = rep(1:4, 10),
                k = rpois(n, 2), b1 = rbinom(n, 1, 0.5), w = runif(n, 0.5, 2),
                wneg = c(-1, runif(n - 1)), wzero = 0, stringsAsFactors = FALSE)
d$yna <- d$y
d$yna[3] <- NA
d$wna <- d$w
d$wna[5] <- NA
d$fails <- 3 - d$b1
DBI::dbWriteTable(con, "d", d)

# LM memory
expect_error(design_lm(y ~ x, d, weights = wneg),
             "missing or negative weights not allowed", fixed = TRUE)
expect_error(design_lm(y ~ x, d, weights = wna, na.action = na.pass),
             "missing or negative weights not allowed", fixed = TRUE)
expect_error(design_lm(y ~ x, d, weights = ch),
             "'weights' must be a numeric vector", fixed = TRUE)
expect_error(design_lm(y ~ x, d, weights = wzero),
             "0 (non-NA) cases", fixed = TRUE)
expect_error(design_lm(y ~ x, d, offset = 1:3),
             "variable lengths differ (found for '(offset)')", fixed = TRUE)
expect_error(design_lm(~ x, d, offset = z),
             "The formula has no response.", fixed = TRUE)
expect_error(design_lm(yna ~ x, d, na.action = na.pass),
             "NA/NaN/Inf in 'y'", fixed = TRUE)
expect_error(design_lm(yna ~ x, d, na.action = na.fail),
             "missing values in object", fixed = TRUE)
expect_error(design_lm(ch ~ x, d),
             "The response must be numeric.", fixed = TRUE)
expect_error(design_lm(cbind(y, x) ~ z, d),
             "Multivariate responses are not supported.", fixed = TRUE)
expect_error(design_lm(~ x, d),
             "The formula has no response.", fixed = TRUE)
expect_error(design_lm(y ~ x, d, subset = x > 100),
             "0 (non-NA) cases", fixed = TRUE)
expect_error(design_lm(y ~ x, d, subset = nosuch > 0),
             "object 'nosuch' not found", fixed = TRUE)
expect_error(design_lm(y ~ x, d, foo = 1),
             "Unknown arguments: foo", fixed = TRUE)
expect_error(design_lm(y ~ x, d, foo = 1, weights = wneg),
             "Unknown arguments: foo", fixed = TRUE)
expect_error(design_lm(ch ~ x, d, weights = wneg),
             "The response must be numeric.", fixed = TRUE)
expect_error(design_lm(cbind(y, x) ~ z, d, weights = ch),
             "Multivariate responses are not supported.", fixed = TRUE)
expect_error(design_lm(y ~ x, d, vcov = "hc1", compress = "no"),
             "`vcov` other than \"iid\" requires fixed effects (y ~ x | fe).", fixed = TRUE)
expect_error(design_lm(y ~ x, d, y = 1, compress = "no"),
             "`y` applies only when the first argument is a design matrix.", fixed = TRUE)
expect_error(design_lm(y ~ x, d, compress = "no", tol = -1),
             "`compress` must be \"auto\", TRUE or FALSE.", fixed = TRUE)
expect_error(design_lm(y ~ x, d, weights = wneg, subset = x > 100),
             "0 (non-NA) cases", fixed = TRUE)
expect_error(design_lm(y ~ x, d, na.action = function(x) stop("custom na")),
             "custom na", fixed = TRUE)

# LM remote
expect_error(design_lm(y ~ x, con, table = "d", weights = wneg),
             "missing or negative weights not allowed", fixed = TRUE)
expect_error(design_lm(y ~ x, con, table = "d", weights = wna, na.action = na.pass),
             "missing or negative weights not allowed", fixed = TRUE)
expect_error(design_lm(y ~ x, con, table = "d", weights = ch),
             "`weights` has no SQL translation: column `ch` is not numeric.", fixed = TRUE)
expect_error(design_lm(y ~ x, con, table = "d", weights = wzero),
             "0 (non-NA) cases", fixed = TRUE)
expect_error(design_lm(~ x, con, table = "d", offset = z),
             "The formula has no response.", fixed = TRUE)
expect_error(design_lm(yna ~ x, con, table = "d", na.action = na.pass),
             "NA/NaN/Inf in 'y'", fixed = TRUE)
expect_error(design_lm(yna ~ x, con, table = "d", na.action = na.fail),
             "missing values in object", fixed = TRUE)
expect_error(design_lm(ch ~ x, con, table = "d"),
             "`response` has no SQL translation: column `ch` is not numeric.", fixed = TRUE)
expect_error(design_lm(cbind(y, x) ~ z, con, table = "d"),
             "`cbind(y, x)` gives 80 values for 40 rows.", fixed = TRUE)
expect_error(design_lm(~ x, con, table = "d"),
             "The formula has no response.", fixed = TRUE)
expect_error(design_lm(y ~ x, con, table = "d", subset = x > 100),
             "The remote data has no complete rows.", fixed = TRUE)
expect_error(design_lm(y ~ x, con, table = "d", subset = nosuch > 0),
             "object 'nosuch' not found", fixed = TRUE)
expect_error(design_lm(y ~ x, con, table = "d", foo = 1),
             "Unknown arguments: foo", fixed = TRUE)
expect_error(design_lm(y ~ x, con, table = "d", foo = 1, weights = wneg),
             "Unknown arguments: foo", fixed = TRUE)
expect_error(design_lm(y ~ x, con, table = "d", na.action = function(x) x),
             "For remote data, `na.action` must be na.omit, na.exclude, na.fail, na.pass.", fixed = TRUE)
expect_error(design_lm(y ~ x, con, table = "d", na.action = function(x) x, foo = 1),
             "For remote data, `na.action` must be na.omit, na.exclude, na.fail, na.pass.", fixed = TRUE)
expect_error(design_lm(~ x, con, table = "d", foo = 1),
             "The formula has no response.", fixed = TRUE)
expect_error(design_lm(y ~ x, con, table = "d", weights = wneg, subset = x > 100),
             "The remote data has no complete rows.", fixed = TRUE)
expect_error(design_lm(y ~ x, con, table = "d", vcov = "hc1", compress = "no"),
             "`vcov` other than \"iid\" requires fixed effects (y ~ x | fe).", fixed = TRUE)

# GLM memory
expect_error(design_glm(b1 ~ x, binomial, d, weights = wneg),
             "negative weights not allowed", fixed = TRUE)
expect_error(design_glm(y ~ x, gaussian, d, weights = wna, na.action = na.fail),
             "missing values in object", fixed = TRUE)
expect_error(design_glm(y ~ x, gaussian, d, weights = ch),
             "'weights' must be a numeric vector", fixed = TRUE)
expect_error(suppressWarnings(design_glm(y ~ x, gaussian, d, weights = wzero)),
             "The fit failed: no valid coefficients were found.", fixed = TRUE)
expect_warning(design_glm(y ~ x, gaussian, d, weights = wzero),
               "no observations informative at iteration 1", fixed = TRUE)
expect_error(design_glm(y ~ x, gaussian, d, offset = 1:3),
             "variable lengths differ (found for '(offset)')", fixed = TRUE)
expect_error(design_glm(~ x, gaussian, d, offset = z),
             "The formula has no response.", fixed = TRUE)
expect_error(design_glm(~ x, gaussian, d),
             "The formula has no response.", fixed = TRUE)
expect_error(design_glm(yna ~ x, gaussian, d, na.action = na.pass),
             "`na.action` must be na.omit, na.exclude, na.fail.", fixed = TRUE)
expect_error(design_glm(yna ~ x, gaussian, d, na.action = na.fail),
             "missing values in object", fixed = TRUE)
expect_error(design_glm(f ~ x, poisson, d),
             "A factor response needs a binomial family.", fixed = TRUE)
expect_error(design_glm(design_matrix(~ x, d), poisson, y = d$f),
             "A factor response needs a binomial family.", fixed = TRUE)
expect_error(design_glm(ch ~ x, poisson, d),
             "A factor response needs a binomial family.", fixed = TRUE)
expect_error(design_glm(cbind(b1, fails, k) ~ x, binomial, d),
             "A two-column response needs two columns (successes and failures).", fixed = TRUE)
expect_error(design_glm(cbind(b1, fails) ~ x, poisson, d),
             "A two-column response (successes and failures) needs a binomial family.", fixed = TRUE)
expect_error(design_glm(as.Date(k) ~ x, gaussian, d),
             "The response must be numeric, logical or a factor.", fixed = TRUE)
expect_error(design_glm(y ~ x, gaussian, d, subset = x > 0),
             "object 'x' not found", fixed = TRUE)
expect_error(design_glm(y ~ x, gaussian, d, foo = 1),
             "Unknown arguments: foo", fixed = TRUE)
expect_error(design_glm(y ~ x, gaussian, d, foo = 1, weights = wneg),
             "Unknown arguments: foo", fixed = TRUE)
expect_error(design_glm(cbind(b1, fails, k) ~ x, binomial, d, weights = ch),
             "'weights' must be a numeric vector", fixed = TRUE)
expect_error(design_glm(f ~ x, poisson, d, weights = wneg),
             "A factor response needs a binomial family.", fixed = TRUE)
expect_error(design_glm(y ~ x, gaussian, d, vcov = "hc1", compress = "no"),
             "`compress` must be \"auto\", TRUE or FALSE.", fixed = TRUE)
expect_error(design_glm(y ~ x, gaussian, d, y = 1, na.action = na.pass),
             "`y` applies only when the first argument is a design matrix.", fixed = TRUE)
expect_error(design_glm(y ~ x, gaussian, d, na.action = na.pass, foo = 1),
             "`na.action` must be na.omit, na.exclude, na.fail.", fixed = TRUE)
expect_error(design_glm(y ~ x, "nofamily", d),
             "object 'nofamily' of mode 'function' was not found", fixed = TRUE)
expect_error(design_glm(y ~ x, gaussian, d, start = 1),
             "length of 'start' should equal 2 and correspond to initial coefs for c(\"(Intercept)\", \"x\")", fixed = TRUE)

# GLM remote
expect_error(design_glm(b1 ~ x, binomial, con, table = "d", weights = wneg),
             "negative weights not allowed", fixed = TRUE)
expect_error(design_glm(y ~ x, gaussian, con, table = "d", weights = wna, na.action = na.fail),
             "missing values in object", fixed = TRUE)
expect_error(design_glm(y ~ x, gaussian, con, table = "d", weights = ch),
             "`weights` has no SQL translation: column `ch` is not numeric.", fixed = TRUE)
expect_error(suppressWarnings(design_glm(y ~ x, gaussian, con, table = "d", weights = wzero)),
             "The fit failed: no valid coefficients were found.", fixed = TRUE)
expect_warning(design_glm(y ~ x, gaussian, con, table = "d", weights = wzero),
               "no observations informative at iteration 1", fixed = TRUE)
expect_error(design_glm(~ x, gaussian, con, table = "d", offset = z),
             "The formula has no response.", fixed = TRUE)
expect_error(design_glm(~ x, gaussian, con, table = "d"),
             "The formula has no response.", fixed = TRUE)
expect_error(design_glm(yna ~ x, gaussian, con, table = "d", na.action = na.pass),
             "`na.action` must be na.omit, na.exclude, na.fail.", fixed = TRUE)
expect_error(design_glm(yna ~ x, gaussian, con, table = "d", na.action = na.fail),
             "missing values in object", fixed = TRUE)
expect_error(design_glm(f ~ x, poisson, con, table = "d"),
             "A factor response needs a binomial family.", fixed = TRUE)
expect_error(design_glm(ch ~ x, poisson, con, table = "d"),
             "A factor response needs a binomial family.", fixed = TRUE)
expect_inherits(design_glm(f ~ x, binomial, con, table = "d"), "design_glm")
expect_error(design_glm(design_matrix(~ x, con, table = "d"), poisson, y = "f"),
             "A factor response needs a binomial family.", fixed = TRUE)
expect_error(design_glm(cbind(b1, fails, k) ~ x, binomial, con, table = "d"),
             "A two-column response needs two columns (successes and failures).", fixed = TRUE)
expect_error(design_glm(cbind(b1, fails) ~ x, poisson, con, table = "d"),
             "A two-column response (successes and failures) needs a binomial family.", fixed = TRUE)
expect_error(design_glm(y ~ x, gaussian, con, table = "d", subset = x > 0),
             "object 'x' not found", fixed = TRUE)
expect_error(design_glm(y ~ x, gaussian, con, table = "d", foo = 1),
             "Unknown arguments: foo", fixed = TRUE)
expect_error(design_glm(y ~ x, gaussian, con, table = "d", foo = 1, weights = wneg),
             "Unknown arguments: foo", fixed = TRUE)
expect_error(design_glm(cbind(b1, fails, k) ~ x, binomial, con, table = "d", weights = ch),
             "A two-column response needs two columns (successes and failures).", fixed = TRUE)
expect_error(design_glm(f ~ x, poisson, con, table = "d", weights = wneg),
             "A factor response needs a binomial family.", fixed = TRUE)
expect_error(design_glm(y ~ x, gaussian, con, table = "d", na.action = na.pass, foo = 1),
             "`na.action` must be na.omit, na.exclude, na.fail.", fixed = TRUE)
expect_error(design_glm(y ~ x, gaussian, con, table = "d", start = 1),
             "length of 'start' should equal 2 and correspond to initial coefs for c(\"(Intercept)\", \"x\")", fixed = TRUE)

# FE LM memory
expect_error(design_lm(y ~ x | id, d, weights = wneg),
             "missing or negative weights not allowed", fixed = TRUE)
expect_error(design_lm(y ~ x | id, d, weights = ch),
             "'weights' must be a numeric vector", fixed = TRUE)
expect_error(design_lm(y ~ x | id, d, offset = 1:3),
             "variable lengths differ (found for '(offset)')", fixed = TRUE)
expect_error(design_lm(~ x | id, d, offset = z),
             "The formula has no response.", fixed = TRUE)
expect_error(design_lm(~ x | id, d),
             "The formula has no response.", fixed = TRUE)
expect_error(design_lm(ch ~ x | id, d),
             "The response must be numeric.", fixed = TRUE)
expect_error(design_lm(cbind(y, x) ~ z | id, d),
             "Multivariate responses are not supported.", fixed = TRUE)
expect_error(design_lm(yna ~ x | id, d, na.action = na.fail),
             "missing values in object", fixed = TRUE)
expect_inherits(design_lm(yna ~ x | id, d, na.action = na.pass), "design_lm")
expect_error(design_lm(y ~ x | id, d, subset = x > 0),
             "`subset` is not supported with fixed effects.", fixed = TRUE)
expect_error(design_lm(y ~ x | id, d, compress = TRUE),
             "`compress = TRUE` is not supported with fixed effects.", fixed = TRUE)
expect_error(design_lm(y ~ x | id, d, foo = 1),
             "Unknown arguments: foo", fixed = TRUE)
expect_error(design_lm(y ~ x | id, d, foo = 1, weights = wneg),
             "Unknown arguments: foo", fixed = TRUE)
expect_error(design_lm(ch ~ x | id, d, weights = wneg),
             "The response must be numeric.", fixed = TRUE)
expect_error(design_lm(ch ~ x | id, d, weights = ch),
             "'weights' must be a numeric vector", fixed = TRUE)
expect_error(design_lm(y ~ x | id, d, subset = x > 0, compress = TRUE),
             "`compress = TRUE` is not supported with fixed effects.", fixed = TRUE)
expect_error(design_lm(y ~ x | id, d, vcov = ~nosuch),
             "object 'nosuch' not found", fixed = TRUE)
expect_error(design_lm(y ~ x | id, d, wzero),
             "`subset` is not supported with fixed effects.", fixed = TRUE)

# FE LM remote
expect_error(design_lm(y ~ x | id, con, table = "d", weights = wneg),
             "missing or negative weights not allowed", fixed = TRUE)
expect_error(design_lm(y ~ x | id, con, table = "d", weights = ch),
             "`weights` has no SQL translation: column `ch` is not numeric.", fixed = TRUE)
expect_error(design_lm(~ x | id, con, table = "d"),
             "The formula has no response.", fixed = TRUE)
expect_error(design_lm(ch ~ x | id, con, table = "d"),
             "`response` has no SQL translation: column `ch` is not numeric.", fixed = TRUE)
expect_error(design_lm(yna ~ x | id, con, table = "d", na.action = na.fail),
             "For remote data, `na.action` must be na.omit, na.exclude.", fixed = TRUE)
expect_error(design_lm(y ~ x | id, con, table = "d", na.action = na.fail, foo = 1),
             "For remote data, `na.action` must be na.omit, na.exclude.", fixed = TRUE)
expect_error(design_lm(y ~ x | id, con, table = "d", foo = 1),
             "Unknown arguments: foo", fixed = TRUE)
expect_error(design_lm(y ~ x | id, con, table = "d", foo = 1, weights = wneg),
             "Unknown arguments: foo", fixed = TRUE)
expect_error(design_lm(y ~ x | id, con, table = "d", subset = x > 0),
             "`subset` is not supported with fixed effects.", fixed = TRUE)
expect_error(design_lm(y ~ x | nosuch, con, table = "d"),
             "No column `nosuch` in the remote data.", fixed = TRUE)
expect_error(design_lm(y ~ x | id, con, table = "d", vcov = ~nosuch),
             "No column `nosuch` in the remote data.", fixed = TRUE)

# FE GLM memory
expect_error(design_glm(b1 ~ x | id, binomial, d, weights = wneg),
             "negative weights not allowed", fixed = TRUE)
expect_error(design_glm(k ~ x | id, poisson, d, weights = ch),
             "'weights' must be a numeric vector", fixed = TRUE)
expect_error(design_glm(k ~ x | id, poisson, d, offset = 1:3),
             "variable lengths differ (found for '(offset)')", fixed = TRUE)
expect_error(design_glm(~ x | id, poisson, d),
             "The formula has no response.", fixed = TRUE)
expect_error(design_glm(f ~ x | id, poisson, d),
             "A factor response needs a binomial family.", fixed = TRUE)
expect_error(design_glm(cbind(b1, fails, k) ~ x | id, binomial, d),
             "A two-column response needs two columns (successes and failures).", fixed = TRUE)
expect_error(design_glm(k ~ x | id, poisson, d, na.action = na.fail),
             "With fixed effects, `na.action` must be na.omit or na.exclude.", fixed = TRUE)
expect_error(design_glm(k ~ x | id, poisson, d, na.action = na.pass),
             "`na.action` must be na.omit, na.exclude, na.fail.", fixed = TRUE)
expect_error(design_glm(k ~ x | id, poisson, d, start = 1),
             "`start` is not supported with fixed effects.", fixed = TRUE)
expect_error(design_glm(k ~ x | id, poisson, d, foo = 1),
             "Unknown arguments: foo", fixed = TRUE)
expect_error(design_glm(k ~ x | id, poisson, d, subset = x > 0),
             "object 'x' not found", fixed = TRUE)
expect_error(design_glm(k ~ x | id, poisson, d, compress = TRUE),
             "`compress = TRUE` is not supported with fixed effects.", fixed = TRUE)
expect_error(design_glm(k ~ x | id, poisson, d, foo = 1, weights = wneg),
             "Unknown arguments: foo", fixed = TRUE)
expect_error(design_glm(f ~ x | id, poisson, d, weights = ch),
             "'weights' must be a numeric vector", fixed = TRUE)

# FE GLM remote
expect_error(design_glm(b1 ~ x | id, binomial, con, table = "d", weights = wneg),
             "negative weights not allowed", fixed = TRUE)
expect_error(design_glm(k ~ x | id, poisson, con, table = "d", weights = ch),
             "`weights` has no SQL translation: column `ch` is not numeric.", fixed = TRUE)
expect_error(design_glm(~ x | id, poisson, con, table = "d"),
             "The formula has no response.", fixed = TRUE)
expect_error(design_glm(f ~ x | id, poisson, con, table = "d"),
             "A factor response needs a binomial family.", fixed = TRUE)
expect_error(design_glm(cbind(b1, fails, k) ~ x | id, binomial, con, table = "d"),
             "A two-column response needs two columns (successes and failures).", fixed = TRUE)
expect_error(design_glm(k ~ x | id, poisson, con, table = "d", na.action = na.fail),
             "With fixed effects, `na.action` must be na.omit or na.exclude.", fixed = TRUE)
expect_error(design_glm(k ~ x | id, poisson, con, table = "d", foo = 1),
             "Unknown arguments: foo", fixed = TRUE)
expect_error(design_glm(k ~ x | id, poisson, con, table = "d", subset = x > 0),
             "object 'x' not found", fixed = TRUE)
expect_error(design_glm(k ~ x | id, poisson, con, table = "d", foo = 1, weights = wneg),
             "Unknown arguments: foo", fixed = TRUE)
expect_error(design_glm(f ~ x | id, poisson, con, table = "d", weights = ch),
             "`weights` has no SQL translation: column `ch` is not numeric.", fixed = TRUE)
# Negative weights on a design matrix: the wording of lm() or glm()
Dm <- design_matrix(~ x, d)
Dr <- design_matrix(~ x, con, table = "d")
expect_error(design_lm(Dm, y = d$y, weights = d$wneg),
             "missing or negative weights not allowed", fixed = TRUE)
expect_error(design_lm(Dr, y = "y", weights = wneg),
             "missing or negative weights not allowed", fixed = TRUE)
expect_error(design_glm(Dm, binomial, y = d$b1, weights = d$wneg),
             "negative weights not allowed", fixed = TRUE)
expect_error(design_glm(Dr, binomial, y = "b1", weights = wneg),
             "negative weights not allowed", fixed = TRUE)
expect_false(grepl("missing", tryCatch(design_glm(Dm, binomial, y = d$b1, weights = d$wneg),
                                       error = conditionMessage)))
expect_false(grepl("missing", tryCatch(design_glm(Dr, binomial, y = "b1", weights = wneg),
                                       error = conditionMessage)))

# Offsets of the wrong length (a matrix offset), two problems, na.action,
# contrasts, non-numeric offsets and subsets, zero weights
expect_error(design_lm(y ~ x, d, offset = cbind(z, z)),
             "number of offsets is 80, should equal 40 (number of observations)", fixed = TRUE)
expect_error(design_glm(y ~ x, gaussian, d, offset = cbind(z, z)),
             "number of offsets is 80 should equal 40 (number of observations)", fixed = TRUE)
expect_error(design_lm(y ~ x | id, d, offset = cbind(z, z)),
             "number of offsets is 80 should equal 40 (number of observations)", fixed = TRUE)
expect_error(design_glm(k ~ x | id, poisson, d, offset = cbind(z, z)),
             "number of offsets is 80 should equal 40 (number of observations)", fixed = TRUE)
expect_error(design_lm(y ~ x, con, table = "d", offset = cbind(z, z)),
             "number of offsets is 80, should equal 40 (number of observations)", fixed = TRUE)
expect_error(design_lm(y ~ x + offset(cbind(z, z)), con, table = "d"),
             "number of offsets is 80, should equal 40 (number of observations)", fixed = TRUE)
expect_error(design_lm(y ~ x, con, table = "d", subset = x > 0, offset = cbind(z, z)),
             sprintf("number of offsets is %d, should equal %d (number of observations)",
                     2L * sum(d$x > 0), sum(d$x > 0)), fixed = TRUE)
expect_error(design_lm(y ~ x, d, subset = x > 0, offset = cbind(z, z)),
             sprintf("number of offsets is %d, should equal %d (number of observations)",
                     2L * sum(d$x > 0), sum(d$x > 0)), fixed = TRUE)
expect_error(design_glm(y ~ x, gaussian, con, table = "d", offset = cbind(z, z)),
             "number of offsets is 80 should equal 40 (number of observations)", fixed = TRUE)
expect_error(design_lm(y ~ x | id, con, table = "d", offset = cbind(z, z)),
             "number of offsets is 80 should equal 40 (number of observations)", fixed = TRUE)
expect_error(design_glm(k ~ x | id, poisson, con, table = "d", offset = cbind(z, z)),
             "number of offsets is 80 should equal 40 (number of observations)", fixed = TRUE)
expect_error(design_glm(cbind(b1, fails) ~ x, binomial, con, table = "d", offset = cbind(z, z)),
             "number of offsets is 80 should equal 40 (number of observations)", fixed = TRUE)
expect_error(design_lm(y ~ x, d, offset = cbind(z, z), weights = ch),
             "'weights' must be a numeric vector", fixed = TRUE)
expect_error(design_glm(y ~ x, gaussian, d, offset = cbind(z, z), weights = ch),
             "'weights' must be a numeric vector", fixed = TRUE)
expect_error(design_glm(cbind(b1, fails, k) ~ x, binomial, d, offset = cbind(z, z)),
             "number of offsets is 80 should equal 40 (number of observations)", fixed = TRUE)
expect_error(design_lm(ch ~ x, d, offset = cbind(z, z)),
             "The response must be numeric.", fixed = TRUE)
expect_error(design_lm(ch ~ x | id, d, offset = cbind(z, z)),
             "number of offsets is 80 should equal 40 (number of observations)", fixed = TRUE)
expect_error(design_lm(y ~ x, d, weights = wneg, offset = cbind(z, z)),
             "number of offsets is 80, should equal 40 (number of observations)", fixed = TRUE)
expect_inherits(design_lm(y ~ x, d, na.action = na.exclude), "design_lm")
expect_inherits(design_lm(y ~ x, con, table = "d", na.action = "na.omit"), "design_lm")
expect_inherits(suppressWarnings(design_lm(y ~ x, d, contrasts = list(f = "nosuch"))), "design_lm")
expect_warning(design_lm(y ~ x, d, contrasts = list(f = "nosuch")),
               "variable 'f' is absent, its contrast will be ignored", fixed = TRUE)
expect_inherits(suppressWarnings(design_lm(y ~ x, con, table = "d", contrasts = list(nosuch = "contr.sum"))), "design_lm")
expect_warning(design_lm(y ~ x, con, table = "d", contrasts = list(nosuch = "contr.sum")),
               "variable 'nosuch' is absent, its contrast will be ignored", fixed = TRUE)
expect_inherits(suppressWarnings(design_glm(y ~ x, gaussian, d, contrasts = list(nosuch = "contr.sum"))), "design_glm")
expect_warning(design_glm(y ~ x, gaussian, d, contrasts = list(nosuch = "contr.sum")),
               "variable 'nosuch' is absent, its contrast will be ignored", fixed = TRUE)
expect_inherits(suppressWarnings(design_glm(y ~ x, gaussian, con, table = "d", contrasts = list(nosuch = "contr.sum"))), "design_glm")
expect_warning(design_glm(y ~ x, gaussian, con, table = "d", contrasts = list(nosuch = "contr.sum")),
               "variable 'nosuch' is absent, its contrast will be ignored", fixed = TRUE)
expect_inherits(suppressWarnings(design_lm(y ~ x | id, d, contrasts = list(nosuch = "contr.sum"))), "design_lm")
expect_warning(design_lm(y ~ x | id, d, contrasts = list(nosuch = "contr.sum")),
               "variable 'nosuch' is absent, its contrast will be ignored", fixed = TRUE)
expect_inherits(suppressWarnings(design_lm(y ~ x | id, con, table = "d", contrasts = list(nosuch = "contr.sum"))), "design_lm")
expect_warning(design_lm(y ~ x | id, con, table = "d", contrasts = list(nosuch = "contr.sum")),
               "variable 'nosuch' is absent, its contrast will be ignored", fixed = TRUE)
expect_error(design_lm(y ~ x, d, offset = ch),
             "'offset' must be numeric", fixed = TRUE)
expect_error(design_glm(y ~ x, gaussian, d, offset = ch),
             "'offset' must be numeric", fixed = TRUE)
expect_error(design_lm(y ~ x, con, table = "d", offset = ch),
             "`offset` has no SQL translation: column `ch` is not numeric.", fixed = TRUE)
expect_error(design_lm(y ~ x | id, d, weights = wzero),
             "0 (non-NA) cases", fixed = TRUE)
expect_error(design_lm(y ~ x | id, con, table = "d", weights = wzero),
             "0 (non-NA) cases", fixed = TRUE)
expect_error(design_lm(y ~ x | id, d, na.action = na.exclude, foo = 1),
             "Unknown arguments: foo", fixed = TRUE)
expect_error(design_glm(k ~ x | id, poisson, d, weights = wzero),
             "0 (non-NA) cases", fixed = TRUE)

# What a fit keeps from its input, per path (names, not order).
d2 <- d
d2$k[4] <- NA
DBI::dbWriteTable(con, "d2", d2)
fit_names <- function(fit) sort(names(fit))
f <- design_lm(yna ~ x + f, d2, weights = w, offset = z, na.action = na.exclude)
expect_identical(fit_names(f), sort(c("coefficients", "rank", "pivot", "R", "tol", "rss", "mss", "intercept",
    "df.residual", "sigma", "nobs", "sumlogw", "weighted", "compressed", "patterns", "rows", "spec", "call",
    "terms", "contrasts", "xlevels", "na.action", "offset_call", "env", "model", "data")))
expect_inherits(f$na.action, "exclude")
expect_identical(f$offset_call, quote(z))
expect_identical(f$env, environment())
expect_null(f$method)
expect_identical(names(f$data), c("design", "y", "weights", "offset"))
expect_inherits(f$model, "data.frame")
f <- design_lm(yna ~ x + f, d2, weights = w, offset = z, na.action = na.exclude, model = FALSE)
expect_null(f$data)
expect_null(f$model)
f <- design_lm(yna ~ x + f, con, table = "d2", weights = w, offset = z, na.action = na.exclude)
expect_identical(fit_names(f), sort(c("coefficients", "rank", "pivot", "R", "tol", "rss", "mss", "intercept",
    "df.residual", "sigma", "nobs", "sumlogw", "weighted", "compressed", "patterns", "rows", "spec", "method",
    "contrasts", "xlevels", "data", "na.action", "offset_call", "env", "call", "terms")))
expect_inherits(f$na.action, "dm_remote_exclude")
expect_identical(f$offset_call, quote(z))
expect_identical(f$method, "SQL moments")
expect_identical(names(f$data), c("design", "y", "weights", "offset"))
expect_identical(f$spec, f$data$design$plan)
f <- design_lm(yna ~ f, con, table = "d2", model = FALSE)
expect_identical(f$method, "patterns")
expect_false(is.null(f$data))

f <- design_glm(k ~ x + f, poisson, d2, weights = w, offset = log(z), na.action = na.exclude)
expect_identical(fit_names(f), sort(c("coefficients", "rank", "pivot", "R", "tol", "deviance", "null.deviance",
    "aic", "iter", "converged", "boundary", "df.residual", "df.null", "nobs", "rows", "pearson", "zero_weights",
    "family", "intercept", "compressed", "patterns", "control", "spec", "method", "data", "call", "terms",
    "formula", "contrasts", "xlevels", "na.action", "offset_call", "env", "model")))
expect_inherits(f$na.action, "exclude")
expect_identical(f$offset_call, quote(log(z)))
expect_identical(f$formula, k ~ x + f)
f <- design_glm(k ~ x + f, poisson, d2, model = FALSE)
expect_null(f$data)
expect_null(f$model)
f <- design_glm(k ~ x + f, poisson, con, table = "d2", weights = w, offset = log(z), na.action = na.exclude)
expect_identical(fit_names(f), sort(c("coefficients", "rank", "pivot", "R", "tol", "deviance", "null.deviance",
    "aic", "iter", "converged", "boundary", "df.residual", "df.null", "nobs", "rows", "pearson", "zero_weights",
    "family", "intercept", "compressed", "patterns", "control", "spec", "method", "data", "na.action", "call",
    "terms", "formula", "contrasts", "xlevels", "offset_call", "env")))
expect_inherits(f$na.action, "dm_remote_exclude")
expect_identical(f$offset_call, quote(log(z)))
expect_identical(f$formula, k ~ x + f)

f <- design_lm(yna ~ x | id, d2, weights = w, offset = z, na.action = na.exclude)
expect_identical(fit_names(f), sort(c("coefficients", "vcov", "vcov_type", "sigma", "rss", "df.residual",
    "nobs", "rank", "fixef", "fe_levels", "fe_vars", "fe_sizes", "fe_dof", "fe_dof_exact", "nsingletons",
    "cov.unscaled", "tss", "tss_within", "compressed", "patterns", "rows", "fe_info", "intercept", "weighted",
    "sumlogw", "x", "cols", "offset_call", "env", "terms", "xlevels", "spec", "na.action", "model", "call")))
expect_inherits(f$na.action, "exclude")
f <- design_lm(yna ~ x | id, con, table = "d2", weights = w, offset = z, na.action = na.exclude)
expect_identical(fit_names(f), sort(c("coefficients", "vcov", "vcov_type", "sigma", "rss", "df.residual",
    "nobs", "rank", "fixef", "fe_levels", "fe_vars", "fe_sizes", "fe_dof", "fe_dof_exact", "nsingletons",
    "cov.unscaled", "tss", "tss_within", "compressed", "patterns", "rows", "fe_info", "intercept", "weighted",
    "sumlogw", "x", "cols", "offset_call", "env", "terms", "xlevels", "spec", "na.action", "method", "call")))
expect_inherits(f$na.action, "dm_remote_exclude")
f <- design_glm(k ~ x | id, poisson, d2, weights = w, offset = log(z), na.action = na.exclude)
expect_inherits(f$na.action, "exclude")
expect_identical(f$rows, 39L)
f <- design_glm(k ~ x | id, poisson, con, table = "d2", weights = w, offset = log(z), na.action = na.exclude)
expect_inherits(f$na.action, "dm_remote_exclude")

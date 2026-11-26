# The data of the fixed-effects benchmark (see ../fixed_effects.R), shared by
# simulate.R (which saves each dataset once, as Parquet), cell.R (timed fits)
# and reference.R (tight-tolerance coefficients).
#
# The data generating process of Kyle Butts' fixest_benchmarks, as used by
# pyfixest: a panel of individuals over 10 years, about 23 individuals per
# firm. "simple": each row's firm is drawn at random. "difficult": firms are
# assigned round robin, so individuals move between firms every year and the
# fixed effects are poorly separated (slow convergence for demeaning).
dgp <- function(n, type, k, seed = 42) {
    set.seed(seed)
    nb_year <- 10
    nb_indiv <- round(n / nb_year)
    nb_firm <- round(nb_indiv / 23)
    n <- nb_indiv * nb_year
    indiv_id <- rep(seq_len(nb_indiv), each = nb_year)
    year <- rep(seq_len(nb_year), nb_indiv)
    firm_id <- if (type == "simple") sample.int(nb_firm, n, TRUE) else rep_len(seq_len(nb_firm), n)
    # x1..xk ~ N(0, 1), column by column (same draws as an n by k matrix)
    d <- data.frame(y = numeric(n))
    y <- numeric(n)
    for (j in seq_len(k)) {
        x <- rnorm(n)
        d[[paste0("x", j)]] <- x
        y <- y + x / j
    }
    firm_fe <- rnorm(nb_firm)
    unit_fe <- rnorm(nb_indiv)
    year_fe <- rnorm(nb_year)
    d$y <- y + firm_fe[firm_id] + unit_fe[indiv_id] + year_fe[year] + rnorm(n)
    d$indiv_id <- indiv_id
    d$year <- year
    d$firm_id <- firm_id
    d
}

# A Parquet file of simulate.R holds the data of one design and size for
# several k: x1..x_max(k), firm_id and the other ids are the same draws
# whatever k (they come first in the random stream), and the response for k
# regressors is column y_<k>. This query gives exactly dgp(n, type, k): y,
# x1..xk, indiv_id, year, firm_id, in that order.
fe_sql <- function(path, k) {
    sprintf("SELECT y_%d AS y, %s, indiv_id, year, firm_id FROM read_parquet('%s')",
            as.integer(k), paste0("x", seq_len(k), collapse = ", "), gsub("'", "''", path))
}

# The data for k regressors as a data frame (in-memory fits)
fe_read <- function(path, k) {
    con <- DBI::dbConnect(duckdb::duckdb())
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
    DBI::dbExecute(con, "SET threads = 1")
    DBI::dbGetQuery(con, fe_sql(path, k))
}

fe_formula <- function(k, nfe) {
    xs <- paste0("x", seq_len(k))
    fe <- c("indiv_id", "year", "firm_id")[seq_len(nfe)]
    f <- as.formula(paste("y ~", paste(xs, collapse = " + "), "|", paste(fe, collapse = " + ")))
    environment(f) <- globalenv()
    f
}

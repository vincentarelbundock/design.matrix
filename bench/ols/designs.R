## Benchmark designs: three families, each a sweep of widths with a fixed
## structure, plus one rank-deficient design. n rows; data built for one
## design at a time.
widths <- list(numeric = c(5, 20, 50, 100), mixed = c(25, 50, 100),
               factor = c(5, 20, 50, 100), interaction = c(6, 20, 50, 100), aliased = 49)

design_spec <- function(family, width) {
    switch(family,
        ## y ~ x1 + ... + xk: k = width - 1 numeric columns
        numeric = list(k = width - 1L, levels = NA),
        ## y ~ x1 + ... + xk + f, f with 20 levels: k = width - 20
        mixed = list(k = width - 20L, levels = 20L),
        ## y ~ x + f: 1 + 1 + (levels - 1) columns
        factor = list(k = 1L, levels = width - 1L),
        ## y ~ x * f: 2 * levels columns
        interaction = list(k = 1L, levels = max(2L, round(width / 2))),
        ## y ~ x + f1 * f3 (12 x 4 levels) with empty cells: 49 columns,
        ## 6 aliased
        aliased = list(k = 1L, levels = NA),
        stop("unknown family"))
}

make_data <- function(family, width, n) {
    s <- design_spec(family, width)
    set.seed(1)
    d <- data.frame(y = rnorm(n))
    if (family %in% c("numeric", "mixed")) {
        for (j in seq_len(s$k)) d[[paste0("x", j)]] <- rnorm(n)
        rhs <- paste0("x", seq_len(s$k))
        if (family == "mixed") {
            d$f <- factor(sample(sprintf("L%03d", seq_len(s$levels)), n, TRUE))
            rhs <- c(rhs, "f")
        }
        fo <- reformulate(rhs, "y")
    } else if (family %in% c("factor", "interaction")) {
        d$x <- rnorm(n)
        d$f <- factor(sample(sprintf("L%03d", seq_len(s$levels)), n, TRUE))
        fo <- if (family == "factor") y ~ x + f else y ~ x * f
    } else {
        d$x <- rnorm(n)
        d$f1 <- factor(sample(sprintf("a%02d", 1:12), n, TRUE))
        d$f3 <- factor(sample(1:4, n, TRUE))
        d$f3[d$f1 %in% c("a01", "a02")] <- "1"
        fo <- y ~ x + f1 * f3
    }
    environment(fo) <- globalenv()
    list(data = d, formula = fo)
}

# The C kernels of the fixed-effect fits (src/fe.c) against the R code
# they replace.
library(design.matrix)
ns <- asNamespace("design.matrix")

# ---- FE codes (dm_fe_code()): direct indexing or hashing, rows per level
# the codes as sort(unique()) and match() give them (factors: level order)
code_ref <- function(v) {
    if (is.factor(v)) {
        used <- tabulate(as.integer(v), nlevels(v)) > 0L
        return(list(codes = cumsum(used)[as.integer(v)], G = sum(used), levels = levels(v)[used]))
    }
    if (is.logical(v)) v <- as.integer(v)
    u <- sort(unique(v), method = "radix")
    list(codes = match(v, u), G = length(u), levels = u)
}
set.seed(1)
n <- 5000
cases <- list(
    int = sample(1:50, n, TRUE),
    int_wide = sample(c(1L, 1000000000L, -5L, 7L), n, TRUE),
    dbl = sample(c(0.5, -0, 0, 3, 1e300, -2.25), n, TRUE),
    dbl_int = as.double(sample(-3:30, n, TRUE)),
    dbl_huge = sample(c(2^60, -2^60, 1), n, TRUE),
    chr = sample(c("b", "a", "Z", "B", "aa", ""), n, TRUE),
    chr_runs = rep(sprintf("id%03d", sample(200)), each = 25),
    fac = factor(sample(letters, n, TRUE), levels = rev(letters)),
    fac_unused = factor(sample(c("x", "y"), n, TRUE), levels = c("y", "q", "x")),
    lgl = sample(c(TRUE, FALSE), n, TRUE),
    const = rep(7, n),
    many = sample(1e5, n, TRUE) + 0.5)
for (nm in names(cases)) {
    v <- cases[[nm]]
    a <- ns$dm_fe_code(v)
    b <- code_ref(v)
    expect_identical(a$codes, b$codes, info = nm)
    expect_identical(a$G, b$G, info = nm)
    expect_identical(a$levels, b$levels, info = nm)
    expect_identical(a$n, tabulate(b$codes, b$G), info = nm)
}
# a string in two encodings: one level, as match() has it
lat <- rawToChar(as.raw(c(0x63, 0x61, 0x66, 0xe9)))
Encoding(lat) <- "latin1"
enc <- c(enc2utf8(lat), lat, "tea")[sample(3, 300, TRUE)]
expect_identical(ns$dm_fe_code(enc)[c("codes", "G", "levels")], code_ref(enc))
# interactions: the combined codes and their levels
d <- data.frame(a = sample(c("u", "v", "w"), n, TRUE), b = sample(c(10.5, 2, 3), n, TRUE))
ct <- ns$dm_fe_code_term(d, c("a", "b"))
key <- paste(d$a, d$b)
lv <- unique(d[order(d$a, d$b, method = "radix"), ])
expect_identical(ct$levels$a[ct$codes], d$a)
expect_identical(ct$levels$b[ct$codes], d$b)
expect_identical(nrow(unique(ct$levels)), ct$G)
expect_identical(ct$n, tabulate(ct$codes, ct$G))
expect_equal(unname(as.list(ct$levels)), unname(as.list(lv)), check.attributes = FALSE)

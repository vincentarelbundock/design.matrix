## Column-structured algebra on a plain design (see src/structured.c).
##
## Term t of the design has X_t[i, colmap_t[u, k]] = N_t(i)[u] C_t[F_t(i), k],
## with F_t(i) the combined level of its factors, N_t(i) the products of its
## numerics and C_t the Kronecker product of its contrasts. Products with the
## design are sums of N's within levels (computed in C), followed by small
## dense products with the C_t (here). Used only when every input is finite
## and no product can overflow; the results then agree with the dense
## computation up to rounding. Option `design.matrix.structured = FALSE`
## turns these kernels off (the row kernels are used instead).

## Is the structured path available for this plain design? Cached.
dm_struct_ok <- function(d) {
    if (!isTRUE(getOption("design.matrix.structured", TRUE))) return(FALSE)
    if (dm_is_comb(d) || !kernel_ok(d)) return(FALSE)
    cache <- d$cache
    if (!is.null(cache) && !is.null(cache$struct_ok)) return(cache$struct_ok)
    ok <- .Call(C_dm_struct_ok, d$setup, d$vars, 2^23)
    if (!is.null(cache)) cache$struct_ok <- ok
    ok
}

## Term descriptions: off, colmap (0-based), fvar (the variable of each
## factor), fcm (their contrast matrices; C_t is their Kronecker product, the
## first factor varying fastest), cols (1-based global columns, ns x kc) and
## L (the combined levels). Cached.
dm_struct_info <- function(d) {
    cache <- d$cache
    if (!is.null(cache) && !is.null(cache$struct_info)) return(cache$struct_info)
    info <- lapply(.Call(C_dm_struct_info, d$setup, d$vars), function(t) {
        names(t) <- c("off", "colmap", "fvar", "fcm")
        t$cols <- t$off + t$colmap + 1L
        t$L <- prod(vapply(t$fcm, nrow, 1L))
        t
    })
    if (!is.null(cache)) cache$struct_info <- info
    info
}

## Products with Kronecker products of contrasts, one factor (mode) at a
## time, instead of with the dense L x kc matrix C_t: X holds an array whose
## first dims are the factors' levels; mats[[f]] is list(A) (contract the
## mode with A: sum_a X[a, ...] A[a, i]) or list(A, B) (a factor two terms
## share: sum_a X[a, ...] A[a, i] B[a, j]). Returns the result with the
## remaining dims first, then the new dims in order. The cost is the size
## of the moments times the non-zero entries per level, not L x kc.
dm_kron_modes <- function(X, mats) {
    for (m in mats) X <- .Call(C_dm_kron_mode, X, m[[1L]], if (length(m) > 1L) m[[2L]])
    X
}

## t(X) %*% diag(w) %*% X, w NULL or finite. The moments are by cell of
## the union of each pair of terms' factors and are contracted with the
## factors' contrasts one at a time, in C (gram_block() in src/structured.c).
dm_struct_gram <- function(d, w = NULL) {
    .Call(C_dm_struct_gram, d$setup, d$vars, if (is.null(w)) NULL else as.double(w), design_threads())
}

## t(X) %*% Y, Y an n x K finite matrix.
dm_struct_xty <- function(d, Y) {
    info <- dm_struct_info(d)
    mom <- .Call(C_dm_struct_xty, d$setup, d$vars, Y, design_threads())
    K <- ncol(Y)
    out <- matrix(0, ncol(d), K)
    for (t in seq_along(info)) {
        Tt <- info[[t]]
        out[Tt$cols, ] <- dm_struct_contract(mom[[t]], Tt, K)
    }
    out
}

## Moments of one term, [L, ns, K], contracted with its contrasts: the
## (ns kc) x K matrix whose rows follow the term's columns (Tt$cols).
dm_struct_contract <- function(M, Tt, K) {
    nf <- length(Tt$fcm)
    A <- dm_kron_modes(M, lapply(Tt$fcm, list))
    if (nf == 0L) return(A)
    dim(A) <- c(nrow(Tt$cols), K, vapply(Tt$fcm, ncol, 1L))
    A <- aperm(A, c(1L, 2L + seq_len(nf), 2L))
    dim(A) <- c(length(Tt$cols), K)
    A
}

## X %*% B, B a p x K finite matrix.
dm_struct_mult <- function(d, B) {
    info <- dm_struct_info(d)
    K <- ncol(B)
    P <- lapply(info, function(Tt) {
        ## P_t = C_t B_t as [L, ns, K]: B_t is [ns, k_1, .., k_m, K]
        ns <- nrow(Tt$cols)
        nf <- length(Tt$fcm)
        Bt <- B[Tt$cols, , drop = FALSE]
        if (nf == 0L) return(Bt)
        dim(Bt) <- c(ns, vapply(Tt$fcm, ncol, 1L), K)
        Bt <- aperm(Bt, c(1L + seq_len(nf), 1L, nf + 2L))
        A <- dm_kron_modes(Bt, lapply(Tt$fcm, function(cm) list(t(cm))))
        dim(A) <- c(ns * K, Tt$L)
        t(A)
    })
    .Call(C_dm_struct_mult, d$setup, d$vars, P, K, design_threads())
}

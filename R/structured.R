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

## Term descriptions: off, C, colmap (0-based). Cached.
dm_struct_info <- function(d) {
    cache <- d$cache
    if (!is.null(cache) && !is.null(cache$struct_info)) return(cache$struct_info)
    info <- lapply(.Call(C_dm_struct_info, d$setup, d$vars), function(t) {
        names(t) <- c("off", "C", "colmap")
        t$cols <- t$off + t$colmap + 1L # 1-based global columns, ns x kc
        t
    })
    if (!is.null(cache)) cache$struct_info <- info
    info
}

## t(X) %*% diag(w) %*% X, w NULL or finite.
dm_struct_gram <- function(d, w = NULL) {
    info <- dm_struct_info(d)
    mom <- .Call(C_dm_struct_gram, d$setup, d$vars, if (is.null(w)) NULL else as.double(w))
    p <- ncol(d)
    G <- matrix(0, p, p)
    q <- 0L
    for (t in seq_along(info)) {
        Tt <- info[[t]]
        Lt <- nrow(Tt$C)
        nt <- nrow(Tt$cols)
        for (s in seq_len(t)) {
            q <- q + 1L
            Ts <- info[[s]]
            Ls <- nrow(Ts$C)
            ns <- nrow(Ts$cols)
            M <- mom[[q]]
            if (s == t) {
                dim(M) <- c(Ls, ns, ns)
                for (u in seq_len(ns)) for (v in seq_len(ns)) {
                    blk <- crossprod(Ts$C, M[, u, v] * Ts$C)
                    G[Ts$cols[u, ], Ts$cols[v, ]] <- G[Ts$cols[u, ], Ts$cols[v, ]] + blk
                }
            } else {
                dim(M) <- c(Ls, Lt, ns, nt)
                for (u in seq_len(ns)) for (v in seq_len(nt)) {
                    Muv <- M[, , u, v]
                    dim(Muv) <- c(Ls, Lt)
                    blk <- crossprod(Ts$C, Muv %*% Tt$C)
                    G[Ts$cols[u, ], Tt$cols[v, ]] <- G[Ts$cols[u, ], Tt$cols[v, ]] + blk
                    G[Tt$cols[v, ], Ts$cols[u, ]] <- G[Tt$cols[v, ], Ts$cols[u, ]] + t(blk)
                }
            }
        }
    }
    G
}

## t(X) %*% Y, Y an n x K finite matrix.
dm_struct_xty <- function(d, Y) {
    info <- dm_struct_info(d)
    mom <- .Call(C_dm_struct_xty, d$setup, d$vars, Y)
    K <- ncol(Y)
    out <- matrix(0, ncol(d), K)
    for (t in seq_along(info)) {
        Tt <- info[[t]]
        L <- nrow(Tt$C)
        ns <- nrow(Tt$cols)
        M <- mom[[t]]
        dim(M) <- c(L, ns, K)
        for (u in seq_len(ns)) {
            Mu <- M[, u, ]
            dim(Mu) <- c(L, K)
            out[Tt$cols[u, ], ] <- crossprod(Tt$C, Mu)
        }
    }
    out
}

## X %*% B, B a p x K finite matrix.
dm_struct_mult <- function(d, B) {
    info <- dm_struct_info(d)
    K <- ncol(B)
    P <- lapply(info, function(Tt) {
        L <- nrow(Tt$C)
        ns <- nrow(Tt$cols)
        out <- array(0, c(L, ns, K))
        for (u in seq_len(ns)) out[, u, ] <- Tt$C %*% B[Tt$cols[u, ], , drop = FALSE]
        out
    })
    .Call(C_dm_struct_mult, d$setup, d$vars, P, K)
}

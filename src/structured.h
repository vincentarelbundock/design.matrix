/*
 *  design.matrix: the term structure of a plain design (structured.c),
 *  shared with the fixed-effect kernels (fe.c).
 *
 *  Copyright (C) 2026 Vincent Arel-Bundock.  GPL (>= 2).
 */
#ifndef DM_STRUCTURED_H
#define DM_STRUCTURED_H
#include <Rinternals.h>

#define SMAXC 32
#define BLOCK 2048

typedef struct {
    int L, ns, nf, nn, kc, count, off;
    const int *code[SMAXC];
    int fvar[SMAXC];     /* index of each factor's variable (shared across terms) */
    SEXP fcm[SMAXC];     /* each factor's contrast matrix */
    int adj[SMAXC], lstride[SMAXC], nrc[SMAXC], fnc[SMAXC], fstride[SMAXC];
    const double *cm[SMAXC];
    const double *num[SMAXC];
    int nnc[SMAXC], nstride[SMAXC];
} sterm;

typedef struct {
    R_xlen_t n;
    int p, nt;
    sterm *t;
} sctx;

void build_sctx(sctx *S, SEXP setup, SEXP vars);
void term_block(const sterm *T, R_xlen_t n, R_xlen_t i0, int nb, int *F, double *N);

#endif

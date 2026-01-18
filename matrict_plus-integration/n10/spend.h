#ifndef _SPENT_H
#define _SPENT_H

#include "param.h"
#include "poly_q.h"
#include <stdint.h>

typedef struct {
  POLY_Q pk[N_BAR];
  POLY_Q cn[N + 1];
} ACT;

typedef struct {
  POLY_R sk[N_BAR + K];
  unsigned char sn[CRYPTO_BYTES];
  POLY_R ck[N + K + 1];
  uint64_t amt_in;
} ASK;

typedef struct {
  unsigned char sn[CRYPTO_BYTES];
} SN_OUT;

typedef struct {
  POLY_Q pk[N_BAR];
} PK_OUT;

typedef struct {
  POLY_Q cn[N + 1];
} CN_OUT;

typedef struct {
  POLY_R ck[N + K + 1];
} CK_OUT;

typedef struct {
  POLY_Q c;
  POLY_Q v[N];
  POLY_Q v1;
  POLY_Q v2;
  POLY_QHAT b[N_HAT];
  unsigned char alpha[CRYPTO_BYTES];
  POLY_R x;
  POLY_R z_out[S][N + K + 1];
  POLY_R z0[N + K + 1];
  POLY_R z1[N + K + 3];
  POLY_R z2[N + K + 3];
  POLY_R f1[N_SPENT + 1];
  POLY_R zb[N_HAT + K_HAT];
  POLY_R z[N_BAR + K];
} SPEND_OUT;

void spend(SPEND_OUT *out, CN_OUT *cn_out, CK_OUT *ck_out, SN_OUT *sn_out,
           const ACT a_in[][M], const uint64_t l, const ASK *ask,
           const PK_OUT *pk_out, const uint64_t *amt_out,
           const POLY_Q mat1[][N], const POLY_Q mat2[][3],
           const POLY_QHAT mat3[][N_HAT], unsigned char *pp);

#endif

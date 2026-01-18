#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "fastrandombytes.h"
#include "keygen.h"
#include "param.h"
#include "poly_mult.h"
#include "poly_q.h"
#include "poly_red.h"
#include "sammat.h"
#include "sampleb.h"

void keygen(POLY_Q *pk, POLY_R *sk, unsigned char *sn, const POLY_Q mat1[][N]) {
  POLY_Q r[N_BAR + K];
  POLY_Q g_gamma[N_BAR];

  uint64_t i, j;

  /* gamma <-- {0, 1}^256 */
  fastrandombytes_prv(sn, CRYPTO_BYTES);
  /* G(gamma) */
  sample_g(g_gamma, sn);

  /* r <-- S_1^{N_BAR + K} */
  sample_1(sk, N_BAR + K);

  /* ntt */
  for (i = 0; i < N_BAR + K; i++) {
    for (j = 0; j < D; j++) {
      r[i].poly[j] = con_add(sk[i].poly[j], Q);
    }

    ntt_q(r + i);
  }

  /* c = G * r + G(gamma) */
  for (i = 0; i < N_BAR; i++) {
    memcpy(pk + i, r + i, sizeof(POLY_Q));
  }
  for (i = 0; i < K; i++) {
    for (j = 0; j < N_BAR; j++) {
      mult_plus_rq(pk + j, mat1[i] + j + (N - N_BAR), r + N_BAR + i);
    }
  }
  for (i = 0; i < N_BAR; i++) {
    for (j = 0; j < D; j++) {
      pk[i].poly[j] = con_sub(pk[i].poly[j] + g_gamma[i].poly[j], Q);
    }
  }
}

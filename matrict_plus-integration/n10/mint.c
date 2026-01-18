#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "mint.h"
#include "param.h"
#include "poly_mult.h"
#include "poly_q.h"
#include "poly_red.h"
#include "sampleb.h"

void mint(POLY_Q *cn, POLY_R *ck, uint64_t *a_hat, const uint64_t amt,
          const POLY_Q mat1[][N], const POLY_Q mat2[][3]) {
  uint64_t i, j;
  POLY_Q r[N + K + 1];

  /* a_hat = [a[0],..., a[r - 1]] */
  for (i = 0; i < R; i++) {
    a_hat[i] = (amt >> i) & 0x1;
  }

  /* r <-- S_1^{N + K + 1} */
  sample_1(ck, N + K + 1);

  /* ntt */
  for (i = 0; i < N + K + 1; i++) {
    for (j = 0; j < D; j++) {
      r[i].poly[j] = con_add(ck[i].poly[j], Q);
    }

    ntt_q(r + i);
  }

  /* t = B * r */
  for (i = 0; i < N + 1; i++) {
    memcpy(cn + i, r + i, sizeof(POLY_Q));
  }
  for (i = 0; i < K + 1; i++) {
    for (j = 0; j < N; j++) {
      mult_plus_rq(cn + j, mat1[i] + j, r + N + i);
    }
  }

  /* t_1 = <b_1, r> + a_hat */
  for (i = 0; i < K - 2; i++) {
    mult_plus_rq(cn + N, mat2[i], r + N + 3 + i);
  }
  for (i = 0; i < R; i++) {
    cn[N].poly[i * SUBRING_SIZE] = con_sub(
        cn[N].poly[i * SUBRING_SIZE] + ((-a_hat[i]) & MONTGOMERY_32_Q), Q);
  }
}

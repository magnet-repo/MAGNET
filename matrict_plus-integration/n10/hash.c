#include "hash.h"
#include "littleendian.h"
#include "param.h"
#include "poly_q.h"
#include "spend.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <libXKCP.a.headers/SimpleFIPS202.h>

#define HASH_ALPHA_LEN                                                         \
  (CRYPTO_BYTES + N_SPENT * M * (N_BAR * D * Q_BYTE + (N + 1) * D * Q_BYTE) +  \
   M * CRYPTO_BYTES + S * N_BAR * D * Q_BYTE + S * (N + 1) * D * Q_BYTE +      \
   S * N * D * Q_BYTE + N * D * Q_BYTE + N * D * Q_BYTE + D * Q_BYTE +         \
   D * Q_BYTE + D * Q_BYTE + N * D * Q_BYTE + N * D * Q_BYTE)

#define HASH_X_INPUT_LEN                                                       \
  (CRYPTO_BYTES + D * Q_BYTE + D * Q_BYTE + N_HAT * D * QHAT_BYTE +            \
   N_HAT * D * QHAT_BYTE + N_BAR * D * Q_BYTE)
#define HASH_X_OUTPUT_LEN 137LL
#define HASH_X_COEFF_LEN 6LL
#define HASH_X_GROUP 4LL
#define HASH_X_W (W / HASH_X_GROUP)
#define HASH_X_MASK (D / HASH_X_GROUP - 1)

void hash_alpha(unsigned char *out, const unsigned char *pp,
                const ACT a_in[][M], const SN_OUT *sn_out, const PK_OUT *pk_out,
                const CN_OUT *cn_out, const POLY_Q w_out[][N], const POLY_Q *w0,
                const POLY_Q *v, const POLY_Q *v2, const POLY_Q *c,
                const POLY_Q *w, const POLY_Q *w1, const POLY_Q *w2) {
  static unsigned char hash_input[HASH_ALPHA_LEN];
  unsigned char *hash_input_head;

  uint64_t i, j, k;

  /* pp */
  memcpy(hash_input, pp, CRYPTO_BYTES);
  hash_input_head = hash_input + CRYPTO_BYTES;

  /* A_in */
  for (i = 0; i < N_SPENT; i++) {
    for (j = 0; j < N_BAR; j++) {
      for (k = 0; k < D; k++) {
        STORE_Q(hash_input_head, a_in[i][0].pk[j].poly[k]);
        hash_input_head += Q_BYTE;
      }
    }
    for (j = 0; j < N + 1; j++) {
      for (k = 0; k < D; k++) {
        STORE_Q(hash_input_head, a_in[i][0].cn[j].poly[k]);
        hash_input_head += Q_BYTE;
      }
    }

    for (j = 0; j < N_BAR; j++) {
      for (k = 0; k < D; k++) {
        STORE_Q(hash_input_head, a_in[i][1].pk[j].poly[k]);
        hash_input_head += Q_BYTE;
      }
    }
    for (j = 0; j < N + 1; j++) {
      for (k = 0; k < D; k++) {
        STORE_Q(hash_input_head, a_in[i][1].cn[j].poly[k]);
        hash_input_head += Q_BYTE;
      }
    }
  }

  /* SN */
  memcpy(hash_input_head, sn_out[0].sn, CRYPTO_BYTES);
  hash_input_head += CRYPTO_BYTES;
  memcpy(hash_input_head, sn_out[1].sn, CRYPTO_BYTES);
  hash_input_head += CRYPTO_BYTES;

  /* PK_out */
  for (i = 0; i < N_BAR; i++) {
    for (j = 0; j < D; j++) {
      STORE_Q(hash_input_head, pk_out[0].pk[i].poly[j]);
      hash_input_head += Q_BYTE;
    }
  }
  for (i = 0; i < N_BAR; i++) {
    for (j = 0; j < D; j++) {
      STORE_Q(hash_input_head, pk_out[1].pk[i].poly[j]);
      hash_input_head += Q_BYTE;
    }
  }

  /* CN_out */
  for (i = 0; i < N + 1; i++) {
    for (j = 0; j < D; j++) {
      STORE_Q(hash_input_head, cn_out[0].cn[i].poly[j]);
      hash_input_head += Q_BYTE;
    }
  }
  for (i = 0; i < N + 1; i++) {
    for (j = 0; j < D; j++) {
      STORE_Q(hash_input_head, cn_out[1].cn[i].poly[j]);
      hash_input_head += Q_BYTE;
    }
  }

  /* w_out */
  for (i = 0; i < N; i++) {
    for (j = 0; j < D; j++) {
      STORE_Q(hash_input_head, w_out[0][i].poly[j]);
      hash_input_head += Q_BYTE;
    }
  }
  for (i = 0; i < N; i++) {
    for (j = 0; j < D; j++) {
      STORE_Q(hash_input_head, w_out[1][i].poly[j]);
      hash_input_head += Q_BYTE;
    }
  }

  /* w_0 */
  for (i = 0; i < N; i++) {
    for (j = 0; j < D; j++) {
      STORE_Q(hash_input_head, w0[i].poly[j]);
      hash_input_head += Q_BYTE;
    }
  }

  /* v */
  for (i = 0; i < N; i++) {
    for (j = 0; j < D; j++) {
      STORE_Q(hash_input_head, v[i].poly[j]);
      hash_input_head += Q_BYTE;
    }
  }

  /* v_2 */
  for (i = 0; i < D; i++) {
    STORE_Q(hash_input_head, v2->poly[i]);
    hash_input_head += Q_BYTE;
  }

  /* C */
  for (i = 0; i < D; i++) {
    STORE_Q(hash_input_head, c->poly[i]);
    hash_input_head += Q_BYTE;
  }

  /* w */
  for (i = 0; i < D; i++) {
    STORE_Q(hash_input_head, w->poly[i]);
    hash_input_head += Q_BYTE;
  }

  /* w_1 */
  for (i = 0; i < N; i++) {
    for (j = 0; j < D; j++) {
      STORE_Q(hash_input_head, w1[i].poly[j]);
      hash_input_head += Q_BYTE;
    }
  }

  /* w_2 */
  for (i = 0; i < N; i++) {
    for (j = 0; j < D; j++) {
      STORE_Q(hash_input_head, w2[i].poly[j]);
      hash_input_head += Q_BYTE;
    }
  }

  SHAKE256(out, CRYPTO_BYTES, hash_input, HASH_ALPHA_LEN);
}

void hash_x(POLY_R *out, const unsigned char *alpha, const POLY_Q *v0,
            const POLY_Q *v1, const POLY_QHAT *a, const POLY_QHAT *b,
            const POLY_Q *e) {
  static unsigned char hash_input[HASH_X_INPUT_LEN];
  unsigned char hash_output[HASH_X_OUTPUT_LEN];

  unsigned char *hash_input_head;

  uint64_t i, j, boo, x_pos, k;
  uint64_t coeff, count = 0;

  unsigned char *hash_output_head = hash_output + HASH_X_COEFF_LEN;
  uint64_t nonzero_pos[HASH_X_GROUP][HASH_X_W];

  /* alpha */
  memcpy(hash_input, alpha, CRYPTO_BYTES);
  hash_input_head = hash_input + CRYPTO_BYTES;

  /* v_0 */
  for (i = 0; i < D; i++) {
    STORE_Q(hash_input_head, v0->poly[i]);
    hash_input_head += Q_BYTE;
  }

  /* v_1 */
  for (i = 0; i < D; i++) {
    STORE_Q(hash_input_head, v1->poly[i]);
    hash_input_head += Q_BYTE;
  }

  /* A */
  for (i = 0; i < N_HAT; i++) {
    for (j = 0; j < D; j++) {
      STORE_QHAT(hash_input_head, a[i].poly[j]);
      hash_input_head += QHAT_BYTE;
    }
  }

  /* B */
  for (i = 0; i < N_HAT; i++) {
    for (j = 0; j < D; j++) {
      STORE_QHAT(hash_input_head, b[i].poly[j]);
      hash_input_head += QHAT_BYTE;
    }
  }

  /* E */
  for (i = 0; i < N_BAR; i++) {
    for (j = 0; j < D; j++) {
      STORE_Q(hash_input_head, e[i].poly[j]);
      hash_input_head += Q_BYTE;
    }
  }

  SHAKE256(hash_output, HASH_X_OUTPUT_LEN, hash_input, HASH_X_INPUT_LEN);

  /* sample and fill the nonzero positions from the hash output
   * since it may generates duplicate positions, we need to do a rejection
   * sampling here */
  memset(out, 0, sizeof(POLY_R));
  coeff = load_48(hash_output);
  for (i = 0; i < HASH_X_GROUP; i++) {
    nonzero_pos[i][0] = (*(hash_output_head++)) & HASH_X_MASK;
    out->poly[nonzero_pos[i][0] * HASH_X_GROUP + i] =
        1 - 2 * ((coeff >> (count++)) & 0x1);
  }
  for (k = 0; k < HASH_X_GROUP; k++) {
    for (i = 1; i < HASH_X_W; i++) {
      do {
        boo = 0;
        x_pos = (*(hash_output_head++)) & HASH_X_MASK;

        for (j = 0; j < i; j++) {
          if (x_pos == nonzero_pos[k][j]) {
            boo = 1;
            break;
          }
        }
      } while (boo);

      nonzero_pos[k][i] = x_pos;
      out->poly[x_pos * HASH_X_GROUP + k] =
          1 - 2 * ((coeff >> (count++)) & 0x1);
    }
  }
}

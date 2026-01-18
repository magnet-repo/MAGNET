#include <stdint.h>

#include "fastrandombytes.h"
#include "littleendian.h"
#include "param.h"
#include "poly_q.h"
#include "sammat.h"

#define SAMMAT_Q_BYTES 4LL
#define SAMMAT_Q_BOUND 4194256025LL

#define SAMMAT_MAT1_LEN 7648LL
#define SAMMAT_MAT1_KEYGEN_LEN 4440LL
#define SAMMAT_MAT2_LEN 3366LL

#define SAMMAT_QHAT_BYTES 5LL
#define SAMMAT_QHAT_BOUND 1099243232955LL
#define SAMMAT_MAT3_LEN 37232LL

/* --------------------------
 * |I_n       | Mat_1       |
 * |------------------------|
 * |0_{3 * n} | I_3 | Mat_2 |
 * --------------------------
 *
 * where Mat_1 is n * (k + 3), and Mat_2 is 3 * k */

/* sample Mat_1 matrix */
static inline void sample_mat_g(POLY_Q out[][N], const uint64_t m,
                                const uint64_t rej_len) {
  static unsigned char r[SAMMAT_MAT1_LEN * SAMMAT_Q_BYTES];
  uint64_t i, j, k, x;

  unsigned char *r_head = r;

  fastrandombytes_setiv_mat1();

  fastrandombytes_pub(r, SAMMAT_Q_BYTES * rej_len);

  for (i = 0; i < m; i++) {
    for (j = 0; j < N; j++) {
      for (k = 0; k < D; k++) {
        do {
          x = load_32(r_head);
          r_head += SAMMAT_Q_BYTES;
        } while (x >= SAMMAT_Q_BOUND);

        out[i][j].poly[k] = x % Q;
      }
    }
  }
}

void sample_mat1(POLY_Q out[][N]) { sample_mat_g(out, K + 3, SAMMAT_MAT1_LEN); }

void sample_mat1_keygen(POLY_Q out[][N]) {
  sample_mat_g(out, K, SAMMAT_MAT1_KEYGEN_LEN);
}

/* sample Mat_2 matrix */
void sample_mat2(POLY_Q out[][3]) {
  static unsigned char r[SAMMAT_MAT2_LEN * SAMMAT_Q_BYTES];
  uint64_t i, j, k, x;

  unsigned char *r_head = r;

  fastrandombytes_setiv_mat2();

  fastrandombytes_pub(r, SAMMAT_MAT2_LEN * SAMMAT_Q_BYTES);

  for (i = 0; i < K; i++) {
    for (j = 0; j < 2; j++) {
      for (k = 0; k < D; k++) {
        do {
          x = load_32(r_head);
          r_head += SAMMAT_Q_BYTES;
        } while (x >= SAMMAT_Q_BOUND);

        out[i][j].poly[k] = x % Q;
      }
    }
  }

  /* sample a_3 to have 0 in first CRT coefficients of both halves */
  for (i = 0; i < K; i++) {
    for (j = 0; j < SUBRING_SIZE; j++) {
      out[i][2].poly[j] = 0;
      out[i][2].poly[j + (D >> 1)] = 0;
    }

    for (k = SUBRING_SIZE; k < (D >> 1); k++) {
      do {
        x = load_32(r_head);
        r_head += SAMMAT_Q_BYTES;
      } while (x >= SAMMAT_Q_BOUND);

      out[i][2].poly[k] = x % Q;

      do {
        x = load_32(r_head);
        r_head += SAMMAT_Q_BYTES;
      } while (x >= SAMMAT_Q_BOUND);

      out[i][2].poly[(D >> 1) + k] = x % Q;
    }
  }
}

/* --------------------------
 * |I_\hat{n} | Mat_3       |
 * |------------------------|
 *
 * where Mat_3 is \hat{n} * (\hat{k} + 2 * N_SPENT + 4) */
/* sample Mat_3 matrix */
void sample_mat3(POLY_QHAT out[][N_HAT]) {
  static unsigned char r[SAMMAT_MAT3_LEN * SAMMAT_QHAT_BYTES];
  uint64_t i, j, k, x;

  unsigned char *r_head = r;

  fastrandombytes_setiv_mat3();

  fastrandombytes_pub(r, SAMMAT_QHAT_BYTES * SAMMAT_MAT3_LEN);

  for (i = 0; i < K_HAT + 2 * N_SPENT + 4; i++) {
    for (j = 0; j < N_HAT; j++) {
      for (k = 0; k < D; k++) {
        do {
          x = load_40(r_head);
          r_head += SAMMAT_QHAT_BYTES;
        } while (x >= SAMMAT_QHAT_BOUND);

        out[i][j].poly[k] = x % QHAT;
      }
    }
  }
}

#include "sampleb.h"
#include "comp.h"
#include "fastrandombytes.h"
#include "littleendian.h"
#include "param.h"
#include "poly_q.h"
#include "poly_red.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define SAMPLE_1_MAX_M (N + K + 3)

#define SAMPLE_G_LEN 1197LL
#define SAMPLE_G_BYTE 4LL
#define SAMPLE_G_BOUND 4194256025LL
#define BARRETT_FACTOR_G 25LL
#define BARRETT_SHIFT_G 32LL

#define SAMPLE_ALPHA_LEN 643LL
#define SAMPLE_ALPHA_BYTE 4LL
#define SAMPLE_ALPHA_BOUND 4194256025LL
#define BARRETT_FACTOR_ALPHA 25LL
#define BARRETT_SHIFT_ALPHA 32LL

#define SAMPLE_ALPHA_PRIME_LEN 127LL
#define SAMPLE_ALPHA_PRIME_COEFF_LEN 5LL

/* out <-- S_1^{m*d} */
void sample_1(POLY_R *out, const uint64_t m) {
  unsigned char r[SAMPLE_1_MAX_M * D / 2];

  uint64_t i, j, x1, x2, r_head = 0;

  fastrandombytes_prv(r, m * D / 2);

  for (i = 0; i < m; i++) {
    for (j = 0; j < D; j += 2) {
      x1 = (r[r_head] & 0x1) + ((r[r_head] >> 1) & 0x1) -
           ((r[r_head] >> 2) & 0x1) - ((r[r_head] >> 3) & 0x1);
      x2 = ((r[r_head] >> 4) & 0x1) + ((r[r_head] >> 5) & 0x1) -
           ((r[r_head] >> 6) & 0x1) - ((r[r_head] >> 7) & 0x1);

      x1 += ct_eq(x1, -2) - ct_eq(x1, 2);
      x2 += ct_eq(x2, -2) - ct_eq(x2, 2);

      out[i].poly[j] = x1;
      out[i].poly[j + 1] = x2;

      r_head++;
    }
  }
}

/* out <-- R_{q}^{\hat{n}} */
void sample_g(POLY_Q *out, const unsigned char *seed) {
  unsigned char r[SAMPLE_G_LEN * SAMPLE_G_BYTE];

  uint64_t i, j, x;
  unsigned char *r_head = r;

  fastrandombytes_setseed_tmp(seed);
  fastrandombytes_tmp(r, SAMPLE_G_LEN * SAMPLE_G_BYTE);

  for (i = 0; i < N_BAR; i++) {
    for (j = 0; j < D; j++) {
      do {
        x = load_32(r_head);
        r_head += SAMPLE_G_BYTE;
      } while (ct_ge(x, SAMPLE_G_BOUND));

      out[i].poly[j] = barrett(x, Q, BARRETT_FACTOR_G, BARRETT_SHIFT_G);
    }
  }
}

/* out <-- R_{q}^{S} */
void sample_alpha(POLY_Q *out) {
  unsigned char r[SAMPLE_ALPHA_LEN * SAMPLE_ALPHA_BYTE];

  uint64_t i, j, x;
  unsigned char *r_head = r;

  fastrandombytes_ch(r, SAMPLE_ALPHA_LEN * SAMPLE_ALPHA_BYTE);

  for (i = 0; i < S; i++) {
    for (j = 0; j < D; j++) {
      do {
        x = load_32(r_head);
        r_head += SAMPLE_ALPHA_BYTE;
      } while (ct_ge(x, SAMPLE_ALPHA_BOUND));

      out[i].poly[j] = barrett(x, Q, BARRETT_FACTOR_ALPHA, BARRETT_SHIFT_ALPHA);
    }
  }
}

/* out \in C' */
void sample_alpha_prime(POLY_R *out) {
  unsigned char hash_output[SAMPLE_ALPHA_PRIME_LEN];

  uint64_t i, j, boo, x_pos;
  uint64_t coeff;

  unsigned char *hash_output_head = hash_output + SAMPLE_ALPHA_PRIME_COEFF_LEN;
  uint64_t nonzero_pos[WA];

  fastrandombytes_ch(hash_output, SAMPLE_ALPHA_PRIME_LEN);

  /* sample and fill the nonzero positions from the hash output
   * since it may generates duplicate positions, we need to do a rejection
   * sampling here */
  memset(out, 0, sizeof(POLY_R));
  coeff = load_40(hash_output);
  nonzero_pos[0] = *(hash_output_head++);
  out->poly[nonzero_pos[0]] = 1 - 2 * (coeff & 0x1);
  for (i = 1; i < WA; i++) {
    do {
      boo = 0;
      x_pos = *(hash_output_head++);

      for (j = 0; j < i; j++) {
        if (x_pos == nonzero_pos[j]) {
          boo = 1;
          break;
        }
      }
    } while (boo);

    nonzero_pos[i] = x_pos;
    out->poly[x_pos] = 1 - 2 * ((coeff >> i) & 0x1);
  }
}

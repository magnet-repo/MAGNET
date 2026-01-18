#include "poly_q.h"
#include "comp.h"
#include "poly_param.h"
#include "poly_red.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static const uint64_t omega_q[R] = {
    0,           25298509LL,  16654445LL,  105313793LL, 123808285LL,
    52751491LL,  98449108LL,  15315657LL,  7313337LL,   159716729LL,
    19434864LL,  64510017LL,  157082068LL, 138478139LL, 153167634LL,
    162331005LL, 93193167LL,  88710966LL,  142824326LL, 68055267LL,
    143926375LL, 2409315LL,   55525686LL,  1320151LL,   132830388LL,
    127359629LL, 77393886LL,  66119492LL,  115649674LL, 74322832LL,
    38402839LL,  133169436LL, 88335286LL,  99920198LL,  164068287LL,
    29239969LL,  68974339LL,  111338276LL, 117209794LL, 20738257LL,
    74667747LL,  127912654LL, 140345503LL, 44063037LL,  118114731LL,
    88157867LL,  112738182LL, 67712443LL,  123200480LL, 75774403LL,
    3140516LL,   12139007LL,  125420424LL, 140275625LL, 108225137LL,
    56438098LL,  128481306LL, 89488910LL,  165713817LL, 101098166LL,
    89492208LL,  117624798LL, 129032747LL, 56613722LL};

static const uint64_t perm_ntt[R] = {
    0,  32, 8,  44, 5,  39, 13, 40, 2,  35, 10, 47, 7,  37, 15, 43,
    1,  33, 9,  45, 4,  38, 12, 41, 3,  34, 11, 46, 6,  36, 14, 42,
    63, 31, 55, 19, 58, 24, 50, 23, 61, 28, 53, 16, 56, 26, 48, 20,
    62, 30, 54, 18, 59, 25, 51, 22, 60, 29, 52, 17, 57, 27, 49, 21};

static const uint64_t perm_inv_aut[D] = {
    0,   3,   6,   9,   12,  15,  18,  21,  24,  27,  30,  33,  36,  39,  42,
    45,  48,  51,  54,  57,  60,  63,  66,  69,  72,  75,  78,  81,  84,  87,
    90,  93,  96,  99,  102, 105, 108, 111, 114, 117, 120, 123, 126, 129, 132,
    135, 138, 141, 144, 147, 150, 153, 156, 159, 162, 165, 168, 171, 174, 177,
    180, 183, 186, 189, 192, 195, 198, 201, 204, 207, 210, 213, 216, 219, 222,
    225, 228, 231, 234, 237, 240, 243, 246, 249, 252, 255, 2,   5,   8,   11,
    14,  17,  20,  23,  26,  29,  32,  35,  38,  41,  44,  47,  50,  53,  56,
    59,  62,  65,  68,  71,  74,  77,  80,  83,  86,  89,  92,  95,  98,  101,
    104, 107, 110, 113, 116, 119, 122, 125, 128, 131, 134, 137, 140, 143, 146,
    149, 152, 155, 158, 161, 164, 167, 170, 173, 176, 179, 182, 185, 188, 191,
    194, 197, 200, 203, 206, 209, 212, 215, 218, 221, 224, 227, 230, 233, 236,
    239, 242, 245, 248, 251, 254, 1,   4,   7,   10,  13,  16,  19,  22,  25,
    28,  31,  34,  37,  40,  43,  46,  49,  52,  55,  58,  61,  64,  67,  70,
    73,  76,  79,  82,  85,  88,  91,  94,  97,  100, 103, 106, 109, 112, 115,
    118, 121, 124, 127, 130, 133, 136, 139, 142, 145, 148, 151, 154, 157, 160,
    163, 166, 169, 172, 175, 178, 181, 184, 187, 190, 193, 196, 199, 202, 205,
    208, 211, 214, 217, 220, 223, 226, 229, 232, 235, 238, 241, 244, 247, 250,
    253};

static const uint64_t omega_qhat[R] = {0,
                                       1783426224LL,
                                       2351192601LL,
                                       1746590616LL,
                                       3398309680LL,
                                       1748669978LL,
                                       1231514205LL,
                                       3013816622LL,
                                       1109357322LL,
                                       1001963140LL,
                                       119423006LL,
                                       2779991497LL,
                                       2908661274LL,
                                       992186722LL,
                                       167818709LL,
                                       1498269039LL,
                                       658934918LL,
                                       2914796113LL,
                                       3206411471LL,
                                       798669888LL,
                                       2269350985LL,
                                       1435531734LL,
                                       490445684LL,
                                       1342098185LL,
                                       1034782091LL,
                                       694861507LL,
                                       2287372836LL,
                                       3001658505LL,
                                       3120609358LL,
                                       2280115831LL,
                                       41486045LL,
                                       3177221765LL,
                                       2577100799LL,
                                       2405279510LL,
                                       2975400305LL,
                                       1030567863LL,
                                       2388041701LL,
                                       1001362987LL,
                                       1320506885LL,
                                       1880364562LL,
                                       280026387LL,
                                       166102516LL,
                                       3065396503LL,
                                       306400228LL,
                                       826469296LL,
                                       2257662681LL,
                                       946639311LL,
                                       3367646933LL,
                                       1105343673LL,
                                       1152204420LL,
                                       558218157LL,
                                       1054335247LL,
                                       1120818485LL,
                                       1201796820LL,
                                       1601846056LL,
                                       2569540653LL,
                                       38135503LL,
                                       362293383LL,
                                       1000961578LL,
                                       930299718LL,
                                       975114724LL,
                                       2593543626LL,
                                       789382750LL,
                                       2762722401LL};

/* in bit-reversed order */
static inline void ntt_q_bitrev(POLY_Q *a) {
  uint64_t t = D >> 1, m, i, j, j_first, j_last;
  uint64_t u, v;

  /* implicitly convert to Montgomery form and perform the butterfly operations
   */
  for (i = 0; i < D >> 1; i++) {
    u = montgomery_q(a->poly[i], MONTGOMERY_CONVERT_FACTOR_Q);
    v = montgomery_q(a->poly[i + (D >> 1)], omega_q[1]);

    a->poly[i] = con_sub(u + v, Q);
    a->poly[i + (D >> 1)] = con_add(u - v, Q);
  }

  for (m = 2; m < R; m <<= 1) {
    t >>= 1;
    for (i = 0; i < m; i++) {
      j_first = 2 * i * t;
      j_last = j_first + t;
      /* butterfly: (u,v)-->(u+v*omega,u-v*omega) */
      for (j = j_first; j < j_last; j++) {
        u = a->poly[j];
        v = montgomery_q(a->poly[j + t], omega_q[m + i]);

        a->poly[j] = con_sub(u + v, Q);
        a->poly[j + t] = con_add(u - v, Q);
      }
    }
  }
}

/* in the ordering (X - w^{(-1)^s * 3^r % 128}) for s = {0..1}, r = {0..31} */
void ntt_q(POLY_Q *a) {
  POLY_Q tmp;

  uint64_t i, j;

  memcpy(&tmp, a, sizeof(POLY_Q));
  ntt_q_bitrev(&tmp);

  for (i = 0; i < R; i++) {
    for (j = 0; j < SUBRING_SIZE; j++) {
      a->poly[i * SUBRING_SIZE + j] = tmp.poly[perm_ntt[i] * SUBRING_SIZE + j];
    }
  }
}

static inline void mult_4_plus_q(uint64_t *out, const uint64_t *a,
                                 const uint64_t *b, const uint64_t omega) {
  out[0] = red_short_q(montgomery_q(a[0], b[0]) +
                       montgomery_q(montgomery_q(a[1], b[3]) +
                                        montgomery_q(a[2], b[2]) +
                                        montgomery_q(a[3], b[1]),
                                    omega));
  out[1] = red_short_q(
      montgomery_q(a[0], b[1]) + montgomery_q(a[1], b[0]) +
      montgomery_q(montgomery_q(a[2], b[3]) + montgomery_q(a[3], b[2]), omega));
  out[2] = red_short_q(montgomery_q(a[0], b[2]) + montgomery_q(a[1], b[1]) +
                       montgomery_q(a[2], b[0]) +
                       montgomery_q(montgomery_q(a[3], b[3]), omega));
  out[3] = red_short_q(montgomery_q(a[0], b[3]) + montgomery_q(a[1], b[2]) +
                       montgomery_q(a[2], b[1]) + montgomery_q(a[3], b[0]));
}

static inline void mult_4_minus_q(uint64_t *out, const uint64_t *a,
                                  const uint64_t *b, const uint64_t omega) {
  out[0] = red_short_q(montgomery_q(a[0], b[0]) + Q -
                       montgomery_q(montgomery_q(a[1], b[3]) +
                                        montgomery_q(a[2], b[2]) +
                                        montgomery_q(a[3], b[1]),
                                    omega));
  out[1] = red_short_q(
      montgomery_q(a[0], b[1]) + montgomery_q(a[1], b[0]) + Q -
      montgomery_q(montgomery_q(a[2], b[3]) + montgomery_q(a[3], b[2]), omega));
  out[2] = red_short_q(montgomery_q(a[0], b[2]) + montgomery_q(a[1], b[1]) +
                       montgomery_q(a[2], b[0]) + Q -
                       montgomery_q(montgomery_q(a[3], b[3]), omega));
  out[3] = red_short_q(montgomery_q(a[0], b[3]) + montgomery_q(a[1], b[2]) +
                       montgomery_q(a[2], b[1]) + montgomery_q(a[3], b[0]));
}

/* in the ordering (X - w^{(-1)^s * 3^r % 128}) for s = {0..1}, r = {0..31} */
void mult_rq(POLY_Q *out, const POLY_Q *a, const POLY_Q *b) {
  mult_4_plus_q(out->poly, a->poly, b->poly, omega_q[32]);
  mult_4_plus_q(out->poly + SUBRING_SIZE, a->poly + SUBRING_SIZE,
                b->poly + SUBRING_SIZE, omega_q[48]);
  mult_4_plus_q(out->poly + 2 * SUBRING_SIZE, a->poly + 2 * SUBRING_SIZE,
                b->poly + 2 * SUBRING_SIZE, omega_q[36]);
  mult_4_plus_q(out->poly + 3 * SUBRING_SIZE, a->poly + 3 * SUBRING_SIZE,
                b->poly + 3 * SUBRING_SIZE, omega_q[54]);
  mult_4_minus_q(out->poly + 4 * SUBRING_SIZE, a->poly + 4 * SUBRING_SIZE,
                 b->poly + 4 * SUBRING_SIZE, omega_q[34]);
  mult_4_minus_q(out->poly + 5 * SUBRING_SIZE, a->poly + 5 * SUBRING_SIZE,
                 b->poly + 5 * SUBRING_SIZE, omega_q[51]);
  mult_4_minus_q(out->poly + 6 * SUBRING_SIZE, a->poly + 6 * SUBRING_SIZE,
                 b->poly + 6 * SUBRING_SIZE, omega_q[38]);
  mult_4_plus_q(out->poly + 7 * SUBRING_SIZE, a->poly + 7 * SUBRING_SIZE,
                b->poly + 7 * SUBRING_SIZE, omega_q[52]);
  mult_4_plus_q(out->poly + 8 * SUBRING_SIZE, a->poly + 8 * SUBRING_SIZE,
                b->poly + 8 * SUBRING_SIZE, omega_q[33]);
  mult_4_minus_q(out->poly + 9 * SUBRING_SIZE, a->poly + 9 * SUBRING_SIZE,
                 b->poly + 9 * SUBRING_SIZE, omega_q[49]);
  mult_4_plus_q(out->poly + 10 * SUBRING_SIZE, a->poly + 10 * SUBRING_SIZE,
                b->poly + 10 * SUBRING_SIZE, omega_q[37]);
  mult_4_minus_q(out->poly + 11 * SUBRING_SIZE, a->poly + 11 * SUBRING_SIZE,
                 b->poly + 11 * SUBRING_SIZE, omega_q[55]);
  mult_4_minus_q(out->poly + 12 * SUBRING_SIZE, a->poly + 12 * SUBRING_SIZE,
                 b->poly + 12 * SUBRING_SIZE, omega_q[35]);
  mult_4_minus_q(out->poly + 13 * SUBRING_SIZE, a->poly + 13 * SUBRING_SIZE,
                 b->poly + 13 * SUBRING_SIZE, omega_q[50]);
  mult_4_minus_q(out->poly + 14 * SUBRING_SIZE, a->poly + 14 * SUBRING_SIZE,
                 b->poly + 14 * SUBRING_SIZE, omega_q[39]);
  mult_4_minus_q(out->poly + 15 * SUBRING_SIZE, a->poly + 15 * SUBRING_SIZE,
                 b->poly + 15 * SUBRING_SIZE, omega_q[53]);
  mult_4_minus_q(out->poly + 16 * SUBRING_SIZE, a->poly + 16 * SUBRING_SIZE,
                 b->poly + 16 * SUBRING_SIZE, omega_q[32]);
  mult_4_minus_q(out->poly + 17 * SUBRING_SIZE, a->poly + 17 * SUBRING_SIZE,
                 b->poly + 17 * SUBRING_SIZE, omega_q[48]);
  mult_4_minus_q(out->poly + 18 * SUBRING_SIZE, a->poly + 18 * SUBRING_SIZE,
                 b->poly + 18 * SUBRING_SIZE, omega_q[36]);
  mult_4_minus_q(out->poly + 19 * SUBRING_SIZE, a->poly + 19 * SUBRING_SIZE,
                 b->poly + 19 * SUBRING_SIZE, omega_q[54]);
  mult_4_plus_q(out->poly + 20 * SUBRING_SIZE, a->poly + 20 * SUBRING_SIZE,
                b->poly + 20 * SUBRING_SIZE, omega_q[34]);
  mult_4_plus_q(out->poly + 21 * SUBRING_SIZE, a->poly + 21 * SUBRING_SIZE,
                b->poly + 21 * SUBRING_SIZE, omega_q[51]);
  mult_4_plus_q(out->poly + 22 * SUBRING_SIZE, a->poly + 22 * SUBRING_SIZE,
                b->poly + 22 * SUBRING_SIZE, omega_q[38]);
  mult_4_minus_q(out->poly + 23 * SUBRING_SIZE, a->poly + 23 * SUBRING_SIZE,
                 b->poly + 23 * SUBRING_SIZE, omega_q[52]);
  mult_4_minus_q(out->poly + 24 * SUBRING_SIZE, a->poly + 24 * SUBRING_SIZE,
                 b->poly + 24 * SUBRING_SIZE, omega_q[33]);
  mult_4_plus_q(out->poly + 25 * SUBRING_SIZE, a->poly + 25 * SUBRING_SIZE,
                b->poly + 25 * SUBRING_SIZE, omega_q[49]);
  mult_4_minus_q(out->poly + 26 * SUBRING_SIZE, a->poly + 26 * SUBRING_SIZE,
                 b->poly + 26 * SUBRING_SIZE, omega_q[37]);
  mult_4_plus_q(out->poly + 27 * SUBRING_SIZE, a->poly + 27 * SUBRING_SIZE,
                b->poly + 27 * SUBRING_SIZE, omega_q[55]);
  mult_4_plus_q(out->poly + 28 * SUBRING_SIZE, a->poly + 28 * SUBRING_SIZE,
                b->poly + 28 * SUBRING_SIZE, omega_q[35]);
  mult_4_plus_q(out->poly + 29 * SUBRING_SIZE, a->poly + 29 * SUBRING_SIZE,
                b->poly + 29 * SUBRING_SIZE, omega_q[50]);
  mult_4_plus_q(out->poly + 30 * SUBRING_SIZE, a->poly + 30 * SUBRING_SIZE,
                b->poly + 30 * SUBRING_SIZE, omega_q[39]);
  mult_4_plus_q(out->poly + 31 * SUBRING_SIZE, a->poly + 31 * SUBRING_SIZE,
                b->poly + 31 * SUBRING_SIZE, omega_q[53]);
  mult_4_minus_q(out->poly + 32 * SUBRING_SIZE, a->poly + 32 * SUBRING_SIZE,
                 b->poly + 32 * SUBRING_SIZE, omega_q[63]);
  mult_4_minus_q(out->poly + 33 * SUBRING_SIZE, a->poly + 33 * SUBRING_SIZE,
                 b->poly + 33 * SUBRING_SIZE, omega_q[47]);
  mult_4_minus_q(out->poly + 34 * SUBRING_SIZE, a->poly + 34 * SUBRING_SIZE,
                 b->poly + 34 * SUBRING_SIZE, omega_q[59]);
  mult_4_minus_q(out->poly + 35 * SUBRING_SIZE, a->poly + 35 * SUBRING_SIZE,
                 b->poly + 35 * SUBRING_SIZE, omega_q[41]);
  mult_4_plus_q(out->poly + 36 * SUBRING_SIZE, a->poly + 36 * SUBRING_SIZE,
                b->poly + 36 * SUBRING_SIZE, omega_q[61]);
  mult_4_plus_q(out->poly + 37 * SUBRING_SIZE, a->poly + 37 * SUBRING_SIZE,
                b->poly + 37 * SUBRING_SIZE, omega_q[44]);
  mult_4_plus_q(out->poly + 38 * SUBRING_SIZE, a->poly + 38 * SUBRING_SIZE,
                b->poly + 38 * SUBRING_SIZE, omega_q[57]);
  mult_4_minus_q(out->poly + 39 * SUBRING_SIZE, a->poly + 39 * SUBRING_SIZE,
                 b->poly + 39 * SUBRING_SIZE, omega_q[43]);
  mult_4_minus_q(out->poly + 40 * SUBRING_SIZE, a->poly + 40 * SUBRING_SIZE,
                 b->poly + 40 * SUBRING_SIZE, omega_q[62]);
  mult_4_plus_q(out->poly + 41 * SUBRING_SIZE, a->poly + 41 * SUBRING_SIZE,
                b->poly + 41 * SUBRING_SIZE, omega_q[46]);
  mult_4_minus_q(out->poly + 42 * SUBRING_SIZE, a->poly + 42 * SUBRING_SIZE,
                 b->poly + 42 * SUBRING_SIZE, omega_q[58]);
  mult_4_plus_q(out->poly + 43 * SUBRING_SIZE, a->poly + 43 * SUBRING_SIZE,
                b->poly + 43 * SUBRING_SIZE, omega_q[40]);
  mult_4_plus_q(out->poly + 44 * SUBRING_SIZE, a->poly + 44 * SUBRING_SIZE,
                b->poly + 44 * SUBRING_SIZE, omega_q[60]);
  mult_4_plus_q(out->poly + 45 * SUBRING_SIZE, a->poly + 45 * SUBRING_SIZE,
                b->poly + 45 * SUBRING_SIZE, omega_q[45]);
  mult_4_plus_q(out->poly + 46 * SUBRING_SIZE, a->poly + 46 * SUBRING_SIZE,
                b->poly + 46 * SUBRING_SIZE, omega_q[56]);
  mult_4_plus_q(out->poly + 47 * SUBRING_SIZE, a->poly + 47 * SUBRING_SIZE,
                b->poly + 47 * SUBRING_SIZE, omega_q[42]);
  mult_4_plus_q(out->poly + 48 * SUBRING_SIZE, a->poly + 48 * SUBRING_SIZE,
                b->poly + 48 * SUBRING_SIZE, omega_q[63]);
  mult_4_plus_q(out->poly + 49 * SUBRING_SIZE, a->poly + 49 * SUBRING_SIZE,
                b->poly + 49 * SUBRING_SIZE, omega_q[47]);
  mult_4_plus_q(out->poly + 50 * SUBRING_SIZE, a->poly + 50 * SUBRING_SIZE,
                b->poly + 50 * SUBRING_SIZE, omega_q[59]);
  mult_4_plus_q(out->poly + 51 * SUBRING_SIZE, a->poly + 51 * SUBRING_SIZE,
                b->poly + 51 * SUBRING_SIZE, omega_q[41]);
  mult_4_minus_q(out->poly + 52 * SUBRING_SIZE, a->poly + 52 * SUBRING_SIZE,
                 b->poly + 52 * SUBRING_SIZE, omega_q[61]);
  mult_4_minus_q(out->poly + 53 * SUBRING_SIZE, a->poly + 53 * SUBRING_SIZE,
                 b->poly + 53 * SUBRING_SIZE, omega_q[44]);
  mult_4_minus_q(out->poly + 54 * SUBRING_SIZE, a->poly + 54 * SUBRING_SIZE,
                 b->poly + 54 * SUBRING_SIZE, omega_q[57]);
  mult_4_plus_q(out->poly + 55 * SUBRING_SIZE, a->poly + 55 * SUBRING_SIZE,
                b->poly + 55 * SUBRING_SIZE, omega_q[43]);
  mult_4_plus_q(out->poly + 56 * SUBRING_SIZE, a->poly + 56 * SUBRING_SIZE,
                b->poly + 56 * SUBRING_SIZE, omega_q[62]);
  mult_4_minus_q(out->poly + 57 * SUBRING_SIZE, a->poly + 57 * SUBRING_SIZE,
                 b->poly + 57 * SUBRING_SIZE, omega_q[46]);
  mult_4_plus_q(out->poly + 58 * SUBRING_SIZE, a->poly + 58 * SUBRING_SIZE,
                b->poly + 58 * SUBRING_SIZE, omega_q[58]);
  mult_4_minus_q(out->poly + 59 * SUBRING_SIZE, a->poly + 59 * SUBRING_SIZE,
                 b->poly + 59 * SUBRING_SIZE, omega_q[40]);
  mult_4_minus_q(out->poly + 60 * SUBRING_SIZE, a->poly + 60 * SUBRING_SIZE,
                 b->poly + 60 * SUBRING_SIZE, omega_q[60]);
  mult_4_minus_q(out->poly + 61 * SUBRING_SIZE, a->poly + 61 * SUBRING_SIZE,
                 b->poly + 61 * SUBRING_SIZE, omega_q[45]);
  mult_4_minus_q(out->poly + 62 * SUBRING_SIZE, a->poly + 62 * SUBRING_SIZE,
                 b->poly + 62 * SUBRING_SIZE, omega_q[56]);
  mult_4_minus_q(out->poly + 63 * SUBRING_SIZE, a->poly + 63 * SUBRING_SIZE,
                 b->poly + 63 * SUBRING_SIZE, omega_q[42]);
}

void aut_ntt_plus(uint64_t *out, const uint64_t *a, const uint64_t omega_pos) {
  out[0] = a[0];
  out[1] = montgomery_q(
      a[3], con_add((1 - 2 * (omega_pos & 0x1)) * omega_q[omega_pos >> 1], Q));
  out[2] = montgomery_q(a[2], omega_q[omega_pos]);
  out[3] = a[1];
}

void aut_ntt_minus(uint64_t *out, const uint64_t *a, const uint64_t omega_pos) {
  out[0] = a[0];
  out[1] = montgomery_q(
      a[3], con_add((1 - 2 * (omega_pos & 0x1)) * omega_q[omega_pos >> 1], Q));
  out[2] = Q - montgomery_q(a[2], omega_q[omega_pos]);
  out[3] = a[1];
}

/* automorphism: \sigma_3(a) (on NTT domain) */
void aut(POLY_Q *out, const POLY_Q *a) {
  POLY_Q tmp;

  uint64_t i;

  for (i = 0; i < (D >> 1) - SUBRING_SIZE; i++) {
    tmp.poly[i] = (a->poly)[i + SUBRING_SIZE];
    tmp.poly[i + (D >> 1)] = (a->poly)[i + (D >> 1) + SUBRING_SIZE];
  }
  for (i = 0; i < SUBRING_SIZE; i++) {
    tmp.poly[(D >> 1) - SUBRING_SIZE + i] = (a->poly)[i];
    tmp.poly[D - SUBRING_SIZE + i] = (a->poly)[i + (D >> 1)];
  }

  aut_ntt_plus(out->poly, tmp.poly, 32);
  aut_ntt_plus(out->poly + SUBRING_SIZE, tmp.poly + SUBRING_SIZE, 48);
  aut_ntt_plus(out->poly + 2 * SUBRING_SIZE, tmp.poly + 2 * SUBRING_SIZE, 36);
  aut_ntt_plus(out->poly + 3 * SUBRING_SIZE, tmp.poly + 3 * SUBRING_SIZE, 54);
  aut_ntt_minus(out->poly + 4 * SUBRING_SIZE, tmp.poly + 4 * SUBRING_SIZE, 34);
  aut_ntt_minus(out->poly + 5 * SUBRING_SIZE, tmp.poly + 5 * SUBRING_SIZE, 51);
  aut_ntt_minus(out->poly + 6 * SUBRING_SIZE, tmp.poly + 6 * SUBRING_SIZE, 38);
  aut_ntt_plus(out->poly + 7 * SUBRING_SIZE, tmp.poly + 7 * SUBRING_SIZE, 52);
  aut_ntt_plus(out->poly + 8 * SUBRING_SIZE, tmp.poly + 8 * SUBRING_SIZE, 33);
  aut_ntt_minus(out->poly + 9 * SUBRING_SIZE, tmp.poly + 9 * SUBRING_SIZE, 49);
  aut_ntt_plus(out->poly + 10 * SUBRING_SIZE, tmp.poly + 10 * SUBRING_SIZE, 37);
  aut_ntt_minus(out->poly + 11 * SUBRING_SIZE, tmp.poly + 11 * SUBRING_SIZE,
                55);
  aut_ntt_minus(out->poly + 12 * SUBRING_SIZE, tmp.poly + 12 * SUBRING_SIZE,
                35);
  aut_ntt_minus(out->poly + 13 * SUBRING_SIZE, tmp.poly + 13 * SUBRING_SIZE,
                50);
  aut_ntt_minus(out->poly + 14 * SUBRING_SIZE, tmp.poly + 14 * SUBRING_SIZE,
                39);
  aut_ntt_minus(out->poly + 15 * SUBRING_SIZE, tmp.poly + 15 * SUBRING_SIZE,
                53);
  aut_ntt_minus(out->poly + 16 * SUBRING_SIZE, tmp.poly + 16 * SUBRING_SIZE,
                32);
  aut_ntt_minus(out->poly + 17 * SUBRING_SIZE, tmp.poly + 17 * SUBRING_SIZE,
                48);
  aut_ntt_minus(out->poly + 18 * SUBRING_SIZE, tmp.poly + 18 * SUBRING_SIZE,
                36);
  aut_ntt_minus(out->poly + 19 * SUBRING_SIZE, tmp.poly + 19 * SUBRING_SIZE,
                54);
  aut_ntt_plus(out->poly + 20 * SUBRING_SIZE, tmp.poly + 20 * SUBRING_SIZE, 34);
  aut_ntt_plus(out->poly + 21 * SUBRING_SIZE, tmp.poly + 21 * SUBRING_SIZE, 51);
  aut_ntt_plus(out->poly + 22 * SUBRING_SIZE, tmp.poly + 22 * SUBRING_SIZE, 38);
  aut_ntt_minus(out->poly + 23 * SUBRING_SIZE, tmp.poly + 23 * SUBRING_SIZE,
                52);
  aut_ntt_minus(out->poly + 24 * SUBRING_SIZE, tmp.poly + 24 * SUBRING_SIZE,
                33);
  aut_ntt_plus(out->poly + 25 * SUBRING_SIZE, tmp.poly + 25 * SUBRING_SIZE, 49);
  aut_ntt_minus(out->poly + 26 * SUBRING_SIZE, tmp.poly + 26 * SUBRING_SIZE,
                37);
  aut_ntt_plus(out->poly + 27 * SUBRING_SIZE, tmp.poly + 27 * SUBRING_SIZE, 55);
  aut_ntt_plus(out->poly + 28 * SUBRING_SIZE, tmp.poly + 28 * SUBRING_SIZE, 35);
  aut_ntt_plus(out->poly + 29 * SUBRING_SIZE, tmp.poly + 29 * SUBRING_SIZE, 50);
  aut_ntt_plus(out->poly + 30 * SUBRING_SIZE, tmp.poly + 30 * SUBRING_SIZE, 39);
  aut_ntt_plus(out->poly + 31 * SUBRING_SIZE, tmp.poly + 31 * SUBRING_SIZE, 53);
  aut_ntt_minus(out->poly + 32 * SUBRING_SIZE, tmp.poly + 32 * SUBRING_SIZE,
                63);
  aut_ntt_minus(out->poly + 33 * SUBRING_SIZE, tmp.poly + 33 * SUBRING_SIZE,
                47);
  aut_ntt_minus(out->poly + 34 * SUBRING_SIZE, tmp.poly + 34 * SUBRING_SIZE,
                59);
  aut_ntt_minus(out->poly + 35 * SUBRING_SIZE, tmp.poly + 35 * SUBRING_SIZE,
                41);
  aut_ntt_plus(out->poly + 36 * SUBRING_SIZE, tmp.poly + 36 * SUBRING_SIZE, 61);
  aut_ntt_plus(out->poly + 37 * SUBRING_SIZE, tmp.poly + 37 * SUBRING_SIZE, 44);
  aut_ntt_plus(out->poly + 38 * SUBRING_SIZE, tmp.poly + 38 * SUBRING_SIZE, 57);
  aut_ntt_minus(out->poly + 39 * SUBRING_SIZE, tmp.poly + 39 * SUBRING_SIZE,
                43);
  aut_ntt_minus(out->poly + 40 * SUBRING_SIZE, tmp.poly + 40 * SUBRING_SIZE,
                62);
  aut_ntt_plus(out->poly + 41 * SUBRING_SIZE, tmp.poly + 41 * SUBRING_SIZE, 46);
  aut_ntt_minus(out->poly + 42 * SUBRING_SIZE, tmp.poly + 42 * SUBRING_SIZE,
                58);
  aut_ntt_plus(out->poly + 43 * SUBRING_SIZE, tmp.poly + 43 * SUBRING_SIZE, 40);
  aut_ntt_plus(out->poly + 44 * SUBRING_SIZE, tmp.poly + 44 * SUBRING_SIZE, 60);
  aut_ntt_plus(out->poly + 45 * SUBRING_SIZE, tmp.poly + 45 * SUBRING_SIZE, 45);
  aut_ntt_plus(out->poly + 46 * SUBRING_SIZE, tmp.poly + 46 * SUBRING_SIZE, 56);
  aut_ntt_plus(out->poly + 47 * SUBRING_SIZE, tmp.poly + 47 * SUBRING_SIZE, 42);
  aut_ntt_plus(out->poly + 48 * SUBRING_SIZE, tmp.poly + 48 * SUBRING_SIZE, 63);
  aut_ntt_plus(out->poly + 49 * SUBRING_SIZE, tmp.poly + 49 * SUBRING_SIZE, 47);
  aut_ntt_plus(out->poly + 50 * SUBRING_SIZE, tmp.poly + 50 * SUBRING_SIZE, 59);
  aut_ntt_plus(out->poly + 51 * SUBRING_SIZE, tmp.poly + 51 * SUBRING_SIZE, 41);
  aut_ntt_minus(out->poly + 52 * SUBRING_SIZE, tmp.poly + 52 * SUBRING_SIZE,
                61);
  aut_ntt_minus(out->poly + 53 * SUBRING_SIZE, tmp.poly + 53 * SUBRING_SIZE,
                44);
  aut_ntt_minus(out->poly + 54 * SUBRING_SIZE, tmp.poly + 54 * SUBRING_SIZE,
                57);
  aut_ntt_plus(out->poly + 55 * SUBRING_SIZE, tmp.poly + 55 * SUBRING_SIZE, 43);
  aut_ntt_plus(out->poly + 56 * SUBRING_SIZE, tmp.poly + 56 * SUBRING_SIZE, 62);
  aut_ntt_minus(out->poly + 57 * SUBRING_SIZE, tmp.poly + 57 * SUBRING_SIZE,
                46);
  aut_ntt_plus(out->poly + 58 * SUBRING_SIZE, tmp.poly + 58 * SUBRING_SIZE, 58);
  aut_ntt_minus(out->poly + 59 * SUBRING_SIZE, tmp.poly + 59 * SUBRING_SIZE,
                40);
  aut_ntt_minus(out->poly + 60 * SUBRING_SIZE, tmp.poly + 60 * SUBRING_SIZE,
                60);
  aut_ntt_minus(out->poly + 61 * SUBRING_SIZE, tmp.poly + 61 * SUBRING_SIZE,
                45);
  aut_ntt_minus(out->poly + 62 * SUBRING_SIZE, tmp.poly + 62 * SUBRING_SIZE,
                56);
  aut_ntt_minus(out->poly + 63 * SUBRING_SIZE, tmp.poly + 63 * SUBRING_SIZE,
                42);
}

/* automorphism inverse: \sigma_{171}(a) (in R) */
void iaut_r(POLY_R *out, const POLY_R *a) {
  uint64_t i;

  for (i = 0; i < D; i++) {
    (out->poly)[i] = (a->poly)[perm_inv_aut[i]];
  }
  for (i = 86; i < 171; i++) {
    (out->poly)[i] = -(out->poly)[i];
  }
}

/* in bit-reversed order (seems there is no need to permute) */
void ntt_qhat(POLY_QHAT *a) {
  uint64_t t = D >> 1, m, i, j, j_first, j_last;
  uint64_t u, v;

  /* implicitly convert to Montgomery form and perform the butterfly operations
   */
  for (i = 0; i < D >> 1; i++) {
    u = montgomery_qhat(a->poly[i], MONTGOMERY_CONVERT_FACTOR_QHAT);
    v = montgomery_qhat(a->poly[i + (D >> 1)], omega_qhat[1]);

    a->poly[i] = con_sub(u + v, QHAT);
    a->poly[i + (D >> 1)] = con_add(u - v, QHAT);
  }

  for (m = 2; m < R; m <<= 1) {
    t >>= 1;
    for (i = 0; i < m; i++) {
      j_first = 2 * i * t;
      j_last = j_first + t;
      /* butterfly: (u,v)-->(u+v*omega,u-v*omega) */
      for (j = j_first; j < j_last; j++) {
        u = a->poly[j];
        v = montgomery_qhat(a->poly[j + t], omega_qhat[m + i]);

        a->poly[j] = con_sub(u + v, QHAT);
        a->poly[j + t] = con_add(u - v, QHAT);
      }
    }
  }
}

static inline void mult_4_plus_qhat(uint64_t *out, const uint64_t *a,
                                    const uint64_t *b, const uint64_t omega) {
  out[0] = red_short_qhat(
      montgomery_qhat(a[0], b[0]) +
      montgomery_qhat(red_short_qhat(montgomery_qhat(a[1], b[3]) +
                                     montgomery_qhat(a[2], b[2]) +
                                     montgomery_qhat(a[3], b[1])),
                      omega));
  out[1] = red_short_qhat(
      montgomery_qhat(a[0], b[1]) + montgomery_qhat(a[1], b[0]) +
      montgomery_qhat(red_short_qhat(montgomery_qhat(a[2], b[3]) +
                                     montgomery_qhat(a[3], b[2])),
                      omega));
  out[2] =
      red_short_qhat(montgomery_qhat(a[0], b[2]) + montgomery_qhat(a[1], b[1]) +
                     montgomery_qhat(a[2], b[0]) +
                     montgomery_qhat(montgomery_qhat(a[3], b[3]), omega));
  out[3] =
      red_short_qhat(montgomery_qhat(a[0], b[3]) + montgomery_qhat(a[1], b[2]) +
                     montgomery_qhat(a[2], b[1]) + montgomery_qhat(a[3], b[0]));
}

static inline void mult_4_minus_qhat(uint64_t *out, const uint64_t *a,
                                     const uint64_t *b, const uint64_t omega) {
  out[0] = red_short_qhat(
      montgomery_qhat(a[0], b[0]) + QHAT -
      montgomery_qhat(red_short_qhat(montgomery_qhat(a[1], b[3]) +
                                     montgomery_qhat(a[2], b[2]) +
                                     montgomery_qhat(a[3], b[1])),
                      omega));
  out[1] = red_short_qhat(
      montgomery_qhat(a[0], b[1]) + montgomery_qhat(a[1], b[0]) + QHAT -
      montgomery_qhat(red_short_qhat(montgomery_qhat(a[2], b[3]) +
                                     montgomery_qhat(a[3], b[2])),
                      omega));
  out[2] =
      red_short_qhat(montgomery_qhat(a[0], b[2]) + montgomery_qhat(a[1], b[1]) +
                     montgomery_qhat(a[2], b[0]) + QHAT -
                     montgomery_qhat(montgomery_qhat(a[3], b[3]), omega));
  out[3] =
      red_short_qhat(montgomery_qhat(a[0], b[3]) + montgomery_qhat(a[1], b[2]) +
                     montgomery_qhat(a[2], b[1]) + montgomery_qhat(a[3], b[0]));
}

/* in bit-reversed order */
void mult_rqhat(POLY_QHAT *out, const POLY_QHAT *a, const POLY_QHAT *b) {
  uint64_t i;

  for (i = 0; i < (R >> 1); i++) {
    mult_4_plus_qhat(out->poly + i * SUBRING_SIZE * 2,
                     a->poly + i * SUBRING_SIZE * 2,
                     b->poly + i * SUBRING_SIZE * 2, omega_qhat[(R >> 1) + i]);
    mult_4_minus_qhat(out->poly + i * SUBRING_SIZE * 2 + SUBRING_SIZE,
                      a->poly + i * SUBRING_SIZE * 2 + SUBRING_SIZE,
                      b->poly + i * SUBRING_SIZE * 2 + SUBRING_SIZE,
                      omega_qhat[(R >> 1) + i]);
  }
}

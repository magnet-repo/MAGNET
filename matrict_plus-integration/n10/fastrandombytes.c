/* Adapted from Intel® Advanced Encryption Standard (Intel® AES) Instructions
 * Set - Rev 3.01
 * https://software.intel.com/sites/default/files/article/165683/aes-wp-2012-09-22-v01.pdf
 */

#include "fastrandombytes.h"
#include <string.h>
#include <x86intrin.h>

static __m128i round_key_pub[15];
static __m128i round_key_prv[15];
static __m128i round_key_tmp[15];
static __m128i round_key_ch[15];

static __m128i iv_pub;
static __m128i iv_prv;
static __m128i iv_tmp;
static __m128i iv_ch;

static const __m128i ONE = {1, 0};

/* we put the input string in one part of the block, and the counter in the
 * other part */
static const __m128i IV_Mat1 = {0, 0x4D617431};
static const __m128i IV_Mat2 = {0, 0x4D617432};
static const __m128i IV_Mat3 = {0, 0x4D617433};
static const __m128i IV_Ch = {0, 0x00004368};

static inline void KEY_256_ASSIST_1(__m128i *temp1, __m128i *temp2) {
  __m128i temp4;
  *temp2 = _mm_shuffle_epi32(*temp2, 0xff);
  temp4 = _mm_slli_si128(*temp1, 0x4);
  *temp1 = _mm_xor_si128(*temp1, temp4);
  temp4 = _mm_slli_si128(temp4, 0x4);
  *temp1 = _mm_xor_si128(*temp1, temp4);
  temp4 = _mm_slli_si128(temp4, 0x4);
  *temp1 = _mm_xor_si128(*temp1, temp4);
  *temp1 = _mm_xor_si128(*temp1, *temp2);
}

static inline void KEY_256_ASSIST_2(__m128i *temp1, __m128i *temp3) {
  __m128i temp2, temp4;
  temp4 = _mm_aeskeygenassist_si128(*temp1, 0x0);
  temp2 = _mm_shuffle_epi32(temp4, 0xaa);
  temp4 = _mm_slli_si128(*temp3, 0x4);
  *temp3 = _mm_xor_si128(*temp3, temp4);
  temp4 = _mm_slli_si128(temp4, 0x4);
  *temp3 = _mm_xor_si128(*temp3, temp4);
  temp4 = _mm_slli_si128(temp4, 0x4);
  *temp3 = _mm_xor_si128(*temp3, temp4);
  *temp3 = _mm_xor_si128(*temp3, temp2);
}

/* round_key <-- aes256_key_expansion(randomness) */
static inline void fastrandombytes_setseed(__m128i *round_key,
                                           const unsigned char *randomness) {
  __m128i temp1, temp2, temp3;

  temp1 = _mm_loadu_si128((__m128i *)randomness);
  temp3 = _mm_loadu_si128((__m128i *)(randomness + 16));
  round_key[0] = temp1;
  round_key[1] = temp3;
  temp2 = _mm_aeskeygenassist_si128(temp3, 0x01);
  KEY_256_ASSIST_1(&temp1, &temp2);
  round_key[2] = temp1;
  KEY_256_ASSIST_2(&temp1, &temp3);
  round_key[3] = temp3;
  temp2 = _mm_aeskeygenassist_si128(temp3, 0x02);
  KEY_256_ASSIST_1(&temp1, &temp2);
  round_key[4] = temp1;
  KEY_256_ASSIST_2(&temp1, &temp3);
  round_key[5] = temp3;
  temp2 = _mm_aeskeygenassist_si128(temp3, 0x04);
  KEY_256_ASSIST_1(&temp1, &temp2);
  round_key[6] = temp1;
  KEY_256_ASSIST_2(&temp1, &temp3);
  round_key[7] = temp3;
  temp2 = _mm_aeskeygenassist_si128(temp3, 0x08);
  KEY_256_ASSIST_1(&temp1, &temp2);
  round_key[8] = temp1;
  KEY_256_ASSIST_2(&temp1, &temp3);
  round_key[9] = temp3;
  temp2 = _mm_aeskeygenassist_si128(temp3, 0x10);
  KEY_256_ASSIST_1(&temp1, &temp2);
  round_key[10] = temp1;
  KEY_256_ASSIST_2(&temp1, &temp3);
  round_key[11] = temp3;
  temp2 = _mm_aeskeygenassist_si128(temp3, 0x20);
  KEY_256_ASSIST_1(&temp1, &temp2);
  round_key[12] = temp1;
  KEY_256_ASSIST_2(&temp1, &temp3);
  round_key[13] = temp3;
  temp2 = _mm_aeskeygenassist_si128(temp3, 0x40);
  KEY_256_ASSIST_1(&temp1, &temp2);
  round_key[14] = temp1;
}

void fastrandombytes_setiv_mat1() { iv_pub = IV_Mat1; }

void fastrandombytes_setiv_mat2() { iv_pub = IV_Mat2; }

void fastrandombytes_setiv_mat3() { iv_pub = IV_Mat3; }

static inline void AES_ctr_round(__m128i *round_key, __m128i *iv,
                                 unsigned char *out) {
  __m128i tmp;

  tmp = _mm_xor_si128(*iv, round_key[0]);
  tmp = _mm_aesenc_si128(tmp, round_key[1]);
  tmp = _mm_aesenc_si128(tmp, round_key[2]);
  tmp = _mm_aesenc_si128(tmp, round_key[3]);
  tmp = _mm_aesenc_si128(tmp, round_key[4]);
  tmp = _mm_aesenc_si128(tmp, round_key[5]);
  tmp = _mm_aesenc_si128(tmp, round_key[6]);
  tmp = _mm_aesenc_si128(tmp, round_key[7]);
  tmp = _mm_aesenc_si128(tmp, round_key[8]);
  tmp = _mm_aesenc_si128(tmp, round_key[9]);
  tmp = _mm_aesenc_si128(tmp, round_key[10]);
  tmp = _mm_aesenc_si128(tmp, round_key[11]);
  tmp = _mm_aesenc_si128(tmp, round_key[12]);
  tmp = _mm_aesenc_si128(tmp, round_key[13]);
  tmp = _mm_aesenclast_si128(tmp, round_key[14]);
  _mm_storeu_si128((__m128i *)out, tmp);

  *iv = _mm_add_epi32(*iv, ONE);
}

/* r <-- aes256_ctr(round_key, iv, rlen) */
static inline void fastrandombytes(__m128i *round_key, __m128i *iv,
                                   unsigned char *r, unsigned long long rlen) {
  unsigned char ct[16];
  unsigned long long num_of_blocks = rlen >> 4;
  unsigned long long i;

  for (i = 0; i < num_of_blocks; i++) {
    AES_ctr_round(round_key, iv, r + (i << 4));
  }

  if (rlen & 0x0f) {
    AES_ctr_round(round_key, iv, ct);

    memcpy(r + (i << 4), ct, rlen & 0x0f);
  }
}

void fastrandombytes_setseed_pub(const unsigned char *randomness) {
  fastrandombytes_setseed(round_key_pub, randomness);
}

void fastrandombytes_setseed_prv(const unsigned char *randomness) {
  fastrandombytes_setseed(round_key_prv, randomness);
  iv_prv = _mm_setzero_si128();
}

void fastrandombytes_setseed_tmp(const unsigned char *randomness) {
  fastrandombytes_setseed(round_key_tmp, randomness);
  iv_tmp = _mm_setzero_si128();
}

void fastrandombytes_setseed_ch(const unsigned char *randomness) {
  fastrandombytes_setseed(round_key_tmp, randomness);
  iv_ch = IV_Ch;
}

void fastrandombytes_pub(unsigned char *r, unsigned long long rlen) {
  fastrandombytes(round_key_pub, &iv_pub, r, rlen);
}

void fastrandombytes_prv(unsigned char *r, unsigned long long rlen) {
  fastrandombytes(round_key_prv, &iv_prv, r, rlen);
}

void fastrandombytes_tmp(unsigned char *r, unsigned long long rlen) {
  fastrandombytes(round_key_tmp, &iv_tmp, r, rlen);
}

void fastrandombytes_ch(unsigned char *r, unsigned long long rlen) {
  fastrandombytes(round_key_ch, &iv_ch, r, rlen);
}

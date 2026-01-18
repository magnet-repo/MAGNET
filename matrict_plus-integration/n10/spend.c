#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "comp.h"
#include "fastrandombytes.h"
#include "gaussian_avx.h"
#include "hash.h"
#include "mint.h"
#include "param.h"
#include "poly_mult.h"
#include "poly_q.h"
#include "poly_red.h"
#include "sampleb.h"
#include "spend.h"

#include "masked_algorithms.h"

#include "params-GR19.h"
#include "utils-GR19.h"
#include "base_gadgets.h"
#include "sign_gadgets.h"

#define CDT_SIZE 52694

/* for M = 2, S = 2 */
void spend(SPEND_OUT *out, CN_OUT *cn_out, CK_OUT *ck_out, SN_OUT *sn_out,
           const ACT a_in[][M], const uint64_t l, const ASK *ask,
           const PK_OUT *pk_out, const uint64_t *amt_out,
           const POLY_Q mat1[][N], const POLY_Q mat2[][3],
           const POLY_QHAT mat3[][N_HAT], unsigned char *pp) {
  uint64_t i, j, k, rej, rej1;

  static uint64_t a_hat[S][R];
  static POLY_R y_out[S][N + K + 1];
  static POLY_Q y_out_ntt[S][N + K + 1];
  static POLY_Q w_out[S][N];
  static uint64_t c[R];
  static POLY_R y0[N + K + 1];
  static POLY_Q y0_ntt[N + K + 1];
  static POLY_Q w0[N];
  static POLY_R a[N_SPENT + 2];
  static POLY_Q a_ntt[N_SPENT + 2];
  static POLY_QHAT a_ntt_qhat[N_SPENT + 2];
  static POLY_R tmp_r;
  static POLY_Q t_in_sum[N_SPENT][N + 1];
  static POLY_R rho[N + K + 3];
  static POLY_Q rho_ntt[N + K + 3];
  static POLY_R y1[N + K + 3], y2[N + K + 3];
  static POLY_Q y1_ntt[N + K + 3], y2_ntt[N + K + 3];
  static POLY_Q w1[N], w2[N];
  static POLY_Q w;
  static POLY_Q tmp_q;
  static POLY_Q a3y1, a3y1_2;
  static POLY_Q b1yout[S];
  static POLY_Q alpha[S];
  static POLY_R alpha_prime_0;
  static POLY_Q v0;
  static POLY_Q alpha_b1yout[S];
  static POLY_Q pk_hat[N_SPENT][M][N_BAR];
  static POLY_Q g_gamma[N_BAR];
  static POLY_Q p[N_SPENT][N_BAR];
  static POLY_Q alpha_prime_0_ntt;
  static POLY_R s[N_BAR + K];
  static POLY_QHAT a2_ntt[N_SPENT + 2];
  static uint64_t delta[N_SPENT + 2];
  static POLY_R rb[N_HAT + K_HAT], ra[N_HAT + K_HAT];
  static POLY_QHAT rb_ntt[N_HAT + K_HAT], ra_ntt[N_HAT + K_HAT];
  static POLY_QHAT a_rs[N_HAT];
  static POLY_R y[N_BAR + K];
  static POLY_Q y_ntt[N_BAR + K];
  static POLY_Q e[N_BAR];
  static POLY_R g[N_SPENT + 2];
  static POLY_R f0;
  static POLY_Q tmp_aut;
  static POLY_R tmp_aut_r;
  static POLY_R x_delta[N_SPENT + 1];
  static POLY_R x_s[N_BAR + K];
  static POLY_R x_rb[N_HAT + K_HAT];
  static POLY_R x_rout[S][N + K + 1];
  static POLY_R x_ri[N + K + 1];
  static POLY_R x_rho[N + K + 3];
  static POLY_R x_aut_rho[N + K + 3];

  /* SN = (sn_0,...,sn_{M - 1}) */
  memcpy(sn_out->sn, ask->sn, CRYPTO_BYTES);
  memcpy((sn_out + 1)->sn, (ask + 1)->sn, CRYPTO_BYTES);

  /* (cn_{out, i}, ck_{out, i}, \hat{a}_i) <-- Mint(amt_{out, i}) */
  mint(cn_out->cn, ck_out->ck, a_hat[0], amt_out[0], mat1, mat2);
  mint((cn_out + 1)->cn, (ck_out + 1)->ck, a_hat[1], amt_out[1], mat1, mat2);

  /* c_{i + 1} = (c_i + \sum_{j = 0}^{S - 1} amt_{out, j}[i] - \sum_{j = 0}^{M -
   * 1} amt_{in, j}[i]) / 2 */
  c[0] = 0;
  for (i = 0; i < R - 1; i++) {
    c[i + 1] =
        ((int64_t)(c[i] + ((amt_out[0] >> i) & 0x1) +
                   ((amt_out[1] >> i) & 0x1) - ((ask[0].amt_in >> i) & 0x1) -
                   ((ask[1].amt_in >> i) & 0x1))) >>
        1;
  }

  /* \sum_{i = 0}^{M - 1} t_{in, j}^(i) */
  for (i = 0; i < N_SPENT; i++) {
    for (j = 0; j < N + 1; j++) {
      for (k = 0; k < D; k++) {
        t_in_sum[i][j].poly[k] =
            con_sub(a_in[i][0].cn[j].poly[k] + a_in[i][1].cn[j].poly[k], Q);
      }
    }
  }

  /* \hat{pk}_{i, j} = pk_{i, j} - G(sn_i) */
  sample_g(g_gamma, ask[0].sn);
  for (i = 0; i < N_SPENT; i++) {
    for (j = 0; j < N_BAR; j++) {
      for (k = 0; k < D; k++) {
        pk_hat[i][0][j].poly[k] =
            con_add(a_in[i][0].pk[j].poly[k] - g_gamma[j].poly[k], Q);
      }
    }
  }
  sample_g(g_gamma, ask[1].sn);
  for (i = 0; i < N_SPENT; i++) {
    for (j = 0; j < N_BAR; j++) {
      for (k = 0; k < D; k++) {
        pk_hat[i][1][j].poly[k] =
            con_add(a_in[i][1].pk[j].poly[k] - g_gamma[j].poly[k], Q);
      }
    }
  }

  /* b = (\delta_{l, 0},...,\delta_{l, N_SPENT - 1}, b_0, b_1) s.t. b_0 - b_1 =
   * c_32 */
  for (i = 0; i < N_SPENT; i++) {
    delta[i] = ct_eq(i, l);
  }
  delta[N_SPENT] = (-c[32]) >> 63;
  delta[N_SPENT + 1] = c[32] >> 63;

  /* c = [c_0,...,c_{r' - 1}, 0, c_{r' + 1},...,c_{R - 1}] */
  c[32] = 0;

  do {
    rej = 0;

    /* a_1,...,a_{N_SPENT + 1} <-- D_{10 * B_a} */
    sample_ba_n1(a + 1);
    for (i = 1; i < N_SPENT + 2; i++) {
      for (j = 0; j < D; j++) {
        a_ntt[i].poly[j] = con_add(a[i].poly[j], Q);
        a_ntt_qhat[i].poly[j] = con_add(a[i].poly[j], QHAT);
      }

      ntt_q(a_ntt + i);
      ntt_qhat(a_ntt_qhat + i);
    }

    /* a_0 = -\sum_{j = 1}^{N_SPENT - 1} a_j */
    memset(a, 0, sizeof(POLY_R));
    for (i = 1; i < N_SPENT; i++) {
      for (j = 0; j < D; j++) {
        a[0].poly[j] -= a[i].poly[j];
      }
    }
    for (j = 0; j < D; j++) {
      a_ntt[0].poly[j] = con_add(a[0].poly[j], Q);
      a_ntt_qhat[0].poly[j] = con_add(a[0].poly[j], QHAT);
    }
    ntt_q(a_ntt);
    ntt_qhat(a_ntt_qhat);

    /* y_{out, i} <-- D_{B}^{N + K + 1} */
    /* y_0 <-- D_B^{N + K + 1} */
    
    int64_t *samps = (int64_t*)calloc(NN, sizeof(int64_t));
    uint32_t ct = 0;

#ifdef MAGNET_    
    MAGNET(samps);
#endif

    sample_b_nk1(y_out[0], samps, &ct);
    sample_b_nk1(y_out[1], samps, &ct);
    sample_b_nk1(y0, samps, &ct);
    for (i = 0; i < N + K + 1; i++) {
      for (j = 0; j < D; j++) {
        y_out_ntt[0][i].poly[j] = con_add(y_out[0][i].poly[j], Q);
        y_out_ntt[1][i].poly[j] = con_add(y_out[1][i].poly[j], Q);
        y0_ntt[i].poly[j] = con_add(y0[i].poly[j], Q);
      }

      ntt_q(y_out_ntt[0] + i);
      ntt_q(y_out_ntt[1] + i);
      ntt_q(y0_ntt + i);
    }

    /* w_{out, i} <-- B * y_{out, i} */
    /* B * y_0 */
    for (i = 0; i < N; i++) {
      memcpy(w_out[0] + i, y_out_ntt[0] + i, sizeof(POLY_Q));
      memcpy(w_out[1] + i, y_out_ntt[1] + i, sizeof(POLY_Q));
      memcpy(w0 + i, y0_ntt + i, sizeof(POLY_Q));
    }
    for (i = 0; i < K + 1; i++) {
      for (j = 0; j < N; j++) {
        mult_plus_rq(w_out[0] + j, mat1[i] + j, y_out_ntt[0] + N + i);
        mult_plus_rq(w_out[1] + j, mat1[i] + j, y_out_ntt[1] + N + i);
        mult_plus_rq(w0 + j, mat1[i] + j, y0_ntt + N + i);
      }
    }

    /* w_0 = B * y_0 - \sum_{j = 0}^{N_SPENT - 1} a_j * \sum_{i = 0}^{M - 1}
     * t_{in, j}^(i) */
    for (i = 0; i < N_SPENT; i++) {
      for (j = 0; j < N; j++) {
        mult_minus_rq(w0 + j, a_ntt + i, t_in_sum[i] + j);
      }
    }

    /* rho <-- S_1^{N + K + 3} */
    /* y_1, y_2 <-- D_{B}^{N + K + 3} */
    sample_1(rho, N + K + 3);
    
    sample_b_nk3(y1, samps, &ct);
    sample_b_nk3(y2, samps, &ct);
    for (i = 0; i < N + K + 3; i++) {
      for (j = 0; j < D; j++) {
        rho_ntt[i].poly[j] = con_add(rho[i].poly[j], Q);
        y1_ntt[i].poly[j] = con_add(y1[i].poly[j], Q);
        y2_ntt[i].poly[j] = con_add(y2[i].poly[j], Q);
      }

      ntt_q(rho_ntt + i);
      ntt_q(y1_ntt + i);
      ntt_q(y2_ntt + i);
    }

    /* v = A * \rho */
    /* w_i = A * y_i */
    for (i = 0; i < N; i++) {
      memcpy((out->v) + i, rho_ntt + i, sizeof(POLY_Q));
      memcpy(w1 + i, y1_ntt + i, sizeof(POLY_Q));
      memcpy(w2 + i, y2_ntt + i, sizeof(POLY_Q));
    }
    for (i = 0; i < K + 3; i++) {
      for (j = 0; j < N; j++) {
        mult_plus_rq((out->v) + j, mat1[i] + j, rho_ntt + N + i);
        mult_plus_rq(w1 + j, mat1[i] + j, y1_ntt + N + i);
        mult_plus_rq(w2 + j, mat1[i] + j, y2_ntt + N + i);
      }
    }

    /* C = <a_3', \rho> + c */
    memcpy(&(out->c), rho_ntt + N + 2, sizeof(POLY_Q));
    for (i = 0; i < SUBRING_SIZE; i++) {
      (out->c).poly[i] = 0;
      (out->c).poly[i + (D >> 1)] = 0;
    }
    for (i = 0; i < K; i++) {
      mult_plus_rq(&(out->c), mat2[i] + 2, rho_ntt + N + 3 + i);
    }
    for (i = 0; i < R; i++) {
      (out->c).poly[i * SUBRING_SIZE] = red_short_q(
          (out->c).poly[i * SUBRING_SIZE] + (c[i] * MONTGOMERY_32_Q) + Q);
    }

    /* w = <a_3', y_1> - 2\sigma(<a_3', y_2>) + <b_1, \sum_{i = 0}^{S - 1}
     * y_{out_i} - y_0> + \sum_{j = 0}^{N_SPENT - 1} (a_j * \sum_{i = 0}^{M - 1}
     * v_{in, j}^(i) - (a_{N_SPENT} - a_{N_SPENT + 1}) * [0^{r' - 1}, -2, 1,
     * 0^{r' - 1}] */
    /* <a_3', y_1> */
    memcpy(&a3y1, y1_ntt + N + 2, sizeof(POLY_Q));
    for (i = 0; i < SUBRING_SIZE; i++) {
      a3y1.poly[i] = 0;
      a3y1.poly[i + (D >> 1)] = 0;
    }
    for (i = 0; i < K; i++) {
      mult_plus_rq(&a3y1, mat2[i] + 2, y1_ntt + N + 3 + i);
    }

    /* <b_1, y_{out, i}> */
    memcpy(b1yout, y_out_ntt[0] + N, sizeof(POLY_Q));
    memcpy(b1yout + 1, y_out_ntt[1] + N, sizeof(POLY_Q));
    for (i = 0; i < K - 2; i++) {
      mult_plus_rq(b1yout, mat2[i], y_out_ntt[0] + N + 3 + i);
      mult_plus_rq(b1yout + 1, mat2[i], y_out_ntt[1] + N + 3 + i);
    }

    /* <a_3', y_2> */
    memcpy(&tmp_q, y2_ntt + N + 2, sizeof(POLY_Q));
    for (i = 0; i < SUBRING_SIZE; i++) {
      tmp_q.poly[i] = 0;
      tmp_q.poly[i + (D >> 1)] = 0;
    }
    for (i = 0; i < K; i++) {
      mult_plus_rq(&tmp_q, mat2[i] + 2, y2_ntt + N + 3 + i);
    }

    /* \sigma(<a_3', y_2>) */
    aut(&tmp_aut, &tmp_q);

    /* <a_3', y_1> - 2\sigma(<a_3', y_2>) + <b_1, \sum_{i = 0}^{S - 1} y_{out_i}
     * - y_0> */
    for (i = 0; i < D; i++) {
      w.poly[i] = red_short_q(con_add(a3y1.poly[i] - 2 * tmp_aut.poly[i] +
                                          b1yout[0].poly[i] +
                                          b1yout[1].poly[i] - y0_ntt[N].poly[i],
                                      Q * 3));
    }
    for (i = 0; i < K - 2; i++) {
      mult_minus_rq(&w, mat2[i], y0_ntt + N + 3 + i);
    }

    /* <a_3', y_1> - 2\sigma(<a_3', y_2>) + <b_1, \sum_{i = 0}^{S - 1} y_{out_i}
     * - y_0> + \sum_{j = 0}^{N_SPENT - 1} (a_j * \sum_{i = 0}^{M - 1} v_{in,
     * j}^(i) */
    for (i = 0; i < N_SPENT; i++) {
      mult_plus_rq(&w, a_ntt + i, t_in_sum[i] + N);
    }

    /* <a_3', y_1> - 2\sigma(<a_3', y_2>) + <b_1, \sum_{i = 0}^{S - 1} y_{out_i}
     * - y_0> + \sum_{j = 0}^{N_SPENT - 1} (a_j * \sum_{i = 0}^{M - 1} v_{in,
     * j}^(i) - (a_{N_SPENT} - a_{N_SPENT + 1}) * [0^{r' - 1}, -2, 1, 0^{r' -
     * 1}] */
    for (i = 0; i < SUBRING_SIZE; i++) {
      w.poly[31 * SUBRING_SIZE + i] = red_short_q(
          con_add(w.poly[31 * SUBRING_SIZE + i] +
                      2 * (a_ntt[N_SPENT].poly[31 * SUBRING_SIZE + i] -
                           a_ntt[N_SPENT + 1].poly[31 * SUBRING_SIZE + i]),
                  Q << 1));
      w.poly[32 * SUBRING_SIZE + i] =
          red_short_q(w.poly[32 * SUBRING_SIZE + i] + Q -
                      a_ntt[N_SPENT].poly[32 * SUBRING_SIZE + i] +
                      a_ntt[N_SPENT + 1].poly[32 * SUBRING_SIZE + i]);
    }

    /* v_2 = <a_2, \rho> + g_2,
     * g_2 = (3c^2 - 1) * <a_3', y_1> */
    /* <a_2, \rho> */
    memcpy(&(out->v2), rho_ntt + N + 1, sizeof(POLY_Q));
    for (i = 0; i < K; i++) {
      mult_plus_rq(&(out->v2), mat2[i] + 1, rho_ntt + N + 3 + i);
    }

    /* <a_2, \rho> + (3c^2 - 1) * <a_3', y_1> */
    for (i = 0; i < R; i++) {
      for (j = 0; j < SUBRING_SIZE; j++) {
        (out->v2).poly[i * SUBRING_SIZE + j] = red_short_q(
            (out->v2).poly[i * SUBRING_SIZE + j] + Q +
            (3 * c[i] * c[i] - 1) * a3y1.poly[i * SUBRING_SIZE + j]);
      }
    }

    /* \alpha <-- H'(IN, COM) */
    hash_alpha(out->alpha, pp, a_in, sn_out, pk_out, cn_out, w_out, w0, out->v,
               &(out->v2), &(out->c), &w, w1, w2);

    /* \alpha_0,...,\alpha_{S - 1}, \alpha_0',...,\alpha_{M - 2}' <--
     * Expand(\alpha, "Ch") */
    fastrandombytes_setseed_ch(out->alpha);
    sample_alpha(alpha);
    sample_alpha_prime(&alpha_prime_0);
    for (i = 0; i < D; i++) {
      alpha_prime_0_ntt.poly[i] = con_add(alpha_prime_0.poly[i], Q);
    }
    ntt_q(&alpha_prime_0_ntt);

    /* check invertibility of \alpha_0' */
    for (i = 0; i < R; i++) {
      rej1 = 1;
      for (j = 0; j < SUBRING_SIZE; j++) {
        rej1 &= ct_eq(alpha_prime_0_ntt.poly[i * SUBRING_SIZE + j], 0);
      }
      rej |= rej1;
    }
    if (__builtin_expect((rej), 0)) {
      continue;
    }

    /* v_0 = <a_1, y_1> + \hat{g}_0,
     * \hat{g}_0 = \sum_{i = 0}^{S - 1} (\alpha_i * g_0^(i)) + g_0^(S)
     * g_0^(i) = <b_1, y_{out, i}>^2
     * g_0^(S) = <a_3', y_1>^3 */
    /* \alpha_i * <b_1, y_{out, i}> */
    mult_rq(alpha_b1yout, alpha, b1yout);
    mult_rq(alpha_b1yout + 1, alpha + 1, b1yout + 1);

    /* <a_1, y_1> */
    memcpy(&v0, y1_ntt + N, sizeof(POLY_Q));
    for (i = 0; i < K; i++) {
      mult_plus_rq(&v0, mat2[i], y1_ntt + N + 3 + i);
    }

    /* <a_1, y_1> + <a_3', y_1>^3 */
    mult_rq(&a3y1_2, &a3y1, &a3y1);
    mult_plus_rq(&v0, &a3y1_2, &a3y1);

    /* <a_1, y_1> + \sum_{i = 0}^{S - 1} (\alpha_i * <b_1, y_{out, i}>^2) +
     * <a_3', y_1>^3 */
    mult_plus_rq(&v0, alpha_b1yout, b1yout);
    mult_plus_rq(&v0, alpha_b1yout + 1, b1yout + 1);

    /* v_1 = <a_1, \rho> + <a_2, y_1> + \hat{g}_1,
     * \hat{g}_1 = \sum_{i = 0}^{S - 1} (\alpha_i * g_1^(i)) + g_1^(S)
     * g_1^(i) = (-2 \hat{a}_i + 1) * <b_1, y_{out, i}>
     * g_1^(S) = (-3c) * <a_3', y_1>^2 */
    /* <a_1, \rho> + <a_2, y_1> */
    for (i = 0; i < D; i++) {
      (out->v1).poly[i] =
          con_sub(rho_ntt[N].poly[i] + y1_ntt[N + 1].poly[i], Q);
    }
    for (i = 0; i < K; i++) {
      mult_plus_rq(&(out->v1), mat2[i], rho_ntt + N + 3 + i);
      mult_plus_rq(&(out->v1), mat2[i] + 1, y1_ntt + N + 3 + i);
    }

    /* <a_1, \rho> + <a_2, y_1> + \sum_{i = 0}^{S - 1} (\alpha_i * (-2 \hat{a}_i
     * + 1) * <b_1, y_{out, i}>) + (-3c) * <a_3', y_1>^2 */
    for (i = 0; i < R; i++) {
      for (j = 0; j < SUBRING_SIZE; j++) {
        (out->v1).poly[i * SUBRING_SIZE + j] =
            red_short_q(red_short_q(con_add(
                            (out->v1).poly[i * SUBRING_SIZE + j] +
                                (1 - 2 * a_hat[0][i]) *
                                    alpha_b1yout[0].poly[i * SUBRING_SIZE + j] +
                                (1 - 2 * a_hat[1][i]) *
                                    alpha_b1yout[1].poly[i * SUBRING_SIZE + j],
                            Q << 1)) +
                        3 * (Q - c[i] * a3y1_2.poly[i * SUBRING_SIZE + j]));
      }
    }

    /* P_j = \sum_{i = 0}^{M - 2} (a_i' * \hat{pk}_{i, j}) + \hat{pk}_{M - 1, j}
     */
    for (i = 0; i < N_SPENT; i++) {
      for (j = 0; j < N_BAR; j++) {
        mult_rq(p[i] + j, &alpha_prime_0_ntt, pk_hat[i][0] + j);

        for (k = 0; k < D; k++) {
          p[i][j].poly[k] =
              con_sub(p[i][j].poly[k] + pk_hat[i][1][j].poly[k], Q);
        }
      }
    }

    /* s = \sum_{i = 0}^{M - 2} (a_i' * sk_i) + sk_{M - 1} */
    for (i = 0; i < N_BAR + K; i++) {
      mult_r(s + i, &alpha_prime_0, (ask[0].sk) + i);

      for (j = 0; j < D; j++) {
        s[i].poly[j] += ask[1].sk[i].poly[j];
      }
    }

    /* (a_0^2,...,a_{N_SPENT + 1}^2) */
    for (i = 0; i < N_SPENT + 2; i++) {
      mult_rqhat(a2_ntt + i, a_ntt_qhat + i, a_ntt_qhat + i);
    }

    /* rb <-- S_1^{\hat{n} + \hat{k}} */
    /* ra <-- D_B^{\hat{n} + \hat{k}} */
    sample_1(rb, N_HAT + K_HAT);
    sample_b_nhatkhat(ra, samps, &ct);
    for (i = 0; i < N_HAT + K_HAT; i++) {
      for (j = 0; j < D; j++) {
        rb_ntt[i].poly[j] = con_add(rb[i].poly[j], QHAT);
        ra_ntt[i].poly[j] = con_add(ra[i].poly[j], QHAT);
      }

      ntt_qhat(rb_ntt + i);
      ntt_qhat(ra_ntt + i);
    }

    /* B = G' * (r_b, b, c) */
    /* A = G' * (r_a, a, d) */
    for (i = 0; i < N_HAT; i++) {
      memcpy((out->b) + i, rb_ntt + i, sizeof(POLY_QHAT));
      memcpy(a_rs + i, ra_ntt + i, sizeof(POLY_QHAT));
    }
    for (i = 0; i < K_HAT; i++) {
      for (j = 0; j < N_HAT; j++) {
        mult_plus_rqhat((out->b) + j, mat3[i] + j, rb_ntt + N_HAT + i);
        mult_plus_rqhat(a_rs + j, mat3[i] + j, ra_ntt + N_HAT + i);
      }
    }
    for (i = 0; i < N_SPENT + 2; i++) {
      for (j = 0; j < N_HAT; j++) {
        for (k = 0; k < D; k++) {
          (out->b)[j].poly[k] += (-delta[i]) & mat3[K + i][j].poly[k];
        }

        mult_plus_rqhat(a_rs + j, mat3[K + i] + j, a_ntt_qhat + i);
      }
    }
    for (i = 0; i < N_SPENT + 2; i++) {
      for (j = 0; j < N_HAT; j++) {
        mult_plus_rqhat_pm((out->b) + j, mat3[K + N_SPENT + 2 + i] + j,
                           a_ntt_qhat + i, delta[i]);
        mult_minus_rqhat(a_rs + j, mat3[K + N_SPENT + 2 + i] + j, a2_ntt + i);
      }
    }

    /* y <-- D_{13 * B_{rs}}^{\bar{n} + k} */
    sample_brs_13_nbark(y);
    for (i = 0; i < N_BAR + K; i++) {
      for (j = 0; j < D; j++) {
        y_ntt[i].poly[j] = con_add(y[i].poly[j], Q);
      }

      ntt_q(y_ntt + i);
    }

    /* E = G * y - \sum_{j = 0}^{N_SPENT - 1} a_j * P_j */
    for (i = 0; i < N_BAR; i++) {
      memcpy(e + i, y_ntt + i, sizeof(POLY_Q));
    }
    for (i = 0; i < K; i++) {
      for (j = 0; j < N_BAR; j++) {
        mult_plus_rq(e + j, mat1[i] + j + (N - N_BAR), y_ntt + N_BAR + i);
      }
    }
    for (i = 0; i < N_SPENT; i++) {
      for (j = 0; j < N_BAR; j++) {
        mult_minus_rq(e + j, a_ntt + i, p[i] + j);
      }
    }

    /* x <-- H(\alpha, v_0, v_1, A, B, E) */
    hash_x(&(out->x), out->alpha, &v0, &(out->v1), a_rs, out->b, e);

    /* f_j = x * b_j + a_j */
    for (i = 1; i < N_SPENT + 2; i++) {
      for (j = 0; j < D; j++) {
        x_delta[i - 1].poly[j] = (-delta[i]) & (out->x).poly[j];
        (out->f1)[i - 1].poly[j] = x_delta[i - 1].poly[j] + a[i].poly[j];
      }
    }
    /* Rej(f_1, x * delta, \phi_a, B_a) */
    rej = rej_f(out->f1, x_delta);
    if (rej) {
      continue;
    }

    /* f_0 */
    for (i = 0; i < D; i++) {
      f0.poly[i] = ((-delta[0]) & (out->x).poly[i]) + a[0].poly[i];
    }

    /* g_0 = f_0(x - f_0) */
    mult_fxf(g, &f0, &(out->x));

    /* check g_0 norm here */
    for (j = 0; j < D; j++) {
      rej |= ct_lt(BG_0, g[0].poly[j]) | ct_lt(g[0].poly[j], -BG_0);
    }
    if (rej) {
      continue;
    }

    /* g_1 = (f_1(x - f_1),...,f_{N_SPENT + 1}(X - f_{N_SPENT + 1})) */
    for (i = 1; i < N_SPENT + 2; i++) {
      mult_fxf(g + i, (out->f1) + i - 1, &(out->x));
    }

    /* check g_1 norm here */
    for (i = 1; i < N_SPENT + 2; i++) {
      for (j = 0; j < D; j++) {
        rej |= ct_lt(BG_1, g[i].poly[j]) | ct_lt(g[i].poly[j], -BG_1);
      }
    }
    if (rej) {
      continue;
    }

    /* z = y + x * s */
    for (i = 0; i < N_BAR + K; i++) {
      mult_r(x_s + i, &(out->x), s + i);

      for (j = 0; j < D; j++) {
        (out->z)[i].poly[j] = x_s[i].poly[j] + y[i].poly[j];
      }
    }
    /* Rej(z, x * s, \phi_rs, B_rs) */
    rej = rej_z(out->z, x_s);
    if (rej) {
      continue;
    }

    /* z_b = r_a + x * r_b */
    for (i = 0; i < N_HAT + K_HAT; i++) {
      mult_r(x_rb + i, &(out->x), rb + i);

      for (j = 0; j < D; j++) {
        (out->zb)[i].poly[j] = x_rb[i].poly[j] + ra[i].poly[j];
      }
    }

    /* z_{out, i} = y_{out, i} + x * ck_{out, i} */
    for (i = 0; i < N + K + 1; i++) {
      mult_r(x_rout[0] + i, &(out->x), (ck_out->ck) + i);

      for (j = 0; j < D; j++) {
        (out->z_out)[0][i].poly[j] = x_rout[0][i].poly[j] + y_out[0][i].poly[j];
      }
    }
    for (i = 0; i < N + K + 1; i++) {
      mult_r(x_rout[1] + i, &(out->x), ((ck_out + 1)->ck) + i);

      for (j = 0; j < D; j++) {
        (out->z_out)[1][i].poly[j] = x_rout[1][i].poly[j] + y_out[1][i].poly[j];
      }
    }

    /* z_0 = y_0 + x * \sum_{i = 0}^{M - 1} ck_i */
    for (i = 0; i < N + K + 1; i++) {
      mult_r(x_ri + i, &(out->x), (ask->ck) + i);
      mult_r(&tmp_r, &(out->x), ((ask + 1)->ck) + i);

      for (j = 0; j < D; j++) {
        x_ri[i].poly[j] += tmp_r.poly[j];

        (out->z0)[i].poly[j] = x_ri[i].poly[j] + y0[i].poly[j];
      }
    }

    /* z_1 = y_1 + x * \rho */
    for (i = 0; i < N + K + 3; i++) {
      mult_r(x_rho + i, &(out->x), rho + i);

      for (j = 0; j < D; j++) {
        (out->z1)[i].poly[j] = x_rho[i].poly[j] + y1[i].poly[j];
      }
    }

    /* sigma^{-1}(x) */
    iaut_r(&tmp_aut_r, &(out->x));

    /* z_2 = y_2 + \sigma^{-1}(x) * \rho */
    for (i = 0; i < N + K + 3; i++) {
      mult_r(x_aut_rho + i, &tmp_aut_r, rho + i);

      for (j = 0; j < D; j++) {
        (out->z2)[i].poly[j] = x_aut_rho[i].poly[j] + y2[i].poly[j];
      }
    }

    /* RejOp((z_{out, i}, z_1, z_2, z_0, z_b), (x * ck_{out, i}, x * \rho,
     * \sigma^{-1}(x) * \rho, x * \sum_{i = 0}^{M - 1} ck_i, x * r_b), 1, B) */
    rej = rej_op(out->z_out, out->z1, out->z2, out->z0, out->zb, x_rout, x_rho,
                 x_aut_rho, x_ri, x_rb);
    
    //printf("rej_op: %d\n", rej);
    free(samps);            
  } while (rej);
}

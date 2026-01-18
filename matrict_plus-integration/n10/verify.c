#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "fastrandombytes.h"
#include "hash.h"
#include "param.h"
#include "poly_mult.h"
#include "poly_q.h"
#include "poly_red.h"
#include "sampleb.h"
#include "spend.h"
#include "verify.h"

#define Z_NORM 3161803581581766LL
#define Z_OUT_NORM 135519727172LL
#define Z_0_NORM 135519727172LL
#define Z_I_NORM 165635222099LL

uint64_t verify(const ACT a_in[][M], const PK_OUT *pk_out, const CN_OUT *cn_out,
                const SN_OUT *sn_out, const SPEND_OUT *in,
                const POLY_Q mat1[][N], const POLY_Q mat2[][3],
                const POLY_QHAT mat3[][N_HAT], unsigned char *pp) {
  uint64_t i, j, k;
  uint64_t tmp;

  static POLY_R f0;
  static POLY_R g[N_SPENT + 2];
  static POLY_Q c;
  static POLY_Q x_ntt;
  static POLY_QHAT x_ntt_qhat;
  static POLY_Q w_out[S][N];
  static POLY_Q z_out_ntt[S][N + K + 1];
  static POLY_Q t_in_sum[N_SPENT][N + 1];
  static POLY_Q z0_ntt[N + K + 1];
  static POLY_Q f_ntt[N_SPENT + 2];
  static POLY_QHAT f_ntt_qhat[N_SPENT + 2];
  static POLY_Q w0[N];
  static POLY_Q z1_ntt[N + K + 3];
  static POLY_Q w1[N];
  static POLY_Q tmp_aut;
  static POLY_Q z2_ntt[N + K + 3];
  static POLY_Q w2[N];
  static POLY_Q a3z1;
  static POLY_Q tmp_q;
  static POLY_Q w;
  static unsigned char alpha_verify[CRYPTO_BYTES];
  static POLY_Q alpha[S];
  static POLY_R alpha_prime_0;
  static POLY_Q h[S + 1];
  static POLY_Q v0;
  static POLY_Q hx[S + 1];
  static POLY_QHAT a[N_HAT];
  static POLY_QHAT zb_ntt[N_HAT + K_HAT];
  static POLY_QHAT g_ntt[N_SPENT + 2];
  static POLY_Q pk_hat[N_SPENT][M][N_BAR];
  static POLY_Q g_gamma[N_BAR];
  static POLY_Q p[N_SPENT][N_BAR];
  static POLY_Q e[N_BAR];
  static POLY_Q z_ntt[N_BAR + K];
  static POLY_Q alpha_prime_0_ntt;
  static POLY_R x;
  static POLY_R tmp_aut_r;

  /* check z_{out, i} norm */
  tmp = 0;
  for (i = 0; i < N + K + 1; i++) {
    for (j = 0; j < D; j++) {
      tmp += (int64_t)((in->z_out)[0][i].poly[j]) *
             (int64_t)((in->z_out)[0][i].poly[j]);
      if (((int64_t)((in->z_out)[0][i].poly[j]) > B_6) ||
          ((int64_t)((in->z_out)[0][i].poly[j]) < -B_6)) {
        return 0;
      }
    }
  }
  if (tmp > Z_OUT_NORM) {
    return 0;
  }
  tmp = 0;
  for (i = 0; i < N + K + 1; i++) {
    for (j = 0; j < D; j++) {
      tmp += (int64_t)((in->z_out)[1][i].poly[j]) *
             (int64_t)((in->z_out)[1][i].poly[j]);
      if (((int64_t)((in->z_out)[1][i].poly[j]) > B_6) ||
          ((int64_t)((in->z_out)[1][i].poly[j]) < -B_6)) {
        return 0;
      }
    }
  }
  if (tmp > Z_OUT_NORM) {
    return 0;
  }

  /* check z_0 norm */
  tmp = 0;
  for (i = 0; i < N + K + 1; i++) {
    for (j = 0; j < D; j++) {
      tmp += (int64_t)((in->z0)[i].poly[j]) * (int64_t)((in->z0)[i].poly[j]);
      if (((int64_t)((in->z0)[i].poly[j]) > B_6) ||
          ((int64_t)((in->z0)[i].poly[j]) < -B_6)) {
        return 0;
      }
    }
  }
  if (tmp > Z_0_NORM) {
    return 0;
  }

  /* check z_i norm */
  tmp = 0;
  for (i = 0; i < N + K + 3; i++) {
    for (j = 0; j < D; j++) {
      tmp += (int64_t)((in->z1)[i].poly[j]) * (int64_t)((in->z1)[i].poly[j]);
      if (((int64_t)((in->z1)[i].poly[j]) > B_6) ||
          ((int64_t)((in->z1)[i].poly[j]) < -B_6)) {
        return 0;
      }
    }
  }
  if (tmp > Z_I_NORM) {
    return 0;
  }
  tmp = 0;
  for (i = 0; i < N + K + 3; i++) {
    for (j = 0; j < D; j++) {
      tmp += (int64_t)((in->z2)[i].poly[j]) * (int64_t)((in->z2)[i].poly[j]);
      if (((int64_t)((in->z2)[i].poly[j]) > B_6) ||
          ((int64_t)((in->z2)[i].poly[j]) < -B_6)) {
        return 0;
      }
    }
  }
  if (tmp > Z_I_NORM) {
    return 0;
  }

  /* check z_b norm */
  for (i = 0; i < N_HAT + K_HAT; i++) {
    for (j = 0; j < D; j++) {
      if (((int64_t)((in->zb)[i].poly[j]) > B_6) ||
          ((int64_t)((in->zb)[i].poly[j]) < -B_6)) {
        return 0;
      }
    }
  }

  /* check f_1 norm */
  for (i = 0; i < N_SPENT + 1; i++) {
    for (j = 0; j < D; j++) {
      if (((int64_t)((in->f1)[i].poly[j]) > BA_10_6) ||
          ((int64_t)((in->f1)[i].poly[j]) < -BA_10_6)) {
        return 0;
      }
    }
  }

  /* check z norm */
  tmp = 0;
  for (i = 0; i < N_BAR + K; i++) {
    for (j = 0; j < D; j++) {
      tmp += (int64_t)((in->z)[i].poly[j]) * (int64_t)((in->z)[i].poly[j]);
      if (((int64_t)((in->z)[i].poly[j]) > BRS_13_6) ||
          ((int64_t)((in->z)[i].poly[j]) < -BRS_13_6)) {
        return 0;
      }
    }
  }
  if (tmp > Z_NORM) {
    return 0;
  }

  /* f_0 = x - \sum_{j = 0}^{N_SPENT - 1} f_j */
  memcpy(&f0, &(in->x), sizeof(POLY_Q));
  for (i = 0; i < N_SPENT - 1; i++) {
    for (j = 0; j < D; j++) {
      f0.poly[j] -= (in->f1)[i].poly[j];
    }
  }

  /* g_0 = f_0(x - f_0) */
  mult_fxf(g, &f0, &(in->x));

  /* check g_0 norm */
  for (j = 0; j < D; j++) {
    if ((((int64_t)g[0].poly[j]) > BG_0) || (((int64_t)g[0].poly[j]) < -BG_0)) {
      return 0;
    }
  }

  /* g_1 = (f_1(x - f_1),...,f_{N_SPENT + 1}(X - f_{N_SPENT + 1})) */
  for (i = 1; i < N_SPENT + 2; i++) {
    mult_fxf(g + i, (in->f1) + i - 1, &(in->x));
  }

  /* check g_1 norm */
  for (i = 1; i < N_SPENT + 2; i++) {
    for (j = 0; j < D; j++) {
      if ((((int64_t)g[i].poly[j]) > BG_1) ||
          (((int64_t)g[i].poly[j]) < -BG_1)) {
        return 0;
      }
    }
  }

  /* C' = C * [0, 1^{r' - 1}, 0, 1^{r' - 1}] */
  memcpy(&c, &(in->c), sizeof(POLY_Q));
  for (i = 0; i < SUBRING_SIZE; i++) {
    c.poly[i] = 0;
    c.poly[(D >> 1) + i] = 0;
  }

  /* \sigma^{-1}(x) */
  iaut_r(&tmp_aut_r, &(in->x));

  /* ntt(x) */
  for (i = 0; i < D; i++) {
    x_ntt.poly[i] = con_add((in->x).poly[i], Q);
    x_ntt_qhat.poly[i] = con_add((in->x).poly[i], QHAT);

    tmp_aut.poly[i] = con_add(tmp_aut_r.poly[i], Q);
  }
  ntt_q(&x_ntt);
  ntt_qhat(&x_ntt_qhat);
  ntt_q(&tmp_aut);

  /* ntt(z_{out, i}), ntt(z_0) */
  for (i = 0; i < N + K + 1; i++) {
    for (j = 0; j < D; j++) {
      z_out_ntt[0][i].poly[j] = con_add((in->z_out)[0][i].poly[j], Q);
      z_out_ntt[1][i].poly[j] = con_add((in->z_out)[1][i].poly[j], Q);
      z0_ntt[i].poly[j] = con_add((in->z0)[i].poly[j], Q);
    }

    ntt_q(z_out_ntt[0] + i);
    ntt_q(z_out_ntt[1] + i);
    ntt_q(z0_ntt + i);
  }

  /* w_{out, i} = B * z_{out, i} - x * t^(i) */
  /* B * z_0 */
  for (i = 0; i < N; i++) {
    memcpy(w_out[0] + i, z_out_ntt[0] + i, sizeof(POLY_Q));
    memcpy(w_out[1] + i, z_out_ntt[1] + i, sizeof(POLY_Q));
    memcpy(w0 + i, z0_ntt + i, sizeof(POLY_Q));
  }
  for (i = 0; i < K + 1; i++) {
    for (j = 0; j < N; j++) {
      mult_plus_rq(w_out[0] + j, mat1[i] + j, z_out_ntt[0] + N + i);
      mult_plus_rq(w_out[1] + j, mat1[i] + j, z_out_ntt[1] + N + i);
      mult_plus_rq(w0 + j, mat1[i] + j, z0_ntt + N + i);
    }
  }
  for (i = 0; i < N; i++) {
    mult_minus_rq(w_out[0] + i, &x_ntt, (cn_out->cn) + i);
    mult_minus_rq(w_out[1] + i, &x_ntt, ((cn_out + 1)->cn) + i);
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

  /* ntt(f) */
  for (i = 0; i < D; i++) {
    f_ntt[0].poly[i] = con_add(f0.poly[i], Q);
    f_ntt_qhat[0].poly[i] = con_add(f0.poly[i], QHAT);
  }
  ntt_q(f_ntt);
  ntt_qhat(f_ntt_qhat);
  for (i = 1; i < N_SPENT + 2; i++) {
    for (j = 0; j < D; j++) {
      f_ntt[i].poly[j] = con_add((in->f1)[i - 1].poly[j], Q);
      f_ntt_qhat[i].poly[j] = con_add((in->f1)[i - 1].poly[j], QHAT);
    }

    ntt_q(f_ntt + i);
    ntt_qhat(f_ntt_qhat + i);
  }

  /* w_0 = B * z_0 - \sum_{j = 0}^{N_SPENT - 1} (f_j * \sum_{i = 0}^{M - 1}
   * t_{in, j}^(i) */
  for (i = 0; i < N_SPENT; i++) {
    for (j = 0; j < N; j++) {
      mult_minus_rq(w0 + j, f_ntt + i, t_in_sum[i] + j);
    }
  }

  /* ntt(z_1), ntt(z_2) */
  for (i = 0; i < N + K + 3; i++) {
    for (j = 0; j < D; j++) {
      z1_ntt[i].poly[j] = con_add((in->z1)[i].poly[j], Q);
      z2_ntt[i].poly[j] = con_add((in->z2)[i].poly[j], Q);
    }

    ntt_q(z1_ntt + i);
    ntt_q(z2_ntt + i);
  }

  /* w_1 = A * z_1 - x * v */
  /* w_2 = A * z_2 - \sigma^{-1}(x) * v */
  for (i = 0; i < N; i++) {
    memcpy(w1 + i, z1_ntt + i, sizeof(POLY_Q));
    memcpy(w2 + i, z2_ntt + i, sizeof(POLY_Q));
  }
  for (i = 0; i < K + 3; i++) {
    for (j = 0; j < N; j++) {
      mult_plus_rq(w1 + j, mat1[i] + j, z1_ntt + N + i);
      mult_plus_rq(w2 + j, mat1[i] + j, z2_ntt + N + i);
    }
  }
  for (i = 0; i < N; i++) {
    mult_minus_rq(w1 + i, &x_ntt, (in->v) + i);
    mult_minus_rq(w2 + i, &tmp_aut, (in->v) + i);
  }

  /* w = <a_3', z_1> - 2\sigma(<a_3', z_2>) + <b_1, \sum_{i = 0}^{S - 1}
   * z_{out_i} - z_0> - x * \sum_{i = 0}^{S - 1} t_1^(i) + \sum_{j = 0}^{N_SPENT
   * - 1} (f_j * \sum_{i = 0}^{M - 1} v_{in, j}^(i) - xC' - (f_{N_SPENT} -
   * f_{N_SPENT + 1}) * [0^{r' - 1}, -2, 1, 0^{r' - 1}] + x * 2\sigma(C') */
  /* <a_3', z_1> */
  memcpy(&a3z1, z1_ntt + N + 2, sizeof(POLY_Q));
  for (i = 0; i < SUBRING_SIZE; i++) {
    a3z1.poly[i] = 0;
    a3z1.poly[i + (D >> 1)] = 0;
  }
  for (i = 0; i < K; i++) {
    mult_plus_rq(&a3z1, mat2[i] + 2, z1_ntt + N + 3 + i);
  }

  /* <a_3', z_2> */
  memcpy(&tmp_q, z2_ntt + N + 2, sizeof(POLY_Q));
  for (i = 0; i < SUBRING_SIZE; i++) {
    tmp_q.poly[i] = 0;
    tmp_q.poly[i + (D >> 1)] = 0;
  }
  for (i = 0; i < K; i++) {
    mult_plus_rq(&tmp_q, mat2[i] + 2, z2_ntt + N + 3 + i);
  }

  /* \sigma(<a_3', z_2>) */
  aut(&tmp_aut, &tmp_q);

  /* <a_3', z_1> - 2\sigma(<a_3', z_2>) + <b_1, \sum_{i = 0}^{S - 1} z_{out_i} -
   * z_0> */
  for (i = 0; i < D; i++) {
    w.poly[i] = red_short_q(
        con_add(a3z1.poly[i] - 2 * tmp_aut.poly[i] + z_out_ntt[0][N].poly[i] +
                    z_out_ntt[1][N].poly[i] - z0_ntt[N].poly[i],
                Q * 3));
  }
  for (i = 0; i < K - 2; i++) {
    for (j = 0; j < D; j++) {
      tmp_q.poly[j] = red_short_q(z_out_ntt[0][N + 3 + i].poly[j] +
                                  z_out_ntt[1][N + 3 + i].poly[j] + Q -
                                  z0_ntt[N + 3 + i].poly[j]);
    }

    mult_plus_rq(&w, mat2[i], &tmp_q);
  }

  /* -\sum_{i = 0}^{S - 1} t_1^(i) - C' + 2\sigma(C') */
  aut(&tmp_aut, &c);
  for (i = 0; i < D; i++) {
    tmp_q.poly[i] =
        red_short_q(con_add(2 * tmp_aut.poly[i] - (cn_out->cn)[N].poly[i] -
                                ((cn_out + 1)->cn)[N].poly[i] - c.poly[i],
                            Q * 3));
  }

  /* <a_3', z_1> - 2\sigma(<a_3', z_2>) + <b_1, \sum_{i = 0}^{S - 1} z_{out_i} -
   * z_0> - x * \sum_{i = 0}^{S - 1} t_1^(i) - xC' + x * 2\sigma(C') */
  mult_plus_rq(&w, &x_ntt, &tmp_q);

  /* <a_3', z_1> - 2\sigma(<a_3', z_2>) + <b_1, \sum_{i = 0}^{S - 1} z_{out_i} -
   * z_0> - x * \sum_{i = 0}^{S - 1} t_1^(i) - xC' + x * 2\sigma(C') + \sum_{j =
   * 0}^{N_SPENT - 1} (f_j * \sum_{i = 0}^{M - 1} v_{in, j}^(i) */
  for (i = 0; i < N_SPENT; i++) {
    mult_plus_rq(&w, f_ntt + i, t_in_sum[i] + N);
  }

  /* <a_3', z_1> - 2\sigma(<a_3', z_2>) + <b_1, \sum_{i = 0}^{S - 1} z_{out_i} -
   * z_0> - x * \sum_{i = 0}^{S - 1} t_1^(i) - xC' + x * 2\sigma(C') + \sum_{j =
   * 0}^{N_SPENT - 1} (f_j * \sum_{i = 0}^{M - 1} v_{in, j}^(i) - (f_{N_SPENT} -
   * f_{N_SPENT + 1}) * [0^{r' - 1}, -2, 1, 0^{r' - 1}] */
  for (i = 0; i < SUBRING_SIZE; i++) {
    w.poly[31 * SUBRING_SIZE + i] = red_short_q(
        con_add(w.poly[31 * SUBRING_SIZE + i] +
                    2 * (f_ntt[N_SPENT].poly[31 * SUBRING_SIZE + i] -
                         f_ntt[N_SPENT + 1].poly[31 * SUBRING_SIZE + i]),
                Q << 1));
    w.poly[32 * SUBRING_SIZE + i] =
        red_short_q(w.poly[32 * SUBRING_SIZE + i] + Q -
                    f_ntt[N_SPENT].poly[32 * SUBRING_SIZE + i] +
                    f_ntt[N_SPENT + 1].poly[32 * SUBRING_SIZE + i]);
  }

  /* \alpha <-- H'(IN, COM) */
  hash_alpha(alpha_verify, pp, a_in, sn_out, pk_out, cn_out, w_out, w0, in->v,
             &(in->v2), &(in->c), &w, w1, w2);
  for (i = 0; i < CRYPTO_BYTES; i++) {
    if (alpha_verify[i] != (in->alpha)[i]) {
      return 0;
    }
  }

  /* \alpha_0,...,\alpha_{S - 1}, \alpha_0',...,\alpha_{M - 2}' <--
   * Expand(\alpha, "Ch") */
  fastrandombytes_setseed_ch(alpha_verify);
  sample_alpha(alpha);
  sample_alpha_prime(&alpha_prime_0);

  /* h_i = <b_1, z_{out, i} - x * t_1^(i) */
  memcpy(h, z_out_ntt[0] + N, sizeof(POLY_Q));
  memcpy(h + 1, z_out_ntt[1] + N, sizeof(POLY_Q));
  for (i = 0; i < K - 2; i++) {
    mult_plus_rq(h, mat2[i], z_out_ntt[0] + N + 3 + i);
    mult_plus_rq(h + 1, mat2[i], z_out_ntt[1] + N + 3 + i);
  }
  mult_minus_rq(h, &x_ntt, (cn_out->cn) + N);
  mult_minus_rq(h + 1, &x_ntt, ((cn_out + 1)->cn) + N);

  /* h_S = <a_3', z_1> - x * C' */
  mult_rq(&tmp_q, &x_ntt, &c);
  for (i = 0; i < D; i++) {
    h[S].poly[i] = con_add(a3z1.poly[i] - tmp_q.poly[i], Q);
  }

  /* v_0 = h_S * (h_S + x) * (h_S - x) + \sum_{i = 0}^{S - 1} (\alpha_i * h_i *
   * (h_i + x)) - x * v_1 - x^2 * v_2 + <a_1, z_1> + x * <a_2, z_2> */
  /* h * (h + x) */
  for (i = 0; i < S + 1; i++) {
    for (j = 0; j < D; j++) {
      tmp_q.poly[j] = con_sub(h[i].poly[j] + x_ntt.poly[j], Q);
    }

    mult_rq(hx + i, h + i, &tmp_q);
  }

  /* h_S - x */
  for (i = 0; i < D; i++) {
    tmp_q.poly[i] = con_add(h[S].poly[i] - x_ntt.poly[i], Q);
  }

  /* h_S * (h_S + x) * (h_S - x) */
  mult_rq(&v0, hx + S, &tmp_q);

  /* h_S * (h_S + x) * (h_S - x) + \sum_{i = 0}^{S - 1} (\alpha_i * h_i * (h_i +
   * x)) */
  mult_plus_rq(&v0, alpha, hx);
  mult_plus_rq(&v0, alpha + 1, hx + 1);

  /* -v_1 - x * v_2 + <a_2, z_1> */
  memcpy(&tmp_q, z1_ntt + N + 1, sizeof(POLY_Q));
  for (i = 0; i < K; i++) {
    mult_plus_rq(&tmp_q, mat2[i] + 1, z1_ntt + N + 3 + i);
  }
  mult_minus_rq(&tmp_q, &x_ntt, &(in->v2));
  for (i = 0; i < D; i++) {
    tmp_q.poly[i] = con_add(tmp_q.poly[i] - (in->v1).poly[i], Q);
  }

  /* h_S * (h_S + x) * (h_S - x) + \sum_{i = 0}^{S - 1} (\alpha_i * h_i * (h_i +
   * x)) - x * v_1 - x^2 * v_2 + x * <a_2, z_1> */
  mult_plus_rq(&v0, &x_ntt, &tmp_q);

  /* h_S * (h_S + x) * (h_S - x) + \sum_{i = 0}^{S - 1} (\alpha_i * h_i * (h_i +
   * x)) - x * v_1 - x^2 * v_2 + x * <a_2, z_1> + <a_1, z_1> */
  for (i = 0; i < D; i++) {
    v0.poly[i] = con_sub(v0.poly[i] + z1_ntt[N].poly[i], Q);
  }
  for (i = 0; i < K; i++) {
    mult_plus_rq(&v0, mat2[i], z1_ntt + N + 3 + i);
  }

  /* ntt(z_b) */
  for (i = 0; i < N_HAT + K_HAT; i++) {
    for (j = 0; j < D; j++) {
      zb_ntt[i].poly[j] = con_add((in->zb)[i].poly[j], QHAT);
    }

    ntt_qhat(zb_ntt + i);
  }

  /* ntt(g) */
  for (i = 0; i < N_SPENT + 2; i++) {
    for (j = 0; j < D; j++) {
      tmp_q.poly[j] = con_add(x_ntt_qhat.poly[j] - f_ntt_qhat[i].poly[j], QHAT);
    }

    mult_rqhat(g_ntt + i, f_ntt_qhat + i, &tmp_q);
  }

  /* A = G' * (z_b, f, g) - x * B */
  for (i = 0; i < N_HAT; i++) {
    memcpy(a + i, zb_ntt + i, sizeof(POLY_QHAT));
  }
  for (i = 0; i < K_HAT; i++) {
    for (j = 0; j < N_HAT; j++) {
      mult_plus_rqhat(a + j, mat3[i] + j, zb_ntt + N_HAT + i);
    }
  }
  for (i = 0; i < N_SPENT + 2; i++) {
    for (j = 0; j < N_HAT; j++) {
      mult_plus_rqhat(a + j, mat3[K + i] + j, f_ntt_qhat + i);
    }
  }
  for (i = 0; i < N_SPENT + 2; i++) {
    for (j = 0; j < N_HAT; j++) {
      mult_plus_rqhat(a + j, mat3[K + N_SPENT + 2 + i] + j, g_ntt + i);
    }
  }
  for (i = 0; i < N_HAT; i++) {
    mult_minus_rqhat(a + i, &x_ntt_qhat, (in->b) + i);
  }

  /* \hat{pk}_{i, j} = pk_{i, j} - G(sn_i) */
  sample_g(g_gamma, sn_out[0].sn);
  for (i = 0; i < N_SPENT; i++) {
    for (j = 0; j < N_BAR; j++) {
      for (k = 0; k < D; k++) {
        pk_hat[i][0][j].poly[k] =
            con_add(a_in[i][0].pk[j].poly[k] - g_gamma[j].poly[k], Q);
      }
    }
  }
  sample_g(g_gamma, sn_out[1].sn);
  for (i = 0; i < N_SPENT; i++) {
    for (j = 0; j < N_BAR; j++) {
      for (k = 0; k < D; k++) {
        pk_hat[i][1][j].poly[k] =
            con_add(a_in[i][1].pk[j].poly[k] - g_gamma[j].poly[k], Q);
      }
    }
  }

  /* ntt(\alpha'_0) */
  for (i = 0; i < D; i++) {
    alpha_prime_0_ntt.poly[i] = con_add(alpha_prime_0.poly[i], Q);
  }
  ntt_q(&alpha_prime_0_ntt);

  /* P_j = \sum_{i = 0}^{M - 2} (a_i' * \hat{pk}_{i, j}) + \hat{pk}_{M - 1, j}
   */
  for (i = 0; i < N_SPENT; i++) {
    for (j = 0; j < N_BAR; j++) {
      mult_rq(p[i] + j, &alpha_prime_0_ntt, pk_hat[i][0] + j);

      for (k = 0; k < D; k++) {
        p[i][j].poly[k] = con_sub(p[i][j].poly[k] + pk_hat[i][1][j].poly[k], Q);
      }
    }
  }

  /* ntt(z) */
  for (i = 0; i < N_BAR + K; i++) {
    for (j = 0; j < D; j++) {
      z_ntt[i].poly[j] = con_add((in->z)[i].poly[j], Q);
    }

    ntt_q(z_ntt + i);
  }

  /* E = G * z - \sum_{j = 0}^{N_SPENT - 1} (f_j * P_j) */
  for (i = 0; i < N_BAR; i++) {
    memcpy(e + i, z_ntt + i, sizeof(POLY_Q));
  }
  for (i = 0; i < K; i++) {
    for (j = 0; j < N_BAR; j++) {
      mult_plus_rq(e + j, mat1[i] + j + (N - N_BAR), z_ntt + N_BAR + i);
    }
  }
  for (i = 0; i < N_SPENT; i++) {
    for (j = 0; j < N_BAR; j++) {
      mult_minus_rq(e + j, f_ntt + i, p[i] + j);
    }
  }

  /* x <-- H(\alpha, v_0, v_1, A, B, E) */
  hash_x(&x, alpha_verify, &v0, &(in->v1), a, in->b, e);
  for (i = 0; i < D; i++) {
    if (x.poly[i] != (in->x).poly[i]) {
      return 0;
    }
  }

  return 1;
}

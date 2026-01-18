#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "poly_mult.h"
#include "poly_param.h"
#include "poly_q.h"
#include "poly_red.h"

/* c=c+a*b on R_q */
void mult_plus_rq(POLY_Q *c, const POLY_Q *a, const POLY_Q *b) {
  POLY_Q tmp;
  uint64_t i;

  mult_rq(&tmp, a, b);

  for (i = 0; i < D; i++) {
    c->poly[i] = con_sub(c->poly[i] + tmp.poly[i], Q);
  }
}

/* c=c-a*b on R_q */
void mult_minus_rq(POLY_Q *c, const POLY_Q *a, const POLY_Q *b) {
  POLY_Q tmp;
  uint64_t i;

  mult_rq(&tmp, a, b);

  for (i = 0; i < D; i++) {
    c->poly[i] = con_add(c->poly[i] - tmp.poly[i], Q);
  }
}

/* c=c+a*b on R_\hat{q} */
void mult_plus_rqhat(POLY_QHAT *c, const POLY_QHAT *a, const POLY_QHAT *b) {
  POLY_QHAT tmp;
  uint64_t i;

  mult_rqhat(&tmp, a, b);

  for (i = 0; i < D; i++) {
    c->poly[i] = con_sub(c->poly[i] + tmp.poly[i], QHAT);
  }
}

/* c=c-a*b on R_{\hat{q}} */
void mult_minus_rqhat(POLY_QHAT *c, const POLY_QHAT *a, const POLY_QHAT *b) {
  POLY_QHAT tmp;
  uint64_t i;

  mult_rqhat(&tmp, a, b);

  for (i = 0; i < D; i++) {
    c->poly[i] = con_add(c->poly[i] - tmp.poly[i], QHAT);
  }
}

/* c=c+a*b when pm=0 or c=c-a*b when pm=1 on R_{\hat{q}} */
void mult_plus_rqhat_pm(POLY_QHAT *c, const POLY_QHAT *a, const POLY_QHAT *b,
                        const uint64_t pm) {
  POLY_QHAT tmp;
  uint64_t i;
  int64_t pm1 = 1 - (pm << 1);

  mult_rqhat(&tmp, a, b);

  for (i = 0; i < D; i++) {
    c->poly[i] = con_sub(con_add(c->poly[i] + pm1 * tmp.poly[i], QHAT), QHAT);
  }
}

/* out=a*b on R (no modulo reduction) */
void mult_r(POLY_R *out, const POLY_R *a, const POLY_R *b) {
  uint64_t tmp[D << 1];

  uint64_t i, j;

  memset(tmp, 0, sizeof(tmp));
  for (i = 0; i < D; i++) {
    for (j = 0; j < D; j++) {
      tmp[i + j] += ((int64_t)a->poly[i]) * ((int64_t)b->poly[j]);
    }
  }

  for (i = 0; i < D; i++) {
    out->poly[i] = tmp[i] - tmp[i + D];
  }
}

/* out=f*(x-f) on R (no modulo reduction) */
void mult_fxf(POLY_R *out, const POLY_R *f, const POLY_R *x) {
  POLY_R fx;

  uint64_t i;

  for (i = 0; i < D; i++) {
    fx.poly[i] = x->poly[i] - f->poly[i];
  }

  mult_r(out, f, &fx);
}

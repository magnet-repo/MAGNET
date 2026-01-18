#ifndef _POLY_MULT_H
#define _POLY_MULT_H

#include "poly_q.h"
#include <stdint.h>

void mult_plus_rq(POLY_Q *c, const POLY_Q *a, const POLY_Q *b);
void mult_minus_rq(POLY_Q *c, const POLY_Q *a, const POLY_Q *b);

void mult_plus_rqhat(POLY_QHAT *c, const POLY_QHAT *a, const POLY_QHAT *b);
void mult_minus_rqhat(POLY_QHAT *c, const POLY_QHAT *a, const POLY_QHAT *b);
void mult_plus_rqhat_pm(POLY_QHAT *c, const POLY_QHAT *a, const POLY_QHAT *b,
                        const uint64_t pm);

void mult_r(POLY_R *out, const POLY_R *a, const POLY_R *b);
void mult_fxf(POLY_R *out, const POLY_R *f, const POLY_R *x);

#endif

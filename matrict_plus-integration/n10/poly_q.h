#ifndef _POLY_Q
#define _POLY_Q

#include "poly_param.h"
#include <stdint.h>

typedef struct {
  uint64_t poly[D];
} POLY_Q;

typedef POLY_Q POLY_R;
typedef POLY_Q POLY_QHAT;

void ntt_q(POLY_Q *a);
void mult_rq(POLY_Q *out, const POLY_Q *a, const POLY_Q *b);

void ntt_qhat(POLY_QHAT *a);
void mult_rqhat(POLY_QHAT *out, const POLY_QHAT *a, const POLY_QHAT *b);

void aut(POLY_Q *out, const POLY_Q *a);
void iaut_r(POLY_R *out, const POLY_R *a);

#endif

#ifndef _HASH_H
#define _HASH_H

#include "param.h"
#include "poly_q.h"
#include "spend.h"

void hash_alpha(unsigned char *out, const unsigned char *pp,
                const ACT a_in[][M], const SN_OUT *sn_out, const PK_OUT *pk_out,
                const CN_OUT *cn_out, const POLY_Q w_out[][N], const POLY_Q *w0,
                const POLY_Q *v, const POLY_Q *v2, const POLY_Q *c,
                const POLY_Q *w, const POLY_Q *w1, const POLY_Q *w2);
void hash_x(POLY_R *out, const unsigned char *alpha, const POLY_Q *v0,
            const POLY_Q *v1, const POLY_QHAT *a, const POLY_QHAT *b,
            const POLY_Q *e);

#endif

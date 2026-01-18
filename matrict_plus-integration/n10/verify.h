#ifndef _VERIFY_H
#define _VERIFY_H

#include "param.h"
#include "poly_q.h"
#include "spend.h"
#include <stdint.h>

uint64_t verify(const ACT a_in[][M], const PK_OUT *pk_out, const CN_OUT *cn_out,
                const SN_OUT *sn_out, const SPEND_OUT *in,
                const POLY_Q mat1[][N], const POLY_Q mat2[][3],
                const POLY_QHAT mat3[][N_HAT], unsigned char *pp);

#endif

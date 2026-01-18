#ifndef _MINT_H
#define _MINT_H

#include "param.h"
#include "poly_q.h"
#include <stdint.h>

void mint(POLY_Q *cn, POLY_R *ck, uint64_t *a_hat, const uint64_t amt,
          const POLY_Q mat1[][N], const POLY_Q mat2[][3]);

#endif

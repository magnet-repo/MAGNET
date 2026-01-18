#ifndef _KEYGEN_H
#define _KEYGEN_H

#include "param.h"
#include "poly_q.h"

void keygen(POLY_Q *pk, POLY_R *sk, unsigned char *sn, const POLY_Q mat1[][N]);

#endif

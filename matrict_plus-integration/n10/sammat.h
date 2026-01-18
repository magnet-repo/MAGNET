#ifndef _SAMMAT_H
#define _SAMMAT_H

#include "param.h"
#include "poly_q.h"

void sample_mat1(POLY_Q out[][N]);
void sample_mat1_keygen(POLY_Q out[][N]);
void sample_mat2(POLY_Q out[][3]);
void sample_mat3(POLY_QHAT out[][N_HAT]);

#endif

#ifndef _GAUSSIAN_AVX_H
#define _GAUSSIAN_AVX_H

#include "param.h"
#include "poly_q.h"
#include <stdint.h>

void sample_b_nk1(POLY_R *out, int64_t *samps, uint32_t *ct);
void sample_b_nk3(POLY_R *out, int64_t *samps, uint32_t *ct);
void sample_b_nhatkhat(POLY_R *out, int64_t *samps, uint32_t *ct);
void sample_brs_13_nbark(POLY_R *out);
void sample_ba_n1(POLY_R *out);

uint64_t rej_f(const POLY_R *z, const POLY_R *c);
uint64_t rej_z(const POLY_R *z, const POLY_R *c);

uint64_t rej_op(const POLY_R z_out[][N + K + 1], const POLY_R *z1,
                const POLY_R *z2, const POLY_R *z0, const POLY_R *zb,
                const POLY_R x_rout[][N + K + 1], const POLY_R *x_rho,
                const POLY_R *x_aut_rho, const POLY_R *x_ri,
                const POLY_R *x_rb);

#endif

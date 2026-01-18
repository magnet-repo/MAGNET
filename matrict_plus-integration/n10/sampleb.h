#ifndef _SAMPLEB_H
#define _SAMPLEB_H

#include "poly_q.h"
#include <stdint.h>

void sample_1(POLY_R *out, const uint64_t m);
void sample_g(POLY_Q *out, const unsigned char *seed);
void sample_alpha(POLY_Q *out);
void sample_alpha_prime(POLY_R *out);

#endif

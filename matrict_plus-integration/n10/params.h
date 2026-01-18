#ifndef PARAMS_H_
#define PARAMS_H_

#include <stdint.h>
#include <stdlib.h>

#define NN 256 * 59 // 59 = 3*(N+K+1) + 2*(N+K+3) + (N_HAT+K_HAT)
#define LMD 64
#define PS 0.76

#define ITER 5

#define NUM_SHARES 2
#define MASK_ORDER NUM_SHARES - 1
#define ww 32
#define WW 5 // ceil(log2(w-1))

#define rand_uint32() xoshiro_next()

//#define VERIFY
//#define VERIFY2

extern int kap, l, mu, m, SIG;
extern double r, p, dlap_bias;
extern uint32_t *p_lap, *p_geo, *p_exp;
extern uint32_t one[NUM_SHARES], mSIG[NUM_SHARES];

void init_parameters();

#endif /* PARAMS_H_ */

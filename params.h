#ifndef PARAMS_H_
#define PARAMS_H_

#define N 10000
#define LMD 64
#define PS 0.76

#define ITER 25

#define NUM_SHARES 2
#define MASK_ORDER NUM_SHARES - 1
#define w 32
#define W 5 // ceil(log2(w-1))

#define rand_uint32() xoshiro_next()

//#define VERIFY

extern int kap, l, mu, m, SIG;
extern double r, p, dlap_bias;

void init_parameters(int j);

#endif /* PARAMS_H_ */

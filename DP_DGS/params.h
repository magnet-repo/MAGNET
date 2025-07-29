#ifndef PARAMS_H_
#define PARAMS_H_

#define N 10000
#define LMD 64
#define PS 0.76

#define ITER 25

#define rand_uint32() xoshiro_next()

extern int kap, l, mu, m, SIG;
extern double r, p, dlap_bias;

//#define VERIFY

void init_parameters(int j);

#endif /* PARAMS_H_ */

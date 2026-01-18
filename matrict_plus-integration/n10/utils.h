#ifndef UTILS_H_
#define UTILS_H_

#include <stdint.h>

void GenUniRandBit(uint32_t* x);
void print_shares(uint32_t* x);
void print_value(uint32_t* x);
void bin_rep(uint32_t *b, double p);
void bin_rep_arr(uint32_t *b, int n, double *arr);
int cmp_double(const void *a, const void *b);
double calc_median(double values[ITER]);

void seed_xoshiro(void);
uint32_t xoshiro_next(void);

#endif /* UTILS_H_ */

#ifndef UTILS_H_
#define UTILS_H_

#include <stdint.h>

uint32_t* bin_rep(double p);
void bin_rep_arr(uint32_t *b_arr, int n, double *arr);
int cmp_double(const void *a, const void *b);
double calc_median(double values[ITER]);

void seed_xoshiro(void);
uint32_t xoshiro_next(void);

#endif /* UTILS_H_ */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#include "params.h"
#include "utils.h"

uint32_t* bin_rep(double p) {
	uint32_t *b = (uint32_t*)malloc(mu * sizeof(uint32_t));
	double t = p;
	for(int i = 0; i < mu; i++) {
		t *= 2;
		b[i] = (uint32_t)t;
		t -= b[i];
	}
	return b;
}

void bin_rep_arr(uint32_t *b_arr, int n, double *arr) {
	for(int i = 0; i < n; ++i) {
		double t = arr[i];
		uint32_t b = 0;
		for(int j = 0; j < mu; j++) {
			t *= 2;
			b = (uint32_t)t;
			t -= b;
			b_arr[i * mu + j] = b;
		}
	}
}

int cmp_double(const void *a, const void *b) {
    double da = *(const double*)a;
    double db = *(const double*)b;
    if (da < db) return -1;
    if (da > db) return  1;
    return 0;
}

double calc_median(double values[ITER]) {
    qsort(values, ITER, sizeof(values[0]), cmp_double);

    double mid1 = values[ITER/2 - 1];
    double mid2 = values[ITER/2];
    return (mid1 + mid2) / 2.0;
}

/* ------------------------------------------- */

// http://xoshiro.di.unimi.it/xoshiro128starstar.c

/*  Written in 2018 by David Blackman and Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

/* This is xoshiro128** 1.0, our 32-bit all-purpose, rock-solid generator. It
   has excellent (sub-ns) speed, a state size (128 bits) that is large
   enough for mild parallelism, and it passes all tests we are aware of.

   For generating just single-precision (i.e., 32-bit) floating-point
   numbers, xoshiro128+ is even faster.

   The state must be seeded so that it is not everywhere zero. */

static uint32_t s[4];

static inline uint32_t rotl(const uint32_t x, int k) {
	return (x << k) | (x >> (32 - k));
}


uint32_t xoshiro_next(void) {
	const uint32_t result_starstar = rotl(s[0] * 5, 7) * 9;

	const uint32_t t = s[1] << 9;

	s[2] ^= s[0];
	s[3] ^= s[1];
	s[1] ^= s[2];
	s[0] ^= s[3];

	s[2] ^= t;

	s[3] = rotl(s[3], 11);

	return result_starstar;
}

void seed_xoshiro(void){
    srand(time(NULL));
    s[0] = rand();
    s[1] = rand();
    s[2] = rand();
    s[3] = rand();
}

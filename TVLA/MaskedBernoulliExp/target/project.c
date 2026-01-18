#include "hal.h"
#include "simpleserial.h"
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define NUM_SHARES 2

#define N 10
#define SIG 1
#define LMD 32
#define PS 0.76
#define mu 2
#define l 2
#define NUM_SHARES 2
#define rand_uint32() xoshiro_next()

#define BER_RAND mu
#define BEXP_RAND l * (BER_RAND + 3) + 20

static inline uint32_t rotl(const uint32_t x, int k) {
	return (x << k) | (x >> (32 - k));
}

static uint32_t s[4];

void seed_xoshiro(uint8_t* pt){
    s[0] = ((uint32_t)pt[11] << 24) | ((uint32_t)pt[10] << 16) | ((uint32_t)pt[9] <<  8) | ((uint32_t)pt[8] );
    s[1] = ((uint32_t)pt[15] << 24) | ((uint32_t)pt[14] << 16) | ((uint32_t)pt[13] <<  8) | ((uint32_t)pt[12] );
    s[2] = ((uint32_t)pt[11] << 24) | ((uint32_t)pt[10] << 16) | ((uint32_t)pt[9] <<  8) | ((uint32_t)pt[8] );
    s[3] = ((uint32_t)pt[15] << 24) | ((uint32_t)pt[14] << 16) | ((uint32_t)pt[13] <<  8) | ((uint32_t)pt[12] );
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

void GenUniRandBit(uint32_t* x) {
	for(int i = 0; i < NUM_SHARES; ++i) x[i] = rand_uint32() & 1;
}

void Refresh(uint32_t* x) {
	for(int i = 0; i < NUM_SHARES; ++i) {
		for(int j = i + 1; j < NUM_SHARES; ++j) {
			uint32_t r = rand_uint32();
			x[i] ^= r;
			x[j] ^= r;
		}
	}
}

void bin_rep_arr(uint32_t *b_msk, int n, double *arr) {
	for(int i = 0; i < n; ++i) {
		double t = arr[i];
		uint8_t b = 0;
		for(int j = 0; j < mu; ++j) {
			t *= 2;
			b = (uint8_t)t;
			t -= b;
			int idx = i * (mu * NUM_SHARES) + j * NUM_SHARES;
			b_msk[idx] = b;
			for(int k = 1; k < NUM_SHARES; k++) {
				b_msk[idx + k] = 0;
			}
			Refresh(&b_msk[idx]);
		}
	}
}

extern void MaskedBernoulliExp_m4(uint32_t *b, uint32_t *u, uint32_t *p_exp, uint32_t *one, uint32_t *rs, uint32_t *pool);

uint8_t get_pt(uint8_t* pt, uint8_t len) {
    seed_xoshiro(pt);
    double r = 2 * SIG * SIG;

    double Bexp_bias[l];
    for(int i = 0; i < l; ++i) Bexp_bias[i] = exp(-(1 << i) / r);
    
    uint32_t *p_exp = (uint32_t*)malloc((l * mu * NUM_SHARES) * sizeof(uint32_t));    
    bin_rep_arr(p_exp, l, Bexp_bias);
    
    uint32_t one[NUM_SHARES] = {0}, u[NUM_SHARES] = {0};
    one[0] = 1;
	Refresh(one);

    u[0] = ((uint32_t)pt[3] << 24) | ((uint32_t)pt[2] << 16) | ((uint32_t)pt[1] << 8) | pt[0];
    u[1] = ((uint32_t)pt[7] << 24) | ((uint32_t)pt[6] << 16) | ((uint32_t)pt[5] << 8) | pt[4];

    uint8_t res[4] = {0};
	uint32_t rr[NUM_SHARES] = {0}, zp[NUM_SHARES] = {0}; int j = 0;
    uint32_t rs[NUM_SHARES * l * mu] = {0}, pool[BEXP_RAND] = {0};
	for(int i = 0; i < l * mu; ++i) {
	    GenUniRandBit(rr);
	    Refresh(rr);
        rs[j] = rr[0]; j++;
        rs[j] = rr[1]; j++;
	}
	for(int i = 0; i < BEXP_RAND; ++i) {
        pool[i] = rand_uint32();
    }
     
    for (volatile int k = 0; k < 1000; k++) {;} // to clean the power trace

    trigger_high();

    MaskedBernoulliExp_m4(zp, u, p_exp, one, rs, pool);   

    trigger_low();

    for (volatile int k = 0; k < 100; k++) {;}  // to clean the power trace

    uint32_t ans = (zp[0] ^ zp[1]);

    res[0] = ans & 0xff;
    simpleserial_put('r', 1, res);

    free(p_exp);

    return 0x00;
}

int main(void) {
    platform_init();
    init_uart();
    trigger_setup();

    simpleserial_init();

    simpleserial_addcmd('p', 16,  get_pt);
        
    while(1)
        simpleserial_get();
}


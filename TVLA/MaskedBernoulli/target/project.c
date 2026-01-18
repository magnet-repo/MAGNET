#include "hal.h"
#include "simpleserial.h"
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define NUM_SHARES 2

#define N 100
#define SIG 16
#define LMD 16
#define PS 0.76

#define BER_RAND 20

#define rand_uint32() xoshiro_next()

int kap, l, mu, m;
double N0, t, r, p, p0, dlap_bias, k1, k2;

void init_parameters(void) {
	N0 = sqrt(2 * log(2) * (LMD + 2 + log2(N))) * SIG;
	kap = (int)ceil(log2(N0 - 1));
	l = 2 * kap;
	t = (SIG * SIG) / round(SIG);
	r = 2 * SIG * SIG;
	p = exp(-1 / t);
	p0 = PS - pow(2, -LMD);
	mu = 3; //(int)ceil(LMD + 2 + log2((2.0 * (kap + 1) + l) * N / p0));
	dlap_bias = 0.5; //(1 - p) / (1 + p - 2 * exp(-(pow(2, kap) + 1) / t));
	k1 = N / p0;
	k2 = ((LMD + 2) * log(2)) / (2 * p0 * p0);
	m = (int)ceil(k1 + k2 / 2 + sqrt(k2 * k2 / 4 + k1 * k2));
}

static inline uint32_t rotl(const uint32_t x, int k) {
	return (x << k) | (x >> (32 - k));
}

static uint32_t s[4];

void seed_xoshiro(uint8_t* pt){
    s[0] = ((uint32_t)pt[3] << 24) | ((uint32_t)pt[2] << 16) | ((uint32_t)pt[1] <<  8) | ((uint32_t)pt[0] );
    s[1] = ((uint32_t)pt[7] << 24) | ((uint32_t)pt[6] << 16) | ((uint32_t)pt[5] <<  8) | ((uint32_t)pt[4] );
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

void Refresh(uint32_t* x) { // [ISW03]
	for(int i = 0; i < NUM_SHARES; ++i) {
		for(int j = i + 1; j < NUM_SHARES; ++j) {
			uint32_t r = rand_uint32();
			x[i] ^= r;
			x[j] ^= r;
		}
	}
}

void bin_rep(uint32_t *b_msk, double p) {
	double t = p;
	uint8_t b = 0;
	for(int i = 0; i < mu; ++i) {
		t *= 2;
		b = (uint8_t)t;
		t -= b;
		b_msk[i * NUM_SHARES] = (uint32_t)b;
		for(int j = 1; j < NUM_SHARES; ++j) {
			b_msk[(i * NUM_SHARES) + j] = 0;
		}
		Refresh(&b_msk[i * NUM_SHARES]);
	}
}

extern void MaskedBernoulli_m4(uint32_t *zp, uint32_t *p_lap, uint32_t *one, uint32_t *rs, uint32_t *pool);

uint8_t get_pt(uint8_t* pt, uint8_t len) {
    seed_xoshiro(pt);
    init_parameters();
    
    uint32_t *p_lap = (uint32_t*)malloc((mu * NUM_SHARES) * sizeof(uint32_t));    
    bin_rep(p_lap, dlap_bias);
    
    uint32_t one[NUM_SHARES] = {0};
    one[0] = 1;
	Refresh(one);    
    
    uint8_t res[4];
	uint32_t r[NUM_SHARES] = {0}, zp[NUM_SHARES] = {0}; int j = 0;
    uint32_t pool[BER_RAND], rs[NUM_SHARES * mu];
	for(int i = 0; i < mu; ++i) {
	    GenUniRandBit(r);
	    Refresh(r);
        rs[j] = r[0]; j++;
        rs[j] = r[1]; j++;
	}

	for(int i = 0; i < BER_RAND; ++i) {	
	    pool[i] = rand_uint32();
	}
        
    for (volatile int k = 0; k < 1000; k++) {;} // to clean the power trace

    trigger_high();

    MaskedBernoulli_m4(zp, p_lap, one, rs, pool);   

    trigger_low();

    for (volatile int k = 0; k < 100; k++) {;}  // to clean the power trace

    uint32_t ans = (zp[0] ^ zp[1]);

    res[0] = ans & 0xff;
    simpleserial_put('r', 1, res);

    free(p_lap);

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


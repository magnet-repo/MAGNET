#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <zephyr/kernel.h>
#include <zephyr/drivers/entropy.h>
#include <zephyr/arch/cpu.h>

int kap, l, mu, m;
double r, p, dlap_bias;
uint32_t s[4]; 

#define N_ 100
#define LMD 64
#define PS 0.76
#define SIG 256
#define NUM_SHARES 2 

#define rand_uint32() xoshiro_next_c()

static void enable_cycle_counter(void) {
    CoreDebug->DEMCR |= CoreDebug_DEMCR_TRCENA_Msk;
    DWT->CYCCNT = 0;
    DWT->CTRL |= DWT_CTRL_CYCCNTENA_Msk; 
}

static inline uint32_t get_cycle_count(void) {
    return DWT->CYCCNT;
}

void init_parameters() {
    double N0 = sqrt(2 * log(2) * (LMD + 2 + log2(N_))) * SIG;
    kap = (int)ceil(log2(N0 - 1));
    l = 2 * kap;
    double t_val = (double)(SIG * SIG) / round(SIG);
    r = 2 * SIG * SIG;
    p = exp(-1 / t_val);
    double p0 = PS - pow(2, -LMD);
    mu = (int)ceil(LMD + 2 + log2((2.0 * (kap + 1) + l) * N_ / p0));
    dlap_bias = (1 - p) / (1 + p - 2 * exp(-(pow(2, kap) + 1) / t_val));
    double k1 = N_ / p0;
    double k2 = ((LMD + 2) * log(2)) / (2 * p0 * p0);
    m = (int)ceil(k1 + k2 / 2 + sqrt(k2 * k2 / 4 + k1 * k2));
}

void seed_xoshiro(void) {
    const struct device *entropy_dev = DEVICE_DT_GET_ONE(nordic_nrf_rng);    
    if(!device_is_ready(entropy_dev)) {
        printf("TRNG device not ready\n");
        return;
    }
    entropy_get_entropy(entropy_dev, (uint8_t*)s, sizeof(s));
}

static inline uint32_t rotl_c(const uint32_t x, int k) {
	return (x << k) | (x >> (32 - k));
}

uint32_t xoshiro_next_c(void) {
	const uint32_t result_starstar = rotl_c(s[1] * 5, 7) * 9;

	const uint32_t t = s[1] << 9;

	s[2] ^= s[0];
	s[3] ^= s[1];
	s[1] ^= s[2];
	s[0] ^= s[3];

	s[2] ^= t;

	s[3] = rotl_c(s[3], 11);

	return result_starstar;
}

void Refresh_c(uint32_t* x) {
	for(int i = 0; i < NUM_SHARES; ++i) {
		for(int j = i + 1; j < NUM_SHARES; ++j) {
			uint32_t r = rand_uint32();
			x[i] ^= r;
			x[j] ^= r;
		}
	}
}

void FullRefresh(uint32_t* z, const uint32_t* x) {
	for(int i = 0; i < NUM_SHARES; ++i) z[i] = x[i];
	for(int i = 0; i < NUM_SHARES; ++i) {
		for(int j = 1; j < NUM_SHARES; ++j) {
			uint32_t r = rand_uint32();
			z[i] ^= r;
			z[j] ^= r;
		}
	}
}

uint32_t FullXOR(uint32_t* x) {
	uint32_t z[NUM_SHARES] = {0};
	FullRefresh(z, x);
	uint32_t val = z[0];
	for(int i = 1; i < NUM_SHARES; ++i) val ^= z[i];

	return val;
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
		Refresh_c(&b_msk[i * NUM_SHARES]);
	}
}

void bin_rep_arr(uint32_t *b_msk, int n, double *arr) {
    for(int i = 0; i < n; ++i) {
        double t_v = arr[i];
        uint8_t b = 0;
        for(int j = 0; j < mu; ++j) {
            t_v *= 2; b = (uint8_t)t_v; t_v -= b;
            int idx = i * (mu * NUM_SHARES) + j * NUM_SHARES;
            b_msk[idx] = b;
            b_msk[idx + 1] = 0;
        }
    }
}

extern void MAGNET_m4(uint32_t *samps, uint32_t *mSIG, uint32_t *one, uint32_t *p_geo, uint32_t *p_exp, uint32_t *p_lap);

int main(void) {
    seed_xoshiro();
    init_parameters();
    
    printf("\nRunning MAGNET First-Order Cortex-M4 Implementation with N=%d, SIG=%d\n", N_, SIG);
    printf("Params: m=%d, mu=%d, kap=%d, l=%d\n", m, mu, kap, l);
    
    double geo_bias[kap], Bexp_bias[l];
    for(int i = 0; i < kap; ++i) geo_bias[i] = pow(p, (1 << i)) / (1 + pow(p, (1 << i)));
    for(int i = 0; i < l; ++i) Bexp_bias[i] = exp(-(1 << i) / r);
      
    uint32_t *p_lap = (uint32_t*)malloc((mu * NUM_SHARES) * sizeof(uint32_t));
    uint32_t *p_geo = (uint32_t*)malloc((kap * mu * NUM_SHARES) * sizeof(uint32_t));
    uint32_t *p_exp = (uint32_t*)malloc((l * mu * NUM_SHARES) * sizeof(uint32_t));
    uint32_t *samps = (uint32_t*)malloc((m * NUM_SHARES) * sizeof(uint32_t));
    
    if (!p_lap || !p_geo || !p_exp || !samps) {
        printf("Error: OOM\n"); return -1;
    }

    memset(samps, 0, (m * NUM_SHARES) * sizeof(uint32_t));
    bin_rep(p_lap, dlap_bias);
    bin_rep_arr(p_geo, kap, geo_bias);
    bin_rep_arr(p_exp, l, Bexp_bias);

    uint32_t one[NUM_SHARES] = {0}, mSIG[NUM_SHARES] = {0};
    one[0] = 1; 
    mSIG[0] = -SIG; 
    Refresh_c(one); 
    Refresh_c(mSIG);

    enable_cycle_counter();

    uint32_t start_cycles = get_cycle_count();

    MAGNET_m4(samps, mSIG, one, p_geo, p_exp, p_lap);

    uint32_t end_cycles = get_cycle_count();
    uint32_t total_cycles = end_cycles - start_cycles;
        
    printf("--- Benchmark Results ---\n");
    printf("Total Cycles: %u\n", total_cycles);
    printf("Cycles per Sample (m=%d): %u\n", m, total_cycles / m);
    
    free(samps); 
    free(p_exp); 
    free(p_geo); 
    free(p_lap);       
                   
    return 0;
}



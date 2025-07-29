#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "params.h"
#include "utils.h"
#include "algorithms.h"

uint32_t SamplingCircuit_B(uint32_t *b) { // Algo 1
	uint32_t x = 1; // BIT(1)
	for(int i = mu - 1; i >= 0; i--) {
		uint32_t r = rand_uint32() & 1; // URBIT
		uint32_t t = b[i] ^ r; // XOR
		uint32_t t1 = x ^ !r;
		uint32_t t2 = t & t1;
		x = x ^ t2;
	}
	return x;
}

uint32_t SamplingCircuit_G(uint32_t *p_geo) { // Algo 2
	uint32_t *b = (uint32_t*)malloc(kap * sizeof(uint32_t));
	for(int i = 0; i < kap; i++) {
		b[i] = SamplingCircuit_B(&p_geo[i * mu]);
	}
	uint32_t x = 0;
	for(int i = 0; i < kap; i++) { // COMPOSE
		x |= (b[i] << i);
	}
	return x;
}

void SamplingCircuit_L(uint32_t *v, uint32_t *p_geo, uint32_t *p_lap) { // Algo 3
	uint32_t b = SamplingCircuit_B(p_lap);
	uint32_t x = SamplingCircuit_G(p_geo);
	uint32_t t1 = -(!b);
	uint32_t t2 = x + 1;
	*v = t1 & t2;
}

void SamplingCircuit_Bexp(uint32_t *res, uint32_t *p_exp, uint32_t u) { // Algo 4
	uint32_t b = 1;
	for(int i = 0; i < l; i++) {
		uint32_t bd = SamplingCircuit_B(&p_exp[i * mu]);
		uint32_t ui = (u >> i) & 1;
		uint32_t v = (!ui) | bd; // OR
		b &= v; // AND
	}
	*res = b;
}

void SamplingCircuit_N(uint32_t *samps, uint32_t *p_geo, uint32_t *p_exp, uint32_t *p_lap) { // Algo 5
	int j = 0;
	uint32_t u, v, s, b, sign;
	for(int i = 0; i < m; i++) {
		SamplingCircuit_L(&s, p_geo, p_lap);
		v = s - round(SIG);
		u = v * v;

		SamplingCircuit_Bexp(&b, p_exp, u);

		sign = rand_uint32() & 1; // URBIT
		s ^= -sign;
		s += sign;

		uint32_t bitmask = -b;
		samps[j] ^= (samps[j] ^ s) & bitmask;
		j += b;
	}

#ifdef VERIFY
	FILE *outfile = fopen("samples.txt", "w");

	for (int i = 0; i < N; i++) {
		fprintf(outfile, "%d ", (int32_t)samps[i]);
	}
	fprintf(outfile, "%d ", (int32_t)(SIG * SIG));
	fprintf(outfile, "%d ", N);

	fclose(outfile);

	exit(system("python3 discretegauss.py"));
#endif
}

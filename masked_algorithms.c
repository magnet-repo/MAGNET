#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "masked_algorithms.h"
#include "gadgets.h"

void MaskedBernoulli(uint32_t *x, uint32_t *p, uint32_t *one) {
	uint32_t r[NUM_SHARES] = {0}, r_prm[NUM_SHARES] = {0}, t[NUM_SHARES] = {0}, u[NUM_SHARES] = {0}, v[NUM_SHARES] = {0};
	
	for(int i = 0; i < NUM_SHARES; ++i) x[i] = one[i];

	for(int j = mu - 1; j >= 0; j--) {
        GenUniRandBit(r);
        Refresh(r);

        for(int i = 0; i < NUM_SHARES; ++i) t[i] = p[(j * NUM_SHARES) + i] ^ r[i];

        r_prm[0] = ~r[0];
	    for(int i = 1; i < NUM_SHARES; ++i) r_prm[i] = r[i];

        for(int i = 0; i < NUM_SHARES; ++i) u[i] = x[i] ^ r_prm[i];
        
		SecAND(v, t, u);
		
		for(int i = 0; i < NUM_SHARES; ++i) x[i] ^= v[i];
		
	}
}

void MaskedGeometric(uint32_t *x, uint32_t *p_geo, uint32_t *one) {
	for(int j = 0; j <= kap - 1; ++j) {
		uint32_t b[NUM_SHARES] = {0}, t[NUM_SHARES] = {0}, tmp[NUM_SHARES] = {0};
		
		MaskedBernoulli(b, &p_geo[j * (mu * NUM_SHARES)], one);

		for(int i = 0; i < NUM_SHARES; ++i) t[i] = b[i] << j;
		
		SecOR(tmp, x, t);
		for(int i = 0; i < NUM_SHARES; ++i) x[i] = tmp[i];
		
		Refresh(x);
	}
}


void MaskedLaplace(uint32_t *v, uint32_t *p_lap, uint32_t *p_geo, uint32_t *one) {
	uint32_t b[NUM_SHARES] = {0}, x[NUM_SHARES] = {0}, b_prm[NUM_SHARES] = {0}, t[NUM_SHARES] = {0}, u[NUM_SHARES] = {0};
	MaskedBernoulli(b, p_lap, one);
	MaskedGeometric(x, p_geo, one);

    b_prm[0] = ~b[0];
    for(int j = 1; j < NUM_SHARES; ++j) b_prm[j] = b[j];
    
    SecNEG(t, b_prm);
	SecINC(u, x, one);
	SecAND(v, t, u);
}

void MaskedBernoulliExp(uint32_t *b, uint32_t *u, uint32_t *p_exp, uint32_t *one) {
	uint32_t b_prm[NUM_SHARES] = {0}, uj[NUM_SHARES] = {0}, uj_til[NUM_SHARES] = {0}, v[NUM_SHARES] = {0}, tmp[NUM_SHARES] = {0};
	
	for(int i = 0; i < NUM_SHARES; ++i) b[i] = one[i];

	for(int j = 0; j <= l - 1; ++j) {
		MaskedBernoulli(b_prm, &p_exp[j * (mu * NUM_SHARES)], one);

		for(int i = 0; i < NUM_SHARES; ++i) uj[i] = (u[i] >> j);
		Refresh(uj);

	    uj_til[0] = ~uj[0];
        for(int i = 1; i < NUM_SHARES; ++i) uj_til[i] = uj[i];
		SecOR(v, uj_til, b_prm);
		SecAND(tmp, b, v);
		
		for(int i = 0; i < NUM_SHARES; ++i) b[i] = tmp[i];
	}
}

void MAGNET(uint32_t *samps, uint32_t *mSIG, uint32_t *one, uint32_t *p_geo,  uint32_t *p_exp, uint32_t *p_lap) {
	uint32_t sj[NUM_SHARES] = {0}, b[NUM_SHARES] = {0}, V[NUM_SHARES] = {0}, v[NUM_SHARES] = {0}, U[NUM_SHARES] = {0}, u[NUM_SHARES] = {0};
	uint32_t t[NUM_SHARES] = {0}, tmp[NUM_SHARES] = {0}, b_prm[NUM_SHARES] = {0}, sign[NUM_SHARES] = {0}, k = 0;
	
	uint32_t sj_til[NUM_SHARES] = {0}, diffmask[NUM_SHARES] = {0},  signmask[NUM_SHARES] = {0}, selmask[NUM_SHARES] = {0};
	
	for(int j = 0; j <= m - 1; ++j) {
		MaskedLaplace(sj, p_lap, p_geo, one);

		SecADD(v, sj, mSIG);

		SecB2A(V, v, NUM_SHARES); // [BCZ18]
		SecSQU(U, V);
		SecA2B(U, u, w); // [CGTZ23]

		MaskedBernoulliExp(b, u, p_exp, one);

		GenUniRandBit(sign);
		Refresh(sign);

        SecNEG(signmask, sign);

        for(int i = 0; i < NUM_SHARES; ++i) sj[i] ^= signmask[i];
	    SecINC(sj_til, sj, one);
        
        for(int i = 0; i < NUM_SHARES; ++i) diffmask[i] = sj_til[i] ^ sj[i];
        SecAND(selmask, diffmask, signmask);

        for(int i = 0; i < NUM_SHARES; ++i) t[i] = sj[i] ^ selmask[i];
        SecNEG(b_prm, b);
		SecAND(tmp, b_prm, t);
		
		for(int i = 0; i < NUM_SHARES; ++i) samps[(k * NUM_SHARES) + i] ^= tmp[i];

		k += FullXOR(b);
	}

#ifdef VERIFY
	FILE *outfile = fopen("samples.txt", "w");

	for (int i = 0; i < N; i++) {
		fprintf(outfile, "%d ", (int32_t)FullXOR(&samps[i * NUM_SHARES]));
	}
	fprintf(outfile, "%d ", (int32_t)(SIG * SIG));
	fprintf(outfile, "%d ", N);

	fclose(outfile);

	exit(system("python3 discretegauss.py")); // plotting samples using the code from "The Discrete Gaussian for Differential Privacy [CKS20]"
#endif
}



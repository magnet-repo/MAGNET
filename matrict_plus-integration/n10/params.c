#include <math.h>

#include "params.h"
#include "utils.h"
#include "gadgets.h"

int kap, l, mu, m, SIG;
double r, p, dlap_bias;
uint32_t one[NUM_SHARES] = {0}, mSIG[NUM_SHARES] = {0};
uint32_t *p_lap = NULL;
uint32_t *p_geo = NULL;
uint32_t *p_exp = NULL;

void init_parameters() {
	SIG = 6390;
	
	double N0 = sqrt(2 * log(2) * (LMD + 2 + log2(NN))) * SIG;
	kap = (int)ceil(log2(N0 - 1));
	l = 2 * kap;
	double t = (SIG * SIG) / round(SIG);
	r = 2 * SIG * SIG;
	p = exp(-1 / t);
	double p0 = PS - pow(2, -LMD);
	mu = (int)ceil(LMD + 2 + log2((2.0 * (kap + 1) + l) * NN / p0));
	dlap_bias = (1 - p) / (1 + p - 2 * exp(-(pow(2, kap) + 1) / t));
	double k1 = NN / p0;
	double k2 = ((LMD + 2) * log(2)) / (2 * p0 * p0);
	m = (int)ceil(k1 + k2 / 2 + sqrt(k2 * k2 / 4 + k1 * k2));
	
    double geo_bias[kap], Bexp_bias[l];
    for(int i = 0; i < kap; ++i) geo_bias[i] = pow(p, (1 << i)) / (1 + pow(p, (1 << i)));
    for(int i = 0; i < l; ++i) Bexp_bias[i] = exp(-(1 << i) / r);

    p_lap = (uint32_t*)malloc((mu * NUM_SHARES) * sizeof(uint32_t));
    p_geo = (uint32_t*)malloc((kap * mu * NUM_SHARES) * sizeof(uint32_t));
    p_exp = (uint32_t*)malloc((l * mu * NUM_SHARES) * sizeof(uint32_t));

    bin_rep(p_lap, dlap_bias);
    bin_rep_arr(p_geo, kap, geo_bias);
    bin_rep_arr(p_exp, l, Bexp_bias);

    one[0] = 1;
    mSIG[0] = -round(SIG);
    Refresh(one);
    Refresh(mSIG);	
}


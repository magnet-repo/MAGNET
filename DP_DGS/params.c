#include <math.h>

#include "params.h"

int kap, l, mu, m, SIG;
double r, p, dlap_bias;

int SIGs[11] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};

void init_parameters(int j) {
	SIG = SIGs[j];
	double N0 = sqrt(2 * log(2) * (LMD + 2 + log2(N))) * SIG;
	kap = (int)ceil(log2(N0 - 1));
	l = 2 * kap;
	double t = (SIG * SIG) / round(SIG);
	r = 2 * SIG * SIG;
	p = exp(-1 / t);
	double p0 = PS - pow(2, -LMD);
	mu = (int)ceil(LMD + 2 + log2((2.0 * (kap + 1) + l) * N / p0));
	dlap_bias = (1 - p) / (1 + p - 2 * exp(-(pow(2, kap) + 1) / t));
	double k1 = N / p0;
	double k2 = ((LMD + 2) * log(2)) / (2 * p0 * p0);
	m = (int)ceil(k1 + k2 / 2 + sqrt(k2 * k2 / 4 + k1 * k2));
}


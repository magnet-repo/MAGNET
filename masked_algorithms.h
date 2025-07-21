#ifndef MASKED_ALGORITHMS_H_
#define MASKED_ALGORITHMS_H_

#include <stdint.h>
#include "params.h"

void MaskedBernoulli(uint32_t *x, uint32_t *b, uint32_t *one);
void MaskedGeometric(uint32_t *x, uint32_t *p_geo, uint32_t *one);
void MaskedLaplace(uint32_t *v, uint32_t *p_geo, uint32_t *p_lap, uint32_t *one);
void MaskedBernoulliExp(uint32_t *b, uint32_t *p_exp, uint32_t *u, uint32_t *one);
void MAGNET(uint32_t *samps, uint32_t *mSIG, uint32_t *one, uint32_t *p_geo,  uint32_t *p_exp, uint32_t *p_lap);

#endif /* MASKED_ALGORITHMS_H_ */

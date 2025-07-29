#ifndef ALGORITHMS_H_
#define ALGORITHMS_H_

#include <stdint.h>

uint32_t SamplingCircuit_B(uint32_t *b);
uint32_t SamplingCircuit_G(uint32_t *p_geo);
void SamplingCircuit_L(uint32_t *v, uint32_t *p_geo, uint32_t *p_lap);
void SamplingCircuit_Bexp(uint32_t *res, uint32_t *p_exp, uint32_t u);
void SamplingCircuit_N(uint32_t *samps, uint32_t *p_geo, uint32_t *p_exp, uint32_t *p_lap);

#endif /* ALGORITHMS_H_ */

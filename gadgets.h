#ifndef GADGETS_H_
#define GADGETS_H_

#include <stdint.h>

void Refresh(uint32_t* x);
void FullRefresh(uint32_t* z, const uint32_t* x);
uint32_t FullXOR(uint32_t* x);

uint32_t SecAND_Order1(uint32_t x, uint32_t y, uint32_t s, uint32_t t, uint32_t u);
void SecADD_Order1(uint32_t* z, const uint32_t* x, const uint32_t* y);
void SecINC_Order1(uint32_t *z, const uint32_t *x, const uint32_t *one);

void SecNEG(uint32_t* z, const uint32_t* x);
void SecAND(uint32_t* z, const uint32_t* x, const uint32_t* y);
void SecOR(uint32_t* z, const uint32_t* x, const uint32_t* y);
void SecINC(uint32_t *z, const uint32_t *x, const uint32_t *one);
void SecADD(uint32_t* z, const uint32_t* x, const uint32_t* y);
void SecSQU(uint32_t* C, const uint32_t* A);

void refreshBool(uint32_t a[], int n);
uint32_t Psi(uint32_t x, uint32_t y);
uint32_t Psi0(uint32_t x, uint32_t y, int n);
void copy(uint32_t *x, uint32_t *y, int n);
void SecB2A(uint32_t *D_, uint32_t *x, int n);
uint32_t GoubinAB(uint32_t A, uint32_t r, int k);
void refreshArith(uint32_t a[], int l, int n);
void bool2ArithSPOG(uint32_t *x, uint32_t *y, int k, int n);
void shift3(uint32_t *x, uint32_t *a, int k, int n);
void SecA2B(uint32_t *x, uint32_t *y, int k);

#endif /* GADGETS_H_ */

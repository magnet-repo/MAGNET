#ifndef BASE_GADGETS_H
#define BASE_GADGETS_H

#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <stdint.h>

#include "params-GR19.h"

typedef __int128_t int128_t;
typedef int masked[MASKING_ORDER+1];
typedef __int128_t masked128[MASKING_ORDER+1];

typedef int32_t poly[PARAM_N];
typedef poly masked_poly[MASKING_ORDER+1];
typedef int16_t masked_small_poly[MASKING_ORDER+1][PARAM_N];

void full_add(masked_poly, poly res);
void full_add_small(masked_poly, uint16_t* res);
int full_add_coef(int* masked);
int full_xor(int* x);
void refresh(int* x, int* res);
void refresh_vs(int* x, int* res, const int N_);
void refresh_masks_n(int* x, int* y, const int N_);
void expand(int* x, int* out, const int N_);
void full_refresh(int* x, int* res);
void full_refresh_arith(int* x, int* res);
void sec_and(int* x, int* y, int* res);
void sec_and_vs(int* x, int* y, int* res, const int N_);
void sec_add(int* x, int* y, int* z);
void sec_add_vs(int* x, int* y, int* z, const int N_);
void sec_and_const(int* x, int* y, int* res);
void sec_add_const(int* x, int* y, int* z);

void sec_add128(__int128_t* x, __int128_t* y, __int128_t* z);
void sec_and128(__int128_t* x, __int128_t* y, __int128_t* res);
void refresh128(__int128_t* x, __int128_t* res);

void sec_arith_bool_mod_p(int* a, int* a_prime);
void sec_bool_arith(int* x_bool, int* x_arith);
void convert_A_B(int* arith_x, int* bool_x, const int N_);

void goubin_bool_arith(int* bool_x, int* arith_x);
void goubin_arith_bool(int* arith_x, int* bool_x);
void order_1_add(int* x, int* y, int* z);
int order_1_AND(int x, int y, int s, int t, int u);
void HO_bool_arith(int* bool_x, int* arith_x, const int N_);
void seed_xoshiro_(void);
uint32_t xoshiro_next_(void);
#endif

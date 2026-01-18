#ifndef UTILS_H
#define UTILS_H

#include "params-GR19.h"
#include "sign_gadgets.h"
#include <stdio.h>

int get_CDT_SIZE(void);
void create_CDT(int64_t *cdt_v, int size);
int cmp_double_(const void *a, const void *b);
double calc_median_(double values[ITER_]);

int mod_q(int a);
int mod_q128(__int128_t a);
void print_bytes(unsigned char* b, int len);
void print_bits(int x);
void print_bits128(__int128_t x);
void print_shares_(int* x);
void print_shares_vs(int* x, const int N_);
void print_shares_bits(int* x);

void print_shares_bits_vs(int* x, const int N_);
void print_full_bits(int x);
void print_full_shares_bits(int* x);

void print_masked_poly(masked_poly p);
void print_poly_py(poly p);
void print_small_poly(int16_t* p);
void print_poly(poly p);
void print_small_masked_poly(masked_small_poly p);
void print_poly_f(poly p);

#endif

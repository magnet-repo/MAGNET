#ifndef SIGN_GADGETS_H
#define SIGN_GADGETS_H

#include "base_gadgets.h"
#include "params.h"

void gaussian(int* a, int64_t *cdt_v, int size);
void masked_gaussian_poly(masked_poly s, int64_t *cdt_v, int size);
void masked_sign_choice(masked_poly p);

#endif

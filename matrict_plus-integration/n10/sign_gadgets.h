#ifndef SIGN_GADGETS_H
#define SIGN_GADGETS_H

#include "base_gadgets.h"
#include "params-GR19.h"

void gaussian(int* a, int size);
void masked_gaussian_poly(masked_poly s, int size);
void masked_sign_choice(masked_poly p);

#endif

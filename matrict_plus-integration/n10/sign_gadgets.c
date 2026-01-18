#include "sign_gadgets.h"
#include "base_gadgets.h"

#include "utils-GR19.h"
#include "params-GR19.h"

#include "CDT_6391.h"

void gaussian(int *a, int size) {
	int i,j;
	masked128 r, delta;
	masked b, b_p;
	for (i=0; i < N_SHARES; ++i){
		r[i]   = rand();
		r[i] <<= 32;
		r[i]  ^= rand();
	}
	masked x = {0};

	for(j=1; j < size; ++j) {
		//masked128 k = {(int128_t)-cdt_v[j]};
		masked128 k = {(int128_t)-CDT_v[j]};
		
		sec_add128(r,k, delta);

		for(i=0; i < N_SHARES; ++i) {
			b[i] = (int)(delta[i] >> 127);
		}
		b[0] = ~b[0];
		masked J = {j-1};
		b_p[0] = ~b[0];
		for(i=1; i < N_SHARES; ++i) b_p[i] = b[i];

		sec_and(J, b, b);
		sec_and(x, b_p, b_p);
		for(i=0; i < N_SHARES; ++i) x[i] = b_p[i] ^ b[i];
	}
	for(i=0; i < N_SHARES; ++i) a[i] = x[i];
}

void masked_gaussian_poly(masked_poly s, int size) {
  masked v;
  int i,j;
  for(j=0; j < PARAM_N; ++j){
    gaussian(v, size);
    for(i=0; i < N_SHARES; ++i) s[i][j] = v[i];
  }
}

void masked_sign_choice(masked_poly p) {
  /*
    Assign a random sign to each coefficient of the masked polynomial
    This function is needed because the table is use to sample only half of the gaussian
  */
  int i,j;
  masked mask, sign, temp;

  for(j=0; j < PARAM_N; ++j){
    for(i=0; i < N_SHARES; ++i) temp[i] = p[i][j];
    for(i=0; i < N_SHARES; ++i) sign[i] = (rand() & 1); 	// 0 or 1
    for(i=0; i < N_SHARES; ++i) mask[i] = -sign[i]; 		// 0x00000000 if sign=0, 0x11111111 if sign=1
    for(i=0; i < N_SHARES; ++i) temp[i] ^= mask[i];
    sec_add(temp, sign, temp);
    for(i=0; i < N_SHARES; ++i) p[i][j] = temp[i];
  }
}


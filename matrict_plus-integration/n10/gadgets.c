#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

#include "params.h"
#include "gadgets.h"
#include "utils.h"     

void Refresh(uint32_t* x) { // [ISW03]
	for(int i = 0; i < NUM_SHARES; ++i) {
		for(int j = i + 1; j < NUM_SHARES; ++j) {
			uint32_t r = rand_uint32();
			x[i] ^= r;
			x[j] ^= r;
		}
	}
}

void FullRefresh(uint32_t* z, const uint32_t* x) { // [Cor14]
	for(int i = 0; i < NUM_SHARES; ++i) z[i] = x[i];
	for(int i = 0; i < NUM_SHARES; ++i) {
		for(int j = 1; j < NUM_SHARES; ++j) {
			uint32_t r = rand_uint32();
			z[i] ^= r;
			z[j] ^= r;
		}
	}
}

uint32_t FullXOR(uint32_t* x) { // [CGV14]
	uint32_t z[NUM_SHARES] = {0};
	FullRefresh(z, x);
	uint32_t val = z[0];
	for(int i = 1; i < NUM_SHARES; ++i) val ^= z[i];

	return val;
}

void SecNEG(uint32_t* z, const uint32_t* x) {
	for(int j = 0; j < NUM_SHARES; ++j) z[j] = -(x[j] & 1U);
    Refresh(z);
}

void SecAND(uint32_t* z, const uint32_t* x, const uint32_t* y) { // [BBE+18]
#if MASK_ORDER == 1
	int u = rand_uint32();
	int t;
	t = u ^ (x[0] & y[0]);
	t = t ^ (x[0] & y[1]);
	t = t ^ (x[1] & y[0]);
	t = t ^ (x[1] & y[1]);
	z[0] = t;
	z[1] = u;
#else
	for(int i = 0; i < NUM_SHARES; ++i) z[i] = x[i] & y[i];

	for(int i = 0; i < NUM_SHARES; ++i) {
		for(int j = i + 1; j < NUM_SHARES; ++j) {
			uint32_t tmp = rand_uint32();
			uint32_t tmp2 = (tmp ^ (x[i] & y[j])) ^ (x[j] & y[i]);
			z[i] ^= tmp;
			z[j] ^= tmp2;
		}
	}
#endif
}

void SecOR(uint32_t* z, const uint32_t* x, const uint32_t* y) { // [CC24]
	uint32_t t1[NUM_SHARES] = {0}, t2[NUM_SHARES] = {0}, t3[NUM_SHARES] = {0};
	t1[0] = ~x[0];
	for(int i = 1; i < NUM_SHARES; ++i) t1[i] = x[i];
	
	t2[0] = ~y[0];
	for(int i = 1; i < NUM_SHARES; ++i) t2[i] = y[i];
	SecAND(t3, t1, t2);
	
	z[0] = ~t3[0];
	for(int i = 1; i < NUM_SHARES; ++i) z[i] = t3[i];
}

void SecINC_Order1(uint32_t *z, const uint32_t *x, const uint32_t *one) {
    uint32_t g0 = x[0];
    uint32_t g1 = x[1];
    
    for(int j = 1; j <= WW; ++j) {
        uint32_t pw = 1 << (j - 1);
        uint32_t high = 0xFFFFFFFF << pw;
        uint32_t low  = ~high;

        uint32_t a0 = g0 << pw;
        uint32_t a1 = g1 << pw;
  
      	uint32_t r = rand_uint32(); 
        
        uint32_t a_prm0 = SecAND_Order1(a0, g0, a1, g1, r);
        uint32_t a_prm1 = r;
        
        g0 = (a_prm0 & high) ^ (g0 & low);
        g1 = (a_prm1 & high) ^ (g1 & low); 
    }
    z[0] = x[0] ^ one[0] ^ (g0 << 1);
    z[1] = x[1] ^ one[1] ^ (g1 << 1);
}

void SecINC(uint32_t *z, const uint32_t *x, const uint32_t *one) {
	#if MASK_ORDER == 1
		SecINC_Order1(z, x, one);
	#else
	    uint32_t g[NUM_SHARES] = {0}, a[NUM_SHARES] = {0}, a_prm[NUM_SHARES] = {0};
        for(int i = 0; i < NUM_SHARES; ++i) g[i] = x[i];

        for(int j = 1; j <= WW; ++j) {
        	uint32_t pw = 1 << (j - 1);
            const uint32_t high = 0xFFFFFFFF << pw;
            const uint32_t low  = ~high;
            
            for(int i = 0; i < NUM_SHARES; ++i) a[i] = g[i] << pw;
            Refresh(a);
            
            SecAND(a_prm, a, g);

            for(int i = 0; i < NUM_SHARES; ++i) g[i] = (a_prm[i] & high) ^ (g[i] & low);
        }
        for(int i = 0; i < NUM_SHARES; ++i) z[i] = x[i] ^ one[i] ^ (g[i] << 1);
    #endif
}


static uint32_t SecXOR_Order1(uint32_t x, uint32_t y, uint32_t u) { // [CGVT15]
	return x ^ y ^ u;
}

static uint32_t SecSHIFT_Order1(uint32_t x, uint32_t s, uint32_t t, uint32_t j) { // [CGVT15]
	uint32_t y;
	y = t ^ (x << j);
	y = y ^ (s << j);
	return y;
}
uint32_t SecAND_Order1(uint32_t x, uint32_t y, uint32_t s, uint32_t t, uint32_t u) { // [CGVT15]
	uint32_t z;
	z = u ^ (x & y);
	z = z ^ (x & t);
	z = z ^ (s & y);
	z = z ^ (s & t);
	return z;
}

void SecADD_Order1(uint32_t* z, const uint32_t* x, const uint32_t* y) {  // [CGVT15]
	int t, u, x_p, y_p, z_p, r, s, P, G, H, U;
	x_p = x[0]; s = x[1];
	y_p = y[0]; r = y[1];
	t = rand_uint32(); u = rand_uint32();

	P = SecXOR_Order1(x_p, y_p, r);
	G = SecAND_Order1(x_p, y_p, s, r, u);
	G = G ^ s;
	G = G ^ u;
	for(int i = 1; i <= WW-1; ++i){
		H = SecSHIFT_Order1(G, s, t, 1<<(i-1));
		U = SecAND_Order1(P,H,s,t,u);
		G = SecXOR_Order1(G, U, u);
		H = SecSHIFT_Order1(P, s, t, 1<<(i-1));
		P = SecAND_Order1(P, H, s, t, u);
		P = P ^ s;
		P = P ^ u;
	}

	H = SecSHIFT_Order1(G,s,t, (1<<(WW-1)));
	U = SecAND_Order1(P,H,s,t,u);
	G = SecXOR_Order1(G, U, u);
	z_p = SecXOR_Order1(y_p, x_p, s);
	z_p = z_p ^ (2*G);
	z_p = z_p ^ (2*s);
	z[0] = z_p;
	z[1] = r;
}

void SecADD(uint32_t* z, const uint32_t* x, const uint32_t* y) { // [BBE+18]
	#if MASK_ORDER == 1
		SecADD_Order1(z, x, y);
	#else
		uint32_t p[NUM_SHARES] = {0}, p_[NUM_SHARES] = {0}, g[NUM_SHARES] = {0}, a[NUM_SHARES] = {0}, a_[NUM_SHARES] = {0}, ap[NUM_SHARES] = {0};
		for(int i = 0; i < NUM_SHARES; ++i) p[i] = x[i] ^ y[i];
		SecAND(g, x, y);
		
		for(int j = 1; j <= WW - 1; ++j) {
			uint32_t pw = 1 << (j - 1);
			for(int i = 0; i < NUM_SHARES; ++i) a[i] = g[i] << pw;
			SecAND(a_, a, p);
			
			for(int i = 0; i < NUM_SHARES; ++i) g[i] ^= a_[i];

			for(int i = 0; i < NUM_SHARES; ++i) ap[i] = p[i] << pw;
			Refresh(ap);
			SecAND(p_, p, ap);
			
			for(int i = 0; i < NUM_SHARES; ++i) p[i] = p_[i];
		}
		for(int i = 0; i < NUM_SHARES; ++i) a[i] = g[i] << (1 << (WW - 1));
		SecAND(a_, a, p);
		
		for(int i = 0; i < NUM_SHARES; ++i) g[i] ^= a_[i];

		for(int i = 0; i < NUM_SHARES; ++i) z[i] = x[i] ^ y[i] ^ (g[i] << 1);
	#endif
}

void SecSQU(uint32_t* C, const uint32_t* A) {
	for(int i = 0; i < NUM_SHARES; ++i) C[i] = ((uint64_t)A[i] * A[i]);

	for(int i = 0; i < NUM_SHARES; ++i) {
		for(int j = i + 1; j < NUM_SHARES; ++j) {
			uint32_t r = rand_uint32();
			uint32_t t1 = (uint64_t)A[i] * A[j];
			uint32_t t2 = r + (t1 << 1);
			C[i] -= r;
			C[j] += t2;
		}
	}
}

/* ---------------------------------------------------------------- */

void refreshBool(uint32_t a[], int n) { // auxiliary gadget used in SecB2A
	for(int i = 0; i < n - 1; i++) {
		uint32_t tmp = rand_uint32();
		a[n-1] = a[n-1] ^ tmp;
		a[i] = a[i] ^ tmp;
	}
}

uint32_t Psi(uint32_t x, uint32_t y) { // auxiliary gadget used in SecB2A
	return (x ^ y) - y;
}

uint32_t Psi0(uint32_t x, uint32_t y, int n) { // auxiliary gadget used in SecB2A
	return Psi(x, y) ^ ((~n & 1) * x);
}

void copy(uint32_t *x, uint32_t *y, int n) { // auxiliary gadget used in SecB2A
	for(int i = 0; i < n; i++) x[i] = y[i];
}

static void impconvBA_rec(uint32_t *D_, uint32_t *x, int n) { // auxiliary gadget used in SecB2A
	if(n == 2) { // // Goubin's first order conversion from Boolean to arithmetic [Gou01]
		uint32_t r1 = rand_uint32();
		uint32_t r2 = rand_uint32();
		uint32_t y0 = (x[0] ^ r1) ^ r2;
		uint32_t y1 = x[1] ^ r1;
		uint32_t y2 = x[2] ^ r2;

		uint32_t z0 = y0 ^ Psi(y0, y1);
		uint32_t z1 = Psi(y0, y2);

		D_[0] = y1 ^ y2;
		D_[1] = z0 ^ z1;

		return;
	}

	uint32_t y[n+1];
	copy(y, x, n+1);

	refreshBool(y, n+1);

	uint32_t z[n];
	z[0] = Psi0(y[0], y[1], n);
	
	for(int i = 1; i < n; i++) z[i] = Psi(y[0], y[i+1]);

	uint32_t A[n-1], B[n-1];
	impconvBA_rec(A, y+1, n-1);
	impconvBA_rec(B, z, n-1);

	for(int i = 0; i < n - 2; i++) D_[i] = A[i] + B[i];

	D_[n-2] = A[n-2];
	D_[n-1] = B[n-2];
}

void SecB2A(uint32_t *D_, uint32_t *x, int n) { //  Boolean-to-Arithmetic conversion algorithm of [BCZ18]
	uint32_t x_ext[n+1];
	copy(x_ext, x, n);
	x_ext[n] = 0;
	impconvBA_rec(D_, x_ext, n);
}

uint32_t GoubinAB(uint32_t A, uint32_t r, int k) { // Goubin's first order conversion from arithmetic to Boolean [Gou01]
	uint32_t G = rand_uint32();
	uint32_t T = G << 1;
	uint32_t x = G ^ r;
	uint32_t O = G & x;
	x = T ^ A;
	G = G ^ x;
	G = G & r;
	O = O ^ G;
	G = T & A;
	O = O ^ G;
	for(int i=1;i<k;i++) {
	    G = T & r;
	    G = G ^ O;
	    T = T & A;
	    G = G ^ T;
	    T = G << 1;
	}
	
	x = x ^ T;
	return x;
}

void refreshArith(uint32_t a[],int l,int n) { // auxiliary gadget used in SecA2B
	uint32_t tmp;
	uint32_t ma = (1 << l) - 1;
	for(int i = 0; i < n - 1; i++) {
		tmp = rand_uint32();
		a[n-1] = (a[n-1] + tmp) & ma;
		a[i] = (a[i] - tmp) & ma;
	}
}

void bool2ArithSPOG(uint32_t *x, uint32_t *y, int k, int n) { // auxiliary gadget used in SecA2B
	y[0] = x[0];
	for(int i = 1; i < n; i++) {
		y[i] = 0;
		refreshArith(y, k, i+1);

		for(int j = 0; j < (i+1); j++) y[j] = y[j] * (1 - (x[i] << 1));
		y[0] = y[0] + x[i];
	}
	refreshArith(y, k, n);
}

void shift3(uint32_t *x, uint32_t *a, int k, int n) { // auxiliary gadget used in SecA2B
	uint32_t b[n];
	for(int i = 0; i < n; i++) b[i] = x[i] & 1;

	uint32_t y[n];
	bool2ArithSPOG(b, y, k, n);

	uint32_t u[n];
	for(int i = 0; i < n; i++) u[i] = x[i] - y[i];

	for(int i = 0; i < n; i++) {
		u[n-1] += (u[i] & 1);
		u[i] -= (u[i] & 1);
	}
	for(int i = 0; i < n; i++) a[i] = u[i] >> 1;
}

void SecA2B(uint32_t *x, uint32_t *y, int k) { // Arithmetic-to-Boolean conversion algorithm of [CGTZ23]
	#if MASK_ORDER == 1
		y[0] = GoubinAB(x[0], x[1], k);
		y[1] = x[1];
	#else
	uint32_t a[NUM_SHARES];
	for(int i = 0; i < NUM_SHARES; i++) a[i] = x[i];

	for(int i = 0; i < NUM_SHARES; i++) y[i] = 0;

	for(int j = 0; j < k; j++) {
		for(int i = 0; i < NUM_SHARES; i++) y[i] += ((a[i] & 1) << j);
		    uint32_t a2[NUM_SHARES];
		    if(j < k - 1) {
			    shift3(a, a2, k-j, NUM_SHARES);
			    for(int i = 0; i < NUM_SHARES; i++) a[i] = a2[i];
		    }
	}
	#endif
}


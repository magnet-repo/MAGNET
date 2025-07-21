#include <math.h>
#include "params.h"
#include "utils.h"

int get_CDT_SIZE() {
    const long double PI  = acosl(-1.0L);
    const long double inv = 1.0L / (SIG * sqrtl(2.0L * PI));

    long double cdf = inv;
    int n = 2;

    for (int i = 1; ; ++i) {
        long double w = inv * expl(-(long double)(i * i) / (2.0L * SIG * SIG));
        if (w * SCALE < 1.0L) {
        	++n;
        	break;
        }
        cdf += 2.0L * w;
        ++n;
    }
    return n;
}

void create_CDT(int64_t *cdt_v, int size) {
	const long double PI   = acosl(-1.0L);
	const long double inv  = 1.0L / (SIG * sqrtl(2.0L * PI));
	cdt_v[0] = 0;

	long double cdf = inv;
	cdt_v[1] = (int64_t)(cdf * SCALE + 0.5L);

    int pos = 2;
	for (int i = 1; pos < size - 1; ++i) {
		long double w = inv * expl(-(long double)(i * i) / (2.0L * SIG * SIG));
		long double next_cdf = cdf + 2.0L * w;
		cdt_v[pos++] = (int64_t)(next_cdf * SCALE + 0.5L);
		cdf = next_cdf;
	}
	cdt_v[pos] = (int64_t)SCALE;
}

int cmp_double(const void *a, const void *b) {
    double da = *(const double*)a;
    double db = *(const double*)b;
    if (da < db) return -1;
    if (da > db) return  1;
    return 0;
}

double calc_median(double values[ITER]) {
    qsort(values, ITER, sizeof(values[0]), cmp_double);

    double mid1 = values[ITER/2 - 1];
    double mid2 = values[ITER/2];
    return (mid1 + mid2) / 2.0;
}

/* ---------------------------------- */

void print_poly(poly p){
    int i;
    for(i=0; i < 10; ++i)
        printf("%i ",p[i]);
    printf("\n");
}

void print_small_poly(int16_t* p){
    int i;
    for(i=0; i < 10; ++i)
        printf("%i ",p[i]);
    printf("\n");
}

void print_poly_py(poly p){
    int i;
    printf("[");
    for(i=0; i < PARAM_N-1; ++i)
        printf("%i, ",p[i]);
    printf("%i]\n",p[PARAM_N-1]);    
}

void print_poly_f(poly p){
    int i;
    for(i=0; i < PARAM_N; ++i)
        printf("%i ",p[i]);
    printf("\n");
}

void print_small_masked_poly(masked_small_poly p){
    int i,j;
    int val=0;
    for(i=0; i < 10; ++i){
        for(j=0; j <= MASKING_ORDER; ++ j){
            printf("%i |",p[j][i]);
            val += p[j][i];
        }
        printf(" = %i = %i\n", val, val&(PARAM_Q-1));
        val = 0;
    }
}

void print_masked_poly(masked_poly p){
    int i,j;
    int val=0;
    for(i=0; i < 30; ++i){
        for(j=0; j <= MASKING_ORDER; ++ j){
            printf("%i |",p[j][i]);
            val += p[j][i];
        }
        printf(" = %i = %i\n", val, val&(PARAM_Q-1));
        val = 0;
    }
}

int mod_q(int a){
    return a&(PARAM_Q-1);
}

int mod_q128(__int128_t a){
    return a&(PARAM_Q-1);
}


void print_bytes(unsigned char* b, int len){
  for(int i = 0; i < len; ++i)
    printf("%X ", b[i]);
  printf("\n");

}

void print_bits(int x){
    //int i;
    //for(i=sizeof(int)*8-1; i >= 0 ; --i) printf("%i",(x>>i)&1);
    //for(i=26; i >= 0 ; --i) printf("%i",(x>>i)&1);
    printf("0x%X",x);
}

void print_bits128(__int128_t x){
    printf("0x%08X", (int32_t)(x>>96));
    printf("%08X", (int32_t)(x>>64));
    printf("%08X", (int32_t)(x>>32));
    printf("%X", (int32_t)x);
}

void print_shares(int* x){
    int i;
    int val=0;
    for(i=0; i < N_SHARES-1; ++i){
        printf("%i | ", x[i]);
        val += x[i];
    }
    printf("%i = ", x[N_SHARES-1]);
    val += x[N_SHARES-1];
    printf("%i = ", val);
    
    printf("%i\n", mod_q(val));
}


void print_shares_vs(int* x, const int N){
    int i;
    int val=0;
    for(i=0; i < N-1; ++i){
        printf("%i | ", x[i]);
        val += x[i];
    }
    printf("%i = ", x[N-1]);
    val += x[N-1];
    printf("%i = ", val);
    
    printf("%i\n", mod_q(val));
}

void print_shares_bits(int* x){
    int i;
    int val=0;
    for(i=0; i < MASKING_ORDER; ++i){
        print_bits(x[i]); printf(" | ");
        val ^= x[i];
    }
    print_bits(x[MASKING_ORDER]);
    val ^= x[MASKING_ORDER];
    printf(" = ");
    print_bits(val);
    printf(" = %i = %i = ",val, mod_q(val));
    print_bits(mod_q(val));
    printf("\n");
}

void print_shares_bits_vs(int* x, int N){
    int i;
    int val=0;
    for(i=0; i < N-1; ++i){
        print_bits(x[i]); printf(" | ");
        val ^= x[i];
    }
    print_bits(x[N-1]);
    val ^= x[N-1];
    printf(" = ");
    print_bits(val);
    printf(" = %i = %i = ",val, mod_q(val));
    print_bits(mod_q(val));
    printf("\n");
}

void print_full_bits(int x){
    int i;
    for(i=sizeof(int)*8-1; i >= 0 ; --i) printf("%i",(x>>i)&1);
    //for(i=20; i >= 0 ; --i) printf("%i",(x>>i)&1);
}


void print_full_shares_bits(int* x){
    int i;
    int val=0;
    for(i=0; i < MASKING_ORDER; ++i){
        print_full_bits(x[i]); printf(" | ");
        val ^= x[i];
    }
    print_full_bits(x[MASKING_ORDER]);
    val ^= x[MASKING_ORDER];
    printf(" = ");
    print_full_bits(val);
    printf(" = %i\n",val);
}


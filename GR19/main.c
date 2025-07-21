#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#include "params.h"
#include "utils.h"
#include "base_gadgets.h"
#include "sign_gadgets.h"

//#define VERIFY

int main(void) {
	seed_xoshiro();
	int res[PARAM_N] = {0};
	
    printf("\nRunning with Masking Order-%d, N = %d, lambda = 64, varying SIGMAs [2^0,...,2^10]:\n\n", MASKING_ORDER, PARAM_N);	
	for(int ct = 0; ct < 11; ct++) {
        init_parameters(ct);
        printf("SIGMA = %d\n", SIG);
	    
	    int CDT_SIZE = get_CDT_SIZE();

        int64_t *cdt_v = (int64_t*)malloc(CDT_SIZE * sizeof(int64_t));
        create_CDT(cdt_v, CDT_SIZE);

        clock_t start_t, end_t;
        double runs[ITER];
        
        for(int k = 0; k < ITER; ++k) {
        
            start_t = clock();
            
            masked_poly t;

            masked_gaussian_poly(t, cdt_v, CDT_SIZE);
            masked_sign_choice(t);

            int i, j;
            masked temp;
            for(i = 0; i < PARAM_N; i++) {
               for(j = 0; j < N_SHARES; j++) {
            	   temp[j] = t[j][i];
               }
               res[i] = full_xor(temp);
            }    
               
            end_t = clock();
            
            runs[k] = (double)(end_t - start_t) / CLOCKS_PER_SEC;
            
            #ifdef VERIFY    
                FILE *outfile = fopen("samples.txt", "w");

	            for (int i = 0; i < PARAM_N; i++) {
		            fprintf(outfile, "%d ", res[i]);
	            }
	            fprintf(outfile, "%d ", (int)(SIG * SIG));
	            fprintf(outfile, "%d ", PARAM_N);

	            fclose(outfile);

	            exit(system("python3 discretegauss.py")); // plotting samples using "The Discrete Gaussian for Differential Privacy [CKS20]"
            #endif   
        }
        double med = calc_median(runs);
        printf("Median of %d runs: %.2f s\n\n", ITER, med);
        
    }
 
	return 0;
}

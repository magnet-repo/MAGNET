#include <math.h>

#include "params-GR19.h"
#include "utils-GR19.h"

int SIG_, CDT_SIZE;
int64_t *cdt_v = NULL;

void init_parameters_() {
    SIG_ = 6390; 
    
    //CDT_SIZE = get_CDT_SIZE();
    //cdt_v = (int64_t*)malloc(CDT_SIZE * sizeof(int64_t));
    //create_CDT(cdt_v, CDT_SIZE);
}

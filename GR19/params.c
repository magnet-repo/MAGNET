#include <math.h>

#include "params.h"

int SIG;

int SIGs[11] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};

void init_parameters(int ct) {
    SIG = SIGs[ct];
}

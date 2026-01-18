#include "hal.h"
#include "simpleserial.h"
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define NUM_SHARES 2
#define MASK_ORDER NUM_SHARES - 1
#define w 32
#define W 5 // ceil(log2(w-1))

#define PARAMS_LOGQ 16
#define PARAMS_Q (1 << PARAMS_LOGQ)

#define N 10
#define SIG 10
#define LMD 16
#define PS 0.76

/* ---------------------------------------------------------------------- */

uint8_t get_pt(uint8_t* pt, uint8_t len) {
    volatile uint32_t zp[NUM_SHARES], xp[NUM_SHARES], yp[NUM_SHARES], up[1];
    uint8_t res[4];
    
    xp[0] = ((uint32_t)pt[3] << 24) | ((uint32_t)pt[2] << 16) | ((uint32_t)pt[1] << 8) | pt[0];
    xp[1] = ((uint32_t)pt[7] << 24) | ((uint32_t)pt[6] << 16) | ((uint32_t)pt[5] << 8) | pt[4];
    
    yp[0] = ((uint32_t)pt[11] << 24) | ((uint32_t)pt[10] << 16) | ((uint32_t)pt[9] << 8) | pt[8];
    yp[1] = ((uint32_t)pt[15] << 24) | ((uint32_t)pt[14] << 16) | ((uint32_t)pt[13] << 8) | pt[12];

    up[0] = ((uint32_t)pt[19] << 24) | ((uint32_t)pt[18] << 16) | ((uint32_t)pt[17] << 8) | pt[16];

    for (volatile int k = 0; k < 1000; k++) {;} // to clean the power trace

    trigger_high();

    SecAND_m4(zp, xp, yp, up);

    trigger_low();

    for (volatile int k = 0; k < 100; k++) {;}  // to clean the power trace

    uint32_t ans = (zp[0] ^ zp[1]);

    res[0] = ans & 0xff;
    res[1] = (ans & 0xffff) >> 8;
    res[2] = (ans & 0xffffff) >> 16;
    res[3] = ans >> 24;
    
    simpleserial_put('r', 4, res);

    return 0x00;
}

int main(void) {
    platform_init();
    init_uart();
    trigger_setup();

    simpleserial_init();

    simpleserial_addcmd('p', 20,  get_pt);
        
    while(1)
        simpleserial_get();
}


#include "hal.h"
#include "simpleserial.h"
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>

uint8_t get_pt(uint8_t* pt, uint8_t len) {

    volatile uint32_t zp[NUM_SHARES], xp[NUM_SHARES], yp[NUM_SHARES], up[1];
    uint8_t res[4];
    
    xp[0] = ((uint32_t)pt[3] << 24) | ((uint32_t)pt[2] << 16) | ((uint32_t)pt[1] << 8) | pt[0];
    xp[1] = ((uint32_t)pt[7] << 24) | ((uint32_t)pt[6] << 16) | ((uint32_t)pt[5] << 8) | pt[4];
    
    up[0] = ((uint32_t)pt[11] << 24) | ((uint32_t)pt[10] << 16) | ((uint32_t)pt[9] << 8) | pt[8];

    for (volatile int k = 0; k < 1000; k++) {;} // to clean the power trace

    trigger_high();

    SecBA_m4(zp, xp, up);

    trigger_low();

    for (volatile int k = 0; k < 100; k++) {;}  // to clean the power trace

    uint32_t ans = (zp[0] + zp[1]);

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

    simpleserial_addcmd('p', 12,  get_pt);
        
    while(1)
        simpleserial_get();
}


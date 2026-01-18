.section .text
.syntax unified
.cpu cortex-m4
.thumb

.global SecAND_m4
SecAND_m4:
    PUSH    {r4-r11, lr}

    @ Register aliases
    pz      .req r0      @ Output pointer
    py      .req r1      @ Input Y pointer
    px      .req r2      @ Input X pointer
    pool    .req r3      @ Scratch memory / pool pointer
    t1      .req r4
    rx0     .req r5
    rx1     .req r12
    rz0     .req r7
    rz1     .req r8
    t2      .req r9
    ry0     .req r10
    ry1     .req r11
    r       .req r6

    LDR     rx0, [px], #4
    LDR     r,   [pool]
    LDR     rx1, [px], #4
    NOP

    LDR     ry0, [py], #4
    NOP
    STR     t1,  [pool]	
    NOP
    LDR     ry1, [py], #4
    NOP

    AND     rz0, rx0, ry0
    EOR     t1,  t1, t1
    AND     rz1, rx1, ry1

    EOR     px,  px, px
    AND     t1,  rx0, ry1

    EOR     px,  px, px
    AND     t2,  rx1, ry0

    NOP
    EOR     px,  px, px
    EOR     t1,  r, t1

    NOP                  
    EOR     px,  px, px 
    EOR     t1,  t2, t1

    EOR     px,  px, px
    EOR     rz0, rz0, r

    NOP
    EOR     px,  px, px
    EOR     rz1, rz1, t1

    EOR     px,  px, px
    NOP
    STR     rz0, [pz], #4
    NOP
    STR     px,  [pool]
    NOP
    STR     rz1, [pz], #4
          
    @ ---- Clean up register aliases ----
    .unreq  pz
    .unreq  py
    .unreq  px
    .unreq  pool
    .unreq  t1
    .unreq  rx0
    .unreq  rx1
    .unreq  rz0
    .unreq  rz1
    .unreq  t2
    .unreq  ry0
    .unreq  ry1
    .unreq  r

    POP     {r4-r11, pc}

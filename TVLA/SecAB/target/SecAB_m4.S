.section .text
.syntax unified
.cpu cortex-m4
.thumb

.global SecAB_m4
SecAB_m4:

    PUSH {r4-r11,r14}

    py          .req r0
    px          .req r1
    pool        .req r2
    T           .req r3
    rx0         .req r4
    rx1         .req r5
    ry0         .req r6
    ry1         .req r7
    G           .req r8
    O           .req r9
    t1          .req r10
    t2          .req r11
    t3          .req r12

    LDR     rx0, [px], #4       @ Load x[0]
    
    LDR     G,   [pool], #4     @ Load Random
    
    LDR     rx1, [px], #4       @ Load x[1]

    LSL     T, G, #1            @ T = 2*G
    AND     O, G, rx1           @ O = G & x1
    EOR     O, O, G             @ O = (G & x1) ^ G = G & (~x1)

    NOP
    NOP

    EOR     ry0, T, rx0
    AND     t2, T, rx0
    EOR     G, G, ry0

    NOP
    EOR     t1, t1, t1  

    AND     G, G, rx1

    NOP
    EOR     t1, t1, t1  
 
    EOR     O, O, G

    NOP
    EOR     t1, t1, t1    
    
    EOR     O, O, t2

    NOP
    EOR     t2, t2, t2

    .rept 31
        NOP
        EOR     t2, t2, t2    
        
        NOP
        EOR     t3, t3, t3
    
        AND     G, T, rx1

        EOR     G, G, O
        AND     T, T, rx0
        EOR     G, G, T
        
        NOP
        EOR     t1, t1, t1        
        
        LSL     T, G, #1
    .endr

    NOP
    EOR     t3, t3, t3
    EOR     ry0, ry0, T

    NOP
    EOR     px, px, px
    NOP
    STR     t2, [pool]          @ Dummy Store
    NOP
    STR     ry0, [py], #4       @ Store y[0]
    NOP
    STR     t2, [pool]          @ Dummy Store
    NOP
    STR     rx1, [py], #4       @ Store y[1]

    .unreq py
    .unreq px
    .unreq pool
    .unreq T    
    .unreq rx0    
    .unreq rx1 
    .unreq ry0 
    .unreq ry1  
    .unreq G    
    .unreq O   
    .unreq t1  
    .unreq t2
    .unreq t3

    POP {r4-r11,pc}

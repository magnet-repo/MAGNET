.section .text
.syntax unified
.cpu cortex-m4
.thumb

.global SecBA_m4
SecBA_m4:

    PUSH {r4-r11,r14}

    @ Register aliases
    px          .req r1
    py          .req r0
    pool        .req r2
    T           .req r3
    rx0         .req r4
    rx1         .req r5
    ry0         .req r6
    ry1         .req r7
    G           .req r8
    t1          .req r9
    t2          .req r10
    t3          .req r11
    t4          .req r12

    @ 1. Process Share 1 (rx1) and Gamma
    LDR G, [pool]           @ Load 32-bit Random (Gamma)
    LDR rx1, [px, #4]       @ Load 32-bit input share x1

    EOR G, G, rx1           @ G = Gamma ^ x1
    STR rx1, [py, #4]       @ Store output share y1 (unchanged x1)
    
    EOR rx1, t4, t4         @ Scrub rx1 register

    @ 2. Process Share 0 (rx0)
    LDR t1, [pool]          @ Reload Gamma
    LDR rx0, [px]           @ Load 32-bit input share x0

    EOR T, rx0, t1          @ T = x0 ^ Gamma
    
    SUB T, T, t1            @ T = (x0 ^ Gamma) - Gamma (Mod 2^32 implicit)

    EOR T, T, rx0           @ T = ((x0 ^ Gamma) - Gamma) ^ x0
    
    @ 3. Cross Combination
    EOR ry0, rx0, G         @ ry0 = x0 ^ (Gamma ^ x1) = x ^ Gamma
    EOR rx0, t4, t4         @ Scrub rx0 register
    
    SUB ry0, ry0, G         @ ry0 = (x ^ Gamma) - (Gamma ^ x1)
    
    @ 4. Final Combination
    EOR t3, T, ry0          @ Final Merge
    STR t3, [py]            @ Store output share y0

    @ Cleanup aliases
    .unreq py            
    .unreq px            
    .unreq pool     
    .unreq T            
    .unreq rx0             
    .unreq rx1          
    .unreq ry0          
    .unreq ry1             
    .unreq G             
    .unreq t1            
    .unreq t2             
    .unreq t3            
    .unreq t4            

    POP {r4-r11,pc}

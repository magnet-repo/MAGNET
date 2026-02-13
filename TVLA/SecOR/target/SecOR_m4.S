.section .text
.syntax unified
.cpu cortex-m4
.thumb

.global SecOR_m4
.extern SecAND_m4

SecOR_m4:
    PUSH    {r4-r11, lr}
    
    @ Arguments: r0=z, r1=x, r2=y, r3=pool
    MOV     r4, r0              @ r4 = z
    MOV     r5, r1              @ r5 = x
    MOV     r6, r2              @ r6 = y
    MOV     r7, r3              @ r7 = pool
    
    @ -------------------------------------------------------------------------
    @ STACK ALLOCATION
    @ -------------------------------------------------------------------------
    SUB     sp, sp, #32
    MOV     r8, sp
    
    @ Initialize Scrub Slot
    LDR     r12, [r7]
    STR     r12, [r8, #24]

    @ -------------------------------------------------------------------------
    @ 1. PREPARE INPUTS (~x)
    @ -------------------------------------------------------------------------
    
    @ --- Share 0 ---
    LDR     r0, [r5]            @ Load x[0]
    LDR     r12, [r7]
    MVN     r0, r0              @ ~x[0]
    STR     r0, [r8]            @ Store t1[0]  <-- Write Bus: ~x0

    @ SCRUB STACK WRITE BUS
    LDR     r12, [r7]           @ Load Random
    STR     r12, [r8, #24]      @ Scrub Stack address (sp+24)

    LDR     r12, [r7]

    @ --- Share 1 ---
    LDR     r1, [r5, #4]        @ Load x[1]
    LDR     r12, [r7]
    STR     r1, [r8, #4]        @ Store t1[1]

    LDR     r12, [r7]

    @ -------------------------------------------------------------------------
    @ 2. PREPARE INPUTS (~y)
    @ -------------------------------------------------------------------------
    
    @ --- Share 0 ---
    LDR     r0, [r6]            @ Load y[0]
    MVN     r0, r0              @ ~y[0]
    STR     r0, [r8, #8]        @ Store t2[0]

    @ SCRUB STACK WRITE BUS
    LDR     r12, [r7]           @ Load Random
    STR     r12, [r8, #24]      @ Scrub Stack address (sp+24)

    LDR     r12, [r7]

    @ --- Share 1 ---
    LDR     r1, [r6, #4]        @ Load y[1]
    LDR     r12, [r7]
    STR     r1, [r8, #12]       @ Store t2[1]

    @ -------------------------------------------------------------------------
    @ 3. CALL secAND(t3, t1, t2)
    @ -------------------------------------------------------------------------
    
    LDR     r12, [r7]
    
    ADD     r0, r8, #16         @ r0 = t3 (Output)
    MOV     r2, r8              @ r2 = t1 (Input X)
    ADD     r1, r8, #8          @ r1 = t2 (Input Y)
    
    @ CRITICAL: PASS REAL POOL (r7)
    @ We rely on the modified secAND (below) to scrub the read bus.
    MOV     r3, r7              @ r3 = pool

    LDR     r12, [r7]

    BL      SecAND_m4
    
    @ -------------------------------------------------------------------------
    @ 4. FINALIZE (z = ~t3)
    @ -------------------------------------------------------------------------

    LDR     r12, [r7]

    @ --- Process Share 1 ---
    LDR     r1, [r8, #20]       @ Load t3[1]
    STR     r1, [r4, #4]        @ Store z[1]   <-- Output Bus: z[1]

    @ SCRUB OUTPUT WRITE BUS
    LDR     r12, [r7]           @ Load Random
    STR     r12, [r4]           @ Scrub Target Address (r4) directly

    @ --- Process Share 0 ---
    EOR     r0, r0, r0
    EOR     r1, r1, r1    
    
    LDR     r0, [r8, #16]       @ Load t3[0]
 
    LDR     r12, [r7]
 
    MVN     r0, r0              @ ~t3[0]

    STR     r12, [r4]           @ Scrub Target Address (r4) directly

    LDR     r12, [r7]

    STR     r0, [r4]            @ Store z[0]
    
    @ Cleanup
    ADD     sp, sp, #32
    POP     {r4-r11, pc}

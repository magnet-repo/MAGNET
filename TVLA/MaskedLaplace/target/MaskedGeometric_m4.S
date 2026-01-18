.section .text
.syntax unified
.cpu cortex-m4
.thumb

.global MaskedGeometric_m4

@ -----------------------------------------------------------------------------
@ Constants
@ -----------------------------------------------------------------------------
.equ KAP, 2
.equ MU, 2
.equ P_GEO_STEP, 16        @ mu * NUM_SHARES * 4 = 16 bytes
.equ RS_STEP, 16           @ mu * NUM_SHARES * 4 = 16 bytes
.equ POOL_STEP_BERN, 8     @ MaskedBernoulli consumes 8 bytes (2 words)

@ -----------------------------------------------------------------------------
@ MaskedGeometric_m4
@ Parameters:
@   R0 = x        (uint32_t*) Output state (Accumulator)
@   R1 = p_geo    (uint32_t*) Geometric probabilities
@   R2 = one      (uint32_t*) Constant one
@   R3 = rs       (uint32_t*) Randomness for p
@   [SP] = pool   (uint32_t*) Fresh randomness stream
@ -----------------------------------------------------------------------------
MaskedGeometric_m4:
    PUSH    {r4-r11, lr}

    @ --- 1. Load Arguments & Setup Registers ---
    LDR     r7, [sp, #36]       @ R7 = pool
    
    MOV     r4, r0              @ R4 = x
    MOV     r5, r1              @ R5 = p_geo
    MOV     r6, r3              @ R6 = rs
    MOV     r10, r2             @ R10 = one

    @ --- 2. Stack Layout ---
    @ We only need space for MaskedBernoulli output 'b' (2 words)
    SUB     sp, sp, #32         @ 32 bytes for safety/alignment
    MOV     r9, sp              @ R9 = Stack Base

    @ --- 3. Loop Initialization ---
    MOV     r8, #0              @ j = 0

loop_geo:
    CMP     r8, #KAP
    BGE     end_geo

    @ -------------------------------------------------------------------------
    @ Step A: MaskedBernoulli
    @ -------------------------------------------------------------------------
    ADD     r0, r9, #4          @ R0 = &b (Output to Stack+4)
    MOV     r1, r5              @ R1 = p_geo
    MOV     r2, r10             @ R2 = one
    MOV     r3, r6              @ R3 = rs
    STR     r7, [sp]            @ [SP] = pool

    BL      MaskedBernoulli_m4

    @ Update Pointers
    ADD     r5, r5, #P_GEO_STEP
    ADD     r6, r6, #RS_STEP
    ADD     r7, r7, #POOL_STEP_BERN

    @ -------------------------------------------------------------------------
    @ Step B & C Merged: Shift & XOR (Accumulate)
    @ x[i] = x[i] ^ (b[i] << j)
    @ -------------------------------------------------------------------------
    
    @ --- Process Share 0 ---
    LDR     r0, [r9, #4]        @ Load b[0]
    LSL     r0, r0, r8          @ r0 = b[0] << j  (This is t[0])
    
    LDR     r1, [r4]            @ Load x[0]
    EOR     r1, r1, r0          @ x[0] = x[0] ^ t[0]
 
    NOP
    EOR     r11, r11, r11
 
    STR     r1, [r4]            @ Store x[0]


    EOR     r0, r0, r0          @ Clear t[0]
    EOR     r1, r1, r1          @ Clear x[0]

    @ --- [HARDENING] BUS SCRUB ---
    @ Context: We just finished Share 0. Bus holds x[0] or t[0].
    LDR     r12, [r7]           @ Dummy Load Randomness

    @ --- Process Share 1 ---
    LDR     r2, [r9, #8]        @ Load b[1]
    LSL     r2, r2, r8          @ r2 = b[1] << j (This is t[1])

    LDR     r12, [r7] 
    
    LDR     r3, [r4, #4]        @ Load x[1]
    EOR     r3, r3, r2          @ x[1] = x[1] ^ t[1]

    NOP
    EOR     r11, r11, r11    

    LDR     r12, [r7] 

    STR     r3, [r4, #4]        @ Store x[1]

    EOR     r2, r2, r2
    EOR     r3, r3, r3

    @ --- [HARDENING] BUS SCRUB ---
    @ Context: Finished Share 1. Scrub before Refresh function call.
    LDR     r12, [r7]           

    @ -------------------------------------------------------------------------
    @ Step E: Refresh
    @ -------------------------------------------------------------------------
    @ Refresh is still needed to prevent "horizontal" leakage accumulation
    MOV     r0, r4              @ R0 = x
    MOV     r1, r7              @ R1 = pool

    BL      Refresh_m4
    
    ADD     r7, r7, #4          @ Refresh consumes 1 word

    @ -------------------------------------------------------------------------
    @ Loop Control
    @ -------------------------------------------------------------------------
    ADD     r8, r8, #1          @ j++
    B       loop_geo

end_geo:
    ADD     sp, sp, #32
    POP     {r4-r11, pc}

.section .text
.syntax unified
.cpu cortex-m4
.thumb

.global MaskedLaplace_m4
.extern MaskedBernoulli_m4
.extern MaskedGeometric_m4
.extern SecNEG_m4
.extern SecINC_m4
.extern SecAND_m4

@ -----------------------------------------------------------------------------
@ Constants
@ -----------------------------------------------------------------------------
.equ INC_RC_BERN, 16
.equ INC_PC_BERN, 32
.equ INC_PC_GEO, 40
.equ INC_PC_NEG, 20
.equ INC_PC_INC, 20

@ Stack Layout
@ 12 Bytes total:
@ [SP + 0] : Arg4 slot for sub-calls (Pool pointer)
@ [SP + 4] : Storage for t[0] 
@ [SP + 8] : Storage for t[1] 
.equ STACK_ALLOC, 12

MaskedLaplace_m4:
    PUSH    {r4-r11, lr}
    SUB     sp, sp, #STACK_ALLOC

    @ Load Stack Args
    LDR     r8, [sp, #48]       @ r8 = rs
    LDR     r7, [sp, #52]       @ r7 = pool

    MOV     r4, r0              @ r4 = v
    MOV     r5, r1              @ r5 = p_lap
    MOV     r6, r2              @ r6 = p_geo
    MOV     r10, r3             @ r10 = one

    @ 1. MaskedBernoulli
    MOV     r0, r4
    MOV     r1, r5
    MOV     r2, r10
    MOV     r3, r8
    STR     r7, [sp]
    BL      MaskedBernoulli_m4

    ADD     r8, r8, #INC_RC_BERN
    ADD     r7, r7, #INC_PC_BERN

    LDR     r12, [r7]

    @ 2. Save 'b' & Clear v
    LDR     r9, [r4]            @ r9 = b[0]

    LDR     r12, [r7], #4

    LDR     r11, [r4, #4]       @ r11 = b[1]

    @ [HARDENING] Bus Scrub
    LDR     r12, [r7]
    EOR     r12, r12, r12

    MOV     r0, #0
    STR     r0, [r4]            @ v[0] = 0
    STR     r0, [r4, #4]        @ v[1] = 0

    @ 3. MaskedGeometric
    MOV     r0, r4
    MOV     r1, r6
    MOV     r2, r10
    MOV     r3, r8
    STR     r7, [sp]
    BL      MaskedGeometric_m4

    ADD     r7, r7, #INC_PC_GEO

    LDR     r12, [r7]

    @ 4. Save 'x' to Registers
    LDR     r5, [r4]            @ r5 = x[0]

    NOP
    LDR     r12, [r7]
    
    LDR     r6, [r4, #4]        @ r6 = x[1]

    @ [HARDENING] Bus Scrub
    LDR     r12, [r7, #4]

    @ 5. SecNEG (Compute t directly into Stack)
    MVN     r9, r9              @ ~b[0]
    NOP
    EOR     r12, r12, r12

    STR     r9, [sp, #4]        @ Stack[0] = ~b[0]

    NOP
    LDR     r12, [r7]

    STR     r11, [sp, #8]       @ Stack[1] = b[1]
    NOP
    EOR     r12, r12, r12
    LDR     r12, [r7]
    
    ADD     r0, sp, #4
    ADD     r1, sp, #4
    MOV     r2, r7
    BL      SecNEG_m4

    ADD     r7, r7, #INC_PC_NEG

    LDR     r12, [r7]

    @ 6. SecINC (Compute u directly into v)
    STR     r5, [r4]            @ v[0] = x[0]
 
    NOP
    LDR     r12, [r7]
 
    STR     r6, [r4, #4]        @ v[1] = x[1]

    NOP
    LDR     r12, [r7]

    @ --- FIX START ---
    MOV     r0, r4              @ R0 = v (Output)
    MOV     r1, r4              @ R1 = v (Input x)
    MOV     r2, r10             @ R2 = one 
    MOV     r3, r7              @ R3 = pool (Correct placement)
    @ --- FIX END ---
    
    BL      SecINC_m4

    ADD     r7, r7, #INC_PC_INC

    LDR     r12, [r7]
    
    @ -------------------------------------------------------------------------
    @ 7. SecAND
    @ -------------------------------------------------------------------------
    MOV     r0, r4              @ Arg0: Output v
    ADD     r1, sp, #4          @ Arg1: Input t (Stack)
    MOV     r2, r4              @ Arg2: Input u (v)
    MOV     r3, r7              @ Arg3: Pool
    BL      SecAND_m4

    @ Stack Scrubbing (Crucial for final clean)
    MOV     r0, #0
    STR     r0, [sp, #4]
    STR     r0, [sp, #8]

    @ Cleanup
    EOR     r12, r12, r12
    LDR     r12, [r7]
    ADD     sp, sp, #STACK_ALLOC
    POP     {r4-r11, pc}

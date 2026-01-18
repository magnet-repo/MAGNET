.section .text
.syntax unified
.cpu cortex-m4
.thumb

.global SecINC_m4

@ -----------------------------------------------------------------------------
@ External Functions
@ -----------------------------------------------------------------------------
@ SecBA_m4 = Boolean to Arithmetic (B2A)
@ SecAB_m4 = Arithmetic to Boolean (A2B)
.extern SecAB_m4
.extern SecBA_m4

@ -----------------------------------------------------------------------------
@ void SecINC_m4(uint32_t *z, const uint32_t *x, const uint32_t *one, const uint32_t *pool)
@ -----------------------------------------------------------------------------
SecINC_m4:
    @ 1. Prologue
    PUSH    {r4-r6, lr}
    SUB     sp, sp, #8      @ Allocate 8 bytes for intermediate Arithmetic Shares

    @ Save pointers
    MOV     r4, r0          @ z_ptr (Output)
    MOV     r5, r1          @ x_ptr (Input)
    MOV     r6, r3          @ pool_ptr

    @ ---------------------------------------------------------
    @ STEP 1: Boolean to Arithmetic (B2A)
    @ ---------------------------------------------------------
    MOV     r0, sp          @ Arg0: dst (Stack -> Arithmetic Temp)
    MOV     r1, r5          @ Arg1: src (x_ptr -> Boolean Input)
    MOV     r2, r6          @ Arg2: pool
    BL      SecBA_m4        

    @ Advance Pool
    ADD     r6, r6, #4

    LDR     r12, [r6]   @ Dummy Random Load

    NOP
    EOR     r12, r12, r12

    @ ---------------------------------------------------------
    @ STEP 2: Increment (Arithmetic Domain)
    @ A = A + 1
    @ ---------------------------------------------------------
    LDR     r0, [sp, #0]    @ Load Arithmetic Share 0
 
    LDR     r12, [r6]   @ Dummy Random Load
 
    ADD     r0, r0, #1      @ Increment
    STR     r0, [sp, #0]    @ Store back

    LDR     r12, [r6]   @ Dummy Random Load

    @ ---------------------------------------------------------
    @ STEP 3: Arithmetic to Boolean (A2B)
    @ ---------------------------------------------------------
    MOV     r0, r4          @ Arg0: dst (z_ptr -> Boolean Output)
    MOV     r1, sp          @ Arg1: src (Stack -> Arithmetic Temp)
    MOV     r2, r6          @ Arg2: pool
 
    BL      SecAB_m4        

    LDR     r12, [r6]   @ Dummy Random Load

    @ ---------------------------------------------------------
    @ Epilogue
    @ ---------------------------------------------------------
    @ Scrub stack
    MOV     r0, #0
    STR     r0, [sp, #0]
    STR     r0, [sp, #4]

    ADD     sp, sp, #8
    POP     {r4-r6, pc}

.section .text
.syntax unified
.cpu cortex-m4
.thumb

.global SecADD_m4

@ -----------------------------------------------------------------------------
@ External Functions
@ -----------------------------------------------------------------------------
.extern SecAB_m4
.extern SecBA_m4

@ -----------------------------------------------------------------------------
@ void SecADD_m4(uint32_t *z, const uint32_t *x, const uint32_t *y, const uint32_t *pool)
@ -----------------------------------------------------------------------------
SecADD_m4:
    @ 1. Prologue
    @ Push 6 registers (24 bytes) to maintain 8-byte stack alignment
    PUSH    {r4-r8, lr}
    
    @ Allocate 16 bytes for intermediate variables:
    @ [sp, #0]: xa (2 words) / za (2 words)
    @ [sp, #8]: ya (2 words)
    SUB     sp, sp, #16

    @ Move arguments to preserved registers
    MOV     r4, r0          @ z_ptr
    MOV     r5, r1          @ x_ptr
    MOV     r6, r2          @ y_ptr
    MOV     r7, r3          @ pool_ptr

    @ ---------------------------------------------------------
    @ STEP 1: xa = BooleanToArithmetic(x)
    @ ---------------------------------------------------------
    MOV     r0, sp          @ Arg0: dst (storage for xa)
    MOV     r1, r5          @ Arg1: src (x_ptr)
    MOV     r2, r7          @ Arg2: pool
    BL      SecBA_m4

    @ Advance pool pointer (assuming SecBA used 1 random word)
    ADD     r7, r7, #4

    @ [SCRUB 1] Bus likely holds a share of X. Scrub before touching Y.
    LDR     r12, [r7, #-4]  @ Load "used" random word to clear bus

    @ ---------------------------------------------------------
    @ STEP 2: ya = BooleanToArithmetic(y)
    @ ---------------------------------------------------------
    ADD     r0, sp, #8      @ Arg0: dst (storage for ya at sp+8)
    MOV     r1, r6          @ Arg1: src (y_ptr)
    MOV     r2, r7          @ Arg2: pool
    BL      SecBA_m4

    @ Advance pool pointer
    ADD     r7, r7, #4

    @ [SCRUB 2] Bus likely holds a share of Y. Scrub before Math.
    LDR     r12, [r7, #-4]  

    @ ---------------------------------------------------------
    @ STEP 3: Arithmetic Addition (Component-wise)
    @ ---------------------------------------------------------
   
    @ --- za[0] = xa[0] + ya[0] ---
    LDR     r0, [sp, #0]    @ Load xa[0]
    LDR     r1, [sp, #8]    @ Load ya[0]

    LDR     r12, [r7, #-4]  @ Dummy load random value
    
    NOP
    EOR     r12, r12, r12

    ADD     r0, r0, r1      @ Add

    NOP
    EOR     r12, r12, r12

    STR     r0, [sp, #0]    @ Store result back to stack (za[0])

    @ [SCRUB 3] CRITICAL!
    @ Bus now holds za[0]. If we immediately load xa[1], 
    @ we leak HW(za[0] ^ xa[1]).
    LDR     r12, [r7, #-4]  @ Dummy load random value
    
    
    @ --- za[1] = xa[1] + ya[1] ---
    LDR     r0, [sp, #4]    @ Load xa[1]
    LDR     r1, [sp, #12]   @ Load ya[1]

    LDR     r12, [r7, #-4]  @ Dummy load random value
    
    NOP
    EOR     r12, r12, r12

    ADD     r0, r0, r1      @ Add

    NOP
    EOR     r12, r12, r12

    LDR     r12, [r7, #-4]  @ Dummy load random value

    STR     r0, [sp, #4]    @ Store result back to stack (za[1])

    @ [SCRUB 4] Scrub before calling A2B
    LDR     r12, [r7, #-4]

    @ ---------------------------------------------------------
    @ STEP 4: z = ArithmeticToBoolean(za)
    @ ---------------------------------------------------------
    MOV     r0, r4          @ Arg0: dst (z_ptr)
    MOV     r1, sp          @ Arg1: src (points to za at sp)
    MOV     r2, r7          @ Arg2: pool
    BL      SecAB_m4

    @ ---------------------------------------------------------
    @ Epilogue
    @ ---------------------------------------------------------
    @ Optional: Zero out stack to clear secrets
    MOV     r0, #0
    STR     r0, [sp, #0]
    STR     r0, [sp, #4]
    STR     r0, [sp, #8]
    STR     r0, [sp, #12]

    ADD     sp, sp, #16     @ Deallocate stack
    POP     {r4-r8, pc}     @ Restore regs and return

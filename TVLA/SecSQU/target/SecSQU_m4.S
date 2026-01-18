.section .text
.syntax unified
.cpu cortex-m4
.thumb

.global SecSQU_m4

@ -----------------------------------------------------------------------------
@ void SecSQU_m4(uint32_t* C, const uint32_t* A, const uint32_t* pool)
@ -----------------------------------------------------------------------------
SecSQU_m4:
    PUSH    {r4-r10, lr}
    SUB     sp, sp, #8          @ Stack for dummy writes

    @ ---------------------------------------------------------
    @ 1. PREPARE BLINDING FACTOR
    @ ---------------------------------------------------------
    @ We need two random values: 
    @   r7 = r (for output masking)
    @   r9 = rho (for input blinding)
    LDR     r7, [r2]            @ Load r (pool[0])
    LDR     r9, [r2, #4]        @ Load rho (pool[1])
    LDR     r10, [r2, #8]        @ t
    
    @ ---------------------------------------------------------
    @ 2. LOAD & BLIND INPUTS
    @ ---------------------------------------------------------
    @ FIX: Randomize shares BEFORE they meet in the ALU.
    @ Prevents HD(A0, A1) leakage which reveals the secret.

    LDR     r3, [r1]            @ Load A[0]
    SUB     r3, r3, r9          @ A[0]' = A[0] - rho

    @ [HARDENING] Scrub Read Bus (Transition A[0] -> A[1])
    LDR     r12, [r2]           @ Dummy Load
    EOR     r12, r12, r12       @ Clear reg
    
    LDR     r4, [r1, #4]        @ Load A[1]
    ADD     r4, r4, r9          @ A[1]' = A[1] + rho

    @ Now r3 and r4 are freshly randomized. 
    @ HD(r3, r4) is random and does not leak the secret.

    
    EOR     r12, r12, r12       @ Clear history
    
    MUL     r12, r7, r7         @ Scratch = r * r
    
    @ ---------------------------------------------------------
    @ 3. COMPUTE SQUARES
    @ ---------------------------------------------------------
    MUL     r5, r3, r3          @ r5 = (A[0]')^2
    
    
    EOR     r12, r12, r12       @ Clear history
    
    MUL     r12, r7, r7         @ Scratch = r * r
    

    MUL     r6, r4, r4          @ r6 = (A[1]')^2

    @ ---------------------------------------------------------
    @ 4. COMPUTE CROSS TERMS
    @ ---------------------------------------------------------
    @ [HARDENING] Dummy Multiplication to pre-charge ALU
    
    EOR     r12, r12, r12       @ Clear history
    
    MUL     r12, r7, r7         @ Scratch = r * r
    
    
    @ Compute Cross Term: 2 * A[0]' * A[1]'
    
    ADD     r11, r3, r10            @ u  = r3 + t

    
    EOR     r12, r12, r12       @ Clear history

    MUL     r8,  r11, r4            @ m1 = u * r4
   
    
    EOR     r12, r12, r12       @ Clear history
    
    
    MUL     r12, r10, r4            @ m2 = t * r4 

    
    EOR     r9, r9, r9          @ Clear history

    SUB     r8,  r8,  r12           @ r8 = r3 * r4

    @MUL     r8, r3, r4          @ r8 = A[0]' * A[1]'

    
    EOR     r12, r12, r12       @ Clear history
    
    MUL     r12, r7, r7         @ Scratch = r * r
    

    LSL     r8, r8, #1          @ r8 = 2 * A[0]' * A[1]'

    
    EOR     r12, r12, r12       @ Clear history
    
    MUL     r12, r7, r7         @ Scratch = r * r
    

    @ Apply Output Mask
    ADD     r8, r7, r8          @ r8 = r + 2*A[0]'*A[1]'

    
    EOR     r12, r12, r12       @ Clear history
    
    MUL     r12, r7, r7         @ Scratch = r * r
    

    @ ---------------------------------------------------------
    @ 5. UPDATE SHARES
    @ ---------------------------------------------------------
    SUB     r5, r5, r7          @ C[0] = (A[0]')^2 - r
    
    
    EOR     r12, r12, r12       @ Clear history
    
    MUL     r12, r7, r7         @ Scratch = r * r
    

    ADD     r6, r6, r8          @ C[1] = (A[1]')^2 + 2*A[0]'*A[1]' + r

    
    EOR     r9, r9, r9          @ Clear history

    @ ---------------------------------------------------------
    @ 6. STORE RESULTS (HARDENED)
    @ ---------------------------------------------------------
    STR     r5, [r0]            @ Store C[0]

    @ [HARDENING] Scrub Write Bus (Transition C[0] -> C[1])
    LDR     r12, [r2]           @ Load Random
    STR     r12, [sp]           @ Dummy Store to Stack

    STR     r6, [r0, #4]        @ Store C[1]

    @ Cleanup
    ADD     sp, sp, #8
    POP     {r4-r10, pc}

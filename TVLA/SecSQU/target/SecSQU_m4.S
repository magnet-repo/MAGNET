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

    LDR     r7, [r2]            @ Load r (pool[0])
    LDR     r9, [r2, #4]        @ Load rho (pool[1])
    LDR     r10, [r2, #8]       @ t
    
    LDR     r3, [r1]            @ Load A[0]
    SUB     r3, r3, r9          @ A[0]' = A[0] - rho

    LDR     r12, [r2]           @ Dummy Load
    EOR     r12, r12, r12       @ Clear reg
    
    LDR     r4, [r1, #4]        @ Load A[1]
    ADD     r4, r4, r9          @ A[1]' = A[1] + rho

    NOP
    EOR     r12, r12, r12       @ Clear history
    NOP
    MUL     r12, r7, r7         @ Scratch = r * r
    NOP

    MUL     r5, r3, r3          @ r5 = (A[0]')^2
    
    NOP
    EOR     r12, r12, r12       @ Clear history
    NOP
    MUL     r12, r7, r7         @ Scratch = r * r
    NOP

    MUL     r6, r4, r4          @ r6 = (A[1]')^2

    NOP
    EOR     r12, r12, r12       @ Clear history
    NOP
    MUL     r12, r7, r7         @ Scratch = r * r
    NOP
    
    ADD     r11, r3, r10            @ u  = r3 + t

    NOP
    EOR     r12, r12, r12       @ Clear history
    NOP
    MUL     r12, r7, r7         @ Dummy MUL to pre-charge ALU
    NOP

    MUL     r8,  r11, r4        @ m1 = u * r4         

    NOP
    NOP
    EOR     r12, r12, r12       @ Clear history
    NOP
    MUL     r12, r7, r7         @ Dummy MUL to pre-charge ALU
    NOP
    NOP
    
    MUL     r12, r10, r4        @ m2 = t * r4

    NOP
    NOP
    EOR     r9, r9, r9          @ Clear history
    NOP
    MUL     r9, r7, r7          @ Dummy MUL to pre-charge ALU
    NOP
    NOP

    SUB     r8,  r8,  r12       @ r8 = r3 * r4

    NOP
    EOR     r12, r12, r12       @ Clear history
    NOP
    MUL     r12, r7, r7         @ Scratch = r * r
    NOP

    LSL     r8, r8, #1          @ r8 = 2 * A[0]' * A[1]'

    NOP
    EOR     r12, r12, r12       @ Clear history
    NOP
    MUL     r12, r7, r7         @ Scratch = r * r
    NOP

    ADD     r8, r7, r8          @ r8 = r + 2*A[0]'*A[1]'

    NOP
    EOR     r12, r12, r12       @ Clear history
    NOP
    MUL     r12, r7, r7         @ Scratch = r * r
    NOP

    SUB     r5, r5, r7          @ C[0] = (A[0]')^2 - r
    
    NOP
    EOR     r12, r12, r12       @ Clear history
    NOP
    MUL     r12, r7, r7         @ Scratch = r * r
    NOP

    ADD     r6, r6, r8          @ C[1] = (A[1]')^2 + 2*A[0]'*A[1]' + r

    NOP
    EOR     r9, r9, r9          @ Clear history

    LDR     r12, [r2]           @ Load Random

    STR     r5, [r0]            @ Store C[0]
    NOP

    LDR     r12, [r2]           @ Load Random
    
    NOP
    STR     r6, [r0, #4]        @ Store C[1]

    LDR     r12, [r2]           @ Load Random

    POP     {r4-r10, pc}

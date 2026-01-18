.section .text
.syntax unified
.cpu cortex-m4
.thumb

.global Refresh_m4

@ -----------------------------------------------------------------------------
@ void Refresh_m4(uint32_t* x, uint32_t* pool)
@ Arguments:
@   r0 = x    (Pointer to 2 shares)
@   r1 = pool (Pointer to Randomness Stream)
@ -----------------------------------------------------------------------------
Refresh_m4:
    PUSH    {r4-r6, lr}
    SUB     sp, sp, #8          @ Allocate Scrub Slot on Stack

    @ 1. Load the Mask 'r'
    LDR     r2, [r1], #4        @ r = pool[next]

    @ =========================================================
    @ UPDATE SHARE 0
    @ =========================================================
    LDR     r3, [r0]            @ Load x[0]
    EOR     r3, r3, r2          @ x[0] ^= r
    STR     r3, [r0]            @ Store x[0]
    
    @ Clear register to prevent ALU leakage
    MOV     r3, #0
    
    LDR     r5, [r1]            @ Random Load
    EOR     r5, r5, r5          @ Dummy xor
    
    @ =========================================================
    @ UPDATE SHARE 1
    @ =========================================================
    LDR     r4, [r0, #4]        @ Load x[1]
    EOR     r4, r4, r2          @ x[1] ^= r
    STR     r4, [r0, #4]        @ Store x[1]
    
    @ Clear register
    MOV     r4, #0

    LDR     r5, [r1]            @ Random Load

    @ Cleanup
    ADD     sp, sp, #8
    POP     {r4-r6, pc}

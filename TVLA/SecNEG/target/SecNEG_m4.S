.section .text
.syntax unified
.cpu cortex-m4
.thumb

.global SecNEG_m4
.extern Refresh_m4

@ -----------------------------------------------------------------------------
@ void SecNEG(uint32_t* z, const uint32_t* x, const uint32_t* pool)
@ Arguments:
@   R0 = z    (Output pointer)
@   R1 = x    (Input pointer)
@   R2 = pool (Randomness pointer)
@ -----------------------------------------------------------------------------
SecNEG_m4:
    PUSH    {r4-r6, lr}
    
    @ Allocate Scrub Slot
    SUB     sp, sp, #8

    @ Preserve pointers
    MOV     r4, r0              @ r4 = z
    MOV     r5, r2              @ r5 = pool

    @ =========================================================
    @ PROCESS SHARE 0: z[0] = -(x[0] & 1)
    @ =========================================================
    LDR     r3, [r1]            @ Load x[0]
    AND     r3, r3, #1          @ Mask LSB
    RSB     r3, r3, #0          @ Negate: 0 - r3 (0 -> 0, 1 -> -1)
    STR     r3, [r4]            @ Store z[0]

    @ =========================================================
    @ SCRUBBING (Break correlation between z[0] and z[1])
    @ =========================================================
    LDR     r6, [r5]            @ Load Random (Peek)
    
    @ Scrub Write Bus
    STR     r6, [sp]            @ Dummy Store to Stack
    @DSB

    LDR     r6, [r5]

    @ =========================================================
    @ PROCESS SHARE 1: z[1] = -(x[1] & 1)
    @ =========================================================
    LDR     r3, [r1, #4]        @ Load x[1]
    AND     r3, r3, #1          @ Mask LSB
    RSB     r3, r3, #0          @ Negate
    STR     r3, [r4, #4]        @ Store z[1]

    MOV     r3, #0
    @ =========================================================
    @ CALL REFRESH
    @ =========================================================
    @ Refresh_m4(z, pool)
    MOV     r0, r4              @ R0 = z
    MOV     r1, r5              @ R1 = pool
    
    BL      Refresh_m4

    @ Cleanup
    ADD     sp, sp, #8
    POP     {r4-r6, pc}
    

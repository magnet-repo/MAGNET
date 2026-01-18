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

.equ STACK_ALLOC, 16

MaskedLaplace_m4:
    PUSH    {r4-r11, lr}
    SUB     sp, sp, #STACK_ALLOC

    MOV     r4, r0
    MOV     r5, r1
    MOV     r6, r2
    MOV     r10, r3

    MOV     r0, r4
    MOV     r1, r5
    MOV     r2, r10
    BL      MaskedBernoulli_m4

    LDR     r0, [r4]
    STR     r0, [sp, #8]
    LDR     r1, [r4, #4]
    STR     r1, [sp, #12]

    MOV     r0, #0
    STR     r0, [r4]
    STR     r0, [r4, #4]

    MOV     r0, r4
    MOV     r1, r6
    MOV     r2, r10
    BL      MaskedGeometric_m4

    LDR     r0, [sp, #8]
    MVN     r0, r0
    STR     r0, [sp]
    LDR     r1, [sp, #12]
    STR     r1, [sp, #4]

    ADD     r0, sp, #8
    MOV     r1, sp
    BL      SecNEG_m4

    MOV     r0, r4
    MOV     r1, r4
    MOV     r2, r10
    BL      SecINC_m4

    MOV     r0, r4
    ADD     r1, sp, #8
    MOV     r2, r4
    BL      SecAND_m4

    ADD     sp, sp, #STACK_ALLOC
    POP     {r4-r11, pc}

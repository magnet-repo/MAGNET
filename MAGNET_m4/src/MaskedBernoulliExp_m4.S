.section .text
.syntax unified
.cpu cortex-m4
.thumb

.global MaskedBernoulliExp_m4
.extern MaskedBernoulli_m4
.extern Refresh_m4
.extern SecOR_m4
.extern SecAND_m4
.extern l
.extern mu

MaskedBernoulliExp_m4:
    PUSH    {r4-r11, lr}
    MOV     r4, r0
    MOV     r5, r1
    MOV     r6, r2
    MOV     r10, r3

    SUB     sp, sp, #32

    LDR     r0, [r10]
    STR     r0, [r4]
    LDR     r1, [r10, #4]
    STR     r1, [r4, #4]

    LDR     r0, =l
    LDR     r9, [r0]
    MOV     r11, #0

loop_exp:
    CMP     r11, r9
    BGE     end_exp

    ADD     r0, sp, #4
    MOV     r1, r6
    MOV     r2, r10
    BL      MaskedBernoulli_m4

    LDR     r0, =mu
    LDR     r1, [r0]
    LSL     r1, r1, #3
    ADD     r6, r6, r1

    LDR     r0, [r5]
    LSR     r0, r0, r11
    STR     r0, [sp, #12]
    LDR     r1, [r5, #4]
    LSR     r1, r1, r11
    STR     r1, [sp, #16]

    ADD     r0, sp, #12
    BL      Refresh_m4

    LDR     r0, [sp, #12]
    MVN     r0, r0
    STR     r0, [sp, #12]

    ADD     r0, sp, #20
    ADD     r1, sp, #12
    ADD     r2, sp, #4
    BL      SecOR_m4

    MOV     r0, r4
    MOV     r1, r4
    ADD     r2, sp, #20
    BL      SecAND_m4

    ADD     r11, r11, #1
    B       loop_exp

end_exp:
    ADD     sp, sp, #32
    POP     {r4-r11, pc}

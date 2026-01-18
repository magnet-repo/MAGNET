.section .text
.syntax unified
.cpu cortex-m4
.thumb

.global MaskedBernoulli_m4
.extern SecAND_m4
.extern Refresh_m4
.extern xoshiro_next
.extern s
.extern mu

MaskedBernoulli_m4:
    PUSH    {r4-r11, lr}
    MOV     r4, r0
    MOV     r5, r1
    MOV     r6, r2
    
    SUB     sp, sp, #32

    LDR     r10, [r6]
    LDR     r11, [r6, #4]
    STR     r10, [r4]
    STR     r11, [r4, #4]

    LDR     r0, =mu
    LDR     r8, [r0]
    SUB     r8, r8, #1

loop_bern:
    CMP     r8, #0
    BLT     end_bern

    PUSH    {r0-r3}
    LDR     r0, =s
    BL      xoshiro_next
    MOV     r10, r0
    POP     {r0-r3}
    AND     r10, r10, #1
    STR     r10, [sp]

    PUSH    {r0-r3}
    LDR     r0, =s
    BL      xoshiro_next
    MOV     r11, r0
    POP     {r0-r3}
    AND     r11, r11, #1
    STR     r11, [sp, #4]

    MOV     r0, sp
    BL      Refresh_m4

    LDR     r10, [sp]
    LDR     r11, [sp, #4]

    LSL     r12, r8, #3
    ADD     r12, r5, r12

    LDR     r0, [r12]
    EOR     r0, r0, r10
    STR     r0, [sp, #24]

    LDR     r1, [r12, #4]
    EOR     r1, r1, r11
    STR     r1, [sp, #28]

    MVN     r10, r10

    LDR     r0, [r4]
    EOR     r0, r0, r10
    STR     r0, [sp, #8]

    LDR     r1, [r4, #4]
    EOR     r1, r1, r11
    STR     r1, [sp, #12]

    ADD     r0, sp, #16
    ADD     r1, sp, #24
    ADD     r2, sp, #8
    BL      SecAND_m4

    LDR     r0, [sp, #16]
    LDR     r10, [r4]
    EOR     r10, r10, r0
    STR     r10, [r4]

    LDR     r1, [sp, #20]
    LDR     r11, [r4, #4]
    EOR     r11, r11, r1
    STR     r11, [r4, #4]

    SUB     r8, r8, #1
    B       loop_bern

end_bern:
    ADD     sp, sp, #32
    POP     {r4-r11, pc}

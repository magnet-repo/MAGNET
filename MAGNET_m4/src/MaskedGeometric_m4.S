.section .text
.syntax unified
.cpu cortex-m4
.thumb

.global MaskedGeometric_m4
.extern MaskedBernoulli_m4
.extern Refresh_m4
.extern SecOR_m4
.extern kap
.extern mu

MaskedGeometric_m4:
    PUSH    {r4-r11, lr}
    MOV     r4, r0
    MOV     r5, r1
    MOV     r6, r2
    
    SUB     sp, sp, #32

    LDR     r0, =kap
    LDR     r10, [r0]
    MOV     r8, #0

loop_geo:
    CMP     r8, r10
    BGE     end_geo

    MOV     r0, sp
    MOV     r1, r5
    MOV     r2, r6
    BL      MaskedBernoulli_m4

    LDR     r0, =mu
    LDR     r1, [r0]
    LSL     r1, r1, #3
    ADD     r5, r5, r1

    LDR     r0, [sp]
    LSL     r0, r0, r8
    STR     r0, [sp, #8]
    LDR     r1, [sp, #4]
    LSL     r1, r1, r8
    STR     r1, [sp, #12]

    ADD     r0, sp, #16
    MOV     r1, r4
    ADD     r2, sp, #8
    BL      SecOR_m4

    LDR     r0, [sp, #16]
    STR     r0, [r4]
    LDR     r1, [sp, #20]
    STR     r1, [r4, #4]

    MOV     r0, r4
    BL      Refresh_m4

    ADD     r8, r8, #1
    B       loop_geo

end_geo:
    ADD     sp, sp, #32
    POP     {r4-r11, pc}

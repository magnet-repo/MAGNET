.syntax unified
.cpu cortex-m4
.thumb

.global xoshiro_next
.type xoshiro_next, %function

xoshiro_next:
    PUSH    {r4-r7, lr}

    LDR     r1, [r0, #0]
    LDR     r2, [r0, #4]
    LDR     r3, [r0, #8]
    LDR     r4, [r0, #12]

    ADD     r5, r2, r2, LSL #2

    ROR     r5, r5, #25

    ADD     r5, r5, r5, LSL #3

    LSL     r6, r2, #9

    EOR     r3, r3, r1
    EOR     r4, r4, r2
    EOR     r2, r2, r3
    EOR     r1, r1, r4
    EOR     r3, r3, r6
    ROR     r4, r4, #21

    STR     r1, [r0, #0]
    STR     r2, [r0, #4]
    STR     r3, [r0, #8]
    STR     r4, [r0, #12]

    MOV     r0, r5

    POP     {r4-r7, pc}

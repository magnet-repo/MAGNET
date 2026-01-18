.section .text
.syntax unified
.cpu cortex-m4
.thumb

.global SecOR_m4
.extern SecAND_m4

SecOR_m4:
    PUSH    {r4-r11, lr}
    MOV     r4, r0
    MOV     r5, r1
    MOV     r6, r2
    
    SUB     sp, sp, #24
    MOV     r8, sp

    LDR     r0, [r5]
    MVN     r0, r0
    STR     r0, [r8]
    LDR     r1, [r5, #4]
    STR     r1, [r8, #4]

    LDR     r0, [r6]
    MVN     r0, r0
    STR     r0, [r8, #8]
    LDR     r1, [r6, #4]
    STR     r1, [r8, #12]

    ADD     r0, r8, #16
    MOV     r1, r8
    ADD     r2, r8, #8
    BL      SecAND_m4

    LDR     r1, [r8, #20]
    STR     r1, [r4, #4]
    LDR     r0, [r8, #16]
    MVN     r0, r0
    STR     r0, [r4]

    ADD     sp, sp, #24
    POP     {r4-r11, pc}

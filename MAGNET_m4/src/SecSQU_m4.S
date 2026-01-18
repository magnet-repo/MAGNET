.section .text
.syntax unified
.cpu cortex-m4
.thumb

.global SecSQU_m4
.extern xoshiro_next
.extern s

SecSQU_m4:
    PUSH    {r4-r8, lr}
    
    LDR     r3, [r1]
    LDR     r4, [r1, #4]

    MUL     r5, r3, r3
    MUL     r6, r4, r4

    PUSH {r0-r3}
    LDR r0, =s
    BL xoshiro_next
    MOV r7, r0
    POP {r0-r3}

    MUL     r8, r3, r4
    LSL     r8, r8, #1
    ADD     r8, r7, r8
    SUB     r5, r5, r7
    ADD     r6, r6, r8

    STR     r5, [r0]
    STR     r6, [r0, #4]
    
    POP     {r4-r8, pc}

.section .text
.syntax unified
.cpu cortex-m4
.thumb

.global SecAB_m4
.extern xoshiro_next
.extern s

SecAB_m4:
    PUSH {r4-r11, r14}

    LDR r4, [r1], #4
    LDR r5, [r1], #4

    PUSH {r0-r3}
    LDR r0, =s
    BL xoshiro_next
    MOV r8, r0
    POP {r0-r3}

    LSL r3, r8, #1
    AND r9, r8, r5
    EOR r9, r9, r8

    EOR r6, r3, r4
    AND r11, r3, r4
    EOR r8, r8, r6

    AND r8, r8, r5
    EOR r9, r9, r8
    EOR r9, r9, r11

    .rept 31
        AND r8, r3, r5
        EOR r8, r8, r9
        AND r3, r3, r4
        EOR r8, r8, r3
        LSL r3, r8, #1
    .endr

    EOR r6, r6, r3

    STR r6, [r0], #4
    STR r5, [r0], #4

    POP {r4-r11, pc}

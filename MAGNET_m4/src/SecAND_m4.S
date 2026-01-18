.section .text
.syntax unified
.cpu cortex-m4
.thumb

.global SecAND_m4
.extern xoshiro_next
.extern s

SecAND_m4:
    PUSH {r4-r11, lr}

    PUSH {r0-r3}
    LDR r0, =s
    BL xoshiro_next
    MOV r6, r0
    POP {r0-r3}

    LDR r5, [r2], #4
    LDR r12, [r2], #4

    LDR r10, [r1], #4
    LDR r11, [r1], #4

    AND r7, r5, r10
    EOR r4, r4, r4
    AND r8, r12, r11

    AND r4, r5, r11
    AND r9, r12, r10

    EOR r4, r6, r4
    EOR r4, r9, r4

    EOR r7, r7, r6
    EOR r8, r8, r4

    STR r7, [r0], #4
    STR r8, [r0], #4

    POP {r4-r11, pc}

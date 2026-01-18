.section .text
.syntax unified
.cpu cortex-m4
.thumb

.global SecNEG_m4
.extern Refresh_m4

SecNEG_m4:
    PUSH {r4-r5, lr}

    MOV r4, r0

    LDR r3, [r1]
    AND r3, r3, #1
    RSB r3, r3, #0
    STR r3, [r4]

    LDR r3, [r1, #4]
    AND r3, r3, #1
    RSB r3, r3, #0
    STR r3, [r4, #4]

    MOV r0, r4
    BL Refresh_m4

    POP {r4-r5, pc}

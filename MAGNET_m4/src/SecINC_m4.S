.section .text
.syntax unified
.cpu cortex-m4
.thumb

.global SecINC_m4
.extern SecBA_m4
.extern SecAB_m4

SecINC_m4:
    PUSH {r4-r6, lr}
    SUB sp, sp, #8

    MOV r4, r0
    MOV r5, r1
    MOV r6, r2

    MOV r0, sp
    MOV r1, r5
    BL SecBA_m4

    LDR r0, [sp]
    ADD r0, r0, #1
    STR r0, [sp]

    MOV r0, r4
    MOV r1, sp
    BL SecAB_m4

    ADD sp, sp, #8
    POP {r4-r6, pc}

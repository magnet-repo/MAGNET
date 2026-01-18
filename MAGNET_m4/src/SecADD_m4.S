.section .text
.syntax unified
.cpu cortex-m4
.thumb

.global SecADD_m4
.extern SecBA_m4
.extern SecAB_m4

SecADD_m4:
    PUSH {r4-r8, lr}
    SUB sp, sp, #16

    MOV r4, r0
    MOV r5, r1
    MOV r6, r2

    MOV r0, sp
    MOV r1, r5
    BL SecBA_m4

    ADD r0, sp, #8
    MOV r1, r6
    BL SecBA_m4

    LDR r0, [sp]
    LDR r1, [sp, #8]
    ADD r0, r0, r1
    STR r0, [sp]

    LDR r0, [sp, #4]
    LDR r1, [sp, #12]
    ADD r0, r0, r1
    STR r0, [sp, #4]

    MOV r0, r4
    MOV r1, sp
    BL SecAB_m4

    ADD sp, sp, #16
    POP {r4-r8, pc}

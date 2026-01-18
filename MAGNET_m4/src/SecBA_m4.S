.section .text
.syntax unified
.cpu cortex-m4
.thumb

.global SecBA_m4
.extern xoshiro_next
.extern s

SecBA_m4:
    PUSH {r4-r6, lr}

    PUSH {r0-r3}
    LDR r0, =s
    BL xoshiro_next
    MOV r12, r0
    
    POP {r0}
    POP {r1-r3}
    
    MOV r3, r12

    LDR r4, [r1, #4]

    EOR r5, r3, r4
    STR r4, [r0, #4]

    LDR r6, [r1]

    EOR r12, r6, r3
    SUB r12, r12, r3
    EOR r12, r12, r6

    EOR r6, r6, r5
    SUB r6, r6, r5

    EOR r6, r12, r6
    STR r6, [r0]

    POP {r4-r6, pc}

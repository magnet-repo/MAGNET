.section .text
.syntax unified
.cpu cortex-m4
.thumb

.global Refresh_m4
.extern xoshiro_next
.extern s

Refresh_m4:
    PUSH {r4, lr}

    PUSH {r0-r3}
    LDR r0, =s
    BL xoshiro_next
    MOV r2, r0
    POP {r0-r3}

    LDR r3, [r0]
    EOR r3, r3, r2
    STR r3, [r0]

    LDR r4, [r0, #4]
    EOR r4, r4, r2
    STR r4, [r0, #4]

    POP {r4, pc}

.section .text
.syntax unified
.cpu cortex-m4
.thumb

.global MAGNET_m4

.extern m
.extern s
.extern xoshiro_next

.extern MaskedLaplace_m4
.extern MaskedBernoulliExp_m4
.extern SecADD_m4
.extern SecBA_m4
.extern SecAB_m4
.extern SecSQU_m4
.extern SecNEG_m4
.extern SecINC_m4
.extern SecAND_m4
.extern Refresh_m4

.equ SJ,        0
.equ V_SML,     8
.equ V_BIG,     16
.equ U_BIG,     24
.equ U_SML,     32
.equ B_BIT,     40
.equ SIGN,      48
.equ SIGN_MSK,  56
.equ SJ_TIL,    64
.equ DIFF_MSK,  72
.equ SEL_MSK,   80
.equ T_VAL,     88
.equ B_PRM,     96
.equ TMP_VAL,   104
.equ STACK_SZ,  128

MAGNET_m4:
    PUSH    {r4-r11, lr}
    
    MOV     r4, r0
    MOV     r5, r1
    MOV     r6, r2
    MOV     r7, r3
    LDR     r8, [sp, #36]
    LDR     r9, [sp, #40]

    SUB     sp, sp, #STACK_SZ

    LDR     r0, =m
    LDR     r11, [r0]
    MOV     r10, #0
    
    MOV     r12, #0
    STR     r12, [sp, #120]

loop_magnet:
    LDR     r12, [sp, #120]
    CMP     r12, r11
    BGE     end_magnet

    ADD     r0, sp, #SJ
    MOV     r1, r9
    MOV     r2, r7
    MOV     r3, r6
    BL      MaskedLaplace_m4

    ADD     r0, sp, #V_SML
    ADD     r1, sp, #SJ
    MOV     r2, r5
    BL      SecADD_m4

    ADD     r0, sp, #V_BIG
    ADD     r1, sp, #V_SML
    BL      SecBA_m4

    ADD     r0, sp, #U_BIG
    ADD     r1, sp, #V_BIG
    BL      SecSQU_m4

    ADD     r0, sp, #U_SML
    ADD     r1, sp, #U_BIG
    BL      SecAB_m4

    ADD     r0, sp, #B_BIT
    ADD     r1, sp, #U_SML
    MOV     r2, r8
    MOV     r3, r6
    BL      MaskedBernoulliExp_m4

    LDR     r0, =s
    BL      xoshiro_next
    AND     r0, r0, #1
    STR     r0, [sp, #SIGN]
    
    LDR     r0, =s
    BL      xoshiro_next
    AND     r0, r0, #1
    STR     r0, [sp, #SIGN+4]

    ADD     r0, sp, #SIGN
    BL      Refresh_m4

    ADD     r0, sp, #SIGN_MSK
    ADD     r1, sp, #SIGN
    BL      SecNEG_m4

    LDR     r0, [sp, #SJ]
    LDR     r1, [sp, #SIGN_MSK]
    EOR     r0, r0, r1
    STR     r0, [sp, #SJ]
    LDR     r0, [sp, #SJ+4]
    LDR     r1, [sp, #SIGN_MSK+4]
    EOR     r0, r0, r1
    STR     r0, [sp, #SJ+4]

    ADD     r0, sp, #SJ_TIL
    ADD     r1, sp, #SJ
    MOV     r2, r6
    BL      SecINC_m4

    LDR     r0, [sp, #SJ_TIL]
    LDR     r1, [sp, #SJ]
    EOR     r0, r0, r1
    STR     r0, [sp, #DIFF_MSK]
    LDR     r0, [sp, #SJ_TIL+4]
    LDR     r1, [sp, #SJ+4]
    EOR     r0, r0, r1
    STR     r0, [sp, #DIFF_MSK+4]

    ADD     r0, sp, #SEL_MSK
    ADD     r1, sp, #DIFF_MSK
    ADD     r2, sp, #SIGN_MSK
    BL      SecAND_m4

    LDR     r0, [sp, #SJ]
    LDR     r1, [sp, #SEL_MSK]
    EOR     r0, r0, r1
    STR     r0, [sp, #T_VAL]
    LDR     r0, [sp, #SJ+4]
    LDR     r1, [sp, #SEL_MSK+4]
    EOR     r0, r0, r1
    STR     r0, [sp, #T_VAL+4]

    ADD     r0, sp, #B_PRM
    ADD     r1, sp, #B_BIT
    BL      SecNEG_m4

    ADD     r0, sp, #TMP_VAL
    ADD     r1, sp, #B_PRM
    ADD     r2, sp, #T_VAL
    BL      SecAND_m4

    LSL     r3, r10, #3
    ADD     r3, r4, r3

    LDR     r0, [r3]
    LDR     r1, [sp, #TMP_VAL]
    EOR     r0, r0, r1
    STR     r0, [r3]

    LDR     r0, [r3, #4]
    LDR     r1, [sp, #TMP_VAL+4]
    EOR     r0, r0, r1
    STR     r0, [r3, #4]

    LDR     r0, [sp, #B_BIT]
    LDR     r1, [sp, #B_BIT+4]
    EOR     r0, r0, r1
    ADD     r10, r10, r0

    LDR     r12, [sp, #120]
    ADD     r12, r12, #1
    STR     r12, [sp, #120]
    B       loop_magnet

end_magnet:
    ADD     sp, sp, #STACK_SZ
    POP     {r4-r11, pc}

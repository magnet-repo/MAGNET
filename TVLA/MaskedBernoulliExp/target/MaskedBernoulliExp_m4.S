.section .text
.syntax unified
.cpu cortex-m4
.thumb

.global MaskedBernoulliExp_m4
.extern MaskedBernoulli_m4
.extern Refresh_m4
.extern SecOR_m4
.extern SecAND_m4

@ -----------------------------------------------------------------------------
@ Constants
@ -----------------------------------------------------------------------------
.equ L, 2
.equ P_EXP_STEP, 16
.equ RS_STEP, 16
.equ POOL_BERN_STEP, 8

MaskedBernoulliExp_m4:
    PUSH    {r4-r11, lr}

    @ --- Load Stack Arguments ---
    LDR     r8, [sp, #36]       @ rs
    LDR     r7, [sp, #40]       @ pool

    MOV     r4, r0              @ b (accumulator)
    MOV     r5, r1              @ u (input)
    MOV     r6, r2              @ p_exp
    MOV     r10, r3             @ one

    SUB     sp, sp, #32         @ Stack Buffer

    @ --- Initialize b = one ---
    LDR     r0, [r10]
    STR     r0, [r4]
    LDR     r1, [r10, #4]
    STR     r1, [r4, #4]

    MOV     r11, #0             @ Loop counter j = 0

loop_exp:
    CMP     r11, #L
    BGE     end_exp

    @ -------------------------------------------------------------------------
    @ Step 1: MaskedBernoulli
    @ -------------------------------------------------------------------------
    ADD     r0, sp, #4          @ Output: b_prm
    MOV     r1, r6              @ p_exp
    MOV     r2, r10             @ one
    MOV     r3, r8              @ rs
    STR     r7, [sp]            @ Pass pool
    BL      MaskedBernoulli_m4

    @ Update pointers
    ADD     r6, r6, #P_EXP_STEP
    ADD     r8, r8, #RS_STEP
    ADD     r7, r7, #POOL_BERN_STEP

    @ -------------------------------------------------------------------------
    @ Step 2: Shift u (uj = u >> j)
    @ -------------------------------------------------------------------------
    @ [HARDENING] Dummy load to break correlation with previous step
    LDR     r12, [r7]

    @ --- Share 0 ---
    LDR     r0, [r5]            @ Load u[0]
    LSR     r0, r0, r11         @ Shift
 
    NOP
    EOR     r12, r12, r12
 
    STR     r0, [sp, #12]       @ Store uj[0]

    @ [HARDENING] Scrub Register & Bus before loading Share 1
    @ Prevents HD(u[0], u[1]) leakage
    MOV     r0, #0
    LDR     r12, [r7, #4]       @ Random Dummy Load

    NOP
    EOR     r12, r12, r12

    @ --- Share 1 ---
    LDR     r1, [r5, #4]        @ Load u[1]
    LSR     r1, r1, r11         @ Shift

    NOP
    EOR     r12, r12, r12

    STR     r1, [sp, #16]       @ Store uj[1]

    @ -------------------------------------------------------------------------
    @ Step 3: Refresh(uj)
    @ -------------------------------------------------------------------------
    ADD     r0, sp, #12         @ uj ptr
    MOV     r1, r7              @ pool
    BL      Refresh_m4

    ADD     r7, r7, #4          @ Refresh consumes 1 word

    LDR     r12, [r7]

    @ -------------------------------------------------------------------------
    @ Step 4: Invert uj (uj_til = ~uj)
    @ -------------------------------------------------------------------------
    @ Only invert Share 0 for Boolean masking
    LDR     r0, [sp, #12]       @ Load uj[0]
    MVN     r0, r0              @ Invert

    NOP
    EOR     r12, r12, r12

    STR     r0, [sp, #12]       @ Store back

    @ [HARDENING] Clear register containing sensitive ~uj[0]
    MOV     r0, #0
    LDR     r12, [r7]           @ Random Dummy Load

    @ -------------------------------------------------------------------------
    @ Step 5: SecOR(v, uj_til, b_prm)
    @ -------------------------------------------------------------------------
    ADD     r0, sp, #20         @ Output: v
    ADD     r1, sp, #12         @ Input: uj_til
    ADD     r2, sp, #4          @ Input: b_prm
    MOV     r3, r7              @ Pool
    BL      SecOR_m4

    ADD     r7, r7, #4          @ SecOR consumes 1 word (check implementation)

    LDR     r12, [r7] 

    @ -------------------------------------------------------------------------
    @ Step 6: SecAND(b, b, v)
    @ -------------------------------------------------------------------------
    MOV     r0, r4              @ Output: b (Accumulate)
    MOV     r1, r4              @ Input: b
    ADD     r2, sp, #20         @ Input: v
    MOV     r3, r7              @ Pool
    BL      SecAND_m4

    ADD     r7, r7, #4          @ SecAND consumes 1 word

    LDR     r12, [r7]

    NOP
    EOR     r12, r12, r12

    ADD     r11, r11, #1
    B       loop_exp

end_exp:
    ADD     sp, sp, #32
    POP     {r4-r11, pc}

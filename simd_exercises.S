.section .data
.align 2

.equ n, 16

a_clamp:  .word 100, 300, 255, 0, 512, 128, 256, 77, 999, 1, 254, 400, 50, 600, 200, 255
c_clamp:  .space 64

a_neg: .word 5, -3, 0, -1, 7, -9, 4, -2, 6, -8, 3, -4, 1, -6, 0, -7
result_neg: .word 0

a_dot: .word 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16
b_dot: .word 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16
result_dot: .word 0


.section .text
.globl _start

_start:
    call ex_clamp
    call ex_neg_count
    call ex_dot_product

    li   a7, 93
    li   a0, 0
    ecall


# exercise 1 - clamp values to max 255
ex_clamp:
    addi sp, sp, -16
    sw   ra, 0(sp)

    la   a0, a_clamp
    la   a1, c_clamp
    li   a2, n
    li   t0, 255

clamp_loop:
    vsetvli t1, a2, e32, m1, ta, ma

    vle32.v  v2, (a0)           # load values into v2

    vmsgt.vx v0, v2, t0         # mask: 1 where value > 255

    vmerge.vxm v2, v2, t0, v0  # replace with 255 where mask=1

    vse32.v  v2, (a1)

    slli  t2, t1, 2
    add   a0, a0, t2
    add   a1, a1, t2
    sub   a2, a2, t1
    bnez  a2, clamp_loop

    lw   ra, 0(sp)
    addi sp, sp, 16
    ret


# exercise 2 - count negatives in array
ex_neg_count:
    addi sp, sp, -16
    sw   ra, 0(sp)

    la   a0, a_neg
    li   a2, n
    li   s0, 0          # count

neg_loop:
    vsetvli t1, a2, e32, m1, ta, ma

    vle32.v  v0, (a0)

    vmslt.vx v1, v0, zero   # mask where element < 0

    vcpop.m  t2, v1         # count the 1s in mask

    add   s0, s0, t2

    slli  t3, t1, 2
    add   a0, a0, t3
    sub   a2, a2, t1
    bnez  a2, neg_loop

    la   t0, result_neg
    sw   s0, 0(t0)

    lw   ra, 0(sp)
    addi sp, sp, 16
    ret


# exercise 3 - dot product of two arrays
ex_dot_product:
    addi sp, sp, -16
    sw   ra, 0(sp)

    la   a0, a_dot
    la   a1, b_dot
    li   a2, n

    # initialize accumulator to 0
    vsetvli t0, zero, e32, m1, ta, ma
    vmv.v.x  v3, zero

dot_loop:
    vsetvli t1, a2, e32, m1, ta, ma

    vle32.v  v0, (a0)
    vle32.v  v1, (a1)

    vmul.vv  v2, v0, v1

    vredsum.vs v3, v2, v3   # add into accumulator

    slli  t2, t1, 2
    add   a0, a0, t2
    add   a1, a1, t2
    sub   a2, a2, t1
    bnez  a2, dot_loop

    vmv.x.s  s0, v3         # get scalar result

    la   t0, result_dot
    sw   s0, 0(t0)

    lw   ra, 0(sp)
    addi sp, sp, 16
    ret

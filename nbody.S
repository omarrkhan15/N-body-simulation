.section .data
.align 2

.equ n,          100
.equ dt,         8640
.equ G_int,      0
.equ soft,       1000000000
.equ float_size, 4
.equ arr_bytes,  n * float_size

# hex representations of float constants
# got these by converting manually
.equ G_f,        0x3A8EEFA5
.equ dt_f,       0x46070000
.equ soft_f,     0x4E6E6B28
.equ half,       0x3F000000

filename: .string "shared.mem"

.include "mydata.S"

.section .bss
.align 2
ax: .space arr_bytes
ay: .space arr_bytes
az: .space arr_bytes

.section .text
.globl _start

_start:
    # open the file
    li   a7, 56
    li   a0, -100
    la   a1, filename
    li   a2, 66
    li   a3, 0644
    ecall
    mv   s0, a0
    li   t0, -1
    beq  s0, t0, exit_err

    # mmap it so python can access it
    li   a7, 192
    li   a0, 0
    li   a1, 16384
    li   a2, 3
    li   a3, 1
    mv   a4, s0
    li   a5, 0
    ecall
    mv   s1, a0
    li   t0, -1
    beq  s1, t0, exit_err

main_loop:
    addi sp, sp, -16
    sw   s1, 0(sp)

    # kick drift kick (leapfrog)
    call calc_kick
    call calc_drift
    call calc_acc
    call calc_kick

    lw   s1, 0(sp)
    addi sp, sp, 16

    li   t0, 1
    sw   t0, 0(s1)      # tell python we are writing now

    la   t1, px
    li   t2, 0
    li   t3, arr_bytes
write_px:
    bge  t2, t3, done_px
    add  t4, t1, t2
    lw   t5, 0(t4)
    add  t6, s1, t2
    addi t6, t6, 4
    sw   t5, 0(t6)
    addi t2, t2, 4
    j    write_px
done_px:

    la   t1, py
    li   t2, 0
write_py:
    bge  t2, t3, done_py
    add  t4, t1, t2
    lw   t5, 0(t4)
    add  t6, s1, t2
    addi t6, t6, arr_bytes
    addi t6, t6, 4
    sw   t5, 0(t6)
    addi t2, t2, 4
    j    write_py
done_py:

    la   t1, pz
    li   t2, 0
write_pz:
    bge  t2, t3, done_pz
    add  t4, t1, t2
    lw   t5, 0(t4)
    add  t6, s1, t2
    li   t0, arr_bytes
    slli t0, t0, 1
    add  t6, t6, t0
    addi t6, t6, 4
    sw   t5, 0(t6)
    addi t2, t2, 4
    j    write_pz
done_pz:

    sw   zero, 0(s1)    # done writing, python can read

    j main_loop

exit_err:
    li   a7, 93
    li   a0, 1
    ecall
    j .

exit_ok:
    li   a7, 93
    li   a0, 0
    ecall


# v += 0.5 * a * dt
calc_kick:
    addi sp, sp, -48
    sw   ra,  0(sp)
    sw   s2,  4(sp)
    sw   s3,  8(sp)
    sw   s4, 12(sp)
    sw   s5, 16(sp)
    sw   s6, 20(sp)
    sw   s7, 24(sp)
    sw   s8, 28(sp)
    sw   s9, 32(sp)

    la   s2, vx
    la   s3, vy
    la   s4, vz
    la   s5, ax
    la   s6, ay
    la   s7, az
    li   s8, 0
    li   s9, n

    lui  t0, %hi(dt_f)
    addi t0, t0, %lo(dt_f)
    fmv.w.x ft0, t0
    lui  t0, %hi(half)
    addi t0, t0, %lo(half)
    fmv.w.x ft1, t0
    fmul.s ft2, ft0, ft1    # half_dt

kick_loop:
    bge  s8, s9, kick_done
    slli t0, s8, 2

    add  t1, s2, t0
    flw  ft3, 0(t1)
    add  t2, s5, t0
    flw  ft4, 0(t2)
    fmadd.s ft3, ft4, ft2, ft3
    fsw  ft3, 0(t1)

    add  t1, s3, t0
    flw  ft3, 0(t1)
    add  t2, s6, t0
    flw  ft4, 0(t2)
    fmadd.s ft3, ft4, ft2, ft3
    fsw  ft3, 0(t1)

    add  t1, s4, t0
    flw  ft3, 0(t1)
    add  t2, s7, t0
    flw  ft4, 0(t2)
    fmadd.s ft3, ft4, ft2, ft3
    fsw  ft3, 0(t1)

    addi s8, s8, 1
    j    kick_loop
kick_done:
    lw   ra,  0(sp)
    lw   s2,  4(sp)
    lw   s3,  8(sp)
    lw   s4, 12(sp)
    lw   s5, 16(sp)
    lw   s6, 20(sp)
    lw   s7, 24(sp)
    lw   s8, 28(sp)
    lw   s9, 32(sp)
    addi sp, sp, 48
    ret


# p += v * dt
calc_drift:
    addi sp, sp, -48
    sw   ra,  0(sp)
    sw   s2,  4(sp)
    sw   s3,  8(sp)
    sw   s4, 12(sp)
    sw   s5, 16(sp)
    sw   s6, 20(sp)
    sw   s7, 24(sp)
    sw   s8, 28(sp)
    sw   s9, 32(sp)

    la   s2, px
    la   s3, py
    la   s4, pz
    la   s5, vx
    la   s6, vy
    la   s7, vz
    li   s8, 0
    li   s9, n

    lui  t0, %hi(dt_f)
    addi t0, t0, %lo(dt_f)
    fmv.w.x ft0, t0

drift_loop:
    bge  s8, s9, drift_done
    slli t0, s8, 2

    add  t1, s2, t0
    flw  ft1, 0(t1)
    add  t2, s5, t0
    flw  ft2, 0(t2)
    fmadd.s ft1, ft2, ft0, ft1  # px += vx*dt
    fsw  ft1, 0(t1)

    add  t1, s3, t0
    flw  ft1, 0(t1)
    add  t2, s6, t0
    flw  ft2, 0(t2)
    fmadd.s ft1, ft2, ft0, ft1
    fsw  ft1, 0(t1)

    add  t1, s4, t0
    flw  ft1, 0(t1)
    add  t2, s7, t0
    flw  ft2, 0(t2)
    fmadd.s ft1, ft2, ft0, ft1
    fsw  ft1, 0(t1)

    addi s8, s8, 1
    j    drift_loop
drift_done:
    lw   ra,  0(sp)
    lw   s2,  4(sp)
    lw   s3,  8(sp)
    lw   s4, 12(sp)
    lw   s5, 16(sp)
    lw   s6, 20(sp)
    lw   s7, 24(sp)
    lw   s8, 28(sp)
    lw   s9, 32(sp)
    addi sp, sp, 48
    ret


# naive n^2 force calculation
calc_acc:
    addi sp, sp, -48
    sw   ra,  0(sp)
    sw   s0,  4(sp)
    sw   s1,  8(sp)
    sw   s2, 12(sp)
    sw   s3, 16(sp)
    sw   s4, 20(sp)
    sw   s5, 24(sp)
    sw   s6, 28(sp)
    sw   s7, 32(sp)
    sw   s8, 36(sp)
    sw   s9, 40(sp)

    la   s2, ax
    la   s3, ay
    la   s4, az
    li   s8, 0
    li   s9, n

    fcvt.s.w ft11, zero     # 0.0 for zeroing arrays

zero_acc:
    bge  s8, s9, zero_done
    slli t0, s8, 2
    add  t1, s2, t0
    fsw  ft11, 0(t1)
    add  t1, s3, t0
    fsw  ft11, 0(t1)
    add  t1, s4, t0
    fsw  ft11, 0(t1)
    addi s8, s8, 1
    j    zero_acc
zero_done:

    li   s8, 0

outer_loop:
    bge  s8, s9, acc_done
    slli t0, s8, 2

    la   t1, px
    add  t1, t1, t0
    flw  fs0, 0(t1)
    la   t1, py
    add  t1, t1, t0
    flw  fs1, 0(t1)
    la   t1, pz
    add  t1, t1, t0
    flw  fs2, 0(t1)

    la   s5, ax
    add  s5, s5, t0
    flw  fs3, 0(s5)
    la   s6, ay
    add  s6, s6, t0
    flw  fs4, 0(s6)
    la   s7, az
    add  s7, s7, t0
    flw  fs5, 0(s7)

    li   s0, 0

inner_loop:
    bge  s0, s9, inner_done
    beq  s0, s8, skip_self
    slli t0, s0, 2

    la   t1, px
    add  t1, t1, t0
    flw  ft0, 0(t1)
    fsub.s ft0, ft0, fs0    # dx

    la   t1, py
    add  t1, t1, t0
    flw  ft1, 0(t1)
    fsub.s ft1, ft1, fs1    # dy

    la   t1, pz
    add  t1, t1, t0
    flw  ft2, 0(t1)
    fsub.s ft2, ft2, fs2    # dz

    fmul.s ft3, ft0, ft0
    fmadd.s ft3, ft1, ft1, ft3
    fmadd.s ft3, ft2, ft2, ft3  # r^2

    lui  t1, %hi(soft_f)
    addi t1, t1, %lo(soft_f)
    fmv.w.x ft4, t1
    fmadd.s ft3, ft4, ft4, ft3  # r^2 + soft^2

    fsqrt.s ft5, ft3
    fmul.s  ft6, ft5, ft3       # r_soft^3

    la   t1, mass
    add  t1, t1, t0
    flw  ft7, 0(t1)

    lui  t1, %hi(G_f)
    addi t1, t1, %lo(G_f)
    fmv.w.x ft8, t1
    fmul.s  ft7, ft8, ft7
    fdiv.s  ft7, ft7, ft6       # G*mj/r^3

    fmadd.s fs3, ft7, ft0, fs3
    fmadd.s fs4, ft7, ft1, fs4
    fmadd.s fs5, ft7, ft2, fs5

skip_self:
    addi s0, s0, 1
    j    inner_loop

inner_done:
    fsw  fs3, 0(s5)
    fsw  fs4, 0(s6)
    fsw  fs5, 0(s7)

    addi s8, s8, 1
    j    outer_loop

acc_done:
    lw   ra,  0(sp)
    lw   s0,  4(sp)
    lw   s1,  8(sp)
    lw   s2, 12(sp)
    lw   s3, 16(sp)
    lw   s4, 20(sp)
    lw   s5, 24(sp)
    lw   s6, 28(sp)
    lw   s7, 32(sp)
    lw   s8, 36(sp)
    lw   s9, 40(sp)
    addi sp, sp, 48
    ret

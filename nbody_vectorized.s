# nbody_vectorized.s
# Vectorized RISC-V N-body naive simulation using RVV
#
# kick_step   : half-update velocities  (vectorized)
# drift_step  : update positions        (vectorized)
# accel_step  : compute accelerations   (vectorized, bitmask for i!=j)
#
# All float arrays are f32, 4-byte aligned.

    .section .data
    .align 2

    .equ N, 1000

    px:       .rept N  .float 0.0  .endr
    py:       .rept N  .float 0.0  .endr
    pz:       .rept N  .float 0.0  .endr

    vx:       .rept N  .float 0.0  .endr
    vy:       .rept N  .float 0.0  .endr
    vz:       .rept N  .float 0.0  .endr

    Ax:       .rept N  .float 0.0  .endr
    Ay:       .rept N  .float 0.0  .endr
    Az:       .rept N  .float 0.0  .endr

    mass:     .rept N  .float 1.0  .endr

    dt_val:   .float 0.01
    softening:.float 0.01

    .section .text
    .globl _start
    .globl kick_step
    .globl drift_step
    .globl accel_step

# kick_step(vx, vy, vz, Ax, Ay, Az, N, dt)
#
# vx[i] += Ax[i] * (dt * 0.5)  for all i  (likewise vy, vz)
#
# a0=&vx  a1=&vy  a2=&vz  a3=&Ax  a4=&Ay  a5=&Az  a6=N  fa0=dt
kick_step:
    fmv.s    ft1, fa0
    lui      t5, %hi(half_dt_const)
    flw      ft2, %lo(half_dt_const)(t5)
    fmul.s   ft1, ft1, ft2          # ft1 = half_dt

    mv       t0, a6
    mv       t1, a0;  mv t2, a1;  mv t3, a2
    mv       s1, a3;  mv s2, a4;  mv s3, a5

kick_loop:
    beqz     t0, kick_done
    vsetvli  t4, t0, e32, m1, ta, ma

    vle32.v  v0, (t1);  vle32.v  v1, (s1);  vfmacc.vf v0, ft1, v1;  vse32.v v0, (t1)
    vle32.v  v0, (t2);  vle32.v  v1, (s2);  vfmacc.vf v0, ft1, v1;  vse32.v v0, (t2)
    vle32.v  v0, (t3);  vle32.v  v1, (s3);  vfmacc.vf v0, ft1, v1;  vse32.v v0, (t3)

    slli     t5, t4, 2
    add      t1, t1, t5;  add t2, t2, t5;  add t3, t3, t5
    add      s1, s1, t5;  add s2, s2, t5;  add s3, s3, t5

    sub      t0, t0, t4
    j        kick_loop

kick_done:
    ret

    .section .rodata
    half_dt_const: .float 0.5

    .section .text

# drift_step(px, py, pz, vx, vy, vz, N, dt)
#
# px[i] += vx[i] * dt  for all i  (likewise py, pz)
#
# a0=&px  a1=&py  a2=&pz  a3=&vx  a4=&vy  a5=&vz  a6=N  fa0=dt
drift_step:
    fmv.s    ft1, fa0

    mv       t0, a6
    mv       t1, a0;  mv t2, a1;  mv t3, a2
    mv       s1, a3;  mv s2, a4;  mv s3, a5

drift_loop:
    beqz     t0, drift_done
    vsetvli  t4, t0, e32, m1, ta, ma

    vle32.v  v0, (t1);  vle32.v  v1, (s1);  vfmacc.vf v0, ft1, v1;  vse32.v v0, (t1)
    vle32.v  v0, (t2);  vle32.v  v1, (s2);  vfmacc.vf v0, ft1, v1;  vse32.v v0, (t2)
    vle32.v  v0, (t3);  vle32.v  v1, (s3);  vfmacc.vf v0, ft1, v1;  vse32.v v0, (t3)

    slli     t5, t4, 2
    add      t1, t1, t5;  add t2, t2, t5;  add t3, t3, t5
    add      s1, s1, t5;  add s2, s2, t5;  add s3, s3, t5

    sub      t0, t0, t4
    j        drift_loop

drift_done:
    ret

# accel_step(px, py, pz, mass, Ax, Ay, Az, N, softening)
#
# For each i, accumulate gravitational acceleration from all j != i.
# The inner loop over j is vectorized; i==j is masked out with vmsne.
#
# a0=&px  a1=&py  a2=&pz  a3=&mass  a4=&Ax  a5=&Ay  a6=&Az  a7=N  fa0=softening
accel_step:
    addi     sp, sp, -48
    sw       ra, 44(sp)
    sw       s0, 40(sp);  sw s1, 36(sp);  sw s2, 32(sp);  sw s3, 28(sp)
    sw       s4, 24(sp);  sw s5, 20(sp);  sw s6, 16(sp);  sw s7, 12(sp)

    mv       s0, a0;  mv s1, a1;  mv s2, a2;  mv s3, a3
    mv       s4, a4;  mv s5, a5;  mv s6, a6;  mv s7, a7
    fmv.s    fs0, fa0

    # zero Ax/Ay/Az
    mv       t0, s7;  mv t1, s4;  mv t2, s5;  mv t3, s6
zero_accel_loop:
    beqz     t0, zero_accel_done
    vsetvli  t4, t0, e32, m1, ta, ma
    vmv.v.i  v0, 0
    vse32.v  v0, (t1);  vse32.v v0, (t2);  vse32.v v0, (t3)
    slli     t5, t4, 2
    add      t1, t1, t5;  add t2, t2, t5;  add t3, t3, t5
    sub      t0, t0, t4
    j        zero_accel_loop
zero_accel_done:

    li       s8, 0                  # i = 0

outer_loop:
    bge      s8, s7, accel_done

    slli     t0, s8, 2
    add      t0, s0, t0;  flw ft0, 0(t0)   # px[i]
    slli     t0, s8, 2
    add      t0, s1, t0;  flw ft1, 0(t0)   # py[i]
    slli     t0, s8, 2
    add      t0, s2, t0;  flw ft2, 0(t0)   # pz[i]

    fmv.s.x  ft3, zero              # acc_x
    fmv.s.x  ft4, zero              # acc_y
    fmv.s.x  ft5, zero              # acc_z

    li       t6, 0                  # j_base
    mv       t0, s7

inner_loop:
    beqz     t0, inner_done
    vsetvli  t4, t0, e32, m1, ta, ma

    vid.v    v8
    vmv.v.x  v9, t6
    vadd.vv  v8, v8, v9             # v8 = [j_base, j_base+1, ...]

    vmv.v.x  v10, s8
    vmsne.vv v0, v8, v10            # mask: j != i

    vsll.vi  v11, v8, 2             # byte offsets

    vloxei32.v v2, (s0), v11        # px[j]
    vloxei32.v v3, (s1), v11        # py[j]
    vloxei32.v v4, (s2), v11        # pz[j]
    vloxei32.v v5, (s3), v11        # mass[j]

    vfsub.vf v2, v2, ft0            # dx
    vfsub.vf v3, v3, ft1            # dy
    vfsub.vf v4, v4, ft2            # dz

    vfmul.vv  v6, v2, v2
    vfmacc.vv v6, v3, v3
    vfmacc.vv v6, v4, v4
    vfadd.vf  v6, v6, fs0           # dist_sq

    vfsqrt.v  v7, v6
    lui       t5, %hi(one_float)
    flw       ft6, %lo(one_float)(t5)
    vfrdiv.vf v7, v7, ft6           # inv_dist

    vfmul.vv  v6, v7, v7
    vfmul.vv  v6, v6, v7            # inv_dist^3

    vfmul.vv  v6, v6, v5            # mass[j] * inv_dist^3

    vmv.v.i   v12, 0
    vmerge.vvm v6, v12, v6, v0      # zero the i==j lane

    vfmul.vv  v13, v2, v6
    vfmul.vv  v14, v3, v6
    vfmul.vv  v15, v4, v6

    vfredusum.vs v16, v13, v0_zero
    vfmv.f.s     ft6, v16;  fadd.s ft3, ft3, ft6

    vfredusum.vs v16, v14, v0_zero
    vfmv.f.s     ft6, v16;  fadd.s ft4, ft4, ft6

    vfredusum.vs v16, v15, v0_zero
    vfmv.f.s     ft6, v16;  fadd.s ft5, ft5, ft6

    add      t6, t6, t4
    sub      t0, t0, t4
    j        inner_loop

inner_done:
    slli     t0, s8, 2
    add      t1, s4, t0;  fsw ft3, 0(t1)
    add      t1, s5, t0;  fsw ft4, 0(t1)
    add      t1, s6, t0;  fsw ft5, 0(t1)

    addi     s8, s8, 1
    j        outer_loop

accel_done:
    lw       ra, 44(sp)
    lw       s0, 40(sp);  lw s1, 36(sp);  lw s2, 32(sp);  lw s3, 28(sp)
    lw       s4, 24(sp);  lw s5, 20(sp);  lw s6, 16(sp);  lw s7, 12(sp)
    addi     sp, sp, 48
    ret

    .section .rodata
    one_float: .float 1.0

    .section .text

_start:
    la       a0, vx;  la a1, vy;  la a2, vz
    la       a3, Ax;  la a4, Ay;  la a5, Az
    li       a6, N
    la       t0, dt_val;  flw fa0, 0(t0)
    call     kick_step

    la       a0, px;  la a1, py;  la a2, pz
    la       a3, vx;  la a4, vy;  la a5, vz
    li       a6, N
    la       t0, dt_val;  flw fa0, 0(t0)
    call     drift_step

    la       a0, px;  la a1, py;  la a2, pz
    la       a3, mass
    la       a4, Ax;  la a5, Ay;  la a6, Az
    li       a7, N
    la       t0, softening;  flw fa0, 0(t0)
    call     accel_step

    la       a0, vx;  la a1, vy;  la a2, vz
    la       a3, Ax;  la a4, Ay;  la a5, Az
    li       a6, N
    la       t0, dt_val;  flw fa0, 0(t0)
    call     kick_step

    li       a7, 93
    li       a0, 0
    ecall

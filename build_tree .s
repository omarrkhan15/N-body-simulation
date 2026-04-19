# build_tree.s
# RISC-V scalar — Barnes-Hut octree build phase
#
# Functions: reset_tree, allocate_node, get_octant,
#            set_child_bbox, insert_particle, build_tree
#
# Child k of node i lives at node_child[i*8 + k].

    .section .data
    .align 2

    .equ  N,         100
    .equ  MAX_NODES, 900        # 8*N + 100

    .globl px
    .globl py
    .globl pz
    .globl mass

    px:   .rept N  .float 0.0  .endr
    py:   .rept N  .float 0.0  .endr
    pz:   .rept N  .float 0.0  .endr
    mass: .rept N  .float 1.0  .endr

    .globl node_cx
    .globl node_cy
    .globl node_cz
    .globl node_mass
    .globl node_child
    .globl node_particle
    .globl node_is_leaf
    .globl node_xmin
    .globl node_xmax
    .globl node_ymin
    .globl node_ymax
    .globl node_zmin
    .globl node_zmax
    .globl node_count

    node_cx:       .rept MAX_NODES  .float 0.0  .endr
    node_cy:       .rept MAX_NODES  .float 0.0  .endr
    node_cz:       .rept MAX_NODES  .float 0.0  .endr
    node_mass:     .rept MAX_NODES  .float 0.0  .endr

    node_child:
        .rept MAX_NODES * 8
        .word -1
        .endr

    node_particle: .rept MAX_NODES  .word -1  .endr
    node_is_leaf:  .rept MAX_NODES  .word  1  .endr

    node_xmin:     .rept MAX_NODES  .float 0.0  .endr
    node_xmax:     .rept MAX_NODES  .float 0.0  .endr
    node_ymin:     .rept MAX_NODES  .float 0.0  .endr
    node_ymax:     .rept MAX_NODES  .float 0.0  .endr
    node_zmin:     .rept MAX_NODES  .float 0.0  .endr
    node_zmax:     .rept MAX_NODES  .float 0.0  .endr

    node_count:    .word 0

    float_half:    .float 0.5
    float_zero:    .float 0.0
    int_neg1:      .word  -1
    int_max_nodes: .word  MAX_NODES

    .equ  STACK_SIZE, 64
    ip_stack:      .rept STACK_SIZE  .word 0  .endr
    ip_stack_part: .rept STACK_SIZE  .word 0  .endr

    .section .text
    .globl allocate_node
    .globl get_octant
    .globl set_child_bbox
    .globl insert_particle
    .globl build_tree
    .globl reset_tree

# reset_tree()
# Zero the node pool so the tree can be rebuilt from scratch.
reset_tree:
    addi  sp, sp, -8
    sw    ra, 4(sp)
    sw    s0, 0(sp)

    la    t0, node_count
    sw    zero, 0(t0)

    li    s0, 0
    li    t6, MAX_NODES

reset_loop:
    bge   s0, t6, reset_done

    slli  t0, s0, 2

    la    t1, node_cx;       add t1, t1, t0;  sw zero, 0(t1)
    la    t1, node_cy;       add t1, t1, t0;  sw zero, 0(t1)
    la    t1, node_cz;       add t1, t1, t0;  sw zero, 0(t1)
    la    t1, node_mass;     add t1, t1, t0;  sw zero, 0(t1)

    la    t1, node_particle; add t1, t1, t0;  li t2, -1;  sw t2, 0(t1)
    la    t1, node_is_leaf;  add t1, t1, t0;  li t2,  1;  sw t2, 0(t1)

    la    t1, node_xmin;  add t1, t1, t0;  sw zero, 0(t1)
    la    t1, node_xmax;  add t1, t1, t0;  sw zero, 0(t1)
    la    t1, node_ymin;  add t1, t1, t0;  sw zero, 0(t1)
    la    t1, node_ymax;  add t1, t1, t0;  sw zero, 0(t1)
    la    t1, node_zmin;  add t1, t1, t0;  sw zero, 0(t1)
    la    t1, node_zmax;  add t1, t1, t0;  sw zero, 0(t1)

    slli  t0, s0, 5         # i * 32  (8 children * 4 bytes)
    la    t1, node_child
    add   t1, t1, t0
    li    t2, -1
    sw    t2,  0(t1)
    sw    t2,  4(t1)
    sw    t2,  8(t1)
    sw    t2, 12(t1)
    sw    t2, 16(t1)
    sw    t2, 20(t1)
    sw    t2, 24(t1)
    sw    t2, 28(t1)

    addi  s0, s0, 1
    j     reset_loop

reset_done:
    lw    ra, 4(sp)
    lw    s0, 0(sp)
    addi  sp, sp, 8
    ret

# allocate_node() -> a0 = new node index
allocate_node:
    la    t0, node_count
    lw    a0, 0(t0)

    li    t1, MAX_NODES
    blt   a0, t1, alloc_ok
    li    a7, 93
    li    a0, 1
    ecall

alloc_ok:
    addi  t1, a0, 1
    sw    t1, 0(t0)
    ret

# get_octant(a0=node_idx, a1=part_idx) -> a0 = octant (0-7)
#
# Bit encoding: bit0=x, bit1=y, bit2=z (set if particle >= midpoint)
get_octant:
    slli  t0, a0, 2

    la    t1, node_xmin;  add t1, t1, t0;  flw ft0, 0(t1)
    la    t1, node_xmax;  add t1, t1, t0;  flw ft1, 0(t1)
    la    t1, node_ymin;  add t1, t1, t0;  flw ft2, 0(t1)
    la    t1, node_ymax;  add t1, t1, t0;  flw ft3, 0(t1)
    la    t1, node_zmin;  add t1, t1, t0;  flw ft4, 0(t1)
    la    t1, node_zmax;  add t1, t1, t0;  flw ft5, 0(t1)

    la    t2, float_half
    flw   ft6, 0(t2)

    fadd.s ft0, ft0, ft1;  fmul.s ft0, ft0, ft6   # mid_x
    fadd.s ft2, ft2, ft3;  fmul.s ft2, ft2, ft6   # mid_y
    fadd.s ft4, ft4, ft5;  fmul.s ft4, ft4, ft6   # mid_z

    slli  t0, a1, 2
    la    t1, px;  add t1, t1, t0;  flw ft1, 0(t1)
    la    t1, py;  add t1, t1, t0;  flw ft3, 0(t1)
    la    t1, pz;  add t1, t1, t0;  flw ft5, 0(t1)

    li    a0, 0

    flt.s t0, ft0, ft1;  or a0, a0, t0          # bit 0
    flt.s t0, ft2, ft3;  slli t0, t0, 1;  or a0, a0, t0   # bit 1
    flt.s t0, ft4, ft5;  slli t0, t0, 2;  or a0, a0, t0   # bit 2

    ret

# set_child_bbox(a0=parent, a1=child, a2=octant)
set_child_bbox:
    addi  sp, sp, -16
    sw    ra,  12(sp)
    sw    s1,   8(sp)
    sw    s2,   4(sp)
    sw    s3,   0(sp)

    mv    s1, a0
    mv    s2, a1
    mv    s3, a2

    slli  t0, s1, 2
    la    t1, node_xmin;  add t1, t1, t0;  flw ft0, 0(t1)
    la    t1, node_xmax;  add t1, t1, t0;  flw ft1, 0(t1)
    la    t1, node_ymin;  add t1, t1, t0;  flw ft2, 0(t1)
    la    t1, node_ymax;  add t1, t1, t0;  flw ft3, 0(t1)
    la    t1, node_zmin;  add t1, t1, t0;  flw ft4, 0(t1)
    la    t1, node_zmax;  add t1, t1, t0;  flw ft5, 0(t1)

    la    t2, float_half
    flw   ft6, 0(t2)

    fadd.s ft7, ft0, ft1;  fmul.s ft7, ft7, ft6
    addi  sp, sp, -12
    fsw   ft7, 8(sp)        # mid_x

    fadd.s ft7, ft2, ft3;  fmul.s ft7, ft7, ft6
    fsw   ft7, 4(sp)        # mid_y

    fadd.s ft7, ft4, ft5;  fmul.s ft7, ft7, ft6
    fsw   ft7, 0(sp)        # mid_z

    slli  t0, s2, 2

    andi  t1, s3, 1
    beqz  t1, set_bbox_x_lo
    flw   ft7, 8(sp)
    la    t2, node_xmin;  add t2, t2, t0;  fsw ft7, 0(t2)
    la    t2, node_xmax;  add t2, t2, t0;  fsw ft1, 0(t2)
    j     set_bbox_y

set_bbox_x_lo:
    la    t2, node_xmin;  add t2, t2, t0;  fsw ft0, 0(t2)
    flw   ft7, 8(sp)
    la    t2, node_xmax;  add t2, t2, t0;  fsw ft7, 0(t2)

set_bbox_y:
    andi  t1, s3, 2
    beqz  t1, set_bbox_y_lo
    flw   ft7, 4(sp)
    la    t2, node_ymin;  add t2, t2, t0;  fsw ft7, 0(t2)
    la    t2, node_ymax;  add t2, t2, t0;  fsw ft3, 0(t2)
    j     set_bbox_z

set_bbox_y_lo:
    la    t2, node_ymin;  add t2, t2, t0;  fsw ft2, 0(t2)
    flw   ft7, 4(sp)
    la    t2, node_ymax;  add t2, t2, t0;  fsw ft7, 0(t2)

set_bbox_z:
    andi  t1, s3, 4
    beqz  t1, set_bbox_z_lo
    flw   ft7, 0(sp)
    la    t2, node_zmin;  add t2, t2, t0;  fsw ft7, 0(t2)
    la    t2, node_zmax;  add t2, t2, t0;  fsw ft5, 0(t2)
    j     set_bbox_done

set_bbox_z_lo:
    la    t2, node_zmin;  add t2, t2, t0;  fsw ft4, 0(t2)
    flw   ft7, 0(sp)
    la    t2, node_zmax;  add t2, t2, t0;  fsw ft7, 0(t2)

set_bbox_done:
    addi  sp, sp, 12
    lw    ra,  12(sp)
    lw    s1,   8(sp)
    lw    s2,   4(sp)
    lw    s3,   0(sp)
    addi  sp, sp, 16
    ret

# insert_particle(a0=node_idx, a1=part_idx)
#
# Iterative insertion using a pre-allocated explicit stack in .data.
# s0 = explicit stack pointer (index, not hardware sp)
# s1 = current node, s2 = current particle
insert_particle:
    addi  sp, sp, -36
    sw    ra,  32(sp)
    sw    s0,  28(sp)
    sw    s1,  24(sp)
    sw    s2,  20(sp)
    sw    s3,  16(sp)
    sw    s4,  12(sp)
    sw    s5,   8(sp)
    sw    s6,   4(sp)
    sw    s7,   0(sp)

    li    s0, 0
    la    t0, ip_stack;      sw a0, 0(t0)
    la    t0, ip_stack_part; sw a1, 0(t0)
    addi  s0, s0, 1

ip_main_loop:
    beqz  s0, ip_done

    addi  s0, s0, -1
    slli  t0, s0, 2

    la    t1, ip_stack;      add t1, t1, t0;  lw s1, 0(t1)
    la    t1, ip_stack_part; add t1, t1, t0;  lw s2, 0(t1)

    slli  t0, s1, 2
    la    t1, node_is_leaf
    add   t1, t1, t0
    lw    t2, 0(t1)

    beqz  t2, ip_internal

    # --- leaf ---
    la    t1, node_particle
    add   t1, t1, t0
    lw    s3, 0(t1)

    li    t2, -1
    beq   s3, t2, ip_empty_leaf

    # occupied leaf: mark internal, push both particles back
    la    t1, node_is_leaf;  slli t0, s1, 2;  add t1, t1, t0;  sw zero, 0(t1)
    la    t1, node_particle; add t1, t1, t0;  li t2, -1;  sw t2, 0(t1)

    slli  t0, s0, 2
    la    t1, ip_stack;      add t1, t1, t0;  sw s1, 0(t1)
    la    t1, ip_stack_part; add t1, t1, t0;  sw s3, 0(t1)
    addi  s0, s0, 1

    slli  t0, s0, 2
    la    t1, ip_stack;      add t1, t1, t0;  sw s1, 0(t1)
    la    t1, ip_stack_part; add t1, t1, t0;  sw s2, 0(t1)
    addi  s0, s0, 1

    j     ip_main_loop

ip_empty_leaf:
    slli  t0, s1, 2
    la    t1, node_particle; add t1, t1, t0;  sw s2, 0(t1)

    slli  t0, s2, 2
    la    t1, px;   add t1, t1, t0;  flw ft0, 0(t1)
    la    t1, py;   add t1, t1, t0;  flw ft1, 0(t1)
    la    t1, pz;   add t1, t1, t0;  flw ft2, 0(t1)
    la    t1, mass; add t1, t1, t0;  flw ft3, 0(t1)

    slli  t0, s1, 2
    la    t1, node_cx;   add t1, t1, t0;  fsw ft0, 0(t1)
    la    t1, node_cy;   add t1, t1, t0;  fsw ft1, 0(t1)
    la    t1, node_cz;   add t1, t1, t0;  fsw ft2, 0(t1)
    la    t1, node_mass; add t1, t1, t0;  fsw ft3, 0(t1)

    j     ip_main_loop

    # --- internal node ---
ip_internal:
    mv    a0, s1
    mv    a1, s2
    call  get_octant
    mv    s3, a0

    slli  t0, s1, 3;  add t0, t0, s3;  slli t0, t0, 2
    la    t1, node_child;  add t1, t1, t0
    lw    s4, 0(t1)

    li    t2, -1
    bne   s4, t2, ip_child_exists

    call  allocate_node
    mv    s4, a0

    mv    a0, s1;  mv a1, s4;  mv a2, s3
    call  set_child_bbox

    slli  t0, s1, 3;  add t0, t0, s3;  slli t0, t0, 2
    la    t1, node_child;  add t1, t1, t0
    sw    s4, 0(t1)

ip_child_exists:
    slli  t0, s0, 2
    la    t1, ip_stack;      add t1, t1, t0;  sw s4, 0(t1)
    la    t1, ip_stack_part; add t1, t1, t0;  sw s2, 0(t1)
    addi  s0, s0, 1

    # incremental CoM update for this internal node
    slli  t0, s1, 2
    la    t1, node_cx;   add t1, t1, t0;  flw ft1, 0(t1)
    la    t1, node_cy;   add t1, t1, t0;  flw ft2, 0(t1)
    la    t1, node_cz;   add t1, t1, t0;  flw ft3, 0(t1)
    la    t1, node_mass; add t1, t1, t0;  flw ft0, 0(t1)

    slli  t0, s2, 2
    la    t1, mass; add t1, t1, t0;  flw ft4, 0(t1)
    la    t1, px;   add t1, t1, t0;  flw ft5, 0(t1)
    la    t1, py;   add t1, t1, t0;  flw ft6, 0(t1)
    la    t1, pz;   add t1, t1, t0;  flw ft7, 0(t1)

    fadd.s  fa0, ft0, ft4

    fmul.s  fa1, ft1, ft0;  fmul.s fa2, ft5, ft4;  fadd.s fa1, fa1, fa2;  fdiv.s fa1, fa1, fa0
    fmul.s  fa2, ft2, ft0;  fmul.s fa3, ft6, ft4;  fadd.s fa2, fa2, fa3;  fdiv.s fa2, fa2, fa0
    fmul.s  fa3, ft3, ft0;  fmul.s fa4, ft7, ft4;  fadd.s fa3, fa3, fa4;  fdiv.s fa3, fa3, fa0

    slli  t0, s1, 2
    la    t1, node_cx;   add t1, t1, t0;  fsw fa1, 0(t1)
    la    t1, node_cy;   add t1, t1, t0;  fsw fa2, 0(t1)
    la    t1, node_cz;   add t1, t1, t0;  fsw fa3, 0(t1)
    la    t1, node_mass; add t1, t1, t0;  fsw fa0, 0(t1)

    j     ip_main_loop

ip_done:
    lw    ra,  32(sp)
    lw    s0,  28(sp)
    lw    s1,  24(sp)
    lw    s2,  20(sp)
    lw    s3,  16(sp)
    lw    s4,  12(sp)
    lw    s5,   8(sp)
    lw    s6,   4(sp)
    lw    s7,   0(sp)
    addi  sp, sp, 36
    ret

# build_tree() -> a0 = root node index
build_tree:
    addi  sp, sp, -16
    sw    ra, 12(sp)
    sw    s0,  8(sp)
    sw    s1,  4(sp)
    sw    s2,  0(sp)

    call  reset_tree
    call  allocate_node
    mv    s0, a0               # root index

    # scan px/py/pz to find axis-aligned bounding box
    la    t0, px;  flw fs0, 0(t0);  fmv.s fs1, fs0
    la    t0, py;  flw fs2, 0(t0);  fmv.s fs3, fs2
    la    t0, pz;  flw fs4, 0(t0);  fmv.s fs5, fs4

    li    t0, 1
    li    t6, N

scan_loop:
    bge   t0, t6, scan_done
    slli  t1, t0, 2

    la    t2, px;  add t2, t2, t1;  flw ft0, 0(t2)
    fmin.s fs0, fs0, ft0;  fmax.s fs1, fs1, ft0

    la    t2, py;  add t2, t2, t1;  flw ft0, 0(t2)
    fmin.s fs2, fs2, ft0;  fmax.s fs3, fs3, ft0

    la    t2, pz;  add t2, t2, t1;  flw ft0, 0(t2)
    fmin.s fs4, fs4, ft0;  fmax.s fs5, fs5, ft0

    addi  t0, t0, 1
    j     scan_loop

scan_done:
    la    t0, bbox_pad
    flw   ft1, 0(t0)

    fsub.s fs0, fs0, ft1;  fadd.s fs1, fs1, ft1
    fsub.s fs2, fs2, ft1;  fadd.s fs3, fs3, ft1
    fsub.s fs4, fs4, ft1;  fadd.s fs5, fs5, ft1

    slli  t0, s0, 2
    la    t1, node_xmin;  add t1, t1, t0;  fsw fs0, 0(t1)
    la    t1, node_xmax;  add t1, t1, t0;  fsw fs1, 0(t1)
    la    t1, node_ymin;  add t1, t1, t0;  fsw fs2, 0(t1)
    la    t1, node_ymax;  add t1, t1, t0;  fsw fs3, 0(t1)
    la    t1, node_zmin;  add t1, t1, t0;  fsw fs4, 0(t1)
    la    t1, node_zmax;  add t1, t1, t0;  fsw fs5, 0(t1)

    li    s1, 0
    li    t6, N

insert_loop:
    bge   s1, t6, insert_done
    mv    a0, s0
    mv    a1, s1
    call  insert_particle
    addi  s1, s1, 1
    j     insert_loop

insert_done:
    mv    a0, s0

    lw    ra, 12(sp)
    lw    s0,  8(sp)
    lw    s1,  4(sp)
    lw    s2,  0(sp)
    addi  sp, sp, 16
    ret

    .section .rodata
    bbox_pad: .float 0.0001

    .section .text

    .globl _start
_start:
    call  build_tree

    la    t0, node_mass
    flw   ft0, 0(t0)

    la    t0, node_count
    lw    t1, 0(t0)

    li    a7, 93
    li    a0, 0
    ecall

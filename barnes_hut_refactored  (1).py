import math
import random

N         = 100
MAX_NODES = 8 * N + 100
DT        = 0.01
SOFTENING = 0.01
THETA     = 0.5

px   = [0.0] * N
py   = [0.0] * N
pz   = [0.0] * N
vx   = [0.0] * N
vy   = [0.0] * N
vz   = [0.0] * N
mass = [1.0] * N

node_cx       = [0.0] * MAX_NODES
node_cy       = [0.0] * MAX_NODES
node_cz       = [0.0] * MAX_NODES
node_mass     = [0.0] * MAX_NODES
node_child    = [-1]  * (MAX_NODES * 8)
node_particle = [-1]  * MAX_NODES
node_is_leaf  = [1]   * MAX_NODES
node_xmin     = [0.0] * MAX_NODES
node_xmax     = [0.0] * MAX_NODES
node_ymin     = [0.0] * MAX_NODES
node_ymax     = [0.0] * MAX_NODES
node_zmin     = [0.0] * MAX_NODES
node_zmax     = [0.0] * MAX_NODES
node_count    = [0]


def reset_tree():
    for i in range(MAX_NODES):
        node_cx[i]       = 0.0
        node_cy[i]       = 0.0
        node_cz[i]       = 0.0
        node_mass[i]     = 0.0
        node_particle[i] = -1
        node_is_leaf[i]  = 1
        node_xmin[i]     = 0.0
        node_xmax[i]     = 0.0
        node_ymin[i]     = 0.0
        node_ymax[i]     = 0.0
        node_zmin[i]     = 0.0
        node_zmax[i]     = 0.0
        for k in range(8):
            node_child[i * 8 + k] = -1
    node_count[0] = 0


def allocate_node():
    idx = node_count[0]
    if idx >= MAX_NODES:
        raise RuntimeError("Node pool exhausted")
    node_count[0] += 1
    return idx


def get_octant(node_idx, part_idx):
    mid_x = 0.5 * (node_xmin[node_idx] + node_xmax[node_idx])
    mid_y = 0.5 * (node_ymin[node_idx] + node_ymax[node_idx])
    mid_z = 0.5 * (node_zmin[node_idx] + node_zmax[node_idx])

    octant = 0
    if px[part_idx] >= mid_x:
        octant |= 1
    if py[part_idx] >= mid_y:
        octant |= 2
    if pz[part_idx] >= mid_z:
        octant |= 4
    return octant


def set_child_bbox(parent, child_idx, octant):
    mid_x = 0.5 * (node_xmin[parent] + node_xmax[parent])
    mid_y = 0.5 * (node_ymin[parent] + node_ymax[parent])
    mid_z = 0.5 * (node_zmin[parent] + node_zmax[parent])

    if octant & 1:
        node_xmin[child_idx] = mid_x
        node_xmax[child_idx] = node_xmax[parent]
    else:
        node_xmin[child_idx] = node_xmin[parent]
        node_xmax[child_idx] = mid_x

    if octant & 2:
        node_ymin[child_idx] = mid_y
        node_ymax[child_idx] = node_ymax[parent]
    else:
        node_ymin[child_idx] = node_ymin[parent]
        node_ymax[child_idx] = mid_y

    if octant & 4:
        node_zmin[child_idx] = mid_z
        node_zmax[child_idx] = node_zmax[parent]
    else:
        node_zmin[child_idx] = node_zmin[parent]
        node_zmax[child_idx] = mid_z


def insert_particle(node_idx, part_idx):
    STACK_SIZE = 64
    stack      = [-1] * STACK_SIZE
    stack_part = [-1] * STACK_SIZE
    sp = 0

    stack[sp]      = node_idx
    stack_part[sp] = part_idx
    sp += 1

    while sp > 0:
        sp -= 1
        cur  = stack[sp]
        pidx = stack_part[sp]

        if node_is_leaf[cur]:
            existing = node_particle[cur]

            if existing == -1:
                node_particle[cur] = pidx
                node_cx[cur]       = px[pidx]
                node_cy[cur]       = py[pidx]
                node_cz[cur]       = pz[pidx]
                node_mass[cur]     = mass[pidx]
            else:
                node_is_leaf[cur]  = 0
                node_particle[cur] = -1

                stack[sp]      = cur
                stack_part[sp] = existing
                sp += 1

                stack[sp]      = cur
                stack_part[sp] = pidx
                sp += 1
        else:
            oct = get_octant(cur, pidx)
            child_slot = cur * 8 + oct
            child = node_child[child_slot]

            if child == -1:
                child = allocate_node()
                set_child_bbox(cur, child, oct)
                node_child[child_slot] = child

            stack[sp]      = child
            stack_part[sp] = pidx
            sp += 1


def propagate_com(node_idx):
    if node_idx == -1:
        return
    if node_is_leaf[node_idx]:
        return

    total_m = 0.0
    cx = cy = cz = 0.0
    for k in range(8):
        child = node_child[node_idx * 8 + k]
        if child == -1:
            continue
        propagate_com(child)
        m = node_mass[child]
        total_m += m
        cx += node_cx[child] * m
        cy += node_cy[child] * m
        cz += node_cz[child] * m

    node_mass[node_idx] = total_m
    if total_m > 0.0:
        node_cx[node_idx] = cx / total_m
        node_cy[node_idx] = cy / total_m
        node_cz[node_idx] = cz / total_m


def build_tree():
    reset_tree()
    root = allocate_node()

    xmin = min(px); xmax = max(px)
    ymin = min(py); ymax = max(py)
    zmin = min(pz); zmax = max(pz)

    pad = 1e-4
    node_xmin[root] = xmin - pad; node_xmax[root] = xmax + pad
    node_ymin[root] = ymin - pad; node_ymax[root] = ymax + pad
    node_zmin[root] = zmin - pad; node_zmax[root] = zmax + pad

    for i in range(N):
        insert_particle(root, i)

    propagate_com(root)
    return root


def compute_force_bh(i, node_idx):
    if node_idx == -1:
        return 0.0, 0.0, 0.0

    if node_is_leaf[node_idx]:
        pidx = node_particle[node_idx]
        if pidx == -1 or pidx == i:
            return 0.0, 0.0, 0.0
        dx = px[pidx] - px[i]
        dy = py[pidx] - py[i]
        dz = pz[pidx] - pz[i]
        dist_sq = dx*dx + dy*dy + dz*dz + SOFTENING
        inv_d   = 1.0 / math.sqrt(dist_sq)
        inv_d3  = inv_d ** 3
        f = mass[pidx] * inv_d3
        return dx*f, dy*f, dz*f

    node_size = node_xmax[node_idx] - node_xmin[node_idx]
    dx_com = node_cx[node_idx] - px[i]
    dy_com = node_cy[node_idx] - py[i]
    dz_com = node_cz[node_idx] - pz[i]
    dist = math.sqrt(dx_com*dx_com + dy_com*dy_com + dz_com*dz_com + SOFTENING)

    if node_size / dist < THETA:
        inv_d  = 1.0 / dist
        inv_d3 = inv_d ** 3
        f = node_mass[node_idx] * inv_d3
        return dx_com*f, dy_com*f, dz_com*f

    ax = ay = az = 0.0
    for k in range(8):
        child = node_child[node_idx * 8 + k]
        if child != -1:
            dax, day, daz = compute_force_bh(i, child)
            ax += dax; ay += day; az += daz
    return ax, ay, az


def accel_bh():
    root = build_tree()
    ax_arr = [0.0] * N
    ay_arr = [0.0] * N
    az_arr = [0.0] * N
    for i in range(N):
        ax_arr[i], ay_arr[i], az_arr[i] = compute_force_bh(i, root)
    return ax_arr, ay_arr, az_arr


def accel_naive():
    ax_arr = [0.0] * N
    ay_arr = [0.0] * N
    az_arr = [0.0] * N
    for i in range(N):
        for j in range(N):
            if i == j:
                continue
            dx = px[j] - px[i]
            dy = py[j] - py[i]
            dz = pz[j] - pz[i]
            dist_sq = dx*dx + dy*dy + dz*dz + SOFTENING
            inv_d   = 1.0 / math.sqrt(dist_sq)
            inv_d3  = inv_d ** 3
            f = mass[j] * inv_d3
            ax_arr[i] += dx * f
            ay_arr[i] += dy * f
            az_arr[i] += dz * f
    return ax_arr, ay_arr, az_arr


def step(use_bh=True):
    global px, py, pz, vx, vy, vz

    if use_bh:
        Ax, Ay, Az = accel_bh()
    else:
        Ax, Ay, Az = accel_naive()

    half_dt = 0.5 * DT

    for i in range(N):
        vx[i] += Ax[i] * half_dt
        vy[i] += Ay[i] * half_dt
        vz[i] += Az[i] * half_dt

    for i in range(N):
        px[i] += vx[i] * DT
        py[i] += vy[i] * DT
        pz[i] += vz[i] * DT

    if use_bh:
        Ax, Ay, Az = accel_bh()
    else:
        Ax, Ay, Az = accel_naive()

    for i in range(N):
        vx[i] += Ax[i] * half_dt
        vy[i] += Ay[i] * half_dt
        vz[i] += Az[i] * half_dt


def init_random(seed=42):
    random.seed(seed)
    for i in range(N):
        px[i]   = random.uniform(-1.0, 1.0)
        py[i]   = random.uniform(-1.0, 1.0)
        pz[i]   = random.uniform(-1.0, 1.0)
        vx[i]   = random.uniform(-0.1, 0.1)
        vy[i]   = random.uniform(-0.1, 0.1)
        vz[i]   = random.uniform(-0.1, 0.1)
        mass[i] = random.uniform(0.5, 1.5)


if __name__ == "__main__":
    init_random()
    root = build_tree()
    print(f"N={N}, nodes used={node_count[0]}")
    print(f"Root mass={node_mass[root]:.4f}, sum={sum(mass):.4f}")

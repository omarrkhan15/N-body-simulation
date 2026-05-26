import numpy as np
import time

# Tree parameters
N = 100
MAX_NODES = 50000  # sized large enough for any N used in benchmarks
THETA = 0.5        # opening criterion
SOFTENING = 0.01   # gravitational softening

# Particle arrays
px_arr = np.random.uniform(-1.0, 1.0, N).astype(np.float32)
py_arr = np.random.uniform(-1.0, 1.0, N).astype(np.float32)
pz_arr = np.random.uniform(-1.0, 1.0, N).astype(np.float32)
mass_arr = np.ones(N, dtype=np.float32)

# Node arrays (matches build_tree.s layout)
node_cx   = np.zeros(MAX_NODES, dtype=np.float32)
node_cy   = np.zeros(MAX_NODES, dtype=np.float32)
node_cz   = np.zeros(MAX_NODES, dtype=np.float32)
node_mass = np.zeros(MAX_NODES, dtype=np.float32)
node_size = np.zeros(MAX_NODES, dtype=np.float32)  # side length of node box
node_child    = np.full((MAX_NODES, 8), -1, dtype=np.int32)
node_particle = np.full(MAX_NODES, -1, dtype=np.int32)
node_is_leaf  = np.ones(MAX_NODES, dtype=np.int32)
node_xmin = np.zeros(MAX_NODES, dtype=np.float32)
node_xmax = np.zeros(MAX_NODES, dtype=np.float32)
node_ymin = np.zeros(MAX_NODES, dtype=np.float32)
node_ymax = np.zeros(MAX_NODES, dtype=np.float32)
node_zmin = np.zeros(MAX_NODES, dtype=np.float32)
node_zmax = np.zeros(MAX_NODES, dtype=np.float32)
node_count = [0]


def reset_tree():
    node_cx[:] = 0.0
    node_cy[:] = 0.0
    node_cz[:] = 0.0
    node_mass[:] = 0.0
    node_size[:] = 0.0
    node_child[:] = -1
    node_particle[:] = -1
    node_is_leaf[:] = 1
    node_xmin[:] = 0.0
    node_xmax[:] = 0.0
    node_ymin[:] = 0.0
    node_ymax[:] = 0.0
    node_zmin[:] = 0.0
    node_zmax[:] = 0.0
    node_count[0] = 0


def allocate_node():
    idx = node_count[0]
    node_count[0] += 1
    return idx


def get_octant(node_idx, part_idx):
    mid_x = (node_xmin[node_idx] + node_xmax[node_idx]) * 0.5
    mid_y = (node_ymin[node_idx] + node_ymax[node_idx]) * 0.5
    mid_z = (node_zmin[node_idx] + node_zmax[node_idx]) * 0.5
    oct = 0
    if px_arr[part_idx] >= mid_x:
        oct |= 1
    if py_arr[part_idx] >= mid_y:
        oct |= 2
    if pz_arr[part_idx] >= mid_z:
        oct |= 4
    return oct


def set_child_bbox(parent, child, octant):
    mid_x = (node_xmin[parent] + node_xmax[parent]) * 0.5
    mid_y = (node_ymin[parent] + node_ymax[parent]) * 0.5
    mid_z = (node_zmin[parent] + node_zmax[parent]) * 0.5

    node_xmin[child] = mid_x if (octant & 1) else node_xmin[parent]
    node_xmax[child] = node_xmax[parent] if (octant & 1) else mid_x
    node_ymin[child] = mid_y if (octant & 2) else node_ymin[parent]
    node_ymax[child] = node_ymax[parent] if (octant & 2) else mid_y
    node_zmin[child] = mid_z if (octant & 4) else node_zmin[parent]
    node_zmax[child] = node_zmax[parent] if (octant & 4) else mid_z

    node_size[child] = node_xmax[child] - node_xmin[child]


def insert_particle(root, part_idx):
    # iterative insertion using an explicit stack
    stack = [(root, part_idx)]
    while stack:
        node_idx, pid = stack.pop()
        if node_is_leaf[node_idx]:
            existing = node_particle[node_idx]
            if existing == -1:
                # empty leaf: place particle here
                node_particle[node_idx] = pid
                node_cx[node_idx] = px_arr[pid]
                node_cy[node_idx] = py_arr[pid]
                node_cz[node_idx] = pz_arr[pid]
                node_mass[node_idx] = mass_arr[pid]
            else:
                # occupied leaf: subdivide and re-insert both
                node_is_leaf[node_idx] = 0
                node_particle[node_idx] = -1
                stack.append((node_idx, existing))
                stack.append((node_idx, pid))
        else:
            # internal node: update centre of mass and recurse into child
            m_old = node_mass[node_idx]
            m_new = mass_arr[pid]
            m_total = m_old + m_new
            if m_total > 0:
                node_cx[node_idx] = (node_cx[node_idx] * m_old + px_arr[pid] * m_new) / m_total
                node_cy[node_idx] = (node_cy[node_idx] * m_old + py_arr[pid] * m_new) / m_total
                node_cz[node_idx] = (node_cz[node_idx] * m_old + pz_arr[pid] * m_new) / m_total
                node_mass[node_idx] = m_total

            oct = get_octant(node_idx, pid)
            child_idx = node_child[node_idx, oct]
            if child_idx == -1:
                child_idx = allocate_node()
                set_child_bbox(node_idx, child_idx, oct)
                node_child[node_idx, oct] = child_idx
            stack.append((child_idx, pid))


def build_tree():
    reset_tree()
    root = allocate_node()

    xmin = float(np.min(px_arr)) - 0.0001
    xmax = float(np.max(px_arr)) + 0.0001
    ymin = float(np.min(py_arr)) - 0.0001
    ymax = float(np.max(py_arr)) + 0.0001
    zmin = float(np.min(pz_arr)) - 0.0001
    zmax = float(np.max(pz_arr)) + 0.0001

    node_xmin[root] = xmin
    node_xmax[root] = xmax
    node_ymin[root] = ymin
    node_ymax[root] = ymax
    node_zmin[root] = zmin
    node_zmax[root] = zmax
    node_size[root] = xmax - xmin

    for i in range(N):
        insert_particle(root, i)

    return root


# SCALAR force calculation (reference implementation)

def build_interaction_list_scalar(i, root):
    """Traverse the tree and return a list of node indices for particle i."""
    interaction = []
    stack = [root]
    while stack:
        node_idx = stack.pop()
        if node_idx == -1:
            continue
        if node_is_leaf[node_idx]:
            p = node_particle[node_idx]
            if p != -1 and p != i:
                interaction.append(node_idx)
        else:
            dx = node_cx[node_idx] - px_arr[i]
            dy = node_cy[node_idx] - py_arr[i]
            dz = node_cz[node_idx] - pz_arr[i]
            dist = np.sqrt(dx*dx + dy*dy + dz*dz)
            s = node_size[node_idx]
            if dist > 0 and s / dist < THETA:
                # treat node as a single pseudo-particle
                interaction.append(node_idx)
            else:
                # open the node
                for k in range(8):
                    c = node_child[node_idx, k]
                    if c != -1:
                        stack.append(c)
    return interaction


def compute_force_scalar(i, interaction_list):
    """Scalar force accumulation loop (from Milestone 3 template)."""
    ax = ay = az = 0.0
    for node in interaction_list:
        dx = node_cx[node] - px_arr[i]
        dy = node_cy[node] - py_arr[i]
        dz = node_cz[node] - pz_arr[i]
        dist_sq = dx*dx + dy*dy + dz*dz + SOFTENING
        inv_dist3 = dist_sq ** -1.5
        f = node_mass[node] * inv_dist3
        ax += f * dx
        ay += f * dy
        az += f * dz
    return ax, ay, az


# NUMPY vectorised force calculation

def compute_force_numpy(i, interaction_list):
    """
    NumPy vectorised force kernel.
    The for-loop over interaction_list is eliminated entirely.
    All arithmetic runs on whole arrays at once.
    """
    nodes = np.array(interaction_list, dtype=np.int32)

    # gather node positions and masses using array indexing
    dx = node_cx[nodes] - px_arr[i]
    dy = node_cy[nodes] - py_arr[i]
    dz = node_cz[nodes] - pz_arr[i]

    dist_sq = dx*dx + dy*dy + dz*dz + SOFTENING

    # single vectorised call: replaces the scalar loop's per-element power
    inv_dist3 = dist_sq ** -1.5

    f = node_mass[nodes] * inv_dist3

    ax = np.sum(f * dx)
    ay = np.sum(f * dy)
    az = np.sum(f * dz)
    return float(ax), float(ay), float(az)


def compute_all_forces_scalar(root):
    ax_all = np.zeros(N, dtype=np.float64)
    ay_all = np.zeros(N, dtype=np.float64)
    az_all = np.zeros(N, dtype=np.float64)
    for i in range(N):
        ilist = build_interaction_list_scalar(i, root)
        ax_all[i], ay_all[i], az_all[i] = compute_force_scalar(i, ilist)
    return ax_all, ay_all, az_all


def compute_all_forces_numpy(root):
    ax_all = np.zeros(N, dtype=np.float64)
    ay_all = np.zeros(N, dtype=np.float64)
    az_all = np.zeros(N, dtype=np.float64)
    for i in range(N):
        ilist = build_interaction_list_scalar(i, root)
        ax_all[i], ay_all[i], az_all[i] = compute_force_numpy(i, ilist)
    return ax_all, ay_all, az_all


# Correctness test

def test_correctness(root):
    print("--- Correctness Test ---")
    ax_s, ay_s, az_s = compute_all_forces_scalar(root)
    ax_n, ay_n, az_n = compute_all_forces_numpy(root)

    max_err_x = np.max(np.abs(ax_s - ax_n))
    max_err_y = np.max(np.abs(ay_s - ay_n))
    max_err_z = np.max(np.abs(az_s - az_n))
    max_err = max(max_err_x, max_err_y, max_err_z)

    print(f"  Max absolute error in accelerations: {max_err:.6e}")
    if max_err < 1e-4:
        print("  PASS: NumPy output matches scalar reference.")
    else:
        print("  FAIL: outputs differ beyond tolerance.")
    return ax_n, ay_n, az_n


# Performance benchmark

def benchmark(n_particles, repeats=3):
    global N, px_arr, py_arr, pz_arr, mass_arr

    N = n_particles
    px_arr = np.random.uniform(-1.0, 1.0, N).astype(np.float32)
    py_arr = np.random.uniform(-1.0, 1.0, N).astype(np.float32)
    pz_arr = np.random.uniform(-1.0, 1.0, N).astype(np.float32)
    mass_arr = np.ones(N, dtype=np.float32)

    root = build_tree()

    # build all interaction lists once (same work for both)
    ilists = [build_interaction_list_scalar(i, root) for i in range(N)]

    # time scalar
    t_start = time.perf_counter()
    for _ in range(repeats):
        for i in range(N):
            compute_force_scalar(i, ilists[i])
    t_scalar = (time.perf_counter() - t_start) / repeats

    # time numpy
    t_start = time.perf_counter()
    for _ in range(repeats):
        for i in range(N):
            compute_force_numpy(i, ilists[i])
    t_numpy = (time.perf_counter() - t_start) / repeats

    return t_scalar, t_numpy


# Main


if __name__ == "__main__":
    # reset global N for the correctness test
    N = 100
    px_arr = np.random.uniform(-1.0, 1.0, N).astype(np.float32)
    py_arr = np.random.uniform(-1.0, 1.0, N).astype(np.float32)
    pz_arr = np.random.uniform(-1.0, 1.0, N).astype(np.float32)
    mass_arr = np.ones(N, dtype=np.float32)

    root = build_tree()
    test_correctness(root)

    print("\n--- Performance Benchmark ---")
    print(f"{'N':>6}  {'Scalar (s)':>12}  {'NumPy (s)':>12}  {'Speedup':>8}")
    for n in [100, 500, 1000, 5000]:
        ts, tn = benchmark(n)
        speedup = ts / tn if tn > 0 else float('inf')
        print(f"{n:>6}  {ts:>12.4f}  {tn:>12.4f}  {speedup:>7.1f}x")

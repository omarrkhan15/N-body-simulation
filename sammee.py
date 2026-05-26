import sys
import os
import time
import csv
import numpy as np

# -----------------------------------------------------------------------------
# Path setup
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_dir, "..", "Milestone_3", "barnes_hut_refactor"))
from bh_refactored import (
    reset_pool, allocate_node, insert_particle,
    node_mass, node_child, node_is_leaf, node_particle, node_size,
    node_com_x, node_com_y, node_com_z
)

# -----------------------------------------------------------------------------
# Scalar baseline (identical)
def force_scalar(i, ilist, px, py, pz, mass, eps):
    ax = ay = az = 0.0
    for n in ilist:
        dx = node_com_x[n] - px[i]
        dy = node_com_y[n] - py[i]
        dz = node_com_z[n] - pz[i]
        r2 = dx*dx + dy*dy + dz*dz + eps
        inv_r3 = r2 ** -1.5
        f = node_mass[n] * inv_r3
        ax += f * dx; ay += f * dy; az += f * dz
    return ax, ay, az

# -----------------------------------------------------------------------------
# Vectorised kernel using .sum() method (explicit temporary array)
def force_numpy(i, ilist, px, py, pz, mass, eps,
                cx_arr, cy_arr, cz_arr, m_arr):
    if len(ilist) == 0:
        return 0.0, 0.0, 0.0
    # gather
    nodes = np.asarray(ilist, dtype=np.int32)
    dx = cx_arr[nodes] - px[i]
    dy = cy_arr[nodes] - py[i]
    dz = cz_arr[nodes] - pz[i]
    # arithmetic
    r2 = dx*dx + dy*dy + dz*dz + eps
    inv_r3 = r2 ** -1.5
    f = m_arr[nodes] * inv_r3
    # reduction using .sum() on temporary product array
    ax = (f * dx).sum()
    ay = (f * dy).sum()
    az = (f * dz).sum()
    return ax, ay, az

# -----------------------------------------------------------------------------
# Interaction list builder (explicit stack, natural child order 0..7)
# IMPORTANT: This version uses ITERATION, NOT recursion.
def build_ilist_iterative(p_idx, px, py, pz, theta=0.7):
    ilist = []
    stack = [0]                      # root node
    tx, ty, tz = px[p_idx], py[p_idx], pz[p_idx]

    while stack:
        node = stack.pop()
        if node == -1:
            continue

        dx = node_com_x[node] - tx
        dy = node_com_y[node] - ty
        dz = node_com_z[node] - tz
        d2 = dx*dx + dy*dy + dz*dz
        if d2 == 0.0:
            continue
        d = d2 ** 0.5

        if node_is_leaf[node]:
            p = node_particle[node]
            if p != p_idx and p != -1:
                ilist.append(node)
        else:
            if (node_size[node] / d) < theta:
                ilist.append(node)
            else:
                base = node * 8
                # push children in natural order 0..7
                # Because stack is LIFO, pushing 0..7 results in 7 being processed first.
                # That is fine – the final list order does not matter (addition is commutative).
                for k in range(8):
                    child = node_child[base + k]
                    if child != -1:
                        stack.append(child)
    return ilist

# -----------------------------------------------------------------------------
# Benchmark pipeline (identical to version A, but using the new functions)
def main():
    G = 6.67430e-11
    eps = 1e9
    theta = 0.7
    datasets = {100: "stable_random_system100.csv", 500: "cluster500.csv",
                1000: "cluster1000.csv", 5000: "cluster5000.csv"}
    data_dir = os.path.join(current_dir, "..", "BARNES-HUT", "data")

    print(f"\n{'N':<8} | {'Scalar (s)':<12} | {'NumPy (s)':<12} | {'Speedup':<10} | {'Max Error'}")
    print("-" * 65)

    for N, fname in datasets.items():
        fpath = os.path.join(data_dir, fname)
        if not os.path.exists(fpath):
            print(f"Skipping N={N}: {fname} not found")
            continue

        px, py, pz, mass = [], [], [], []
        with open(fpath, 'r') as f:
            rdr = csv.DictReader(f)
            for row in rdr:
                mass.append(float(row['mass']))
                px.append(float(row['distanceX']))
                py.append(float(row['distanceY']))
                pz.append(float(row['distanceZ']))
        px, py, pz, mass = map(np.array, (px, py, pz, mass))
        n = len(mass)

        # build tree
        reset_pool()
        xmin, xmax = px.min(), px.max()
        ymin, ymax = py.min(), py.max()
        zmin, zmax = pz.min(), pz.max()
        root_cx = (xmin + xmax) * 0.5
        root_cy = (ymin + ymax) * 0.5
        root_cz = (zmin + zmax) * 0.5
        root_sz = max(xmax-xmin, ymax-ymin, zmax-zmin) * 1.1
        root = allocate_node(root_cx, root_cy, root_cz, root_sz)
        for i in range(n):
            insert_particle(root, i, px, py, pz, mass)

        # NumPy views
        np_cx = np.array(node_com_x, dtype=np.float64)
        np_cy = np.array(node_com_y, dtype=np.float64)
        np_cz = np.array(node_com_z, dtype=np.float64)
        np_mass = np.array(node_mass, dtype=np.float64)

        # build interaction lists (iterative, non‑recursive)
        raw = [build_ilist_iterative(i, px, py, pz, theta) for i in range(n)]
        np_lists = [np.array(lst, dtype=np.int32) for lst in raw]

        # scalar timing
        t0 = time.perf_counter()
        sf = []
        for i in range(n):
            fx, fy, fz = force_scalar(i, raw[i], px, py, pz, mass, eps)
            sf.append((fx*G, fy*G, fz*G))
        t_scalar = time.perf_counter() - t0

        # NumPy (.sum) timing
        t0 = time.perf_counter()
        nf = []
        for i in range(n):
            fx, fy, fz = force_numpy(i, np_lists[i], px, py, pz, mass, eps,
                                     np_cx, np_cy, np_cz, np_mass)
            nf.append((fx*G, fy*G, fz*G))
        t_numpy = time.perf_counter() - t0

        speedup = t_scalar / t_numpy
        max_err = max(abs(sf[i][c] - nf[i][c]) for i in range(n) for c in range(3))
        print(f"{n:<8} | {t_scalar:<12.5f} | {t_numpy:<12.5f} | {speedup:<10.2f} | {max_err:.2e}")

if __name__ == "__main__":
    main()
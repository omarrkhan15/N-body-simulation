import math
import copy
import sys

from barnes_hut_refactored import (
    N, MAX_NODES,
    px, py, pz, vx, vy, vz, mass,
    node_cx, node_cy, node_cz, node_mass,
    node_child, node_particle, node_is_leaf,
    node_xmin, node_xmax, node_ymin, node_ymax, node_zmin, node_zmax,
    node_count,
    build_tree, init_random, step, accel_naive, accel_bh,
    DT
)


def collect_all_particles(root):
    found = set()
    stack = [root]
    while stack:
        node = stack.pop()
        if node == -1:
            continue
        if node_is_leaf[node]:
            pidx = node_particle[node]
            if pidx != -1:
                found.add(pidx)
        else:
            for k in range(8):
                child = node_child[node * 8 + k]
                if child != -1:
                    stack.append(child)
    return found


def find_leaf_for_particle(root, part_idx):
    stack = [root]
    while stack:
        node = stack.pop()
        if node == -1:
            continue
        if node_is_leaf[node]:
            if node_particle[node] == part_idx:
                return node
        else:
            for k in range(8):
                child = node_child[node * 8 + k]
                if child != -1:
                    stack.append(child)
    return -1


def test_mass_conservation(root):
    total_mass = sum(mass)
    root_mass  = node_mass[root]
    diff = abs(root_mass - total_mass)
    tol  = 1e-4 * total_mass

    passed = diff <= tol
    status = "PASS" if passed else "FAIL"
    print(f"[{status}] Test 1 — Mass conservation")
    print(f"         root_mass={root_mass:.6f}  sum_mass={total_mass:.6f}  diff={diff:.2e}  tol={tol:.2e}")
    return passed


def test_particle_count(root):
    found = collect_all_particles(root)
    expected = set(range(N))

    missing = expected - found
    extra   = found - expected

    passed = (len(missing) == 0) and (len(extra) == 0)
    status = "PASS" if passed else "FAIL"
    print(f"[{status}] Test 2 — Particle count")
    print(f"         found={len(found)}/{N}  missing={sorted(missing)[:5]}  extra={sorted(extra)[:5]}")
    return passed


def test_centre_of_mass(root):
    total_m = sum(mass)
    com_x = sum(mass[i] * px[i] for i in range(N)) / total_m
    com_y = sum(mass[i] * py[i] for i in range(N)) / total_m
    com_z = sum(mass[i] * pz[i] for i in range(N)) / total_m

    dx = abs(node_cx[root] - com_x)
    dy = abs(node_cy[root] - com_y)
    dz = abs(node_cz[root] - com_z)
    tol = 1e-4

    passed = (dx <= tol) and (dy <= tol) and (dz <= tol)
    status = "PASS" if passed else "FAIL"
    print(f"[{status}] Test 3 — Centre of mass")
    print(f"         computed=({com_x:.5f},{com_y:.5f},{com_z:.5f})")
    print(f"         root_com=({node_cx[root]:.5f},{node_cy[root]:.5f},{node_cz[root]:.5f})")
    print(f"         diff=({dx:.2e},{dy:.2e},{dz:.2e})  tol={tol:.2e}")
    return passed


def test_bounding_boxes(root):
    failures = []
    for i in range(N):
        leaf = find_leaf_for_particle(root, i)
        if leaf == -1:
            failures.append((i, "not found in tree"))
            continue

        ok = (
            node_xmin[leaf] <= px[i] <= node_xmax[leaf] and
            node_ymin[leaf] <= py[i] <= node_ymax[leaf] and
            node_zmin[leaf] <= pz[i] <= node_zmax[leaf]
        )
        if not ok:
            failures.append((i, f"pos=({px[i]:.4f},{py[i]:.4f},{pz[i]:.4f}) "
                                 f"box=([{node_xmin[leaf]:.4f},{node_xmax[leaf]:.4f}],"
                                 f"[{node_ymin[leaf]:.4f},{node_ymax[leaf]:.4f}],"
                                 f"[{node_zmin[leaf]:.4f},{node_zmax[leaf]:.4f}])"))

    passed = len(failures) == 0
    status = "PASS" if passed else "FAIL"
    print(f"[{status}] Test 4 — Bounding boxes  ({N - len(failures)}/{N} particles OK)")
    if failures:
        for pidx, reason in failures[:3]:
            print(f"         particle {pidx}: {reason}")
    return passed


def test_regression():
    import barnes_hut_refactored as bh

    N_test = 20
    TOL    = 0.05

    orig_N = bh.N
    bh.N = N_test

    test_px   = [0.0] * N_test
    test_py   = [0.0] * N_test
    test_pz   = [0.0] * N_test
    test_vx   = [0.0] * N_test
    test_vy   = [0.0] * N_test
    test_vz   = [0.0] * N_test
    test_mass = [1.0] * N_test

    import random
    random.seed(77)
    for i in range(N_test):
        test_px[i]   = random.uniform(-1, 1)
        test_py[i]   = random.uniform(-1, 1)
        test_pz[i]   = random.uniform(-1, 1)
        test_vx[i]   = random.uniform(-0.1, 0.1)
        test_vy[i]   = random.uniform(-0.1, 0.1)
        test_vz[i]   = random.uniform(-0.1, 0.1)
        test_mass[i] = random.uniform(0.5, 1.5)

    naive_px = test_px[:]
    naive_py = test_py[:]
    naive_pz = test_pz[:]
    naive_vx = test_vx[:]
    naive_vy = test_vy[:]
    naive_vz = test_vz[:]

    def naive_accel(lpx, lpy, lpz, lmass, n):
        ax = [0.0]*n; ay = [0.0]*n; az = [0.0]*n
        for i in range(n):
            for j in range(n):
                if i == j: continue
                dx = lpx[j]-lpx[i]; dy = lpy[j]-lpy[i]; dz = lpz[j]-lpz[i]
                d2 = dx*dx+dy*dy+dz*dz + bh.SOFTENING
                inv_d = 1/math.sqrt(d2); inv_d3 = inv_d**3
                f = lmass[j]*inv_d3
                ax[i]+=dx*f; ay[i]+=dy*f; az[i]+=dz*f
        return ax, ay, az

    def naive_step(lpx,lpy,lpz,lvx,lvy,lvz,lmass,n):
        Ax,Ay,Az = naive_accel(lpx,lpy,lpz,lmass,n)
        hdt = 0.5*bh.DT
        for i in range(n):
            lvx[i]+=Ax[i]*hdt; lvy[i]+=Ay[i]*hdt; lvz[i]+=Az[i]*hdt
        for i in range(n):
            lpx[i]+=lvx[i]*bh.DT; lpy[i]+=lvy[i]*bh.DT; lpz[i]+=lvz[i]*bh.DT
        Ax,Ay,Az = naive_accel(lpx,lpy,lpz,lmass,n)
        for i in range(n):
            lvx[i]+=Ax[i]*hdt; lvy[i]+=Ay[i]*hdt; lvz[i]+=Az[i]*hdt

    naive_step(naive_px,naive_py,naive_pz,naive_vx,naive_vy,naive_vz,test_mass,N_test)

    for i in range(N_test):
        bh.px[i]=test_px[i]; bh.py[i]=test_py[i]; bh.pz[i]=test_pz[i]
        bh.vx[i]=test_vx[i]; bh.vy[i]=test_vy[i]; bh.vz[i]=test_vz[i]
        bh.mass[i]=test_mass[i]

    bh.step(use_bh=True)
    bh_px = bh.px[:N_test]
    bh_py = bh.py[:N_test]
    bh_pz = bh.pz[:N_test]

    max_rel_err = 0.0
    for i in range(N_test):
        ref_d = math.sqrt(naive_px[i]**2 + naive_py[i]**2 + naive_pz[i]**2) + 1e-10
        err   = math.sqrt((bh_px[i]-naive_px[i])**2 +
                          (bh_py[i]-naive_py[i])**2 +
                          (bh_pz[i]-naive_pz[i])**2)
        rel   = err / ref_d
        if rel > max_rel_err:
            max_rel_err = rel

    bh.N = orig_N

    passed = max_rel_err <= TOL
    status = "PASS" if passed else "FAIL"
    print(f"[{status}] Test 5 — Regression (N={N_test}, 1 timestep)")
    print(f"         max relative position error = {max_rel_err:.4f}  tol={TOL}")
    return passed


def run_all_tests():
    print("=" * 60)
    print("Barnes-Hut Refactored Implementation — Test Suite")
    print("=" * 60)

    init_random(seed=42)
    root = build_tree()
    print(f"Tree built: {node_count[0]} nodes for {N} particles\n")

    results = []
    results.append(test_mass_conservation(root));   print()
    results.append(test_particle_count(root));       print()
    results.append(test_centre_of_mass(root));       print()
    results.append(test_bounding_boxes(root));       print()
    results.append(test_regression());               print()

    passed = sum(results)
    total  = len(results)
    print("=" * 60)
    print(f"Results: {passed}/{total} tests passed")
    print("=" * 60)

    if passed == total:
        print("ALL TESTS PASSED ")
        return 0
    else:
        print("SOME TESTS FAILED ")
        return 1


if __name__ == "__main__":
    sys.exit(run_all_tests())

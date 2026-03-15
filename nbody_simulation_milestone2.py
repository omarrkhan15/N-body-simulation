import math
import csv
from nbody_visualizer import draw_gui

G = 6.67430e-11
dt = 8640
softening = 1e9
N = 100

# using fixed size arrays now, no append
# this is so it maps better to assembly
m   = [0.0] * N
p_x = [0.0] * N
p_y = [0.0] * N
p_z = [0.0] * N
v_x = [0.0] * N
v_y = [0.0] * N
v_z = [0.0] * N
a_x = [0.0] * N
a_y = [0.0] * N
a_z = [0.0] * N

i = 0
with open('stable_random_system100.csv', mode='r') as file:
    reader = csv.DictReader(file)
    for row in reader:
        m[i]   = float(row["mass"])
        p_x[i] = float(row["distanceX"])
        p_y[i] = float(row["distanceY"])
        p_z[i] = float(row["distanceZ"])
        v_x[i] = float(row["velocityX"])
        v_y[i] = float(row["velocityY"])
        v_z[i] = float(row["velocityZ"])
        i += 1


def calculate_acceleration():
    for i in range(N):
        a_x[i] = 0.0
        a_y[i] = 0.0
        a_z[i] = 0.0

    for i in range(N):
        for j in range(N):
            if i != j:
                dx = p_x[j] - p_x[i]
                dy = p_y[j] - p_y[i]
                dz = p_z[j] - p_z[i]

                r_sq = dx*dx + dy*dy + dz*dz
                r_soft_sq = r_sq + softening*softening
                r_soft = math.sqrt(r_soft_sq)
                r_soft_cubed = r_soft * r_soft_sq

                f = G * m[j] / r_soft_cubed

                a_x[i] += f * dx
                a_y[i] += f * dy
                a_z[i] += f * dz


def kick():
    for i in range(N):
        v_x[i] += 0.5 * a_x[i] * dt
        v_y[i] += 0.5 * a_y[i] * dt
        v_z[i] += 0.5 * a_z[i] * dt


def drift():
    for i in range(N):
        p_x[i] += v_x[i] * dt
        p_y[i] += v_y[i] * dt
        p_z[i] += v_z[i] * dt


calculate_acceleration()

step = 0
while draw_gui(p_x, p_y, p_z):
    kick()
    drift()
    calculate_acceleration()
    kick()
    step += 1
    if step % 100 == 0:
        print(f"step {step}")

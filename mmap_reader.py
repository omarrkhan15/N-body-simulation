import mmap
import struct
import time
import os
from nbody_visualizer import draw_gui

N = 100
ARR_BYTES = N * 4
# layout: [flag(4 bytes)][px][py][pz]
TOTAL_SIZE = 4 + 3 * ARR_BYTES
SHARED_FILE = "shared.mem"

# make the file if it doesnt exist yet
if not os.path.exists(SHARED_FILE):
    with open(SHARED_FILE, 'wb') as f:
        f.write(b'\x00' * TOTAL_SIZE)


def read_positions(mm):
    # spin until assembly says its done writing (flag goes to 0)
    while True:
        mm.seek(0)
        flag = struct.unpack('I', mm.read(4))[0]
        if flag == 0:
            break
        time.sleep(0.0001)

    mm.seek(4)
    px = list(struct.unpack(f'{N}f', mm.read(ARR_BYTES)))

    mm.seek(4 + ARR_BYTES)
    py = list(struct.unpack(f'{N}f', mm.read(ARR_BYTES)))

    mm.seek(4 + 2 * ARR_BYTES)
    pz = list(struct.unpack(f'{N}f', mm.read(ARR_BYTES)))

    return px, py, pz


with open(SHARED_FILE, 'r+b') as f:
    mm = mmap.mmap(f.fileno(), TOTAL_SIZE)

    px = [0.0] * N
    py = [0.0] * N
    pz = [0.0] * N

    while draw_gui(px, py, pz):
        try:
            px, py, pz = read_positions(mm)
        except Exception as e:
            print(f"error reading shared mem: {e}")

    mm.close()

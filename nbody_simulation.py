import math
import csv
from nbody_visualizer import draw_gui

# Simulation parameters
G = 6.67430e-11  
dt = 8640
softening = 1e9
METHOD = "naive"  # Change to "barnes-hut" to switch algorithms
THETA = 0.5  # Barnes-Hut accuracy parameter (lower = more accurate, slower)

# Particle data
m = []      # masses
v_x = []    # x velocities
v_y = []    # y velocities
v_z = []    # z velocities
p_x = []    # x positions
p_y = []    # y positions
p_z = []    # z positions
a_x = []    # x accelerations
a_y = []    # y accelerations
a_z = []    # z accelerations

# Load initial conditions
with open('stable_random_system100.csv', mode='r') as file:
    reader = csv.DictReader(file)
    for row in reader:
        m.append(float(row["mass"]))
        p_x.append(float(row["distanceX"]))
        p_y.append(float(row["distanceY"]))
        p_z.append(float(row["distanceZ"]))
        v_x.append(float(row["velocityX"]))
        v_y.append(float(row["velocityY"]))
        v_z.append(float(row["velocityZ"]))
        a_x.append(0.0)
        a_y.append(0.0)
        a_z.append(0.0)

N = len(m)

#  NAIVE ALGORITHM 
def calculate_acceleration_naive():
    
    # Reset all accelerations to zero
    for i in range(N):
        a_x[i] = 0.0
        a_y[i] = 0.0
        a_z[i] = 0.0
    
    # Calculate pairwise forces
    for i in range(N):
        for j in range(N):
            if i != j:  
                # Calculate distance vector from body i to body j
                dx = p_x[j] - p_x[i]
                dy = p_y[j] - p_y[i]
                dz = p_z[j] - p_z[i]
                
                # Calculate distance with softening to avoid singularities
                r_squared = dx*dx + dy*dy + dz*dz
                r_soft_squared = r_squared + softening*softening
                r_soft = math.sqrt(r_soft_squared)
                r_soft_cubed = r_soft * r_soft_squared
                
                # Calculate force magnitude divided by mass of body 
                force_per_mass = G * m[j] / r_soft_cubed
                
                # Add acceleration contribution in each direction
                a_x[i] += force_per_mass * dx
                a_y[i] += force_per_mass * dy
                a_z[i] += force_per_mass * dz


# BARNES-HUT ALGORITHM 
class OctreeNode:

    def __init__(self, x_min, x_max, y_min, y_max, z_min, z_max):
        
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.z_min = z_min
        self.z_max = z_max
        
        # Center of mass coordinates and total mass
        self.com_x = 0.0
        self.com_y = 0.0
        self.com_z = 0.0
        self.total_mass = 0.0
        
        # Eight children (octants)
        self.children = [None] * 8
        
        # If this is a leaf node, stores the body index
        self.body_index = None
        self.is_leaf = True
        self.is_empty = True


def get_octant(node, x, y, z):

    mid_x = (node.x_min + node.x_max) / 2
    mid_y = (node.y_min + node.y_max) / 2
    mid_z = (node.z_min + node.z_max) / 2
    
    octant = 0
    if x >= mid_x: octant += 1
    if y >= mid_y: octant += 2
    if z >= mid_z: octant += 4
    
    return octant


def get_octant_bounds(node, octant):

    mid_x = (node.x_min + node.x_max) / 2
    mid_y = (node.y_min + node.y_max) / 2
    mid_z = (node.z_min + node.z_max) / 2
    
    # Define bounds for each of the 8 octants
    bounds = [
        [node.x_min, mid_x, node.y_min, mid_y, node.z_min, mid_z],  # 0: (-, -, -)
        [mid_x, node.x_max, node.y_min, mid_y, node.z_min, mid_z],  # 1: (+, -, -)
        [node.x_min, mid_x, mid_y, node.y_max, node.z_min, mid_z],  # 2: (-, +, -)
        [mid_x, node.x_max, mid_y, node.y_max, node.z_min, mid_z],  # 3: (+, +, -)
        [node.x_min, mid_x, node.y_min, mid_y, mid_z, node.z_max],  # 4: (-, -, +)
        [mid_x, node.x_max, node.y_min, mid_y, mid_z, node.z_max],  # 5: (+, -, +)
        [node.x_min, mid_x, mid_y, node.y_max, mid_z, node.z_max],  # 6: (-, +, +)
        [mid_x, node.x_max, mid_y, node.y_max, mid_z, node.z_max],  # 7: (+, +, +)
    ]
    
    return bounds[octant]


def insert_body(node, body_idx):
  
    x, y, z = p_x[body_idx], p_y[body_idx], p_z[body_idx]
    mass = m[body_idx]
    
    # Case 1: Node is empty - just place the body here
    if node.is_empty:
        node.body_index = body_idx
        node.is_empty = False
        node.is_leaf = True
        node.total_mass = mass
        node.com_x = x
        node.com_y = y
        node.com_z = z
        return
    
    # Case 2: Node contains one body - need to subdivide
    if node.is_leaf:
        # Save the existing body
        old_body_idx = node.body_index
        old_x, old_y, old_z = p_x[old_body_idx], p_y[old_body_idx], p_z[old_body_idx]
        
        # Convert to internal node
        node.is_leaf = False
        node.body_index = None
        
        # Re-insert the old body into appropriate child
        octant_old = get_octant(node, old_x, old_y, old_z)
        bounds = get_octant_bounds(node, octant_old)
        node.children[octant_old] = OctreeNode(*bounds)
        insert_body(node.children[octant_old], old_body_idx)
    
    # Case 3: Node is internal - insert into appropriate child
    octant = get_octant(node, x, y, z)
    if node.children[octant] is None:
        bounds = get_octant_bounds(node, octant)
        node.children[octant] = OctreeNode(*bounds)
    
    insert_body(node.children[octant], body_idx)
    
    # Update center of mass using weighted average
    total_mass_new = node.total_mass + mass
    node.com_x = (node.com_x * node.total_mass + x * mass) / total_mass_new
    node.com_y = (node.com_y * node.total_mass + y * mass) / total_mass_new
    node.com_z = (node.com_z * node.total_mass + z * mass) / total_mass_new
    node.total_mass = total_mass_new


def build_octree():
    # Find bounding box containing all particles
    x_min = min(p_x)
    x_max = max(p_x)
    y_min = min(p_y)
    y_max = max(p_y)
    z_min = min(p_z)
    z_max = max(p_z)
    
    # Add padding to avoid particles exactly on boundaries
    padding = max(x_max - x_min, y_max - y_min, z_max - z_min) * 0.1
    x_min -= padding
    x_max += padding
    y_min -= padding
    y_max += padding
    z_min -= padding
    z_max += padding
    
    # Create root node spanning entire space
    root = OctreeNode(x_min, x_max, y_min, y_max, z_min, z_max)
    
    # Insert all bodies into the tree
    for i in range(N):
        insert_body(root, i)
    
    return root


def calculate_force_from_node(body_idx, node, ax, ay, az):

    # Skip empty nodes
    if node is None or node.is_empty:
        return ax, ay, az
    
    # Don't calculate force from a body on itself
    if node.is_leaf and node.body_index == body_idx:
        return ax, ay, az
    
    # Calculate distance to node's center of mass
    dx = node.com_x - p_x[body_idx]
    dy = node.com_y - p_y[body_idx]
    dz = node.com_z - p_z[body_idx]
    
    distance = math.sqrt(dx*dx + dy*dy + dz*dz)
    
    # Calculate size of the node (maximum dimension)
    size = max(node.x_max - node.x_min, 
               node.y_max - node.y_min,
               node.z_max - node.z_min)
    
    # If node is far enough or is a leaf, treat as single body
    if node.is_leaf or (distance > 0 and size / distance < THETA):
        # Apply softening parameter
        r_soft_squared = distance*distance + softening*softening
        r_soft = math.sqrt(r_soft_squared)
        r_soft_cubed = r_soft * r_soft_squared
        
        # Calculate force magnitude per unit mass
        force_per_mass = G * node.total_mass / r_soft_cubed
        
        # Add acceleration contribution
        ax += force_per_mass * dx
        ay += force_per_mass * dy
        az += force_per_mass * dz
    else:
        # Node is too close - recurse into children for better accuracy
        for child in node.children:
            if child is not None:
                ax, ay, az = calculate_force_from_node(body_idx, child, ax, ay, az)
    
    return ax, ay, az


def calculate_acceleration_barnes_hut():

    # Reset all accelerations
    for i in range(N):
        a_x[i] = 0.0
        a_y[i] = 0.0
        a_z[i] = 0.0
    
    # Build octree from current positions
    tree = build_octree()
    
    # Calculate force on each body using the tree
    for i in range(N):
        a_x[i], a_y[i], a_z[i] = calculate_force_from_node(i, tree, 0.0, 0.0, 0.0)


# INTEGRATION FUNCTIONS
def calculate_acceleration():
    
    if METHOD == "naive":
        calculate_acceleration_naive()
    elif METHOD == "barnes-hut":
        calculate_acceleration_barnes_hut()
    else:
        raise ValueError(f"Unknown method: {METHOD}. Use 'naive' or 'barnes-hut'")


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


# MAIN SIMULATION LOOP
print(f"=== N-Body Simulation ===")
print(f"Number of bodies: {N}")
print(f"Method: {METHOD}")
print(f"Timestep: {dt} seconds")
print(f"Softening parameter: {softening}")
if METHOD == "barnes-hut":
    print(f"Barnes-Hut theta: {THETA}")
print(f"Starting simulation...")

step = 0
while draw_gui(p_x, p_y, p_z):
    # (Kick-Drift-Kick)
    kick()                      
    drift()                     
    calculate_acceleration() # Calculate forces at new positions
    kick()                     
    
    step += 1
    if step % 100 == 0:
        print(f"Step {step} complete")

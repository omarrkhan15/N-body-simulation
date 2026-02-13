# N-body-simulation
# N-Body Simulation Project

## What This Does
Simulates how stars and planets move in space due to gravity.

## Team name
mubarak paan shop

## Team Members
- muhammad omar khan 31599
- sammee mudassar 30765
- sajjad ahmed 30406

## How It Works
We built two different algorithms:

### 1. Naive Algorithm
- Checks every object against every other object
- Simple but slow
- Good for small simulations (< 500 bodies)

### 2. Barnes-Hut Algorithm  
- Uses a tree to group far-away objects
- Much faster for large simulations
- Can handle thousands of bodies

## Running the Code

```bash
python nbody_simulation.py
```

Change the `METHOD` variable at the top to switch between algorithms:
- `METHOD = "naive"` - Simple direct calculation
- `METHOD = "barnes-hut"` - Fast tree-based method

## Files
- `nbody_simulation.py` - Main simulation code
- `nbody_visualizer.py` - Graphics/GUI (provided by TA)
- Use the CSV files provided by the TA for initial conditions

## What We Learned
- Leapfrog integration is way better than Euler for keeping energy stable
- Barnes-Hut is ~10x faster than naive for 500 bodies
- Softening parameter prevents crashes when objects get too close
- Smart algorithms make a huge difference in performance

## Requirements
- Python 3.x
- numpy
- (whatever graphics library the visualizer uses)

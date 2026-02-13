import pygame
import math
import random

# Global state for the GUI
_gui_state = None

class _GUIState:
    """Internal state management for the GUI"""
    def __init__(self, grid_size, trail_length, num_particles):
        self.grid_size = grid_size
        self.trail_length = trail_length
        self.num_particles = num_particles
        
        # Visual properties
        self.radii = [5]
        self.radii += [random.uniform(0.5, 2.5) for _ in range(num_particles - 1)]
        
        self.colors = [(255, 255, 0)]
        self.colors += [(random.randint(100, 255), random.randint(100, 255), 
                        random.randint(100, 255)) for _ in range(num_particles - 1)]
        
        # Camera controls
        self.angle_y = 0
        self.angle_x = 0
        self.view_scale = 1.0
        
        # Trail management
        self.trails = [[] for _ in range(num_particles)]
        
        # Pygame setup
        pygame.init()
        self.screen = pygame.display.set_mode((grid_size, grid_size))
        pygame.display.set_caption("N-Body Simulation")
        self.running = True
        
    def update_camera(self):
        """Handle keyboard input for camera rotation"""
        keys = pygame.key.get_pressed()
        if keys[pygame.K_LEFT]:
            self.angle_y += 0.02
        if keys[pygame.K_RIGHT]:
            self.angle_y -= 0.02
        if keys[pygame.K_UP]:
            self.angle_x += 0.02
        if keys[pygame.K_DOWN]:
            self.angle_x -= 0.02
    
    def update_view_scale(self, positions):
        """Calculate appropriate view scale based on particle positions"""
        max_val = 0
        for particle in positions:
            max_val = max(max_val, abs(particle[0]), abs(particle[1]), abs(particle[2]))
        self.view_scale = max(max_val * 0.5, 1e-10)  # Avoid division by zero
    
    def project_3d_to_2d(self, x, y, z):
        """Project 3D coordinates to 2D screen coordinates"""
        # Rotation on y axis
        curr_x = x * math.cos(self.angle_y) + z * math.sin(self.angle_y)
        curr_z = -x * math.sin(self.angle_y) + z * math.cos(self.angle_y)
        
        # Rotation on x axis
        curr_y = y * math.cos(self.angle_x) - curr_z * math.sin(self.angle_x)
        final_z = y * math.sin(self.angle_x) + curr_z * math.cos(self.angle_x)
        
        # Perspective projection
        z_offset = (final_z / self.view_scale) + 2
        factor = 1 / z_offset * 2
        pos_x = (curr_x * factor / self.view_scale + 1) / 2 * self.grid_size
        pos_y = (curr_y * factor / self.view_scale + 1) / 2 * self.grid_size
        
        return pos_x, pos_y, final_z, factor


def draw_gui(p_x, p_y, p_z, grid_size=800, trail_length=100):
    """
    Draw particles on the GUI. Call this function once per frame in your main loop.
    
    Args:
        p_x: List of x-coordinates for all particles
        p_y: List of y-coordinates for all particles
        p_z: List of z-coordinates for all particles
        grid_size: Size of the window (default: 800)
        trail_length: Number of past positions to show as trails (default: 100)
    
    Returns:
        True if the window is still open, False if user closed it
    
    Example usage:
        while draw_gui(p_x, p_y, p_z):
            # Your physics code here
            # Update p_x, p_y, p_z
    """
    global _gui_state
    
    # Initialize on first call
    if _gui_state is None:
        num_particles = len(p_x)
        _gui_state = _GUIState(grid_size, trail_length, num_particles)
    
    # Handle events
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            _gui_state.running = False
            pygame.quit()
            return False
    
    if not _gui_state.running:
        return False
    
    # Update camera
    _gui_state.update_camera()
    
    # Combine coordinates into positions
    positions = [[p_x[i], p_y[i], p_z[i]] for i in range(len(p_x))]
    
    # Update view scale
    _gui_state.update_view_scale(positions)
    
    # Update trails
    for i in range(len(positions)):
        _gui_state.trails[i].append(tuple(positions[i]))
        if len(_gui_state.trails[i]) > _gui_state.trail_length:
            _gui_state.trails[i].pop(0)
    
    # Clear screen
    _gui_state.screen.fill((0, 0, 0))
    
    # Create draw list with depth sorting
    draw_list = []
    for i, particle in enumerate(positions):
        x, y, z = particle[0], particle[1], particle[2]
        pos_x, pos_y, final_z, factor = _gui_state.project_3d_to_2d(x, y, z)
        draw_list.append((final_z, pos_x, pos_y, factor, i))
    
    # Sort by depth (furthest first)
    draw_list.sort(key=lambda x: x[0], reverse=True)
    
    # Draw trails and particles
    for z_depth, x, y, f, i in draw_list:
        # Draw trail
        if len(_gui_state.trails[i]) > 1:
            projected_trail = []
            for pt in _gui_state.trails[i]:
                tx, ty, _, _ = _gui_state.project_3d_to_2d(pt[0], pt[1], pt[2])
                projected_trail.append((tx, ty))
            
            pygame.draw.lines(_gui_state.screen, (70, 70, 70), False, projected_trail, 1)
        
        # Draw particle
        r_size = max(1, int(_gui_state.radii[i] * f * 2))
        pygame.draw.circle(_gui_state.screen, _gui_state.colors[i], (int(x), int(y)), r_size)
    
    # Update display
    pygame.display.flip()
    
    return True

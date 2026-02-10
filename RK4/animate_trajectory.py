"""
Animation script for RK4 projectile simulation visualization.
Displays the trajectory, velocity direction (speed versor), angular velocity (spin vector),
and body axis. Highlights frames where the magnetic field is active.

Features:
- Cylindrical target visualization
- Fixed scenario with configurable starting point
- Collision detection to determine hit/miss
- Visual feedback on whether magnetic field successfully deflected the bullet

Usage:
    python animate_trajectory.py [trajectory.csv] [--save animation.mp4] [--fps 30] [--skip N]
    python animate_trajectory.py [trajectory.csv] --scenario  # Run fixed scenario with target

Dependencies:
    pip install numpy pandas matplotlib
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection
import matplotlib.animation as animation
import argparse
import sys


# ============================================================================
# SCENARIO CONFIGURATION
# ============================================================================
# Fixed scenario: Bullet fired at a cylindrical target, magnetic field attempts deflection

SCENARIO_CONFIG = {
    # Target cylinder parameters (in meters)
    'target_center': np.array([5.0, -0.75, 0.0]),  # Center of cylinder base (5m away, centered on trajectory)
    'target_radius': 0.5,                           # Cylinder radius (50 cm diameter = 1m wide)
    'target_height': 1.5,                           # Cylinder height (1.5 m)
    'target_axis': np.array([0.0, 1.0, 0.0]),      # Cylinder axis direction (vertical)
    
    # Bullet starting point (should match simulation initial conditions)
    'bullet_start': np.array([0.0, 0.0, 0.0]),
    
    # Visualization settings
    'cylinder_color': 'orange',
    'cylinder_alpha': 0.4,
    'hit_color': 'red',
    'miss_color': 'green',
}


def create_cylinder_mesh(center, radius, height, axis, n_segments=30):
    """
    Create mesh vertices for a 3D cylinder.
    
    Parameters:
    -----------
    center : array-like
        Center of the cylinder base
    radius : float
        Cylinder radius
    height : float
        Cylinder height
    axis : array-like
        Unit vector along cylinder axis
    n_segments : int
        Number of segments for circular cross-section
        
    Returns:
    --------
    vertices : list of arrays
        List of polygon vertices for cylinder surface
    """
    center = np.array(center)
    axis = np.array(axis)
    axis = axis / np.linalg.norm(axis)
    
    # Find two perpendicular vectors to the axis
    if abs(axis[0]) < 0.9:
        perp1 = np.cross(axis, np.array([1, 0, 0]))
    else:
        perp1 = np.cross(axis, np.array([0, 1, 0]))
    perp1 = perp1 / np.linalg.norm(perp1)
    perp2 = np.cross(axis, perp1)
    
    # Generate circle points
    theta = np.linspace(0, 2 * np.pi, n_segments + 1)
    circle_points = np.array([
        radius * (np.cos(t) * perp1 + np.sin(t) * perp2) for t in theta
    ])
    
    # Bottom and top circles
    bottom = center + circle_points
    top = center + axis * height + circle_points
    
    vertices = []
    
    # Side faces (quads)
    for i in range(n_segments):
        quad = [bottom[i], bottom[i+1], top[i+1], top[i]]
        vertices.append(quad)
    
    # Bottom cap (as triangle fan)
    for i in range(n_segments):
        tri = [center, bottom[i], bottom[i+1]]
        vertices.append(tri)
    
    # Top cap
    top_center = center + axis * height
    for i in range(n_segments):
        tri = [top_center, top[i+1], top[i]]
        vertices.append(tri)
    
    return vertices


def point_to_cylinder_distance(point, center, radius, height, axis):
    """
    Compute the shortest distance from a point to a cylinder.
    
    Returns:
    --------
    distance : float
        Signed distance (negative if inside cylinder)
    is_inside : bool
        True if point is inside the cylinder
    """
    center = np.array(center)
    axis = np.array(axis)
    axis = axis / np.linalg.norm(axis)
    point = np.array(point)
    
    # Vector from center to point
    v = point - center
    
    # Project onto axis
    along_axis = np.dot(v, axis)
    
    # Radial component
    radial_vec = v - along_axis * axis
    radial_dist = np.linalg.norm(radial_vec)
    
    # Check if inside
    inside_height = 0 <= along_axis <= height
    inside_radius = radial_dist <= radius
    
    is_inside = inside_height and inside_radius
    
    if is_inside:
        # Distance to nearest surface
        dist_to_side = radius - radial_dist
        dist_to_bottom = along_axis
        dist_to_top = height - along_axis
        distance = -min(dist_to_side, dist_to_bottom, dist_to_top)
    else:
        # Distance from outside
        if inside_height:
            distance = radial_dist - radius
        elif inside_radius:
            if along_axis < 0:
                distance = -along_axis
            else:
                distance = along_axis - height
        else:
            # Corner case
            if along_axis < 0:
                corner = center
            else:
                corner = center + axis * height
            to_corner = point - corner
            along = np.dot(to_corner, axis) * axis
            radial = to_corner - along
            radial_norm = np.linalg.norm(radial)
            if radial_norm > 1e-10:
                radial_dir = radial / radial_norm
            else:
                # Point is on axis, pick arbitrary perpendicular direction
                if abs(axis[0]) < 0.9:
                    radial_dir = np.cross(axis, np.array([1, 0, 0]))
                else:
                    radial_dir = np.cross(axis, np.array([0, 1, 0]))
                radial_dir = radial_dir / np.linalg.norm(radial_dir)
            nearest_on_rim = corner + radial_dir * radius
            distance = np.linalg.norm(point - nearest_on_rim)
    
    return distance, is_inside


def check_trajectory_hits_cylinder(positions, center, radius, height, axis):
    """
    Check if any point in the trajectory hits the cylinder.
    
    Parameters:
    -----------
    positions : array (N, 3)
        Trajectory positions
    center, radius, height, axis : cylinder parameters
    
    Returns:
    --------
    hit : bool
        True if trajectory intersects cylinder
    hit_index : int or None
        Index of first hit point (or None if no hit)
    min_distance : float
        Minimum distance from trajectory to cylinder surface
    """
    min_distance = float('inf')
    hit_index = None
    
    for i, pos in enumerate(positions):
        dist, is_inside = point_to_cylinder_distance(pos, center, radius, height, axis)
        if dist < min_distance:
            min_distance = dist
        if is_inside and hit_index is None:
            hit_index = i
    
    hit = hit_index is not None
    return hit, hit_index, min_distance
# Optional progress bar for long operations (animation save, plotting)
try:
    from tqdm import tqdm
except Exception:
    tqdm = None


def load_trajectory(filename):
    """Load trajectory data from CSV file."""
    df = pd.read_csv(filename)
    return df


def determine_scale(values, threshold_ratio=1e-3):
    """
    Determine appropriate scale and unit for a set of values.
    
    Parameters:
    -----------
    values : array-like
        The values to analyze
    threshold_ratio : float
        If range/max_abs < threshold_ratio, consider values as negligible
    
    Returns:
    --------
    scale : float
        Multiplier to convert to the chosen unit
    unit : str
        Unit label (m, mm, μm, nm)
    is_negligible : bool
        True if the variation is negligible
    """
    values = np.asarray(values)
    val_min, val_max = np.min(values), np.max(values)
    val_range = val_max - val_min
    max_abs = max(abs(val_min), abs(val_max), 1e-15)
    
    # Check if negligible variation
    is_negligible = val_range / max_abs < threshold_ratio if max_abs > 1e-15 else True
    
    # Determine unit based on the magnitude of the values
    if max_abs >= 1:
        return 1.0, 'm', is_negligible
    elif max_abs >= 1e-3:
        return 1e3, 'mm', is_negligible
    elif max_abs >= 1e-6:
        return 1e6, 'μm', is_negligible
    else:
        return 1e9, 'nm', is_negligible


def determine_common_scale(x, y, z):
    """
    Determine a common scale for 3D trajectory data.
    Uses the maximum range across all axes to pick appropriate units.
    
    Returns:
    --------
    scale : float
    unit : str
    negligible_axes : list of str ('x', 'y', 'z' that have negligible variation)
    """
    ranges = {
        'x': (np.min(x), np.max(x)),
        'y': (np.min(y), np.max(y)),
        'z': (np.min(z), np.max(z))
    }
    
    # Find max absolute value across all data
    all_vals = np.concatenate([x, y, z])
    max_abs = np.max(np.abs(all_vals))
    
    # Determine scale based on max absolute value
    if max_abs >= 1:
        scale, unit = 1.0, 'm'
    elif max_abs >= 1e-3:
        scale, unit = 1e3, 'mm'
    elif max_abs >= 1e-6:
        scale, unit = 1e6, 'μm'
    else:
        scale, unit = 1e9, 'nm'
    
    # Identify negligible axes (variation << max range)
    max_range = max(ranges['x'][1] - ranges['x'][0],
                   ranges['y'][1] - ranges['y'][0],
                   ranges['z'][1] - ranges['z'][0])
    max_range = max(max_range, 1e-15)
    
    negligible_axes = []
    threshold = 1e-3  # 0.1% of max range
    for axis, (vmin, vmax) in ranges.items():
        if (vmax - vmin) / max_range < threshold:
            negligible_axes.append(axis)
    
    return scale, unit, negligible_axes


def create_animation(df, skip=1, fps=30, save_path=None, scenario_mode=False):
    """
    Create 3D animation of the projectile motion.
    
    Parameters:
    -----------
    df : DataFrame
        Trajectory data with columns: t, x, y, z, vhat_x, vhat_y, vhat_z, 
        omega_x, omega_y, omega_z, ex, ey, ez, impulse
    skip : int
        Frame skip factor for faster animation
    fps : int
        Frames per second for saved animation
    save_path : str or None
        Path to save animation (MP4 format), or None for interactive display
    scenario_mode : bool
        If True, show cylinder target and hit/miss detection
    """
    
    # Subsample data for animation
    df_anim = df.iloc[::skip].reset_index(drop=True)
    n_frames = len(df_anim)
    
    # Extract data arrays
    t = df_anim['t'].values
    pos_raw = df_anim[['x', 'y', 'z']].values
    vhat = df_anim[['vhat_x', 'vhat_y', 'vhat_z']].values
    omega = df_anim[['omega_x', 'omega_y', 'omega_z']].values
    e_axis = df_anim[['ex', 'ey', 'ez']].values
    impulse = df_anim['impulse'].values
    
    # Determine common scale for all position axes
    scale, unit, negligible_axes = determine_common_scale(
        pos_raw[:, 0], pos_raw[:, 1], pos_raw[:, 2]
    )
    
    # Scale position data
    pos = pos_raw * scale
    
    # Normalize omega for visualization (scale to reasonable arrow length)
    omega_mag = np.linalg.norm(omega, axis=1, keepdims=True)
    omega_mag_max = np.max(omega_mag) if np.max(omega_mag) > 0 else 1.0
    omega_normalized = omega / omega_mag_max

    # If scenario_mode, swap Y and Z for plotting so cylinder axis appears vertical
    if scenario_mode:
        # Re-map coordinates: (x, y, z) -> (x, z, y)
        pos = pos[:, [0, 2, 1]]
        vhat = vhat[:, [0, 2, 1]]
        omega_normalized = omega_normalized[:, [0, 2, 1]]
        e_axis = e_axis[:, [0, 2, 1]]
    
    # Compute trajectory bounds for axis limits
    margin_factor = 0.1
    x_min, x_max = pos[:, 0].min(), pos[:, 0].max()
    y_min, y_max = pos[:, 1].min(), pos[:, 1].max()
    z_min, z_max = pos[:, 2].min(), pos[:, 2].max()
    
    # Add margins
    x_range = max(x_max - x_min, 1e-10)
    y_range = max(y_max - y_min, 1e-10)
    z_range = max(z_max - z_min, 1e-10)
    
    x_min -= x_range * margin_factor
    x_max += x_range * margin_factor
    y_min -= y_range * margin_factor
    y_max += y_range * margin_factor
    z_min -= z_range * margin_factor
    z_max += z_range * margin_factor
    
    # Make axes roughly equal scale for proper 3D visualization
    max_range = max(x_max - x_min, y_max - y_min, z_max - z_min)
    x_center = (x_max + x_min) / 2
    y_center = (y_max + y_min) / 2
    z_center = (z_max + z_min) / 2
    
    # Vector arrow scale (relative to trajectory size)
    arrow_scale = max_range * 0.15
    
    # Create figure and 3D axis
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Initialize plot elements
    # Trajectory line (will be updated)
    traj_line, = ax.plot([], [], [], 'b-', linewidth=1, alpha=0.6, label='Trajectory')
    
    # Current position marker
    pos_marker, = ax.plot([], [], [], 'ko', markersize=8, label='Position')
    
    # Velocity direction arrow (green)
    vhat_quiver = ax.quiver(0, 0, 0, 0, 0, 0, color='green', linewidth=2, 
                            arrow_length_ratio=0.2, label='Velocity direction')
    
    # Angular velocity arrow (red)
    omega_quiver = ax.quiver(0, 0, 0, 0, 0, 0, color='red', linewidth=2,
                             arrow_length_ratio=0.2, label='Angular velocity')
    
    # Body axis arrow (blue)
    e_quiver = ax.quiver(0, 0, 0, 0, 0, 0, color='blue', linewidth=2,
                         arrow_length_ratio=0.2, label='Body axis')
    
    # Impulse indicator (background color changes)
    impulse_text = ax.text2D(0.02, 0.98, '', transform=ax.transAxes, fontsize=12,
                             verticalalignment='top', fontweight='bold')
    
    # Time display
    time_text = ax.text2D(0.02, 0.92, '', transform=ax.transAxes, fontsize=10,
                          verticalalignment='top')
    
    # Negligible axes info
    if negligible_axes:
        negl_text = f"Note: {', '.join(negligible_axes).upper()} axes have negligible variation"
        ax.text2D(0.02, 0.02, negl_text, transform=ax.transAxes, fontsize=9,
                  color='gray', verticalalignment='bottom')
    
    # Set axis limits
    ax.set_xlim(x_center - max_range/2, x_center + max_range/2)
    ax.set_ylim(y_center - max_range/2, y_center + max_range/2)
    ax.set_zlim(z_center - max_range/2, z_center + max_range/2)
    
    # Set axis labels with consistent units
    ax.set_xlabel(f'X ({unit})')
    ax.set_ylabel(f'Y ({unit})')
    ax.set_zlabel(f'Z ({unit})')
    
    # Scenario mode: add cylinder target and check for collision
    hit_detected = False
    hit_index = None
    min_distance = None
    
    if scenario_mode:
        # Get scenario configuration
        cfg = SCENARIO_CONFIG
        target_center = cfg['target_center'] * scale
        target_radius = cfg['target_radius'] * scale
        target_height = cfg['target_height'] * scale
        target_axis = cfg['target_axis']
        
        # Check if trajectory hits the cylinder (use raw positions for accuracy)
        hit_detected, hit_index, min_distance = check_trajectory_hits_cylinder(
            pos_raw, cfg['target_center'], cfg['target_radius'], 
            cfg['target_height'], cfg['target_axis']
        )
        
        # Create and add cylinder mesh (mapped to plotting coordinates if swapped)
        if scenario_mode:
            # Map target center and axis to plotting coords: (x,y,z)->(x,z,y)
            plot_center = np.array([target_center[0], target_center[2], target_center[1]])
            plot_axis = np.array([target_axis[0], target_axis[2], target_axis[1]])
            cylinder_verts = create_cylinder_mesh(
                plot_center, target_radius, target_height, plot_axis
            )
        else:
            cylinder_verts = create_cylinder_mesh(
                target_center, target_radius, target_height, target_axis
            )
        
        # Determine cylinder color based on hit/miss
        if hit_detected:
            cyl_color = cfg['hit_color']
            result_text = 'TARGET HIT!'
            result_color = 'red'
        else:
            cyl_color = cfg['miss_color']
            result_text = f'TARGET MISSED (min dist: {min_distance*1000:.1f} mm)'
            result_color = 'green'
        
        # Add cylinder to plot
        cylinder_collection = Poly3DCollection(
            cylinder_verts, 
            alpha=cfg['cylinder_alpha'],
            facecolor=cyl_color,
            edgecolor='darkgray',
            linewidth=0.5
        )
        ax.add_collection3d(cylinder_collection)
        
        # In scenario mode: Focus view on the region around the target
        # Set scene bounds around plotted target (use plot_center when coordinates are swapped)
        scene_x_min = -1.0 * scale  # A bit before start
        scene_x_max = (cfg['target_center'][0] + cfg['target_radius'] + 2.0) * scale  # Past target

        # If plotting coords are swapped, plot_center was computed above
        if scenario_mode:
            pc = plot_center
            # pc[1] is plotted Y (original Z), pc[2] is plotted Z (original Y)
            scene_y_min = (pc[1] - target_radius - 1.5 * scale)
            scene_y_max = (pc[1] + target_radius + 1.5 * scale)
            scene_z_min = (pc[2] - target_height/2 - 1.5 * scale)
            scene_z_max = (pc[2] + target_height/2 + 1.5 * scale)
        else:
            # Set Y and Z range to show the cylinder with margin (original coords)
            cyl_half_height = cfg['target_height'] / 2
            cyl_center_y = cfg['target_center'][1] + cyl_half_height
            scene_y_min = (cyl_center_y - cyl_half_height - 1.5) * scale
            scene_y_max = (cyl_center_y + cyl_half_height + 1.5) * scale

            scene_z_min = (cfg['target_center'][2] - cfg['target_radius'] - 1.5) * scale
            scene_z_max = (cfg['target_center'][2] + cfg['target_radius'] + 1.5) * scale
        
        # Make the view roughly cubic for proper proportions
        x_span = scene_x_max - scene_x_min
        y_span = scene_y_max - scene_y_min
        z_span = scene_z_max - scene_z_min
        max_span = max(x_span, y_span, z_span)
        
        # Center each axis and expand to max_span
        x_center = (scene_x_max + scene_x_min) / 2
        y_center = (scene_y_max + scene_y_min) / 2
        z_center = (scene_z_max + scene_z_min) / 2
        
        ax.set_xlim(x_center - max_span/2 * 1.1, x_center + max_span/2 * 1.1)
        ax.set_ylim(y_center - max_span/2 * 1.1, y_center + max_span/2 * 1.1)
        ax.set_zlim(z_center - max_span/2 * 1.1, z_center + max_span/2 * 1.1)
        
        # Update arrow scale for the zoomed view
        arrow_scale = max_span * 0.08
        
        # Update max_range for later use
        max_range = max_span
        
        # Add result text
        result_display = ax.text2D(0.5, 0.95, result_text, transform=ax.transAxes, 
                                   fontsize=16, fontweight='bold', color=result_color,
                                   ha='center', va='top',
                                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Add scenario info
        scenario_info = f"Target: {cfg['target_radius']*2:.1f}m diameter cylinder at x={cfg['target_center'][0]:.0f}m"
        ax.text2D(0.5, 0.88, scenario_info, transform=ax.transAxes, fontsize=10, 
                  ha='center', va='top', color='gray')
        
        ax.set_title('Projectile vs Magnetic Deflection Scenario')
    else:
        ax.set_title('Projectile Motion with Magnetic Impulse')
    
    # Legend
    ax.legend(loc='upper right')
    
    def init():
        """Initialize animation."""
        traj_line.set_data([], [])
        traj_line.set_3d_properties([])
        pos_marker.set_data([], [])
        pos_marker.set_3d_properties([])
        impulse_text.set_text('')
        time_text.set_text('')
        return traj_line, pos_marker, impulse_text, time_text
    
    def update(frame):
        """Update animation frame."""
        nonlocal vhat_quiver, omega_quiver, e_quiver
        
        # Update trajectory line (show path up to current frame)
        traj_line.set_data(pos[:frame+1, 0], pos[:frame+1, 1])
        traj_line.set_3d_properties(pos[:frame+1, 2])
        
        # Update position marker - change color at hit point in scenario mode
        marker_color = 'ko'
        if scenario_mode and hit_detected and hit_index is not None:
            # Scale hit_index for subsampled data
            hit_frame_approx = hit_index // skip
            if frame >= hit_frame_approx:
                marker_color = 'ro'  # Red marker after hit
        
        pos_marker.set_data([pos[frame, 0]], [pos[frame, 1]])
        pos_marker.set_3d_properties([pos[frame, 2]])
        pos_marker.set_color('red' if 'r' in marker_color else 'black')
        
        # Remove old quivers
        vhat_quiver.remove()
        omega_quiver.remove()
        e_quiver.remove()
        
        # Current position
        px, py, pz = pos[frame]
        
        # Velocity direction arrow
        vhat_quiver = ax.quiver(px, py, pz, 
                                vhat[frame, 0] * arrow_scale,
                                vhat[frame, 1] * arrow_scale,
                                vhat[frame, 2] * arrow_scale,
                                color='green', linewidth=2, arrow_length_ratio=0.2)
        
        # Angular velocity arrow (scaled)
        omega_quiver = ax.quiver(px, py, pz,
                                 omega_normalized[frame, 0] * arrow_scale,
                                 omega_normalized[frame, 1] * arrow_scale,
                                 omega_normalized[frame, 2] * arrow_scale,
                                 color='red', linewidth=2, arrow_length_ratio=0.2)
        
        # Body axis arrow
        e_quiver = ax.quiver(px, py, pz,
                             e_axis[frame, 0] * arrow_scale,
                             e_axis[frame, 1] * arrow_scale,
                             e_axis[frame, 2] * arrow_scale,
                             color='blue', linewidth=2, arrow_length_ratio=0.2)
        
        # Update impulse indicator
        if impulse[frame]:
            impulse_text.set_text('MAGNETIC FIELD ACTIVE')
            impulse_text.set_color('red')
            ax.set_facecolor('#ffe6e6')  # Light red background
        else:
            impulse_text.set_text('Magnetic field inactive')
            impulse_text.set_color('gray')
            ax.set_facecolor('white')
        
        # Update time display
        time_text.set_text(f't = {t[frame]:.3f} s')
        
        return traj_line, pos_marker, impulse_text, time_text
    
    # Create animation
    anim = animation.FuncAnimation(fig, update, frames=n_frames,
                                   init_func=init, blit=False,
                                   interval=1000/fps, repeat=True)
    
    # Print scenario results
    if scenario_mode:
        print("\n" + "="*50)
        print("SCENARIO RESULTS")
        print("="*50)
        if hit_detected:
            hit_time = t[hit_index // skip] if hit_index // skip < len(t) else t[-1]
            print(f"  Result: TARGET HIT!")
            print(f"  Hit at approximately t = {hit_time:.4f} s")
        else:
            print(f"  Result: TARGET MISSED")
            print(f"  Minimum distance to target: {min_distance*1000:.2f} mm")
        print(f"  Magnetic field: {'deflected bullet' if not hit_detected else 'insufficient deflection'}")
        print("="*50 + "\n")
    
    if save_path:
        print(f"Saving animation to {save_path}...")
        # Try to use ffmpeg writer, fall back to pillow for gif
        if save_path.endswith('.gif'):
            writer = animation.PillowWriter(fps=fps)
        else:
            writer = animation.FFMpegWriter(fps=fps, bitrate=2000)

        # Use tqdm progress bar if available and matplotlib supports progress_callback
        if tqdm is not None:
            try:
                with tqdm(total=n_frames, desc='Saving animation', unit='frames') as pbar:
                    def _progress_callback(i, n):
                        # i is the current frame index (1-based in some mpl versions)
                        # Ensure we don't over-update the bar
                        to_update = int(i) - pbar.n
                        if to_update > 0:
                            pbar.update(to_update)

                    anim.save(save_path, writer=writer, progress_callback=_progress_callback)
            except TypeError:
                # Older matplotlib may not support progress_callback
                print('matplotlib progress_callback not supported; saving without tqdm')
                anim.save(save_path, writer=writer)
        else:
            # No tqdm installed
            anim.save(save_path, writer=writer)

        print(f"Animation saved to {save_path}")
    else:
        plt.show()
    
    return anim


def plot_static_analysis(df):
    """
    Create static plots for trajectory analysis.
    Uses appropriate unit scaling for all axes.
    """
    fig = plt.figure(figsize=(14, 10))
    
    t = df['t'].values
    pos_raw = df[['x', 'y', 'z']].values
    omega = df[['omega_x', 'omega_y', 'omega_z']].values
    impulse = df['impulse'].values
    
    # Determine position scale
    pos_scale, pos_unit, negligible_pos = determine_common_scale(
        pos_raw[:, 0], pos_raw[:, 1], pos_raw[:, 2]
    )
    pos = pos_raw * pos_scale
    
    # Determine omega scale (rad/s is usually fine, but check magnitude)
    omega_mag_max = np.max(np.abs(omega))
    if omega_mag_max >= 1e3:
        omega_scale, omega_unit = 1e-3, 'krad/s'
    elif omega_mag_max >= 1:
        omega_scale, omega_unit = 1, 'rad/s'
    else:
        omega_scale, omega_unit = 1e3, 'mrad/s'
    omega_scaled = omega * omega_scale
    
    # Find impulse region for highlighting
    impulse_mask = impulse == 1
    t_impulse_start = t[impulse_mask][0] if impulse_mask.any() else None
    t_impulse_end = t[impulse_mask][-1] if impulse_mask.any() else None
    
    # Plot 1: 3D trajectory
    ax1 = fig.add_subplot(2, 2, 1, projection='3d')
    
    # Color trajectory by impulse state
    iter_range = range(len(t) - 1)
    if tqdm is not None:
        iter_range = tqdm(iter_range, desc='Plotting trajectory segments')
    for i in iter_range:
        color = 'red' if impulse[i] else 'blue'
        ax1.plot(pos[i:i+2, 0], pos[i:i+2, 1], pos[i:i+2, 2], color=color, linewidth=1)
    
    ax1.set_xlabel(f'X ({pos_unit})')
    ax1.set_ylabel(f'Y ({pos_unit})')
    ax1.set_zlabel(f'Z ({pos_unit})')
    ax1.set_title('3D Trajectory (blue=normal, red=impulse)')
    
    # Add note about negligible axes
    if negligible_pos:
        ax1.text2D(0.02, 0.02, f"Negligible: {', '.join(negligible_pos).upper()}", 
                   transform=ax1.transAxes, fontsize=8, color='gray')
    
    # Make 3D plot axes equal scale
    max_range = max(pos[:, 0].max() - pos[:, 0].min(),
                   pos[:, 1].max() - pos[:, 1].min(),
                   pos[:, 2].max() - pos[:, 2].min())
    if max_range > 0:
        mid_x = (pos[:, 0].max() + pos[:, 0].min()) / 2
        mid_y = (pos[:, 1].max() + pos[:, 1].min()) / 2
        mid_z = (pos[:, 2].max() + pos[:, 2].min()) / 2
        ax1.set_xlim(mid_x - max_range/2, mid_x + max_range/2)
        ax1.set_ylim(mid_y - max_range/2, mid_y + max_range/2)
        ax1.set_zlim(mid_z - max_range/2, mid_z + max_range/2)
    
    # Plot 2: Position vs time
    ax2 = fig.add_subplot(2, 2, 2)
    
    # Determine which position components to plot (skip negligible ones or scale them appropriately)
    labels = ['x', 'y', 'z']
    colors = ['tab:blue', 'tab:orange', 'tab:green']
    
    for i, (lbl, col) in enumerate(zip(labels, colors)):
        val_range = pos[:, i].max() - pos[:, i].min()
        max_val = max(abs(pos[:, i].max()), abs(pos[:, i].min()))
        
        # Only plot if there's meaningful variation
        if val_range > 1e-10 * max_val or max_val < 1e-15:
            ax2.plot(t, pos[:, i], label=lbl, color=col)
        else:
            # Negligible variation - plot with note
            ax2.plot(t, pos[:, i], label=f'{lbl} (const)', color=col, alpha=0.3, linestyle='--')
    
    if t_impulse_start:
        ax2.axvspan(t_impulse_start, t_impulse_end, alpha=0.3, color='red', label='Impulse')
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel(f'Position ({pos_unit})')
    ax2.set_title('Position vs Time')
    ax2.legend()
    ax2.grid(True)
    
    # Plot 3: Angular velocity magnitude vs time
    ax3 = fig.add_subplot(2, 2, 3)
    omega_mag = np.linalg.norm(omega_scaled, axis=1)
    ax3.plot(t, omega_mag, 'r-')
    if t_impulse_start:
        ax3.axvspan(t_impulse_start, t_impulse_end, alpha=0.3, color='red', label='Impulse')
    ax3.set_xlabel('Time (s)')
    ax3.set_ylabel(f'|ω| ({omega_unit})')
    ax3.set_title('Angular Velocity Magnitude vs Time')
    ax3.grid(True)
    
    # Plot 4: Angular velocity components
    ax4 = fig.add_subplot(2, 2, 4)
    omega_labels = ['ωx', 'ωy', 'ωz']
    omega_colors = ['tab:blue', 'tab:orange', 'tab:green']
    
    for i, (lbl, col) in enumerate(zip(omega_labels, omega_colors)):
        val_range = omega_scaled[:, i].max() - omega_scaled[:, i].min()
        max_val = max(abs(omega_scaled[:, i].max()), abs(omega_scaled[:, i].min()))
        
        if val_range > 1e-10 * max_val or max_val < 1e-15:
            ax4.plot(t, omega_scaled[:, i], label=lbl, color=col)
        else:
            ax4.plot(t, omega_scaled[:, i], label=f'{lbl} (const)', color=col, alpha=0.3, linestyle='--')
    
    if t_impulse_start:
        ax4.axvspan(t_impulse_start, t_impulse_end, alpha=0.3, color='red', label='Impulse')
    ax4.set_xlabel('Time (s)')
    ax4.set_ylabel(f'ω ({omega_unit})')
    ax4.set_title('Angular Velocity Components vs Time')
    ax4.legend()
    ax4.grid(True)
    
    plt.tight_layout()
    plt.savefig('trajectory_analysis.png', dpi=150)
    print("Static analysis saved to trajectory_analysis.png")
    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description='Animate RK4 projectile simulation trajectory'
    )
    parser.add_argument('trajectory', nargs='?', default='output.csv',
                        help='Path to trajectory CSV file (default: output.csv). '
                             'Must match the --out path used with sim.exe')
    parser.add_argument('--save', type=str, default=None,
                        help='Save animation to file (e.g., animation.mp4 or animation.gif)')
    parser.add_argument('--fps', type=int, default=30,
                        help='Frames per second (default: 30)')
    parser.add_argument('--skip', type=int, default=10,
                        help='Frame skip factor for faster playback (default: 10)')
    parser.add_argument('--static', action='store_true',
                        help='Generate static analysis plots instead of animation')
    parser.add_argument('--scenario', action='store_true',
                        help='Run fixed scenario mode with target cylinder to test magnetic deflection')
    
    args = parser.parse_args()
    
    # Load trajectory data
    print(f"Loading trajectory from {args.trajectory}...")
    try:
        df = load_trajectory(args.trajectory)
    except FileNotFoundError:
        print(f"Error: Could not find file {args.trajectory}")
        print("Run the simulation first with: ./sim --out traj.csv")
        sys.exit(1)
    
    print(f"Loaded {len(df)} data points")
    print(f"Time range: {df['t'].min():.3f} - {df['t'].max():.3f} s")
    
    # Check for impulse frames
    impulse_frames = df['impulse'].sum()
    print(f"Frames with magnetic field active: {impulse_frames}")
    
    if args.scenario:
        print("\n--- SCENARIO MODE ---")
        print(f"Target cylinder at x={SCENARIO_CONFIG['target_center'][0]}m")
        print(f"Cylinder radius: {SCENARIO_CONFIG['target_radius']}m, height: {SCENARIO_CONFIG['target_height']}m")
        print("Testing if magnetic field can deflect bullet from target...\n")
    
    if args.static:
        plot_static_analysis(df)
    else:
        create_animation(df, skip=args.skip, fps=args.fps, save_path=args.save, 
                         scenario_mode=args.scenario)


if __name__ == '__main__':
    main()
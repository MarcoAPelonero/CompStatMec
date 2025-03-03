import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

# Parameters for the big cube and grid subdivision
L = 10       # side length of the big cube
M = 2        # side length of each small cube
n = int(L / M)  # number of small cubes along each axis

# Function to compute the 8 vertices of a small cube given its lower corner and side length
def cube_vertices(origin, M):
    x, y, z = origin
    return np.array([
        [x,     y,     z],
        [x+M,   y,     z],
        [x+M,   y+M,   z],
        [x,     y+M,   z],
        [x,     y,     z+M],
        [x+M,   y,     z+M],
        [x+M,   y+M,   z+M],
        [x,     y+M,   z+M],
    ])

# Function to get the faces (each face is a list of vertices) of a cube
def cube_faces(origin, M):
    v = cube_vertices(origin, M)
    faces = [
        [v[0], v[1], v[2], v[3]],  # bottom
        [v[4], v[5], v[6], v[7]],  # top
        [v[0], v[1], v[5], v[4]],  # front
        [v[2], v[3], v[7], v[6]],  # back
        [v[1], v[2], v[6], v[5]],  # right
        [v[0], v[3], v[7], v[4]]   # left
    ]
    return faces

# Function to draw a cube: either just its edges or filled with a transparent color.
def draw_cube(ax, origin, M, color='k', lw=0.5, fill=False, facecolor=None, alpha=1.0):
    faces = cube_faces(origin, M)
    if fill:
        poly3d = Poly3DCollection(faces, facecolors=facecolor, edgecolors=color, linewidths=lw, alpha=alpha)
        ax.add_collection3d(poly3d)
    else:
        # Plot each face's outline
        for face in faces:
            xs = [vertex[0] for vertex in face] + [face[0][0]]
            ys = [vertex[1] for vertex in face] + [face[0][1]]
            zs = [vertex[2] for vertex in face] + [face[0][2]]
            ax.plot(xs, ys, zs, color=color, linewidth=lw)

# Create the figure and 3D axis.
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection='3d')

# Draw the entire grid of small cubes in light gray.
for i in range(n):
    for j in range(n):
        for k in range(n):
            origin = (i * M, j * M, k * M)
            draw_cube(ax, origin, M, color='gray', lw=0.5, fill=False)

# Define three points inside the big cube. 
# (We choose points that are not too near the boundary so that the full 3×3×3 neighborhood exists.)
points = [
    {'pt': (3, 3, 3), 'color': 'red'},    # Lies in cube with grid index (1,1,1)
    {'pt': (7, 7, 3), 'color': 'green'},    # Lies in cube with grid index (3,3,1)
    {'pt': (7, 3, 7), 'color': 'blue'}      # Lies in cube with grid index (3,1,3)
]

# For each point, highlight the 27 (3×3×3) adjacent cubes and plot the point.
for item in points:
    pt = item['pt']
    color = item['color']
    # Plot the point
    ax.scatter(pt[0], pt[1], pt[2], color=color, s=50)
    
    # Determine the grid cell (small cube) in which the point resides.
    grid_idx = (int(pt[0] // M), int(pt[1] // M), int(pt[2] // M))
    ix, iy, iz = grid_idx
    
    # Determine indices for the 3×3×3 neighborhood (taking care of boundary conditions).
    ix_range = range(max(0, ix-1), min(n, ix+2))
    iy_range = range(max(0, iy-1), min(n, iy+2))
    iz_range = range(max(0, iz-1), min(n, iz+2))
    
    # Highlight each small cube in the neighborhood
    for i_n in ix_range:
        for j_n in iy_range:
            for k_n in iz_range:
                origin = (i_n * M, j_n * M, k_n * M)
                draw_cube(ax, origin, M, color=color, lw=2, fill=True, facecolor=color, alpha=0.3)

# Set axis limits and labels.
ax.set_xlim(0, L)
ax.set_ylim(0, L)
ax.set_zlim(0, L)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Cube subdivided into small cubes with highlighted 3x3x3 neighborhoods')

plt.show()

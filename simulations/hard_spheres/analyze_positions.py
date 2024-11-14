import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

OUTPUT_INTERVAL = 1000  # Should match the value in your C++ code
BOX_SIZE = 10.0         # Should match BOX_SIZE in your C++ code

def read_positions(filename):
    positions = []
    with open(filename, 'r') as f:
        timestep = []
        for line in f:
            stripped_line = line.strip()
            if stripped_line == '':
                if timestep:
                    positions.append(np.array(timestep))
                    timestep = []
            else:
                coords = list(map(float, stripped_line.split()))
                timestep.append(coords)
        if timestep:
            positions.append(np.array(timestep))
    return positions

def animate_evolution(positions, interval=OUTPUT_INTERVAL):
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(0, BOX_SIZE)
    ax.set_ylim(0, BOX_SIZE)
    ax.set_zlim(0, BOX_SIZE)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Particle Evolution')

    scatters = [ax.scatter([], [], [], s=20, c='blue', alpha=0.6) for _ in range(len(positions[0]))]

    def update(frame):
        for scatter, pos in zip(scatters, positions[frame]):
            scatter._offsets3d = (pos[:, 0], pos[:, 1], pos[:, 2])
        return scatters

    ani = FuncAnimation(fig, update, frames=len(positions), interval=interval, blit=False)
    plt.show()

def main():
    positions = read_positions('positions.dat')
    animate_evolution(positions)

if __name__ == '__main__':
    main()

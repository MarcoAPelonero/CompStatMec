import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

def run_simulation(output_file, density, num_particles, dt, num_steps):
    """
    Runs the simulation with the given parameters.
    """
    try:
        subprocess.run(["./build/main.exe", output_file, str(density), str(num_particles), str(dt), str(num_steps)], check=True)
    except Exception as e:
        print(f"Error running simulation: {e}")
        raise

def read_simulation_data(filename):
    """
    Reads the simulation data from the output file.
    """
    try:
        with open(filename, 'r') as f:
            lines = f.read().strip().split('\n\n')
        header = lines[0].strip().split()
        data_blocks = [block.strip().split('\n') for block in lines[1:]]
        raw_data = []
        for block in data_blocks:
            snapshot_data = [list(map(float, line.split())) for line in block]
            raw_data.append(snapshot_data)
        data_reshaped = np.array(raw_data)
        return header, data_reshaped
    except Exception as e:
        print(f"Error reading data: {e}")
        raise

def animate_simulation(data, boxSize):
    fig = plt.figure(figsize=(12, 12))
    ax3d = fig.add_subplot(221, projection='3d')
    ax3d.set_title("Particle Positions")
    ax3d.set_xlim([0, boxSize])
    ax3d.set_ylim([0, boxSize])
    ax3d.set_zlim([0, boxSize])

    ax_energy = fig.add_subplot(222)
    ax_energy.set_title("Energy Evolution")
    ax_energy.set_xlabel("Snapshot")
    ax_energy.set_ylabel("Energy")

    ax_velocity = fig.add_subplot(223)
    ax_velocity.set_title("Velocity Evolution")
    ax_velocity.set_xlabel("Snapshot")
    ax_velocity.set_ylabel("Velocity")

    particles, = ax3d.plot([], [], [], 'bo', markersize=2)
    energy_line, = ax_energy.plot([], [], 'r-')
    velocity_line, = ax_velocity.plot([], [], 'b-')

    energies = []
    velocities = []

    def update(frame):
        positions = data[frame, :, 0:3]
        particles.set_data(positions[:, 0], positions[:, 1])
        particles.set_3d_properties(positions[:, 2])

        energy = data[frame, :, 6].mean()
        energies.append(energy)
        energy_line.set_data(range(len(energies)), energies)
        ax_energy.set_xlim([0, len(energies)])
        ax_energy.set_ylim([min(energies), max(energies)])

        velocity = np.linalg.norm(data[frame, :, 3:6], axis=1).mean()
        velocities.append(velocity)
        velocity_line.set_data(range(len(velocities)), velocities)
        ax_velocity.set_xlim([0, len(velocities)])
        ax_velocity.set_ylim([min(velocities), max(velocities)])

        return particles, energy_line, velocity_line

    ani = animation.FuncAnimation(fig, update, frames=data.shape[0], interval=1, blit=True)
    plt.show()

if __name__ == "__main__":
    output_file = "output.txt"
    density = 0.5
    num_particles = 100
    dt = 0.001
    num_steps = 10000
    boxSize = np.cbrt(num_particles/density)
    run_simulation(output_file, density, num_particles, dt, num_steps)
    header, data = read_simulation_data(output_file)
    animate_simulation(data, boxSize)
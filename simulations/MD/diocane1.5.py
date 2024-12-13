import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt

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

if __name__ == "__main__":
    output_file = "output.txt"
    density = 0.5
    num_particles = 100
    dt = 0.001
    num_steps = 10000

    run_simulation(output_file, density, num_particles, dt, num_steps)
    header, data = read_simulation_data(output_file)

    # Calculate average position
    avg_positions = np.mean(data[:, :, 0:3], axis=1)  # Shape: (num_snapshots, 3)

    # Calculate average speed
    speeds = np.linalg.norm(data[:, :, 3:6], axis=2)  # Shape: (num_snapshots, num_particles)
    avg_speeds = np.mean(speeds, axis=1)  # Shape: (num_snapshots,)

    # Calculate average energy
    energies = data[:, :, 6]  # Assuming energy is in column 6
    avg_energies = np.mean(energies, axis=1)  # Shape: (num_snapshots,)

    # Plotting
    fig = plt.figure(figsize=(15, 5))

    # Average Position (3D)
    ax1 = fig.add_subplot(131, projection='3d')
    ax1.set_title('Average Position (3D)')
    ax1.plot(avg_positions[:, 0], avg_positions[:, 1], avg_positions[:, 2], 'o-')
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z')

    # Average Speed
    ax2 = fig.add_subplot(132)
    ax2.set_title('Average Speed')
    ax2.plot(avg_speeds, 'r-')
    ax2.set_xlabel('Snapshot')
    ax2.set_ylabel('Average Speed')

    # Average Energy
    ax3 = fig.add_subplot(133)
    ax3.set_title('Average Energy')
    ax3.plot(avg_energies, 'b-')
    ax3.set_xlabel('Snapshot')
    ax3.set_ylabel('Average Energy')

    plt.tight_layout()
    plt.show()
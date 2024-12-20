import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from tqdm import tqdm

def read_ensemble_data(filename):
    """
    Reads the ensemble data file generated by the simulation.
    The file contains multiple snapshots separated by blank lines.
    The first line is a header: x y z vx vy vz energy
    """
    print(f"Opening file: {filename}")
    try:
        with open(filename, 'r') as f:
            lines = f.read().strip().split('\n\n')  # Split by blank lines

        print("File read successfully.")
        header = lines[0].strip().split()  # First line is the header

        data_blocks = [block.strip().split('\n') for block in lines[1:]]
        print(f"Number of snapshots detected: {len(data_blocks)}")

        numParticles = len(data_blocks[0])  # Number of lines in the first snapshot
        print(f"Number of particles per snapshot: {numParticles}")

        raw_data = []
        for block in tqdm(data_blocks, desc="Processing snapshots"):
            snapshot_data = [list(map(float, line.split())) for line in block]
            raw_data.append(snapshot_data)

        data_reshaped = np.array(raw_data)  # Convert to numpy array
        print(f"Data shape: {data_reshaped.shape}")
        return header, data_reshaped

    except Exception as e:
        print(f"Error reading data: {e}")
        raise

def animate_particles(data):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    boxSize = 6.93361
    ax.set_xlim([0, boxSize])
    ax.set_ylim([0, boxSize])
    ax.set_zlim([0, boxSize])

    particles, = ax.plot([], [], [], 'bo', markersize=2)

    def update(frame):
        positions = data[frame, :, 0:3]
        particles.set_data(positions[:, 0], positions[:, 1])
        particles.set_3d_properties(positions[:, 2])
        return particles,

    ani = animation.FuncAnimation(fig, update, frames=data.shape[0], interval=50, blit=True)
    plt.show()

def plot_energies(data):
    try:
        energies = data[:, :, 6].mean(axis=1)
        plt.plot(energies, label="Average Energy")
        plt.xlabel("Snapshot")
        plt.ylabel("Energy")
        plt.title("Energy Evolution")
        plt.legend()
        plt.show()
    except Exception as e:
        print(f"Error plotting energies: {e}")

def plot_velocities(data):
    try:
        velocities = np.linalg.norm(data[:, :, 3:6], axis=2).mean(axis=1)
        plt.plot(velocities, label="Average Velocity")
        plt.xlabel("Snapshot")
        plt.ylabel("Velocity")
        plt.title("Velocity Evolution")
        plt.legend()
        plt.show()
    except Exception as e:
        print(f"Error plotting velocities: {e}")

def live_plot(data):
    try:
        fig, axes = plt.subplots(4, 4, figsize=(12, 12))
        fig.suptitle("Simulation Live Metrics")

        ax3d = fig.add_subplot(4, 4, 1, projection='3d')
        ax3d.set_xlim([0, 100])
        ax3d.set_ylim([0, 100])
        ax3d.set_zlim([0, 100])

        particles, = ax3d.plot([], [], [], 'bo', markersize=2)
        avg_energy_ax = axes[1, 0]
        avg_vel_ax = axes[2, 0]
        avg_pos_ax = axes[3, 0]

        avg_energy_data = []
        avg_velocity_data = []
        avg_positions_data = []

        def update(frame):
            positions = data[frame, :, 0:3]
            particles.set_data(positions[:, 0], positions[:, 1])
            particles.set_3d_properties(positions[:, 2])

            avg_energy_data.append(data[frame, :, 6].mean())
            avg_velocity_data.append(np.linalg.norm(data[frame, :, 3:6], axis=1).mean())
            avg_positions_data.append(np.linalg.norm(positions, axis=1).mean())

            avg_energy_ax.plot(avg_energy_data, color='red')
            avg_vel_ax.plot(avg_velocity_data, color='blue')
            avg_pos_ax.plot(avg_positions_data, color='green')

            return particles,

        ani = animation.FuncAnimation(fig, update, frames=data.shape[0], interval=20, blit=False)
        plt.show()
    except Exception as e:
        print(f"Error in live plotting: {e}")

if __name__ == "__main__":
    try:
        header, data = read_ensemble_data("ensemble_data.txt")
        print("Header:", header)
        print("Data shape:", data.shape)

        while True:
            try:
                print("\nOptions:")
                print("1. Animate particles")
                print("2. Plot energy evolution")
                print("3. Plot velocity evolution")
                print("4. Live plot with metrics")
                print("5. Exit")

                choice = input("Enter your choice: ")
                if choice == "1":
                    animate_particles(data)
                elif choice == "2":
                    plot_energies(data)
                elif choice == "3":
                    plot_velocities(data)
                elif choice == "4":
                    live_plot(data)
                elif choice == "5":
                    break
                else:
                    print("Invalid choice. Please try again.")
            except Exception as e:
                print(f"Error during choice execution: {e}")
    except Exception as e:
        print(f"Fatal error: {e}")

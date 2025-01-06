import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from tqdm import tqdm

# Load data function
def load_data(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    particle_data = []
    step_data = []
    L_values = []
    N_particles = 0
    skip_next = False

    for line in tqdm(lines, desc='Loading data', total=len(lines)):
        stripped = line.strip()
        if skip_next:
            skip_next = False
            continue
        if stripped:
            # Assuming L is a single float value
            try:
                L = float(stripped)
                L_values.append(L)
                skip_next = True  # Skip the next line (snapshot data)
            except ValueError:
                numbers = list(map(float, stripped.split()))
                step_data.append(numbers)
        else:
            if step_data:
                particle_data.append(step_data)
                if not N_particles:
                    N_particles = len(step_data)
                step_data = []

    if step_data:
        particle_data.append(step_data)
        if not N_particles:
            N_particles = len(step_data)

    particle_data = np.array(particle_data)
    return particle_data, N_particles, L_values

# Visualization mode
def visualize_3d_animation(data, N_particles, L_values):
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    scatter = ax.scatter([], [], [], c='b', s=20)

    def update(frame):
        L = L_values[frame]
        ax.set_xlim(0, L)
        ax.set_ylim(0, L)
        ax.set_zlim(0, L)
        ax.set_title(f"Time step: {frame}")
        positions = data[frame, :, :3]
        scatter._offsets3d = (positions[:, 0], positions[:, 1], positions[:, 2])

    ani = FuncAnimation(fig, update, frames=len(data), interval=20)
    plt.show()

# Subplots for energy evolution
def visualize_energy(data, N_particles):
    steps = len(data)

    avg_potential = []
    avg_kinetic = []
    avg_mechanical = []

    for step in range(steps):
        potential = data[step, :, 6]
        kinetic = data[step, :, 7]
        mechanical = data[step, :, 8]

        avg_potential.append(np.mean(potential))
        avg_kinetic.append(np.mean(kinetic))
        avg_mechanical.append(np.mean(mechanical))

    # Calculate average potential energy after 10% of steps
    avg_potential_after_10pct = np.mean(avg_potential[int(0.2 * steps):])

    time = np.arange(steps)

    fig, axs = plt.subplots(2, 2, figsize=(10, 8))

    axs[0, 0].plot(time, avg_potential, label='Avg Potential Energy')
    axs[0, 0].set_title(f'Avg Potential Energy (After 20%: {avg_potential_after_10pct:.2f})')

    axs[0, 1].plot(time, avg_kinetic, label='Avg Kinetic Energy')
    axs[0, 1].set_title('Avg Kinetic Energy')

    axs[1, 0].plot(time, avg_mechanical, label='Avg Mechanical Energy')
    axs[1, 0].set_title('Avg Mechanical Energy')

    for ax in axs.flat:
        ax.set_xlabel('Time Step')
        ax.set_ylabel('Energy')
        ax.legend()

    plt.tight_layout()
    plt.show()

# Combined mode: 3D animation and live energy updates
def visualize_combined(data, N_particles, L_values):
    from mpl_toolkits.mplot3d import Axes3D

    steps = len(data)
    fig = plt.figure(figsize=(14, 8))

    # 3D subplot
    ax3d = fig.add_subplot(121, projection='3d')
    scatter = ax3d.scatter([], [], [], c='b', s=20)

    # Energy subplot
    ax_energy = fig.add_subplot(122)
    ax_energy.set_title('Energy Evolution')

    potential_line, = ax_energy.plot([], [], label='Avg Potential Energy', color='r')
    kinetic_line, = ax_energy.plot([], [], label='Avg Kinetic Energy', color='g')
    mechanical_line, = ax_energy.plot([], [], label='Avg Mechanical Energy', color='b')
    ax_energy.legend()
    ax_energy.set_xlabel('Time Step')
    ax_energy.set_ylabel('Energy')

    avg_potential = []
    avg_kinetic = []
    avg_mechanical = []

    def update(frame):
        L = L_values[frame]
        ax3d.set_xlim(0, L)
        ax3d.set_ylim(0, L)
        ax3d.set_zlim(0, L)

        positions = data[frame, :, :3]
        scatter._offsets3d = (positions[:, 0], positions[:, 1], positions[:, 2])
        ax3d.set_title(f"Time step: {frame}")

        potential = data[frame, :, 6]
        kinetic = data[frame, :, 7]
        mechanical = data[frame, :, 8]

        avg_potential.append(np.mean(potential))
        avg_kinetic.append(np.mean(kinetic))
        avg_mechanical.append(np.mean(mechanical))

        potential_line.set_data(range(len(avg_potential)), avg_potential)
        kinetic_line.set_data(range(len(avg_kinetic)), avg_kinetic)
        mechanical_line.set_data(range(len(avg_mechanical)), avg_mechanical)

        ax_energy.set_xlim(0, len(avg_potential))
        ax_energy.set_ylim(min(min(avg_potential), min(avg_kinetic), min(avg_mechanical)), max(max(avg_potential), max(avg_kinetic), max(avg_mechanical)) * 1.1)

    ani = FuncAnimation(fig, update, frames=steps, interval=8)
    plt.show()

# Main script
def main():
    filename = 'trajectory.dat'
    data, N_particles, L_values = load_data(filename)
    print(f"Loaded data with {N_particles} particles per step.")

    while True:
        print("\nSelect mode:")
        print("1: 3D Animation")
        print("2: Energy Visualization")
        print("3: Combined (Animation + Live Energy)")
        print("4: Exit")

        mode = input("Enter mode (1/2/3/4): ").strip()

        if mode == '1':
            visualize_3d_animation(data, N_particles, L_values)
        elif mode == '2':
            visualize_energy(data, N_particles)
        elif mode == '3':
            visualize_combined(data, N_particles, L_values)
        elif mode == '4':
            print("Exiting...")
            break
        else:
            print("Invalid mode. Please enter 1, 2, 3, or 4.")

if __name__ == "__main__":
    main()

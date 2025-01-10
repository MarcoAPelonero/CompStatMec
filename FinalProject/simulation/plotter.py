import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from tqdm import tqdm
import seaborn as sns  # Import Seaborn for enhanced aesthetics

# Set Seaborn style for better aesthetics
sns.set(style="whitegrid")
sns.set_palette("deep")  # You can choose other palettes like "muted", "pastel", "bright", "dark", "colorblind"

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

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    scatter = ax.scatter([], [], [], c='purple', s=50, alpha=0.6, edgecolors='w', linewidth=0.5)

    def update(frame):
        L = L_values[frame]
        ax.set_xlim(0, L)
        ax.set_ylim(0, L)
        ax.set_zlim(0, L)
        ax.set_title(f"Time Step: {frame + 1}/{len(data)}", fontsize=14, color='purple')
        positions = data[frame, :, :3]
        scatter._offsets3d = (positions[:, 0], positions[:, 1], positions[:, 2])

    ani = FuncAnimation(fig, update, frames=len(data), interval=50, blit=False)
    plt.tight_layout()
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

    # Calculate average potential energy after 20% of steps
    cutoff = int(0.2 * steps)
    avg_potential_after_20pct = np.mean(avg_potential[cutoff:])

    time_steps = np.arange(steps)

    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Energy Evolution', fontsize=16, color='purple')

    # Avg Potential Energy Plot
    axs[0, 0].plot(time_steps, avg_potential, label='Avg Potential Energy', color='purple', linewidth=2)
    axs[0, 0].set_title(f'Avg Potential Energy (After 20%: {avg_potential_after_20pct:.2f})', fontsize=12)
    axs[0, 0].set_xlabel('Time Step')
    axs[0, 0].set_ylabel('Energy')
    axs[0, 0].legend()

    # Avg Kinetic Energy Plot
    axs[0, 1].plot(time_steps, avg_kinetic, label='Avg Kinetic Energy', color='teal', linewidth=2)
    axs[0, 1].set_title('Avg Kinetic Energy', fontsize=12)
    axs[0, 1].set_xlabel('Time Step')
    axs[0, 1].set_ylabel('Energy')
    axs[0, 1].legend()

    # Avg Mechanical Energy Plot
    axs[1, 0].plot(time_steps, avg_mechanical, label='Avg Mechanical Energy', color='orange', linewidth=2)
    axs[1, 0].set_title('Avg Mechanical Energy', fontsize=12)
    axs[1, 0].set_xlabel('Time Step')
    axs[1, 0].set_ylabel('Energy')
    axs[1, 0].legend()

    # Hide the 4th subplot
    axs[1, 1].axis('off')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

# Combined mode: 3D animation and live energy updates
def visualize_combined(data, N_particles, L_values):
    from mpl_toolkits.mplot3d import Axes3D

    steps = len(data)
    fig = plt.figure(figsize=(16, 8))
    fig.suptitle('Combined Visualization: 3D Animation & Energy Evolution', fontsize=16, color='purple')

    # 3D subplot
    ax3d = fig.add_subplot(121, projection='3d')
    scatter = ax3d.scatter([], [], [], c='purple', s=50, alpha=0.6, edgecolors='w', linewidth=0.5)

    # Energy subplot
    ax_energy = fig.add_subplot(122)
    ax_energy.set_title('Energy Evolution', fontsize=14, color='purple')
    potential_line, = ax_energy.plot([], [], label='Avg Potential Energy', color='purple', linewidth=2)
    kinetic_line, = ax_energy.plot([], [], label='Avg Kinetic Energy', color='teal', linewidth=2)
    mechanical_line, = ax_energy.plot([], [], label='Avg Mechanical Energy', color='orange', linewidth=2)
    ax_energy.legend()
    ax_energy.set_xlabel('Time Step')
    ax_energy.set_ylabel('Energy')
    ax_energy.grid(True)

    avg_potential = []
    avg_kinetic = []
    avg_mechanical = []

    def update(frame):
        # Update 3D animation
        L = L_values[frame]
        ax3d.set_xlim(0, L)
        ax3d.set_ylim(0, L)
        ax3d.set_zlim(0, L)
        ax3d.set_title(f"Time Step: {frame + 1}/{steps}", fontsize=14, color='purple')
        positions = data[frame, :, :3]
        scatter._offsets3d = (positions[:, 0], positions[:, 1], positions[:, 2])

        # Update Energy plots
        potential = data[frame, :, 6]
        kinetic = data[frame, :, 7]
        mechanical = data[frame, :, 8]

        avg_potential.append(np.mean(potential))
        avg_kinetic.append(np.mean(kinetic))
        avg_mechanical.append(np.mean(mechanical))

        potential_line.set_data(range(len(avg_potential)), avg_potential)
        kinetic_line.set_data(range(len(avg_kinetic)), avg_kinetic)
        mechanical_line.set_data(range(len(avg_mechanical)), avg_mechanical)

        ax_energy.set_xlim(0, len(avg_potential) + 1)
        current_min = min(min(avg_potential, default=0),
                          min(avg_kinetic, default=0),
                          min(avg_mechanical, default=0))
        current_max = max(max(avg_potential, default=1),
                          max(avg_kinetic, default=1),
                          max(avg_mechanical, default=1)) * 1.1
        ax_energy.set_ylim(current_min, current_max)

        return scatter, potential_line, kinetic_line, mechanical_line

    ani = FuncAnimation(fig, update, frames=steps, interval=50, blit=False)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
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
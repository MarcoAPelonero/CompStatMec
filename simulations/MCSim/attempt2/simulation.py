import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

def read_positions(file_path):
    positions = {}
    with open(file_path, 'r') as f:
        current_step = None
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            if line.startswith('Step'):
                current_step = int(line.strip().split()[1])
                positions[current_step] = []
            else:
                positions[current_step].append([float(x) for x in line.strip().split()])
    return positions

def determine_saving_delta(positions):
    steps = sorted(positions.keys())
    if len(steps) > 1:
        return steps[1] - steps[0]
    return None

def read_energies(file_path):
    energies = {}
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            step, energy = line.strip().split()
            energies[int(step)] = float(energy)
    return energies

def animate(i):
    global energy_steps, energy_values
    if i == 0:
        energy_steps = []
        energy_values = []
    
    step = steps[i]
    pos_array = np.array(positions[step])
    
    ax1.cla()
    ax1.scatter(pos_array[:, 0], pos_array[:, 1], pos_array[:, 2])
    ax1.set_title(f'Particle Positions at Step {step}')
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z')
    
    energy_steps.append(step)
    energy_values.append(energies.get(step, 0))
    ax2.cla()
    ax2.plot(energy_steps, energy_values, marker='o')
    ax2.set_title('Energy vs Step')
    ax2.set_xlabel('Step')
    ax2.set_ylabel('Energy')
    ax2.grid(True)

def main():
    global positions, energies, steps, ax1, ax2, energy_steps, energy_values
    positions = read_positions('data\positions_1.txt')
    energies = read_energies('data\energy_1.txt')
    steps = sorted(positions.keys())
    
    energy_steps = []
    energy_values = []
    
    fig = plt.figure(figsize=(12, 6))
    ax1 = fig.add_subplot(121, projection='3d')
    ax2 = fig.add_subplot(122)
    
    ani = FuncAnimation(fig, animate, frames=len(steps), interval=50, repeat=True)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
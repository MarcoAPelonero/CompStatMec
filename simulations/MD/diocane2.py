import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt

def run_simulation(output_file, density, num_particles, dt, num_steps):
    try:
        subprocess.run(["./build/main.exe", output_file, str(density), str(num_particles), str(dt), str(num_steps)], check=True)
    except Exception as e:
        print(f"Error running simulation: {e}")
        raise

def read_simulation_data(filename):
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

def calculate_total_energy(data):
    total_energy = []
    for snapshot in data:
        energy = np.sum(snapshot[:, -1])  # Assuming the last column is the energy
        total_energy.append(energy)
    return np.array(total_energy)

if __name__ == "__main__":
    output_file = "output.txt"
    density = 0.5
    num_particles = 100
    deltaTime = np.linspace(0.001, 0.005, 5)
    num_steps = 10000

    errors = []
    for dt in deltaTime:
        run_simulation(output_file, density, num_particles, dt, num_steps)
        header, data = read_simulation_data(output_file)
        total_energy = calculate_total_energy(data)
        # Calculate error with respect to initial total energy
        energy_error = np.abs(total_energy - total_energy[0])
        # Use the maximum error over all snapshots
        max_error = np.max(energy_error)
        errors.append(max_error)
    
    # Plotting on a log-log scale
    plt.plot(deltaTime, errors, 'o-')
    plt.xlabel('Time Step Size (dt)')
    plt.ylabel('Maximum Energy Error')
    plt.title('Error vs Time Step Size')
    plt.grid(True, which="both", ls="--")
    plt.show()
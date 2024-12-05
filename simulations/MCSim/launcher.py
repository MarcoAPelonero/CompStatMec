import subprocess
import os
import sys
import matplotlib.pyplot as plt
import numpy as np

def compile_simulation():
    print("Compiling the simulation using Makefile...")
    result = subprocess.run(['make'], capture_output=True, text=True)
    if result.returncode != 0:
        print("Compilation failed:")
        print(result.stderr)
        subprocess.run(['make', 'clean'], capture_output=True, text=True)
        sys.exit(1)
    print("Compilation successful.")

def run_simulation(args=None):
    print("Running the simulation...")
    run_command = ['make', 'run']
    if args:
        run_command += args
    result = subprocess.run(run_command, capture_output=True, text=True)
    if result.returncode != 0:
        print("Simulation failed:")
        print(result.stderr)
        sys.exit(1)
    print("Simulation completed.")

def read_energies(file_path):
    if not os.path.exists(file_path):
        print(f"Error: {file_path} does not exist.")
        sys.exit(1)
    
    steps = []
    energies = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                steps.append(int(parts[0]))
                energies.append(float(parts[1]))
    return np.array(steps), np.array(energies)

def plot_energies(steps, energies):
    plt.figure(figsize=(10, 6))
    plt.plot(steps, energies, label='Energy per Step')
    plt.xlabel('Step')
    plt.ylabel('Energy')
    plt.title('Monte Carlo Simulation Energy')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('energy_plot.png')
    plt.show()

def find_stabilization_point(energies, window=100, threshold=1e-5):
    """
    Finds the stabilization point where the moving average change is below the threshold.
    """
    moving_avg = np.convolve(energies, np.ones(window)/window, mode='valid')
    diffs = np.abs(np.diff(moving_avg))
    for i, diff in enumerate(diffs):
        if diff < threshold:
            stabilization_step = i + window
            print(f"Stabilization point detected at step {stabilization_step}.")
            return stabilization_step
    print("No stabilization point detected within the simulation steps.")
    return None

def main():
    compile_simulation()

    run_simulation(['200', '1.0', '6.0', '100']) # 200 particles, 2.0 temperature, 15.0 L, 50000 steps
    run_simulation()

    steps, energies = read_energies('energies.dat')

    plot_energies(steps, energies)

    find_stabilization_point(energies)

if __name__ == "__main__":
    main()

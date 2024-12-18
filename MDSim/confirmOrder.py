import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.stats import linregress

# -----------------------------
# Configuration Parameters
# -----------------------------

# Number of simulation runs per dt to ensure statistical robustness
NUM_RUNS = 1

# Simulation parameters
SIM_TIME = 2500         # Total simulation time
numParticles = 5       # Number of particles
density = 0.1         # Density of the system
temperature = 300      # Temperature in Kelvin
methods = ["EulerCromer"]  # List of numerical methods to evaluate

# Time steps to evaluate
dts = [0.001, 0.003, 0.005,0.007]

# Directory to store simulation data
DATA_DIR = "data"

# Percentage of data to exclude (e.g., first 30%)
EXCLUSION_PERCENTAGE = 0.4

# -----------------------------
# Function Definitions
# -----------------------------

def run_simulation(numParticles, density, temperature=300, dt=0.008, numSteps=10000, 
                   filename="trajectory.dat", method="EulerCromer", pressure=1.0, taup=0.1):
    """
    Executes the simulation by running the external executable with specified parameters.

    Parameters:
        numParticles (int): Number of particles.
        density (float): Density of the system.
        temperature (float): Temperature in Kelvin.
        dt (float): Time step.
        numSteps (int): Number of simulation steps.
        filename (str): Output filename for trajectory data.
        method (str): Numerical integration method.
        pressure (float): Pressure of the system.
        taup (float): Tau parameter for pressure coupling.
    """
    # Construct the command as a list
    command = [
        './build/main.exe',       # Path to the executable
        str(numParticles),        # Number of particles
        str(density),             # Density
        str(dt),                  # Time step
        str(numSteps),            # Number of steps
        filename,                 # Output filename
        method,                   # Numerical method
        str(temperature),         # Temperature
        str(pressure),            # Pressure
        str(taup)                 # Taup parameter
    ]
    
    # Run the simulation
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # Check for errors
    if result.returncode != 0:
        print(f"Simulation failed for dt={dt}. Error message:")
        print(result.stderr.decode())
    else:
        print(f"Simulation completed for dt={dt}.")

def load_data(filename):
    """
    Loads the simulation data from the specified file.

    Returns:
        snapshots (np.ndarray): Array of shape (#snapshots, #particles, 9), where each
                                particle line contains:
                                x, y, z, vx, vy, vz, PE, KE, TE
        L_values (list of float): List of box lengths for each snapshot.
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    snapshots = []
    current_snapshot = []
    L_values = []
    N_particles = None

    i = 0
    total_lines = len(lines)
    while i < total_lines:
        line = lines[i].strip()
        if not line:
            # Empty line indicates end of snapshot
            if current_snapshot:
                snapshots.append(current_snapshot)
                if N_particles is None:
                    N_particles = len(current_snapshot)
                current_snapshot = []
            i += 1
            continue

        # Try reading boxLength
        try:
            L = float(line)
            L_values.append(L)
            i += 1  # Move to the next line to read particle data
        except ValueError:
            # If it's not a boxLength line, it must be particle data
            numbers = list(map(float, line.split()))
            if len(numbers) < 9:
                raise ValueError(f"Expected 9 values per particle, got {len(numbers)} in line: {line}")
            current_snapshot.append(numbers)
            i += 1

    # Add the last snapshot if not already added
    if current_snapshot:
        snapshots.append(current_snapshot)
        if N_particles is None:
            N_particles = len(current_snapshot)

    # Convert to numpy arrays
    snapshots = np.array(snapshots)
    if len(snapshots) == 0:
        raise ValueError("No snapshots found in file.")

    return snapshots, L_values

def compute_average_mechanical_energy(snapshots):
    """
    Computes the average mechanical energy per particle at each snapshot.

    Parameters:
        snapshots (np.ndarray): Array of shape (#snapshots, #particles, 9).

    Returns:
        avg_energy_per_snapshot (np.ndarray): 1D array of average mechanical energy per particle.
    """
    # Extract total energy (TE) which is at index 8
    TE = snapshots[:, :, 8]  # shape (#snapshots, #particles)

    # Sum over all particles to get total mechanical energy per snapshot
    total_energy_per_snapshot = np.sum(TE, axis=1)

    # Compute average per particle
    N = snapshots.shape[1]
    avg_energy_per_snapshot = total_energy_per_snapshot / N
    return avg_energy_per_snapshot

def plot_energy_evolutions(energy_evolutions, dt, exclusion_steps, method):
    """
    Plots the energy evolutions for multiple runs, highlighting excluded and included regions.

    Parameters:
        energy_evolutions (list of np.ndarray): List of average energy arrays for each run.
        dt (float): Time step used in the simulation.
        exclusion_steps (int): Number of steps excluded from analysis.
        method (str): Numerical integration method.
    """
    plt.figure(figsize=(12, 8))
    # Assuming all energy arrays have the same length
    # Use the length of the first run to define time
    time = np.arange(len(energy_evolutions[0])) * dt

    for idx, energy in enumerate(energy_evolutions):
        if len(energy) != len(time):
            print(f"Run {idx+1}: Mismatch in energy and time array lengths ({len(energy)} vs {len(time)}). Skipping plot for this run.")
            continue
        plt.plot(time, energy, label=f'Run {idx+1}')
    
    # Highlight excluded and included regions once
    plt.axvspan(0, exclusion_steps * dt, color='red', alpha=0.1, label='Excluded Region (First 30%)')
    plt.axvspan(exclusion_steps * dt, time[-1], color='green', alpha=0.1, label='Included Region (Last 70%)')

    plt.xlabel("Time", fontsize=14)
    plt.ylabel("Average Mechanical Energy per Particle", fontsize=14)
    plt.title(f"Energy Evolution for dt={dt} using {method}", fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def analyze_method(method, dts, numParticles, density, temperature, SIM_TIME, NUM_RUNS):
    """
    Runs simulations for each dt and computes energy differences to estimate error scaling.

    Parameters:
        method (str): Numerical integration method.
        dts (list of float): List of time steps to evaluate.
        numParticles (int): Number of particles.
        density (float): Density of the system.
        temperature (float): Temperature in Kelvin.
        SIM_TIME (float): Total simulation time.
        NUM_RUNS (int): Number of simulation runs per dt.
    """
    DATA_DIR = "data"
    if not os.path.exists(DATA_DIR):
        os.makedirs(DATA_DIR)

    final_energy_differences = []
    for dt in dts:
        filename = f"{DATA_DIR}/data_{method}_{dt}.txt"
        simSteps = int(SIM_TIME / dt)

        run_differences = []
        energy_evolutions = []  # To store energy data for plotting

        for run in range(NUM_RUNS):
            print(f"Running {method} with dt={dt}, run {run+1}/{NUM_RUNS}")
            run_simulation(
                numParticles=numParticles,
                density=density,
                temperature=temperature,
                dt=dt,
                numSteps=simSteps,
                filename=filename,
                method=method,
                pressure=1.0,
                taup=0.1
            )

            try:
                snapshots, L_values = load_data(filename)
            except Exception as e:
                print(f"Error loading data for dt={dt}, run={run+1}: {e}")
                continue

            avg_energy = compute_average_mechanical_energy(snapshots)

            # Exclude the first 30% of the data
            exclusion_steps = int(EXCLUSION_PERCENTAGE * len(avg_energy))
            if exclusion_steps >= len(avg_energy):
                print(f"Run {run+1}: Exclusion steps ({exclusion_steps}) exceed total steps ({len(avg_energy)}). Skipping this run.")
                continue
            included_energy = avg_energy[exclusion_steps:]

            # Store energy evolution for plotting
            energy_evolutions.append(avg_energy)

            # Measure the difference between final and initial average energy after exclusion
            initial_energy = included_energy[0]
            final_energy = included_energy[-1]
            energy_diff = final_energy - initial_energy
            run_differences.append(energy_diff)

        if not run_differences:
            print(f"No successful runs for dt={dt}. Skipping.")
            continue

        # Average and standard deviation of final energy differences for this dt
        mean_diff = np.mean(run_differences)
        std_diff = np.std(run_differences)
        final_energy_differences.append((dt, mean_diff, std_diff))
        print(f"dt={dt}: Mean energy diff={mean_diff:.6e}, Std={std_diff:.6e}")

        # Plot energy evolutions for this dt
        plot_energy_evolutions(energy_evolutions, dt, exclusion_steps, method)

    # Estimating the order of the algorithm by fitting energy differences vs dt
    dts_values = np.array([x[0] for x in final_energy_differences])
    mean_diffs = np.array([x[1] for x in final_energy_differences])

    if len(dts_values) > 1:
        # Use absolute values to avoid issues with negative differences
        log_dt = np.log(dts_values)
        log_diff = np.log(np.abs(mean_diffs))
        slope, intercept, r_value, p_value, std_err = linregress(log_dt, log_diff)
        print(f"\nMethod: {method}")
        print(f"Estimated Order: {slope:.2f}")
        print(f"R-squared: {r_value**2:.4f}")

        # Plotting the error scaling
        plt.figure(figsize=(10, 7))
        plt.loglog(dts_values, np.abs(mean_diffs), '-o', label=f"{method} (order ≈ {slope:.2f})")
        plt.xlabel("Time Step (dt)", fontsize=14)
        plt.ylabel("Absolute Energy Drift", fontsize=14)
        plt.title("Error Scaling with Time Step", fontsize=16)
        plt.grid(True, which="both", ls="--", lw=0.5)
        plt.legend(fontsize=12)
        plt.tight_layout()
        plt.show()
    else:
        print("Not enough dt values to estimate order.")

def main():
    """
    Main function to perform simulations, process data, estimate algorithm order, and plot results.
    """
    # Ensure the data directory exists
    if not os.path.exists(DATA_DIR):
        os.makedirs(DATA_DIR)

    for method in methods:
        print(f"\nAnalyzing method: {method}")
        analyze_method(method, dts, numParticles, density, temperature, SIM_TIME, NUM_RUNS)

# -----------------------------
# Entry Point
# -----------------------------

if __name__ == "__main__":
    main()

import matplotlib.pyplot as plt
import numpy as np
import subprocess
import os
from scipy.optimize import curve_fit

global TIME 
TIME = 1

def runSimulation(numParticles, density, timeStep, numSteps, filename, simIndex):
    """
    Runs the simulation with the given parameters.
    """
    try:
        print(f"Running simulation {simIndex}...")
        print(f"Parameters: {numParticles} particles, density {density}, time step {timeStep}, {numSteps} steps.")
        subprocess.run(["./build/main.exe", str(numParticles), str(density), str(timeStep), str(numSteps), filename], check=True)
        print("Simulation complete.")
    except Exception as e:
        print(f"Error running simulation: {e}")
        raise

def readEnsembleData(filename):
    """
    Reads the ensemble data file generated by the simulation.
    The file contains multiple snapshots separated by blank lines.
    The first line should be a header: x y z vx vy vz energy.
    """
    print(f"Opening file: {filename}")
    try:
        with open(filename, 'r') as f:
            # Read the file content
            content = f.read().strip()

        # Split into lines
        lines = content.split('\n')

        # Extract the header
        header = lines[0].strip().split()
        print(f"Header: {header}")

        data_lines = lines[1:]
        # Split snapshots by blank lines
        snapshot_blocks = []
        current_block = []
        for line in data_lines:
            line = line.strip()
            if line == '' and current_block:
                snapshot_blocks.append(current_block)
                current_block = []
            elif line != '':
                current_block.append(line)
        # Add the last block if not empty
        if current_block:
            snapshot_blocks.append(current_block)

        # Parse each snapshot into numerical data
        raw_data = []
        for block in snapshot_blocks:
            snapshot_data = [list(map(float, l.split())) for l in block]
            raw_data.append(snapshot_data)

        data_reshaped = np.array(raw_data)
        print(f"Data shape: {data_reshaped.shape}")

    except Exception as e:
        print(f"Error reading file: {e}")
        raise

    return header, data_reshaped

def computeOrder(minDelta, maxDelta, numSimulations, numParticles=100, density=0.5, numSteps=1000):
    deltaTime = np.linspace(minDelta, maxDelta, numSimulations)
    errors = []
    for i, dt in enumerate(deltaTime):
        filename = f"data/ensemble_{i+1}.dat"
        steps = int(TIME / dt)
        print(f"Steps for simulation {i+1}: {steps}")
        runSimulation(numParticles, density, dt, 10000, filename, i+1)
        header, data = readEnsembleData(filename)

        # Extract the energy column index from the header
        energy_index = header.index('energy')

        # Use all data without skipping snapshots
        polished_data = data

        # Compute total energy for each snapshot
        total_energies = np.sum(polished_data[:, :, energy_index], axis=1)


        # Compute energy errors during the stabilized period
        energy_errors = np.std(total_energies[0] - total_energies)
        errors.append((dt, energy_errors))
        print(f"Timestep: {dt}, Average Energy Error (stabilized period): { energy_errors }")

    # Analyze scaling behavior
    if errors:
        dts, errs = zip(*errors)

        # Fit a quadratic curve
        def quadratic_fit(dt, C):
            return C * dt**2

        popt, _ = curve_fit(quadratic_fit, dts, errs)

        # Plot the results
        plt.figure()
        plt.plot(dts, errs, 'bo', label='Energy error')
        plt.plot(dts, quadratic_fit(np.array(dts), *popt), 'r-', label='Quadratic fit')
        plt.xlabel('Timestep (\u0394t)')
        plt.ylabel('Average Energy Error')
        plt.legend()
        plt.show()

    return errors

def find_stabilization_index(total_energies, window_size=10, std_threshold=1e-3):
    for i in range(window_size, len(total_energies)):
        window = total_energies[i-window_size:i]
        mean = np.mean(window)
        std = np.std(window)
        if mean != 0 and (std / mean) < std_threshold:
            return i - window_size
    return None

if __name__ == "__main__":
    if not os.path.exists("data"):
        os.system("mkdir data")

    numParticles = 20
    density = 0.5
    numSteps = 10000
    minDelta = 0.0001
    maxDelta = 0.01
    numSimulations = 5

    errors = computeOrder(minDelta, maxDelta, numSimulations, numParticles, density, numSteps)

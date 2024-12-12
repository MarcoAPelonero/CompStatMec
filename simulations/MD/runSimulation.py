import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import subprocess
import os

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

def computeOrder(minDelta, maxDelta, numSimulations, numParticles=100, density=0.5, numSteps=1000, skipSnapShots=100):
    deltaTime = np.linspace(minDelta, maxDelta, numSimulations)
    errors = []

    for i, dt in enumerate(deltaTime):
        filename = f"data/ensemble_{i+1}.dat"
        runSimulation(numParticles, density, dt, numSteps, filename, i+1)
        header, data = readEnsembleData(filename)

        # Extract the energy column index from the header
        energy_index = header.index('energy')

        # Always skip the first few snapshots
        if data.shape[0] <= skipSnapShots:
            print("Not enough snapshots to skip the desired number. Consider increasing the number of steps.")
            continue

        polished_data = data[skipSnapShots:]  # Skip the initial snapshots

        # Compute total energy for each snapshot in the polished data
        total_energies = np.sum(polished_data[:, :, energy_index], axis=1)

        # Compute the drift as the average difference between consecutive snapshots
        energy_diffs = np.diff(total_energies)
        average_drift = np.mean(np.abs(energy_diffs))
        errors.append((dt, average_drift))
        print(f"Timestep: {dt}, Average Drift: {average_drift}")

    # Check if the error scales with dt by fitting a line
    if errors:
        dts, errs = zip(*errors)
        coeffs = np.polyfit(dts, errs, 1)
        print(f"Linear fit coefficients: {coeffs}")

        plt.figure()
        plt.plot(dts, errs, 'bo', label='Energy drift')
        plt.plot(dts, np.polyval(coeffs, dts), 'r-', label='Linear fit')
        plt.xlabel('Timestep')
        plt.ylabel('Average Drift')
        plt.legend()
        plt.show()

    return errors

if "__main__" == __name__:
    if not os.path.exists("data"):
        os.system("mkdir data")

    numParticles = 50
    density = 0.5
    numSteps = 1000
    minDelta = 0.001
    maxDelta = 0.01
    numSimulations = 5
    skipSnapShots = 100

    errors = computeOrder(minDelta, maxDelta, numSimulations, numParticles, density, numSteps, skipSnapShots)
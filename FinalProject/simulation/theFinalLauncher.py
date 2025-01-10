import os
import subprocess
from pathlib import Path
from itertools import product
from tqdm import tqdm

def run_simulation(rho, T, simulation_executable, output_folder, N=300, dt=0.005, numSteps=2500, p=1.0, taup=2.0, method="MLSpeedVerlet"):
    """
    Runs a single simulation with given rho and T, saving the results to the output_folder.
    
    Parameters:
    - rho (float): Density value.
    - T (float): Temperature value.
    - simulation_executable (Path): Path to the simulation executable.
    - output_folder (Path): Directory to save the simulation results.
    - N (int): Number of particles. Default is 300.
    - dt (float): Time step. Default is 0.005.
    - numSteps (int): Number of simulation steps. Default is 2000.
    - p (float): Pressure. Default is 1.0.
    - taup (float): Pressure relaxation time. Default is 1.0.
    - method (str): Integration method. Default is "MLSpeedVerlet".
    """
    output_folder.mkdir(parents=True, exist_ok=True)
    ml_comparison_path = output_folder / f'mlUnFairComparison_rho_{rho}_T_{T}.dat'
    trajectory_path = output_folder / f'trajectory_rho_{rho}_T_{T}.dat'
    rdf_path = output_folder / f'radialDistribution_rho_{rho}_T_{T}.dat'

    cmd = [
        str(simulation_executable),
        str(N),
        str(rho),
        str(dt),
        str(numSteps),
        str(trajectory_path),
        str(rdf_path),
        method,
        str(T),
        str(p),
        str(taup),
        str(ml_comparison_path)
    ]
    
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return (rho, T, True, "")
    except subprocess.CalledProcessError as e:
        return (rho, T, False, e.stderr.decode())

def main():
    simulation_executable = Path('./build/main.exe')  # Update the path if necessary
    if not simulation_executable.exists():
        print(f"Simulation executable not found at {simulation_executable}")
        return
    
    ml_comparisons_folder = Path('mlcomparisons')
    ml_comparisons_folder.mkdir(exist_ok=True)
    
    # Define ranges for rho and T
    rho_values = [0.2, 0.9]  # Example values, update as needed
    T_values = [0.7, 2.0]    # Example values, update as needed
    
    parameter_grid = list(product(rho_values, T_values))
    
    with tqdm(total=len(parameter_grid), desc='Running Simulations') as pbar:
        for rho, T in parameter_grid:
            rho, T, success, message = run_simulation(rho, T, simulation_executable, ml_comparisons_folder)
            if success:
                tqdm.write(f"Simulation completed for rho={rho}, T={T}")
            else:
                tqdm.write(f"Simulation failed for rho={rho}, T={T}. Error: {message}")
            pbar.update(1)

if __name__ == "__main__":
    main()

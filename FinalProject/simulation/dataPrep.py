import os
import subprocess
import itertools
from tqdm import tqdm
import numpy as np
import concurrent.futures
from pathlib import Path
import logging
import sys
import shutil

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.FileHandler('simulation.log'),
        logging.StreamHandler(sys.stdout)
    ]
)

def read_rdf_file(file_path):
    """
    Reads an RDF file and extracts g(r) data and simulation metadata.

    Parameters
    ----------
    file_path : str or Path
        Path to the RDF file.

    Returns
    -------
    r_values : list of float
        Radii values corresponding to g(r).
    rdfs : numpy.ndarray
        Array of shape (n_steps, n_r), where each row is the g(r) for a step.
    steps : list of int
        List of step indices.
    temperatures : list of float
        Temperature for each step block.
    pressures : list of float
        Pressure for each step block.
    pepp : list of float
        Potential energy per particle for each step block.
    potential_energies : list of float
        Total potential energy for each step block.
    """
    steps = []
    rdfs = []
    temperatures = []
    pressures = []
    pepp = []
    potential_energies = []
    r_values = None

    with open(file_path, 'r') as f:
        current_step_data = []
        current_step = None
        for line in tqdm(f, desc=f'Reading RDF file {file_path}'):
            line = line.strip()
            if line.startswith('# Step'):
                if current_step_data and current_step is not None:
                    rdfs.append([g for _, g in current_step_data])
                    steps.append(current_step)

                parts = line.split()
                if len(parts) < 3:
                    logging.warning(f"Malformed step header in {file_path}: {line}")
                    continue
                try:
                    current_step = int(parts[2])
                except ValueError:
                    logging.warning(f"Invalid step number in {file_path}: {parts[2]}")
                    continue

                try:
                    temp_line = next(f).strip()
                except StopIteration:
                    logging.warning(f"Unexpected EOF after step {current_step} in {file_path}.")
                    break
                temp_parts = temp_line.split()
                if len(temp_parts) < 4:
                    logging.warning(f"Insufficient data after step {current_step} in {file_path}.")
                    break
                try:
                    temperature = float(temp_parts[0])
                    potential_energy_per_particle = float(temp_parts[1])
                    potential_energy = float(temp_parts[2])
                    pressure = float(temp_parts[3])
                except ValueError:
                    logging.warning(
                        f"Non-numeric data in temperature line after step {current_step} in {file_path}."
                    )
                    break
                temperatures.append(temperature)
                pressures.append(pressure)
                pepp.append(potential_energy_per_particle)
                potential_energies.append(potential_energy)
                current_step_data = []
            else:
                if line:
                    vals = line.split()
                    if len(vals) < 2:
                        continue
                    try:
                        r = float(vals[0])
                        g = float(vals[1])
                        current_step_data.append((r, g))
                    except ValueError:
                        logging.warning(f"Non-numeric data encountered in {file_path}: {line}")
                        continue

        if current_step_data and current_step is not None:
            rdfs.append([g for _, g in current_step_data])
            steps.append(current_step)

    if rdfs:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('# Step'):
                    try:
                        next(f)
                    except StopIteration:
                        break
                    break
            r_values = []
            for line in f:
                line = line.strip()
                if line.startswith('# Step'):
                    break
                if line:
                    vals = line.split()
                    if len(vals) >= 2:
                        try:
                            r_val = float(vals[0])
                            r_values.append(r_val)
                        except ValueError:
                            logging.warning(f"Non-numeric r value in {file_path}: {line}")
                            continue
        if not r_values:
            logging.warning(f"No r_values extracted from {file_path}.")
    else:
        r_values = []

    rdfs = np.array(rdfs)
    return r_values, rdfs, steps, temperatures, pressures, pepp, potential_energies


def run_simulation(rho, T, simulation_executable, output_dir, dataset_dir):
    """
    Runs a single simulation for given density (rho) and temperature (T),
    then immediately reads, processes, and stores the RDF data.
    Afterward, removes the raw data files to save space.

    Parameters
    ----------
    rho : float
        Density for the simulation.
    T : float
        Temperature for the simulation.
    simulation_executable : str or Path
        Path to the compiled simulation executable.
    output_dir : str or Path
        Directory to save the intermediate simulation outputs (trajectories, raw RDF).
    dataset_dir : str or Path
        Directory to save the processed RDF data (averages and sampled RDF).
    """
    output_dir = Path(output_dir)
    dataset_dir = Path(dataset_dir)
    dataset_dir.mkdir(parents=True, exist_ok=True)

    rdf_output_file = f'radialDistribution_rho_{rho}_T_{T}.dat'
    rdf_output_path = output_dir / rdf_output_file

    trajectory_filename = f'trajectory_rho_{rho}_T_{T}.dat'
    trajectory_path = output_dir / trajectory_filename

    N = '300'
    rho_str = str(rho)
    dt = '0.008'
    numSteps = '8000'
    fileName = str(trajectory_path)
    rdfName = str(rdf_output_path)
    method = 'RDFSpeedVerlet'
    T_str = str(T)
    p = '0.6'
    taup = '3'

    cmd = [
        str(simulation_executable),
        N,
        rho_str,
        dt,
        numSteps,
        fileName,
        rdfName,
        method,
        T_str,
        p,
        taup
    ]
    logging.info(f"Executing simulation for rho={rho}, T={T} with command: {' '.join(cmd)}")

    try:
        subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        logging.info(f"Simulation completed for rho={rho}, T={T}.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Simulation failed for rho={rho}, T={T}. Error message:\n{e.stderr}")
        return

    if not rdf_output_path.exists():
        logging.error(f"RDF output file {rdf_output_path} not found for rho={rho}, T={T}.")
        return

    try:
        (
            r_values, 
            rdfs, 
            steps, 
            temperatures, 
            pressures, 
            pepp, 
            potential_energies
        ) = read_rdf_file(rdf_output_path)
    except Exception as e:
        logging.error(f"Failed to read RDF file {rdf_output_path} for rho={rho}, T={T}. Error: {e}")
        return

    if rdfs.size == 0:
        logging.error(f"No RDF data found in {rdf_output_path} for rho={rho}, T={T}.")
        return

    total_steps = len(steps)
    if total_steps == 0:
        logging.error(f"No steps found in {rdf_output_path} for rho={rho}, T={T}.")
        return

    selected_indices_avg = np.arange(0, total_steps, 2)
    selected_indices_avg = selected_indices_avg[: total_steps // 2]
    selected_rdfs_avg = rdfs[selected_indices_avg]
    selected_pressures_avg = np.array(pressures)[selected_indices_avg]
    selected_pepp_avg = np.array(pepp)[selected_indices_avg]

    avg_pressure = np.mean(selected_pressures_avg)
    avg_pepp = np.mean(selected_pepp_avg)
    avg_histogram = np.mean(selected_rdfs_avg, axis=0)

    avg_rdf_filename = f'avg_rdf_rho_{rho}_T_{T}.dat'
    avg_rdf_path = dataset_dir / avg_rdf_filename
    try:
        with open(avg_rdf_path, 'w') as f_avg:
            f_avg.write(f"# Average Pressure: {avg_pressure}\n")
            f_avg.write(f"# Average Potential Energy per Particle: {avg_pepp}\n")
            f_avg.write("# r\tg(r)\n")
            for r_val, g_val in zip(r_values, avg_histogram):
                f_avg.write(f"{r_val:.6f}\t{g_val:.6f}\n")
        logging.info(f"Saved average RDF data to {avg_rdf_path}")
    except Exception as e:
        logging.error(f"Failed to save average RDF data to {avg_rdf_path}. Error: {e}")

    remaining_indices = np.setdiff1d(np.arange(total_steps), selected_indices_avg)
    if len(remaining_indices) == 0:
        logging.warning(f"No remaining steps after 50% selection in {rdf_output_path}.")
        return

    remaining_steps = len(remaining_indices)
    desired_samples = remaining_steps // 100
    interval = remaining_steps // desired_samples if desired_samples > 0 else 1
    sampled_remaining_indices = remaining_indices[::interval]

    if sampled_remaining_indices.size == 0:
        logging.warning(f"No data sampled from remaining steps in {rdf_output_path}.")
        return

    sampled_rdfs = rdfs[sampled_remaining_indices]
    sampled_pressures = np.array(pressures)[sampled_remaining_indices]
    sampled_pepp = np.array(pepp)[sampled_remaining_indices]

    rdf_filename = f'rdf_rho_{rho}_T_{T}.dat'
    rdf_path = dataset_dir / rdf_filename
    try:
        with open(rdf_path, 'w') as f_rdf:
            f_rdf.write("# Sampled RDF Data from remaining 50%\n")
            f_rdf.write("# Each block: Pressure, Potential Energy per Particle, then r, g(r)\n\n")
            for idx in range(len(sampled_rdfs)):
                f_rdf.write(f"# Sample {idx + 1}\n")
                f_rdf.write(f"# Pressure: {sampled_pressures[idx]}\n")
                f_rdf.write(f"# Potential Energy per Particle: {sampled_pepp[idx]}\n")
                f_rdf.write("# r\tg(r)\n")
                for r_val, g_val in zip(r_values, sampled_rdfs[idx]):
                    f_rdf.write(f"{r_val:.6f}\t{g_val:.6f}\n")
                f_rdf.write("\n")
        logging.info(f"Saved sampled RDF data to {rdf_path}")
    except Exception as e:
        logging.error(f"Failed to save sampled RDF data to {rdf_path}. Error: {e}")

    # -------------------------------------------------------------------------
    # Remove the raw simulation files to prevent the folder from growing too large
    # -------------------------------------------------------------------------
    try:
        if rdf_output_path.exists():
            rdf_output_path.unlink()  
        if trajectory_path.exists():
            trajectory_path.unlink()  
        logging.info(f"Removed raw data files for rho={rho}, T={T}.")
    except Exception as e:
        logging.warning(f"Could not remove raw data files for rho={rho}, T={T}. Error: {e}")


def main():
    """
    Main entry point of the script. Defines the parameter grid for rho and T,
    runs all simulations (possibly in parallel), and processes each simulation's
    data immediately after completion. Raw data is removed to save disk space.
    """
    num_examples = 50
    rho_values = np.round(np.linspace(0.1, 0.95, num_examples), 2)
    T_values = np.round(np.linspace(0.5, 10, num_examples), 2)

    simulation_executable = Path('./build/main.exe')
    if not simulation_executable.exists():
        logging.error(f"Simulation executable not found at {simulation_executable}.")
        return

    output_dir = Path('rdfData')
    output_dir.mkdir(parents=True, exist_ok=True)

    dataset_dir = 'dataset'

    param_grid = list(itertools.product(rho_values, T_values))
    total_simulations = len(param_grid)

    with concurrent.futures.ProcessPoolExecutor(max_workers=5) as executor:
        futures = []
        for rho, T in param_grid:
            futures.append(
                executor.submit(
                    run_simulation,
                    rho, T,
                    simulation_executable,
                    str(output_dir),
                    dataset_dir
                )
            )

        for _ in tqdm(concurrent.futures.as_completed(futures), total=total_simulations, desc='Running Simulations'):
            pass

    logging.info("All simulations completed.")
    try:
        shutil.rmtree(output_dir)
        logging.info(f"Removed RDF data directory: {output_dir}")
    except Exception as e:
        logging.warning(f"Could not remove folder {output_dir}. Error: {e}")

if __name__ == "__main__":
    main()

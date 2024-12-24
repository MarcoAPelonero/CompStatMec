import os
import subprocess
import itertools
from tqdm import tqdm
import numpy as np
import concurrent.futures
from pathlib import Path
import logging
import sys
import shutil  # Imported for directory deletion

# Setup logging
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
    Reads an RDF file and extracts relevant data.
    
    Parameters:
        file_path (str): Path to the RDF file.
    
    Returns:
        tuple: Contains r_values, rdfs, steps, temperatures, pressures, pepp, potential_energies.
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
                    rdfs.append([g for r, g in current_step_data])
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
                    logging.warning(f"Unexpected end of file after step {current_step} in {file_path}.")
                    break
                temp_parts = temp_line.split()
                if len(temp_parts) < 4:
                    logging.warning(f"Insufficient data in temperature line after step {current_step} in {file_path}.")
                    break
                try:
                    temperature = float(temp_parts[0])
                    potential_energy_per_particle = float(temp_parts[1])
                    potential_energy = float(temp_parts[2])
                    pressure = float(temp_parts[3])
                except ValueError:
                    logging.warning(f"Non-numeric data in temperature line after step {current_step} in {file_path}.")
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

        # Handle the last step
        if current_step_data and current_step is not None:
            rdfs.append([g for r, g in current_step_data])
            steps.append(current_step)

    if rdfs:
        # Extract r_values from the first step
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('# Step'):
                    try:
                        next(f)  # Skip temperature line
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
                            r = float(vals[0])
                            r_values.append(r)
                        except ValueError:
                            logging.warning(f"Non-numeric r value in {file_path}: {line}")
                            continue
        if not r_values:
            logging.warning(f"No r_values extracted from {file_path}.")
    else:
        r_values = []

    rdfs = np.array(rdfs)
    return r_values, rdfs, steps, temperatures, pressures, pepp, potential_energies

def run_simulation(rho, T, simulation_executable, output_dir):
    # Define unique filenames with directory paths
    rdf_output_file = f'radialDistribution_rho_{rho}_T_{T}.dat'
    rdf_output_path = Path(output_dir) / rdf_output_file

    rdf_data_filename = f'rdf_rho_{rho}_T_{T}.dat'
    rdf_data_path = Path(output_dir) / rdf_data_filename

    # Make 'trajectory.dat' unique per simulation
    trajectory_filename = f'trajectory_rho_{rho}_T_{T}.dat'
    trajectory_path = Path(output_dir) / trajectory_filename

    # Define simulation parameters
    N = '300'  # Adjust as needed
    rho_str = str(rho)
    dt = '0.008'
    numSteps = '8000'
    fileName = str(trajectory_path)  # Updated to include directory
    rdfName = str(rdf_output_path)   # Updated to include directory
    method = 'RDFSpeedVerlet'
    T_str = str(T)
    p = '0.6'
    taup = '3'

    # Construct the command to run the simulation
    cmd = [
        simulation_executable,
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

    # Log the command being executed
    logging.info(f"Executing simulation for rho={rho}, T={T} with command: {' '.join(cmd)}")

    # Execute the simulation
    try:
        result = subprocess.run(
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

    # Verify RDF output file existence
    if not rdf_output_path.exists():
        logging.error(f"RDF output file {rdf_output_path} not found for rho={rho}, T={T}.")
        return  
    
    # Parse the RDF output file
    try:
        r_values, rdfs, steps, temperatures, pressures, pepp, potential_energies = read_rdf_file(rdf_output_path)
    except Exception as e:
        logging.error(f"Failed to read RDF file {rdf_output_path} for rho={rho}, T={T}. Error: {e}")
        return 

    # Check if RDF data is valid
    if not r_values or rdfs.size == 0:
        logging.error(f"No RDF data found in {rdf_output_path} for rho={rho}, T={T}.")
        return  

    # Compute average g(r)
    avg_g_r = np.mean(rdfs, axis=0)
    rdf_data = np.column_stack((r_values, avg_g_r))
    header = 'r\tg(r)'

    # Save the averaged RDF data
    try:
        np.savetxt(rdf_data_path, rdf_data, header=header, comments='', fmt='%.6f')
        logging.info(f"Saved averaged RDF data to {rdf_data_path}")
    except Exception as e:
        logging.error(f"Failed to save RDF data for rho={rho}, T={T}. Error: {e}")

def process_rdf_data(output_dir, dataset_dir):
    """
    Processes RDF data by averaging and sampling as specified.

    Parameters:
        output_dir (str): Directory containing the RDF data files.
        dataset_dir (str): Directory to save the processed dataset.
    """
    dataset_path = Path(dataset_dir)
    dataset_path.mkdir(parents=True, exist_ok=True)
    logging.info(f"Created dataset directory at {dataset_path}")

    # Iterate over all radialDistribution files
    rdf_files = list(Path(output_dir).glob('radialDistribution_rho_*_T_*.dat'))
    if not rdf_files:
        logging.warning(f"No RDF data files found in {output_dir}.")
        return

    for rdf_file in tqdm(rdf_files, desc='Processing RDF Data'):
        # Extract rho and T from filename
        filename = rdf_file.stem  # e.g., 'radialDistribution_rho_0.6_T_1.0'
        parts = filename.split('_')
        try:
            rho_index = parts.index('rho') + 1
            T_index = parts.index('T') + 1
            rho = float(parts[rho_index])
            T = float(parts[T_index])
        except (ValueError, IndexError) as e:
            logging.warning(f"Could not parse rho and T from filename {filename}. Skipping file.")
            continue

        # Read RDF data
        try:
            r_values, rdfs, steps, temperatures, pressures, pepp, potential_energies = read_rdf_file(rdf_file)
        except Exception as e:
            logging.error(f"Failed to read RDF file {rdf_file}. Error: {e}")
            continue

        if not rdfs.size:
            logging.warning(f"No RDF data in {rdf_file}. Skipping.")
            continue

        total_steps = len(steps)
        if total_steps == 0:
            logging.warning(f"No steps found in {rdf_file}. Skipping.")
            continue

        # Determine indices to select 50% of the data (e.g., every second step)
        selected_indices_avg = np.arange(0, total_steps, 2)  # Take every second step
        selected_indices_avg = selected_indices_avg[:total_steps // 2]  # Ensure exactly 50%

        # Extract the 50% data for averaging
        selected_rdfs_avg = rdfs[selected_indices_avg]
        selected_pressures_avg = np.array(pressures)[selected_indices_avg]
        selected_pepp_avg = np.array(pepp)[selected_indices_avg]

        # Compute averages for the selected 50%
        avg_pressure = np.mean(selected_pressures_avg)
        avg_pepp = np.mean(selected_pepp_avg)
        avg_histogram = np.mean(selected_rdfs_avg, axis=0)

        # Save to avg_rdf file
        avg_rdf_filename = f'avg_rdf_rho_{rho}_T_{T}.dat'
        avg_rdf_path = dataset_path / avg_rdf_filename
        try:
            with open(avg_rdf_path, 'w') as f_avg:
                f_avg.write(f"# Average Pressure: {avg_pressure}\n")
                f_avg.write(f"# Average Potential Energy per Particle: {avg_pepp}\n")
                f_avg.write("# r\tg(r)\n")
                for r, g in zip(r_values, avg_histogram):
                    f_avg.write(f"{r:.6f}\t{g:.6f}\n")
            logging.info(f"Saved average RDF data to {avg_rdf_path}")
        except Exception as e:
            logging.error(f"Failed to save average RDF data to {avg_rdf_path}. Error: {e}")

        # Determine the remaining 50% indices
        remaining_indices = np.setdiff1d(np.arange(total_steps), selected_indices_avg)
        remaining_steps = len(remaining_indices)
        if remaining_steps == 0:
            logging.warning(f"No remaining steps after selecting 50% in {rdf_file}. Skipping rdf_ file.")
            continue

        # Calculate the step interval to sample approximately every 100 steps
        desired_samples = remaining_steps // 100
        if desired_samples == 0:
            interval = 1  # Sample every step if less than 100 remaining steps
        else:
            interval = remaining_steps // desired_samples

        # Select sampled indices from the remaining data
        sampled_remaining_indices = remaining_indices[::interval]

        if sampled_remaining_indices.size == 0:
            logging.warning(f"No data sampled from remaining steps in {rdf_file}. Skipping rdf_ file.")
            continue

        sampled_rdfs = rdfs[sampled_remaining_indices]
        sampled_pressures = np.array(pressures)[sampled_remaining_indices]
        sampled_pepp = np.array(pepp)[sampled_remaining_indices]

        # Save sampled data to rdf file
        rdf_filename = f'rdf_rho_{rho}_T_{T}.dat'
        rdf_path = dataset_path / rdf_filename
        try:
            with open(rdf_path, 'w') as f_rdf:
                f_rdf.write(f"# Sampled RDF Data from remaining 50%\n")
                f_rdf.write(f"# Each block contains: Pressure, Potential Energy per Particle, followed by r and g(r)\n\n")
                for idx in range(len(sampled_rdfs)):
                    f_rdf.write(f"# Sample {idx + 1}\n")
                    f_rdf.write(f"# Pressure: {sampled_pressures[idx]}\n")
                    f_rdf.write(f"# Potential Energy per Particle: {sampled_pepp[idx]}\n")
                    f_rdf.write("# r\tg(r)\n")
                    for r, g in zip(r_values, sampled_rdfs[idx]):
                        f_rdf.write(f"{r:.6f}\t{g:.6f}\n")
                    f_rdf.write("\n")  # Separate samples with a newline
            logging.info(f"Saved sampled RDF data to {rdf_path}")
        except Exception as e:
            logging.error(f"Failed to save sampled RDF data to {rdf_path}. Error: {e}")

    # After processing all files, delete rdfData folder
    try:
        shutil.rmtree(output_dir)
        logging.info(f"Deleted RDF data directory {output_dir}.")
    except Exception as e:
        logging.error(f"Failed to delete RDF data directory {output_dir}. Error: {e}")

def main():
    # Define the ranges for rho and T
    rho_values = [0.3,0.4,0.5,0.6, 0.7, 0.8,0.9]  
    T_values = [1.0, 1.2, 1.4,1.6,1.8,2.0]  
    
    # Path to the simulation executable
    simulation_executable = Path('./build/main.exe')  # Adjust as needed

    # Check if the simulation executable exists
    if not simulation_executable.exists():
        logging.error(f"Simulation executable not found at {simulation_executable}. Please check the path.")
        return

    # Create output directory if it doesn't exist
    output_dir = Path('rdfData')
    output_dir.mkdir(parents=True, exist_ok=True)

    # Calculate total number of simulations
    total_simulations = len(rho_values) * len(T_values)

    # Use ProcessPoolExecutor for parallel execution
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [
            executor.submit(run_simulation, rho, T, str(simulation_executable), str(output_dir))
            for rho, T in itertools.product(rho_values, T_values)
        ]
        # Track progress with tqdm
        for future in tqdm(concurrent.futures.as_completed(futures), total=total_simulations, desc='Running Simulations'):
            try:
                future.result()
            except Exception as e:
                logging.error(f"An unexpected error occurred during simulation: {e}")

    # After all simulations are complete, process the RDF data
    dataset_dir = 'dataset'
    process_rdf_data(str(output_dir), dataset_dir)
    logging.info("Data processing complete.")

if __name__ == "__main__":
    main()

import os
import subprocess
import itertools
from tqdm import tqdm
import numpy as np

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
        for line in tqdm(f, desc='Reading RDF file'):
            line = line.strip()
            if line.startswith('# Step'):
                if current_step_data and current_step is not None:
                    rdfs.append([g for r, g in current_step_data])
                    steps.append(current_step)
                parts = line.split()
                current_step = int(parts[2])
                try:
                    temp_line = next(f).strip()
                except StopIteration:
                    print(f"Unexpected end of file after step {current_step}.")
                    break
                temp_parts = temp_line.split()
                if len(temp_parts) < 4:
                    print(f"Insufficient data in temperature line after step {current_step}.")
                    break
                temperature = float(temp_parts[0])
                pressure = float(temp_parts[3])
                potential_energy_per_particle = float(temp_parts[1])
                potential_energy = float(temp_parts[2])
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
                        continue  

        # Handle the last step
        if current_step_data and current_step is not None:
            rdfs.append([g for r, g in current_step_data])
            steps.append(current_step)

    if rdfs:
        r_values = []
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('# Step'):
                    try:
                        next(f)
                    except StopIteration:
                        break
                    break
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
                            continue
            if r_values:
                r_values = r_values
            else:
                r_values = []
    else:
        r_values = []

    rdfs = np.array(rdfs)
    return r_values, rdfs, steps, temperatures, pressures, pepp, potential_energies

def main():
    rho_values = [0.6, 0.7, 0.8]  
    T_values = [1.0, 1.2, 1.4]   

    simulation_executable = './build/main.exe'  

    if not os.path.isfile(simulation_executable):
        print(f"Simulation executable not found at {simulation_executable}. Please check the path.")
        return

    output_dir = 'rdfData'
    os.makedirs(output_dir, exist_ok=True)

    total_simulations = len(rho_values) * len(T_values)

    for rho, T in tqdm(itertools.product(rho_values, T_values), total=total_simulations, desc='Running Simulations'):
        rdf_output_file = f'radialDistribution_rho_{rho}_T_{T}.dat'
        rdf_output_path = os.path.join(output_dir, rdf_output_file)

        rdf_data_filename = f'rdf_rho_{rho}_T_{T}.dat'
        rdf_data_path = os.path.join(output_dir, rdf_data_filename)

        # Prepare command-line arguments based on the C++ program's expectations
        # argv[1] = N (default: 800)
        # argv[2] = rho
        # argv[3] = dt (default: 0.008)
        # argv[4] = numSteps (default: 3000)
        # argv[5] = fileName (default: "trajectory.dat")
        # argv[6] = rdfName (unique for each run)
        # argv[7] = method (default: "ThermoEulerCromer")
        # argv[8] = T
        # argv[9] = p (default: 1.011)
        # argv[10] = taup (default: 1)

        N = '800'
        rho_str = str(rho)
        dt = '0.008'
        numSteps = '3000'
        fileName = 'trajectory.dat'
        rdfName = rdf_output_file  
        method = 'ThermoEulerCromer'
        T_str = str(T)
        p = '1.0'
        taup = '3'

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

        try:
            result = subprocess.run(
                cmd,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            print(f"Simulation completed for rho={rho}, T={T}.")
        except subprocess.CalledProcessError as e:
            print(f"Simulation failed for rho={rho}, T={T}. Error message:\n{e.stderr}")
            continue 

        simulation_rdf_path = rdf_output_path

        if not os.path.exists(simulation_rdf_path):
            print(f"RDF output file {simulation_rdf_path} not found for rho={rho}, T={T}.")
            continue  
        try:
            r_values, rdfs, steps, temperatures, pressures, pepp, potential_energies = read_rdf_file(simulation_rdf_path)
        except Exception as e:
            print(f"Failed to read RDF file {simulation_rdf_path} for rho={rho}, T={T}. Error: {e}")
            continue 
        if not r_values or rdfs.size == 0:
            print(f"No RDF data found in {simulation_rdf_path} for rho={rho}, T={T}.")
            continue  

        avg_g_r = np.mean(rdfs, axis=0)

        rdf_data = np.column_stack((r_values, avg_g_r))

        header = 'r\tg(r)'
        try:
            np.savetxt(rdf_data_path, rdf_data, header=header, comments='', fmt='%.6f')
            print(f"Saved RDF data to {rdf_data_path}")
        except Exception as e:
            print(f"Failed to save RDF data for rho={rho}, T={T}. Error: {e}")

if __name__ == "__main__":
    main()

from tqdm import tqdm
import numpy as np

def read_rdf_file(file_path):
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
                    # Read the next line for temperature, pressure, potential energy
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
                        continue  # Skip malformed lines
                    try:
                        r = float(vals[0])
                        g = float(vals[1])
                        current_step_data.append((r, g))
                    except ValueError:
                        continue  # Skip lines with non-numeric data

        # Handle the last step
        if current_step_data and current_step is not None:
            rdfs.append([g for r, g in current_step_data])
            steps.append(current_step)

    # Extract r_values from the first step for consistency
    if rdfs:
        r_values = [r for r, g in current_step_data]  # Assuming r is the same across all steps
    else:
        r_values = []

    rdfs = np.array(rdfs)
    return r_values, rdfs, steps, temperatures, pressures, pepp, potential_energies
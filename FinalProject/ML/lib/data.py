# lib/data.py

from pathlib import Path
from tqdm import tqdm
import os
import numpy as np
import re

datasetFolder = "dataset"

class Sample:
    """
    Represents a single RDF sample with its associated properties.
    """
    def __init__(self, rho, T, pressure, energy, r, rdf):
        self.rho = rho
        self.T = T
        self.pressure = pressure
        self.energy = energy
        self.r = r
        self.rdf = rdf

    def show(self, single=True):
        import matplotlib.pyplot as plt
        plt.plot(self.r, self.rdf, '-')
        if single:
            plt.xlabel('r')
            plt.ylabel('g(r)')
            plt.title(f'RDF for rho={self.rho}, T={self.T}\nPressure={self.pressure}, Energy={self.energy}')
            plt.grid(True)
            plt.show()

class RDFFile:
    """
    Handles reading RDF data from a single file, which may contain multiple samples.
    """
    def __init__(self, filename):
        self.filename = filename

    def parse_filename(self):
        """
        Extracts rho and T values from the filename using regex for robustness.
        Expected filename format: rdf_rho_<rho>_T_<T>.dat
        Example: rdf_rho_0.9_T_1.0.dat
        """
        pattern = r'rdf_rho_(?P<rho>[\d\.]+)_T_(?P<T>[\d\.]+)\.dat'
        match = re.match(pattern, self.filename)
        if not match:
            raise ValueError(f"Filename format incorrect: {self.filename}")
        try:
            rho = float(match.group('rho'))
            T = float(match.group('T'))
        except ValueError as e:
            raise ValueError(f"Could not convert rho or T to float in filename: {self.filename}") from e
        return rho, T

    def readFile(self):
        rho, T = self.parse_filename()

        filepath = os.path.join(datasetFolder, self.filename)
        if not os.path.isfile(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")

        samples = []

        with open(filepath, 'r') as f:
            lines = f.readlines()

        i = 0
        total_lines = len(lines)
        while i < total_lines:
            line = lines[i].strip()
            # Only match lines that start with '# Sample' followed by a number
            if re.match(r'^# Sample\s+\d+', line):
                # Parse Pressure
                if i + 1 >= total_lines:
                    raise ValueError(f"Unexpected end of file after '# Sample' at line {i+1} in {self.filename}.")
                pressure_line = lines[i + 1].strip()
                pressure_pattern = r'#\s*Pressure:\s*([-\d\.]+)'
                pressure_match = re.match(pressure_pattern, pressure_line)
                if not pressure_match:
                    raise ValueError(f"Expected Pressure line after Sample header at line {i+2} in {self.filename}. Found: {pressure_line}")
                try:
                    pressure = float(pressure_match.group(1))
                except ValueError as e:
                    raise ValueError(f"Invalid Pressure value at line {i+2} in {self.filename}.") from e

                # Parse Potential Energy
                if i + 2 >= total_lines:
                    raise ValueError(f"Unexpected end of file after Pressure line at line {i+2} in {self.filename}.")
                energy_line = lines[i + 2].strip()
                energy_pattern = r'#\s*Potential Energy per Particle:\s*([-\d\.]+)'
                energy_match = re.match(energy_pattern, energy_line)
                if not energy_match:
                    raise ValueError(f"Expected Potential Energy line after Pressure line at line {i+3} in {self.filename}. Found: {energy_line}")
                try:
                    energy = float(energy_match.group(1))
                except ValueError as e:
                    raise ValueError(f"Invalid Potential Energy value at line {i+3} in {self.filename}.") from e

                # Expecting header line for r and g(r)
                if i + 3 >= total_lines:
                    raise ValueError(f"Unexpected end of file after Potential Energy line at line {i+3} in {self.filename}.")
                header_line = lines[i + 3].strip()
                header_pattern = r'#\s*r\s+g\(r\)'
                if not re.match(header_pattern, header_line):
                    raise ValueError(f"Expected header line for r and g(r) at line {i+4} in {self.filename}. Found: {header_line}")

                i += 4  # Move to the first data line after header

                r = []
                rdf = []
                # Read data lines until the next sample or end of file
                while i < total_lines:
                    data_line = lines[i].strip()
                    if re.match(r'^# Sample\s+\d+', data_line):
                        break  # Next sample starts
                    if not data_line or data_line.startswith('#'):
                        i += 1
                        continue  # Skip empty or comment lines

                    # Split by any whitespace
                    parts = re.split(r'\s+', data_line)
                    if len(parts) < 2:
                        print(f"Warning: Skipping invalid line {i+1} in {self.filename}: '{data_line}'")
                        i += 1
                        continue

                    try:
                        r_val = float(parts[0])
                        rdf_val = float(parts[1])
                        r.append(r_val)
                        rdf.append(rdf_val)
                    except ValueError:
                        print(f"Warning: Skipping non-numeric line {i+1} in {self.filename}: '{data_line}'")
                    i += 1

                if not r or not rdf:
                    raise ValueError(f"No valid RDF data found for a sample starting at line {i+1} in file: {filepath}")

                sample = Sample(
                    rho=rho,
                    T=T,
                    pressure=pressure,
                    energy=energy,
                    r=np.array(r),
                    rdf=np.array(rdf)
                )
                samples.append(sample)
            else:
                i += 1  # Move to next line if not a sample header

        if not samples:
            raise ValueError(f"No samples found in file: {filepath}")

        return samples

class Dataset:
    def __init__(self, dataset_folder, filetype='rdf_rho_*_T_*.dat'):
        self.filetype = filetype
        self.dataset_folder = dataset_folder
        self.data = []
        self.train = []
        self.vali = []
        self.test = []

    def read_dataset(self):
        files = list(Path(self.dataset_folder).glob(self.filetype))
        if not files:
            raise FileNotFoundError(f"No files found with pattern {self.filetype} in {self.dataset_folder}")

        for rdf_file in tqdm(files, desc='Processing RDF Data'):
            filename = rdf_file.name
            rdf_handler = RDFFile(filename)
            try:
                samples = rdf_handler.readFile()
                self.data.extend(samples)
            except Exception as e:
                print(f"Error reading {filename}: {e}")

        if not self.data:
            raise ValueError("No valid RDF data was loaded.")

        self.data = np.array(self.data, dtype=object)

    def split(self, train=0.8, vali=0.1, test=0.1, random_state=None):
        if not np.isclose(train + vali + test, 1.0):
            raise ValueError("Percentages for train, vali, and test must sum to 1.0")

        if train < 0 or vali < 0 or test < 0:
            raise ValueError("Split proportions must be non-negative")

        data_length = len(self.data)
        indices = np.arange(data_length)
        np.random.seed(random_state)
        np.random.shuffle(indices)

        train_end = int(train * data_length)
        vali_end = train_end + int(vali * data_length)

        self.train = self.data[indices[:train_end]]
        self.vali = self.data[indices[train_end:vali_end]] if vali > 0 else np.array([], dtype=object)
        self.test = self.data[indices[vali_end:]] if test > 0 else np.array([], dtype=object)

        print(f"Dataset split into {len(self.train)} training, {len(self.vali)} validation, and {len(self.test)} testing samples.")

    def get_statistics(self):
        def average_param(data, param):
            if len(data) == 0:
                return float('nan')
            return np.mean([getattr(d, param) for d in data])

        def mean_or_nan(data_list):
            return np.mean(data_list) if data_list else float('nan')

        stats = {
            'train': {
                'rho_mean': average_param(self.train, 'rho'),
                'T_mean': average_param(self.train, 'T'),
                'pressure_mean': mean_or_nan([d.pressure for d in self.train]),
                'energy_mean': mean_or_nan([d.energy for d in self.train]),
            },
            'vali': {
                'rho_mean': average_param(self.vali, 'rho'),
                'T_mean': average_param(self.vali, 'T'),
                'pressure_mean': mean_or_nan([d.pressure for d in self.vali]) if len(self.vali) > 0 else float('nan'),
                'energy_mean': mean_or_nan([d.energy for d in self.vali]) if len(self.vali) > 0 else float('nan'),
            },
            'test': {
                'rho_mean': average_param(self.test, 'rho'),
                'T_mean': average_param(self.test, 'T'),
                'pressure_mean': mean_or_nan([d.pressure for d in self.test]) if len(self.test) > 0 else float('nan'),
                'energy_mean': mean_or_nan([d.energy for d in self.test]) if len(self.test) > 0 else float('nan'),
            },
        }
        return stats

    def plot_sample(self, split='train', index=0):
        split_data = getattr(self, split, None)
        if split_data is None:
            raise ValueError(f"Invalid split name: {split}. Choose from 'train', 'vali', 'test'.")

        if index >= len(split_data):
            raise IndexError(f"Index {index} out of range for split '{split}' with size {len(split_data)}.")

        sample = split_data[index]
        sample.show(single=True)

    def __len__(self):  
        return len(self.data)

if __name__ == '__main__':
    # Test RDF reading and plotting with a specific file
    test_filename = 'rdf_rho_0.9_T_1.0.dat'
    rdf_file = RDFFile(test_filename)
    
    try:
        samples = rdf_file.readFile()
        for idx, sample in enumerate(samples, start=1):
            print(f"Sample {idx}:")
            print(f"  rho: {sample.rho}")
            print(f"  T: {sample.T}")
            print(f"  Pressure: {sample.pressure}")
            print(f"  Energy: {sample.energy}")
            # Uncomment the next line to plot each sample
            # sample.show(single=True)
    except Exception as e:
        print(f"Error: {e}")

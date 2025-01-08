from pathlib import Path
from tqdm import tqdm
import os
import matplotlib.pyplot as plt
import numpy as np
from sklearn.model_selection import train_test_split

datasetFolder = "dataset"

class RDF:
    def __init__(self, filename):
        self.filename = filename 

    def readFile(self):
        parts = self.filename.split('_')
        try:
            rho_index = parts.index('rho') + 1
            T_index = parts.index('T') + 1
        except ValueError as e:
            raise ValueError(f"Filename format incorrect: {self.filename}") from e

        r = []
        rdf = []

        try:
            self.rho = float(parts[rho_index])
            self.T = float(parts[T_index][:3])
        except (IndexError, ValueError) as e:
            raise ValueError(f"Could not parse rho or T from filename: {self.filename}") from e

        filepath = os.path.join(datasetFolder, self.filename)
        if not os.path.isfile(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")

        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('# Average Pressure'):
                    self.avg_pressure = float(line.split()[-1])
                    continue
                if line.startswith('# Average Potential'):
                    self.avg_energy = float(line.split()[-1])
                    continue
                if line.startswith('# r\tg(r)'):
                    continue
                try:
                    r_val, rdf_val = map(float, line.strip().split())
                    r.append(r_val)
                    rdf.append(rdf_val)
                except ValueError:
                    continue

        if not r or not rdf:
            raise ValueError(f"No valid RDF data found in file: {filepath}")

        self.r = np.array(r) 
        self.rdf = np.array(rdf)

    def show(self, single=True):
        plt.plot(self.r, self.rdf, '-') 
        if single:
            plt.xlabel('r')
            plt.ylabel('g(r)')
            plt.title(f'RDF for rho={self.rho}, T={self.T}')
            plt.grid(True)
            plt.show()


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
            element = RDF(filename)
            try:
                element.readFile()
                self.data.append(element)
            except Exception as e:
                print(f"Error reading {filename}: {e}")

        if not self.data:
            raise ValueError("No valid RDF data was loaded.")

        self.data = np.array(self.data, dtype=object) 

    def split(self, train=0.8, vali=0.0, test=0.2, random_state=None):
        if not np.isclose(train + vali + test, 1.0):
            raise ValueError("Percentages for train, vali, and test must sum to 1.0")

        train_data, temp_data = train_test_split(
            self.data, test_size=(vali + test), random_state=random_state
        )
        
        vali_size = vali / (vali + test)  
        vali_data, test_data = train_test_split(
            temp_data, test_size=(test), random_state=random_state
        )

        self.train = train_data
        self.vali = vali_data
        self.test = test_data

        print(f"Dataset split into {len(self.train)} training, {len(self.vali)} validation, and {len(self.test)} testing samples.")

    def get_statistics(self):
        def average_param(data, param):
            return np.mean([getattr(d, param) for d in data])

        stats = {
            'train': {
                'rho_mean': average_param(self.train, 'rho'),
                'T_mean': average_param(self.train, 'T'),
                'avg_pressure_mean': np.mean([d.avg_pressure for d in self.train]),
                'avg_energy_mean': np.mean([d.avg_energy for d in self.train]),
            },
            'vali': {
                'rho_mean': average_param(self.vali, 'rho'),
                'T_mean': average_param(self.vali, 'T'),
                'avg_pressure_mean': np.mean([d.avg_pressure for d in self.vali]),
                'avg_energy_mean': np.mean([d.avg_energy for d in self.vali]),
            },
            'test': {
                'rho_mean': average_param(self.test, 'rho'),
                'T_mean': average_param(self.test, 'T'),
                'avg_pressure_mean': np.mean([d.avg_pressure for d in self.test]),
                'avg_energy_mean': np.mean([d.avg_energy for d in self.test]),
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


if __name__ == "__main__":
    dataset = Dataset(dataset_folder=datasetFolder, filetype='avg_rdf_rho_*_T_*.dat')

    dataset.read_dataset()

    dataset.split(train=0.5, vali=0.0, test=0.5, random_state=42)

    stats = dataset.get_statistics()
    for split, split_stats in stats.items():
        print(f"\n{split.capitalize()} Set Statistics:")
        for key, value in split_stats.items():
            print(f"  {key}: {value:.3f}")

    dataset.plot_sample(split='train', index=0)

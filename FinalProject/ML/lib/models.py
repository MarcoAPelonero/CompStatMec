import numpy as np
import torch
from torch.utils.data import Dataset
import torch.nn as nn

class RDFDensityDataset(Dataset):
    """
    PyTorch Dataset for (RDF, density) -> energy.
    We assume all RDF arrays have been zero-padded or truncated to the same length.
    But typically we pick a fixed bin length, e.g. 1000, across all samples.
    """
    def __init__(self, samples, max_len=1000):
        """
        samples: list of Sample objects (from data.py).
                 Each has .rdf, .rho, .energy
        max_len: how many bins we use for the RDF (truncate or zero-pad).
        """
        self.samples = samples
        self.max_len = max_len

        self.X_rdf, self.X_density, self.y = self._build_dataset()

    def _pad_rdf(self, rdf):
        """Zero-pad or truncate to self.max_len."""
        if len(rdf) >= self.max_len:
            return rdf[:self.max_len]
        else:
            arr = np.zeros(self.max_len, dtype=np.float32)
            arr[:len(rdf)] = rdf
            return arr

    def _build_dataset(self):
        X_rdf_list = []
        X_density_list = []
        y_list = []
        for sample in self.samples:
            rdf_fixed = self._pad_rdf(sample.rdf)
            X_rdf_list.append(rdf_fixed)
            X_density_list.append(sample.rho)
            y_list.append(sample.energy)

        X_rdf_array = np.array(X_rdf_list, dtype=np.float32)    
        X_density_array = np.array(X_density_list, dtype=np.float32).reshape(-1, 1) 
        y_array = np.array(y_list, dtype=np.float32).reshape(-1, 1)
        return X_rdf_array, X_density_array, y_array

    def __len__(self):
        return len(self.samples)

    def __getitem__(self, idx):
        rdf = torch.from_numpy(self.X_rdf[idx])         
        density = torch.from_numpy(self.X_density[idx]) 
        y = torch.from_numpy(self.y[idx])               
        return (rdf, density, y)


class RDFDensityMLP(nn.Module):
    """
    Two-branch MLP:
      - Branch 1: RDF (e.g. 1000 bins)
      - Branch 2: single density
    We merge them and output energy.
    """

    def __init__(self, rdf_size=1000, hidden_dim=128):
        super().__init__()
        self.rdf_branch = nn.Sequential(
            nn.Linear(rdf_size, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU()
        )
        self.density_branch = nn.Sequential(
            nn.Linear(1, hidden_dim // 2),
            nn.ReLU(),
            nn.Linear(hidden_dim // 2, hidden_dim // 2),
            nn.ReLU()
        )
        merged_dim = hidden_dim + (hidden_dim // 2)
        self.merged_head = nn.Sequential(
            nn.Linear(merged_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, 1)  
        )

    def forward(self, x_rdf, x_density):
        out_rdf = self.rdf_branch(x_rdf) 
        out_density = self.density_branch(x_density)  
        merged = torch.cat([out_rdf, out_density], dim=1)
        return self.merged_head(merged)

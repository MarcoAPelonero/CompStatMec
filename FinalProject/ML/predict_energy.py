import torch
import numpy as np
import argparse
import os
from lib.data import RDFFile 
from lib.models import RDFDensityMLP  

def pad_rdf(rdf, max_len=1000):
    """
    Pads or truncates the RDF array to a fixed length.

    Parameters:
    - rdf (list or np.array): The RDF data.
    - max_len (int): The desired fixed length.

    Returns:
    - np.array: Padded or truncated RDF.
    """
    rdf = np.array(rdf, dtype=np.float32)
    if len(rdf) >= max_len:
        return rdf[:max_len]
    else:
        padded = np.zeros(max_len, dtype=np.float32)
        padded[:len(rdf)] = rdf
        return padded

def load_model(model_path, rdf_size=1000, hidden_dim=128, device='cpu'):
    """
    Loads the averaged model from the given path.

    Parameters:
    - model_path (str): Path to the saved averaged model.
    - rdf_size (int): Size of the RDF input.
    - hidden_dim (int): Hidden dimension size used during training.
    - device (str): Device to load the model on ('cpu' or 'cuda').

    Returns:
    - model (RDFDensityMLP): The loaded PyTorch model.
    """
    if not os.path.isfile(model_path):
        raise FileNotFoundError(f"Model file not found: {model_path}")

    model = RDFDensityMLP(rdf_size=rdf_size, hidden_dim=hidden_dim)
    model.load_state_dict(torch.load(model_path, map_location=device))
    model.to(device)
    model.eval()
    return model

def main():
    data_file = "rdf_rho_0.2_T_1.0.dat"
    model_path = "best_model.pth"
    max_len = 1000

    dataset_folder = "dataset"
    data_file_path = os.path.join(dataset_folder, data_file)

    if not os.path.isfile(data_file_path):
        print(f"Data file not found: {data_file_path}")
        return

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    try:
        model = load_model(model_path, rdf_size=max_len, hidden_dim=128, device=device)
        print(f"Loaded model from {model_path}")
    except Exception as e:
        print(f"Error loading model: {e}")
        return

    rdf_file_handler = RDFFile(data_file)
    try:
        samples = rdf_file_handler.readFile()
    except Exception as e:
        print(f"Error reading data file: {e}")
        return

    if not samples:
        print("No samples found in the data file.")
        return

    print(f"Loaded {len(samples)} sample(s) from {data_file}.\n")

    print(f"{'Sample':<10} {'True Energy':<15} {'Predicted Energy':<20}")
    print("-" * 45)
    for idx, sample in enumerate(samples, start=1):
        # Pad or truncate RDF
        rdf_padded = pad_rdf(sample.rdf, max_len=max_len)  
        rho = sample.rho  

        rdf_tensor = torch.from_numpy(rdf_padded).to(device, dtype=torch.float).unsqueeze(0)  
        rho_tensor = torch.tensor([[rho]], dtype=torch.float).to(device)  

        with torch.no_grad():
            predicted_energy = model(rdf_tensor, rho_tensor).item()

        true_energy = sample.energy

        print(f"{idx:<10} {true_energy:<15.4f} {predicted_energy:<20.4f}")

    print("\nPrediction completed.")

if __name__ == "__main__":
    main()
# scripts/predict_energy.py

import torch
import numpy as np
import sys
import json
import os
import time

class RDFDensityMLP(torch.nn.Module):
    """
    Two-branch MLP:
      - Branch 1: RDF (e.g. 1000 bins)
      - Branch 2: single density
    We merge them and output energy.
    """

    def __init__(self, rdf_size=1000, hidden_dim=128):
        super().__init__()
        # 1) RDF branch
        self.rdf_branch = torch.nn.Sequential(
            torch.nn.Linear(rdf_size, hidden_dim),
            torch.nn.ReLU(),
            torch.nn.Linear(hidden_dim, hidden_dim),
            torch.nn.ReLU()
        )
        # 2) Density branch
        self.density_branch = torch.nn.Sequential(
            torch.nn.Linear(1, hidden_dim // 2),
            torch.nn.ReLU(),
            torch.nn.Linear(hidden_dim // 2, hidden_dim // 2),
            torch.nn.ReLU()
        )
        # 3) Merge the two branches
        merged_dim = hidden_dim + (hidden_dim // 2)
        self.merged_head = torch.nn.Sequential(
            torch.nn.Linear(merged_dim, hidden_dim),
            torch.nn.ReLU(),
            torch.nn.Linear(hidden_dim, 1)  # final energy prediction
        )

    def forward(self, x_rdf, x_density):
        # Forward pass RDF
        out_rdf = self.rdf_branch(x_rdf)           # shape: [batch, hidden_dim]
        # Forward pass density
        out_density = self.density_branch(x_density)  # shape: [batch, hidden_dim//2]
        # Concat
        merged = torch.cat([out_rdf, out_density], dim=1)
        # Final
        return self.merged_head(merged)

def main():
    """
    Expects a JSON input from stdin with 'rdf' and 'density' keys.
    Outputs the predicted energy to stdout.
    """
    if len(sys.argv) != 2:
        print("Usage: python predict_energy.py <path_to_model.pth>", file=sys.stderr)
        sys.exit(1)

    model_path = sys.argv[1]

    if not os.path.exists(model_path):
        print(f"Error: Model file '{model_path}' does not exist.", file=sys.stderr)
        sys.exit(1)

    # Load the model
    model = RDFDensityMLP()
    model.load_state_dict(torch.load(model_path, map_location=torch.device('cpu'), weights_only=True))
    model.eval()

    # Read JSON input from stdin
    
    try:
        input_data = sys.stdin.read()
        input_json = json.loads(input_data)
        rdf = input_json['rdf']          # Should be a list of 1000 floats
        density = input_json['density']  # Should be a single float
    except Exception as e:
        print(f"Error reading input data: {e}", file=sys.stderr)
        sys.exit(1)

    # Validate input sizes
    if len(rdf) != 1000:
        print(f"Error: RDF histogram must have 1000 bins, but got {len(rdf)}.", file=sys.stderr)
        sys.exit(1)
    if not isinstance(density, float) and not isinstance(density, int):
        print(f"Error: Density must be a float, but got {type(density)}.", file=sys.stderr)
        sys.exit(1)

    # Convert inputs to tensors
    rdf_tensor = torch.tensor([rdf], dtype=torch.float32)       # Shape: [1, 1000]
    density_tensor = torch.tensor([[density]], dtype=torch.float32)  # Shape: [1, 1]

    start_time = time.perf_counter()
    with torch.no_grad():
        predicted_energy = model(rdf_tensor, density_tensor).item()
    inference_time = time.perf_counter() - start_time

    # Output both energy and timing
    print(f"{predicted_energy} {inference_time}")

if __name__ == "__main__":
    main()

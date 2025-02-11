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
        self.rdf_branch = torch.nn.Sequential(
            torch.nn.Linear(rdf_size, hidden_dim),
            torch.nn.ReLU(),
            torch.nn.Linear(hidden_dim, hidden_dim),
            torch.nn.ReLU()
        )
        self.density_branch = torch.nn.Sequential(
            torch.nn.Linear(1, hidden_dim // 2),
            torch.nn.ReLU(),
            torch.nn.Linear(hidden_dim // 2, hidden_dim // 2),
            torch.nn.ReLU()
        )
        merged_dim = hidden_dim + (hidden_dim // 2)
        self.merged_head = torch.nn.Sequential(
            torch.nn.Linear(merged_dim, hidden_dim),
            torch.nn.ReLU(),
            torch.nn.Linear(hidden_dim, 1) 
        )

    def forward(self, x_rdf, x_density):
        out_rdf = self.rdf_branch(x_rdf)   
        out_density = self.density_branch(x_density)  
      
        merged = torch.cat([out_rdf, out_density], dim=1)
        
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

    model = RDFDensityMLP()
    model.load_state_dict(torch.load(model_path, map_location=torch.device('cpu'), weights_only=True))
    model.eval()

    
    try:
        input_data = sys.stdin.read()
        input_json = json.loads(input_data)
        rdf = input_json['rdf']          
        density = input_json['density']  
    except Exception as e:
        print(f"Error reading input data: {e}", file=sys.stderr)
        sys.exit(1)

    if len(rdf) != 1000:
        print(f"Error: RDF histogram must have 1000 bins, but got {len(rdf)}.", file=sys.stderr)
        sys.exit(1)
    if not isinstance(density, float) and not isinstance(density, int):
        print(f"Error: Density must be a float, but got {type(density)}.", file=sys.stderr)
        sys.exit(1)

    rdf_tensor = torch.tensor([rdf], dtype=torch.float32)      
    density_tensor = torch.tensor([[density]], dtype=torch.float32)  

    start_time = time.perf_counter()
    with torch.no_grad():
        predicted_energy = model(rdf_tensor, density_tensor).item()
    inference_time = time.perf_counter() - start_time

    print(f"{predicted_energy} {inference_time}")

if __name__ == "__main__":
    main()

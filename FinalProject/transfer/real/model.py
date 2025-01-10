# export_rdf_density_model.py

import torch
import numpy as np
import os
import torch.nn as nn

class RDFDensityMLP(nn.Module):
    """
    Two-branch MLP:
      - Branch 1: RDF (e.g. 1000 bins)
      - Branch 2: single density
    We merge them and output energy.
    """

    def __init__(self, rdf_size=1000, hidden_dim=128):
        super().__init__()
        # 1) RDF branch
        self.rdf_branch = nn.Sequential(
            nn.Linear(rdf_size, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU()
        )
        # 2) Density branch
        #   (just a small sub-network; we can keep it smaller or bigger as you wish)
        self.density_branch = nn.Sequential(
            nn.Linear(1, hidden_dim // 2),
            nn.ReLU(),
            nn.Linear(hidden_dim // 2, hidden_dim // 2),
            nn.ReLU()
        )
        # 3) Merge the two branches
        #   RDF branch outputs hidden_dim, density branch outputs hidden_dim//2
        #   so total is hidden_dim + hidden_dim//2
        merged_dim = hidden_dim + (hidden_dim // 2)
        self.merged_head = nn.Sequential(
            nn.Linear(merged_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, 1)  # final energy prediction
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

def save_weights_to_csv(state_dict, layer_name, folder_path):
    """
    Saves weights and biases of a given layer to CSV files within the specified folder.
    """
    weight = state_dict[f"{layer_name}.weight"].cpu().numpy()
    bias = state_dict[f"{layer_name}.bias"].cpu().numpy()

    weight_path = os.path.join(folder_path, f"{layer_name}_weight.csv")
    bias_path = os.path.join(folder_path, f"{layer_name}_bias.csv")

    np.savetxt(weight_path, weight, delimiter=",")
    np.savetxt(bias_path, bias, delimiter=",")
    print(f"Saved {weight_path} and {bias_path}")

def main():
    # Define model dimensions
    rdf_size = 1000
    hidden_dim = 128

    # Initialize the model and load weights
    model = RDFDensityMLP(rdf_size=rdf_size, hidden_dim=hidden_dim)
    model.load_state_dict(torch.load("best_model.pth", map_location=torch.device('cpu')))
    model.eval()

    # Print model summary (optional)
    from torchsummary import summary
    print("===== PyTorch Model Summary =====")
    summary(model, [(rdf_size,), (1,)])

    # Create 'layers' directory if it doesn't exist
    layers_dir = "layers"
    os.makedirs(layers_dir, exist_ok=True)
    print(f"All CSV files will be saved to the '{layers_dir}' directory.\n")

    # Export all weights and biases
    state_dict = model.state_dict()

    # Define layer names corresponding to the model's layers
    layers = [
        "rdf_branch.0",  # First Linear layer in RDF branch
        "rdf_branch.2",  # Second Linear layer in RDF branch
        "density_branch.0",  # First Linear layer in Density branch
        "density_branch.2",  # Second Linear layer in Density branch
        "merged_head.0",  # First Linear layer in Merged head
        "merged_head.2"   # Second Linear layer in Merged head
    ]

    for layer in layers:
        save_weights_to_csv(state_dict, layer, layers_dir)

    # Create example input
    x_rdf = torch.tensor([[0.5]*rdf_size], dtype=torch.float32)  # shape: (1, 1000)
    x_density = torch.tensor([[1.0]], dtype=torch.float32)       # shape: (1, 1)

    # Forward pass
    with torch.no_grad():
        output = model(x_rdf, x_density)
    print("\nPyTorch forward pass output for x_rdf=[0.5,...,0.5], x_density=[1.0]:", output.item())

    print("\nAll weights, biases exported, and example forward pass completed successfully.")

if __name__ == "__main__":
    main()

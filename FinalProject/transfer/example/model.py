#!/usr/bin/env python3

"""
export_model.py

- Defines a simple MLP in PyTorch
- Uses torchsummary to print a summary
- Exports weights to CSV
- Prints an example forward pass output
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torchsummary import summary  # pip install torchsummary
import numpy as np

# 1) Define a simple MLP class in PyTorch
class SimpleMLP(nn.Module):
    def __init__(self, in_dim, hidden_dim, out_dim):
        super(SimpleMLP, self).__init__()
        self.fc1 = nn.Linear(in_dim, hidden_dim)
        self.fc2 = nn.Linear(hidden_dim, out_dim)
    
    def forward(self, x):
        # x: shape (batch_size, in_dim)
        x = self.fc1(x)       # linear
        x = F.relu(x)         # ReLU
        x = self.fc2(x)       # linear
        return x

if __name__ == "__main__":
    # 2) Instantiate the model with fixed dimensions
    in_dim = 4       # Input size
    hidden_dim = 8   # Hidden layer size
    out_dim = 2      # Output size
    model = SimpleMLP(in_dim, hidden_dim, out_dim)

    # 3) Print a summary
    #    We can feed in a "dummy" input shape (batch_size=1, input_dim=4).
    #    torchsummary expects (channels, height, width) for CNNs, but for an MLP,
    #    we can treat channels=1, height=1, width=in_dim, or just do input_dim=4 as "width".
    #    Another approach is to pass (in_dim,) and set "batch_size=1" in summary(...) call.
    print("===== PyTorch Model Summary =====")
    summary(model, (in_dim,))  # (4,) means just 4 features (like a 1D input)

    # 4) Create some example input
    x = torch.tensor([[0.5, -1.0, 2.0, 0.0]], dtype=torch.float32)  # shape: (1,4)

    # 5) Forward pass in PyTorch
    with torch.no_grad():
        pytorch_output = model(x)
    print("\nPyTorch forward pass output for x=[0.5, -1.0, 2.0, 0.0]:", pytorch_output)

    state_dict = model.state_dict()

    np.savetxt("fc1_weight.csv", state_dict["fc1.weight"].numpy(), delimiter=",")
    np.savetxt("fc1_bias.csv",   state_dict["fc1.bias"].numpy()[None], delimiter=",")

    np.savetxt("fc2_weight.csv", state_dict["fc2.weight"].numpy(), delimiter=",")
    np.savetxt("fc2_bias.csv",   state_dict["fc2.bias"].numpy()[None], delimiter=",")

    print("\nExported CSVs: fc1_weight.csv, fc1_bias.csv, fc2_weight.csv, fc2_bias.csv.")
    print("Done.")

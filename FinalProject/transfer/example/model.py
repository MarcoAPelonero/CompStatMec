import torch
import torch.nn as nn
import torch.nn.functional as F
from torchsummary import summary  
import numpy as np

class SimpleMLP(nn.Module):
    def __init__(self, in_dim, hidden_dim, out_dim):
        super(SimpleMLP, self).__init__()
        self.fc1 = nn.Linear(in_dim, hidden_dim)
        self.fc2 = nn.Linear(hidden_dim, out_dim)
    
    def forward(self, x):
        x = self.fc1(x)       
        x = F.relu(x)         
        x = self.fc2(x)       
        return x

if __name__ == "__main__":
    
    in_dim = 4       
    hidden_dim = 8   
    out_dim = 2      
    model = SimpleMLP(in_dim, hidden_dim, out_dim)

    #    We can feed in a "dummy" input shape (batch_size=1, input_dim=4).
    #    torchsummary expects (channels, height, width) for CNNs, but for an MLP,
    #    we can treat channels=1, height=1, width=in_dim, or just do input_dim=4 as "width".
    #    Another approach is to pass (in_dim,) and set "batch_size=1" in summary(...) call.
    print("===== PyTorch Model Summary =====")
    summary(model, (in_dim,))  # (4,) means just 4 features (like a 1D input)

    x = torch.tensor([[0.5, -1.0, 2.0, 0.0]], dtype=torch.float32)  

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

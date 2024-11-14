import numpy as np
import matplotlib.pyplot as plt

# Generate time data
t = np.linspace(0, 10, 1000)

# Define a sigmoid function for smooth transition
def sigmoid(x, x0, k):
    return 1 / (1 + np.exp(-k * (x - x0)))

# Generate data: initial fluctuations and then smooth decay to equilibrium
chaotic_part = 2 + 0.05 * np.random.randn(len(t)) + 2 * np.exp(-0.5 * t)
equilibrium_part = 0.5 + 0.05 * np.random.randn(len(t)) + 0.5 * np.exp(-0.5 * (t - 3))
transition = sigmoid(t, 3, 10)  # Smooth transition around t=3

y = chaotic_part * (1 - transition) + equilibrium_part * transition

# Plot the graph
plt.figure(figsize=(8, 5))
plt.plot(t, y, color='black')
plt.axhline(y=0.5, color='black', linestyle='--', label=r'$U_{ex}$')
plt.axvline(x=3, color='green', linestyle='--')
plt.text(3.1, 1.0, r'$t_{eq}$' '\n' r'(Equilibration' '\n' r'"Time")', color='green', fontsize=10, ha='left')

# Labels and customizations
plt.xlabel(r'$t$', fontsize=14, color='green')
plt.ylabel(r'$Y_t$', fontsize=14, color='green')
plt.xticks(color='green')
plt.yticks(color='green')
plt.title(r'$Y_t$ over time with Equilibration Time $t_{eq}$', fontsize=14, color='green')
plt.grid(True)

# Show plot
plt.show()

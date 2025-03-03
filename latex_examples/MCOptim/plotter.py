#!/usr/bin/env python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Read results CSV
df = pd.read_csv("results.csv")

# Create a seaborn plot
sns.set(style="whitegrid")
plt.figure(figsize=(8,6))
ax = sns.lineplot(data=df, x="N", y="AvgTimePerStep", hue="Mode", marker="o")

ax.set_title("Average Time per Monte Carlo Step")
ax.set_xlabel("Number of Particles (N)")
ax.set_ylabel("Avg Time per Step (s)")

plt.tight_layout()
plt.savefig("simulation_timing.png")
plt.show()

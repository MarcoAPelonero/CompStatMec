import numpy as np
import matplotlib.pyplot as plt

with open("nile_depths.txt", "r") as file:
    actual_avg_depth = float(file.readline().split(",")[1])

depth_data = np.loadtxt("nile_depths.txt", delimiter=",", skiprows=1)
mc_depths = depth_data[:, 0]
mcm_depths = depth_data[:, 1]

steps = np.arange(1, len(mc_depths) + 1)

plt.figure(figsize=(10, 5))
plt.plot(steps, mc_depths, label='Monte Carlo Depth')
plt.plot(steps, mcm_depths, label='Metropolis Depth')
plt.xlabel("Steps")
plt.ylabel("Average Depth")
plt.title("Comparison of Average Depth Estimations Over Time")
plt.legend()
plt.grid()
plt.savefig("depth_comparison.png")
plt.show()

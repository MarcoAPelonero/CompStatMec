import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import os
import shutil

#timesteps = [0.001, 0.003, 0.005, 0.008, 0.01]
T = 1
timesteps = [0.001,0.002,0.003,0.004,0.005,0.006,0.007]
avg_std_standard = []
avg_std_shadow = []

total_runs = len(timesteps) * 3
progress_bar = tqdm(total=total_runs, desc="Simulations")

for dt in timesteps:
    stds_standard = []
    stds_shadow = []
    for run in range(3):
        if os.path.exists("energy.dat"):
            os.remove("energy.dat")
        steps = int(T / dt)
        result = subprocess.run(["./simulation.exe", str(dt), str(steps)], capture_output=True)
        if result.returncode != 0:
            print(f"Simulation failed at dt={dt} on run {run+1}.")
            continue

        try:
            df = pd.read_csv("energy.dat", comment='#', sep='\s+',
                             names=["time", "H_standard", "H_shadow"])
        except Exception as e:
            print(f"Error reading energy.dat: {e}")
            continue
        
        df = df[int(len(df) * 0.2):]
        mean_standard = df["H_standard"].mean()
        mean_shadow = df["H_shadow"].mean()

        regularized_standard = (df["H_standard"] - mean_standard)/mean_standard
        regularized_shadow = (df["H_shadow"] - mean_shadow)/mean_shadow 
        std_standard = regularized_standard.std()
        std_shadow = regularized_shadow.std()
        stds_standard.append(std_standard)
        stds_shadow.append(std_shadow)

        progress_bar.update(1)
    
    if stds_standard and stds_shadow:
        avg_std_standard.append(np.mean(stds_standard))
        avg_std_shadow.append(np.mean(stds_shadow))
    else:
        avg_std_standard.append(np.nan)
        avg_std_shadow.append(np.nan)

progress_bar.close()

sns.set(style="whitegrid")
plt.figure(figsize=(10, 6))
sns.lineplot(x=timesteps, y=avg_std_standard, marker='o', label="Standard Hamiltonian")
sns.lineplot(x=timesteps, y=avg_std_shadow, marker='o', label="Shadow Hamiltonian")
plt.xlabel("Timestep")
plt.ylabel("Average Standard Deviation")
plt.title("Averaged Standard Deviation vs Timestep")
plt.legend()
plt.tight_layout()
plt.savefig("std_plot.png")
plt.show()
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv("energy.dat", comment='#', delim_whitespace=True,
                 names=["time", "H_standard", "H_shadow"])

# Cut the first 20% of the data
df = df[int(len(df) * 0.2):]

sns.set(style="whitegrid")
plt.figure(figsize=(10,6))
sns.lineplot(x="time", y="H_standard", data=df, label="Standard Hamiltonian")
sns.lineplot(x="time", y="H_shadow", data=df, label="Shadow Hamiltonian")
plt.xlabel("Time")
plt.ylabel("Total Energy")
plt.legend()
plt.tight_layout()
plt.savefig("energy_plot.png")
plt.show()
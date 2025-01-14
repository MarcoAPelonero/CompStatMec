import re
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

sns.set_theme(style="darkgrid")
data_dir = Path("mlcomparisons")
files = sorted(data_dir.glob("mlComparison_rho_*_T_*.dat"))
if len(files) != 4:
    raise ValueError("Expected exactly 4 simulations in 'mlcomparisons/'")

cols = sns.color_palette("Purples", 3)
cols = ["#6A0DAD", "#333333", "#00008B"]
pat = re.compile(r"mlComparison_rho_(?P<rho>[\d\.e-]+)_T_(?P<T>[\d\.e-]+)\.dat")
all_data = []

for f in files:
    m = pat.match(f.name)
    arr = np.loadtxt(f) if m else None
    if arr is None or arr.ndim not in [1, 2]: 
        continue
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    if arr.shape[1] != 6:
        continue
    all_data.append({
        "rho": m.group("rho"),
        "T": m.group("T"),
        "t_ana": arr[:, 0],
        "E_ana": arr[:, 1],
        "t_cpp": arr[:, 2],
        "E_cpp": arr[:, 3],
        "t_py":  arr[:, 4],
        "E_py":  arr[:, 5],
    })

tol = 1e-6
ok = all(np.allclose(d["E_cpp"], d["E_py"], atol=tol) for d in all_data)
print("ML C++ and ML Python energies are within tolerance." if ok 
      else "ML C++ and ML Python differ beyond tolerance.")

fig_t, axs_t = plt.subplots(2, 2, figsize=(12, 8))
for i, d in enumerate(all_data):
    ax = axs_t.flat[i]
    steps = range(len(d["t_ana"]))
    ax.plot(steps, d["t_ana"], c=cols[0], label="$t_{ana}$", lw=2)
    ax.plot(steps, d["t_cpp"], c=cols[1], label="$t_{cpp}$", lw=2, ls="--")
    ax.plot(steps, d["t_py"],  c=cols[2], label="$t_{py}$",  lw=2, ls=":")
    ax.set_title(f"Time (rho={d['rho']}, T={d['T']})")
    ax.set_xlabel("Step")
    ax.set_yscale('log')
    ax.set_ylabel("Time [s]")
handles, labels = axs_t[0, 0].get_legend_handles_labels()
fig_t.legend(handles, labels, loc="upper center", ncol=3)
fig_t.suptitle("Time Evolution (4 Simulations)", y=1.03, fontsize=14)
fig_t.tight_layout()
plt.savefig("Time_Evolution_4Simulations.png", dpi=150)

fig_e, axs_e = plt.subplots(2, 2, figsize=(12, 8))
for i, d in enumerate(all_data):
    ax = axs_e.flat[i]
    st = range(len(d["E_ana"]))
    ax.plot(st, d["E_ana"], c=cols[0], label="E_ana", lw=2)
    ax.plot(st, d["E_cpp"], c=cols[1], label="E_cpp", lw=2, ls="--")
    ax.plot(st, d["E_py"],  c=cols[2], label="E_py",  lw=2, ls=":")
    ax.set_title(f"Energy (rho={d['rho']}, T={d['T']})")
    ax.set_xlabel("Step")
    ax.set_ylabel("Energy")
handles, labels = axs_e[0, 0].get_legend_handles_labels()
fig_e.legend(handles, labels, loc="upper center", ncol=3)
fig_e.suptitle("Energy Evolution (4 Simulations)", y=1.03, fontsize=14)
fig_e.tight_layout()
plt.savefig("Energy_Evolution_4Simulations.png", dpi=150)
plt.show()

mE_ana = np.mean([np.mean(d["E_ana"]) for d in all_data])
mE_cpp = np.mean([np.mean(d["E_cpp"]) for d in all_data])
mE_py  = np.mean([np.mean(d["E_py"])  for d in all_data])
mT_ana = np.mean([np.mean(d["t_ana"]) for d in all_data])
mT_cpp = np.mean([np.mean(d["t_cpp"]) for d in all_data])
mT_py  = np.mean([np.mean(d["t_py"])  for d in all_data])

print(rf"""
\begin{{table}}[h!]
\centering
\caption{{Comparison of average results (4 simulations)}}
\begin{{tabular}}{{lcc}}
\toprule
Method & Avg. Energy & Avg. Time (s) \\
\midrule
Analytical & {mE_ana:.6f} & {mT_ana:.6f} \\
ML C++     & {mE_cpp:.6f} & {mT_cpp:.6f} \\
ML Python  & {mE_py:.6f} & {mT_py:.6f} \\
\bottomrule
\end{{tabular}}
\end{{table}}
""")

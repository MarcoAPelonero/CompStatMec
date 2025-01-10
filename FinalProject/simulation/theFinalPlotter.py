import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D  # For custom legend entries

# -----------------------------
# 1) CONFIGURE FOLDERS & STYLING
# -----------------------------
sns.set_theme(style="darkgrid")  # You can choose 'whitegrid', 'darkgrid', 'ticks', etc.

data_folder = Path("mlcomparisons")  # Folder containing the .dat files
output_folder = Path("finalImages")  # Folder to save the plots
output_folder.mkdir(parents=True, exist_ok=True)  # Create the folder if it doesn't exist

# Regex pattern to extract rho and T from filenames
file_pattern = re.compile(
    r"mlComparison_rho_(?P<rho>[\d\.e-]+)_T_(?P<T>[\d\.e-]+)\.dat"
)

# -----------------------------
# 2) DEFINE COLOR PALETTES
# -----------------------------

# Color palette for simulations (34 distinct colors)
# If you have exactly 34 simulations, adjust n_colors accordingly.
# Seaborn's "tab20" has 20 distinct colors; to get 34, we'll cycle through multiple palettes.
simulation_palette = sns.color_palette("tab20", n_colors=20) + sns.color_palette("hsv", n_colors=14)

# Ensure the palette has exactly 34 colors
simulation_palette = simulation_palette[:34]

# Color palette for computation methods
# Assign fixed colors to ensure consistency across all plots
method_colors = {
    "Analytical": "purple",
    "ML C++": "pink",
    "ML Python": "green"
}

# -----------------------------
# 3) LOAD ALL SIMULATIONS
# -----------------------------
all_sims = []
for fpath in data_folder.glob("mlComparison_rho_*_T_*.dat"):
    match = file_pattern.match(fpath.name)
    if not match:
        print(f"Filename '{fpath.name}' does not match the pattern. Skipped.")
        continue

    rho = match.group("rho")
    T = match.group("T")

    try:
        # Load data: (time_ana, E_ana, time_cpp, E_cpp, time_py, E_py)
        data = np.loadtxt(fpath)
        if data.ndim == 1:
            data = data.reshape(1, -1)
    except Exception as e:
        print(f"Error loading '{fpath.name}': {e}. Skipped.")
        continue

    if data.shape[1] != 6:
        print(f"File '{fpath.name}' does not have 6 columns. Skipped.")
        continue

    # Extract columns
    t_ana, E_ana = data[:, 0], data[:, 1]
    t_cpp, E_cpp = data[:, 2], data[:, 3]
    t_py, E_py = data[:, 4], data[:, 5]

    sim_label = f"Sim (rho={rho}, T={T})"

    all_sims.append({
        "filename": fpath.name,
        "rho": rho,
        "T": T,
        "t_ana": t_ana,
        "E_ana": E_ana,
        "t_cpp": t_cpp,
        "E_cpp": E_cpp,
        "t_py":  t_py,
        "E_py":  E_py,
        "label": sim_label
    })

# Sort simulations alphabetically by label
all_sims.sort(key=lambda s: s["label"])

# Check if the number of simulations exceeds the palette size
num_sims = len(all_sims)
if num_sims > len(simulation_palette):
    print(f"Warning: Number of simulations ({num_sims}) exceeds the number of unique colors in the simulation palette ({len(simulation_palette)}). Colors will repeat.")

# -----------------------------
# 4) GENERATE INDIVIDUAL PLOTS
# -----------------------------
for idx, sim in enumerate(all_sims):
    # Assign a unique color to the simulation for potential future use
    # (Not used directly in individual plots since methods have fixed colors)
    simulation_color = simulation_palette[idx % len(simulation_palette)]

    # Extract data
    t_ana = sim["t_ana"]
    t_cpp = sim["t_cpp"]
    t_py  = sim["t_py"]

    # Create a figure with one subplot
    fig, ax = plt.subplots(figsize=(12, 6))

    # Title of the figure
    fig.suptitle(f"Simulation (rho={sim['rho']}, T={sim['T']})", fontsize=16)

    # --- Subplot: Computation Time ---
    ax.set_title("Computation Time Evolution")
    ax.set_xlabel("Step")
    ax.set_ylabel("Time [s]")
    ax.set_yscale('log')  # Set y-axis to log scale

    # Plot computation times for the 3 methods with fixed colors
    ax.plot(np.arange(len(t_ana)), t_ana, label="Analytical", color=method_colors["Analytical"], linestyle="-", linewidth=2)
    ax.plot(np.arange(len(t_cpp)), t_cpp, label="ML C++", color=method_colors["ML C++"], linestyle="--", linewidth=2)
    ax.plot(np.arange(len(t_py)),  t_py,  label="ML Python", color=method_colors["ML Python"], linestyle=":", linewidth=2)

    # Add legend for methods
    ax.legend()

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust layout to make room for the main title

    # Save the figure
    output_filename = f"Time_Evolution_rho_{sim['rho']}_T_{sim['T']}.png"
    output_path = output_folder / output_filename
    plt.savefig(output_path, dpi=150)
    plt.close(fig)

    print(f"Saved individual plot: {output_path}")

# -----------------------------
# 5) GENERATE COMBINED PLOTS
# -----------------------------
# These plots show all simulations in a single plot for comparison

# -----------------------------
# 5a) Combined Energy Evolution Plot
# -----------------------------
fig_energy, ax_energy = plt.subplots(figsize=(14, 8))
fig_energy.suptitle("Comparison of Energy Evolution Across All Simulations", fontsize=18)
ax_energy.set_xlabel("Step")
ax_energy.set_ylabel("Energy per Particle")

for idx, sim in enumerate(all_sims):
    color_i = simulation_palette[idx % len(simulation_palette)]

    # Energy
    E_ana = sim["E_ana"]
    E_cpp = sim["E_cpp"]
    E_py  = sim["E_py"]
    E_min = np.minimum.reduce([E_ana, E_cpp, E_py]) - 0.3
    E_max = np.maximum.reduce([E_ana, E_cpp, E_py]) + 0.3

    # Fill between for energy range
    ax_energy.fill_between(
        np.arange(len(E_ana)),
        E_min,
        E_max,
        color=color_i,
        alpha=0.2
    )
    # Plot the energy curves without adding to the legend
    ax_energy.plot(np.arange(len(E_ana)), E_ana, color=color_i, linestyle="-", linewidth=1, alpha=0.6)
    ax_energy.plot(np.arange(len(E_cpp)), E_cpp, color=color_i, linestyle="--", linewidth=1, alpha=0.6)
    ax_energy.plot(np.arange(len(E_py)),  E_py,  color=color_i, linestyle=":", linewidth=1, alpha=0.6)

# Create custom legend for simulations using colored patches
simulation_patches = [Line2D([0], [0], color=simulation_palette[idx % len(simulation_palette)], lw=4, label=sim["label"]) 
                      for idx, sim in enumerate(all_sims)]

# Due to the large number of simulations, it's practical to place the legend outside the plot
ax_energy.legend(handles=simulation_patches, title="Simulations", bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout(rect=[0, 0, 0.85, 0.95])  # Adjust layout to make room for the legend

# Save the combined energy plot
combined_energy_path = output_folder / "Combined_Energy_Evolution_AllSimulations.png"
fig_energy.savefig(combined_energy_path, dpi=150, bbox_inches='tight')
plt.close(fig_energy)

print(f"Saved combined energy plot: {combined_energy_path}")

# -----------------------------
# 5b) Combined Computation Time Evolution Plot
# -----------------------------
fig_time, ax_time = plt.subplots(figsize=(14, 8))
fig_time.suptitle("Comparison of Computation Time Evolution Across All Simulations", fontsize=18)
ax_time.set_xlabel("Step")
ax_time.set_ylabel("Time [s]")
ax_time.set_yscale('log')  # Log scale for better visualization

for idx, sim in enumerate(all_sims):
    color_i = simulation_palette[idx % len(simulation_palette)]

    # Computation Time
    t_ana = sim["t_ana"]
    t_cpp = sim["t_cpp"]
    t_py  = sim["t_py"]

    # Plot computation times for the 3 methods with fixed colors
    ax_time.plot(np.arange(len(t_ana)), t_ana, label="Analytical" if idx == 0 else None, color=method_colors["Analytical"], linestyle="-", linewidth=1, alpha=0.6)
    ax_time.plot(np.arange(len(t_cpp)), t_cpp, label="ML C++" if idx == 0 else None, color=method_colors["ML C++"], linestyle="--", linewidth=1, alpha=0.6)
    ax_time.plot(np.arange(len(t_py)),  t_py,  label="ML Python" if idx == 0 else None, color=method_colors["ML Python"], linestyle=":", linewidth=1, alpha=0.6)

# Create custom legend for computation methods
method_lines = [
    Line2D([0], [0], color=method_colors["Analytical"], linestyle='-', lw=2, label='Analytical'),
    Line2D([0], [0], color=method_colors["ML C++"], linestyle='--', lw=2, label='ML C++'),
    Line2D([0], [0], color=method_colors["ML Python"], linestyle=':', lw=2, label='ML Python')
]

# Add the legend for methods to the plot
ax_time.legend(handles=method_lines, title="Computation Methods", loc='upper right')

plt.tight_layout()

# Save the combined computation time plot
combined_time_path = output_folder / "Combined_Time_Evolution_AllSimulations.png"
fig_time.savefig(combined_time_path, dpi=150)
plt.close(fig_time)

print(f"Saved combined computation time plot: {combined_time_path}")

print("\nAll plots have been generated and saved in the 'finalImages/' folder.")

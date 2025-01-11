import subprocess
import time
from pathlib import Path
import statistics

# Define simulation parameters
methods = ["RDFSpeedVerlet"]
step_counts = [1000, 5000, 10000]
N = 100
dt = 0.005

# Path to the simulation executable
simulation_executable = Path('./build/main.exe')

# Dictionary to store execution times
execution_times = {method: {steps: [] for steps in step_counts} for method in methods}

# Run simulations and record times
for method in methods:
    for steps in step_counts:
        print(f"Running {method} with {steps} steps...")
        times = []
        for run in range(5):
            print(f"Run {run+1}/5...")
            start_time = time.time()
            
            # Build command arguments
            cmd = [
                str(simulation_executable),
                str(N),
                "0.6",           # rho
                str(dt),
                str(steps),
                "trajectory.dat",
                "radialDistribution.dat",
                method,
                "1.0",           # T
                "1.0",           # p
                "1.0"            # taup
            ]
            
            # Execute the simulation
            subprocess.run(cmd, check=True)
            
            end_time = time.time()
            elapsed_time = end_time - start_time
            times.append(elapsed_time)
            print(f"Completed in {elapsed_time:.2f} seconds.\n")
        
        average_time = statistics.mean(times)
        sigma_time = statistics.stdev(times)
        execution_times[method][steps] = (average_time, sigma_time)
        print(f"{method} with {steps} steps: Avg={average_time:.2f}s, σ={sigma_time:.2f}s\n")

# Generate LaTeX table
latex_table = """
\\begin{tabular}{lccc}
\\hline
Method & 1000 Steps & 5000 Steps & 10000 Steps \\\\
\\hline
"""

for method in methods:
    row = f"{method}"
    for steps in step_counts:
        avg, sigma = execution_times[method][steps]
        row += f" & ${avg:.2f} \\pm {sigma:.2f}$"
    row += " \\\\\n"
    latex_table += row

latex_table += "\\hline\n\\end{tabular}"

print("LaTeX Table of Execution Times:")
print(latex_table)
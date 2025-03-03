#!/usr/bin/env python3
import subprocess
import csv

# List of simulation modes with a label (0: Full, 1: Truncated, 2: Linked Cells, 3: Neighbor List)
modes = {
    0: "Full",
    1: "Truncated",
    2: "Linked Cells",
    3: "Neighbor List"
}

# List of particle numbers to test
N_list = [100, 500, 1000,2000,3000,5000,10000,15000,20000,30000]
steps = 50000  # Number of MC steps for each run

# Path to the compiled simulation executable
executable = "./sim"

results = []  # Will store rows: mode_label, N, avg_time

for mode in modes:
    for N in N_list:
        print(f"Running mode {modes[mode]} with N={N}...")
        # Run the simulation
        proc = subprocess.run([executable, str(mode), str(N), str(steps)], 
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        # The simulation prints: mode,N,steps,accepted,avg_time_per_step
        if proc.returncode != 0:
            print(f"Error running simulation: {proc.stderr}")
            continue
        try:
            output = proc.stdout.strip()
            parts = output.split(',')
            avg_time = float(parts[4])
            results.append([modes[mode], N, avg_time])
            print(f"  Avg time per MC step: {avg_time:.6e} s")
        except Exception as e:
            print("Error parsing output:", e)

# Write results to CSV
csv_file = "results.csv"
with open(csv_file, "w", newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["Mode", "N", "AvgTimePerStep"])
    writer.writerows(results)

print(f"Results written to {csv_file}")

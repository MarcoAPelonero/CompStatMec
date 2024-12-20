import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from tqdm import tqdm

# Input file
rdf_file = "radialDistribution.dat"

def read_rdf_file(file_path):
    steps = []
    rdfs = []
    r_values = None

    with open(file_path, 'r') as f:
        current_step_data = []
        current_step = None
        for line in tqdm(f, desc='Reading RDF file'):
            line = line.strip()
            if line.startswith('# Step'):
                if current_step_data and current_step is not None:
                    rdfs.append([g for r, g in current_step_data])
                    steps.append(current_step)
                current_step = int(line.split()[2])
                current_step_data = []
            else:
                if line:
                    vals = line.split()
                    r = float(vals[0])
                    g = float(vals[1])
                    current_step_data.append((r, g))

        if current_step_data and current_step is not None:
            rdfs.append([g for r, g in current_step_data])
            steps.append(current_step)

    r_values = [d[0] for d in current_step_data]
    rdfs = np.array(rdfs)
    return r_values, rdfs, steps

# Compute stabilization index and averages
def compute_stabilized_rdf(rdfs, threshold_factor=0.1):
    fluctuations = np.std(rdfs, axis=0)
    threshold = np.mean(fluctuations) * threshold_factor
    stabilization_index = None
    for i in tqdm(range(len(rdfs)), desc='Finding stabilization index', total=len(rdfs)):
        if np.mean(np.std(rdfs[i:], axis=0)) < threshold:
            stabilization_index = i
            break
    if stabilization_index is None:
        stabilization_index = len(rdfs) // 2
    stabilized_rdf = np.mean(rdfs[stabilization_index:], axis=0)
    return stabilization_index, stabilized_rdf

# Plotting function
def plot_rdf(r_values, rdfs, steps, stabilized_rdf, final_rdf, stabilization_index):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    plt.style.use('seaborn-darkgrid')

    # Animation axis
    line1, = ax1.plot([], [], lw=2, color='royalblue')
    ax1.set_xlim(0, max(r_values))
    ax1.set_ylim(0, max(final_rdf) * 1.2)
    ax1.set_xlabel('Distance $r$', fontsize=12)
    ax1.set_ylabel('Radial Distribution $g(r)$', fontsize=12)
    ax1.set_title('RDF Evolution Over Time', fontsize=14)
    title = ax1.text(0.5, 1.05, '', transform=ax1.transAxes, ha='center', fontsize=12)

    # Stabilized and final RDFs
    ax2.plot(r_values, final_rdf, label='Final RDF', color='firebrick', lw=2)
    ax2.plot(r_values, stabilized_rdf, label='Stabilized RDF', color='darkgreen', lw=2)
    ax2.axvline(r_values[stabilization_index], color='grey', linestyle='--', label='Stabilization Start')
    ax2.set_xlim(0, max(r_values))
    ax2.set_ylim(0, max(final_rdf) * 1.2)
    ax2.set_xlabel('Distance $r$', fontsize=12)
    ax2.set_ylabel('Radial Distribution $g(r)$', fontsize=12)
    ax2.set_title('Final and Stabilized RDFs', fontsize=14)
    ax2.legend(fontsize=10)

    # Animation initialization
    def init():
        line1.set_data([], [])
        title.set_text('')
        return line1, title

    # Animation update
    def animate(i):
        g = rdfs[i]
        line1.set_data(r_values, g)
        title.set_text(f'Step: {steps[i]}')
        return line1, title

    ani = animation.FuncAnimation(fig, animate, frames=len(steps), init_func=init, interval=100, blit=True)

    plt.tight_layout()
    plt.show()

# Main script
r_values, rdfs, steps = read_rdf_file(rdf_file)
stabilization_index, stabilized_rdf = compute_stabilized_rdf(rdfs)
final_rdf = rdfs[-1]

plot_rdf(r_values, rdfs, steps, stabilized_rdf, final_rdf, stabilization_index)

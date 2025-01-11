import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import seaborn as sns
import sys

class TqdmPillowWriter(animation.PillowWriter):
    def __init__(self, fps=10, codec='gif', bitrate=None, metadata=None, total_frames=0):
        super().__init__(fps=fps, codec=codec, bitrate=bitrate, metadata=metadata)
        self.pbar = None
        self.total_frames = total_frames

    def setup(self, fig, outfile, dpi=None, **kwargs):
        super().setup(fig, outfile, dpi=dpi, **kwargs)
        if self.total_frames > 0:
            self.pbar = tqdm(total=self.total_frames, desc='Saving Animation')

    def grab_frame(self, **savefig_kwargs):
        super().grab_frame(**savefig_kwargs)
        if self.pbar:
            self.pbar.update(1)

    def finish(self):
        super().finish()
        if self.pbar:
            self.pbar.close()

sns.set_style("whitegrid")
sns.set_palette("viridis")

rdf_file = "radialDistribution.dat"

def read_rdf_file(file_path):
    steps = []
    rdfs = []
    temperatures = []
    pressures = []
    pepp = []
    potential_energies = []
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
                parts = line.split()
                current_step = int(parts[2])
                try:
                    temp_line = next(f).strip()
                except StopIteration:
                    break
                temp_parts = temp_line.split()
                if len(temp_parts) < 4:
                    break
                temperature = float(temp_parts[0])
                pressure = float(temp_parts[3])
                potential_energy_per_particle = float(temp_parts[1])
                potential_energy = float(temp_parts[2])
                temperatures.append(temperature)
                pressures.append(pressure)
                pepp.append(potential_energy_per_particle)
                potential_energies.append(potential_energy)
                current_step_data = []
            else:
                if line:
                    vals = line.split()
                    if len(vals) < 2:
                        continue
                    try:
                        r = float(vals[0])
                        g = float(vals[1])
                        current_step_data.append((r, g))
                    except ValueError:
                        continue

        if current_step_data and current_step is not None:
            rdfs.append([g for r, g in current_step_data])
            steps.append(current_step)

    if rdfs:
        r_values = [r for r, g in current_step_data]
    else:
        r_values = []

    rdfs = np.array(rdfs)
    return r_values, rdfs, steps, temperatures, pressures, pepp, potential_energies

def plot_rdf(r_values, rdfs, steps, final_rdf, temperatures, pressures, pepp, potential_energies, skipsteps=1):
    fig, axes = plt.subplots(3, 2, figsize=(18, 12))
    start_index = int(0.2 * len(steps))
    
    avg_pressure = np.mean(pressures[start_index:])
    avg_pepp = np.mean(pepp[start_index:])
    avg_potential_energy = np.mean(potential_energies[start_index:])

    ax1 = axes[0, 0]
    line1, = ax1.plot([], [], lw=2)
    ax1.set_xlim(0, max(r_values) if r_values else 1)
    ax1.set_ylim(0, max(final_rdf) * 1.2 if final_rdf.size > 0 else 1)
    ax1.set_xlabel('Distance $r$', fontsize=14)
    ax1.set_ylabel('Radial Distribution $g(r)$', fontsize=14)
    ax1.set_title('RDF Evolution Over Time', fontsize=16)
    title = ax1.text(0.5, 1.05, '', transform=ax1.transAxes, ha='center', fontsize=12)

    ax2 = axes[0, 1]
    ax2.plot(r_values, final_rdf, label='Final RDF', lw=2)
    ax2.set_xlim(0, max(r_values) if r_values else 1)
    ax2.set_ylim(0, max(final_rdf) * 1.2 if final_rdf.size > 0 else 1)
    ax2.set_xlabel('Distance $r$', fontsize=14)
    ax2.set_ylabel('Radial Distribution $g(r)$', fontsize=14)
    ax2.set_title('Final and Stabilized RDFs', fontsize=16)
    ax2.legend(fontsize=12)

    # **Corrected Section for ax3**
    ax3 = axes[1, 0]
    avg_g = np.mean(rdfs, axis=0)
    ax3.plot(r_values, avg_g, color='skyblue', lw=2)
    ax3.set_xlabel('Distance $r$', fontsize=14)
    ax3.set_ylabel('Average Radial Distribution $g(r)$', fontsize=14)
    ax3.set_title('Average Radial Distribution Over Time', fontsize=16)

    ax4 = axes[1, 1]
    ax4.plot(steps, pressures, color='purple', lw=2)
    ax4.set_xlabel('Step', fontsize=14)
    ax4.set_ylabel('Pressure', fontsize=14)
    ax4.set_title(f'Pressure Over Time (Avg after 20%: {avg_pressure:.2f})', fontsize=16)

    ax5 = axes[2, 0]
    ax5.plot(steps, pepp, color='green', lw=2)
    ax5.set_xlabel('Step', fontsize=14)
    ax5.set_ylabel('Potential Energy per Particle', fontsize=14)
    ax5.set_title(f'Potential Energy per Particle Over Time (Avg after 20%: {avg_pepp:.2f})', fontsize=16)

    ax6 = axes[2, 1]
    ax6.plot(steps, potential_energies, color='brown', lw=2)
    ax6.set_xlabel('Step', fontsize=14)
    ax6.set_ylabel('Potential Energy', fontsize=14)
    ax6.set_title(f'Potential Energy Over Time (Avg after 20%: {avg_potential_energy:.2f})', fontsize=16)

    def init():
        line1.set_data([], [])
        title.set_text('')
        return line1, title

    selected_indices = np.arange(0, len(steps), skipsteps)
    selected_steps = [steps[i] for i in selected_indices]
    selected_rdfs = rdfs[selected_indices]

    def animate(i):
        g = selected_rdfs[i]
        line1.set_data(r_values, g)
        title.set_text(f'Step: {selected_steps[i]}')
        return line1, title

    ani = animation.FuncAnimation(
        fig, animate, frames=len(selected_steps), init_func=init, interval=100, blit=True
    )

    plt.tight_layout()
    writer = TqdmPillowWriter(fps=10, total_frames=len(selected_steps))
    plt.show()


if __name__ == "__main__":
    default_skipsteps = 1

    if len(sys.argv) > 1:
        try:
            skipsteps = int(sys.argv[1])
            if skipsteps < 1:
                skipsteps = default_skipsteps
        except ValueError:
            skipsteps = default_skipsteps
    else:
        skipsteps = default_skipsteps

    print(f"Using skipsteps = {skipsteps}")

    r_values, rdfs, steps, temperatures, pressures, pepp, potential_energies = read_rdf_file(rdf_file)
    
    if rdfs.size == 0:
        print("No RDF data found. Please check the input file.")
    else:
        final_rdf = rdfs[-1]
        plot_rdf(
            r_values,
            rdfs,
            steps,
            final_rdf,
            temperatures,
            pressures,
            pepp,
            potential_energies,
            skipsteps=skipsteps
        )

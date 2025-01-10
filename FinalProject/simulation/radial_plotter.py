import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import seaborn as sns  # Added import
import sys

class TqdmPillowWriter(animation.PillowWriter):
    def __init__(self, fps=10, codec='gif', bitrate=None, metadata=None, total_frames=0):
        super().__init__(fps=fps, codec=codec, bitrate=bitrate, metadata=metadata)
        self.pbar = None
        self.total_frames = total_frames  # Store total frames

    def setup(self, fig, outfile, dpi=None, **kwargs):
        super().setup(fig, outfile, dpi=dpi, **kwargs)
        # Initialize tqdm with the total number of frames if total_frames > 0
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

sns.set_style("darkgrid")

# Input file
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
                    # Read the next line for temperature, pressure, potential energy
                    temp_line = next(f).strip()
                except StopIteration:
                    print(f"Unexpected end of file after step {current_step}.")
                    break
                temp_parts = temp_line.split()
                if len(temp_parts) < 4:
                    print(f"Insufficient data in temperature line after step {current_step}.")
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
                        continue  # Skip malformed lines
                    try:
                        r = float(vals[0])
                        g = float(vals[1])
                        current_step_data.append((r, g))
                    except ValueError:
                        continue  # Skip lines with non-numeric data

        # Handle the last step
        if current_step_data and current_step is not None:
            rdfs.append([g for r, g in current_step_data])
            steps.append(current_step)

    # Extract r_values from the first step for consistency
    if rdfs:
        r_values = [r for r, g in current_step_data]  # Assuming r is the same across all steps
    else:
        r_values = []

    rdfs = np.array(rdfs)
    return r_values, rdfs, steps, temperatures, pressures, pepp, potential_energies

# Plotting function with skipsteps
def plot_rdf(r_values, rdfs, steps, final_rdf, temperatures, pressures, pepp, potential_energies, skipsteps=1):
    fig, axes = plt.subplots(3, 2, figsize=(15, 10))  # Increased figure size for better readability

    # Calculate the index after 20% of the simulation
    start_index = int(0.2 * len(steps))
    
    # Calculate averages after 20%
    avg_temperature = np.mean(temperatures[start_index:])
    avg_pressure = np.mean(pressures[start_index:])
    avg_pepp = np.mean(pepp[start_index:])
    avg_potential_energy = np.mean(potential_energies[start_index:])
    
    # RDF Evolution Over Time
    ax1 = axes[0, 0]
    line1, = ax1.plot([], [], lw=2, color='royalblue')
    ax1.set_xlim(0, max(r_values) if r_values else 1)
    ax1.set_ylim(0, max(final_rdf) * 1.2 if final_rdf.size > 0 else 1)
    ax1.set_xlabel('Distance $r$', fontsize=12)
    ax1.set_ylabel('Radial Distribution $g(r)$', fontsize=12)
    ax1.set_title('RDF Evolution Over Time', fontsize=14)
    title = ax1.text(0.5, 1.05, '', transform=ax1.transAxes, ha='center', fontsize=12)
    
    # Final and Stabilized RDFs
    ax2 = axes[0, 1]
    ax2.plot(r_values, final_rdf, label='Final RDF', color='firebrick', lw=2)
    ax2.set_xlim(0, max(r_values) if r_values else 1)
    ax2.set_ylim(0, max(final_rdf) * 1.2 if final_rdf.size > 0 else 1)
    ax2.set_xlabel('Distance $r$', fontsize=12)
    ax2.set_ylabel('Radial Distribution $g(r)$', fontsize=12)
    ax2.set_title('Final and Stabilized RDFs', fontsize=14)
    ax2.legend(fontsize=10)
    
    # Temperature Plot
    ax3 = axes[1, 0]
    ax3.plot(steps, temperatures, color='orange')
    ax3.set_xlabel('Step', fontsize=12)
    ax3.set_ylabel('Temperature', fontsize=12)
    ax3.set_title(f'Temperature Over Time (Avg after 20%: {avg_temperature:.2f})', fontsize=14)
    
    # Pressure Plot
    ax4 = axes[1, 1]
    ax4.plot(steps, pressures, color='purple')
    ax4.set_xlabel('Step', fontsize=12)
    ax4.set_ylabel('Pressure', fontsize=12)
    ax4.set_title(f'Pressure Over Time (Avg after 20%: {avg_pressure:.2f})', fontsize=14)
    
    # Potential Energy per Particle Plot
    ax5 = axes[2, 0]
    ax5.plot(steps, pepp, color='green')
    ax5.set_xlabel('Step', fontsize=12)
    ax5.set_ylabel('Potential Energy per Particle', fontsize=12)
    ax5.set_title(f'Potential Energy per Particle Over Time (Avg after 20%: {avg_pepp:.2f})', fontsize=14)
    
    # Potential Energy Plot
    ax6 = axes[2, 1]
    ax6.plot(steps, potential_energies, color='brown')
    ax6.set_xlabel('Step', fontsize=12)
    ax6.set_ylabel('Potential Energy', fontsize=12)  # Corrected ylabel
    ax6.set_title(f'Potential Energy Over Time (Avg after 20%: {avg_potential_energy:.2f})', fontsize=14)
    
    # Animation initialization
    def init():
        line1.set_data([], [])
        title.set_text('')
        return line1, title
    
    # Select frames based on skipsteps
    selected_indices = np.arange(0, len(steps), skipsteps)
    selected_steps = [steps[i] for i in selected_indices]
    selected_rdfs = rdfs[selected_indices]
    
    # Animation update
    def animate(i):
        g = selected_rdfs[i]
        line1.set_data(r_values, g)
        title.set_text(f'Step: {selected_steps[i]}')
        return line1, title
    
    ani = animation.FuncAnimation(
        fig, animate, frames=len(selected_steps), init_func=init, interval=100, blit=True
    )
    
    plt.tight_layout()

    # Initialize the writer with total_frames based on skipsteps
    writer = TqdmPillowWriter(fps=10, total_frames=len(selected_steps))  # Adjusted fps to 10 to match interval=100ms

    # Save the animation
    #ani.save('rdf_evolution.gif', writer=writer)
    plt.show()

# Main script
if __name__ == "__main__":
    # Default skipsteps value
    default_skipsteps = 1

    # Parse skipsteps from command-line arguments if provided
    if len(sys.argv) > 1:
        try:
            skipsteps = int(sys.argv[1])
            if skipsteps < 1:
                print("skipsteps must be a positive integer. Using default value of 1.")
                skipsteps = default_skipsteps
        except ValueError:
            print("Invalid skipsteps value provided. Using default value of 1.")
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

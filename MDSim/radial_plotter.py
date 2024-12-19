import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys

# Usage: python animate_rdf.py RDFSpeedVerlet.dat
rdf_file = "radialDistribution.dat"

steps = []
rdfs = []
r_values = None

with open(rdf_file, 'r') as f:
    current_step_data = []
    current_step = None
    for line in f:
        line = line.strip()
        if line.startswith('# Step'):
            # New step
            if current_step_data and current_step is not None:
                # Save previous step data
                rdfs.append(current_step_data)
                steps.append(current_step)
            # Parse step number
            current_step = int(line.split()[2])
            current_step_data = []
        else:
            # r g(r) pairs
            if line:
                vals = line.split()
                r = float(vals[0])
                g = float(vals[1])
                current_step_data.append((r, g))

    # Don't forget the last step
    if current_step_data and current_step is not None:
        rdfs.append(current_step_data)
        steps.append(current_step)

# Assume all steps have the same r values
r_values = [d[0] for d in rdfs[0]]

fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2)
ax.set_xlim(0, max(r_values))
ax.set_ylim(0, 3.0)  # Set according to your expected g(r) range

title = ax.text(0.5, 1.05, '', transform=ax.transAxes, ha='center')
ax.set_xlabel('r')
ax.set_ylabel('g(r)')

def init():
    line.set_data([], [])
    title.set_text('')
    return line, title

def animate(i):
    # i-th frame is i-th step
    rd = rdfs[i]
    r = [x[0] for x in rd]
    g = [x[1] for x in rd]
    line.set_data(r, g)
    title.set_text('Step: {}'.format(steps[i]))
    return line, title

ani = animation.FuncAnimation(fig, animate, frames=len(steps), init_func=init, interval=200, blit=True)
plt.show()

import matplotlib.pyplot as plt

def load_gr(filename):
    r, gr = [], []
    with open(filename, 'r') as file:
        for line in file:
            if line.strip() == "":
                continue
            values = line.strip().split()
            r.append(float(values[0]))
            gr.append(float(values[1]))
    return r, gr

# Load data
r_direct, gr_direct = load_gr('gr_direct.dat')
r_cell, gr_cell = load_gr('gr_cell.dat')

# Plot
plt.figure(figsize=(8,6))
plt.plot(r_direct, gr_direct, label='g(r) Direct', linestyle='-', color='blue')
plt.plot(r_cell, gr_cell, label='g(r) Cell Lists', linestyle='--', color='red')
plt.xlabel('r')
plt.ylabel('g(r)')
plt.title('Radial Distribution Function')
plt.legend()
plt.grid(True)
plt.show()

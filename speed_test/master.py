import subprocess
import matplotlib.pyplot as plt

# Run the C++ benchmark
subprocess.run(["g++", "-o", "cpp_benchmark", "cpp_benchmark.cpp", "-std=c++11", "-Ofast"])
subprocess.run(["./cpp_benchmark"])

# Run the PyTorch and TensorFlow benchmarks
from pytorch_benchmark import benchmark_pytorch
from tensorflow_benchmark import benchmark_tensorflow

benchmark_pytorch(runs=100, filename="pytorch_benchmark.dat")
benchmark_tensorflow(runs=100, filename="tensorflow_benchmark.dat")

def read_benchmark_file(filename):
    times = {}
    with open(filename, 'r') as f:
        for line in f:
            operation, time = line.strip().split(": ")
            times[operation] = float(time)
    return times

pytorch_times = read_benchmark_file("pytorch_benchmark.dat")
tensorflow_times = read_benchmark_file("tensorflow_benchmark.dat")
cpp_times = read_benchmark_file("cpp_benchmark.dat")

# Plot the comparison
labels = list(pytorch_times.keys())
pytorch_values = list(pytorch_times.values())
tensorflow_values = list(tensorflow_times.values())
cpp_values = list(cpp_times.values())

x = range(len(labels))
width = 0.25

plt.figure(figsize=(12, 8))
plt.bar([p - width for p in x], pytorch_values, width=width, label="PyTorch", color='b', align='center')
plt.bar(x, tensorflow_values, width=width, label="TensorFlow", color='g', align='center')
plt.bar([p + width for p in x], cpp_values, width=width, label="C++", color='r', align='center')

plt.xlabel('Operations')
plt.ylabel('Time (microseconds)')
plt.title('Benchmark Comparison')
plt.xticks(x, labels, rotation=45)
plt.legend()
plt.grid(axis='y')

# Add values on top of bars
for i in range(len(labels)):
    plt.text(i - width, pytorch_values[i] + 0.01, f'{pytorch_values[i]:.2f}', ha='center', va='bottom')
    plt.text(i, tensorflow_values[i] + 0.01, f'{tensorflow_values[i]:.2f}', ha='center', va='bottom')
    plt.text(i + width, cpp_values[i] + 0.01, f'{cpp_values[i]:.2f}', ha='center', va='bottom')

plt.tight_layout()
plt.show()

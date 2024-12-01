import subprocess
import time
import numpy as np
import torch
import tensorflow as tf
import matplotlib.pyplot as plt

# Function to run the C++ program and capture its output
def run_cpp():
    # Run the compiled C++ executable and capture its output
    result = subprocess.run(['./cpp_bench.exe'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return result.stdout

# Function to time dot product in PyTorch
def pytorch_dot_product(v1, v2, iterations):
    start_time = time.time()
    for _ in range(iterations):
        torch.dot(v1, v2)
    return time.time() - start_time

# Function to time matrix-vector product in PyTorch
def pytorch_mat_vec_product(mat, vec, iterations):
    start_time = time.time()
    for _ in range(iterations):
        torch.matmul(mat, vec)
    return time.time() - start_time

# Function to time matrix-matrix product in PyTorch
def pytorch_mat_mat_product(mat1, mat2, iterations):
    start_time = time.time()
    for _ in range(iterations):
        torch.matmul(mat1, mat2)
    return time.time() - start_time

# Function to time dot product in TensorFlow
def tensorflow_dot_product(v1, v2, iterations):
    start_time = time.time()
    for _ in range(iterations):
        tf.tensordot(v1, v2, axes=1)
    return time.time() - start_time

# Function to time matrix-vector product in TensorFlow
def tensorflow_mat_vec_product(mat, vec, iterations):
    start_time = time.time()
    for _ in range(iterations):
        tf.linalg.matvec(mat, vec)
    return time.time() - start_time

# Function to time matrix-matrix product in TensorFlow
def tensorflow_mat_mat_product(mat1, mat2, iterations):
    start_time = time.time()
    for _ in range(iterations):
        tf.linalg.matmul(mat1, mat2)
    return time.time() - start_time

# Sizes and data
n = 1000
iterations = 10000
vec1 = torch.ones(n).float()
vec2 = torch.ones(n).float() * 2
mat1 = torch.ones(n, n).float()
mat2 = torch.ones(n, n).float() * 2

# Run C++ program and capture output
cpp_output = run_cpp()
cpp_times = {
    "dot_product": int(cpp_output.split("C++ Dot Product Duration: ")[1].split()[0]),
    "mat_vec": int(cpp_output.split("C++ Matrix-Vector Product Duration: ")[1].split()[0]),
    "mat_mat": int(cpp_output.split("C++ Matrix-Matrix Product Duration: ")[1].split()[0]),
}

# Results for PyTorch
pt_dot = pytorch_dot_product(vec1, vec2, iterations)
pt_mat_vec = pytorch_mat_vec_product(mat1, vec1, iterations)
pt_mat_mat = pytorch_mat_mat_product(mat1, mat2, iterations)

# Results for TensorFlow
tf_dot = tensorflow_dot_product(vec1.numpy(), vec2.numpy(), iterations)
tf_mat_vec = tensorflow_mat_vec_product(mat1.numpy(), vec1.numpy(), iterations)
tf_mat_mat = tensorflow_mat_mat_product(mat1.numpy(), mat2.numpy(), iterations)

# Plotting the comparison
labels = ['Dot Product', 'Matrix-Vector Product', 'Matrix-Matrix Product']
cpp_times_list = [cpp_times["dot_product"], cpp_times["mat_vec"], cpp_times["mat_mat"]]
pytorch_times = [pt_dot, pt_mat_vec, pt_mat_mat]
tensorflow_times = [tf_dot, tf_mat_vec, tf_mat_mat]

x = np.arange(len(labels))

fig, ax = plt.subplots()
ax.bar(x - 0.2, cpp_times_list, 0.4, label='C++')
ax.bar(x, pytorch_times, 0.4, label='PyTorch')
ax.bar(x + 0.2, tensorflow_times, 0.4, label='TensorFlow')

ax.set_ylabel('Time (seconds)')
ax.set_title('Time comparison of operations (10,000 iterations)')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()

plt.show()

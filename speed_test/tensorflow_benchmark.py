import time
import tensorflow as tf
from tqdm import tqdm

def benchmark_tensorflow(runs=100, filename="tensorflow_benchmark.dat"):
    times = {'Scalar Product': 0, 'Matrix-Vector Multiplication': 0, 'Matrix-Matrix Multiplication': 0}
    
    n = 1000
    v1 = tf.random.uniform([n], dtype=tf.float32)
    v2 = tf.random.uniform([n], dtype=tf.float32)
    M = tf.random.uniform([n, n], dtype=tf.float32)

    for _ in tqdm(range(runs), total=runs, desc='TensorFlow Benchmark'):
        start_time = time.time()
        scalar_product = tf.reduce_sum(v1 * v2)
        times['Scalar Product'] += (time.time() - start_time) * 1e6

        start_time = time.time()
        result_vec = tf.matmul(M, tf.expand_dims(v1, -1))
        times['Matrix-Vector Multiplication'] += (time.time() - start_time) * 1e6

        M2 = tf.random.uniform([n, n], dtype=tf.float32)
        start_time = time.time()
        result_mat = tf.matmul(M, M2)
        times['Matrix-Matrix Multiplication'] += (time.time() - start_time) * 1e6

    for key in times:
        times[key] /= runs

    with open(filename, 'w') as f:
        for key, value in times.items():
            f.write(f"{key}: {value}\n")

    return times


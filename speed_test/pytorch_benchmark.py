import time
import torch
from tqdm import tqdm

def benchmark_pytorch(runs=100, filename="pytorch_benchmark.dat"):
    times = {'Scalar Product': 0, 'Matrix-Vector Multiplication': 0, 'Matrix-Matrix Multiplication': 0}
    
    n = 1000
    v1 = torch.rand(n)
    v2 = torch.rand(n)
    M = torch.rand(n, n)

    for _ in tqdm(range(runs), total=runs, desc='PyTorch Benchmark'):
        start_time = time.time()
        scalar_product = torch.dot(v1, v2)
        times['Scalar Product'] += (time.time() - start_time) * 1e6

        start_time = time.time()
        result_vec = torch.matmul(M, v1)
        times['Matrix-Vector Multiplication'] += (time.time() - start_time) * 1e6

        M2 = torch.rand(n, n)
        start_time = time.time()
        result_mat = torch.matmul(M, M2)
        times['Matrix-Matrix Multiplication'] += (time.time() - start_time) * 1e6

    for key in times:
        times[key] /= runs

    with open(filename, 'w') as f:
        for key, value in times.items():
            f.write(f"{key}: {value}\n")

    return times


from joblib import Parallel, delayed
from cmap import cmap, cmap_improve
import os, sys
from cmapPy.pandasGEXpress.parse import parse
import numpy as np
import random
import time
import json

def load_matrix_from_file(file_path):
    with open(file_path, 'r') as file:
        return [[int(value) for value in line.split()] for line in file]

reference_data_file = sys.argv[1]
begin = int(sys.argv[2])
end = int(sys.argv[3])
data = parse(reference_data_file).data_df
data.index = [int(i) for i in data.index]
nearest_neighbor_euclidean = np.loadtxt("../distance_files/nearest_neighbor_euclidean.txt", dtype=int)
nearest_neighbor_cosine = np.loadtxt("../distance_files/nearest_neighbor_cosine.txt", dtype=int)
nearest_neighbor_manhattan = np.loadtxt("../distance_files/nearest_neighbor_cityblock.txt", dtype=int)
result_dir = '../result'
os.makedirs(result_dir, exist_ok=True)

def convert_to_int(arr):
    for i in range(len(arr)):
        arr[i] = int(arr[i])
    return arr

def get_random_gene_list():
    n_up_genes = random.randint(10, 150)
    n_down_genes = random.randint(10, 150)
    rand_genes = random.choices(data.index, k = n_up_genes + n_down_genes)
    return rand_genes[:n_up_genes], rand_genes[n_up_genes:]

def run_test(i):
    up_genes_list, down_genes_list = get_random_gene_list()
    up_genes = set(up_genes_list)
    down_genes = set(down_genes_list)
    start_time = time.time()
    original_result = cmap(up_genes, down_genes, data)
    original_time = time.time()
    euclidean_result = cmap_improve(up_genes, down_genes, data, nearest_neighbor_euclidean)
    euclidean_time = time.time()
    manhattan_result = cmap_improve(up_genes, down_genes, data, nearest_neighbor_manhattan)
    manhattan_time = time.time()
    cosine_result = cmap_improve(up_genes, down_genes, data, nearest_neighbor_cosine)
    cosine_time = time.time()
    result = {}
    result['up_genes'] = convert_to_int(up_genes_list)
    result['down_genes'] = convert_to_int(down_genes_list)
    result['original_mimic_result'] = original_result['expression'].head(100).values.tolist()
    result['original_mimic_score'] = original_result['c_score'].head(100).values.tolist()
    result['original_running_time'] = original_time - start_time
    result['euclidean_mimic_result'] = euclidean_result['expression'].head(100).values.tolist()
    result['euclidean_mimic_score'] = euclidean_result['c_score'].head(100).values.tolist()
    result['euclidean_running_time'] = euclidean_time - original_time
    result['manhattan_mimic_result'] = manhattan_result['expression'].head(100).values.tolist()
    result['manhattan_mimic_score'] = manhattan_result['c_score'].head(100).values.tolist()
    result['manhattan_running_time'] = manhattan_time - euclidean_time
    result['cosine_mimic_result'] = cosine_result['expression'].head(100).values.tolist()
    result['cosine_mimic_score'] = cosine_result['c_score'].head(100).values.tolist()
    result['cosine_running_time'] = cosine_time - manhattan_time
    with open(f'{result_dir}/output_{i}.json', 'w') as file:
        json.dump(result, file)
    print(f'output_{i}.json file created')

Parallel(n_jobs=-1)(
    delayed(run_test)(i) for i in range(begin, end)
)

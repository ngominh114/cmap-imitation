import numpy as np
from cmapPy.pandasGEXpress.parse import parse
from sklearn.metrics import pairwise_distances
from joblib import Parallel, delayed
import os
import sys
from os import listdir
from os.path import isfile, join

data_path = sys.argv[1]

data = parse(data_path).data_df
data = data.astype('float32').to_numpy().transpose()
output_directory = "../distance_files"
tmp_dicectory = "../tmp"
os.makedirs(output_directory, exist_ok=True)
os.makedirs(tmp_dicectory, exist_ok=True)

n_rows = len(data)

metrics = ['euclidean', 'cosine', 'cityblock']

def file_index(file):
    dot_index = file.index(".")
    return int(file[17:dot_index])

def calculate_and_save_distances(i, metric):
    start_row = i * 10000
    end_row = min((i+1)*10000, n_rows)
    subset_data = data[start_row: end_row]
    distances = pairwise_distances(subset_data, data, metric=metric)
    top_indices = np.argsort(distances, axis=1)[:, 1:51]
    np.savetxt(f"{tmp_dicectory}/nearest_neighbor_{i}.txt", top_indices, fmt='%d')

for metric in metrics:
    Parallel(n_jobs=-1)(
        delayed(calculate_and_save_distances)(i, metric) for i in range(n_rows//10000 + 1)
    )

    files = [f for f in listdir(tmp_dicectory) if isfile(join(tmp_dicectory, f))]
    files.sort(key=lambda f: file_index(f))
    nearest_neighbor = []
    for file in files:
        file_dir = tmp_dicectory + "/" + file
        with open(file_dir) as f:
            for line in f:
                nearest_neighbor.append(list(map(int, line.split(" "))))
        os.remove(file_dir)
    np.savetxt(f"{output_directory}/nearest_neighbor_{metric}.txt", nearest_neighbor, fmt='%d')


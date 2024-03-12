from cmap import cmap
from cmapPy.pandasGEXpress.parse import parse
import os, sys, random, multiprocessing
from joblib import Parallel, delayed
import numpy as np

reference_data_file = sys.argv[1]
data = parse(reference_data_file).data_df
data.index = [int(i) for i in data.index]

result_dir = '../original_result'
os.makedirs(result_dir, exist_ok=True)

def get_random_gene_list():
    n_up_genes = random.randint(10, 150)
    n_down_genes = random.randint(10, 150)
    rand_genes = random.choices(data.index, k = n_up_genes + n_down_genes)
    return rand_genes[:n_up_genes], rand_genes[n_up_genes:]

result = multiprocessing.Manager().list()

def run_cmap(result):
    up_genes_list, down_genes_list = get_random_gene_list()
    up_genes = set(up_genes_list)
    down_genes = set(down_genes_list)
    res = cmap(up_genes, down_genes, data)
    print(res)
    result.append(res['c_score'].values.tolist())

Parallel(n_jobs=-1)(
    delayed(run_cmap)(result) for _ in range(1000)
)

print(result)

np.savetxt(f'{result_dir}/original_result.txt', result, fmt='%0.2f')
import time
import datetime
import sys
from cmapPy.pandasGEXpress.parse import parse
from cmap import cmap, cmap_improve
import pandas as pd
import numpy as np

start_time = time.time()
query_name = sys.argv[1]
reference_data_file = sys.argv[2] # reference data in gctx format
gene_info_file = sys.argv[3] # ../data/geneinfo_beta.txt
up_gene_file = sys.argv[4]   # ../data/Example_up_genes.txt
down_gene_file = sys.argv[5] # ../data/Example_down_genes.txt

print("Start calculating query:", query_name, "with the following arguments.")
print("Reference dataset location:", reference_data_file)
print("Gene info file location:", gene_info_file)
print("Up gene file location:", up_gene_file)
print("Down gene file location:", down_gene_file)

data = parse(reference_data_file).data_df
data.index = data.index.astype(int)
reference_gene_list = [gene for gene in data.index]
gene_info = pd.read_csv(gene_info_file, delimiter="\t")
gene_dict = {}
for i in range(len(gene_info.index)):
    gene_dict[str(gene_info["gene_symbol"][i])] = gene_info["gene_id"][i]

def read_genes(file_name):
    with open(file_name) as file:
        return [gene.strip() for gene in file.readlines()]

up_genes = read_genes(up_gene_file)
down_genes = read_genes(down_gene_file)

up_genes = set([gene_dict[i] for i in up_genes if i in gene_dict])
down_genes = set([gene_dict[i] for i in down_genes if i in gene_dict])

nearest_neighbor_euclidean = np.loadtxt("../distance_files/nearest_neighbor_euclidean.txt", dtype=int)
nearest_neighbor_cosine = np.loadtxt("../distance_files/nearest_neighbor_cosine.txt", dtype=int)

result = cmap(up_genes, down_genes, data)
result_euc = cmap_improve(up_genes, down_genes, data, nearest_neighbor_euclidean)
result_cos = cmap_improve(up_genes, down_genes, data, nearest_neighbor_cosine)

end_time = time.time()
print("Total running time is", str(datetime.timedelta(seconds=end_time - start_time)))
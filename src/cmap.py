import sys
from cmapPy.pandasGEXpress.parse import parse
import numpy as np
import pandas as pd
import time
import datetime

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

n_expression = len(data.columns)
n_ref_gene = len(data.index)

def calc_es_score(ref_df, cid, gene_list):
    gene_indexes = ref_df.index.isin(gene_list)
    ns = len(gene_list)
    nr = ref_df.loc[gene_indexes, cid].sum()
    cumsum_score = np.zeros(n_ref_gene)
    cumsum_score[gene_indexes] = ref_df.loc[gene_indexes, cid] / nr
    cumsum_score[~gene_indexes] = -1 / (n_ref_gene - ns)
    scores = np.cumsum(cumsum_score)
    min_score = np.min(scores)
    max_score = np.max(scores)
    return max_score if abs(max_score) > abs(min_score) else min_score

def calc_connectivity_score(up_gene, down_gene, df, cid):
    up_es_score = calc_es_score(df, cid, up_gene)
    down_es_score = calc_es_score(df, cid, down_gene)
    return (up_es_score - down_es_score) / 2 if up_es_score * down_es_score < 0 else 0

def cmap(up_gene_list, down_gene_list):
    
    up_gene_list = [gene_dict[i] for i in up_gene_list if i in gene_dict]
    down_gene_list = [gene_dict[i] for i in down_gene_list if i in gene_dict]
    c_scores = []
    for cid in data.columns:
        df = data[cid]
        df = pd.DataFrame(df.sort_values(ascending=False)).abs()
        c_scores.append(calc_connectivity_score(up_gene_list, down_gene_list, df, cid))
    return pd.DataFrame({"expression": data.columns, "c_score": c_scores})

result = cmap(up_genes, down_genes)
result = result.sort_values(by=["c_score"], ascending=False)
result.to_csv(f"{query_name}_py.csv")

end_time = time.time()
print("Total running time is", str(datetime.timedelta(seconds=end_time - start_time)))

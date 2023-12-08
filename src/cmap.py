import sys
from cmapPy.pandasGEXpress.parse import parse
import numpy as np
import pandas as pd
import time
import datetime

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
    cumsum_score = 0
    ref_df = ref_df.loc[gene_list].sort_values(by=["ind"])
    es_scores = []
    ns = len(gene_list)
    nr = ref_df[cid].sum()
    for i in range(ns):
        if i == 0: pre = -1
        else: pre = ref_df['ind'].iloc[i-1]

        cumsum_score -= (ref_df['ind'].iloc[i] - pre - 1)/(n_ref_gene-ns)
        es_scores.extend([cumsum_score])
        cumsum_score += ref_df[cid].iloc[i]/nr
        es_scores.extend([cumsum_score])
        if i == ns-1:
            cumsum_score -= (n_ref_gene - 1 - ref_df['ind'].iloc[i])/(n_ref_gene-ns)
            es_scores.extend([cumsum_score])
    max_es = max(es_scores)
    min_es = min(es_scores)
    return max_es if abs(max_es) > abs(min_es) else min_es

def cmap(up_gene_list, down_gene_list):
    start_time = time.time()
    up_gene_list = [gene_dict[i] for i in up_gene_list if i in gene_dict]
    down_gene_list = [gene_dict[i] for i in down_gene_list if i in gene_dict]
    c_scores = []
    for cid in data.columns:
        df = data[cid]
        df = pd.DataFrame(df.sort_values(ascending=False)).abs()
        df['ind'] = [i for i in range(n_ref_gene)]
        up_es_score = calc_es_score(df, cid, up_gene_list)
        down_es_score = calc_es_score(df, cid, down_gene_list)
        if up_es_score*down_es_score < 0:
            c_scores.append((up_es_score - down_es_score)/2)
        else:
            c_scores.append(0)
    end_time = time.time()
    print("Total running time is", str(datetime.timedelta(seconds=end_time - start_time)))
    return pd.DataFrame({"expression": data.columns, "c_score": c_scores})

result = cmap(up_genes, down_genes)
result = result.sort_values(by=["c_score"], ascending=False)
result.to_csv(f"{query_name} result.csv")

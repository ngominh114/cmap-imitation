import sys
from cmapPy.pandasGEXpress.parse import parse
import pandas as pd
import time
import datetime

start_time = time.time()
query_name = sys.argv[1]
reference_data_file = sys.argv[2]
gene_info_file = sys.argv[3]
up_gene_file = sys.argv[4]
down_gene_file = sys.argv[5]

print("Start calculating query:", query_name, "with the following arguments.")
print("Reference dataset location:", reference_data_file)
print("Gene info file location:", gene_info_file)
print("Up gene file location:", up_gene_file)
print("Down gene file location:", down_gene_file)

data = parse(reference_data_file).data_df
gene_info = pd.read_csv(gene_info_file, delimiter="\t")
gene_dict = {}
for i in range(len(gene_info.index)):
    gene_dict[str(gene_info["gene_id"][i])] = gene_info["gene_symbol"][i]

def read_genes(file_name):
    with open(file_name) as file:
        return [gene.strip() for gene in file.readlines()]

up_genes = read_genes(up_gene_file)
down_genes = read_genes(down_gene_file)

n_expression = len(data.columns)
n_ref_gene = len(data.index)

def calc_es_score(ref_df, gene_list):
    ref_df = abs(ref_df)
    cumsum_score = 0
    es_score = 0
    n = len(ref_df.index)
    ns = len(gene_list)
    nr = sum(ref_df)
    for i in range(n):
        if gene_dict[ref_df.index[i]] in gene_list:
            cumsum_score += ref_df[i]/nr
        else:
            cumsum_score -= 1/(n-ns)

        if abs(cumsum_score) > abs(es_score): es_score = cumsum_score
    return es_score
def cmap(up_gene_list, down_gene_list):
    up_gene_list = set(up_gene_list)
    down_gene_list = set(down_gene_list)
    c_scores = []
    for cid in data.columns:
        df = data[cid]
        df = df.sort_values(ascending=False)
        up_es_score = calc_es_score(df, up_gene_list)
        down_es_score = calc_es_score(df, down_gene_list)
        if up_es_score*down_es_score < 0:
            c_scores.append((up_es_score - down_es_score)/2)
        else:
            c_scores.append(0)
    return pd.DataFrame({"expression": data.columns, "c_score": c_scores})

result = cmap(up_genes, down_genes)
result.to_csv(f"{query_name} result.csv")
end_time = time.time()

print("Total running time is", str(datetime.timedelta(end_time - start_time)))
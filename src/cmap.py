import numpy as np
import pandas as pd
import random

def calc_es_score(ref_df, cid, gene_list):
    n_genes = len(ref_df.index)
    gene_indexes = ref_df.index.isin(gene_list)
    ns = len(gene_list)
    nr = ref_df.loc[gene_indexes, cid].sum()
    cumsum_score = np.zeros(n_genes)
    cumsum_score[gene_indexes] = ref_df.loc[gene_indexes, cid] / nr
    cumsum_score[~gene_indexes] = -1 / (n_genes - ns)
    scores = np.cumsum(cumsum_score)
    min_score = np.min(scores)
    max_score = np.max(scores)
    return max_score if abs(max_score) > abs(min_score) else min_score

def calc_connectivity_score(up_gene, down_gene, data, cid):
    df = data[cid]
    df = pd.DataFrame(df.sort_values(ascending=False)).abs()
    up_es_score = calc_es_score(df, cid, up_gene)
    down_es_score = calc_es_score(df, cid, down_gene)
    return (up_es_score - down_es_score) / 2 if up_es_score * down_es_score < 0 else 0

def cmap(up_gene_list, down_gene_list, data):
    c_scores = []
    for cid in data.columns:
        c_score = calc_connectivity_score(up_gene_list, down_gene_list, data, cid)
        c_scores.append(c_score)
    
    return pd.DataFrame({"expression": data.columns, "c_score": c_scores}).sort_values(by=["c_score"], ascending=False)

check = set()
def cmap_improve(up_gene_list, down_gene_list, data, nearest_neighbor_info):
    checked = set()
    proceeded = []
    c_scores = []
    for i, cid in enumerate(data.columns):
        if i not in checked:
            c_score = calc_connectivity_score(up_gene_list, down_gene_list, data, cid)
            proceeded.append(cid)
            c_scores.append(c_score)
            if c_score == 0:
                checked.update(nearest_neighbor_info[i,:10])
        
    return pd.DataFrame({"expression": proceeded, "c_score": c_scores}).sort_values(by=["c_score"], ascending=False)

def cmap_random(up_gene_list, down_gene_list, data, indexes):
    checked = set()
    proceeded = []
    c_scores = []
    for i in indexes:
        cid = data.columns[i]
        if i not in checked:
            c_score = calc_connectivity_score(up_gene_list, down_gene_list, data, cid)
            proceeded.append(cid)
            c_scores.append(c_score)
            if c_score == 0:
                checked.update(random.choices(data.index, k=10))

    return pd.DataFrame({"expression": proceeded, "c_score": c_scores}).sort_values(by=["c_score"], ascending=False)
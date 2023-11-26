#!/bin/bash

query_name=$1
reference_data_file=$2
gene_info_file=$3
up_gene_file=$4
down_gene_file=$5

conda activate cmap

python -u ../src/cmap.py $query_name $reference_data_file $gene_info_file $up_gene_file $down_gene_file
#!/bin/bash
#SBATCH --job-name=cmap_run
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=512GB
#SBATCH --time=01:00:00


query_name=$1
reference_data_file=$2
gene_info_file=$3
up_gene_file=$4
down_gene_file=$5

conda activate cmap

python -u ../src/cmap.py $query_name $reference_data_file $gene_info_file $up_gene_file $down_gene_file
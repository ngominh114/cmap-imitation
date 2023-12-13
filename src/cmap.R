if (!require("cmapR", character.only = TRUE)){
  install.packages("cmapR")
}

library(cmapR)

# reference_file = "C:/Users/ngomi/Downloads/Documents/data.gctx"
# up_gene_file = "../data/Example_up_genes.txt"
# down_gene_file = "../data/Example_down_genes.txt"
# gene_file = "../data/geneinfo_beta.txt"
args = commandArgs(trailingOnly=TRUE)
query_name = args[1]
reference_file = args[2]
up_gene_file = args[3]
down_gene_file = args[4]
gene_file = args[5]
col_meta = read_gctx_meta(reference_file, dim="col")
cids = col_meta$id
row_meta = read_gctx_meta(reference_file, dim="row")

data = mat(parse_gctx(reference_file))

up_gene_list = read.table(up_gene_file, header = FALSE)$V1
down_gene_list = read.table(down_gene_file, header = FALSE)$V1
gene_data = read.csv(gene_file, header = TRUE, sep = "\t")

gene_dict = list()
for(i in 1:length(gene_data$gene_id)){
  gene_dict[[gene_data$gene_symbol[i]]] <- gene_data$gene_id[i][1]
}

for(i in 1:length(up_gene_list)){
  if(up_gene_list[i] %in% names(gene_dict)){
    up_gene_list[i] = as.integer(gene_dict[up_gene_list[i]])
  }
}

for(i in 1:length(down_gene_list)){
  if(down_gene_list[i] %in% names(gene_dict)){
    down_gene_list[i] = as.integer(gene_dict[down_gene_list[i]])
  }
}

up_gene_list = as.integer(up_gene_list)
down_gene_list = as.integer(down_gene_list)

calc_es_score = function(reference_df, gene_list){
  gene_indexes <- reference_df$gene %in% gene_dict[gene_list]
  nr <- sum(abs(reference_df[gene_indexes, "expression"]))
  n <- nrow(reference_df)
  ns <- length(gene_list)
  
  cumsum_score <- rep(0, n)
  cumsum_score[gene_indexes] <- abs(reference_df[gene_indexes, "expression"]) / nr
  cumsum_score[!gene_indexes] <- -1 /(n - ns)
  
  es_score <- max(cumsum(cumsum_score))
  
  return(es_score)
}

cmap = function(){
  n = length(cids)
  c_scores = numeric(n)
  genes = as.integer(row_meta$id)
  for(i in 1:n){
    m <- data[, i]
    ordered_index <- order(m, decreasing = TRUE)
    df <- data.frame("gene" = genes, "expression" = m[ordered_index])
    es_up <- calc_es_score(df, up_gene_list)
    es_down <- calc_es_score(df, down_gene_list)
    if (es_up * es_down < 0) {
      c_scores[i] <- (es_up - es_down) / 2
    }
  }
  result_df <- data.frame("expression" = cids, "connectivity_score" = c_scores)
  return(result_df)
}

result = cmap()
write.csv(result, paste(query_name, "result.csv"))
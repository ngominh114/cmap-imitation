library(cmapR)

reference_file = "C:/Users/ngomi/Downloads/Documents/level5_beta_trt_oe_n34171x12328.gctx"
up_gene_file = "../data/Example up genes.txt"
down_gene_file = "../data/Example down genes.txt"
gene_file = "../data/geneinfo_beta.txt"
# args = commandArgs(trailingOnly=TRUE)
# reference_file = args[1]
# up_gene_file = args[2]
# down_gene_file = args[3]
typeof(gene_data$gene_id[1])
gene_dict = list()
for(i in 1:length(gene_data$gene_id)){
  gene_dict[[gene_data$gene_id[i]]] <- gene_data$gene_symbol[i]
}

col_meta = read_gctx_meta(reference_file, dim="col")
cids = col_meta$id
up_gene_list = read.table(up_gene_file, header = FALSE)$V1
down_gene_list = read.table(down_gene_file, header = FALSE)$V1
gene_data = read.table(gene_file, header = TRUE, sep = "\t")

calc_es_score_1 = function(reference_df, gene_list){
  cumsum_arr = c(0)
  cumsum_score = 0
  es_score = 0
  n = nrow(reference_df)
  n_s = length(gene_list)
  
  for(i in 1:n){
    if(gene_dict[as.integer(reference_df$gene[i])] %in% gene_list){
      cumsum_score = cumsum_score + sqrt((n - n_s) / n_s)
    }else{
      cumsum_score = cumsum_score - sqrt(n_s / (n_s + n))
    }
    cumsum_arr = c(cumsum_arr, cumsum_score)
    print(cumsum_score)
    break
    if(cumsum_score < 0){
      es_score = max(es_score, abs(cumsum_score))
    }else{
      es_score = max(es_score, cumsum_score)
    }
  }
  print(n)
  print(length(cumsum_arr))
  plot(1:(n+1), cumsum_arr)
  return(es_score)
}

calc_es_score_2 = function(reference_df, gene_list){
  cumsum_score = 0
  es_score = 0
  n = nrow(reference_df)
  ns = length(gene_list)
  nr = 0
  for(i in 1:n){
    if(gene_dict[as.integer(reference_df$gene[i])] %in% gene_list){
      nr = nr + abs(reference_df$expression[i])
    }
  }

  for(i in 1:n){
    if(gene_dict[as.integer(reference_df$gene[i])] %in% gene_list){
      cumsum_score = cumsum_score + abs(reference_df$expression[i])/nr
    }else{
      cumsum_score = cumsum_score - 1/(n - ns)
    }
    if(abs(cumsum_score) > abs(es_score)){
      es_score = cumsum_score
    }
  }
  return(es_score)
}

cmap = function(){
  n = length(cids)
  c_scores = rep(0, n)
  for(i in 1:n){
    cid = cids[i]
    data = parse_gctx(reference_file, cid = cid)
    m = mat(data)
    ordered_index = order(m[,1], decreasing = TRUE)
    genes = row.names(m)
    df = data.frame("gene" = genes, "expression" = c(m))
    df = df[order(df$expression, decreasing = TRUE), ]
    es_up = calc_es_score_2(df, up_gene_list)
    es_down = calc_es_score_2(df, down_gene_list)
    if(es_up*es_down < 0){
      c_scores[i] = (es_up - es_down) / 2
    }
  }
  data.frame("cid"=cids, "connectivity_score" = c_scores)
}
result = cmap()

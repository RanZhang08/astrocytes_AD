#!/usr/bin/env Rscript
library(data.table)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
outdir = args[1]
outname = args[2]
diff_input = args[3]
contamination_input = args[4]
option = args[5]


## ===================================================
## 1. train - generate gstd
generate_netdiff_gstd <- function(outdir, outname, diff_input, contamination_input){
  diff_gene_mat <- fread(diff_input)
  contamination_genes <- fread(contamination_input)
  
  ## output diff for rank-based enrichment analysis
  diff_gene_mat <- diff_gene_mat[!gene %in% contamination_genes]
  diff_gene_mat <- diff_gene_mat[!is.na(PValue)]
  diff_gene_mat$FDR <- p.adjust(diff_gene_mat$PValue, method = 'BH')
  diff_gene_mat <- unique(diff_gene_mat[, c('gene','PValue','FDR')])
  
  # define a set of genes to be used as positive gstds
  diff_gene_mat_diff_fil_human_genes <- unique(diff_gene_mat[PValue<=0.05]$gene)
  fdr_cut_n <- length(diff_gene_mat_diff_fil_human_genes)
  if (fdr_cut_n < 1000){
    fdr_cut_n <- 1000
  }
  diff_gene_mat_m_top <- diff_gene_mat[order(PValue)][1:fdr_cut_n, ]
  
  # define a set of genes to be used as negative gstds
  diff_gene_mat_m_bottom <- diff_gene_mat[PValue>0.1 & !gene %in% diff_gene_mat_m_top$gene]
  fdr_cut_n_bottom <- nrow(diff_gene_mat_m_bottom)
  
  ## bootstrap: subsample from the overall positive and negative gold-standards gene sets 50 times. In each model, we subsample goldstandards according to a weight derived from p-value.
  for (randi in 1:50){
    diff_gene_mat_diff_fil_human_genes_b <- unique(sample(diff_gene_mat_m_top$gene, fdr_cut_n/3, prob=(diff_gene_mat_m_top$PValue)^(-0.2)))
    diff_gene_mat_nondiff_fil_human_genes_b <- setdiff(unique(sample(diff_gene_mat_m_bottom$gene, fdr_cut_n_bottom/3, prob=diff_gene_mat_m_bottom$PValue)), diff_gene_mat_m_top$human)
    npos <- length(diff_gene_mat_diff_fil_human_genes_b)
    nneg <- length(diff_gene_mat_nondiff_fil_human_genes_b)
    
    diff_gene_mat_svm_b <- rbind(cbind(diff_gene_mat_diff_fil_human_genes_b, rep(1, npos)),
                                    cbind(diff_gene_mat_nondiff_fil_human_genes_b, rep(-1, nneg)))
    
    write.table(diff_gene_mat_svm_b, paste0(outdir, '/gstds/', outname, '_', randi, '.txt'), sep='\t', quote=F, row.names=F, col.names=F)
  }
}


## ===================================================
## 2. generate overall pred through bagging:
summarize_result <- function(ourdir, outname){
  rm(netdiff_pred_mat)
  for (randi in 1:50){
      mat <- fread(outdir, '/gstds/', outname, '_gstd_', randi, '.txt')
      names(mat)[3] <- randi
      if (exists('netdiff_pred_mat')){
        netdiff_pred_mat <- merge(netdiff_pred_mat, mat[, c(1,3)], by='V1')
      }else{
        netdiff_pred_mat <- mat
      }
    }
  netdiff_pred_mat$score <- apply(netdiff_pred_mat[, 3:ncol(netdiff_pred_mat)], 1, mean)
  netdiff_pred_mat <- netdiff_pred_mat[, c('V1','V2','score')]
  colnames(netdiff_pred_mat) <- c('gene','std','score')
  write.table(netdiff_pred_mat, paste0(outdir, outname, '_predidction.txt'), sep='\t', quote=F, row.names=F, col.names=F)
}


if (option=='generate_gstd'){
  generate_netdiff_gstd(outdir, outname, diff_input, contamination_input)
}else if (option=='summarize'){
  summarize_result(outdir, outname)
}


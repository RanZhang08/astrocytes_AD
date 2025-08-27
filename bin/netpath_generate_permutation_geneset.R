#!/usr/bin/env Rscript
library(data.table)

args = commandArgs(trailingOnly=TRUE)
geneset = args[1]
geneset_bg = args[2]

## =========================================
## generate random genes on both query and target genesets for permutation analysis
## =====================================
bc_bg_genes <- fread(geneset_bg, header=F)$V1
regional_genes <- fread(geneset, header=F)$V1
for (seedi in 0:4){
  rand_mat <- c()
  for (randi in 1:2000){
    rand_mat <- rbind(rand_mat, sample(bc_bg_genes, length(regional_genes)))
  }
  write.table(rand_mat, paste0(substr(geneset,1,nchar(geneset)-4),'_rand', seedi), quote=F, row.names=F, sep='\t', col.names=F)
}


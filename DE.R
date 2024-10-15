library(DESeq2)
library(tibble)
library(magrittr)
library(dplyr)


perform_de <- function(counts, metadata){
  dds <- DESeqDataSetFromMatrix(counts, colData = metadata, design = ~ Status)
  
  dds <- DESeq(dds)
  
  res <- results(dds, 
                 alpha = 0.05)
  
  res <- lfcShrink(dds, coef = resultsNames(dds)[2], res = res, type = "apeglm")
  
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble() %>%
    arrange(padj)
  
  sig_genes <- res_tbl %>%
    filter(padj < 0.05) %>%
    select(gene) %>%
    pull(gene)
  
  return(sig_genes)
  
  
  print(paste0("Number of Differentially Expressed Genes: ", length(sig_genes)))
}
  


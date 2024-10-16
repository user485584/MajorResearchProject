# Load necessary libraries
library(DESeq2)
library(tibble)
library(magrittr)
library(dplyr)


perform_de <- function(counts, metadata){
  
  # ----------------------------
  # Step 1: Input Validation
  # ----------------------------
  # Check if 'Status' column exists in metadata
  if (!"Status" %in% colnames(metadata)) {
    stop("Error: Metadata must contain a 'Status' column.")
  }
  
  # Convert 'Status' to a factor and ensure it has exactly two levels
  metadata$Status <- as.factor(metadata$Status)
  
  if (nlevels(metadata$Status) != 2) {
    stop("Error: The 'Status' column must have exactly two levels (e.g., 'Healthy' and 'Diseased').")
  }
  
  
  # ----------------------------
  # Step 2: Create DESeq2 Dataset
  # ----------------------------
  
  dds <- DESeqDataSetFromMatrix(counts, colData = metadata, design = ~ Status)
  
  
  # ----------------------------
  # Step 3: Differential Expression Analysis
  # ----------------------------
  # Run DESeq normalisation and differential expression analysis
  dds <- DESeq(dds)
  
  # Extract results with specified significance level (alpha)
  res <- results(dds, 
                 alpha = 0.05)
  
  # Perform log fold change shrinkage using the 'apeglm' method
  res <- lfcShrink(dds, coef = resultsNames(dds)[2], res = res, type = "apeglm")
  
  
  # ----------------------------
  # Step 4: Convert Results
  # ----------------------------
  
  
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble() %>%
    arrange(padj)
  
  
  # ----------------------------
  # Step 5: Identify Significantly Differentially Expressed Genes
  # ----------------------------
  
  sig_genes <- res_tbl %>%
    filter(padj < 0.05) %>%
    filter(abs(log2FoldChange) > log2(1.5)) %>%
    select(gene) %>%
    pull(gene)
  
  
  # ----------------------------
  # Step 6: Output Results
  # ----------------------------
  
  print(paste0("Number of Differentially Expressed Genes: ", length(sig_genes)))
  
  return(sig_genes)
}
  


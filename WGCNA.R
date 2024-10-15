# Load necessary libraries
library(WGCNA)
library(magrittr)
library(dplyr)
library(DESeq2)

perform_wgcna <- function(counts, metadata, includeCorrelation = FALSE, power = NULL){
  
  # ====================================
  # 1. Data Preparation and Quality Control
  # ====================================
  
  # Transpose counts to have samples as rows and genes as columns
  t_counts <- t(counts)
  
  # Perform quality control to identify good samples and genes
  gsg <- goodSamplesGenes(t_counts, verbose = 3)
  print(summary(gsg))
  
  # Remove bad samples and genes if any are identified
  if (!gsg$allOK) {
    # Filter out bad samples and genes
    if (any(!gsg$goodSamples)) {
      print("Removing bad samples")
      t_counts <- t_counts[gsg$goodSamples, ]
    }
    if (any(!gsg$goodGenes)) {
      print("Removing bad genes")
      t_counts <- t_counts[, gsg$goodGenes]
    }
  }
  
  
  # ====================================
  # 2. Normalisation and Transformation
  # ====================================
  
  #Normalisation
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = metadata,
                                design = ~ 1)
  
  #Filter genes: retain genes with a count of 15 or more in at least 10 samples
  print("Filtering genes")
  dds75 <- dds[rowSums(counts(dds) >= 15) >= 10,]
  print(paste0("Number of remaining Genes: ",nrow(dds75)))
  
  # Apply Variance Stabilising Transformation (VST)
  dds_norm <- vst(dds75)
  
  # Extract normalized counts and transpose to have genes as columns
  norm_counts <- assay(dds_norm) %>% t()
  
  
  # ====================================
  # 3. Network Construction with WGCNA
  # ====================================
  
  # Determine the soft-thresholding power if not provided
  if (is.null(power)) {
    powerOptions <- c(1:10, seq(from = 12, to = 30, by = 2))
    sft <- pickSoftThreshold(norm_counts,
                             powerVector = powerOptions,
                             networkType = "signed",
                             verbose = 5)
    softPower <- sft$powerEstimate
  } else {
    softPower <- power
  }
  
  # Temporarily override the cor function to use WGCNA's cor
  original_cor <- cor
  cor <- WGCNA::cor
  
  
  # Identify modules using blockwiseModules
  bwnet <- blockwiseModules(norm_counts, maxBlockSize = 25000,
                            TOMType = 'signed', power = softPower,
                            mergeCutHeight = 0.25, numericLabels = FALSE,
                            verbose = 3)
  
  # Restore the original cor function
  cor <- original_cor
  
  # Extract module eigengenes
  module_eigengenes <- bwnet$MEs
  
  # ====================================
  # 4. Trait Correlation and Module Selection
  # ====================================
  
  # Prepare trait data: Binary variable indicating disease presence
  traits <- metadata %>% 
    mutate(Disease = ifelse(grepl("ALS|AD", Status), 1, 0)) %>% 
    select(Disease)
  
  nSamples <- nrow(norm_counts)
  
  # Calculate correlation between module eigengenes and traits
  print("Selecting significant modules")
  module_trait_corr <- cor(module_eigengenes, traits, use = 'p')
  module_trait_corr_pvals <- corPvalueStudent(module_trait_corr, nSamples) %>% as.data.frame()
  
  # Identify significant modules with p-value < 0.05
  significant_modules <- module_trait_corr_pvals %>%
    filter(Disease < 0.05) %>%
    rownames() %>%
    gsub("^ME", "", .)
  
  print(paste0("Number of Significant Modules: ", length(significant_modules)))
  
  # Map genes to their respective modules
  module_gene_mapping <- as.data.frame(bwnet$colors, stringsAsFactors = FALSE)
  colnames(module_gene_mapping) <- "Module"
  
  # Extract genes belonging to significant modules
  SigGenesModule <- module_gene_mapping %>%
    filter(Module %in% significant_modules) %>%
    rownames() %>%
    unlist()
  
  print(paste0("Number of Significant Genes inside Modules: ", length(SigGenesModule)))
  
  combinedGenes <- SigGenesModule
  
  # ====================================
  # 5. Optional: Correlation with Trait, Was never run in this research
  # ====================================
  
  #Calculating Correlation of Genes with Trait
  if (includeCorrelation) {
    # Calculate correlation of each gene with the disease trait
    print("Selecting Significantly correlated Genes with Trait")
    gene_signf_corr <- cor(norm.counts, traits$Disease, use = 'p')
    gene_signf_corr_pvals <- corPvalueStudent(gene_signf_corr, nSamples)
    
    # Identify genes with p-value < 0.01
    SigGenesCor <- gene_signf_corr_pvals %>%
      as.data.frame() %>%
      arrange(V1) %>%
      filter(V1 < 0.01) %>% 
      rownames() %>% 
      unlist()
    
    print(paste0("Number of Significantly correlated Genes with Trait: ", length(SigGenesCor)))
    
    # Combine genes from significant modules and significantly correlated genes
    combinedGenes <- unique(c(SigGenesCor, SigGenesModule))
    
    print(paste0("Number of Combined Genes: ", length(combinedGenes)))
  }
  
  #Return WGCNA Genes
  return(combinedGenes)
}
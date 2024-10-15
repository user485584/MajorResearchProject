library(WGCNA)
library(magrittr)
library(dplyr)
library(DESeq2)

perform_wgcna <- function(counts, metadata, includeCorrelation = FALSE, power = NULL){
  
  t_counts <- t(counts)
  
  # Quality control step to check good samples and genes
  gsg <- goodSamplesGenes(t_counts, verbose = 3)
  print(summary(gsg))
  
  # Check all samples and genes are OK
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
  
  
  #Normalisation
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = metadata,
                                design = ~ 1)
  
  #Filtering so that each gene has at least count of 15 in at least 10 samples
  print("Filtering genes")
  dds75 <- dds[rowSums(counts(dds) >= 15) >= 10,]
  print(paste0("Number of remaining Genes: ",nrow(dds75)))
  
  #Variance stabilisation
  dds_norm <- vst(dds75)
  
  norm.counts <- assay(dds_norm) %>% t()
  
  
  # Choose Power for adjacency matrix or use the provided power
  if (is.null(power)) {
    powerOptions <- c(1:10, seq(from = 12, to = 30, by = 2))
    sft <- pickSoftThreshold(norm.counts,
                             powerVector = powerOptions,
                             networkType = "signed",
                             verbose = 5)
    softPower <- sft$powerEstimate
  } else {
    softPower <- power
  }
  
  temp_cor <- cor
  cor <- WGCNA::cor
  
  #Calculating modules
  bwnet <- blockwiseModules(norm.counts, maxBlockSize = 25000,
                            TOMType = 'signed', power = softPower,
                            mergeCutHeight = 0.25, numericLabels = FALSE,
                            randomSeed = 1234, verbose = 3)
  cor <- temp_cor
  
  module_eigengenes <- bwnet$MEs
  
  #preparation for follow-up steps
  traits <- metadata %>% 
    mutate(Disease = ifelse(grepl("ALS|AD", Status), 1, 0)) %>% 
    select(Disease)
  
  nSamples <- nrow(norm.counts)
  nGenes <- ncol(norm.counts)
  
  #Calculating significance of modules
  print("Selecting significant modules")
  module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
  module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples) %>% as.data.frame()
  
  significant_modules <- module.trait.corr.pvals %>%
    filter(Disease < 0.05) %>%
    rownames()
  significant_modules <- gsub("^ME", "", significant_modules)
  
  print(paste0("Number of Significant Modules: ", length(significant_modules)))
  
  #Filtering Significant modules and extracting Genes inside module
  module.gene.mapping <- as.data.frame(bwnet$colors, stringsAsFactors = FALSE)
  colnames(module.gene.mapping) <- "Module"
  
  SigGenesModule <- module.gene.mapping %>%
    filter(Module %in% significant_modules) %>%
    rownames() %>%
    unlist()
  
  print(paste0("Number of Significant Genes inside Modules: ", length(SigGenesModule)))
  
  combinedGenes <- SigGenesModule
  
  #Calculating Correlation of Genes with Trait
  if (includeCorrelation) {
    print("Selecting Significantly correlated Genes with Trait")
    gene.signf.corr <- cor(norm.counts, traits$Disease, use = 'p')
    gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)
    
    SigGenesCor <- gene.signf.corr.pvals %>%
      as.data.frame() %>%
      arrange(V1) %>%
      filter(V1 < 0.01) %>% 
      rownames() %>% 
      unlist()
    
    print(paste0("Number of Significantly correlated Genes with Trait: ", length(SigGenesCor)))
    
    # Combining Genes
    combinedGenes <- unique(c(SigGenesCor, SigGenesModule))
    
    print(paste0("Number of Combined Genes: ", length(combinedGenes)))
  }
  
  return(combinedGenes)
}
library(caret)
library(dplyr)
library(magrittr)
library(Seurat)
library(foreach)
library(doParallel)

# Load external scripts
source("/hpc/hers_en/rballieux/AlzheimersProofConcept/Classifiers/Scripts/CreatePseudobulk.R")
source("/hpc/hers_en/rballieux/AlzheimersProofConcept/Classifiers/Scripts/DE.R")
source("/hpc/hers_en/rballieux/AlzheimersProofConcept/Classifiers/Scripts/WGCNA.R")
source("/hpc/hers_en/rballieux/AlzheimersProofConcept/Classifiers/Scripts/PrepareDfClassifier.R")
source("/hpc/hers_en/rballieux/AlzheimersProofConcept/Classifiers/Scripts/TrainClassifier.R")
source("/hpc/hers_en/rballieux/AlzheimersProofConcept/Classifiers/Scripts/SummariseResults.R")
source("/hpc/hers_en/rballieux/AlzheimersProofConcept/Classifiers/Scripts/SaveObjects.R")

# Register parallel backend to use multiple cores
required_packages <- c("Seurat", "dplyr", "caret", "magrittr", "DESeq2", "tibble", "WGCNA", "ROCR", "glmnet", "randomForest", "GRaNIE", "GRaNPA")
no_cores <- 6 - 1  # Use one less than the total number of cores
registerDoParallel(cores=no_cores)

# Load and convert Data
SeuratObject <- readRDS("/hpc/hers_en/rballieux/ALS_FTLD/Data/SALS_Seurat_Downsampled.rds")
Disease <- "SALS"

SeuratCounts <- CreateSeuratObject(counts = SeuratObject@assays$RNA@layers$counts) %>%
  AddMetaData(metadata = SeuratObject@meta.data %>% select(id, subs, predicted.id, Status))
colnames(SeuratCounts) <- colnames(SeuratObject)
rownames(SeuratCounts) <- rownames(SeuratObject)

metadata <- SeuratObject@meta.data %>% 
  tibble::rownames_to_column("cellbarcode") %>%
  select(cellbarcode, id, subs, predicted.id, Status) 

# Create stratified folds for Cross Validation
folds <- createFolds(metadata$id, k = 5, list = TRUE, returnTrain = TRUE)

save_objects(folds, fold_number = "all")

# Initialize an empty dataframe to store results
results_df <- foreach(i = 1:5, .combine = rbind, .packages = required_packages, 
                      .export = c("SeuratObject", "SeuratCounts", "metadata")) %dopar% {
  print(paste0("Processing fold ", i))
  tryCatch({
    # Create Train and Test Subsets for DE, WGCNA, and GRN
    train_barcodes <- metadata$cellbarcode[folds[[i]]]
    TrainSubset <- SeuratCounts[, train_barcodes]
    
    test_barcodes <- metadata$cellbarcode[-folds[[i]]]
    TestSubset <- SeuratCounts[, test_barcodes]
    Testmetadata <- metadata[metadata$cellbarcode %in% test_barcodes,]
    
    # Create pseudobulk counts and metadata for DE, WGCNA, and GRN analysis
    pseudobulk <- CreatePseudobulk(TrainSubset, metadata, type = "standard")
    pseudobulk_df <- pseudobulk$pseudobulk_df
    pseudobulk_metadata <- pseudobulk$pseudobulk_metadata
    
    
    # Perform WGCNA and DE analysis
    de_features <- perform_de(pseudobulk_df, pseudobulk_metadata)
    wgcna_features <- perform_wgcna(pseudobulk_df, pseudobulk_metadata, disease = "ALS", includeCorrelation = FALSE, power = NULL)
    
    # Combine and deduplicate feature lists
    combined_features <- unique(c(wgcna_features, de_features))
    
    
    save_objects(train_barcodes, test_barcodes, combined_features, fold_number = i)
    
    max_features <- length(combined_features)
  
    
    # Create data frames for Training and Testing Classifier
    TrainData_DEWGCNA <- prepare_df_classifier(TrainSubset, combined_features)
    TestData_DEWGCNA <- prepare_df_classifier(TestSubset, combined_features)
    
    save_objects(TrainData_DEWGCNA, TestData_DEWGCNA, fold_number = i, combined_features)
    
    
    # Training Classifier with DE/WGCNA Features
    fold_results_DEWGCNA <- train_classifier(TrainData_DEWGCNA, TestData_DEWGCNA, Testmetadata, "all", PerCelltype = TRUE, foldNr = i, Analysis = "DEWGCNA")
    fold_results_DEWGCNA$Disease <- Disease
    fold_results_DEWGCNA$FoldNumber <- i
    fold_results_DEWGCNA$NumFeatures <- length(combined_features)
    fold_results_DEWGCNA$Analysis <- "DE/WGCNA"
    
    save_objects(fold_results_DEWGCNA, fold_number = i)
    
    
    fold_results <- fold_results_DEWGCNA %>%
      dplyr::select(Disease, FoldNumber, Analysis, NumFeatures, Classifier, everything())
    
    fold_results  # Return this dataframe for each iteration
    
    save_objects(fold_results, fold_number = i)
    
  }, error = function(e) {
    print(paste0("Error in fold ", i, ": ", e))
    NULL  # Return NULL for this fold in case of error
  })
}

save_objects(results_df, fold_number = "all")


write.csv(results_df, "/hpc/hers_en/rballieux/AlzheimersProofConcept/Classifiers/RunFiles/results_SALS.csv", row.names = FALSE)
print(results_df)
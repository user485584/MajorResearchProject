# Load necessary libraries
library(caret)
library(dplyr)
library(magrittr)
library(Seurat)
library(foreach)
library(doParallel)


# ============================
# Configuration Section
# ============================
# Update these paths as needed to configure the script
config <- list(
  Disease = "C9ALS",
  PathToScriptsFolder = "/hpc/hers_en/rballieux/AlzheimersProofConcept/Classifiers/Scripts/", # Folder containing helper scripts
  PathToDataset = "/hpc/hers_en/rballieux/ALS_FTLD/Data/C9ALS_Seurat_Downsampled.rds",           # Path to the Seurat dataset
  OutputFolder = "/hpc/hers_en/rballieux/ALS_FTLD/C9ALS/DEWGCNA/"        # Base folder for all outputs
)


# ============================
# Sourcing External Scripts
# ============================
# Source the necessary R scripts using the PathToScriptsFolder from config
source(paste0(config$PathToScriptsFolder, "CreatePseudobulk.R"))
source(paste0(config$PathToScriptsFolder, "DE.R"))
source(paste0(config$PathToScriptsFolder, "WGCNA.R"))
source(paste0(config$PathToScriptsFolder, "PrepareDfClassifier.R"))
source(paste0(config$PathToScriptsFolder, "TrainClassifier.R"))
source(paste0(config$PathToScriptsFolder, "SaveObjects.R"))


# ============================
# Parallel Backend Setup
# ============================
# Register parallel backend to utilize multiple cores for processing
required_packages <- c("Seurat", "dplyr", "caret", "magrittr", "DESeq2", "tibble", "WGCNA", "ROCR", "glmnet", "randomForest", "GRaNIE", "GRaNPA")
#no_cores <- detectCores() - 1  # Use one less than the total number of cores
no_cores <- 5                   # Hardcode nr of codes
registerDoParallel(cores=no_cores)


# ============================
# Loading and Preparing Data
# ============================
# Load the Seurat object from the specified dataset path
SeuratObject <- readRDS(config$PathToDataset)

# Create a Seurat object with selected metadata
SeuratCounts <- CreateSeuratObject(counts = SeuratObject@assays$RNA@layers$counts) %>%
  AddMetaData(metadata = SeuratObject@meta.data %>% select(id, subs, predicted.id, Status))
colnames(SeuratCounts) <- colnames(SeuratObject)
rownames(SeuratCounts) <- rownames(SeuratObject)

# Prepare metadata with cell barcodes
metadata <- SeuratObject@meta.data %>% 
  tibble::rownames_to_column("cellbarcode") %>%
  select(cellbarcode, id, subs, predicted.id, Status) 


# ============================
# Creating Cross-Validation Folds
# ============================
# Create stratified folds for cross-validation based on donor IDs
folds <- createFolds(metadata$id, k = 5, list = TRUE, returnTrain = TRUE)

# Save the folds for future reference
save_objects(folds, fold_number = "all", path_prefix = config$OutputFolder)


# ============================
# Feature Processing and Model Training
# ============================
# Initialize an empty dataframe to store results
results_df <- foreach(i = 1:5, .combine = rbind, .packages = required_packages, 
                      .export = c("SeuratObject", "SeuratCounts", "metadata")) %dopar% {
  print(paste0("Processing fold ", i))
  tryCatch({
    # Subset training and testing barcodes for the current fold
    train_barcodes <- metadata$cellbarcode[folds[[i]]]
    TrainSubset <- SeuratCounts[, train_barcodes]
    
    test_barcodes <- metadata$cellbarcode[-folds[[i]]]
    TestSubset <- SeuratCounts[, test_barcodes]
    Testmetadata <- metadata[metadata$cellbarcode %in% test_barcodes,]
    
    # Create pseudobulk counts and metadata for DE and WGCNA analyses
    pseudobulk <- CreatePseudobulk(TrainSubset, metadata)
    pseudobulk_df <- pseudobulk$pseudobulk_df
    pseudobulk_metadata <- pseudobulk$pseudobulk_metadata
    
    # Perform Differential Expression (DE) and Weighted Gene Co-expression Network Analysis (WGCNA)
    de_features <- perform_de(pseudobulk_df, pseudobulk_metadata)
    wgcna_features <- perform_wgcna(pseudobulk_df, pseudobulk_metadata, includeCorrelation = FALSE, power = NULL)
    
    # Combine and deduplicate feature lists
    combined_features <- unique(c(wgcna_features, de_features))
    
    # Save intermediate results
    save_objects(train_barcodes, test_barcodes, de_features, wgcna_features, combined_features, fold_number = i, path_prefix = config$OutputFolder)
  
    # Prepare data frames for training and testing classifiers using DE/WGCNA features
    TrainData_DEWGCNA <- prepare_df_classifier(TrainSubset, combined_features)
    TestData_DEWGCNA <- prepare_df_classifier(TestSubset, combined_features)
    
    # Define the analysis identifier
    Analysis <- "DE/WGCNA"
    
    # Train classifiers and obtain results
    results <- train_classifier(TrainData_DEWGCNA, TestData_DEWGCNA, Testmetadata, "all", PerCelltype = TRUE, foldNr = i, Analysis = Analysis,
                                preloaded_model = NULL, ROCoutput_dir = paste0(config$OutputFolder, "ROCplots/"))
    
    results$Disease <- Disease
    results$FoldNumber <- i
    results$NumFeatures <- length(combined_features)
    results$Analysis <- Analysis
    
    results <- select(results, Disease, FoldNumber, Analysis, NumFeatures, Classifier, CellType, NrOfCells, everything())
    
    save_objects(results, fold_number = i, path_prefix = config$OutputFolder)
    
    
    #Display the results for the current fold for logfile
    results  
  }, error = function(e) {
    print(paste0("Error in fold ", i, ": ", e))
    NULL  # Return NULL for this fold in case of error
  })
}

# ============================
# Combining and Saving All Results
# ============================
# Save the combined results to the specified output folder
save_objects(results_df, fold_number = "all", path_prefix = config$OutputFolder)

# Export all accumulated results to a CSV file in the designated output folder
write.csv(results_df, paste0(config$OutputFolder, Disease, "DEWGCNAresults.csv"), row.names = FALSE)

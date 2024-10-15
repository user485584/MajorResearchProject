# Load necessary libraries
library(caret)
library(dplyr)
library(magrittr)
library(Seurat)


# ============================
# Configuration Section
# ============================
# Update these paths as needed to configure the script
config <- list(
  Disease = "C9ALS",
  PathToScriptsFolder = "/hpc/hers_en/rballieux/AlzheimersProofConcept/Classifiers/Scripts/", # Folder containing helper scripts
  PathToDataset = "/hpc/hers_en/rballieux/ALS_FTLD/Data/C9ALS_Seurat_Downsampled.rds",  # Path to the Seurat dataset
  chunk_size = 3000,                                                                     # Number of features to process per chunk
  OutputFolder = "/hpc/hers_en/rballieux/ALS_FTLD/C9ALS/3000randomfeatures/"           # Base folder for all outputs
)


# ============================
# Sourcing External Scripts
# ============================
# Source the necessary R scripts using the PathToScriptsFolder from config
source(paste0(config$PathToScriptsFolder, "PrepareDfClassifier.R"))
source(paste0(config$PathToScriptsFolder, "TrainClassifier.R"))
source(paste0(config$PathToScriptsFolder, "SaveObjects.R"))


# ============================
# Loading and Preparing Data
# ============================
# Load the Seurat object from the specified dataset path
SeuratObject <- readRDS(config$PathToDataset)

# Initialize the starting index for feature processing
start <- 1

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
# Create cross-validation folds based on donor IDs
folds <- createFolds(metadata$id, k = 5, list = TRUE, returnTrain = TRUE)
i <- 1

# Extract training barcodes and subset the Seurat object for training
train_barcodes <- metadata$cellbarcode[folds[[i]]]
TrainSubset <- SeuratCounts[, train_barcodes]

# Extract testing barcodes and subset the Seurat object for testing
test_barcodes <- metadata$cellbarcode[-folds[[i]]]
TestSubset <- SeuratCounts[, test_barcodes]
Testmetadata <- metadata[metadata$cellbarcode %in% test_barcodes,]


# ============================
# Feature Processing and Model Training
# ============================
# Retrieve all gene features
all_features <- rownames(SeuratCounts)
total_features <- length(all_features)

# Initialize an empty data frame to store all results
all_results <- data.frame()

for (start_index in seq(start, total_features, by = config$chunk_size)) {
  # Determine the end index for the current chunk
  end_index <- min(start_index + config$chunk_size - 1, total_features)
  
  # Subset the features for the current chunk
  feature_subset <- all_features[start_index:end_index]
  
  #Info for Logfile
  cat(paste0("Processing features ", start_index, " to ", end_index, "\n"))
  
  # Prepare training and testing data using the feature subset
  TrainData <- prepare_df_classifier(TrainSubset, feature_subset)
  TestData <- prepare_df_classifier(TestSubset, feature_subset)
  
  # Define the analysis identifier based on feature range
  Analysis <- paste0("Features_", start_index, "_to_", end_index)
  
  # Train classifiers and obtain results
  results <- train_classifier(TrainData, TestData, Testmetadata, "all", PerCelltype = TRUE, foldNr = i, Analysis = Analysis, 
                              preloaded_models = NULL, ROCoutput_dir = paste0(config$OutputFolder, "ROCplots/"))
  
  results$Disease <- Disease
  results$FoldNumber <- paste0(start_index, "-", end_index)
  results$NumFeatures <- config$chunk_size
  results$Analysis <- "3000RandomFeatures"
  
  results <- select(results, Disease, FoldNumber, Analysis, NumFeatures, Classifier, CellType, NrOfCells, everything())
  
  
  # Display the results for the current chunk for logfile
  print(results)
  
  # Save intermediate results to the specified output directory
  save_objects(train_barcodes, test_barcodes, results, fold_number = start_index, path_prefix = config$OutputFolder)
  
  # Append the current results to the cumulative all_results data frame
  all_results <- rbind(all_results, results)
}

# ============================
# Exporting Final Results
# ============================
# Save the combined results to the specified output folder
save_objects(all_results, fold_number = "all", path_prefix = config$OutputFolder)

# Export all accumulated results to a CSV file in the designated output folder
write.csv(all_results, paste0(config$OutputFolder, Disease, "ClassifierResults3000randomfeats.csv"), row.names = FALSE)
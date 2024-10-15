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
  OutputFolder = "/hpc/hers_en/rballieux/ALS_FTLD/C9ALS/MostVariableFeatures/" # Base folder for all outputs
)

# ============================
# Sourcing External Scripts
# ============================
# Source the necessary R scripts using the PathToScriptsFolder from config
source(paste0(config$PathToScriptsFolder, "PrepareDfClassifier.R"))
source(paste0(config$PathToScriptsFolder, "TrainClassifier.R"))
source(paste0(config$PathToScriptsFolder, "SaveObjects.R"))


# ============================
# Parallel Backend Setup
# ============================
# Register parallel backend to utilise multiple cores for processing
required_packages <- c("Seurat", "dplyr", "caret", "magrittr", "foreach", "doParallel")
#no_cores <- detectCores() - 1  # Use one less than the total number of cores
no_cores <- 5                   # Hardcode nr of codes
registerDoParallel(cores=no_cores)


# ============================
# Loading and Preparing Data
# ============================
# Load the Seurat object from the specified dataset path
# Load the Seurat object
SeuratObject <- readRDS(config$PathToDataset)

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

# Initialize an empty dataframe to store results
results_df <- foreach(i = 1:5, .combine = rbind, .packages = required_packages, 
                      .export = c("SeuratObject", "SeuratCounts", "metadata")) %dopar% {
                        print(paste0("Processing fold ", i))
                        tryCatch({
                          # Extract training and testing barcodes for the current fold
                          train_barcodes <- metadata$cellbarcode[folds[[i]]]
                          TrainSubset <- SeuratCounts[, train_barcodes]
                          
                          test_barcodes <- metadata$cellbarcode[-folds[[i]]]
                          TestSubset <- SeuratCounts[, test_barcodes]
                          Testmetadata <- metadata[metadata$cellbarcode %in% test_barcodes,]
                          
                          # Calculate variance for each gene in the training subset
                          gene_variances <- apply(TrainSubset@assays$RNA@layers$counts, 1, var)
                          names(gene_variances) <- rownames(TrainSubset)
                          
                          variance_df <- data.frame(
                            gene = names(gene_variances), 
                            variance = gene_variances,  
                            stringsAsFactors = FALSE  
                          ) %>% arrange(desc(variance))
                        
                          # Save the variance dataframe and barcodes to the specified output folder
                          save_objects(train_barcodes, test_barcodes, variance_df, fold_number = i, path_prefix = config$OutputFolder)
                          
                          # Define the numbers of top variable features to consider
                          top_feature_counts <- c(100, 500, 1000, 1500, 2000)
                          
                          # Initialise an empty dataframe to store all results for the current fold
                          fold_results <- data.frame()
                          
                          # Iterate over each specified number of top features
                          for (top_features in top_feature_counts) {
                            cat(paste0("Selecting top ", top_features, " variable features\n"))
                            
                            # Select the top variable features
                            top_variable_features <- variance_df %>%
                              head(top_features) %>%
                              pull(gene)
                            
                            # Prepare training and testing data using the selected features
                            TrainData <- prepare_df_classifier(TrainSubset, top_variable_features)
                            TestData <- prepare_df_classifier(TestSubset, top_variable_features)
                            
                            # Define the analysis identifier based on the number of features
                            Analysis <- paste0("Top", top_features, "VariableFeatures")
                            
                            # Train classifiers and obtain results
                            results <- train_classifier(TrainData, TestData, Testmetadata, "all", PerCelltype = TRUE, foldNr = i, Analysis = Analysis, 
                                                        preloaded_model = NULL, ROCoutput_dir = paste0(config$OutputFolder, "ROCplots/"))
                            
                            results$Disease <- Disease
                            results$FoldNumber <- i
                            results$NumFeatures <- top_features
                            results$Analysis <- "MostVariableFeatures"
                            
                            results <- select(results, Disease, FoldNumber, Analysis, NumFeatures, Classifier, CellType, NrOfCells, everything())
                           
                            
                            # Display the results for the current feature set for logfile
                            print(results)
                            
                            # Save the results to the specified output folder
                            save_objects(results, fold_number = paste0(i,"_",top_features,"VarFeats"), path_prefix = config$OutputFolder)
                            
                            # Append the current results to the fold_results dataframe
                            fold_results <- rbind(fold_results, results)
                          }
                          
                          # Save the final results for the current fold
                          save_objects(fold_results, fold_number = i, path_prefix = config$OutputFolder)
                          
                          # Display the results for the current fold for logfile
                          fold_results 
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
write.csv(results_df, paste0(config$OutputFolder, Disease, "VariableFeaturesResults.csv"), row.names = FALSE)

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
  OutputFolder = "/hpc/hers_en/rballieux/ALS_FTLD/C9ALS/TopPCs/"                               # Base folder for all outputs
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
SeuratObject <- readRDS(config$PathToDataset)

# Create a Seurat object with selected metadata
SeuratCounts <- CreateSeuratObject(counts = SeuratObject@assays$RNA@layers$counts) %>%
  AddMetaData(metadata = SeuratObject@meta.data %>% select(id, subs, predicted.id, Status))
colnames(SeuratCounts) <- colnames(SeuratObject)
rownames(SeuratCounts) <- rownames(SeuratObject)

# Prepare metadata with cell barcodes
metadata <- SeuratCounts@meta.data %>% 
  tibble::rownames_to_column("cellbarcode") %>%
  select(cellbarcode, id, subs, predicted.id, Status) 


# ============================
# Creating Cross-Validation Folds
# ============================
# Create cross-validation folds based on donor IDs
folds <- createFolds(metadata$id, k = 5, list = TRUE, returnTrain = TRUE)


# ============================
# Feature Processing and Model Training
# ============================
# Initialize an empty dataframe to store results
results_df <- foreach(i = 1:5, .combine = rbind, .packages = required_packages, 
                      .export = c("SeuratObject", "SeuratCounts", "metadata")) %dopar% {
                        print(paste0("Processing fold ", i))
                        tryCatch({
                          # Subset train and test data
                          train_barcodes <- metadata$cellbarcode[folds[[i]]]
                          TrainSubset <- SeuratCounts[, train_barcodes]
                          
                          test_barcodes <- metadata$cellbarcode[-folds[[i]]]
                          TestSubset <- SeuratCounts[, test_barcodes]
                          Testmetadata <- metadata[metadata$cellbarcode %in% test_barcodes,]
                          
                          # Identify non-zero genes in the training subset
                          non_zero_genes <- rownames(TrainSubset)[rowSums(TrainSubset@assays$RNA@layers$counts) != 0]
                          TrainSubset <- subset(TrainSubset, features = non_zero_genes)
                          TestSubset <- subset(TestSubset, features = non_zero_genes)
                          
                          # Identify variable features using Seurat's FindVariableFeatures
                          TrainSubset <- FindVariableFeatures(TrainSubset, selection.method = "vst", nfeatures = 2000)
                          
                          # Normalise the data
                          TrainSubset <- NormalizeData(TrainSubset)
                          TestSubset <- NormalizeData(TestSubset)
                          
                          # Scale the data based on variable features
                          TrainSubset <- ScaleData(TrainSubset, features = VariableFeatures(object = TrainSubset))
                          TestSubset <- ScaleData(TestSubset, features = VariableFeatures(object = TrainSubset))
                          
                          # Run PCA on the training data
                          TrainSubset <- RunPCA(TrainSubset, npcs = 1000, verbose = FALSE, approx = FALSE)
                          
                          # Extract PCA embeddings and rotation matrix for the training data
                          pca_train <- Embeddings(TrainSubset, reduction = "pca")
                          rotation_matrix <- Loadings(TrainSubset, "pca")
                          
                          variable_features <- VariableFeatures(TrainSubset)
                          
                          # Project test data using the same PCA model
                          test_scaled_data <- GetAssayData(TestSubset, slot = "scale.data")[variable_features, ]
                          pca_test <- t(as.matrix(test_scaled_data)) %*% as.matrix(rotation_matrix)
                          
                          # Save the barcodes to the specified output folder
                          save_objects(train_barcodes, test_barcodes, fold_number = i, path_prefix = config$OutputFolder)
                          
                          # Define the numbers of top PCs to consider
                          top_pc_counts <- c(10, 50, 100, 250, 500, 1000)
                          
                          # Initialise an empty dataframe to store all results for the current fold
                          fold_results <- data.frame()
                          
                          # Iterate over each specified number of top PCs
                          for (top_pcs in top_pc_counts) {
                            cat(paste0("Selecting top ", top_pcs, " principal components\n"))
                            
                            # Select the top PCs for training data
                            pca_features_train <- pca_train[, 1:top_pcs]
                            colnames(pca_features_train) <- paste0("PC", 1:top_pcs)
                            
                            # Calculate cumulative variance explained
                            explained_variance <- sum(TrainSubset@reductions$pca@stdev[1:top_pcs]^2)
                            total_variance <- sum(TrainSubset@reductions$pca@stdev^2)
                            cumulative_variance_explained <- explained_variance / total_variance
                            
                            cat(paste0("Cumulative variance explained by the top ", top_pcs, " PCs: ",
                                       round(cumulative_variance_explained * 100, 2), "%\n"))
                            
                            # Prepare Test PCA data
                            pca_features_test <- pca_test[, 1:top_pcs]
                            colnames(pca_features_test) <- paste0("PC", 1:top_pcs)
                            
                            # Create data frames for training and testing classifiers
                            TrainData <- as.data.frame(pca_features_train)
                            TrainData$Status <- TrainSubset@meta.data$Status
                            TrainData$Status <- ifelse(TrainData$Status %in% c("ALS", "AD"), 1, 0)
                            rownames(TrainData) <- colnames(TrainSubset)
                            
                            TestData <- as.data.frame(pca_features_test)
                            TestData$Status <- TestSubset@meta.data$Status
                            TestData$Status <- ifelse(TestData$Status %in% c("ALS", "AD"), 1, 0)
                            rownames(TestData) <- colnames(TestSubset)
                            
                            # Define the analysis identifier based on the number of PCs
                            Analysis <- paste0("Top", top_pcs, "PCs")
                            
                            # Train classifiers and obtain results
                            results <- train_classifier(TrainData, TestData, Testmetadata, "all", PerCelltype = TRUE, foldNr = i, Analysis = Analysis,
                                                        preloaded_model = NULL, ROCoutput_dir = paste0(config$OutputFolder, "ROCplots/"))
                            
                            results$Disease <- Disease
                            results$FoldNumber <- i
                            results$NumFeatures <- top_pcs
                            results$Analysis <- "TopPCs"
                          
                            results <- select(results, Disease, FoldNumber, Analysis, NumFeatures, Classifier, CellType, NrOfCells, everything())
                            
                            # Display the results for the current feature set for logfile
                            print(results)
                            
                            # Save the results to the specified output folder
                            save_objects(results, fold_number = paste0(i,"_Top",top_pcs,"PCs"), path_prefix = config$OutputFolder)
                            
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
write.csv(results_df, paste0(config$OutputFolder, Disease, "TopPCsResults.csv"), row.names = FALSE)

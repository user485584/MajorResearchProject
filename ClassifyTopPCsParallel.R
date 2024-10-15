library(caret)
library(dplyr)
library(magrittr)
library(Seurat)
library(foreach)
library(doParallel)

# Load external scripts
source("/hpc/hers_en/rballieux/AlzheimersProofConcept/Classifiers/Scripts/PrepareDfClassifier.R")
source("/hpc/hers_en/rballieux/AlzheimersProofConcept/Classifiers/Scripts/TrainClassifier.R")
source("/hpc/hers_en/rballieux/AlzheimersProofConcept/Classifiers/Scripts/SaveObjects.R")

# Register parallel backend to use multiple cores
required_packages <- c("Seurat", "dplyr", "caret", "magrittr", "foreach", "doParallel")
no_cores <- 6 - 1  # Use one less than the total number of cores
registerDoParallel(cores=no_cores)

# Load the Seurat object
SeuratObject <- readRDS("/hpc/hers_en/rballieux/ALS_FTLD/Data/SALS_Seurat_Downsampled.rds")

# Create Seurat counts object
SeuratCounts <- CreateSeuratObject(counts = SeuratObject@assays$RNA@layers$counts) %>%
  AddMetaData(metadata = SeuratObject@meta.data %>% select(id, subs, predicted.id, Status))
colnames(SeuratCounts) <- colnames(SeuratObject)
rownames(SeuratCounts) <- rownames(SeuratObject)

# Extract metadata
metadata <- SeuratCounts@meta.data %>% 
  tibble::rownames_to_column("cellbarcode") %>%
  select(cellbarcode, id, subs, predicted.id, Status) 

# Create folds for cross-validation
folds <- createFolds(metadata$id, k = 5, list = TRUE, returnTrain = TRUE)

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
                          
                          non_zero_genes <- rownames(TrainSubset)[rowSums(TrainSubset@assays$RNA@layers$counts) != 0]
                          TrainSubset <- subset(TrainSubset, features = non_zero_genes)
                          TestSubset <- subset(TestSubset, features = non_zero_genes)
                          
                          # Identify variable features
                          TrainSubset <- FindVariableFeatures(TrainSubset, selection.method = "vst", nfeatures = 2000)
                          
                          # Normalize the data
                          TrainSubset <- NormalizeData(TrainSubset)
                          TestSubset <- NormalizeData(TestSubset)
                          
                          # Run PCA using Seurat's RunPCA
                          TrainSubset <- ScaleData(TrainSubset, features = VariableFeatures(object = TrainSubset))
                          TestSubset <- ScaleData(TestSubset, features = VariableFeatures(object = TrainSubset))
                          
                          TrainSubset <- RunPCA(TrainSubset, npcs = 1000, verbose = FALSE, approx = FALSE)
                          
                          # Extract PCA embeddings for training data
                          pca_train <- Embeddings(TrainSubset, reduction = "pca")
                          rotation_matrix <- Loadings(TrainSubset, "pca")
                          
                          variable_features <- VariableFeatures(TrainSubset)
                          
                          # Project test data using the same PCA model
                          test_scaled_data <- GetAssayData(TestSubset, slot = "scale.data")[variable_features, ]
                          
                          pca_test <- t(as.matrix(test_scaled_data)) %*% as.matrix(rotation_matrix)
                          
                          # Define the numbers of top PCs to consider
                          top_pc_counts <- c(10, 50, 100, 250, 500, 1000)
                          
                          # Initialize an empty dataframe to store all results for the current fold
                          fold_results <- data.frame()
                          
                          for (top_pcs in top_pc_counts) {
                            cat(paste0("Selecting top ", top_pcs, " principal components\n"))
                            
                            # Select the top PCs
                            pca_features_train <- pca_train[, 1:top_pcs]
                            colnames(pca_features_train) <- paste0("PC", 1:top_pcs)
                            
                            explained_variance <- sum(TrainSubset@reductions$pca@stdev[1:top_pcs]^2)
                            total_variance <- sum(TrainSubset@reductions$pca@stdev^2)
                            cumulative_variance_explained <- explained_variance / total_variance
                            
                            cat(paste0("Cumulative variance explained by the top ", top_pcs, " PCs: ",
                                       round(cumulative_variance_explained * 100, 2), "%\n"))
                            
                            # Prepare Test PCA data
                            pca_features_test <- pca_test[, 1:top_pcs]
                            colnames(pca_features_test) <- paste0("PC", 1:top_pcs)
                            
                            TrainData <- as.data.frame(pca_features_train)
                            TrainData$Status <- TrainSubset@meta.data$Status
                            TrainData$Status <- ifelse(TrainData$Status == "ALS", 1, 0)
                            rownames(TrainData) <- colnames(TrainSubset)
                            
                            TestData <- as.data.frame(pca_features_test)
                            TestData$Status <- TestSubset@meta.data$Status
                            TestData$Status <- ifelse(TestData$Status == "ALS", 1, 0)
                            rownames(TestData) <- colnames(TestSubset)
                            
                            # Train classifier
                            Analysis <- paste0("Top", top_pcs, "PCs")
                            results <- train_classifier(TrainData, TestData, Testmetadata, "all", PerCelltype = TRUE, foldNr = i, Analysis = Analysis)
                            
                            # Add feature count to results
                            results$NumFeatures <- top_pcs
                            
                            print(results)
                            
                            save_objects(results, fold_number = paste0(i,"_Top",top_pcs,"PCs"), path_prefix = "/hpc/hers_en/rballieux/ALS_FTLD/SALS/TopPCs/")
                            
                            # Append results to fold_results
                            fold_results <- rbind(fold_results, results)
                          }
                          
                          # Save final fold results
                          save_objects(fold_results, fold_number = i, path_prefix = "/hpc/hers_en/rballieux/ALS_FTLD/SALS/TopPCs/")
                          
                          fold_results  # Return this dataframe for each iteration
                        }, error = function(e) {
                          print(paste0("Error in fold ", i, ": ", e))
                          NULL  # Return NULL for this fold in case of error
                        })
                      }

# Combine results from all folds
all_results <- results_df %>%
  dplyr::select(Classifier, NumFeatures, everything())

save_objects(all_results, fold_number = "all", path_prefix = "/hpc/hers_en/rballieux/ALS_FTLD/SALS/TopPCs/")

write.csv(all_results, "/hpc/hers_en/rballieux/ALS_FTLD/SALS/TopPCs/SALSTopPCsResults.csv", row.names = FALSE)

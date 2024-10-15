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

SeuratCounts <- CreateSeuratObject(counts = SeuratObject@assays$RNA@layers$counts) %>%
  AddMetaData(metadata = SeuratObject@meta.data %>% select(id, subs, predicted.id, Status))
colnames(SeuratCounts) <- colnames(SeuratObject)
rownames(SeuratCounts) <- rownames(SeuratObject)

metadata <- SeuratObject@meta.data %>% 
  tibble::rownames_to_column("cellbarcode") %>%
  select(cellbarcode, id, subs, predicted.id, Status)

# Create folds for cross-validation
folds <- createFolds(metadata$id, k = 5, list = TRUE, returnTrain = TRUE)

# Initialize an empty dataframe to store results
results_df <- foreach(i = 1:5, .combine = rbind, .packages = required_packages, 
                      .export = c("SeuratObject", "SeuratCounts", "metadata")) %dopar% {
                        print(paste0("Processing fold ", i))
                        tryCatch({
                          train_barcodes <- metadata$cellbarcode[folds[[i]]]
                          TrainSubset <- SeuratCounts[, train_barcodes]
                          
                          test_barcodes <- metadata$cellbarcode[-folds[[i]]]
                          TestSubset <- SeuratCounts[, test_barcodes]
                          Testmetadata <- metadata[metadata$cellbarcode %in% test_barcodes,]
                          
                          # Calculate variance for each gene
                          gene_variances <- apply(TrainSubset@assays$RNA@layers$counts, 1, var)
                          names(gene_variances) <- rownames(TrainSubset)
                          
                          variance_df <- data.frame(
                            gene = names(gene_variances),  # Gene names
                            variance = gene_variances,  # Variances
                            stringsAsFactors = FALSE  # To prevent automatic conversion to factors
                          ) %>% arrange(desc(variance))
                          
                          variance_df_file <- paste0("/hpc/hers_en/rballieux/ALS_FTLD/C9ALS/MostVariableFeatures/Parallel23-7/variance_df_Fold",i,".rds")
                          variance_df <- readRDS(variance_df_file)
                          
                          save_objects(variance_df, fold_number = i, path_prefix = "/hpc/hers_en/rballieux/ALS_FTLD/C9ALS_SALS_FeatureSwitch/MostVariableFeatures/SALS_C9ALSfeatures/")
                          
                          # Define the numbers of top variable features to consider
                          top_feature_counts <- c(100, 500, 1000, 1500, 2000)
                          # Initialize an empty dataframe to store all results for the current fold
                          fold_results <- data.frame()
                          
                          for (top_features in top_feature_counts) {
                            cat(paste0("Selecting top ", top_features, " variable features\n"))
                            
                            # Select the top variable features
                            top_variable_features <- variance_df %>%
                              head(top_features) %>%
                              pull(gene)
                            
                            # Prepare data for classifier
                            TrainData <- prepare_df_classifier(TrainSubset, top_variable_features)
                            TestData <- prepare_df_classifier(TestSubset, top_variable_features)
                            
                            # Train classifier
                            Analysis = paste0("Top", top_features, "VariableFeatures")
                            results <- train_classifier(TrainData, TestData, Testmetadata, "all", PerCelltype = TRUE, foldNr = i, Analysis = Analysis, preloaded_model = NULL)
                            
                            # Add feature count to results
                            results$NumFeatures <- top_features
                            
                            print(results)
                            
                            save_objects(results, fold_number = paste0(i,"_",top_features,"VarFeats"), path_prefix = "/hpc/hers_en/rballieux/ALS_FTLD/C9ALS_SALS_FeatureSwitch/MostVariableFeatures/SALS_C9ALSfeatures/")
                            
                            # Append results to fold_results
                            fold_results <- rbind(fold_results, results)
                          }
                          
                          # Save final fold results
                          save_objects(fold_results, fold_number = i, path_prefix = "/hpc/hers_en/rballieux/ALS_FTLD/C9ALS_SALS_FeatureSwitch/MostVariableFeatures/SALS_C9ALSfeatures/")
                          
                          fold_results  # Return this dataframe for each iteration
                        }, error = function(e) {
                          print(paste0("Error in fold ", i, ": ", e))
                          NULL  # Return NULL for this fold in case of error
                        })
                      }

# Combine results from all folds
all_results <- results_df %>%
  dplyr::select(Classifier, NumFeatures, everything())

save_objects(all_results, fold_number = "all", path_prefix = "/hpc/hers_en/rballieux/ALS_FTLD/C9ALS_SALS_FeatureSwitch/MostVariableFeatures/SALS/")

write.csv(all_results, "/hpc/hers_en/rballieux/ALS_FTLD/C9ALS_SALS_FeatureSwitch/MostVariableFeatures/SALS_C9ALSfeatures/VariableFeaturesResultsSALS.csv", row.names = FALSE)

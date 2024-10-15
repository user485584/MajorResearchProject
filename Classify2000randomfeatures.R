library(caret)
library(dplyr)
library(magrittr)
library(Seurat)


source("/hpc/hers_en/rballieux/AlzheimersProofConcept/Classifiers/Scripts/PrepareDfClassifier.R")
source("/hpc/hers_en/rballieux/AlzheimersProofConcept/Classifiers/Scripts/TrainClassifierTest2.R")
source("/hpc/hers_en/rballieux/AlzheimersProofConcept/Classifiers/Scripts/SaveObjects.R")


SeuratObject <- readRDS("/hpc/hers_en/rballieux/AlzheimersProofConcept/RObjects/SeuratObject.rds")

SeuratCounts <- CreateSeuratObject(counts = SeuratObject@assays$RNA@layers$counts) %>%
  AddMetaData(metadata = SeuratObject@meta.data %>% select(id, subs, predicted.id, Status))
colnames(SeuratCounts) <- colnames(SeuratObject)
rownames(SeuratCounts) <- rownames(SeuratObject)

metadata <- SeuratObject@meta.data %>% 
  tibble::rownames_to_column("cellbarcode") %>%
  select(cellbarcode, id, subs, predicted.id, Status) 


folds <- createFolds(metadata$id, k = 5, list = TRUE, returnTrain = TRUE)
i <- 1

train_barcodes <- metadata$cellbarcode[folds[[i]]]
TrainSubset <- SeuratCounts[, train_barcodes]

test_barcodes <- metadata$cellbarcode[-folds[[i]]]
TestSubset <- SeuratCounts[, test_barcodes]
Testmetadata <- metadata[metadata$cellbarcode %in% test_barcodes,]

all_features <- rownames(SeuratCounts)
total_features <- length(all_features)
chunk_size <- 2000
all_results <- data.frame()

for (start_index in seq(22001, total_features, by = chunk_size)) {
  end_index <- min(start_index + chunk_size - 1, total_features)
  feature_subset <- all_features[start_index:end_index]
  
  cat(paste0("Processing features ", start_index, " to ", end_index, "\n"))
  
  TrainData <- prepare_df_classifier(TrainSubset, feature_subset)
  TestData <- prepare_df_classifier(TestSubset, feature_subset)
  
  Analysis <- paste0("Features_", start_index, "_to_", end_index)
  results <- train_classifier(TrainData, TestData, Testmetadata, "all", PerCelltype = TRUE, foldNr = i, Analysis = Analysis)
  
  # Add feature range to results
  results$FeatureRange <- paste0(start_index, "-", end_index)
  
  print(results)
  
  save_objects(results, fold_number = start_index, path_prefix = "/hpc/hers_en/rballieux/AlzheimersProofConcept/Classifiers/Baseline/2000randomfeatures/")
  
  # Append results to all_results
  all_results <- rbind(all_results, results)
}



write.csv(all_results, "/hpc/hers_en/rballieux/AlzheimersProofConcept/Classifiers/Baseline/2000randomfeatures/ClassifierResults2000randomfeats.csv", row.names = FALSE)
library(caret)
library(ROCR)
library(glmnet)
library(randomForest)

source("/hpc/hers_en/rballieux/AlzheimersProofConcept/Classifiers/Scripts/PlotROC.R")

train_classifier <- function(TrainData, TestData, metadata, classifierType = "all", PerCelltype = FALSE, foldNr, Analysis, preloaded_models = NULL, ROCoutput_dir = NULL) {
  # Separate features and labels for the entire dataset
  train_labels <- factor(TrainData$Status, levels = c(0, 1))
  train_features <- TrainData[, -ncol(TrainData)]
  
  test_labels <- factor(TestData$Status, levels = c(0, 1))
  test_features <- TestData[, -ncol(TestData)]
  
  # Define scaling process: Min-Max scaling
  preProcValues <- preProcess(train_features, method = 'range')
  train_features <- predict(preProcValues, train_features)
  test_features <- predict(preProcValues, test_features)
  
  # Initialize a list to store results
  results_list <- list()
  
  classifiers <- if (classifierType == "all") c("glm", "randomForest") else c(classifierType)
  
  models <- list()
  
  # Initialize lists to store ROC data
  roc_data_list <- list()
  
  # Train the models if preloaded_models is not provided
  if (is.null(preloaded_models)) {
    for (cls in classifiers) {
      if (cls == "glm") {
        model <- cv.glmnet(as.matrix(train_features), as.numeric(train_labels) - 1, family = "binomial", alpha = 1, nfolds = 3, nlambda = 50)
        models$glm <- model
      
        coef_matrix <- coef(model, s = "lambda.min")
        top_features_glm <- names(sort(abs(coef_matrix[,1]), decreasing = TRUE)[1:3])
        save_objects(coef_matrix, fold_number = foldNr, path_prefix = paste0("/hpc/hers_en/rballieux/AlzheimersProofConcept/Classifiers/RunFiles/", Analysis))
      
      } else if (cls == "randomForest") {
        model <- randomForest(x = train_features, y = train_labels, ntree = 200, mtry = sqrt(ncol(train_features)), nodesize = 1, maxnodes = NULL, replace = FALSE)
        models$randomForest <- model
      
        importance_vals <- round(importance(model), 2)
        top_features_rf <- names(sort(importance_vals[,"MeanDecreaseGini"], decreasing = TRUE)[1:3])
        save_objects(importance_vals, fold_number = foldNr, path_prefix = paste0("/hpc/hers_en/rballieux/AlzheimersProofConcept/Classifiers/RunFiles/", Analysis))
      }
    }
  
    save_objects(models, fold_number = foldNr, path_prefix = paste0("/hpc/hers_en/rballieux/AlzheimersProofConcept/Classifiers/RunFiles/", Analysis))

  } else {
  models <- preloaded_models
  
  coef_matrix <- coef(models$glm, s = "lambda.min")
  top_features_glm <- names(sort(abs(coef_matrix[,1]), decreasing = TRUE)[1:3])
  
  importance_vals <- round(importance(models$randomForest), 2)
  top_features_rf <- names(sort(importance_vals[,"MeanDecreaseGini"], decreasing = TRUE)[1:3])
}
  
  # Determine cell types to loop through
  if (PerCelltype) {
    celltypes <- unique(metadata$predicted.id)
    celltypes <- c("all celltypes", celltypes)
  } else {
    celltypes <- "all celltypes"
  }
  
  for (celltype in celltypes) {
    if (celltype == "all celltypes") {
      subset_TestData <- TestData
      subset_test_features <- test_features
      subset_test_labels <- test_labels
    } else {
      barcodes <- metadata$cellbarcode[metadata$predicted.id == celltype]
      subset_TestData <- TestData[rownames(TestData) %in% barcodes, ]
      subset_test_features <- subset_TestData[, -ncol(subset_TestData)]
      subset_test_labels <- factor(subset_TestData$Status, levels = c(0, 1))
      
      # Scale the features using the same preProcess object
      subset_test_features <- predict(preProcValues, subset_test_features)
    }
    
    for (cls in classifiers) {
      model <- models[[cls]]
      
      if (cls == "glm") {
        train_pred_probs <- predict(model, newx = as.matrix(train_features), s = "lambda.min", type = "response")
        test_pred_probs <- predict(model, newx = as.matrix(subset_test_features), s = "lambda.min", type = "response")
        
        top_features <- top_features_glm
        
      } else if (cls == "randomForest") {
        train_pred_probs <- predict(model, newdata = train_features, type = "prob")[,2]
        test_pred_probs <- predict(model, newdata = subset_test_features, type = "prob")[,2]
        
        top_features <- top_features_rf
      }
      
      # Store the ROC data
      roc_data_list[[paste(cls, celltype, sep = "_")]] <- list(
        classifier = cls,
        celltype = celltype,
        true_labels = subset_test_labels,
        predicted_probs = test_pred_probs
      )
      
      train_predictions <- factor(ifelse(train_pred_probs >= 0.5, 1, 0), levels = c(0, 1))
      test_predictions <- factor(ifelse(test_pred_probs >= 0.5, 1, 0), levels = c(0, 1))
      
      train_confusionMatrix <- confusionMatrix(train_predictions, train_labels)
      test_confusionMatrix <- confusionMatrix(test_predictions, subset_test_labels)
      
      if (is.null(preloaded_models)) {
        train_auc <- ROCR::performance(ROCR::prediction(train_pred_probs, train_labels), "auc")@y.values[[1]]
      } else {
        train_auc <- NA
      }
      test_auc <- ROCR::performance(ROCR::prediction(test_pred_probs, subset_test_labels), "auc")@y.values[[1]]
      
      train_accuracy <- train_confusionMatrix$overall['Accuracy']
      test_accuracy <- test_confusionMatrix$overall['Accuracy']
      
      train_precision <- posPredValue(train_predictions, train_labels)
      test_precision <- posPredValue(test_predictions, subset_test_labels)
      
      train_recall <- sensitivity(train_predictions, train_labels)
      test_recall <- sensitivity(test_predictions, subset_test_labels)
      
      train_f1_score <- (2 * train_precision * train_recall) / (train_precision + train_recall)
      test_f1_score <- (2 * test_precision * test_recall) / (test_precision + test_recall)
      
      # Compile results
      results_row <- data.frame(
        Classifier = cls,
        CellType = celltype,
        NrOfCells = nrow(subset_TestData),
        TrainAUC = train_auc,
        TestAUC = test_auc,
        TrainAccuracy = train_accuracy,
        TestAccuracy = test_accuracy,
        TrainPrecision = train_precision,
        TestPrecision = test_precision,
        TrainRecall = train_recall,
        TestRecall = test_recall,
        TrainF1Score = train_f1_score,
        TestF1Score = test_f1_score,
        MostImportantFeature = top_features[1],
        SecondMostImportantFeature = top_features[2],
        ThirdMostImportantFeature = top_features[3]
      )
      
      results_list[[paste(cls, celltype, sep = "_")]] <- results_row
    }
  }
  
  results_df <- do.call(rbind, results_list)
  row.names(results_df) <- NULL
  
  save_objects(roc_data_list, fold_number = foldNr, path_prefix = paste0(ROCoutput_dir, Analysis))
  
  #Create ROC-curves
  if (!is.null(ROCoutput_dir)) {
    for (roc_item in roc_data_list) {
      plot_individual_roc(roc_item, ROCoutput_dir, Analysis)
    }
    
    for (celltype in celltypes) {
      plot_combined_roc(roc_data_list, ROCoutput_dir, Analysis, celltype)
    }
  }
  
  return(results_df)
}

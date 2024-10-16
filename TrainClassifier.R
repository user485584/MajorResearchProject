# Load necessary libraries
library(caret)
library(ROCR)
library(glmnet)
library(randomForest)


# ============================
# Configuration Section
# ============================
# Update these paths as needed to configure the script
config <- list(
  PathToScriptsFolder = "/hpc/hers_en/rballieux/AlzheimersProofConcept/Classifiers/Scripts/", # Folder containing helper scripts
  OutputFolder = "/hpc/hers_en/rballieux/AlzheimersProofConcept/Classifiers/RunFiles/"        # Base folder for all outputs
)


# ============================
# Sourcing External Scripts
# ============================
# Source the necessary R scripts using the PathToScriptsFolder from config
source(paste0(config$PathToScriptsFolder, "PlotROC.R"))
source(paste0(config$PathToScriptsFolder, "SaveObjects.R"))


# ----------------------------
# Function: train_classifier
# ----------------------------

#' Train and Evaluate Machine Learning Classifiers on scRNA-seq Data
#'
#' This function trains specified machine learning classifiers on training data,
#' evaluates their performance on testing data, and identifies key transcriptomic
#' features contributing to disease classification. It supports classification
#' across all cell types or specific cell types and generates ROC curves for performance assessment.
#'
#' @param TrainData A data frame containing training data, including features and the 'Status' label.
#' @param TestData A data frame containing testing data, structured similarly to TrainData.
#' @param metadata A data frame containing metadata for cell types, including 'predicted.id' and 'cellbarcode'.
#' @param classifierType A string specifying which classifiers to use ("glm", "randomForest", or "all"). Default is "all".
#' @param PerCelltype A boolean indicating whether to perform classification per cell type. Default is FALSE.
#' @param foldNr An integer indicating the current fold number in cross-validation.
#' @param Analysis A string specifying the analysis name for saving outputs.
#' @param preloaded_models An optional list of pre-trained models. Default is NULL.
#' @param ROCoutput_dir An optional directory path for saving ROC curve outputs. Default is NULL.
#'
#' @return A data frame summarizing classifier performance metrics across classifiers and cell types.
train_classifier <- function(TrainData, TestData, metadata, classifierType = "all", PerCelltype = FALSE, foldNr, Analysis, preloaded_models = NULL, ROCoutput_dir = NULL) {
  
  # ----------------------------
  # Step 1: Prepare Labels and Features
  # ----------------------------
  # Convert 'Status' to a factor with levels 0 (Healthy) and 1 (Diseased)
  train_labels <- factor(TrainData$Status, levels = c(0, 1))
  train_features <- TrainData[, -ncol(TrainData)]
  
  test_labels <- factor(TestData$Status, levels = c(0, 1))
  test_features <- TestData[, -ncol(TestData)]
  
  # ----------------------------
  # Step 2: Feature Scaling
  # ----------------------------
  
  # Apply Min-Max scaling to training features
  preProcValues <- preProcess(train_features, method = 'range')
  train_features <- predict(preProcValues, train_features)
  test_features <- predict(preProcValues, test_features)
  
  
  # ----------------------------
  # Step 3: Initialise Storage Structures
  # ----------------------------
  
  results_list <- list()    # To store results for each classifier and cell type
  roc_data_list <- list()   # To store ROC data for each classifier and cell type
  models <- list()          # To store trained models
  
  
  # Determine which classifiers to use
  if (classifierType == "all") {
    classifiers <- c("glm", "randomForest")
  } else {
    classifiers <- c(classifierType)
  }
  
  
  # ----------------------------
  # Step 4: Train or Load Models
  # ----------------------------
  
  
  if (is.null(preloaded_models)) {
    # Train classifiers if preloaded models are not provided
    for (cls in classifiers) {
      if (cls == "glm") {
        # Train Logistic Regression with Lasso (alpha = 1) using cross-validation
        model <- cv.glmnet(as.matrix(train_features), as.numeric(train_labels) - 1, family = "binomial", alpha = 1, nfolds = 3, nlambda = 50)
        models$glm <- model
        
        # Extract coefficients at the lambda that minimizes the cross-validated error
        coef_matrix <- coef(model, s = "lambda.min")
        
        # Identify top 3 features based on absolute coefficient values
        top_features_glm <- names(sort(abs(coef_matrix[,1]), decreasing = TRUE)[1:3])
        
        # Save the coefficient matrix for later reference
        save_objects(coef_matrix, fold_number = foldNr, path_prefix = paste0(config$OutputFolder, Analysis))
      
      } else if (cls == "randomForest") {
        # Train Random Forest classifier
        model <- randomForest(x = train_features, y = train_labels, ntree = 200, mtry = floor(sqrt(ncol(train_features))), nodesize = 1, maxnodes = NULL, replace = FALSE)
        models$randomForest <- model
      
        # Extract importance values
        importance_vals <- round(importance(model), 2)
        
        # Identify top 3 features based on MeanDecreaseGini
        top_features_rf <- names(sort(importance_vals[,"MeanDecreaseGini"], decreasing = TRUE)[1:3])
        
        # Save the importance values for later reference
        save_objects(importance_vals, fold_number = foldNr, path_prefix = paste0(config$OutputFolder, Analysis))
      }
    }
  
    # Save all trained models
    save_objects(models, fold_number = foldNr, path_prefix = paste0(config$OutputFolder, Analysis))

  } else {
  # Use preloaded models
  models <- preloaded_models
  
  # Extract top features from preloaded Logistic Regression model
  coef_matrix <- coef(models$glm, s = "lambda.min")
  top_features_glm <- names(sort(abs(coef_matrix[,1]), decreasing = TRUE)[1:3])
  
  # Extract top features from preloaded Random Forest model
  importance_vals <- round(importance(models$randomForest), 2)
  top_features_rf <- names(sort(importance_vals[,"MeanDecreaseGini"], decreasing = TRUE)[1:3])
}
  
  
  # ----------------------------
  # Step 5: Determine Cell Types for Analysis
  # ----------------------------
  
  
  if (PerCelltype) {
    # Include all unique cell types plus a combined category
    celltypes <- unique(metadata$predicted.id)
    celltypes <- c("all celltypes", celltypes)
  } else {
    # Analyse all cell types together
    celltypes <- "all celltypes"
  }
  
  # ----------------------------
  # Step 6: Loop Through Cell Types and Classifiers
  # ----------------------------
  
  for (celltype in celltypes) {
    # Subset Test Data based on current cell type
    if (celltype == "all celltypes") {
      subset_TestData <- TestData
      subset_test_features <- test_features
      subset_test_labels <- test_labels
    } else {
      # Identify barcodes corresponding to the current cell type
      barcodes <- metadata$cellbarcode[metadata$predicted.id == celltype]
      
      # Subset TestData for the current cell type
      subset_TestData <- TestData[rownames(TestData) %in% barcodes, ]
      subset_test_features <- subset_TestData[, -ncol(subset_TestData)]
      subset_test_labels <- factor(subset_TestData$Status, levels = c(0, 1))
      
      # Scale the subset test features using the same pre-processing as training data
      subset_test_features <- predict(preProcValues, subset_test_features)
    }
    
    # Loop through each classifier to perform predictions and evaluations
    for (cls in classifiers) {
      model <- models[[cls]]
      
      if (cls == "glm") {
        # Predict probabilities using Logistic Regression
        train_pred_probs <- predict(model, newx = as.matrix(train_features), s = "lambda.min", type = "response")
        test_pred_probs <- predict(model, newx = as.matrix(subset_test_features), s = "lambda.min", type = "response")
        
        # Assign top features for Logistic Regression
        top_features <- top_features_glm
        
      } else if (cls == "randomForest") {
        # Predict probabilities using Random Forest
        train_pred_probs <- predict(model, newdata = train_features, type = "prob")[,2]
        test_pred_probs <- predict(model, newdata = subset_test_features, type = "prob")[,2]
        
        # Assign top features for Random Forest
        top_features <- top_features_rf
      }
      
      # Store ROC data for the current classifier and cell type
      roc_data_list[[paste(cls, celltype, sep = "_")]] <- list(
        classifier = cls,
        celltype = celltype,
        true_labels = subset_test_labels,
        predicted_probs = test_pred_probs
      )
      
      # Convert predicted probabilities to binary class predictions using a threshold of 0.5
      train_predictions <- factor(ifelse(train_pred_probs >= 0.5, 1, 0), levels = c(0, 1))
      test_predictions <- factor(ifelse(test_pred_probs >= 0.5, 1, 0), levels = c(0, 1))
      
      # Generate confusion matrices for training and testing data
      train_confusionMatrix <- confusionMatrix(train_predictions, train_labels)
      test_confusionMatrix <- confusionMatrix(test_predictions, subset_test_labels)
      
      # Calculate Area Under the Curve (AUC) for training and testing data
      if (is.null(preloaded_models)) {
        train_auc <- ROCR::performance(ROCR::prediction(train_pred_probs, train_labels), "auc")@y.values[[1]]
      } else {
        train_auc <- NA # AUC not available for training set if using preloaded models
      }
      test_auc <- ROCR::performance(ROCR::prediction(test_pred_probs, subset_test_labels), "auc")@y.values[[1]]
      
      # Extract performance metrics from confusion matrices
      train_accuracy <- train_confusionMatrix$overall['Accuracy']
      test_accuracy <- test_confusionMatrix$overall['Accuracy']
      
      train_precision <- posPredValue(train_predictions, train_labels)
      test_precision <- posPredValue(test_predictions, subset_test_labels)
      
      train_recall <- sensitivity(train_predictions, train_labels)
      test_recall <- sensitivity(test_predictions, subset_test_labels)
      
      train_f1_score <- (2 * train_precision * train_recall) / (train_precision + train_recall)
      test_f1_score <- (2 * test_precision * test_recall) / (test_precision + test_recall)
      
      # Compile results into a data frame row
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
      
      # Add the results row to the results list
      results_list[[paste(cls, celltype, sep = "_")]] <- results_row
    }
  }
  
  # ----------------------------
  # Step 7: Combine All Results
  # ----------------------------
  
  # Combine all individual results into a single data frame
  results_df <- do.call(rbind, results_list)
  row.names(results_df) <- NULL
  
  
  # ----------------------------
  # Step 8: Save ROC Data
  # ----------------------------
  
  if (!is.null(ROCoutput_dir)) {
    # Save ROC data objects for later use
    save_objects(roc_data_list, fold_number = foldNr, 
                 path_prefix = paste0(ROCoutput_dir, Analysis))
  }
  
  # ----------------------------
  # Step 9: Generate ROC Curves
  # ----------------------------
  if (!is.null(ROCoutput_dir)) {
    # Plot individual ROC curves for each classifier and cell type
    for (roc_item in roc_data_list) {
      plot_individual_roc(roc_item, ROCoutput_dir, Analysis)
    }
    
    # Plot combined ROC curves for each cell type across classifiers
    for (celltype in celltypes) {
      plot_combined_roc(roc_data_list, ROCoutput_dir, Analysis, celltype)
    }
  }
  
  # ----------------------------
  # Step 10: Return Compiled Results
  # ----------------------------
  
  return(results_df)
}

# Load Necessary Libraries
library(pROC)

# ----------------------------
# Configuration
# ----------------------------
config <- list(
  OutputFolder = "/hpc/hers_en/rballieux/AlzheimersProofConcept/Classifiers/RunFiles/"        # Base folder for all outputs
)

# ----------------------------
# Function: plot_individual_roc
# ----------------------------

#' Plot and Save Individual ROC Curve
#'
#' This function generates and saves an individual ROC curve plot for a given classifier and cell type.
#' The plot includes the ROC curve with AUC annotation.
#'
#' @param roc_item A list containing the following elements:
#'                 - `classifier`: A string indicating the classifier type (e.g., "glm", "randomForest").
#'                 - `celltype`: A string indicating the cell type.
#'                 - `true_labels`: A numeric or binary vector of true class labels.
#'                 - `predicted_probs`: A numeric vector of predicted probabilities for the positive class.
#' @param output_dir A string specifying the directory path where the ROC plot should be saved.
#'                   Defaults to the `OutputFolder` specified in the `config` list.
#' @param analysis A string specifying the analysis name to include in the plot filename.
#'
#' @return NULL. The function saves the ROC plot as a PNG file in the specified directory.
plot_individual_roc <- function(roc_item, output_dir = config$OutputFolder, analysis) {
  
  # ----------------------------
  # Step 1: Generate ROC Curve
  # ----------------------------
  # Create ROC curve object using pROC
  roc_curve <- pROC::roc(roc_item$true_labels, roc_item$predicted_probs)
  
  # ----------------------------
  # Step 2: Define Plot Filename
  # ----------------------------
  plot_filename <- paste0(output_dir, analysis, "_", roc_item$classifier, "_", gsub(" ", "_", roc_item$celltype), "_ROC.png")
  
  # ----------------------------
  # Step 3: Plot and Save ROC Curve
  # ----------------------------
  # Open a PNG device to save the plot
  png(plot_filename, width = 800, height = 600)
  
  # Set plot parameters for a square plot
  par(pty = "s")
  
  # Plot the ROC curve with AUC annotation
  pROC::plot.roc(roc_curve, 
                 main = paste("ROC Curve -", roc_item$classifier, "-", roc_item$celltype), 
                 col = ifelse(roc_item$classifier == "glm", "blue", "red"), 
                 lwd = 4, 
                 legacy.axes = TRUE, 
                 print.auc=TRUE,
                 xlab = "False Positive Rate (1 - Specificity)", 
                 ylab = "True Positive Rate (Sensitivity)")
  
  # Close the PNG device to save the plot
  dev.off()
  
  # Inform the user that the plot has been saved
  message("Saved individual ROC plot to ", plot_filename)
}


# ----------------------------
# Function: plot_combined_roc
# ----------------------------

#' Plot and Save Combined ROC Curves for a Cell Type
#'
#' This function generates and saves a combined ROC curve plot for multiple classifiers within a specific cell type.
#' The plot includes ROC curves from different classifiers with distinct colors and a legend.
#'
#' @param roc_data_list A named list where each element is a list containing the following:
#'                      - `classifier`: A string indicating the classifier type (e.g., "glm", "randomForest").
#'                      - `celltype`: A string indicating the cell type.
#'                      - `true_labels`: A numeric or binary vector of true class labels.
#'                      - `predicted_probs`: A numeric vector of predicted probabilities for the positive class.
#' @param output_dir A string specifying the directory path where the ROC plots should be saved.
#'                   Defaults to the `OutputFolder` specified in the `config` list.
#' @param analysis A string specifying the analysis name to include in the plot filename.
#' @param celltype A string specifying the cell type for which combined ROC curves are plotted.
#'
#' @return NULL. The function saves the combined ROC plot as a PNG file in the specified directory.
plot_combined_roc <- function(roc_data_list, output_dir = config$OutputFolder, analysis, celltype) {
  
  # ----------------------------
  # Step 1: Define Plot Filename
  # ----------------------------
  plot_filename_combined <- paste0(output_dir, analysis, "_Combined_", gsub(" ", "_", celltype), "_ROC.png")
  
  # ----------------------------
  # Step 2: Plot and Save ROC Curve
  # ----------------------------
  # Open a PNG device to save the plot
  png(plot_filename_combined, width = 800, height = 600)
  
  # Set plot parameters for a square plot
  par(pty = "s")
  
  # Create ROC curve object using pROC for Logistic Regression
  roc_curve_glm <- pROC::roc(roc_data_list[[paste0("glm_", celltype)]]$true_labels, roc_data_list[[paste0("glm_", celltype)]]$predicted_probs)
  
  # Plot the ROC curve with AUC annotation for Logistic Regression
  pROC::plot.roc(roc_curve_glm, 
                 main = paste("Combined ROC Curve -", celltype), 
                 col = "blue", 
                 lwd = 4, 
                 legacy.axes = TRUE, 
                 print.auc = TRUE,
                 xlab = "False Positive Rate (1 - Specificity)", 
                 ylab = "True Positive Rate (Sensitivity)")
  
  # Create ROC curve object using pROC for Random Forest
  roc_curve_rf <- pROC::roc(roc_data_list[[paste0("randomForest_", celltype)]]$true_labels, roc_data_list[[paste0("randomForest_", celltype)]]$predicted_probs)
  
  # Add Random Forest ROC curve to plot above
  pROC::plot.roc(roc_curve_rf, add = TRUE, col = "red", lwd = 4, print.auc = TRUE, print.auc.y = .4)
  
  #Add Legend
  legend("bottomright", legend = c("Glm", "RandomForest"), col = c("blue", "red"), lwd = 4)
  
  # Close the PNG device to save the plot
  dev.off()
  
  # Inform the user that the combined ROC plot has been saved
  message("Saved combined ROC plot to ", plot_filename_combined)
}
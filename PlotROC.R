library(pROC)

# Function to plot individual ROC curves
plot_individual_roc <- function(roc_item, output_dir, analysis) {
  roc_curve <- pROC::roc(roc_item$true_labels, roc_item$predicted_probs)
  plot_filename <- paste0(output_dir, analysis, "_", roc_item$classifier, "_", gsub(" ", "_", roc_item$celltype), "_ROC.png")
  
  png(plot_filename, width = 800, height = 600)
  par(pty = "s")
  pROC::plot.roc(roc_curve, 
                 main = paste("ROC Curve -", roc_item$classifier, "-", roc_item$celltype), 
                 col = ifelse(roc_item$classifier == "glm", "blue", "red"), 
                 lwd = 4, 
                 legacy.axes = TRUE, 
                 print.auc=TRUE,
                 xlab = "False Positive Rate (1 - Specificity)", 
                 ylab = "True Positive Rate (Sensitivity)")
  dev.off()
}

# Function to plot combined ROC curves for each cell type
plot_combined_roc <- function(roc_data_list, output_dir, analysis, celltype) {
  plot_filename_combined <- paste0(output_dir, analysis, "_Combined_", gsub(" ", "_", celltype), "_ROC.png")
  
  png(plot_filename_combined, width = 800, height = 600)
  par(pty = "s")
  
  # Plot GLM ROC curve
  roc_curve_glm <- pROC::roc(roc_data_list[[paste0("glm_", celltype)]]$true_labels, roc_data_list[[paste0("glm_", celltype)]]$predicted_probs)
  pROC::plot.roc(roc_curve_glm, 
                 main = paste("Combined ROC Curve -", celltype), 
                 col = "blue", 
                 lwd = 4, 
                 legacy.axes = TRUE, 
                 print.auc = TRUE,
                 xlab = "False Positive Rate (1 - Specificity)", 
                 ylab = "True Positive Rate (Sensitivity)")
  
  # Add RandomForest ROC curve
  roc_curve_rf <- pROC::roc(roc_data_list[[paste0("randomForest_", celltype)]]$true_labels, roc_data_list[[paste0("randomForest_", celltype)]]$predicted_probs)
  pROC::plot.roc(roc_curve_rf, add = TRUE, col = "red", lwd = 4, print.auc = TRUE, print.auc.y = .4)
  
  legend("bottomright", legend = c("Glm", "RandomForest"), col = c("blue", "red"), lwd = 4)
  dev.off()
}
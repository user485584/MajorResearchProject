library(magrittr)

combine_folds <- function(FilePath, from, to, by, OutputPath) {
  # Initialize an empty list to store the folds
  folds_list <- list()
  
  # Loop over the specified range
  for (i in seq(from, to, by = by)) {
    # Construct the file path
    file_path <- paste0(FilePath, "results_Fold", i, ".rds")
    
    # Read the RDS file and store it in the list
    folds_list[[paste0("fold", i)]] <- readRDS(file_path)
  }
  
  # Combine all folds into a single data frame
  combined_folds <- do.call(rbind, folds_list)
  
  rownames(combined_folds) <- NULL
  
  # Save the combined data frame as an RDS file
  saveRDS(combined_folds, paste0(OutputPath, "resultsFeats", from, "-", to, ".rds"))
  
  # Save the combined data frame as a CSV file
  write.csv(combined_folds, paste0(OutputPath, "resultsFeats", from, "-", to, ".csv"), row.names = FALSE)
  
  return(combined_folds)  # Return the combined data frame if needed
}


combine_folds_coef_importances <- function(Disease, Analysis, FilePath, from, to, by, OutputPath) {
  # Initialize empty lists to store the coefficients and importances for each range of features
  coef_list <- list()
  import_list <- list()
  
  # Loop over the specified range of feature sets (e.g., 1-500, 501-1000, etc.)
  for (j in seq(from, to, by = by)) {
    # Construct the file paths for coefficients and importances
    
    end <- j + by - 1
    if (end > to){
      end <- to
    }
    
    Fold <- paste0(j,"-",end)
    
    coef_path <- paste0(FilePath, "Features_", j, "_to_", end, "coef_matrix_Fold1.rds")
    import_path <- paste0(FilePath, "Features_", j, "_to_", end, "importance_vals_Fold1.rds")
    
    
    # Import the RDS files into separate variables
    coef_df <- readRDS(coef_path) %>% as.data.frame() %>% rownames_to_column()
    
    #modify dfs
    coef_df <- coef_df %>% 
      rename(gene = rowname, importance = s1) %>% 
      mutate(absimportance = abs(importance)) %>%
      mutate(Disease = Disease) %>%
      mutate(Analysis = Analysis) %>%
      mutate(Classifier = "glm") %>%
      mutate(Fold = Fold) %>%
      select(Disease, Analysis, Fold, Classifier, gene, absimportance, importance) %>%
      arrange(-absimportance)
    
    import_df <- readRDS(import_path) %>% as.data.frame() %>% rownames_to_column()
    
    import_df <- import_df %>% 
      rename(gene = rowname, importance = MeanDecreaseGini) %>%
      mutate(absimportance = abs(importance)) %>%
      mutate(Disease = Disease) %>%
      mutate(Analysis = Analysis) %>%
      mutate(Classifier = "RandomForest") %>%
      mutate(Fold = Fold) %>%
      select(Disease, Analysis, Fold, Classifier, gene, absimportance, importance) %>%
      arrange(-absimportance)
    
    
    # Add the manipulated data frames to the lists
    coef_list[[length(coef_list) + 1]] <- coef_df
    import_list[[length(import_list) + 1]] <- import_df
  }
  
  # Combine all coefficients and importances across feature ranges
  combined_coef <- do.call(rbind, coef_list)
  combined_import <- do.call(rbind, import_list)
  
  # Save the combined results as RDS files
  saveRDS(combined_coef, file = paste0(OutputPath, "combined_coef_", from, "_to_", to, ".rds"))
  saveRDS(combined_import, file = paste0(OutputPath, "combined_import_", from, "_to_", to, ".rds"))
  
  write_csv(combined_coef, file = paste0(OutputPath, "combined_coef_", from, "_to_", to, ".csv"))
  write_csv(combined_import, file = paste0(OutputPath, "combined_import_", from, "_to_", to, ".csv"))
  
  # Return the combined coefficients and importances
  return(list(combined_coef = combined_coef, combined_import = combined_import))
}


# Example usage:
# combined_folds <- combine_folds(
#   FilePath = "/Users/rutger/Desktop/Bureaublad/School/Master/Stage/Classifiers/Baseline/500randomfeatures/",
#   from = 1,
#   to = 24501,
#   by = 500,
#   OutputPath = "/Users/rutger/Desktop/Bureaublad/School/Master/Stage/Classifiers/Baseline/500randomfeatures/"
# )


combine_folds_coef_importances_nonrandomfeatures <- function(Disease, Analysis, NrFolds, FilePath, OutputPath){
  
  coef_list <- list()
  import_list <- list()
  
  for (Fold in 1:NrFolds){
    
    if (Analysis == "TopPCs"){
      FeatureSequence <- c("10PCs", "50PCs", "100PCs", "250PCs", "500PCs", "1000PCs")
      
    } else if (Analysis == "VariableFeatures"){
      FeatureSequence <- c(100, 500, 1000, 1500, 2000)
      
    } else if (Analysis == "DEWGCNA"){
      FeatureSequence <- ""
      
    } else {
      stop("Error: Invalid value for 'Analysis'. Please provide one of 'TopPCs', 'VariableFeatures', or 'DEWGCNA'.")
    }
    
    
    for (NrFeatures in FeatureSequence) {
      
      if (Analysis == "TopPCs"){
        prefix <- paste0("Top", NrFeatures)
        
      } else if (Analysis == "VariableFeatures"){
        prefix <- paste0("Top", NrFeatures, Analysis)
        
      } else if (Analysis == "DEWGCNA"){
        prefix <- Analysis
        
      } else {
        stop("Error: Invalid value for 'Analysis'. Please provide one of 'TopPCs', 'VariableFeatures', or 'DEWGCNA'.")
      }
      
      coef_path <- paste0(FilePath, prefix, "coef_matrix_Fold", Fold, ".rds")
      import_path <- paste0(FilePath, prefix, "importance_vals_Fold", Fold, ".rds")
      
      coef_df <- readRDS(coef_path) %>% as.data.frame() %>% rownames_to_column()
      
      #modify dfs
      coef_df <- coef_df %>% 
        rename(gene = rowname, importance = s1) %>% 
        mutate(absimportance = abs(importance)) %>%
        mutate(Disease = Disease) %>%
        mutate(Analysis = Analysis) %>%
        mutate(Classifier = "glm") %>%
        mutate(Fold = paste0(NrFeatures, "Fold", Fold)) %>%
        select(Disease, Analysis, Fold, Classifier, gene, absimportance, importance) %>%
        arrange(-absimportance)
      
      import_df <- readRDS(import_path) %>% as.data.frame() %>% rownames_to_column()
      
      import_df <- import_df %>% 
        rename(gene = rowname, importance = MeanDecreaseGini) %>%
        mutate(absimportance = abs(importance)) %>%
        mutate(Disease = Disease) %>%
        mutate(Analysis = Analysis) %>%
        mutate(Classifier = "RandomForest") %>%
        mutate(Fold = paste0(NrFeatures, "Fold", Fold)) %>%
        select(Disease, Analysis, Fold, Classifier, gene, absimportance, importance) %>%
        arrange(-absimportance)
      
      
      # Add the manipulated data frames to the lists
      coef_list[[length(coef_list) + 1]] <- coef_df
      import_list[[length(import_list) + 1]] <- import_df
    }
  }
  
  # Combine all coefficients and importances across feature ranges
  combined_coef <- do.call(rbind, coef_list)
  combined_import <- do.call(rbind, import_list)
  
  # Save the combined results as RDS files
  saveRDS(combined_coef, file = paste0(OutputPath, "combined_coef_", Analysis, ".rds"))
  saveRDS(combined_import, file = paste0(OutputPath, "combined_import_", Analysis, ".rds"))
  
  write_csv(combined_coef, file = paste0(OutputPath, "combined_coef_", Analysis, ".csv"))
  write_csv(combined_import, file = paste0(OutputPath, "combined_import_", Analysis, ".csv"))
  
  # Return the combined coefficients and importances
  return(list(combined_coef = combined_coef, combined_import = combined_import))
}

# Load Necessary Libraries
library(magrittr)
library(dplyr)


# ----------------------------
# Function: combine_folds
# ----------------------------

#' Combine Multiple Fold Results into Single Files
#'
#' This function reads multiple RDS files corresponding to different folds, combines them into a single data frame,
#' and saves the combined results as both RDS and CSV files.
#'
#' @param FilePath String. Directory path where fold RDS files are located.
#' @param from Integer. Starting fold number.
#' @param to Integer. Ending fold number.
#' @param by Integer. Step size for fold numbers.
#' @param OutputPath String. Directory path to save the combined results.
#'
#' @return Combined data frame of all folds.
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



# ----------------------------
# Function: combine_folds_coef_importances
# ----------------------------

#' Combine Coefficients and Importances Across Feature Ranges
#'
#' This function reads coefficient and importance RDS files for specified feature ranges and combines them into single data frames.
#' The combined results are saved as both RDS and CSV files.
#'
#' @param Disease String. Disease label.
#' @param Analysis String. Type of analysis (in this case a Number of RandomFeatures).
#' @param FilePath String. Directory path where coefficient and importance RDS files are located.
#' @param from Integer. Starting feature range.
#' @param to Integer. Ending feature range.
#' @param by Integer. Step size for feature ranges.
#' @param OutputPath String. Directory path to save the combined results.
#'
#' @return List containing combined coefficients and importances.
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
    
    
    # Read and process coefficients
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
    
    
    # Read and process importances
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



# ----------------------------
# Function: combine_folds_coef_importances_nonrandomfeatures
# ----------------------------

#' Combine Coefficients and Importances for Non-Random Features Across Folds
#'
#' This function processes coefficient and importance RDS files for different analyses and combines them into single data frames.
#' The combined results are saved as both RDS and CSV files.
#'
#' @param Disease String. Disease label.
#' @param Analysis String. Type of analysis (e.g., "TopPCs", "VariableFeatures", "DEWGCNA").
#' @param NrFolds Integer. Number of folds to process.
#' @param FilePath String. Directory path where coefficient and importance RDS files are located.
#' @param OutputPath String. Directory path to save the combined results.
#'
#' @return List containing combined coefficients and importances.
combine_folds_coef_importances_nonrandomfeatures <- function(Disease, Analysis, NrFolds, FilePath, OutputPath){
  
  coef_list <- list()
  import_list <- list()
  
  # Define feature sequences based on analysis type
  FeatureSequence <- switch(
    Analysis,
    "TopPCs" = c("10PCs", "50PCs", "100PCs", "250PCs", "500PCs", "1000PCs"),
    "VariableFeatures" = c(100, 500, 1000, 1500, 2000),
    "DEWGCNA" = "",
    stop("Error: Invalid value for 'Analysis'. Please provide one of 'TopPCs', 'VariableFeatures', or 'DEWGCNA'.")
  )
  
  for (Fold in 1:NrFolds){
    for (NrFeatures in FeatureSequence) {
      prefix <- switch(
        Analysis,
        "TopPCs" = paste0("Top", NrFeatures),
        "VariableFeatures" = paste0("Top", NrFeatures, Analysis),
        "DEWGCNA" = Analysis
      )
      
      coef_path <- paste0(FilePath, prefix, "coef_matrix_Fold", Fold, ".rds")
      import_path <- paste0(FilePath, prefix, "importance_vals_Fold", Fold, ".rds")
      
      
      # Read and process coefficients
      coef_df <- readRDS(coef_path) %>% as.data.frame() %>% rownames_to_column()
      
      coef_df <- coef_df %>% 
        rename(gene = rowname, importance = s1) %>% 
        mutate(absimportance = abs(importance)) %>%
        mutate(Disease = Disease) %>%
        mutate(Analysis = Analysis) %>%
        mutate(Classifier = "glm") %>%
        mutate(Fold = paste0(NrFeatures, "Fold", Fold)) %>%
        select(Disease, Analysis, Fold, Classifier, gene, absimportance, importance) %>%
        arrange(-absimportance)
      
      
      # Read and process importances
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

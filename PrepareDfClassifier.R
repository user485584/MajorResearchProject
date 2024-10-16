# Load Necessary Libraries
library(Seurat)

# ----------------------------
# Function: prepare_df_classifier
# ----------------------------

#' Prepare Data Frame for Classification from Seurat Object
#'
#' This function processes a Seurat object by subsetting it based on specified features,
#' extracting the count data, and preparing a data frame suitable for classification tasks.
#' It also adds a binary 'Status' column indicating disease presence.
#'
#' @param SeuratObject A Seurat object containing single-cell RNA sequencing data.
#' @param features A character vector of gene names or identifiers to subset the Seurat object.
#'
#' @return A data frame where rows represent cells, columns represent gene expression levels,
#'         and an additional 'Status' column indicating disease status (1 for disease, 0 for healthy).
prepare_df_classifier <- function(SeuratObject, features) {
  
  # ----------------------------
  # Step 1: Input Validation
  # ----------------------------
  
  # Check if SeuratObject is a valid Seurat object
  if (!inherits(SeuratObject, "Seurat")) {
    stop("Error: SeuratObject must be a valid Seurat object.")
  }
  
  # Check if 'Status' column exists in metadata
  if (!"Status" %in% colnames(SeuratObject@meta.data)) {
    stop("Error: SeuratObject metadata must contain a 'Status' column.")
  }
  
  # Check if specified features exist in the Seurat object
  missing_features <- setdiff(features, rownames(SeuratObject))
  if (length(missing_features) > 0) {
    warning("Warning: The following features are not present in the Seurat object and will be ignored: ",
            paste(missing_features, collapse = ", "))
    # Remove missing features from the list
    features <- setdiff(features, missing_features)
  }
  
  # Check if there are any features left after removing missing ones
  if (length(features) == 0) {
    stop("Error: No valid features provided for subsetting.")
  }
  
  
  # ----------------------------
  # Step 2: Subset Seurat Object
  # ----------------------------
  # Subset the Seurat object to include only the specified features (genes)
  SubsetData <- SeuratObject[features, ]
  
  # ----------------------------
  # Step 3: Extract Count Data
  # ----------------------------
  # Retrieve the count data from the subsetted Seurat object
  data_matrix <- GetAssayData(object = SubsetData, layer = "counts")
  
  # Convert sparse matrix to a dense matrix and transpose it (cells as rows, genes as columns)
  dense_matrix <- as.matrix(data_matrix)
  transposed_dense_matrix <- t(dense_matrix)
  
  # Convert the transposed matrix to a data frame
  df <- as.data.frame(transposed_dense_matrix)
  
  # ----------------------------
  # Step 4: Add 'Status' Column
  # ----------------------------
  # Add the 'Status' column from metadata
  df$Status <- SubsetData@meta.data$Status
  
  # Convert 'Status' to a binary indicator: 1 for disease (AD, ALS, FTLD), 0 for healthy
  # Adjust the disease labels as necessary based on your specific dataset
  df$Status <- ifelse(df$Status %in% c("AD", "ALS", "FTLD"), 1, 0)
  
  # ----------------------------
  # Step 5: Return Prepared Data Frame
  # ----------------------------
  return(df)
}





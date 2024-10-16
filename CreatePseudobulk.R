# Load Necessary Libraries
library(Seurat)
library(dplyr)

# ----------------------------
# Function: CreatePseudobulk
# ----------------------------

#' Create Pseudobulk Expression Data from Seurat Subset
#'
#' This function aggregates single-cell RNA-seq data into pseudobulk expression profiles by summing
#' expression values for each specified group (e.g., sample or condition). It also prepares the
#' corresponding pseudobulk metadata.
#'
#' @param SeuratSubset A Seurat object containing the subset of single-cell data to be aggregated.
#' @param metadata A data frame containing metadata associated with the cells in SeuratSubset.
#'                 Must include the grouping variable 'id' and the 'Status' column.
#'
#' @return A list containing:
#'         - `pseudobulk_df`: A data frame of aggregated expression values.
#'         - `pseudobulk_metadata`: A data frame of corresponding pseudobulk metadata.

CreatePseudobulk <- function(SeuratSubset, metadata) {
  
  # ----------------------------
  # Step 1: Input Validation
  # ----------------------------
  
  # Check if SeuratSubset is a valid Seurat object
  if (!inherits(SeuratSubset, "Seurat")) {
    stop("Error: SeuratSubset must be a valid Seurat object.")
  }
  
  # Check if metadata is a data frame
  if (!is.data.frame(metadata)) {
    stop("Error: metadata must be a data frame.")
  }
  
  # Check for required columns in metadata
  required_columns <- c("id", "Status")
  missing_columns <- setdiff(required_columns, colnames(metadata))
  if (length(missing_columns) > 0) {
    stop("Error: metadata is missing the following required columns: ", paste(missing_columns, collapse = ", "))
  }
  
  # ----------------------------
  # Step 2: Aggregate Expression Data
  # ----------------------------
  
  # Aggregate expression data by 'id', summing expression values
  pseudobulk_df <- AggregateExpression(object = SeuratSubset, group.by = "id", fun = "sum")
  
  # Convert aggregated RNA assay data to a data frame
  pseudobulk_df <- as.data.frame(as.matrix(pseudobulk_df$RNA))
  
  # ----------------------------
  # Step 3: Clean Column Names
  # ----------------------------
  # Remove leading 'g' from gene names if present
  colnames(pseudobulk_df) <- sub("^g", "", colnames(pseudobulk_df))
  
  # ----------------------------
  # Step 4: Prepare Pseudobulk Metadata
  # ----------------------------
  pseudobulk_metadata <- metadata %>%
    group_by(id) %>%
    summarise(Status = dplyr::first(Status), .groups = 'drop') 
  
  # ----------------------------
  # Step 5: Return Results
  # ----------------------------
  
  # Return a list containing both the pseudobulk data frame and metadata
  return(list(pseudobulk_df = pseudobulk_df, pseudobulk_metadata = pseudobulk_metadata))

}
  
  
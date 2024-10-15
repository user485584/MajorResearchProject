library(Seurat)
library(dplyr)


CreatePseudobulk <- function(SeuratSubset, metadata, type = "standard") {
  if (type == "standard") {
    # Aggregate expression data by 'id', summing expression values
    pseudobulk_df <- AggregateExpression(object = SeuratSubset, group.by = "id", fun = "sum")
    pseudobulk_df <- as.data.frame(as.matrix(pseudobulk_df$RNA))
    
    # Remove 'g' prefix from column names if present
    colnames(pseudobulk_df) <- sub("^g", "", colnames(pseudobulk_df))
    
    # Create pseudobulk metadata
    pseudobulk_metadata <- metadata %>%
      group_by(id) %>%
      summarise(Status = dplyr::first(Status), .groups = 'drop') # Ensuring groups are dropped after summarisation
    
    # Return a list containing both the pseudobulk data frame and metadata
    return(list(pseudobulk_df = pseudobulk_df, pseudobulk_metadata = pseudobulk_metadata))
    
    
  } else if (type == "GRN") {
    SeuratSubset <- subset(SeuratSubset, Status == "AD")
    
    pseudobulk <- AggregateExpression(object = SeuratSubset, group.by = "subs", fun = "sum")
    
    countsPeaks <- as.data.frame(as.matrix(pseudobulk$ATAC)) %>% tibble::rownames_to_column(var = "peakID")
    countsRNA <- as.data.frame(as.matrix(pseudobulk$RNA)) %>% tibble::rownames_to_column(var = "GeneID")
    
    # Remove 'g' prefix from column names if present
    colnames(countsPeaks) <- sub("^g", "", colnames(countsPeaks))
    colnames(countsRNA) <- sub("^g", "", colnames(countsRNA))
    
    pseudobulk_metadata <- metadata %>%
      group_by(subs) %>%
      summarise(Status = dplyr::first(Status), .groups = 'drop')
    
    # Return a placeholder list or some relevant output for now
    return(list(countsPeaks = countsPeaks, countsRNA = countsRNA, metadata = pseudobulk_metadata))
  }
}
  
  
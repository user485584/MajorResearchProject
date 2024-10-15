library(Seurat)

# Define a function to process a Seurat object for specified features and return a data frame
prepare_df_classifier <- function(SeuratObject, features) {
  # Subset Seurat object based on combined features
  SubsetData <- SeuratObject[features, ]
  
  # Get the count data from the subset
  data_matrix <- GetAssayData(object = SubsetData, layer = "counts")
  
  dense_matrix <- as.matrix(data_matrix)
  transposed_dense_matrix <- t(dense_matrix)
  df <- as.data.frame(transposed_dense_matrix)
  
  # Add the 'Status' column from metadata
  df$Status <- SubsetData@meta.data$Status
  
  df$Status <- ifelse(df$Status %in% c("AD", "ALS", "FTLD"), 1, 0)
  
  
  return(df)
}



prepareGRN_df_classifier <- function(SeuratObject, features, target = "separate") {
  
  
  SubsetData <- SeuratObject[features, ]
  
  data_matrix <- GetAssayData(object = SubsetData, layer = "counts")
  
  dense_matrix <- as.matrix(data_matrix)
  transposed_dense_matrix <- t(dense_matrix)
  df <- as.data.frame(transposed_dense_matrix)
  
  if (target == "interaction") {
    df <- df + 1
    
    for (i in 1:nrow(featuresdf)) {
      tf_col <- featuresdf$TF.ENSEMBL[i]
      gene_col <- featuresdf$gene.ENSEMBL[i]
      
      # Check if both columns exist in df
      if (gene_col %in% colnames(df) && tf_col %in% colnames(df)) {
        # Create a new column name based on TF and GeneID
        new_col_name <- paste(featuresdf$TF.name[i], gene_col, tf_col, sep="_")
        
        # Multiply the corresponding columns in df
        df[[new_col_name]] <- df[[gene_col]] * df[[tf_col]]
      }
    }
    
    df <- df[, (length(features) + 1) : ncol(df)]
  }
  
  # Add the Status column
  df$Status <- SubsetData@meta.data$Status
  df$Status <- ifelse(df$Status == "AD", 1, 0)
  
  return(df)

}

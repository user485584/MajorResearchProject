# Load Necessary Libraries
library(magrittr)
library(dplyr)
library(Seurat)


# ----------------------------
# Function: process_metadata
# ----------------------------

#' Process Metadata for a Seurat Object
#'
#' This function renames the 'predicted.id' column to 'orig.substructure' and creates a new
#' 'predicted.id' column based on specified logic.
#'
#' @param seurat_obj A Seurat object whose metadata will be processed.
#'
#' @return The Seurat object with updated metadata.
process_metadata <- function(seurat_obj) {
  # Extract metadata
  metadata <- seurat_obj@meta.data
  
  # Step 1: Rename 'predicted.id' to 'orig.substructure'
  metadata <- metadata %>%
    rename(orig.substructure = predicted.id)
  
  # Step 2: Create a new 'predicted.id' column based on logic
  metadata <- metadata %>%
    mutate(predicted.id = case_when(
      CellClass == "Ex" ~ "Excitatory",
      CellClass == "In" ~ "Inhibitory",
      CellClass == "Vasc" ~ "Vasc",
      CellClass == "Glia" & orig.substructure == "Oligo" ~ "Oligodendrocytes",
      CellClass == "Glia" & orig.substructure == "Astro" ~ "Astrocytes",
      CellClass == "Glia" & orig.substructure == "OPC" ~ "OPCs",
      CellClass == "Glia" & orig.substructure == "Micro" ~ "Microglia",
      TRUE ~ orig.substructure  # Default case to handle any other values
    ))
  
  # Assign modified metadata back to Seurat object
  seurat_obj@meta.data <- metadata
  
  return(seurat_obj)
}


# ----------------------------
# Load Seurat Objects
# ----------------------------
# Load the Seurat objects from the specified dataset paths
C9ALS <- readRDS("/Users/rutger/Desktop/Bureaublad/School/Master/Stage/ALS_FTLD/Data/C9ALS_Seurat_Downsampled.rds")
SALS <- readRDS("/Users/rutger/Desktop/Bureaublad/School/Master/Stage/ALS_FTLD/Data/SALS_Seurat_Downsampled.rds")


# ----------------------------
# Process Metadata for Both Objects
# ----------------------------
# Apply the metadata processing function to both Seurat objects
C9ALS <- process_metadata(C9ALS)
SALS <- process_metadata(SALS)


# ----------------------------
# Verify the Changes
# ----------------------------
# View the first few rows of the updated metadata for verification
head(SALS@meta.data)
head(C9ALS@meta.data)


# ----------------------------
# Save the Modified Seurat Objects
# ----------------------------
# Save the updated Seurat objects back to their respective RDS files
saveRDS(SALS, "/Users/rutger/Desktop/Bureaublad/School/Master/Stage/ALS_FTLD/Data/SALS_Seurat_Downsampled.rds")
saveRDS(C9ALS, "/Users/rutger/Desktop/Bureaublad/School/Master/Stage/ALS_FTLD/Data/C9ALS_Seurat_Downsampled.rds")
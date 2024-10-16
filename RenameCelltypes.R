library(magrittr)
library(dplyr)
library(Seurat)

C9ALS <- readRDS("/Users/rutger/Desktop/Bureaublad/School/Master/Stage/ALS_FTLD/Data/C9ALS_Seurat_Downsampled.rds")


# Step 1: Rename predicted.id to orig.substructure
Sdownmeta <- Sdownmeta %>%
  rename(orig.substructure = predicted.id)

# Step 2: Create a new predicted.id column with the specified logic
Sdownmeta <- Sdownmeta %>%
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

SALS@meta.data <- Sdownmeta

# Verify the changes
head(SALS)
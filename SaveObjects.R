
# ----------------------------
# Configuration
# ----------------------------

# Define configuration settings for script paths and output directories
config <- list(
  OutputFolder = "/hpc/hers_en/rballieux/AlzheimersProofConcept/Classifiers/RunFiles/"        # Base folder for all outputs
)

# ----------------------------
# Function: save_objects
# ----------------------------

#' Save Multiple R Objects to Disk with Specified Fold Number and Path Prefix
#'
#' This function saves multiple R objects to disk as RDS files. Each object is saved with a filename
#' that includes its name and the specified fold number. The function ensures that the output directory
#' exists and provides informative messages about the saving process.
#'
#' @param ... R objects to be saved.
#' @param fold_number Integer indicating the current fold number in cross-validation.
#' @param path_prefix String specifying the directory path where the objects should be saved.
#'                    Defaults to the `OutputFolder` specified in the `config` list.
#'
#' @return NULL. The function saves the objects to disk and prints confirmation messages.

save_objects <- function(..., fold_number, path_prefix = config$OutputFolder) {
  
  # ----------------------------
  # Step 1: Capture Objects and Their Names
  # ----------------------------
  # Capture the list of objects passed to the function
  objects <- list(...)
  
  # Capture the names of the objects using substitute and deparse
  # 'substitute(list(...))[-1]' gets the list of expressions passed to '...'
  # 'deparse' converts them into character strings
  names_list <- sapply(substitute(list(...))[-1], deparse)
  
  # Identify any objects without names
  if (any(is.na(names_list) | names_list == "")) {
    # Find the indices of unnamed or empty-named objects
    empty_names <- which(is.na(names_list) | names_list == "")
    
    # Assign default names to unnamed objects (e.g., object1, object2, ...)
    names_list[empty_names] <- paste0("object", empty_names)
  }
  
  # ----------------------------
  # Step 3: Ensure Output Directory Exists
  # ----------------------------
  
  # Check if the specified output directory exists
  if (!dir.exists(path_prefix)) {
    # Create the directory, including any necessary parent directories
    dir.create(path_prefix, recursive = TRUE)
    
    # Inform the user that the directory was created
    message("Created output directory: ", path_prefix)
  }
  
  # ----------------------------
  # Step 4: Save Each Object to Disk
  # ----------------------------
  
  # Iterate over each object and save it as an RDS file
  for (i in seq_along(objects)) {
    # Construct the file name using the object name and fold number
    file_path <- paste0(path_prefix, names_list[i], "_Fold", fold_number, ".rds")
    
    # Save the object to the specified file path using saveRDS
    saveRDS(objects[[i]], file_path)
    
    # Inform the user that the object was saved successfully
    message("Saved '", names_list[i], "' to ", file_path)
  }
}

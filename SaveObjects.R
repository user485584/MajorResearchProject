


save_objects <- function(..., fold_number, path_prefix = "/hpc/hers_en/rballieux/AlzheimersProofConcept/Classifiers/RunFiles/") {
  # Capture the list of objects
  objects <- list(...)
  
  # Try to capture names using a more explicit method
  names_list <- sapply(substitute(list(...))[-1], deparse)
  
  # Check if any names are NA or missing, and replace them with default names
  if (any(is.na(names_list) | names_list == "")) {
    empty_names <- which(is.na(names_list) | names_list == "")
    names_list[empty_names] <- paste0("object", empty_names)
  }
  
  # Save each object using its name
  for (i in seq_along(objects)) {
    # Construct file path using the names and save the object
    file_path <- paste0(path_prefix, names_list[i], "_Fold", fold_number, ".rds")
    saveRDS(objects[[i]], file_path)
  }
}

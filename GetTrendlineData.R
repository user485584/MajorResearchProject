# Load necessary libraries
library(dplyr)
library(broom)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)



get_slope_intercept_summary <- function(filtered_data, plot = FALSE, combined = FALSE, plot_stats = FALSE) {
  
  # Get a list of unique diseases from the filtered data
  diseases <- unique(filtered_data$Disease)
  
  # Initialize a list to store slopes, intercepts, and model statistics
  slope_intercept_summary_list <- list()
  model_stats_list <- list()
  
  # Initialize a list to store plots
  plot_list <- list()
  
  # Determine the overall y-axis limits by looking at the min and max TestAUC values across all data
  y_min <- min(filtered_data$TestAUC, na.rm = TRUE)
  y_max <- 1
  
  # Loop over each disease
  for (disease_name in diseases) {
    
    # Filter data for the specific disease
    disease_data <- filtered_data %>%
      filter(Disease == disease_name)
    
    # Fit linear models for each classifier
    lm_glm <- lm(TestAUC ~ NumFeatures, data = disease_data %>% filter(Classifier == "glm"))
    lm_rf <- lm(TestAUC ~ NumFeatures, data = disease_data %>% filter(Classifier == "randomForest"))
    
    # Extract slope, p-value, and R^2 for both models (scientific notation for slope)
    glm_params <- broom::tidy(lm_glm) %>%
      filter(term == "NumFeatures") %>%
      mutate(slope = format(estimate, scientific = TRUE), pval = format(p.value, digits = 2)) %>%
      dplyr::select(slope, pval) %>%
      mutate(Classifier = "glm", Disease = disease_name)
    
    rf_params <- broom::tidy(lm_rf) %>%
      filter(term == "NumFeatures") %>%
      mutate(slope = format(estimate, scientific = TRUE), pval = format(p.value, digits = 2)) %>%
      dplyr::select(slope, pval) %>%
      mutate(Classifier = "randomForest", Disease = disease_name)
    
    glm_stats <- broom::glance(lm_glm) %>%
      mutate(Classifier = "glm", Disease = disease_name)
    
    rf_stats <- broom::glance(lm_rf) %>%
      mutate(Classifier = "randomForest", Disease = disease_name)
    
    # Store model stats in the list
    model_stats_list[[disease_name]] <- bind_rows(glm_stats, rf_stats)
    
    # Calculate summary statistics for each classifier
    glm_summary_stats <- disease_data %>%
      filter(Classifier == "glm") %>%
      summarise(
        Mean_TestAUC = mean(TestAUC, na.rm = TRUE),
        Median_TestAUC = median(TestAUC, na.rm = TRUE),
        Max_TestAUC = max(TestAUC, na.rm = TRUE),
        Min_TestAUC = min(TestAUC, na.rm = TRUE),
        Spread = Max_TestAUC - Min_TestAUC
      ) %>%
      mutate(Classifier = "glm", Disease = disease_name)
    
    rf_summary_stats <- disease_data %>%
      filter(Classifier == "randomForest") %>%
      summarise(
        Mean_TestAUC = mean(TestAUC, na.rm = TRUE),
        Median_TestAUC = median(TestAUC, na.rm = TRUE),
        Max_TestAUC = max(TestAUC, na.rm = TRUE),
        Min_TestAUC = min(TestAUC, na.rm = TRUE),
        Spread = Max_TestAUC - Min_TestAUC
      ) %>%
      mutate(Classifier = "randomForest", Disease = disease_name)
    
    glm_combined <- glm_params %>%
      left_join(glm_summary_stats, by = c("Disease", "Classifier"))
    
    rf_combined <- rf_params %>%
      left_join(rf_summary_stats, by = c("Disease", "Classifier"))
    
    # Combine the parameters for this disease
    slope_intercept_summary <- bind_rows(glm_combined, rf_combined)
    
    # Store the result in the list
    slope_intercept_summary_list[[disease_name]] <- slope_intercept_summary
    
    # Prepare data for plotting
    if (plot) {
      # For plotting, use the data for both classifiers
      disease_plot_data <- disease_data %>%
        filter(Classifier %in% c("glm", "randomForest"))
      
      # Get unique values of NumFeatures
      unique_features <- unique(disease_data$NumFeatures)
      
      # Check if number of unique NumFeatures values is less than 7
      if (length(unique_features) < 7) {
        # Set unique NumFeatures values as x-ticks
        x_breaks <- unique_features
      } else {
        # Use the default ggplot breaks
        x_breaks <- waiver()
      }
      
      # Prepare the text for the slope and p-value (scientific notation for slope)
      glm_text <- paste0("GLM: Slope = ", glm_params$slope,
                         ", p = ", glm_params$pval,
                         ", R² = ", round(glm_stats$r.squared, 2))
      
      rf_text <- paste0("RF: Slope = ", rf_params$slope,
                        ", p = ", rf_params$pval,
                        ", R² = ", round(rf_stats$r.squared, 2))
      
      # Generate the plot with enhanced aesthetics and optional statistics overlay
      p <- ggplot(disease_plot_data, aes(x = NumFeatures, y = TestAUC, color = Classifier, shape = Classifier)) +
        geom_point(size = 4) +  # Larger points for better visibility
        geom_smooth(method = "lm", se = FALSE, size = 1.5) +  # Thicker line for publication readiness
        scale_color_brewer(palette = "Set1") +  # Scientific journal color palette
        theme_bw(base_size = 16) +  # Larger base size for readability
        theme(
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # Larger and centered title
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.position = "top",
          legend.key = element_blank()
        ) +
        labs(title = disease_name,
             x = "Number of Features",
             y = "Test AUC",
             color = "Classifier",
             shape = "Classifier") +
        ylim(y_min, y_max) + # Ensure consistent y-axis range 
        scale_x_continuous(breaks = x_breaks)
      
      # Add slope, p-value, and R-squared text to the bottom of the plot if plot_stats is TRUE
      if (plot_stats) {
        p <- p + annotate("text", x = Inf, y = -Inf, label = glm_text, hjust = 1.1, vjust = -2, size = 4, color = "red") +
          annotate("text", x = Inf, y = -Inf, label = rf_text, hjust = 1.1, vjust = -3.5, size = 4, color = "blue")
      }
      
      # Store the plot in the list
      plot_list[[disease_name]] <- p
    }
  }
  
  # Combine all results into a single dataframe
  final_slope_intercept_summary <- bind_rows(slope_intercept_summary_list) %>%
    dplyr::select(Disease, Classifier, Mean_TestAUC, Median_TestAUC, Min_TestAUC, Max_TestAUC, Spread, slope, pval)
  
  # Combine all model statistics
  final_model_stats <- bind_rows(model_stats_list)
  
  # Return the summary and model stats without reshaping
  if (plot) {
    if (combined) {
      # Combine the plots into a single figure if combined == TRUE
      combined_plot <- ggarrange(plot_list[[diseases[1]]], plot_list[[diseases[2]]], plot_list[[diseases[3]]], 
                                 ncol = 3, nrow = 1, common.legend = TRUE, legend = "top",
                                 labels = c("A", "B", "C"))  # Adding subplot labels for clarity
      
      # Add a unified title for the combined plot
      combined_plot <- annotate_figure(combined_plot,
                                       top = text_grob("Test AUC vs Number of Features for Alzheimers, C9ALS, and SALS",
                                                       face = "bold", size = 20))
      
      return(list(summary = final_slope_intercept_summary, model_stats = final_model_stats, combined_plot = combined_plot))
      
    } else {
      # Return the list of individual plots along with the summary and model stats
      return(list(summary = final_slope_intercept_summary, model_stats = final_model_stats, plots = plot_list))
    }
  } else {
    # Return the summary and model stats dataframes
    return(list(summary = final_slope_intercept_summary, model_stats = final_model_stats))
  }
}

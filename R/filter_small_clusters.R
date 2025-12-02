library(dplyr)

#' Filter Small Cluster Observations from a Data Frame
#'
#' This function removes all observations belonging to clusters (e.g., schools)
#' that contain fewer than a specified minimum number of subjects.
#'
#' @param data The input data frame.
#' @param cluster_id_col The name of the column containing the cluster IDs (e.g., school IDs).
#' @param min_subjects The minimum number of subjects required for a cluster to be retained.
#'
#' @return A new data frame with observations from small clusters removed.
#' @export
#'
#' @examples
#' # Filter out clusters with fewer than 20 subjects
#' data_filtered <- filter_small_clusters(
#'   data = score,
#'   cluster_id_col = "IDSCHOOL",
#'   min_subjects = 20
#' )
#'
#' print(table(data_filtered$IDSCHOOL))
filter_small_clusters <- function(data, cluster_id_col, min_subjects = 20) {
  
  # 1. Input Validation and Standardizing Column Name
  # Use rlang::ensym() to safely capture the column name argument for dplyr functions
  id_col_sym <- rlang::ensym(cluster_id_col)
  
  # Check if the cluster ID column exists
  if (!(cluster_id_col %in% names(data))) {
    stop(paste0("Error: Column '", cluster_id_col, "' not found in the data frame."))
  }
  
  # 2. Calculate Cluster Sizes
  # Use dplyr::count to get a data frame of IDs and their counts
  cluster_counts <- data %>%
    dplyr::count(!!id_col_sym, name = "n_subjects")
  
  # 3. Identify Clusters to Keep
  # Filter to find only the IDs that meet the minimum size requirement
  valid_clusters <- cluster_counts %>%
    dplyr::filter(n_subjects >= min_subjects) %>%
    dplyr::pull(!!id_col_sym) # Extract the vector of valid IDs
  
  # Display Summary Statistics
  cat("\n--- Cluster Size Summary ---\n")
  print(summary(cluster_counts$n_subjects))
  cat("Clusters retained:", length(valid_clusters), "\n")
  cat("Total subjects before filtering:", nrow(data), "\n")
  
  # 4. Filter the Original Data
  # Keep only the rows whose cluster ID is in the vector of valid_clusters
  data_filtered <- data %>%
    dplyr::filter(!!id_col_sym %in% valid_clusters)
  
  cat("Total subjects after filtering:", nrow(data_filtered), "\n")
  cat("--------------------------------\n\n")
  
  return(data_filtered)
}
#' Merge Multiple Simulation Results
#'
#' @description
#' Merge multiple simulation results from \code{drop.simulate()} into a single 
#' combined result. This function intelligently handles results with different 
#' experimental conditions by grouping them appropriately when calculating summary 
#' statistics. Each unique combination of (n, p, graph_type, contamination_scenario) 
#' is treated as a separate experimental condition.
#'
#' @param ... One or more list objects returned by \code{drop.simulate()}. 
#'   Results can have different n, p, graph_type, and contamination_scenario.
#' @param recalculate_summary Logical. Whether to recalculate the summary statistics 
#'   from the combined results. Default is TRUE. When TRUE, summary statistics are 
#'   calculated separately for each unique experimental condition (combination of 
#'   n, p, graph_type, contamination_scenario), ensuring valid comparisons.
#'
#' @return A list with the same structure as \code{drop.simulate()} output, containing:
#' \itemize{
#'   \item \strong{results}: Combined data frame with all iteration results plus a "scenario" column
#'   \item \strong{summary}: Aggregated summary statistics grouped by scenario and method
#'   \item \strong{params}: Parameters from the first simulation result (for reference)
#'   \item \strong{graph_info}: Graph information from the first simulation result (for reference)
#'   \item \strong{scenarios}: A data frame describing each unique scenario
#' }
#'
#' @details
#' This function combines multiple simulation results by:
#' \enumerate{
#'   \item Identifying unique scenarios (combinations of n, p, graph_type, contamination_scenario)
#'   \item Adding a "scenario" column to label each result
#'   \item Merging all iteration-level results into a single data frame
#'   \item Calculating summary statistics separately for each scenario-method combination
#'   \item Creating a scenarios data frame describing all experimental conditions
#' }
#'
#' \strong{Statistical Validity:}
#' \itemize{
#'   \item Summary statistics are calculated \strong{separately} for each unique scenario
#'   \item Methods are only compared within the same experimental conditions
#'   \item This keeps scenario-level summaries internally consistent
#'   \item The scenario column allows easy filtering and analysis
#' }
#'
#' \strong{When to use drop.merge():}
#' \itemize{
#'   \item Combining many different experimental conditions into one comprehensive dataset
#'   \item Running different methods separately for parallel processing
#'   \item Creating a single Excel file with results from multiple experiments
#'   \item Increasing simulation_times by combining multiple runs
#' }
#'
#' @examples
#' \dontrun{
#' # Run multiple simulation scenarios
#' result1 <- drop.simulate(
#'   n = 100, p = 10, 
#'   graph_type = "hub",
#'   contamination_scenario = "cauchy",
#'   simulation_times = 3,
#'   methods = c("DROP", "HUGE_Glasso")
#' )
#' 
#' result2 <- drop.simulate(
#'   n = 200, p = 10, 
#'   graph_type = "hub",
#'   contamination_scenario = "cauchy",
#'   simulation_times = 3,
#'   methods = c("DROP", "HUGE_Glasso")
#' )
#' 
#' result3 <- drop.simulate(
#'   n = 100, p = 10, 
#'   graph_type = "cluster",
#'   contamination_scenario = "leverage",
#'   simulation_times = 3,
#'   methods = c("DROP", "HUGE_Glasso")
#' )
#' 
#' # Merge all results
#' merged_results <- drop.merge(result1, result2, result3)
#' 
#' # Export to Excel
#' drop.excel(merged_results, "merged_experiments.xlsx")
#' }
#'
#' @export
drop.merge <- function(..., recalculate_summary = TRUE) {
  
  # Collect all result objects
  results_list <- list(...)
  
  # Validate we have at least one result
  if (length(results_list) == 0) {
    stop("At least one simulation result must be provided")
  }
  
  # Validate all inputs have the required structure
  for (i in seq_along(results_list)) {
    if (!is.list(results_list[[i]]) || 
        !all(c("results", "summary", "params") %in% names(results_list[[i]]))) {
      stop(sprintf("Result %d must be a list returned by drop.simulate() with components: results, summary, params", i))
    }
  }
  
  cat(sprintf("Merging %d simulation result(s)...\n", length(results_list)))
  
  # Step 1: Identify all unique scenarios and add scenario labels
  all_results <- NULL
  scenario_info <- data.frame(
    scenario = character(),
    n = integer(),
    p = integer(),
    graph_type = character(),
    contamination_scenario = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(results_list)) {
    temp_results <- results_list[[i]]$results
    
    # Get unique scenarios from this result's data (handles multi-scenario simulations)
    unique_scenarios_in_result <- unique(temp_results[, c("n", "p", "graph_type", "contamination_scenario")])
    
    # Process each unique scenario in this result
    for (j in seq_len(nrow(unique_scenarios_in_result))) {
      n_val <- unique_scenarios_in_result$n[j]
      p_val <- unique_scenarios_in_result$p[j]
      graph_val <- unique_scenarios_in_result$graph_type[j]
      contam_val <- unique_scenarios_in_result$contamination_scenario[j]
      
      # Create scenario label for this specific combination
      scenario_label <- sprintf("n%d_p%d_%s_%s", 
                               n_val, p_val, graph_val, contam_val)
      
      # Add to scenario info if new
      if (!any(scenario_info$scenario == scenario_label)) {
        scenario_info <- rbind(scenario_info, data.frame(
          scenario = scenario_label,
          n = n_val,
          p = p_val,
          graph_type = graph_val,
          contamination_scenario = contam_val,
          stringsAsFactors = FALSE
        ))
      }
      
      # Add scenario column to matching rows
      matching_rows <- temp_results$n == n_val & 
                       temp_results$p == p_val & 
                       temp_results$graph_type == graph_val & 
                       temp_results$contamination_scenario == contam_val
      temp_results$scenario[matching_rows] <- scenario_label
    }
    
    # Merge with all_results
    if (is.null(all_results)) {
      all_results <- temp_results
    } else {
      all_results <- rbind(all_results, temp_results)
    }
    
    cat(sprintf("  - Added result %d: %d iterations across %d scenario(s)\n", 
                i, nrow(temp_results), nrow(unique_scenarios_in_result)))
  }
  
  # Count unique scenario-method combinations
  unique_combos <- unique(all_results[, c("scenario", "method_name")])
  cat(sprintf("  - Total scenarios: %d\n", nrow(scenario_info)))
  cat(sprintf("  - Total scenario-method combinations: %d\n", nrow(unique_combos)))
  
  # Step 2: Recalculate summary statistics grouped by scenario
  if (recalculate_summary) {
    cat("  - Calculating summary statistics grouped by scenario and method...\n")
    
    summary_list <- list(
      scenario = character(),
      method_name = character(),
      F1_mean = numeric(),
      F1_sd = numeric(),
      F1_median = numeric(),
      Precision_mean = numeric(),
      Precision_sd = numeric(),
      Recall_mean = numeric(),
      Recall_sd = numeric(),
      Specificity_mean = numeric(),
      Specificity_sd = numeric(),
      MCC_mean = numeric(),
      MCC_sd = numeric(),
      Sparsity_mean = numeric(),
      time_mean = numeric(),
      time_sd = numeric(),
      time_median = numeric(),
      success_rate = numeric()
    )
    
    # Calculate statistics for each scenario-method combination
    for (i in seq_len(nrow(unique_combos))) {
      scenario <- unique_combos$scenario[i]
      method <- unique_combos$method_name[i]
      
      # Filter results for this specific scenario and method
      method_results <- all_results[all_results$scenario == scenario & 
                                    all_results$method_name == method, ]
      
      # Calculate statistics (handling NA values)
      # Note: Statistics are calculated only within the same scenario
      summary_list$scenario <- c(summary_list$scenario, scenario)
      summary_list$method_name <- c(summary_list$method_name, method)
      summary_list$F1_mean <- c(summary_list$F1_mean, mean(method_results$F1, na.rm = TRUE))
      summary_list$F1_sd <- c(summary_list$F1_sd, sd(method_results$F1, na.rm = TRUE))
      summary_list$F1_median <- c(summary_list$F1_median, median(method_results$F1, na.rm = TRUE))
      
      summary_list$Precision_mean <- c(summary_list$Precision_mean, mean(method_results$Precision, na.rm = TRUE))
      summary_list$Precision_sd <- c(summary_list$Precision_sd, sd(method_results$Precision, na.rm = TRUE))
      
      summary_list$Recall_mean <- c(summary_list$Recall_mean, mean(method_results$Recall, na.rm = TRUE))
      summary_list$Recall_sd <- c(summary_list$Recall_sd, sd(method_results$Recall, na.rm = TRUE))
      
      summary_list$Specificity_mean <- c(summary_list$Specificity_mean, mean(method_results$Specificity, na.rm = TRUE))
      summary_list$Specificity_sd <- c(summary_list$Specificity_sd, sd(method_results$Specificity, na.rm = TRUE))
      
      summary_list$MCC_mean <- c(summary_list$MCC_mean, mean(method_results$MCC, na.rm = TRUE))
      summary_list$MCC_sd <- c(summary_list$MCC_sd, sd(method_results$MCC, na.rm = TRUE))
      
      summary_list$Sparsity_mean <- c(summary_list$Sparsity_mean, mean(method_results$Sparsity, na.rm = TRUE))
      
      summary_list$time_mean <- c(summary_list$time_mean, mean(method_results$time, na.rm = TRUE))
      summary_list$time_sd <- c(summary_list$time_sd, sd(method_results$time, na.rm = TRUE))
      summary_list$time_median <- c(summary_list$time_median, median(method_results$time, na.rm = TRUE))
      
      # Calculate success rate (non-NA F1 scores)
      success_rate <- sum(!is.na(method_results$F1)) / nrow(method_results)
      summary_list$success_rate <- c(summary_list$success_rate, success_rate)
    }
  } else {
    # Just combine the summaries from individual results (add scenario column)
    summary_list <- list(scenario = character())
    
    for (i in seq_along(results_list)) {
      temp_summary <- results_list[[i]]$summary
      temp_results <- results_list[[i]]$results
      
      # Get unique scenarios from the results
      # Map each summary row to its scenario
      for (j in seq_len(nrow(temp_summary))) {
        method <- temp_summary$method_name[j]
        contam_scenario <- temp_summary$contamination_scenario[j]
        
        # Find matching params from results
        matching_row <- temp_results[temp_results$method_name == method & 
                                     temp_results$contamination_scenario == contam_scenario, ][1, ]
        
        scenario_label <- sprintf("n%d_p%d_%s_%s", 
                                 matching_row$n, matching_row$p, 
                                 matching_row$graph_type, 
                                 matching_row$contamination_scenario)
        
        # Add scenario label
        summary_list$scenario <- c(summary_list$scenario, scenario_label)
        
        # Append other fields
        for (field in names(temp_summary)) {
          if (field != "scenario") {
            if (!(field %in% names(summary_list))) {
              summary_list[[field]] <- temp_summary[[field]][j]
            } else {
              summary_list[[field]] <- c(summary_list[[field]], temp_summary[[field]][j])
            }
          }
        }
      }
    }
  }
  
  # Convert summary_list to data.frame
  summary_df <- as.data.frame(summary_list, stringsAsFactors = FALSE)
  
  # Add individual columns (n, p, graph_type, contamination_scenario) from scenario info
  if (nrow(summary_df) > 0) {
    summary_df$n <- NA
    summary_df$p <- NA
    summary_df$graph_type <- NA
    summary_df$contamination_scenario <- NA
    
    for (i in seq_len(nrow(summary_df))) {
      scenario_match <- scenario_info[scenario_info$scenario == summary_df$scenario[i], ]
      if (nrow(scenario_match) > 0) {
        summary_df$n[i] <- scenario_match$n[1]
        summary_df$p[i] <- scenario_match$p[1]
        summary_df$graph_type[i] <- scenario_match$graph_type[1]
        summary_df$contamination_scenario[i] <- scenario_match$contamination_scenario[1]
      }
    }
  }
  
  # Keep params and graph_info from first result for reference
  params <- results_list[[1]]$params
  graph_info <- if (!is.null(results_list[[1]]$graph_info)) results_list[[1]]$graph_info else NULL
  
  # Create merged result in the same format as drop.simulate()
  merged_result <- list(
    results = all_results,
    summary = summary_df,
    params = params,  # Reference only - actual params may vary by scenario
    graph_info = graph_info,  # Reference only
    scenarios = scenario_info  # NEW: Information about all scenarios
  )
  
  cat("✓ Merge completed successfully!\n")
  cat(sprintf("  - Summary statistics calculated separately for each scenario\n"))
  cat(sprintf("  - Results maintain statistical validity within each scenario\n"))
  
  return(merged_result)
}

#' Export Simulation Results to Excel
#'
#' @description
#' Export simulation results from \code{drop.simulate()} to an Excel file with 
#' multiple sheets including detailed results and summary statistics.
#'
#' @param sim_result A list object returned by \code{drop.simulate()}.
#' @param filename Character string. The output Excel filename. Default is "simulation_results.xlsx".
#' @param include_details Logical. Whether to include a sheet with all iteration details. Default is TRUE.
#' @param include_summary Logical. Whether to include a sheet with summary statistics. Default is TRUE.

#' @param digits Integer. Number of decimal places for formatting. Default is 4.
#' @param overwrite Logical. Whether to overwrite existing file. Default is TRUE.
#'
#' @return Invisibly returns the path to the created Excel file.
#'
#' @details
#' This function creates an Excel workbook with up to two sheets:
#' \itemize{
#'   \item \strong{DETAILS}: All iteration-level results (if \code{include_details = TRUE})
#'   \item \strong{SUMMARY}: Aggregated statistics including mean, SD, and median (if \code{include_summary = TRUE})
#' }
#'
#' \strong{IMPORTANT:} This function REQUIRES the \code{openxlsx} package to be installed.
#' Install it using: \code{install.packages("openxlsx")}
#'
#' The SUMMARY sheet automatically sorts methods by F1 score (descending) and highlights
#' the best performing method with a green background.
#'
#' @examples
#' \dontrun{
#' # Run simulation
#' result <- drop.simulate(
#'   n = 100, p = 10,
#'   graph_type = "hub",
#'   contamination_scenario = "cauchy",
#'   simulation_times = 3,
#'   methods = c("HUGE_Glasso", "NPN", "Kendall")
#' )
#'
#' # Export to Excel (all sheets)
#' drop.excel(result, "my_results.xlsx")
#'
#' # Export only summary
#' drop.excel(result, "summary_only.xlsx",
#'            include_details = FALSE)
#'
#' # Customize decimal places
#' drop.excel(result, "results_2decimals.xlsx", digits = 2)
#' }
#'
#' @export
drop.excel <- function(
    sim_result,
    filename = "simulation_results.xlsx",
    include_details = TRUE,
    include_summary = TRUE,
    digits = 4,
    overwrite = TRUE
) {
  
  # Check if openxlsx is installed (required dependency)
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("\n*** Package 'openxlsx' is REQUIRED for drop.excel() ***\n\n",
         "Please install it by running:\n",
         "  install.packages('openxlsx')\n\n",
         "This function cannot work without openxlsx package.\n",
         call. = FALSE)
  }
  
  # Validate input
  if (!is.list(sim_result) || 
      !all(c("results", "summary", "params") %in% names(sim_result))) {
    stop("sim_result must be a list returned by drop.simulate() with components: results, summary, params")
  }
  
  # Extract components
  results_df <- sim_result$results
  summary_list <- sim_result$summary
  
  # Create workbook
  wb <- openxlsx::createWorkbook()
  
  cat(sprintf("Creating Excel file: %s\n", filename))
  
  # Helper function to format numeric values, handling NA
  format_numeric <- function(x, digits) {
    sapply(x, function(val) {
      if (is.na(val)) {
        "NA"
      } else {
        sprintf(paste0("%.", digits, "f"), val)
      }
    })
  }
  
  # ============================================================
  # Sheet 1: DETAILS (all iterations)
  # ============================================================
  if (include_details && !is.null(results_df) && nrow(results_df) > 0) {
    cat("  - Adding DETAILS sheet...\n")
    
    openxlsx::addWorksheet(wb, "DETAILS")
    
    # Format the details data - handle NA values properly
    # Now results_df already contains n, p, graph_type, contamination_scenario, simulation_times
    details_export <- data.frame(
      n = results_df$n,
      p = results_df$p,
      graph_type = results_df$graph_type,
      contamination_scenario = results_df$contamination_scenario,
      simulation_times = results_df$simulation_times,
      Iteration = results_df$iteration,
      Method = results_df$method_name,
      F1 = format_numeric(results_df$F1, digits),
      Precision = format_numeric(results_df$Precision, digits),
      Recall = format_numeric(results_df$Recall, digits),
      Specificity = format_numeric(results_df$Specificity, digits),
      MCC = format_numeric(results_df$MCC, digits),
      Sparsity = format_numeric(results_df$Sparsity, digits),
      True_Sparsity = format_numeric(results_df$True_Sparsity, digits),
      Time_sec = format_numeric(results_df$time, digits),
      stringsAsFactors = FALSE
    )
    
    openxlsx::writeData(wb, "DETAILS", details_export)
    
    # Format header
    header_style <- openxlsx::createStyle(
      textDecoration = "bold",
      fgFill = "#4F81BD",
      fontColour = "#FFFFFF",
      border = "TopBottomLeftRight"
    )
    openxlsx::addStyle(wb, "DETAILS", header_style, rows = 1,
               cols = seq_len(ncol(details_export)), gridExpand = TRUE)
    
    # Set column widths (specific widths for better display)
    col_widths <- c(
      12,  # n
      12,  # p
      18,  # graph_type
      25,  # contamination_scenario
      18,  # simulation_times
      12,  # Iteration
      15,  # Method
      12,  # F1
      12,  # Precision
      12,  # Recall
      14,  # Specificity
      12,  # MCC
      12,  # Sparsity
      14,  # True_Sparsity
      12   # Time_sec
    )
    openxlsx::setColWidths(wb, "DETAILS", cols = seq_len(ncol(details_export)), widths = col_widths)
    
    # Freeze first row
    openxlsx::freezePane(wb, "DETAILS", firstRow = TRUE)
  }
  
  # ============================================================
  # Sheet 2: SUMMARY (aggregated statistics)
  # ============================================================
  if (include_summary && !is.null(summary_list) && length(summary_list) > 0) {
    cat("  - Adding SUMMARY sheet...\n")
    
    openxlsx::addWorksheet(wb, "SUMMARY")
    
    # Build summary from results_df which now contains all parameter columns
    # For merged results, we need to handle all scenario-method combinations
    if (!is.null(results_df) && nrow(results_df) > 0) {
      summary_export_list <- list()
      
      # Determine number of summary rows
      if (is.data.frame(summary_list)) {
        n_summary_rows <- nrow(summary_list)
      } else if (is.list(summary_list) && !is.null(summary_list$scenario)) {
        n_summary_rows <- length(summary_list$scenario)
      } else {
        n_summary_rows <- 0
      }
      
      # For each summary row, find corresponding parameter info from results_df
      for (i in seq_len(n_summary_rows)) {
        # Get method name from summary
        if (is.data.frame(summary_list)) {
          method <- summary_list$method_name[i]
          scenario <- summary_list$contamination_scenario[i]
          summary_row <- summary_list[i, ]
          
          # For merged results, get params directly from summary if available
          if ("n" %in% colnames(summary_list) && "p" %in% colnames(summary_list) &&
              "graph_type" %in% colnames(summary_list)) {
            param_n <- summary_list$n[i]
            param_p <- summary_list$p[i]
            param_graph_type <- summary_list$graph_type[i]
            param_contamination_scenario <- summary_list$contamination_scenario[i]
            use_summary_params <- TRUE
          } else {
            use_summary_params <- FALSE
          }
        } else {
          # summary_list is a list
          method <- summary_list$method_name[i]
          scenario <- summary_list$contamination_scenario[i]
          use_summary_params <- FALSE
          summary_row <- list(
            F1_mean = summary_list$F1_mean[i],
            F1_sd = summary_list$F1_sd[i],
            F1_median = summary_list$F1_median[i],
            Precision_mean = summary_list$Precision_mean[i],
            Precision_sd = summary_list$Precision_sd[i],
            Recall_mean = summary_list$Recall_mean[i],
            Recall_sd = summary_list$Recall_sd[i],
            Specificity_mean = summary_list$Specificity_mean[i],
            Specificity_sd = summary_list$Specificity_sd[i],
            MCC_mean = summary_list$MCC_mean[i],
            MCC_sd = summary_list$MCC_sd[i],
            Sparsity_mean = summary_list$Sparsity_mean[i],
            time_mean = summary_list$time_mean[i],
            time_sd = summary_list$time_sd[i],
            time_median = summary_list$time_median[i],
            success_rate = summary_list$success_rate[i]
          )
        }
        
        # Find parameter info from results_df or use summary params
        if (use_summary_params) {
          # Use parameters directly from summary (merged results case)
          param_row <- list(
            n = param_n,
            p = param_p,
            graph_type = param_graph_type,
            contamination_scenario = param_contamination_scenario,
            simulation_times = NA  # Not critical for summary
          )
        } else {
          # Match by method AND scenario to get correct parameters for multi-scenario merged results
          if (!is.null(scenario) && "contamination_scenario" %in% colnames(results_df)) {
            # For merged results with multiple (n,p,graph_type) combinations,
            # also match by n, p, graph_type if available in summary
            if (is.data.frame(summary_list) && 
                all(c("n", "p", "graph_type") %in% colnames(summary_list))) {
              # Match by n, p, graph_type, contamination_scenario AND method
              matching_rows <- which(
                results_df$n == summary_list$n[i] &
                results_df$p == summary_list$p[i] &
                results_df$graph_type == summary_list$graph_type[i] &
                results_df$contamination_scenario == scenario & 
                results_df$method_name == method
              )
            } else {
              # Fallback: match by scenario and method only
              matching_rows <- which(results_df$contamination_scenario == scenario & 
                                   results_df$method_name == method)
            }
          } else {
            # Single scenario case, just match by method
            matching_rows <- which(results_df$method_name == method)
          }
          
          if (length(matching_rows) > 0) {
            param_row <- results_df[matching_rows[1], ]
          } else {
            next  # Skip if no matching row found in results
          }
        }
        
        # Combine parameter info with summary statistics
        summary_export_list[[i]] <- data.frame(
          n = param_row$n,
          p = param_row$p,
          graph_type = param_row$graph_type,
          contamination_scenario = param_row$contamination_scenario,
          simulation_times = param_row$simulation_times,
          Method = method,
          F1_mean = format_numeric(summary_row$F1_mean, digits),
          F1_sd = format_numeric(summary_row$F1_sd, digits),
          F1_median = format_numeric(summary_row$F1_median, digits),
          Precision_mean = format_numeric(summary_row$Precision_mean, digits),
          Precision_sd = format_numeric(summary_row$Precision_sd, digits),
          Recall_mean = format_numeric(summary_row$Recall_mean, digits),
          Recall_sd = format_numeric(summary_row$Recall_sd, digits),
          Specificity_mean = format_numeric(summary_row$Specificity_mean, digits),
          Specificity_sd = format_numeric(summary_row$Specificity_sd, digits),
          MCC_mean = format_numeric(summary_row$MCC_mean, digits),
          MCC_sd = format_numeric(summary_row$MCC_sd, digits),
          Sparsity_mean = format_numeric(summary_row$Sparsity_mean, digits),
          Time_mean = format_numeric(summary_row$time_mean, digits),
          Time_sd = format_numeric(summary_row$time_sd, digits),
          Time_median = format_numeric(summary_row$time_median, digits),
          Success_rate = sprintf("%.2f%%", summary_row$success_rate * 100),
          stringsAsFactors = FALSE
        )
      }
      
      summary_export <- do.call(rbind, summary_export_list)
    } else {
      # Fallback: no results_df, should not happen normally
      stop("Cannot create SUMMARY sheet: results data frame is missing")
    }
    
    # Sort by n, p, graph_type, contamination_scenario, then F1_mean (descending)
    f1_values <- suppressWarnings(as.numeric(summary_export$F1_mean))
    f1_values[is.na(f1_values)] <- -Inf  # Put NA at the end
    
    summary_export <- summary_export[order(summary_export$n, summary_export$p, 
                                           summary_export$graph_type, 
                                           summary_export$contamination_scenario, 
                                           -f1_values), ]
    
    openxlsx::writeData(wb, "SUMMARY", summary_export)
    
    # Format header
    header_style <- openxlsx::createStyle(
      textDecoration = "bold",
      fgFill = "#4F81BD",
      fontColour = "#FFFFFF",
      border = "TopBottomLeftRight"
    )
    openxlsx::addStyle(wb, "SUMMARY", header_style, rows = 1,
               cols = seq_len(ncol(summary_export)), gridExpand = TRUE)
    
    # Highlight best method for each scenario (highest F1)
    if (nrow(summary_export) > 0) {
      best_style <- openxlsx::createStyle(
        fgFill = "#C6EFCE",
        border = "TopBottomLeftRight"
      )
      
      # Highlight best method within each scenario
      # Create scenario identifier
      scenario_ids <- paste(summary_export$n, summary_export$p, 
                           summary_export$graph_type, 
                           summary_export$contamination_scenario, sep = "_")
      unique_scenarios <- unique(scenario_ids)
      
      for (scenario in unique_scenarios) {
        scenario_rows <- which(scenario_ids == scenario)
        if (length(scenario_rows) > 0) {
          # First row of each scenario is the best (due to sorting)
          best_row <- scenario_rows[1]
          openxlsx::addStyle(wb, "SUMMARY", best_style, rows = best_row + 1,
                             cols = seq_len(ncol(summary_export)), gridExpand = TRUE)
        }
      }
    }
    
    # Set column widths (specific widths for better display)
    summary_col_widths <- c(
      12,  # n
      12,  # p
      18,  # graph_type
      25,  # contamination_scenario
      18,  # simulation_times
      15,  # Method
      12,  # F1_mean
      12,  # F1_sd
      12,  # F1_median
      15,  # Precision_mean
      12,  # Precision_sd
      12,  # Recall_mean
      12,  # Recall_sd
      16,  # Specificity_mean
      14,  # Specificity_sd
      12,  # MCC_mean
      12,  # MCC_sd
      14,  # Sparsity_mean
      12,  # Time_mean
      12,  # Time_sd
      13,  # Time_median
      14   # Success_rate
    )
    openxlsx::setColWidths(wb, "SUMMARY", cols = seq_len(ncol(summary_export)), widths = summary_col_widths)
    
    # Freeze first row
    openxlsx::freezePane(wb, "SUMMARY", firstRow = TRUE)
  }
  
  # ============================================================
  # Save workbook
  # ============================================================
  openxlsx::saveWorkbook(wb, filename, overwrite = overwrite)
  
  cat(sprintf("✓ Excel file created successfully: %s\n", filename))
  cat(sprintf("  File size: %.1f KB\n", file.info(filename)$size / 1024))
  
  invisible(filename)
}

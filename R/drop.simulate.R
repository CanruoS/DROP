#' Run GGM Simulation Experiment
#'
#' Run simulation experiments for Gaussian Graphical Model estimation with
#' multiple methods, graph types, and contamination scenarios.
#'
#' @param n Integer. Sample size (number of observations).
#' @param p Integer. Dimension (number of variables).
#' @param graph_type Character. Type of graph structure. One of: "band", "hub", 
#'   "cluster", "random", "scale-free".
#' @param contamination_scenario Character or character vector. Type(s) of contamination. 
#'   Can be one or more of: "clean", "cauchy", "leverage". If multiple scenarios are provided,
#'   they will be tested sequentially. Default is "clean".
#' @param contamination_rate Numeric. Proportion of contaminated observations. Default is 0.1.
#' @param skip_failed_methods Logical. When testing multiple scenarios, if TRUE (default), methods 
#'   that fail completely on "clean" scenario (all iterations timeout/error) will be skipped in 
#'   subsequent scenarios to save computation time. If FALSE, all methods will be tested on all 
#'   scenarios regardless of failures. Only applies when multiple scenarios are specified.
#' @param save_matrices Logical. If TRUE, saves the estimated precision matrices and adjacency matrices 
#'   for each method/iteration/scenario. Default is FALSE (to save memory). When enabled, results will 
#'   include \code{estimated_matrices} containing \code{theta} (precision matrices) and \code{adjacency} 
#'   (adjacency matrices) for each combination.
#' @param g Integer. Graph parameter passed to drop.generator(). See \code{?drop.generator} for details. Default is NULL.
#' @param prob Numeric. Connection probability passed to drop.generator(). See \code{?drop.generator} for details. Default is NULL.
#' @param simulation_times Integer. Number of simulation replications. Default is 3.
#' @param methods Either a character vector of pre-defined method names, or a named
#'   list of custom estimation functions. 
#'   
#'   Pre-defined methods: "DROP", "HUGE_Glasso", "HUGE_MB", "NPN", "Scaled_Lasso",
#'   "SCIO", "SCAD", "MCP", "TIGER", "CLIME", "Kendall", "Spearman".
#'   
#'   For custom methods: each list element should be a function that takes a data 
#'   matrix X (n x p) as input and returns a precision matrix estimate (p x p). 
#'   List names will be used as method labels in output.
#' @param seed Integer. Random seed for reproducibility. Default is NULL.
#' @param verbose Logical. If TRUE, print progress information. Default is TRUE.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{results}{Data frame with detailed results for each iteration and method,
#'     including: iteration, method_name, contamination_scenario, F1, Precision, Recall, 
#'     Specificity, MCC, Sparsity, time.}
#'   \item{summary}{Data frame with aggregated statistics (mean, sd, median) for 
#'     each unique combination of method and contamination_scenario. Each row represents
#'     one (method_name, contamination_scenario) combination.}
#'   \item{true_theta}{The true precision matrix used in the simulation.}
#'   \item{graph_info}{List containing graph structure information (adjacency matrix, etc.).}
#'   \item{estimated_matrices}{(Only if save_matrices=TRUE) List with two components: 
#'     \code{theta} and \code{adjacency}, each containing estimated matrices indexed by 
#'     method_name, scenario, and iteration.}
#'   \item{params}{List of simulation parameters (n, p, graph_type, contamination_scenario, 
#'     simulation_times, seed, save_matrices).}
#' }
#'
#' @details
#' This function conducts simulation experiments to evaluate the performance of 
#' multiple precision matrix estimation methods under various conditions.
#' 
#' For each simulation iteration:
#' \enumerate{
#'   \item Generate a graph structure using \code{drop.generator}
#'   \item Generate contaminated data based on specified scenario
#'   \item Apply each estimation method to the data
#'   \item Evaluate performance using \code{drop.evaluate}
#' }
#' 
#' The \code{methods} parameter can be:
#' \itemize{
#'   \item A character vector of pre-defined method names (e.g., c("DROP", "HUGE_Glasso"))
#'   \item A named list of custom functions with signature \code{function(X) { ... }}
#' }
#' 
#' Pre-defined methods require their respective packages to be installed:
#' \itemize{
#'   \item huge: HUGE_Glasso, HUGE_MB, NPN, Kendall, Spearman
#'   \item scalreg: Scaled_Lasso
#'   \item scio: SCIO
#'   \item GGMncv: SCAD, MCP
#'   \item flare: TIGER, CLIME
#' }
#' 
#' \strong{Timeout control:} To prevent excessively long computations, all methods 
#' except "DROP" are subject to a 60-second timeout per iteration. If a method 
#' exceeds this limit, it will be terminated and recorded as NA in the results.
#'
#' @examples
#' \dontrun{
#' # Method 1: Use pre-defined methods (character vector)
#' sim_result <- drop.simulate(
#'   n = 100,
#'   p = 10,
#'   graph_type = "hub",
#'   contamination_scenario = "cauchy",
#'   contamination_rate = 0.1,
#'   simulation_times = 3,
#'   methods = c("DROP", "HUGE_Glasso", "NPN"),
#'   seed = 123,
#'   verbose = TRUE
#' )
#' 
#' # Method 2: Use custom methods (named list)
#' custom_methods <- list(
#'   DROP = function(X) {
#'     fit <- drop(X, verbose = FALSE)
#'     return(fit$est_K)
#'   },
#'   Sample_Cov = function(X) {
#'     S <- cov(X)
#'     K <- tryCatch(solve(S), error = function(e) diag(ncol(X)))
#'     return(K)
#'   }
#' )
#' 
#' sim_result <- drop.simulate(
#'   n = 100,
#'   p = 10,
#'   graph_type = "cluster",
#'   simulation_times = 3,
#'   methods = custom_methods,
#'   seed = 456
#' )
#' 
#' # View summary statistics
#' print(sim_result$summary)
#' 
#' # Access detailed results
#' head(sim_result$results)
#' 
#' # Get top methods by F1 score
#' top_methods <- sim_result$summary[order(-sim_result$summary$F1_mean), ]
#' head(top_methods)
#' 
#' # Method 3: Specify custom graph parameters
#' sim_result <- drop.simulate(
#'   n = 100,
#'   p = 10,
#'   graph_type = "hub",
#'   g = 10,  # Use 10 hubs instead of default
#'   contamination_scenario = "clean",
#'   simulation_times = 3,
#'   methods = c("DROP", "HUGE_Glasso"),
#'   seed = 789
#' )
#' 
#' # Method 4: Test multiple contamination scenarios
#' # By default, methods that fail on clean are skipped in later scenarios
#' sim_result <- drop.simulate(
#'   n = 100,
#'   p = 10,
#'   graph_type = "hub",
#'   contamination_scenario = c("clean", "cauchy", "leverage"),
#'   contamination_rate = 0.1,
#'   simulation_times = 3,
#'   methods = c("DROP", "HUGE_Glasso"),
#'   skip_failed_methods = TRUE,  # Default: skip failed methods
#'   seed = 999
#' )
#' 
#' # To test all methods on all scenarios regardless of failures:
#' sim_result_full <- drop.simulate(
#'   n = 100,
#'   p = 10,
#'   graph_type = "hub",
#'   contamination_scenario = c("clean", "cauchy", "leverage"),
#'   contamination_rate = 0.1,
#'   simulation_times = 3,
#'   methods = c("DROP", "HUGE_Glasso"),
#'   skip_failed_methods = FALSE,  # Test all on all scenarios
#'   seed = 999
#' )
#' 
#' # Results include all scenarios
#' table(sim_result$results$contamination_scenario)
#' 
#' # Method 5: Save estimated matrices for further analysis
#' sim_with_matrices <- drop.simulate(
#'   n = 100,
#'   p = 10,
#'   graph_type = "hub",
#'   contamination_scenario = "clean",
#'   simulation_times = 2,
#'   methods = c("DROP"),
#'   save_matrices = TRUE,  # Save all estimated matrices
#'   seed = 888
#' )
#' 
#' # Access estimated precision matrices
#' # Format: method_scenario_iteration
#' drop_iter1 <- sim_with_matrices$estimated_matrices$theta$DROP_scenario1_iter1
#' drop_adj1 <- sim_with_matrices$estimated_matrices$adjacency$DROP_scenario1_iter1
#' 
#' # Compare with true matrix
#' true_matrix <- sim_with_matrices$true_theta
#' )
#' }
#'
#' @export
drop.simulate <- function(
    n,
    p,
    graph_type = c("band", "hub", "cluster", "random", "scale-free"),
    contamination_scenario = c("clean", "cauchy", "leverage"),
    contamination_rate = 0.1,
    g = NULL,
    prob = NULL,
    simulation_times = 3,
    methods,
    skip_failed_methods = TRUE,
    save_matrices = FALSE,
    seed = NULL,
    verbose = TRUE
) {
  
  # Validate inputs
  graph_type <- match.arg(graph_type)
  
  # Validate contamination scenarios (allow vector input)
  valid_scenarios <- c("clean", "cauchy", "leverage")
  if (!all(contamination_scenario %in% valid_scenarios)) {
    invalid <- setdiff(contamination_scenario, valid_scenarios)
    stop(sprintf("Invalid contamination_scenario: %s. Must be one or more of: %s",
                 paste(invalid, collapse=", "), 
                 paste(valid_scenarios, collapse=", ")))
  }
  
  # Handle methods parameter - can be character vector or named list
  if (is.character(methods)) {
    # Pre-defined methods requested
    method_names <- methods
    
    # Define all available pre-defined methods
    available_methods <- list(
      
      DROP = function(X) {
        fit <- drop(X, verbose = FALSE)
        return(fit$est_K)
      },
      
      HUGE_Glasso = function(X) {
        if (!requireNamespace("huge", quietly = TRUE)) {
          stop("Package 'huge' is required for HUGE_Glasso method")
        }
        data_huge <- huge::huge(X, method = "glasso", verbose = FALSE)
        select <- huge::huge.select(data_huge, criterion = "stars", verbose = FALSE)
        return(select$refit)
      },
      
      HUGE_MB = function(X) {
        if (!requireNamespace("huge", quietly = TRUE)) {
          stop("Package 'huge' is required for HUGE_MB method")
        }
        data_huge <- huge::huge(X, method = "mb", verbose = FALSE)
        select <- huge::huge.select(data_huge, verbose = FALSE)
        return(as.matrix(select$refit))
      },
      
      NPN = function(X) {
        if (!requireNamespace("huge", quietly = TRUE)) {
          stop("Package 'huge' is required for NPN method")
        }
        X_npn <- huge::huge.npn(X, verbose = FALSE)
        data_huge <- huge::huge(X_npn, method = "glasso", verbose = FALSE)
        select <- huge::huge.select(data_huge, criterion = "stars", verbose = FALSE)
        return(select$refit)
      },
      
      Scaled_Lasso = function(X) {
        if (!requireNamespace("scalreg", quietly = TRUE)) {
          stop("Package 'scalreg' is required for Scaled_Lasso method")
        }
        obj <- scalreg::scalreg(X = X, y = NULL, lam0 = "quantile", LSE = FALSE)
        return(obj$precision)
      },
      
      SCIO = function(X) {
        if (!requireNamespace("scio", quietly = TRUE)) {
          stop("Package 'scio' is required for SCIO method")
        }
        scio_fit <- scio::scio.cv(X)
        return(as.matrix(scio_fit$w))
      },
      
      SCAD = function(X) {
        if (!requireNamespace("GGMncv", quietly = TRUE)) {
          stop("Package 'GGMncv' is required for SCAD method")
        }
        corr <- cor(X)
        scad_fit <- GGMncv::ggmncv(corr, n = nrow(X), penalty = "scad")
        return(as.matrix(scad_fit$Theta))
      },
      
      MCP = function(X) {
        if (!requireNamespace("GGMncv", quietly = TRUE)) {
          stop("Package 'GGMncv' is required for MCP method")
        }
        corr <- cor(X)
        mcp_fit <- GGMncv::ggmncv(corr, n = nrow(X), penalty = "mcp")
        return(as.matrix(mcp_fit$Theta))
      },
      
      TIGER = function(X) {
        if (!requireNamespace("flare", quietly = TRUE)) {
          stop("Package 'flare' is required for TIGER method")
        }
        re <- flare::sugm(X, method = "tiger", verbose = FALSE)
        tiger_select <- flare::sugm.select(est = re, criterion = "stars", verbose = FALSE)
        return(as.matrix(tiger_select$opt.icov))
      },
      
      CLIME = function(X) {
        if (!requireNamespace("flare", quietly = TRUE)) {
          stop("Package 'flare' is required for CLIME method")
        }
        re <- flare::sugm(X, method = "clime", verbose = FALSE)
        clime_select <- flare::sugm.select(est = re, criterion = "cv", verbose = FALSE)
        return(as.matrix(clime_select$opt.icov))
      },
      
      Kendall = function(X) {
        if (!requireNamespace("huge", quietly = TRUE)) {
          stop("Package 'huge' is required for Kendall method")
        }
        if (!requireNamespace("Matrix", quietly = TRUE)) {
          stop("Package 'Matrix' is required for Kendall method")
        }
        n <- nrow(X)
        p <- ncol(X)
        
        # Kendall tau correlation
        tau_matrix <- cor(X, method = "kendall")
        
        # Transform according to formula: R_ij^τ = sin(π/2 * τ_ij)
        S_kendall <- sin(pi / 2 * tau_matrix)
        
        # Ensure positive definiteness using nearPD
        S_kendall_pd <- as.matrix(Matrix::nearPD(S_kendall, corr = TRUE)$mat)
        
        # Use huge with transformed correlation
        X_pseudo <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = S_kendall_pd)
        data_huge <- huge::huge(X_pseudo, method = "glasso", verbose = FALSE)
        select <- huge::huge.select(data_huge, criterion = "stars", verbose = FALSE)
        return(select$refit)
      },
      
      Spearman = function(X) {
        if (!requireNamespace("huge", quietly = TRUE)) {
          stop("Package 'huge' is required for Spearman method")
        }
        if (!requireNamespace("Matrix", quietly = TRUE)) {
          stop("Package 'Matrix' is required for Spearman method")
        }
        n <- nrow(X)
        p <- ncol(X)
        
        # Spearman correlation
        rho_matrix <- cor(X, method = "spearman")
        
        # Transform according to formula: R_ij^ρ = 2*sin(π/6 * ρ_ij)
        S_spearman <- 2 * sin(pi / 6 * rho_matrix)
        
        # Ensure positive definiteness using nearPD
        S_spearman_pd <- as.matrix(Matrix::nearPD(S_spearman, corr = TRUE)$mat)
        
        # Use huge with transformed correlation
        X_pseudo <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = S_spearman_pd)
        data_huge <- huge::huge(X_pseudo, method = "glasso", verbose = FALSE)
        select <- huge::huge.select(data_huge, criterion = "stars", verbose = FALSE)
        return(select$refit)
      }
    )
    
    # Check if requested methods are available
    invalid <- setdiff(method_names, names(available_methods))
    if (length(invalid) > 0) {
      stop(sprintf("Invalid method names: %s\nAvailable methods: %s",
                   paste(invalid, collapse = ", "),
                   paste(names(available_methods), collapse = ", ")))
    }
    
    # Extract requested methods
    methods <- available_methods[method_names]
    
  } else if (is.list(methods)) {
    # Custom methods provided
    if (is.null(names(methods))) {
      stop("methods must be a named list when providing custom functions")
    }
  } else {
    stop("methods must be either a character vector or a named list of functions")
  }
  
  if (contamination_rate < 0 || contamination_rate > 1) {
    stop("contamination_rate must be between 0 and 1")
  }
  
  if (simulation_times < 1) {
    stop("simulation_times must be at least 1")
  }
  
  if (verbose) {
    cat("\n========================================\n")
    cat("  DROP Simulation Experiment\n")
    cat("========================================\n")
    if (!is.null(seed)) cat(sprintf("Seed: %d\n", seed))
    cat(sprintf("Dimension (p): %d\n", p))
    cat(sprintf("Sample size (n): %d\n", n))
    cat(sprintf("Graph type: %s\n", graph_type))
    
    # Display contamination scenario(s)
    scenario_str <- paste(contamination_scenario, collapse = ", ")
    cat(sprintf("Contamination: %s", scenario_str))
    if (any(contamination_scenario != "clean")) {
      cat(sprintf(" (rate=%.1f%%)", contamination_rate * 100))
    }
    cat("\n")
    cat(sprintf("Replications: %d\n", simulation_times))
    cat(sprintf("Methods: %d\n", length(methods)))
    cat("========================================\n\n")
  }
  
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Generate graph structure (only once, reuse for all iterations)
  if (verbose) {
    cat("Generating graph structure...\n")
  }
  
  # Build arguments for drop.generator
  gen_args <- list(
    graph_type = graph_type,
    n = n,
    p = p,
    contamination_scenario = "clean",  # Generate clean data first
    g = g,
    prob = prob,
    seed = if (!is.null(seed)) seed else NULL
  )
  
  # Generate graph structure
  graph_data <- do.call(drop.generator, gen_args)
  
  true_theta <- graph_data$omega  # Precision matrix (Omega)
  true_adj <- graph_data$theta    # Adjacency matrix (Theta)
  
  # Save Sigma for data generation in iterations
  true_sigma <- as.matrix(graph_data$sigma) 
  
  # Calculate true sparsity
  true_edges <- sum(true_adj[upper.tri(true_adj)] != 0)
  total_possible_edges <- p * (p - 1) / 2
  true_sparsity <- 1 - (true_edges / total_possible_edges)
  
  if (verbose) {
    cat(sprintf("True edges: %d / %d\n", true_edges, total_possible_edges))
    cat(sprintf("True sparsity: %.4f\n\n", true_sparsity))
    cat("Starting simulation...\n")
  }
  
  # Storage for results
  results_list <- list()
  result_counter <- 1
  
  # Storage for estimated matrices (if save_matrices=TRUE)
  if (save_matrices) {
    estimated_theta_list <- list()
    estimated_adj_list <- list()
  }
  
  # Track which methods failed on clean scenario (for early stopping)
  failed_methods <- character(0)
  
  # Outer loop over contamination scenarios
  for (scenario_idx in seq_along(contamination_scenario)) {
    current_scenario <- contamination_scenario[scenario_idx]
    
    if (verbose && length(contamination_scenario) > 1) {
      cat(sprintf("\n--- Scenario %d/%d: %s ---\n", 
                  scenario_idx, length(contamination_scenario), current_scenario))
    }
    
    # Generate NA results for methods that failed in clean scenario (if skip_failed_methods=TRUE)
    if (skip_failed_methods && current_scenario != "clean" && length(failed_methods) > 0) {
      if (verbose && any(names(methods) %in% failed_methods)) {
        skipped <- names(methods)[names(methods) %in% failed_methods]
        cat(sprintf("  Skipping methods that failed on clean: %s\n", 
                    paste(skipped, collapse=", ")))
      }
      
      # Generate NA results for failed methods
      for (iteration in 1:simulation_times) {
        for (method_name in failed_methods) {
          if (method_name %in% names(methods)) {
            results_list[[result_counter]] <- data.frame(
              n = n,
              p = p,
              graph_type = graph_type,
              contamination_scenario = current_scenario,
              simulation_times = simulation_times,
              iteration = iteration,
              method_name = method_name,
              F1 = NA,
              Precision = NA,
              Recall = NA,
              Specificity = NA,
              MCC = NA,
              Sparsity = NA,
              True_Sparsity = true_sparsity,
              time = NA,
              stringsAsFactors = FALSE
            )
            result_counter <- result_counter + 1
          }
        }
      }
    }
    
    # Iteration loop for current scenario
    for (iteration in 1:simulation_times) {
      if (verbose && (iteration %% 5 == 0 || iteration == 1)) {
        cat(sprintf("  Progress: %d/%d\n", iteration, simulation_times))
        flush.console()
      }
      
      # Generate data for this iteration (reuse same graph structure)
      # Set seed for this iteration
      iter_seed <- if (!is.null(seed)) seed + iteration + (scenario_idx - 1) * 1000 else NULL
      if (!is.null(iter_seed)) set.seed(iter_seed)
      
      # Generate clean data from the same graph
      X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = true_sigma)
      
      # Add contamination if needed
      if (current_scenario == "cauchy") {
        n_cont <- floor(n * contamination_rate)
        if (n_cont > 0) {
          cont_idx <- sample(1:n, n_cont)
          for (i in cont_idx) {
            X[i, ] <- X[i, ] + rcauchy(p, location = 0, scale = 5)
          }
        }
      } else if (current_scenario == "leverage") {
        n_cont <- floor(n * contamination_rate)
        if (n_cont > 0) {
          cont_idx <- sample(1:n, n_cont)
          leverage_cov <- true_sigma * 100
          X[cont_idx, ] <- MASS::mvrnorm(n_cont, mu = rep(0, p), Sigma = leverage_cov)
        }
      }
    
    # Run each method (skip if it failed in clean scenario and skip_failed_methods=TRUE)
    for (method_name in names(methods)) {
      # Skip this method if skip_failed_methods enabled and it failed in clean
      if (skip_failed_methods && current_scenario != "clean" && method_name %in% failed_methods) {
        next
      }
      method_func <- methods[[method_name]]
      
      # Run method and time it
      start_time <- proc.time()
      timeout_reached <- FALSE
      
      # Set timeout for non-DROP methods (60 seconds)
      if (method_name != "DROP") {
        setTimeLimit(cpu = Inf, elapsed = 60, transient = TRUE)
      }
      
      est_theta <- tryCatch({
        method_func(X)
      }, warning = function(w) {
        # Handle warnings (including timeout warnings)
        if (grepl("reached elapsed time limit|reached CPU time limit", w$message, ignore.case = TRUE)) {
          timeout_reached <<- TRUE
          if (verbose) {
            cat(sprintf("    Warning: %s exceeded 60 sec timeout at iteration %d\n", 
                        method_name, iteration))
          }
        }
        NULL
      }, error = function(e) {
        # Handle errors (including timeout errors)
        if (grepl("reached elapsed time limit|reached CPU time limit", e$message, ignore.case = TRUE)) {
          timeout_reached <<- TRUE
          if (verbose) {
            cat(sprintf("    Warning: %s exceeded 60 sec timeout at iteration %d\n", 
                        method_name, iteration))
          }
        } else if (verbose) {
          cat(sprintf("    Warning: %s failed at iteration %d: %s\n", 
                      method_name, iteration, e$message))
        }
        NULL
      }, finally = {
        # Always reset time limit, even if error occurred
        tryCatch({
          setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
        }, error = function(e) {
          # Silently ignore errors when resetting time limit
        })
      })
      
      end_time <- proc.time()
      elapsed_time <- (end_time - start_time)[3]
      
      # If timeout, set elapsed_time to 60
      if (timeout_reached) {
        elapsed_time <- 60
      }
      
      # Evaluate if method succeeded
      if (!is.null(est_theta) && is.matrix(est_theta)) {
        # Ensure dimensions match
        if (nrow(est_theta) != p || ncol(est_theta) != p) {
          if (verbose) {
            cat(sprintf("    Warning: %s returned wrong dimensions at iteration %d\n",
                        method_name, iteration))
          }
          est_theta <- NULL
        }
      }
      
      if (!is.null(est_theta)) {
        # Evaluate performance
        eval_result <- drop.evaluate(est_theta, true_theta, threshold = 1e-6)
        
        # Save estimated matrix if requested
        if (save_matrices) {
          matrix_key <- sprintf("%s_scenario%d_iter%d", method_name, scenario_idx, iteration)
          estimated_theta_list[[matrix_key]] <- est_theta
          # Calculate adjacency matrix from precision matrix
          est_adj <- abs(est_theta) > 1e-6
          diag(est_adj) <- 0
          estimated_adj_list[[matrix_key]] <- est_adj
        }
        
        # Store results
        results_list[[result_counter]] <- data.frame(
          n = n,
          p = p,
          graph_type = graph_type,
          contamination_scenario = current_scenario,
          simulation_times = simulation_times,
          iteration = iteration,
          method_name = method_name,
          F1 = eval_result$F1,
          Precision = eval_result$Precision,
          Recall = eval_result$Recall,
          Specificity = eval_result$Specificity,
          MCC = eval_result$MCC,
          Sparsity = eval_result$estimated_sparsity,
          True_Sparsity = eval_result$true_sparsity,
          time = elapsed_time,
          stringsAsFactors = FALSE
        )
      } else {
        # Method failed, record NA
        # Still calculate true_sparsity from true_theta
        true_edges <- sum(abs(true_theta) > 1e-6) - p
        total_possible_edges <- p * (p - 1) / 2
        true_sparsity <- 1 - (true_edges / total_possible_edges)
        
        results_list[[result_counter]] <- data.frame(
          n = n,
          p = p,
          graph_type = graph_type,
          contamination_scenario = current_scenario,
          simulation_times = simulation_times,
          iteration = iteration,
          method_name = method_name,
          F1 = NA,
          Precision = NA,
          Recall = NA,
          Specificity = NA,
          MCC = NA,
          Sparsity = NA,
          True_Sparsity = true_sparsity,
          time = elapsed_time,
          stringsAsFactors = FALSE
        )
      }
      
      result_counter <- result_counter + 1
    }
  }
  
  # After completing clean scenario, identify which methods failed (if skip_failed_methods=TRUE)
  if (skip_failed_methods && current_scenario == "clean" && scenario_idx < length(contamination_scenario)) {
    # Get results for clean scenario
    clean_results <- do.call(rbind, results_list)
    clean_scenario_results <- clean_results[clean_results$contamination_scenario == "clean", ]
    
    # Identify methods that failed (all iterations have NA F1)
    for (method_name in names(methods)) {
      method_results <- clean_scenario_results[clean_scenario_results$method_name == method_name, ]
      if (all(is.na(method_results$F1))) {
        failed_methods <- c(failed_methods, method_name)
      }
    }
    
    if (length(failed_methods) > 0 && verbose) {
      cat(sprintf("\n  WARNING: The following methods failed on clean scenario: %s\n", 
                  paste(failed_methods, collapse=", ")))
      cat("  These methods will be skipped in remaining scenarios.\n")
    }
  }
  
}  # End of scenario loop
  
  # Combine all results
  results_df <- do.call(rbind, results_list)
  
  if (verbose) {
    cat(sprintf("  Progress: %d/%d\n", simulation_times, simulation_times))
    cat("\nSimulation complete!\n\n")
  }
  
  # Compute summary statistics
  if (verbose) {
    cat("Computing summary statistics...\n")
  }
  
  summary_list <- list()
  summary_counter <- 1
  
  # Get all unique combinations of method and contamination scenario
  unique_combos <- unique(results_df[, c("method_name", "contamination_scenario")])
  
  for (i in seq_len(nrow(unique_combos))) {
    method_name <- unique_combos$method_name[i]
    scenario <- unique_combos$contamination_scenario[i]
    
    # Filter results for this combination
    combo_results <- results_df[results_df$method_name == method_name & 
                                  results_df$contamination_scenario == scenario, ]
    
    # Remove NA values for computation
    metrics <- c("F1", "Precision", "Recall", "Specificity", "MCC", "Sparsity", "time")
    
    summary_row <- data.frame(
      method_name = method_name,
      contamination_scenario = scenario,
      stringsAsFactors = FALSE
    )
    
    for (metric in metrics) {
      values <- combo_results[[metric]]
      values_clean <- values[!is.na(values)]
      
      if (length(values_clean) > 0) {
        summary_row[[paste0(metric, "_mean")]] <- mean(values_clean)
        summary_row[[paste0(metric, "_sd")]] <- sd(values_clean)
        summary_row[[paste0(metric, "_median")]] <- median(values_clean)
      } else {
        summary_row[[paste0(metric, "_mean")]] <- NA
        summary_row[[paste0(metric, "_sd")]] <- NA
        summary_row[[paste0(metric, "_median")]] <- NA
      }
    }
    
    # Add success rate
    summary_row$success_rate <- sum(!is.na(combo_results$F1)) / nrow(combo_results)
    
    summary_list[[summary_counter]] <- summary_row
    summary_counter <- summary_counter + 1
  }
  
  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL
  
  # Sort by contamination_scenario first, then F1 mean (descending)
  summary_df <- summary_df[order(summary_df$contamination_scenario, -summary_df$F1_mean), ]
  
  if (verbose) {
    cat("\n========================================\n")
    cat("  Summary Results by Scenario\n")
    cat("========================================\n")
    
    # Display results grouped by scenario
    for (scenario in unique(summary_df$contamination_scenario)) {
      scenario_df <- summary_df[summary_df$contamination_scenario == scenario, ]
      
      cat(sprintf("\n--- %s scenario ---\n", scenario))
      for (i in seq_len(nrow(scenario_df))) {
        row <- scenario_df[i, ]
        cat(sprintf("%d. %s\n", i, row$method_name))
        cat(sprintf("   F1: %.4f (±%.4f)\n", row$F1_mean, row$F1_sd))
        cat(sprintf("   Precision: %.4f (±%.4f)\n", row$Precision_mean, row$Precision_sd))
        cat(sprintf("   Recall: %.4f (±%.4f)\n", row$Recall_mean, row$Recall_sd))
        cat(sprintf("   MCC: %.4f (±%.4f)\n", row$MCC_mean, row$MCC_sd))
        cat(sprintf("   Time: %.4f sec\n", row$time_mean))
        if (row$success_rate < 1) {
          cat(sprintf("   Success rate: %.1f%%\n", row$success_rate * 100))
        }
        cat("\n")
      }
    }
    cat("========================================\n\n")
  }
  
  # Return results
  result <- list(
    results = results_df,
    summary = summary_df,
    true_theta = true_theta,
    graph_info = list(
      adjacency = true_adj,
      theta = true_theta,
      true_sparsity = true_sparsity,
      true_edges = true_edges,
      total_edges = total_possible_edges
    ),
    params = list(
      n = n,
      p = p,
      graph_type = graph_type,
      contamination_scenario = contamination_scenario,
      contamination_rate = contamination_rate,
      simulation_times = simulation_times,
      save_matrices = save_matrices,
      seed = seed
    )
  )
  
  # Add estimated matrices if saved
  if (save_matrices) {
    result$estimated_matrices <- list(
      theta = estimated_theta_list,
      adjacency = estimated_adj_list
    )
    
    if (verbose) {
      cat(sprintf("Saved %d estimated precision matrices and adjacency matrices.\n", 
                  length(estimated_theta_list)))
    }
  }
  
  return(result)
}

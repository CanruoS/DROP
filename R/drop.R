# ============================================================================
# DROP: Distributionally Robust Optimization for Precision Matrices
# C++ accelerated with parallel computing support
# ============================================================================

#' @useDynLib DROP, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats cor cov median qnorm rcauchy sd
#' @importFrom utils flush.console head
#' @keywords internal
"_PACKAGE"

#' Distributionally Robust Optimization for Precision Matrices
#'
#' Estimates the precision matrix (inverse covariance matrix) for Gaussian graphical
#' models using an adaptive robust approach based on coordinate descent optimization
#' with rank-based inverse normal transformation. Supports C++ acceleration and 
#' parallel computing for improved performance.
#'
#' @param X An n x p data matrix where n is the number of observations and p is
#'   the number of variables.
#' @param nlambda Integer. Number of lambda values to consider. Default is 30.
#' @param gamma Numeric. EBIC tuning parameter. Default is 0.5.
#' @param use_npn Logical. If TRUE (default), applies rank-based inverse normal
#'   transformation (nonparanormal) to the data for robustness. If FALSE, uses raw data.
#' @param verbose Logical. If TRUE, prints progress information. Default is FALSE.
#'
#' @return A list containing:
#'   \item{est_K}{Estimated precision matrix (p x p)}
#'   \item{adj}{Adjacency matrix (0/1) indicating edges}
#'   \item{edges}{Logical matrix indicating edges}
#'   \item{best_lambda}{Optimal lambda value selected by EBIC}
#'   \item{best_idx}{Index of the optimal lambda in the lambda sequence}
#'   \item{lambdas}{Vector of lambda values considered}
#'   \item{ebic_values}{EBIC values for each lambda}
#'   \item{edge_counts}{Number of edges for each lambda}
#'   \item{loss_values}{Loss values for each lambda}
#'   \item{computation_time}{Total computation time in seconds}
#'   \item{n_cores}{Number of cores used}
#'
#' @examples
#' \dontrun{
#' # Generate contaminated data using drop.generator
#' set.seed(123)
#' data_gen <- drop.generator(
#'   graph_type = "hub",
#'   n = 500,
#'   p = 10,
#'   contamination_scenario = "leverage",
#'   contamination_rate = 0.10,
#'   seed = 123
#' )
#'
#' # Estimate precision matrix with DROP
#' result <- drop(X = data_gen$X, verbose = TRUE)
#'
#' # Evaluate performance
#' eval_result <- drop.evaluate(
#'   estimated_theta = result$est_K,
#'   true_theta = data_gen$omega
#' )
#' 
#' # Inspect selected performance metrics
#' cat("F1 Score:", eval_result$F1, "\n")          
#' cat("MCC:", eval_result$MCC, "\n")      
#'
#' }
#'
#' @export
drop <- function(X, nlambda = 30, gamma = 0.5, use_npn = TRUE, verbose = FALSE) {
  n <- nrow(X)
  p <- ncol(X)
  
  # Auto-detect cores for parallel processing
  n_cores <- max(1, parallel::detectCores() - 1)
  use_parallel <- (nlambda >= n_cores && n_cores > 1)
  
  if (verbose) {
    cat(sprintf("[DROP] n=%d, p=%d, cores=%d\n", n, p, n_cores))
  }
  
  # Nonparanormal (rank-based) transformation
  if (use_npn) {
    X_use <- apply(X, 2, function(col) qnorm((rank(col) - 0.5) / length(col)))
    if (verbose) cat("  ✓ Rank transformation\n")
  } else {
    X_use <- X
  }
  
  # Initialize K using CovTools (required)
  cov0 <- CovTools::CovEst.2010OAS(X_use)$S
  K0 <- solve(cov0)
  
  # Lambda sequence
  S <- cov(X_use)
  lambda_max <- max(max(S - diag(p)), -min(S - diag(p)))
  lambda_min <- 0.01 * lambda_max
  lambdas <- exp(seq(log(lambda_max), log(lambda_min), length = nlambda))
  
  if (verbose) {
    cat(sprintf("  λ: [%.2e, %.2e], n=%d\n", lambda_min, lambda_max, nlambda))
  }
  
  start_time <- Sys.time()
  
  # Process lambda sequence
  if (use_parallel) {
    if (verbose) cat(sprintf("  🚀 Parallel mode (%d cores)\n", n_cores))
    
    # Split lambda indices into chunks
    chunks <- split(1:nlambda, cut(1:nlambda, breaks = n_cores, labels = FALSE))
    
    # Create cluster
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    # Export data and functions
    parallel::clusterExport(cl, c("X_use", "K0", "lambdas", "gamma"),
                           envir = environment())
    
    # Load DROP package on each worker
    parallel::clusterEvalQ(cl, {
      library(DROP)
    })
    
    # Process in parallel
    chunk_results <- parallel::parLapply(cl, chunks, function(chunk_indices) {
      # Worker function: process lambda range
      K_current <- K0
      results <- list()
      
      for (i in seq_along(chunk_indices)) {
        j <- chunk_indices[i]
        
        # Use C++ functions
        result <- update_k_cpp(X_use, K_current, lambdas[j], 
                              tol = 1e-4, max_iter = 1000)
        
        K_j <- result$K_update
        diag(K_j) <- diag(K0)
        
        # Calculate EBIC with C++
        ebic_result <- calculate_ebic_cpp(K_j, X_use, gamma)
        
        results[[i]] <- list(
          idx = j,
          K = K_j,
          ebic = ebic_result$ebic,
          loss = ebic_result$loss,
          num_edges = ebic_result$num_edges
        )
        
        K_current <- K_j  # Warm start
      }
      
      results
    })
    
    # Combine results
    all_results <- do.call(c, chunk_results)
    all_results <- all_results[order(sapply(all_results, function(x) x$idx))]
    
    K_list <- lapply(all_results, function(x) x$K)
    ebic_values <- sapply(all_results, function(x) x$ebic)
    loss_values <- sapply(all_results, function(x) x$loss)
    edge_counts <- sapply(all_results, function(x) x$num_edges)
    
  } else {
    # Serial processing
    if (verbose) cat("  → Serial mode\n")
    
    K_current <- K0
    all_results <- list()
    
    for (j in 1:nlambda) {
      if (verbose && (j == 1 || j %% 10 == 0)) {
        elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        cat(sprintf("    λ %d/%d (%.1fs)\n", j, nlambda, elapsed))
      }
      
      # Use C++ functions
      result <- update_k_cpp(X_use, K_current, lambdas[j], 
                            tol = 1e-4, max_iter = 1000)
      
      K_j <- result$K_update
      diag(K_j) <- diag(K0)
      
      # Calculate EBIC with C++
      ebic_result <- calculate_ebic_cpp(K_j, X_use, gamma)
      
      all_results[[j]] <- list(
        idx = j,
        K = K_j,
        ebic = ebic_result$ebic,
        loss = ebic_result$loss,
        num_edges = ebic_result$num_edges
      )
      
      K_current <- K_j  # Warm start
    }
    
    K_list <- lapply(all_results, function(x) x$K)
    ebic_values <- sapply(all_results, function(x) x$ebic)
    loss_values <- sapply(all_results, function(x) x$loss)
    edge_counts <- sapply(all_results, function(x) x$num_edges)
  }
  
  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Select optimal
  best_idx <- which.min(ebic_values)
  best_lambda <- lambdas[best_idx]
  K_est <- K_list[[best_idx]]
  
  edges <- (abs(K_est) > 1e-6)
  diag(edges) <- FALSE
  
  if (verbose) {
    cat(sprintf("\n  ✓ Done: %.2f sec (%.3f sec/λ)\n", total_time, total_time/nlambda))
    cat(sprintf("  Best: λ=%.4e (#%d), %d edges\n",
                best_lambda, best_idx, edge_counts[best_idx]))
  }
  
  list(
    est_K = K_est,
    adj = edges * 1,
    edges = edges,
    best_lambda = best_lambda,
    best_idx = best_idx,
    lambdas = lambdas,
    ebic_values = ebic_values,
    edge_counts = edge_counts,
    loss_values = loss_values,
    computation_time = total_time,
    n_cores = n_cores
  )
}

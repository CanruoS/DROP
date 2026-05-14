# ============================================================================
# Generate Plot Data for GGM Simulations
# ============================================================================

#' Generate Data for Gaussian Graphical Model Simulations
#'
#' Generates graph structures and contaminated data for Gaussian graphical model
#' simulations. Supports various graph types and contamination scenarios with
#' configurable graph and contamination settings.
#'
#' @param graph_type Character. Type of graph structure. One of: "band", "hub", "cluster", "random", "scale-free".
#' @param n Integer. Number of observations (samples). Required.
#' @param p Integer. Number of variables (dimensions). Required.
#' @param contamination_scenario Character. Type of contamination. One of: "clean" (No contamination), "cauchy" (Heavy-tailed Cauchy distribution contamination), "leverage" (Leverage point contamination with large variance).
#' @param contamination_rate Numeric. Proportion of contaminated observations. Default is 0.1.
#' @param g Integer. Graph parameter passed to huge::huge.generator(). For "band": bandwidth. For "hub": number of hubs. For "cluster": number of clusters. Default is NULL (uses: g=1 for band, g=floor(p/20) for hub, g=floor(p/10) for cluster).
#' @param prob Numeric. Connection probability passed to huge::huge.generator(). For "cluster": within-cluster connection probability. For "random": edge connection probability. Default is NULL (uses: prob=0.2 for cluster, prob=0.02 for random).
#' @param cauchy_location Numeric. Location parameter for Cauchy contamination. Default is 0. Only used when contamination_scenario = "cauchy".
#' @param cauchy_scale Numeric. Scale parameter for Cauchy contamination. Default is 5. Only used when contamination_scenario = "cauchy".
#' @param leverage_multiplier Numeric. Variance multiplier for leverage contamination. The contaminated covariance is Sigma * leverage_multiplier. Default is 100. Only used when contamination_scenario = "leverage".
#' @param seed Integer. Random seed for reproducibility. If NULL, uses graph-type-specific default seeds. Default is NULL.
#' @param return_graph_info Logical. If TRUE, returns additional graph structure information (theta, omega, sparsity). Default is TRUE.
#' @param verbose Logical. If TRUE, prints progress information. Default is FALSE.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{X}: Generated data matrix (n x p)
#'   \item \code{graph_type}: Graph type used
#'   \item \code{contamination_scenario}: Contamination scenario used
#'   \item \code{n}: Number of observations
#'   \item \code{p}: Number of variables
#'   \item \code{contamination_rate}: Contamination rate used
#'   \item \code{contaminated_indices}: Vector of contaminated observation indices (NULL if clean)
#'   \item \code{seed}: Seed used for generation
#'   \item \code{theta}: True adjacency matrix (p x p), if return_graph_info = TRUE
#'   \item \code{omega}: True precision matrix (p x p), if return_graph_info = TRUE
#'   \item \code{sigma}: True covariance matrix (p x p), if return_graph_info = TRUE
#'   \item \code{true_edges}: Number of true edges, if return_graph_info = TRUE
#'   \item \code{total_edges}: Total possible edges, if return_graph_info = TRUE
#'   \item \code{sparsity}: True sparsity level (proportion of zero edges = 1 - edges/total), if return_graph_info = TRUE
#' }
#'
#' @examples
#' \dontrun{
#' # Generate clean data with cluster graph
#' data1 <- drop.generator(
#'   graph_type = "cluster",
#'   n = 200,
#'   p = 50,
#'   contamination_scenario = "clean"
#' )
#'
#' # Generate contaminated data with hub graph
#' data2 <- drop.generator(
#'   graph_type = "hub",
#'   n = 200,
#'   p = 50,
#'   contamination_scenario = "cauchy",
#'   contamination_rate = 0.15,
#'   seed = 123
#' )
#'
#' # Access the data and graph structure
#' X <- data2$X
#' true_precision <- data2$omega
#' adjacency <- data2$theta
#' 
#' # Generate data with custom graph parameters
#' data3 <- drop.generator(
#'   graph_type = "hub",
#'   n = 200,
#'   p = 50,
#'   g = 10,  # Use 10 hubs instead of default
#'   contamination_scenario = "clean"
#' )
#' 
#' # Generate random graph with custom connection probability
#' data4 <- drop.generator(
#'   graph_type = "random",
#'   n = 150,
#'   p = 30,
#'   prob = 0.05,  # 5% connection probability
#'   contamination_scenario = "clean"
#' )
#' }
#'
#' @details
#' This function generates synthetic data for Gaussian graphical model simulations.
#' 
#' Default graph parameters when g or prob are NULL:
#' \itemize{
#'   \item Band graph: bandwidth g=1
#'   \item Hub graph: g=floor(p/20) hubs
#'   \item Cluster graph: g=floor(p/10) clusters, connection probability prob=0.2
#'   \item Random graph: connection probability prob=0.02
#'   \item Scale-free graph: uses huge package defaults
#' }
#'
#' When contamination_scenario is not "clean":
#' \itemize{
#'   \item "cauchy": Adds Cauchy(cauchy_location, cauchy_scale) noise to contaminated observations
#'   \item "leverage": Replaces contaminated observations with N(0, leverage_multiplier * Sigma)
#' }
#'
#' @export
drop.generator <- function(graph_type = "cluster",
                               n,
                               p,
                               contamination_scenario = "clean",
                               contamination_rate = 0.1,
                               g = NULL,
                               prob = NULL,
                               cauchy_location = 0,
                               cauchy_scale = 5,
                               leverage_multiplier = 100,
                               seed = NULL,
                               return_graph_info = TRUE,
                               verbose = FALSE) {
  
  # Validate inputs
  valid_graphs <- c("band", "hub", "cluster", "random", "scale-free")
  if (!(graph_type %in% valid_graphs)) {
    stop(sprintf("Invalid graph_type. Must be one of: %s",
                 paste(valid_graphs, collapse = ", ")))
  }
  
  valid_scenarios <- c("clean", "cauchy", "leverage")
  if (!(contamination_scenario %in% valid_scenarios)) {
    stop(sprintf("Invalid contamination_scenario. Must be one of: %s",
                 paste(valid_scenarios, collapse = ", ")))
  }
  
  if (contamination_rate < 0 || contamination_rate > 1) {
    stop("contamination_rate must be between 0 and 1")
  }
  
  if (cauchy_scale <= 0) {
    stop("cauchy_scale must be positive")
  }
  
  if (leverage_multiplier <= 0) {
    stop("leverage_multiplier must be positive")
  }
  
  # Set default graph parameters if not specified
  if (is.null(g)) {
    if (graph_type == "band") {
      g <- 1
    } else if (graph_type == "hub") {
      g <- max(1, floor(p / 20))
    } else if (graph_type == "cluster") {
      g <- max(2, floor(p / 10))
    }
  }
  
  if (is.null(prob)) {
    if (graph_type == "cluster") {
      prob <- 0.2
    } else if (graph_type == "random") {
      prob <- 0.02
    }
  }
  
  # Use graph-type-specific default seeds if not provided
  if (is.null(seed)) {
    optimal_seeds <- list(
      cluster = 123,
      band = 123,
      hub = 123,
      random = 456,
      "scale-free" = 456
    )
    seed <- optimal_seeds[[graph_type]]
  }
  
  if (verbose) {
    cat("\n========================================\n")
    cat("  Generate Plot Data\n")
    cat("========================================\n")
    cat(sprintf("Graph type: %s\n", graph_type))
    if (!is.null(g)) {
      cat(sprintf("  g parameter: %d\n", g))
    }
    if (!is.null(prob)) {
      cat(sprintf("  prob parameter: %.4f\n", prob))
    }
    cat(sprintf("Dimensions: n=%d, p=%d\n", n, p))
    cat(sprintf("Contamination: %s (rate=%.2f)\n", 
                contamination_scenario, contamination_rate))
    if (contamination_scenario == "cauchy") {
      cat(sprintf("  Cauchy params: location=%.2f, scale=%.2f\n", 
                  cauchy_location, cauchy_scale))
    } else if (contamination_scenario == "leverage") {
      cat(sprintf("  Leverage multiplier: %.1f\n", leverage_multiplier))
    }
    cat(sprintf("Seed: %d\n", seed))
    cat("========================================\n\n")
  }
  
  set.seed(seed)
  
  if (!requireNamespace("huge", quietly = TRUE)) {
    stop("Package 'huge' is required but not installed. Please install it.")
  }
  
  if (verbose) cat("Generating graph structure...\n")
  
  base_args <- list(n = n, d = p, graph = graph_type, verbose = FALSE)
  
  if (graph_type == "band") {
    if (!is.null(g)) base_args$g <- g
    graph_data <- do.call(huge::huge.generator, base_args)
    
  } else if (graph_type == "hub") {
    if (!is.null(g)) base_args$g <- g
    graph_data <- do.call(huge::huge.generator, base_args)
    
  } else if (graph_type == "cluster") {
    if (!is.null(g)) base_args$g <- g
    if (!is.null(prob)) base_args$prob <- prob
    graph_data <- do.call(huge::huge.generator, base_args)
    
  } else if (graph_type == "random") {
    if (!is.null(prob)) base_args$prob <- prob
    graph_data <- do.call(huge::huge.generator, base_args)
    
  } else if (graph_type == "scale-free") {
    graph_data <- do.call(huge::huge.generator, base_args)
  }
  
  Omega <- graph_data$omega
  Theta <- graph_data$theta
  Sigma <- graph_data$sigma
  
  true_edges <- sum(Theta[upper.tri(Theta)] != 0)
  total_edges <- p * (p - 1) / 2
  sparsity <- 1 - (true_edges / total_edges)
  
  if (verbose) {
    cat(sprintf("  Edges: %d / %d\n", true_edges, total_edges))
    cat(sprintf("  Sparsity: %.4f\n", sparsity))
  }
  
  if (verbose) cat("Generating samples...\n")
  
  set.seed(seed)
  X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  
  contaminated_idx <- NULL
  if (contamination_scenario != "clean") {
    if (verbose) cat(sprintf("Adding %s contamination...\n", contamination_scenario))
    
    n_contaminated <- floor(n * contamination_rate)
    if (n_contaminated > 0) {
      contaminated_idx <- sample(1:n, n_contaminated)
      
      if (contamination_scenario == "cauchy") {
        # Heavy-tailed contamination (additive)
        for (i in contaminated_idx) {
          X[i, ] <- X[i, ] + rcauchy(p, location = cauchy_location, scale = cauchy_scale)
        }
        
      } else if (contamination_scenario == "leverage") {
        # Leverage points with large variance (replacement)
        leverage_cov <- Sigma * leverage_multiplier
        X[contaminated_idx, ] <- MASS::mvrnorm(
          n_contaminated, 
          mu = rep(0, p), 
          Sigma = leverage_cov
        )
      }
      
      if (verbose) {
        cat(sprintf("  Contaminated %d / %d observations\n", 
                    n_contaminated, n))
        cat("  Contaminated indices:", paste(head(contaminated_idx, 10), collapse=", "))
        if (n_contaminated > 10) {
          cat(sprintf(", ... (%d more)\n", n_contaminated - 10))
        } else {
          cat("\n")
        }
      }
    }
  }
  
  if (verbose) cat("✓ Data generation complete\n\n")
  
  result <- list(
    X = X,
    graph_type = graph_type,
    contamination_scenario = contamination_scenario,
    n = n,
    p = p,
    contamination_rate = contamination_rate,
    contaminated_indices = contaminated_idx,
    seed = seed
  )
  
  if (return_graph_info) {
    result$theta <- Theta
    result$omega <- Omega
    result$sigma <- Sigma
    result$true_edges <- true_edges
    result$total_edges <- total_edges
    result$sparsity <- sparsity
  }
  
  return(result)
}

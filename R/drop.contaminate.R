#' Add Contamination to Existing Data
#'
#' Add various types of contamination (outliers) to an existing data matrix. 
#' Supports Cauchy contamination (additive outliers) and leverage contamination 
#' (high leverage outliers).
#'
#' @param data A numeric matrix (n x p) representing the data to be contaminated. Each row is an observation and each column is a variable.
#' @param contamination_scenario Character string specifying the contamination type. Options: "cauchy" (additive outliers) or "leverage" (high leverage outliers). Default is "cauchy".
#' @param contamination_rate Numeric. Proportion of contaminated observations. Default is 0.1.
#' @param cauchy_location Numeric value for the location parameter of Cauchy distribution. Used for Cauchy contamination. Default is 0.
#' @param cauchy_scale Numeric value for the scale parameter of Cauchy distribution. Used for Cauchy contamination. Default is 5.
#' @param leverage_multiplier Numeric value specifying the factor to multiply the covariance matrix by when generating leverage points. Default is 100.
#' @param seed Optional integer to set random seed for reproducibility. Default is NULL.
#' @param verbose Logical indicating whether to print progress information. Default is FALSE.
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{contaminated_data}: The data matrix with contamination added
#'   \item \code{contaminated_indices}: Indices of the contaminated samples
#'   \item \code{contamination_scenario}: The type of contamination applied
#'   \item \code{contamination_rate}: The contamination rate used
#'   \item \code{n_contaminated}: Number of contaminated samples
#'   \item \code{seed}: The random seed used (if specified)
#' }
#'
#' @details
#' The function supports two types of contamination:
#' \itemize{
#'   \item \strong{Cauchy}: Additive contamination. Adds Cauchy-distributed noise
#'         to the contaminated samples: X[i,] <- X[i,] + rcauchy(p, location, scale)
#'   \item \strong{Leverage}: Replacement contamination. Replaces contaminated samples
#'         with samples from a multivariate normal distribution with an inflated
#'         covariance matrix (multiplier * Sigma).
#' }
#'
#' @examples
#' \dontrun{
#' # Generate clean data using drop.generator
#' set.seed(123)
#' data_gen <- drop.generator(
#'   graph_type = "hub",
#'   n = 100,
#'   p = 10,
#'   contamination_scenario = "clean",
#'   seed = 123
#' )
#' clean_data <- data_gen$X
#'
#' # Add Cauchy contamination
#' result_cauchy <- drop.contaminate(
#'   data = clean_data,
#'   contamination_scenario = "cauchy",
#'   contamination_rate = 0.1,
#'   seed = 123
#' )
#'
#' # Add leverage contamination
#' result_leverage <- drop.contaminate(
#'   data = clean_data,
#'   contamination_scenario = "leverage",
#'   contamination_rate = 0.15,
#'   leverage_multiplier = 100,
#'   seed = 456
#' )
#'
#' # Access results
#' head(result_cauchy$contaminated_data)
#' result_cauchy$contaminated_indices
#' }
#'
#' @export
drop.contaminate <- function(data,
                              contamination_scenario = "cauchy",
                              contamination_rate = 0.1,
                              cauchy_location = 0,
                              cauchy_scale = 5,
                              leverage_multiplier = 100,
                              seed = NULL,
                              verbose = FALSE) {
  
  # ============================================================================
  # Input validation
  # ============================================================================
  
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("'data' must be a matrix or data frame")
  }
  
  # Convert data frame to matrix if needed
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }
  
  if (!is.numeric(data)) {
    stop("'data' must contain numeric values")
  }
  
  if (any(is.na(data))) {
    stop("'data' contains NA values")
  }
  
  n <- nrow(data)
  p <- ncol(data)
  
  if (n < 2 || p < 1) {
    stop("'data' must have at least 2 rows and 1 column")
  }
  
  # Validate contamination scenario
  valid_scenarios <- c("cauchy", "leverage")
  if (!(contamination_scenario %in% valid_scenarios)) {
    stop(sprintf("'contamination_scenario' must be one of: %s", 
                 paste(valid_scenarios, collapse = ", ")))
  }
  
  # Validate contamination rate
  if (!is.numeric(contamination_rate) || length(contamination_rate) != 1) {
    stop("'contamination_rate' must be a single numeric value")
  }
  if (contamination_rate < 0 || contamination_rate > 1) {
    stop("'contamination_rate' must be between 0 and 1")
  }
  
  # Validate other parameters
  if (!is.numeric(cauchy_location) || length(cauchy_location) != 1) {
    stop("'cauchy_location' must be a single numeric value")
  }
  if (!is.numeric(cauchy_scale) || length(cauchy_scale) != 1 || cauchy_scale <= 0) {
    stop("'cauchy_scale' must be a single positive numeric value")
  }
  if (!is.numeric(leverage_multiplier) || length(leverage_multiplier) != 1 || leverage_multiplier <= 0) {
    stop("'leverage_multiplier' must be a single positive numeric value")
  }
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("'verbose' must be a single logical value (TRUE or FALSE)")
  }
  
  # ============================================================================
  # Set random seed if provided
  # ============================================================================
  
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1) {
      stop("'seed' must be a single numeric value or NULL")
    }
    set.seed(seed)
    if (verbose) {
      cat(sprintf("Random seed set to: %d\n", seed))
    }
  }
  
  # ============================================================================
  # Determine contaminated samples
  # ============================================================================
  
  # Calculate number of samples to contaminate
  n_contaminated <- max(1, floor(n * contamination_rate))
  
  if (verbose) {
    cat(sprintf("\n=== Data Contamination ===\n"))
    cat(sprintf("Original data dimensions: n=%d, p=%d\n", n, p))
    cat(sprintf("Contamination scenario: %s\n", contamination_scenario))
    cat(sprintf("Contamination rate: %.1f%% (%d samples)\n", 
                contamination_rate * 100, n_contaminated))
  }
  
  # Randomly select samples to contaminate
  contaminated_indices <- sample(1:n, size = n_contaminated, replace = FALSE)
  contaminated_indices <- sort(contaminated_indices)
  
  # Make a copy of the data
  contaminated_data <- data
  
  # ============================================================================
  # Apply contamination based on scenario
  # ============================================================================
  
  if (contamination_scenario == "cauchy") {
    # Cauchy contamination: additive
    contaminated_data <- add_cauchy_contamination_internal(
      data = contaminated_data,
      indices = contaminated_indices,
      cauchy_location = cauchy_location,
      cauchy_scale = cauchy_scale,
      verbose = verbose
    )
    
  } else if (contamination_scenario == "leverage") {
    # Leverage contamination: replacement
    contaminated_data <- add_leverage_contamination_internal(
      data = contaminated_data,
      indices = contaminated_indices,
      leverage_multiplier = leverage_multiplier,
      verbose = verbose
    )
  }
  
  if (verbose) {
    cat(sprintf("=========================\n\n"))
  }
  
  # ============================================================================
  # Return results
  # ============================================================================
  
  result <- list(
    contaminated_data = contaminated_data,
    contaminated_indices = contaminated_indices,
    contamination_scenario = contamination_scenario,
    contamination_rate = contamination_rate,
    n_contaminated = n_contaminated
  )
  
  # Add seed to output if it was used
  if (!is.null(seed)) {
    result$seed <- seed
  }
  
  return(result)
}


# ==============================================================================
# Internal helper function: Add Cauchy contamination
# ==============================================================================
add_cauchy_contamination_internal <- function(data, indices, 
                                              cauchy_location = 0, 
                                              cauchy_scale = 5, 
                                              verbose = TRUE) {
  p <- ncol(data)
  
  for (i in indices) {
    # Add Cauchy noise to each contaminated sample
    cauchy_noise <- rcauchy(p, location = cauchy_location, scale = cauchy_scale)
    data[i, ] <- data[i, ] + cauchy_noise
  }
  
  if (verbose) {
    cat(sprintf("Added Cauchy contamination (location=%.1f, scale=%.1f) to %d samples\n", 
                cauchy_location, cauchy_scale, length(indices)))
  }
  
  return(data)
}


# ==============================================================================
# Internal helper function: Add leverage contamination
# ==============================================================================
add_leverage_contamination_internal <- function(data, indices, 
                                               leverage_multiplier = 100, 
                                               verbose = TRUE) {
  p <- ncol(data)
  n_contaminated <- length(indices)
  
  # Estimate covariance matrix from original data
  Sigma <- stats::cov(data)
  
  # Create inflated covariance matrix
  Sigma_leverage <- Sigma * leverage_multiplier
  
  # Generate leverage points (sample from inflated covariance)
  tryCatch({
    leverage_samples <- MASS::mvrnorm(n = n_contaminated, 
                                      mu = rep(0, p), 
                                      Sigma = Sigma_leverage)
    
    # Replace original data with leverage samples
    data[indices, ] <- leverage_samples
    
    if (verbose) {
      cat(sprintf("Added leverage contamination (multiplier=%d) to %d samples\n", 
                  leverage_multiplier, n_contaminated))
    }
    
  }, error = function(e) {
    warning(sprintf("Leverage contamination failed: %s. Using fallback method.", e$message))
    # Fallback: simply scale the original samples
    data[indices, ] <- data[indices, ] * sqrt(leverage_multiplier)
  })
  
  return(data)
}

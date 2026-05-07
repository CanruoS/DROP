#' Evaluate Precision Matrix Estimation Accuracy
#'
#' Evaluate the accuracy of an estimated precision matrix by comparing it to the true
#' precision matrix. Computes edge recovery metrics including F1 score, precision, recall,
#' specificity, MCC, and sparsity.
#'
#' @param estimated_theta A numeric matrix representing the estimated precision matrix.
#' @param true_theta A numeric matrix representing the true precision matrix.
#' @param threshold Numeric threshold for determining non-zero edges. Edges with absolute
#'   values greater than this threshold are considered non-zero. Default is 1e-6.
#'
#' @return A list with the following components:
#'   \item{F1}{F1 score (harmonic mean of precision and recall)}
#'   \item{Precision}{Precision (TP / (TP + FP))}
#'   \item{Recall}{Recall/Sensitivity (TP / (TP + FN))}
#'   \item{Specificity}{Specificity (TN / (TN + FP))}
#'   \item{MCC}{Matthews Correlation Coefficient}
#'   \item{TP}{Number of true positives}
#'   \item{FP}{Number of false positives}
#'   \item{FN}{Number of false negatives}
#'   \item{TN}{Number of true negatives}
#'   \item{estimated_edges}{Number of edges in estimated graph}
#'   \item{true_edges}{Number of edges in true graph}
#'   \item{total_possible_edges}{Total number of possible edges (p*(p-1)/2)}
#'   \item{estimated_sparsity}{Sparsity of estimated graph (1 - estimated_edges/total_possible_edges)}
#'   \item{true_sparsity}{Sparsity of true graph (1 - true_edges/total_possible_edges)}
#'
#' @details
#' The function compares the edge structures of two precision matrices by:
#' \enumerate{
#'   \item Converting each matrix to a binary adjacency matrix using the threshold
#'   \item Computing confusion matrix statistics (TP, FP, FN, TN) on the upper triangle
#'   \item Calculating standard classification metrics
#'   \item Computing sparsity as the proportion of zero edges
#' }
#'
#' Edge detection uses the upper triangle of the matrices (excluding diagonal) to avoid
#' counting edges twice in undirected graphs.
#'
#' @examples
#' \dontrun{
#' # Example 1: Evaluate DROP on contaminated data
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
#' # Estimate with DROP
#' drop_result <- drop(X = data_gen$X, verbose = FALSE)
#'
#' # Evaluate performance
#' eval_result <- drop.evaluate(
#'   estimated_theta = drop_result$est_K,
#'   true_theta = data_gen$omega
#' )
#'
#' cat("DROP Performance on Contaminated Data:\n")
#' cat("  F1 Score:", eval_result$F1, "\n")          
#' cat("  MCC:", eval_result$MCC, "\n")             
#'
#' # Example 2: Compare with other methods
#' # HUGE Glasso (non-robust method)
#' huge_fit <- huge::huge(data_gen$X, method = "glasso", verbose = FALSE)
#' huge_select <- huge::huge.select(huge_fit, criterion = "ebic", verbose = FALSE)
#' 
#' huge_eval <- drop.evaluate(
#'   estimated_theta = huge_select$opt.icov,
#'   true_theta = data_gen$omega
#' )
#'
#' cat("\nHUGE Glasso Performance:\n")
#' cat("  F1 Score:", huge_eval$F1, "\n")
#' cat("  MCC:", huge_eval$MCC, "\n")
#'
#' }
#'
#' @export
drop.evaluate <- function(estimated_theta,
                          true_theta,
                          threshold = 1e-6) {
  
  # ============================================================================
  # Input validation
  # ============================================================================
  
  if (!is.matrix(estimated_theta) && !is.data.frame(estimated_theta)) {
    stop("'estimated_theta' must be a matrix or data frame")
  }
  
  if (!is.matrix(true_theta) && !is.data.frame(true_theta)) {
    stop("'true_theta' must be a matrix or data frame")
  }
  
  # Convert to matrix if needed
  if (is.data.frame(estimated_theta)) {
    estimated_theta <- as.matrix(estimated_theta)
  }
  if (is.data.frame(true_theta)) {
    true_theta <- as.matrix(true_theta)
  }
  
  if (!is.numeric(estimated_theta)) {
    stop("'estimated_theta' must contain numeric values")
  }
  if (!is.numeric(true_theta)) {
    stop("'true_theta' must contain numeric values")
  }
  
  # Check dimensions match
  if (!all(dim(estimated_theta) == dim(true_theta))) {
    stop(sprintf("Dimension mismatch: estimated_theta is %dx%d but true_theta is %dx%d",
                 nrow(estimated_theta), ncol(estimated_theta),
                 nrow(true_theta), ncol(true_theta)))
  }
  
  # Check square matrices
  if (nrow(estimated_theta) != ncol(estimated_theta)) {
    stop("'estimated_theta' must be a square matrix")
  }
  if (nrow(true_theta) != ncol(true_theta)) {
    stop("'true_theta' must be a square matrix")
  }
  
  p <- nrow(estimated_theta)
  
  if (p < 2) {
    stop("Matrices must be at least 2x2")
  }
  
  # Validate threshold
  if (!is.numeric(threshold) || length(threshold) != 1 || threshold < 0) {
    stop("'threshold' must be a single non-negative numeric value")
  }
  
  # Handle NULL or all-zero estimated matrix
  if (is.null(estimated_theta) || all(abs(estimated_theta) < threshold)) {
    warning("Estimated matrix is NULL or all zeros. Using identity matrix.")
    estimated_theta <- diag(p)
  }
  
  # ============================================================================
  # Extract upper triangle for edge comparison (exclude diagonal)
  # ============================================================================
  
  upper_est <- estimated_theta[upper.tri(estimated_theta)]
  upper_true <- true_theta[upper.tri(true_theta)]
  
  # Convert to binary edge indicators
  est_edges <- (abs(upper_est) > threshold)
  true_edges <- (abs(upper_true) > threshold)
  
  # ============================================================================
  # Compute confusion matrix
  # ============================================================================
  
  TP <- sum(est_edges & true_edges)
  FP <- sum(est_edges & !true_edges)
  FN <- sum(!est_edges & true_edges)
  TN <- sum(!est_edges & !true_edges)
  
  # ============================================================================
  # Compute classification metrics
  # ============================================================================
  
  # Precision: TP / (TP + FP)
  precision <- ifelse(TP + FP > 0, TP / (TP + FP), 0)
  
  # Recall/Sensitivity: TP / (TP + FN)
  recall <- ifelse(TP + FN > 0, TP / (TP + FN), 0)
  
  # F1 score: harmonic mean of precision and recall
  f1 <- ifelse(precision + recall > 0, 
               2 * precision * recall / (precision + recall), 
               0)
  
  # Specificity: TN / (TN + FP)
  specificity <- ifelse(TN + FP > 0, TN / (TN + FP), 0)
  
  # Matthews Correlation Coefficient (MCC)
  # MCC = (TP*TN - FP*FN) / sqrt((TP+FP)(TP+FN)(TN+FP)(TN+FN))
  mcc_num <- as.numeric(TP) * as.numeric(TN) - as.numeric(FP) * as.numeric(FN)
  mcc_den <- sqrt(as.numeric(TP + FP) * as.numeric(TP + FN) * 
                  as.numeric(TN + FP) * as.numeric(TN + FN))
  mcc <- ifelse(mcc_den > 0, mcc_num / mcc_den, 0)
  
  # ============================================================================
  # Compute sparsity metrics
  # ============================================================================
  
  # Total possible edges (upper triangle, excluding diagonal)
  total_possible_edges <- p * (p - 1) / 2
  
  # Number of edges in each graph
  estimated_edges <- sum(est_edges)
  true_edges_count <- sum(true_edges)
  
  # Sparsity = proportion of zero edges = 1 - (number of edges / total possible edges)
  estimated_sparsity <- 1 - (estimated_edges / total_possible_edges)
  true_sparsity <- 1 - (true_edges_count / total_possible_edges)
  
  # ============================================================================
  # Return results
  # ============================================================================
  
  result <- list(
    F1 = f1,
    Precision = precision,
    Recall = recall,
    Specificity = specificity,
    MCC = mcc,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    estimated_edges = estimated_edges,
    true_edges = true_edges_count,
    total_possible_edges = total_possible_edges,
    estimated_sparsity = estimated_sparsity,
    true_sparsity = true_sparsity
  )
  
  return(result)
}

// ============================================================================
// DRO: C++ core functions for runtime optimization
// Compiled with Rcpp
// ============================================================================

#include <Rcpp.h>
using namespace Rcpp;

// Soft-thresholding operator (inline for speed)
inline double soft_threshold(double x, double lambda) {
  if (x > lambda) return x - lambda;
  else if (x < -lambda) return x + lambda;
  else return 0.0;
}

// [[Rcpp::export]]
NumericVector node_mse_vec_cpp(NumericMatrix X, NumericMatrix K, double eps = 1e-8) {
  int n = X.nrow();
  int p = X.ncol();
  NumericVector mse(p);
  
  for (int i = 0; i < p; i++) {
    double Kii = std::max(K(i, i), eps);
    
    // Compute beta_i = -K[-i, i] / Kii
    NumericVector beta_i(p - 1);
    int idx = 0;
    for (int j = 0; j < p; j++) {
      if (j != i) {
        beta_i[idx++] = -K(j, i) / Kii;
      }
    }
    
    // Compute residuals: X[, i] - X[, -i] %*% beta_i
    NumericVector resid(n);
    for (int row = 0; row < n; row++) {
      double pred = 0.0;
      idx = 0;
      for (int col = 0; col < p; col++) {
        if (col != i) {
          pred += X(row, col) * beta_i[idx++];
        }
      }
      resid[row] = X(row, i) - pred;
    }
    
    // Compute mean squared error
    double sum_sq = 0.0;
    for (int row = 0; row < n; row++) {
      sum_sq += resid[row] * resid[row];
    }
    mse[i] = sum_sq / n;
    
    // Apply eps threshold
    if (mse[i] < eps) mse[i] = eps;
  }
  
  return mse;
}

// [[Rcpp::export]]
List update_k_cpp(NumericMatrix X, NumericMatrix K, double lambda, 
                  double tol = 1e-4, int max_iter = 1000) {
  int n = X.nrow();
  int p = X.ncol();
  
  // Clone K to avoid modifying input
  NumericMatrix K_update = clone(K);
  NumericMatrix K_prev(p, p);
  
  bool converged = false;
  int iteration = 0;
  double eps = 1e-8;
  
  // Pre-compute X'X (crossprod)
  NumericMatrix XtX(p, p);
  for (int i = 0; i < p; i++) {
    for (int j = 0; j < p; j++) {
      double sum = 0.0;
      for (int k = 0; k < n; k++) {
        sum += X(k, i) * X(k, j);
      }
      XtX(i, j) = sum;
    }
  }
  
  while (!converged && iteration < max_iter) {
    // Save previous K
    K_prev = clone(K_update);
    
    // Compute current MSE
    NumericVector mse_vec = node_mse_vec_cpp(X, K_update, eps);
    
    // Update each off-diagonal element
    for (int i = 0; i < p - 1; i++) {
      for (int j = i + 1; j < p; j++) {
        
        double Kii = std::max(K_update(i, i), eps);
        double Kjj = std::max(K_update(j, j), eps);
        
        // Gradient term (from data fitting part)
        double kii_kjj = 1.0/Kii + 1.0/Kjj;
        double term1 = 2.0 * XtX(i, j) * kii_kjj;
        
        // term2: 2/Kii^2 * sum(K_update[i, -c(i, j)] * XtX[j, -c(i, j)])
        double term2_sum = 0.0;
        for (int k = 0; k < p; k++) {
          if (k != i && k != j) {
            term2_sum += K_update(i, k) * XtX(j, k);
          }
        }
        double term2 = 2.0 / (Kii * Kii) * term2_sum;
        
        // term3: 2/Kjj^2 * sum(K_update[j, -c(i, j)] * XtX[i, -c(i, j)])
        double term3_sum = 0.0;
        for (int k = 0; k < p; k++) {
          if (k != i && k != j) {
            term3_sum += K_update(j, k) * XtX(i, k);
          }
        }
        double term3 = 2.0 / (Kjj * Kjj) * term3_sum;
        
        double S0 = term1 + term2 + term3;
        double S1 = 2.0 / (Kii * Kii) * XtX(j, j) + 2.0 / (Kjj * Kjj) * XtX(i, i);
        
        // Adaptive weight wᵢⱼ = √MSE⁽ⁱ⁾/Kᵢᵢ + √MSE⁽ʲ⁾/Kⱼⱼ
        double weight_ij = std::sqrt(mse_vec[i]) / Kii + std::sqrt(mse_vec[j]) / Kjj;
        
        // Soft-thresholding
        double threshold = n * lambda * weight_ij;
        double k_ij = soft_threshold(-S0 / S1, threshold / S1);
        
        // Check for NA or Inf
        if (std::isnan(k_ij) || std::isinf(k_ij)) {
          k_ij = 0.0;
        }
        
        K_update(i, j) = k_ij;
        K_update(j, i) = k_ij;
      }
    }
    
    iteration++;
    
    // Check convergence
    double diff_sum = 0.0;
    for (int i = 0; i < p; i++) {
      for (int j = 0; j < p; j++) {
        double diff = std::abs(K_update(i, j) - K_prev(i, j));
        if (!std::isnan(diff)) {
          diff_sum += diff;
        } else {
          diff_sum = R_PosInf;
          break;
        }
      }
      if (std::isinf(diff_sum)) break;
    }
    
    converged = (diff_sum < tol);
  }
  
  return List::create(
    Named("K_update") = K_update,
    Named("iteration") = iteration
  );
}

// [[Rcpp::export]]
List calculate_ebic_cpp(NumericMatrix K, NumericMatrix X, double gamma = 0.5) {
  int n = X.nrow();
  int p = X.ncol();
  
  // Compute loss
  NumericVector mse_vec = node_mse_vec_cpp(X, K);
  double loss = 0.0;
  for (int i = 0; i < p; i++) {
    loss += mse_vec[i];
  }
  
  // Count number of edges
  double threshold = 1e-6;
  int E = 0;
  for (int i = 0; i < p - 1; i++) {
    for (int j = i + 1; j < p; j++) {
      if (std::abs(K(i, j)) > threshold) {
        E++;
      }
    }
  }
  
  // EBIC = n*log(loss) + log(n)*E + 4*gamma*log(p)*E
  double ebic = n * std::log(loss) + std::log(n) * E + 4.0 * gamma * std::log(p) * E;
  
  return List::create(
    Named("ebic") = ebic,
    Named("loss") = loss,
    Named("num_edges") = E
  );
}

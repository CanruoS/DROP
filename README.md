# DROP: Distributionally Robust Optimization for Multi-task Learning in Graphical Models
This is an R package for estimating Gaussian graphical models.

## Overview

**DROP** implements adaptive distributionally robust optimization (DRO) for estimating high-dimensional Gaussian graphical models. The method achieves robust graph recovery against outliers and heavy-tailed distributions through a rank-based inverse normal transformation, and efficiently solves the adaptive-weighted objective function using a scalable coordinate descent algorithm.


## Installation

You can install **DROP** from GitHub:

```r
# Install remotes if needed
install.packages("remotes")

# Install DROP from GitHub
remotes::install_github("CanruoS/DROP")

# Load the package
library(DROP)
```

Because DROP contains C++ code, users need a working R package compilation
toolchain.

### Platform-Specific Notes

#### Windows Users

**DROP** contains C++ code that must be compiled during installation. Windows requires **Rtools** for this compilation step.

**Important**: 
- **Installing DROP's dependencies** (CovTools, MASS, Rcpp, huge, etc.): NO Rtools needed. These packages download pre-compiled binaries from CRAN automatically.
- **Installing DROP itself**: Rtools IS needed because DROP is distributed as source code (.tar.gz) and must be compiled on your computer.

**One-time Rtools Setup**:
1. Download Rtools: https://cran.r-project.org/bin/windows/Rtools/
2. Install Rtools (follow the installer instructions)
3. Restart R/RStudio
4. Then install DROP

```r
install.packages("remotes")
remotes::install_github("CanruoS/DROP")
library(DROP)
```

#### Mac Users

DROP requires the Xcode Command Line Tools to compile its C++ code. If they are
not already installed, run this once in Terminal:

```bash
xcode-select --install
```

Then install DROP from GitHub:

```r
install.packages("remotes")
remotes::install_github("CanruoS/DROP")
library(DROP)
```

#### Linux Users

Ensure you have R development tools:
```bash
# Ubuntu/Debian
sudo apt-get install r-base-dev

# CentOS/RHEL  
sudo yum install R-devel
```

Then install as normal:
```r
install.packages("remotes")
remotes::install_github("CanruoS/DROP")
library(DROP)
```

### Local Installation

You can also install DROP from a local source archive:

```r
install.packages("DROP_0.1.0.tar.gz", repos = NULL, type = "source")
library(DROP)
```

### Verifying Installation

```r
library(DROP)

# Test basic functionality
data <- drop.generator(graph_type = "hub", n = 50, p = 10)
result <- drop(data$X)
```

## Usage

### Main Function: drop()

```r
library(DROP)

# Generate contaminated data using drop.generator
set.seed(123)
data_gen <- drop.generator(
  graph_type = "hub",
  n = 500,
  p = 10,
  contamination_scenario = "leverage",
  contamination_rate = 0.10,
  seed = 123
)

# Estimate precision matrix with DROP
result <- drop(X = data_gen$X, verbose = TRUE)

# Evaluate performance
eval_result <- drop.evaluate(
  estimated_theta = result$est_K,
  true_theta = data_gen$omega
)

# Inspect selected performance metrics
cat("F1 Score:", eval_result$F1, "\n")  
cat("MCC:", eval_result$MCC, "\n")       

# Access results
precision_matrix <- result$est_K
adjacency_matrix <- result$adj
```

### Data Generation: drop.generator()

Generate data with different graph structures and optional contamination:

```r
# Generate clean data with cluster graph
clean_data <- drop.generator(
  graph_type = "cluster",
  n = 200,
  p = 50,
  contamination_scenario = "clean"
)

# Generate contaminated data with hub graph
contaminated_data <- drop.generator(
  graph_type = "hub",
  n = 200,
  p = 50,
  contamination_scenario = "cauchy",
  contamination_rate = 0.15,
  seed = 123
)

# Access the data and graph structure
X <- contaminated_data$X
true_precision <- contaminated_data$omega
adjacency <- contaminated_data$theta

# Generate data with custom graph parameters
data_custom <- drop.generator(
  graph_type = "hub",
  n = 200,
  p = 50,
  g = 10,  # Use 10 hubs instead of default
  contamination_scenario = "clean"
)

# Supported graph types: "band", "hub", "cluster", "random", "scale-free"
# Contamination scenarios: "clean", "cauchy", "leverage"
```

### Add Contamination: drop.contaminate()

Add contamination to existing data:

```r
# Start with clean data from drop.generator
set.seed(123)
data_clean <- drop.generator(
  graph_type = "hub",
  n = 100,
  p = 10,
  contamination_scenario = "clean",
  seed = 123
)
clean_data <- data_clean$X

# Add Cauchy contamination
result_cauchy <- drop.contaminate(
  data = clean_data,
  contamination_scenario = "cauchy",
  contamination_rate = 0.15,
  cauchy_scale = 5,
  seed = 456
)

# Add leverage contamination
result_leverage <- drop.contaminate(
  data = clean_data,
  contamination_scenario = "leverage",
  contamination_rate = 0.15,
  leverage_multiplier = 100,
  seed = 789
)

# Access contaminated data and indices
contaminated_data <- result_cauchy$contaminated_data
contaminated_indices <- result_cauchy$contaminated_indices
```

### Evaluate Estimation: drop.evaluate()

Evaluate the accuracy of precision matrix estimation:

```r
# Generate data with known structure
set.seed(123)
data_result <- drop.generator(
  graph_type = "hub",
  n = 500,
  p = 10,
  contamination_scenario = "leverage",
  contamination_rate = 0.10,
  seed = 123
)

# Estimate precision matrix using DROP
drop_result <- drop(X = data_result$X, verbose = FALSE)

# Evaluate DROP directly
eval_drop <- drop.evaluate(drop_result$est_K, data_result$omega)

# Compare with HUGE (if installed)
if (requireNamespace("huge", quietly = TRUE)) {
  huge_result <- huge::huge(data_result$X, method = "glasso", verbose = FALSE)
  huge_select <- huge::huge.select(huge_result, criterion = "stars", verbose = FALSE)
  
  # Evaluate both methods
  eval_huge <- drop.evaluate(huge_select$refit, data_result$omega)
  
  cat("DROP F1:", eval_drop$F1, "\n")  
  cat("HUGE F1:", eval_huge$F1, "\n") 
}

# Access detailed metrics
print(eval_drop$Precision)  # Precision
print(eval_drop$Recall)     # Recall
print(eval_drop$MCC)        # Matthews Correlation Coefficient
```

### Run Simulations: drop.simulate()

Run simulation experiments to compare multiple estimation methods:

```r
# Method 1: Use pre-defined methods (character vector)
sim_result <- drop.simulate(
  n = 100,
  p = 10,
  graph_type = "hub",
  contamination_scenario = "cauchy",
  contamination_rate = 0.1,
  simulation_times = 3,
  methods = c("DROP", "HUGE_Glasso", "NPN"),
  seed = 123,
  verbose = TRUE
)

# Available pre-defined methods:
# "DROP", "HUGE_Glasso", "HUGE_MB", "NPN", "Scaled_Lasso", 
# "SCIO", "SCAD", "MCP", "TIGER", "CLIME", "Kendall", "Spearman"

# Method 2: Use custom methods (named list)
custom_methods <- list(
  DROP = function(X) {
    fit <- drop(X, verbose = FALSE)
    return(fit$est_K)
  },
  Sample_Cov = function(X) {
    S <- cov(X)
    K <- tryCatch(solve(S), error = function(e) diag(ncol(X)))
    return(K)
  }
)

sim_result2 <- drop.simulate(
  n = 100,
  p = 10,
  graph_type = "cluster",
  simulation_times = 3,
  methods = custom_methods,
  seed = 456
)

# Method 3: Specify custom graph parameters
sim_result3 <- drop.simulate(
  n = 100,
  p = 10,
  graph_type = "hub",
  g = 10,  # Use 10 hubs instead of default
  contamination_scenario = "clean",
  simulation_times = 3,
  methods = c("DROP", "HUGE_Glasso"),
  seed = 789
)

# View summary statistics
print(sim_result$summary)

# Access detailed results for all iterations
head(sim_result$results)

# Get top methods by F1 score
top_methods <- sim_result$summary[order(-sim_result$summary$F1_mean), ]
print(top_methods)
```

### Visualize Graph Structure: drop.plot()

Visualize the true graph structure from simulation results or directly from an adjacency matrix:

```r
# Method 1: From simulation results
result <- drop.simulate(
  n = 100,
  p = 10,
  graph_type = "hub",
  contamination_scenario = "clean",
  simulation_times = 3,
  methods = c("DROP", "HUGE_Glasso"),
  seed = 123
)
drop.plot(result)

# Method 2: From drop.generator output
data_gen <- drop.generator(
  graph_type = "hub",
  n = 100,
  p = 10,
  contamination_scenario = "clean",
  seed = 456
)
drop.plot(data_gen)

# Method 3: From adjacency matrix generated by drop.generator
data_gen2 <- drop.generator(
  graph_type = "band",
  n = 100,
  p = 10,
  contamination_scenario = "clean",
  seed = 789
)
drop.plot(data_gen2$theta, layout = "spiral")

# Customize appearance
drop.plot(
  data_gen,
  layout = "circle",           # Layout: "fr", "circle", "kk", "star", "grid", "spiral", "auto"
  vertex.size = 5,
  vertex.color = "red",
  edge.color = "gray",
  title = "My Custom Graph"
)

# Save to file
library(ggplot2)
p <- drop.plot(data_gen)
ggsave("true_graph.pdf", p, width = 8, height = 8)
```

**Supported inputs:**
- Simulation results from `drop.simulate()`
- Output from `drop.generator()` (automatically extracts adjacency matrix)
- Binary adjacency matrix (n × n matrix with 0/1 values)
- Sparse matrices (e.g., from Matrix package)

### Export Results to Excel: drop.excel()

Export simulation results to a formatted Excel file:

```r
# Run simulation
result <- drop.simulate(
  n = 100,
  p = 10,
  graph_type = "hub",
  contamination_scenario = "cauchy",
  contamination_rate = 0.1,
  simulation_times = 3,
  methods = c("DROP", "HUGE_Glasso", "NPN"),
  seed = 123
)

# Export to Excel (includes iteration details and summary)
drop.excel(
  result,
  filename = "simulation_results.xlsx",
  overwrite = TRUE
)

# Export with custom decimal formatting
drop.excel(
  result,
  filename = "results_2decimals.xlsx",
  digits = 2,
  overwrite = TRUE
)
```

**Excel Structure:**
- **DETAILS** sheet: All iteration results with complete metrics:
  - n, p, graph_type, contamination_scenario, simulation_times
  - iteration, method_name
  - F1, Precision, Recall, Specificity, MCC
  - Sparsity, True_Sparsity (calculated as 1 - true_edges/total_edges)
  - time (computation time in seconds)
- **SUMMARY** sheet: Aggregated statistics per method:
  - n, p, graph_type, contamination_scenario, simulation_times, method_name
  - Mean, SD, Median for each metric (F1, Precision, Recall, etc.)
  - Sorted by F1_mean (descending)

**⚠️ REQUIRED DEPENDENCY:**  
This function requires the `openxlsx` package to be installed:
```r
install.packages("openxlsx")
```
The function will show a clear error message if `openxlsx` is not installed.

### Merge Results Across Experiments: drop.merge()

Combine several outputs from `drop.simulate()` into one object while keeping
scenario-level summaries separate:

```r
result1 <- drop.simulate(
  n = 100,
  p = 10,
  graph_type = "hub",
  contamination_scenario = "clean",
  simulation_times = 3,
  methods = c("DROP", "HUGE_Glasso"),
  seed = 123
)

result2 <- drop.simulate(
  n = 100,
  p = 10,
  graph_type = "hub",
  contamination_scenario = "cauchy",
  contamination_rate = 0.1,
  simulation_times = 3,
  methods = c("DROP", "HUGE_Glasso"),
  seed = 456
)

merged_results <- drop.merge(result1, result2)
head(merged_results$summary)
```

### Generate LaTeX Tables: drop.latex()

Create LaTeX tables from simulation output stored in memory or on disk:

```r
sim_result <- drop.simulate(
  n = 100,
  p = 10,
  graph_type = "hub",
  contamination_scenario = "clean",
  simulation_times = 3,
  methods = c("DROP", "HUGE_Glasso"),
  seed = 789
)

drop.latex(
  sim_result,
  output_file = "results.tex",
  metrics = c("F1", "MCC")
)
```

## Output

### drop()

The `drop()` function returns a list with:

- `est_K`: Estimated precision matrix
- `adj`: Adjacency matrix (0/1)
- `edges`: Logical edge matrix
- `best_lambda`: Optimal lambda value
- `best_idx`: Index of best lambda
- `lambdas`: Sequence of lambda values
- `ebic_values`: EBIC values for each lambda

### drop.generator()

The `drop.generator()` function returns a list with:

- `X`: Generated data matrix (n x p)
- `graph_type`: Type of graph structure used
- `contamination_scenario`: Type of contamination applied
- `n`: Number of observations
- `p`: Number of variables
- `contamination_rate`: Proportion of contaminated samples
- `contaminated_indices`: Indices of contaminated observations (if contaminated)
- `seed`: Random seed used (if specified)
- `theta`: Adjacency matrix (if return_graph_info = TRUE)
- `omega`: Precision matrix (if return_graph_info = TRUE)
- `sigma`: Covariance matrix (if return_graph_info = TRUE)
- `true_edges`: Number of true edges (if return_graph_info = TRUE)
- `total_edges`: Total possible edges (if return_graph_info = TRUE)
- `sparsity`: Graph sparsity (if return_graph_info = TRUE)

### drop.contaminate()

The `drop.contaminate()` function returns a list with:

- `contaminated_data`: Data matrix with contamination added
- `contaminated_indices`: Indices of contaminated samples
- `contamination_scenario`: Type of contamination applied
- `contamination_rate`: Contamination rate used
- `n_contaminated`: Number of contaminated samples
- `seed`: Random seed used (if specified)

### drop.evaluate()

The `drop.evaluate()` function returns a list with:

- `F1`: F1 score (harmonic mean of precision and recall)
- `Precision`: Precision (TP / (TP + FP))
- `Recall`: Recall/Sensitivity (TP / (TP + FN))
- `Specificity`: Specificity (TN / (TN + FP))
- `MCC`: Matthews Correlation Coefficient
- `TP`: Number of true positives
- `FP`: Number of false positives
- `FN`: Number of false negatives
- `TN`: Number of true negatives
- `estimated_edges`: Number of edges in estimated graph
- `true_edges`: Number of edges in true graph
- `total_possible_edges`: Total possible edges (p*(p-1)/2)
- `estimated_sparsity`: Sparsity of estimated graph (1 - estimated_edges/total_possible_edges)
- `true_sparsity`: Sparsity of true graph (1 - true_edges/total_possible_edges)

### drop.simulate()

The `drop.simulate()` function returns a list with:

- `results`: Data frame with detailed results for each iteration and method:
  - `n`, `p`, `graph_type`, `contamination_scenario`, `simulation_times`
  - `iteration`: Simulation iteration number (1 to simulation_times)
  - `method_name`: Name of estimation method
  - `F1`: F1 score for this iteration
  - `Precision`: Precision for this iteration
  - `Recall`: Recall for this iteration
  - `Specificity`: Specificity for this iteration
  - `MCC`: Matthews Correlation Coefficient for this iteration
  - `Sparsity`: Estimated graph sparsity for this iteration
  - `True_Sparsity`: True graph sparsity for the scenario
  - `time`: Computation time in seconds
- `summary`: Data frame with aggregated statistics for each method:
  - `method_name`: Name of estimation method
  - `contamination_scenario`: Scenario represented by the summary row
  - `*_mean`: Mean value across all iterations (for each metric)
  - `*_sd`: Standard deviation across all iterations (for each metric)
  - `*_median`: Median value across all iterations (for each metric)
  - `success_rate`: Proportion of successful runs (non-NA results)
- `true_theta`: True precision matrix used in simulation
- `graph_info`: List with graph structure information:
  - `adjacency`: True adjacency matrix
  - `theta`: True precision matrix
  - `true_sparsity`: Sparsity of true graph
  - `true_edges`: Number of edges in true graph
  - `total_edges`: Total possible edges
- `params`: List of simulation parameters used (n, p, graph_type, etc.)
- `estimated_matrices`: Optional list returned when `save_matrices = TRUE`

### drop.merge()

The `drop.merge()` function returns a list with:

- `results`: Combined iteration-level results with a `scenario` column
- `summary`: Summary statistics grouped by scenario and method
- `params`: Parameter list copied from the first result for reference
- `graph_info`: Graph information copied from the first result for reference
- `scenarios`: Data frame describing all merged scenarios

### drop.latex()

The `drop.latex()` function invisibly returns the summary data frame used to
generate the LaTeX tables.

## Method

The method minimizes the objective function:

K̂ = argmin Σᵢ{(1/n)||Xᵢ - X₋ᵢβ⁽ⁱ⁾||² + λw⁽ⁱ⁾||K₋ᵢ,ᵢ||₁}

where the adaptive weights are: w⁽ⁱ⁾ = √MSE⁽ⁱ⁾/Kᵢᵢ, wᵢⱼ = w⁽ⁱ⁾ + w⁽ʲ⁾

## Contamination Types

The package supports two types of contamination for robustness testing:

- **Cauchy**: Additive contamination. Adds Cauchy-distributed noise to observations:
  - X[i,] ← X[i,] + rcauchy(p, location, scale)
  - Default: location = 0, scale = 5
  
- **Leverage**: Replacement contamination. Replaces observations with high-leverage points:
  - Generated from mvrnorm with inflated covariance (multiplier × Σ)
  - Default: multiplier = 100

## Dependencies

- MASS (required)
- CovTools (required)
- Rcpp (required)
- parallel (required)
- huge (required for graph generation and several comparison methods)
- openxlsx (optional, for `drop.excel()`)
- readxl (optional, for reading Excel files in `drop.latex()`)
- igraph, ggplot2 (optional, for `drop.plot()`)
- Matrix (optional, for sparse matrix input and rank-based comparison methods)
- scalreg, scio, GGMncv, flare (optional, for additional methods in `drop.simulate()`)

## License

GPL-3

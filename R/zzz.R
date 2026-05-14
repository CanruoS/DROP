#' Package Startup Functions
#'
#' These functions are called when the package is loaded.
#' @keywords internal
#' @noRd

.onAttach <- function(libname, pkgname) {
  # Check and auto-install required packages
  required_packages <- c("MASS", "CovTools", "Rcpp", "parallel", "huge")
  suggested_packages <- c(
    "openxlsx", "igraph", "ggplot2", "readxl", "Matrix",
    "scalreg", "scio", "GGMncv", "flare"
  )
  
  missing_required <- character(0)
  missing_suggested <- character(0)
  
  # Check required packages
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_required <- c(missing_required, pkg)
    }
  }
  
  # Check suggested packages
  for (pkg in suggested_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_suggested <- c(missing_suggested, pkg)
    }
  }
  
  # Show startup message
  packageStartupMessage("DROP: Distributionally Robust Optimization for Precision Matrices")
  packageStartupMessage("Version: ", utils::packageVersion("DROP"))
  
  # Auto-install missing required packages
  if (length(missing_required) > 0) {
    packageStartupMessage("\n*** Auto-installing missing required packages ***")
    packageStartupMessage("Missing: ", paste(missing_required, collapse = ", "))
    
    # Note: Installing from CRAN uses pre-compiled binaries on Windows, no Rtools needed
    for (pkg in missing_required) {
      packageStartupMessage("Installing ", pkg, "...")
      tryCatch({
        utils::install.packages(pkg, repos = "https://cran.r-project.org", 
                                dependencies = TRUE, quiet = FALSE)
        packageStartupMessage("  -> ", pkg, " installed successfully!")
      }, error = function(e) {
        packageStartupMessage("  -> Failed to install ", pkg, ": ", e$message)
      })
    }
    
    # Re-check if installation was successful
    still_missing <- character(0)
    for (pkg in missing_required) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        still_missing <- c(still_missing, pkg)
      }
    }
    
    if (length(still_missing) > 0) {
      packageStartupMessage("\n*** WARNING: Some packages failed to install ***")
      packageStartupMessage("Still missing: ", paste(still_missing, collapse = ", "))
      packageStartupMessage("Please install manually: install.packages(c('", 
                           paste(still_missing, collapse = "', '"), "'))")
    } else {
      packageStartupMessage("\nAll required packages are now installed!")
    }
  }
  
  # Show info about optional packages (don't auto-install)
  if (length(missing_suggested) > 0) {
    packageStartupMessage("\nOptional packages for additional features:")
    packageStartupMessage("  ", paste(missing_suggested, collapse = ", "))
    packageStartupMessage("Install with: install.packages(c('", 
                         paste(missing_suggested, collapse = "', '"), "'))")
  }
  
  if (length(missing_required) == 0 && length(missing_suggested) == 0) {
    packageStartupMessage("\nAll listed dependencies are available.")
  }
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("x", "y", "x1", "y1", "x2", "y2"))
}

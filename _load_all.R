# ============================================================================
# MASTER LOADER: Sources all functions and configuration
# ============================================================================

suppressWarnings({
  
  # LOAD REQUIRED LIBRARIES FIRST
  # ============================================================================
  library(tidyverse, quietly = TRUE)
  library(ggplot2, quietly = TRUE)
  library(dplyr, quietly = TRUE)
  library(readr, quietly = TRUE)
  library(here, quietly = TRUE)
  
  # Bioconductor packages for RNA-seq (optional, load if available)
  if (requireNamespace("DESeq2", quietly = TRUE)) {
    suppressMessages(library(DESeq2, quietly = TRUE))
  }
  if (requireNamespace("edgeR", quietly = TRUE)) {
    suppressMessages(library(edgeR, quietly = TRUE))
  }
  if (requireNamespace("limma", quietly = TRUE)) {
    suppressMessages(library(limma, quietly = TRUE))
  }
  
  # ============================================================================
  # LOAD CONFIGURATION - THIS MUST COME FIRST
  # ============================================================================
  config_path <- here::here("R", "config.R")
  if (!file.exists(config_path)) {
    stop("Cannot find config.R at: ", config_path)
  }
  source(config_path, local = FALSE)  # Use FALSE to make variables global
  
  # ============================================================================
  # LOAD FUNCTION MODULES
  # ============================================================================
  source(here::here("R", "functions_analysis.R"), local = FALSE)      
  source(here::here("R", "functions_data_rnaseq.R"), local = FALSE)
  source(here::here("R", "functions_plots.R"), local = FALSE)         
  
  # ============================================================================
  # CREATE OUTPUT DIRECTORIES (Now that PATH variables are defined)
  # ============================================================================
  # Only create if PATH variables exist
  if (exists("PATH_DATA_RAW")) {
    dir.create(here::here(PATH_DATA_RAW), showWarnings = FALSE, recursive = TRUE)
  }
  if (exists("PATH_DATA_PROCESSED")) {
    dir.create(here::here(PATH_DATA_PROCESSED), showWarnings = FALSE, recursive = TRUE)
  }
  if (exists("PATH_OUTPUT_FIGURES")) {
    dir.create(here::here(PATH_OUTPUT_FIGURES), showWarnings = FALSE, recursive = TRUE)
  }
  if (exists("PATH_OUTPUT_TABLES")) {
    dir.create(here::here(PATH_OUTPUT_TABLES), showWarnings = FALSE, recursive = TRUE)
  }
  
})

cat("\n")
cat("╔════════════════════════��═══════════════════════════════════════════╗\n")
cat("║  ✓ Configuration and all functions loaded successfully            ║\n")
cat("║                                                                    ║\n")
cat("║  Available modules:                                                ║\n")
cat("║    • Analysis functions (general + inflammatory)                  ║\n")
cat("║    • RNA-seq data processing functions                            ║\n")
cat("║    • Plotting & visualization (general + RNA-seq)                 ║\n")
cat("║                                                                    ║\n")
cat("║  Project directories verified                                      ║\n")
cat("╚════════════════════════════════════════════════════════════════════╝\n")
cat("\n")
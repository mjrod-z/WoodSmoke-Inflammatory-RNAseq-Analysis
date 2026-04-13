# ============================================================================
# RNA-SEQ DATA PROCESSING FUNCTIONS
# ============================================================================

# ============================================================================
# 1. NORMALIZE COUNTS (MEDIAN-OF-RATIOS)
# ============================================================================
#' Normalize RNA-seq counts using median-of-ratios method
#'
#' @param count_matrix Matrix with genes as rows, samples as columns
#'
#' @return Normalized count matrix
#'
normalize_counts <- function(count_matrix) {
  
  # Calculate size factors (geometric mean approach)
  size_factors <- colSums(count_matrix) / mean(colSums(count_matrix))
  
  # Normalize
  normalized <- sweep(count_matrix, 2, size_factors, "/")
  
  return(normalized)
}

# ============================================================================
# 2. FILTER GENES (REMOVE LOW-ABUNDANCE)
# ============================================================================
#' Filter out genes with low counts (background genes)
#'
#' @param data Data frame with Geneid column and count columns
#' @param background_threshold Proportion threshold (default 0.20)
#'
#' @return Filtered data frame
#'
filter_genes <- function(data, background_threshold = 0.20) {
  require(tidyverse)
  
  # Extract counts (exclude Geneid column)
  counts <- data %>%
    select(-Geneid) %>%
    as.matrix()
  
  # Calculate median count per gene
  median_counts <- apply(counts, 1, median)
  overall_median <- median(median_counts)
  
  # Filter: keep genes with median > threshold * overall_median in >20% of samples
  keep_genes <- apply(counts > (overall_median * background_threshold), 1, function(x) sum(x) > ncol(counts) * 0.20)
  
  filtered_data <- data[keep_genes, ]
  
  return(filtered_data)
}

# ============================================================================
# 3. PREPARE METADATA FOR RNA-SEQ
# ============================================================================
#' Extract and prepare metadata for samples in RNA-seq data
#'
#' @param metadata_full Full metadata data frame
#' @param sample_ids Vector of sample IDs in RNA-seq dataset
#'
#' @return Filtered metadata matching RNA-seq samples
#'
prepare_metadata_seq <- function(metadata_full, sample_ids) {
  require(tidyverse)
  
  # Match column names (handle different naming conventions)
  id_cols <- c("SAMPLEID", "Sample_ID", "SampleID")
  id_col <- id_cols[which(id_cols %in% names(metadata_full))][1]
  
  metadata_seq <- metadata_full %>%
    filter(!!sym(id_col) %in% sample_ids) %>%
    mutate_if(is.character, as.factor)
  
  return(metadata_seq)
}

# ============================================================================
# 4. CLASSIFY DEGS (UP/DOWN/NO)
# ============================================================================
#' Classify genes as UP, DOWN, or NO based on FC and p-value
#'
#' @param data Data frame with log2FoldChange and padj columns
#' @param fc_cutoff Log2 fold-change threshold (default 1.0)
#' @param p_cutoff Adjusted p-value threshold (default 0.05)
#'
#' @return Data frame with DEG_Status column
#'
classify_degs <- function(data, fc_cutoff = LOG2FC_CUTOFF, p_cutoff = ADJ_P_CUTOFF) {
  data %>%
    mutate(DEG = case_when(
      adj.P.Val <= p_cutoff & log2FC >=  fc_cutoff ~ "UP",
      adj.P.Val <= p_cutoff & log2FC <= -fc_cutoff ~ "DOWN",
      TRUE ~ "NO"
    ))
}

# ============================================================================
# 5. READ GSEA HTML OUTPUT
# ============================================================================
#' Parse GSEA HTML output files
#'
#' @param html_file Path to GSEA HTML report
#'
#' @return Data frame with pathway results
#'
read_gsea_html <- function(html_file) {
  require(rvest)
  
  page <- read_html(html_file)
  
  # Extract table (structure varies by GSEA version)
  tables <- html_table(page)
  
  if (length(tables) > 0) {
    results <- tables[[1]]
  } else {
    results <- data.frame()
  }
  
  return(results)
}

# ============================================================================
# 6. FILTER GSEA RESULTS
# ============================================================================
#' Filter and annotate GSEA results by FDR
#'
#' @param gsea_results Data frame with GSEA pathway results
#' @param fdr_cutoff FDR threshold (default 0.25)
#'
#' @return Filtered results with Significance column
#'
filter_gsea_results <- function(gsea_results, fdr_cutoff = 0.25) {
  
  results_filtered <- gsea_results %>%
    mutate(Significant = ifelse(FDR < fdr_cutoff, "Yes", "No")) %>%
    arrange(FDR)
  
  return(results_filtered)
}

# ============================================================================
# 7. PREPARE DEGS FOR EXPORT (IPA FORMAT)
# ============================================================================
#' Format DEG data for import to IPA or similar tools
#'
#' @param deg_data Data frame with gene names, log2FC, padj
#'
#' @return Data frame in IPA-compatible format
#'
prepare_deg_export <- function(deg_data) {
  require(tidyverse)
  
  export_data <- deg_data %>%
    select(Gene = gene_name, log2FoldChange, padj) %>%
    mutate(
      FC = 2^log2FoldChange,
      log10_padj = -log10(padj)
    ) %>%
    filter(padj < 0.05) %>%
    arrange(desc(abs(log2FoldChange)))
  
  return(export_data)
}

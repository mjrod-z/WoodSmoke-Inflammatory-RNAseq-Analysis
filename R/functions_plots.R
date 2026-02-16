# ============================================================================
# PLOTTING & VISUALIZATION FUNCTIONS
# ============================================================================
# Consolidated module for all plotting functions (general + RNA-seq specific)
# Source via: source("R/_load_all.R")
#
# Sections:
#   1. Theme & General Utilities
#   2. General Plots (Boxplot, PCA, MDS)
#   3. RNA-Seq Plots (Volcano, UpSet, Gene Expression)
#   4. GSEA Plots (Barplot, Dotplot, Heatmap)
#   5. Save Functions & Color Scales

library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)

# ============================================================================
# 1. THEME & GENERAL UTILITIES
# ============================================================================

#' Apply consistent project theme to ggplot
#'
#' @return ggplot theme object
#'
theme_project <- function() {
  theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      strip.text = element_text(size = 12, face = "bold")
    )
}

# ============================================================================
# 2. GENERAL PLOTS
# ============================================================================

#' Create boxplot with jittered points by exposure
#'
#' @param data Data frame with EXPOSURE and measurement column
#' @param measure_col Column name with measurements
#' @param fill_color Color for fill (default automatic)
#'
#' @return ggplot object
#'
plot_exposure_boxplot <- function(data, measure_col, fill_color = NULL) {
  require(ggplot2)
  
  p <- ggplot(data, aes(x = EXPOSURE, y = !!sym(measure_col), fill = EXPOSURE)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 2) +
    theme_project() +
    labs(x = "Exposure", y = measure_col, title = paste("Distribution of", measure_col)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  if (!is.null(fill_color)) {
    p <- p + scale_fill_manual(values = fill_color)
  }
  
  return(p)
}

#' Create PCA scatter plot with ellipses by group
#'
#' @param pca_result prcomp object
#' @param metadata Data frame with group information
#' @param group_var Column name for grouping
#' @param title Plot title
#'
#' @return ggplot object
#'
plot_pca <- function(pca_result, metadata, group_var, title = "PCA") {
  require(ggplot2)
  
  pca_data <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    Group = metadata[[group_var]]
  )
  
  var_exp <- summary(pca_result)$importance[2, 1:2]
  
  p <- ggplot(pca_data, aes(x = PC1, y = PC2, fill = Group, color = Group)) +
    geom_point(size = 3, alpha = 0.7) +
    stat_ellipse(geom = "polygon", alpha = 0.2, show.legend = FALSE) +
    theme_project() +
    labs(
      x = paste0("PC1 (", round(var_exp[1]*100, 1), "%)"),
      y = paste0("PC2 (", round(var_exp[2]*100, 1), "%)"),
      title = title
    )
  
  return(p)
}

#' Create MDS scatter plot
#'
#' @param lcpm Log-CPM matrix (samples as columns)
#' @param group_col Group vector for coloring
#' @param title Plot title
#'
#' @return ggplot object
#'
plot_mds <- function(lcpm, group_col, title = "MDS Plot") {
  require(ggplot2)
  require(limma)
  
  mds_result <- cmdscale(dist(t(lcpm)), k = 2)
  
  mds_data <- data.frame(
    MDS1 = mds_result[, 1],
    MDS2 = mds_result[, 2],
    Group = group_col
  )
  
  p <- ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Group, fill = Group)) +
    geom_point(size = 3, alpha = 0.7) +
    theme_project() +
    labs(x = "MDS Dimension 1", y = "MDS Dimension 2", title = title)
  
  return(p)
}

# ============================================================================
# 3. RNA-SEQ PLOTS
# ============================================================================

#' Create volcano plot for differential expression analysis
#'
#' @param data Data frame with log2FC, adj.P.Val, DEG, and delabel columns
#' @param facet_by Column name for faceting (default "sample_name")
#' @param title Plot title
#' @param x_limits X-axis limits (default c(-10, 10))
#' @param y_limits Y-axis limits (default c(0, 50))
#'
#' @return ggplot object
#'
plot_volcano_deg <- function(data, facet_by = "sample_name", title = NULL,
                             x_limits = c(-10, 10), y_limits = c(0, 50)) {
  require(ggplot2)
  require(ggrepel)
  
  p <- ggplot(data, aes(x = log2FC, y = -log10(adj.P.Val), color = DEG, label = delabel)) +
    geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed', size = 0.5) +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed', size = 0.5) +
    geom_point(size = 2, alpha = 0.7) +
    geom_text_repel(
      data = subset(data, !is.na(delabel)),
      size = 3,
      max.overlaps = 20,
      fontface = "bold",
      box.padding = 0.4,
      segment.color = "black"
    ) +
    scale_color_manual(
      values = c("DOWN" = "#00AFBB", "NO" = "grey", "UP" = "#bb0c00"),
      labels = c("Downregulated", "Not significant", "Upregulated")
    ) +
    coord_cartesian(ylim = y_limits, xlim = x_limits) +
    labs(
      x = expression("log"[2] * "FC"),
      y = expression("-log"[10] * "adj.P.Val"),
      title = title
    ) +
    theme_minimal(base_size = 14) +
    theme(
      strip.text = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_line(color = "gray80"),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "#ffffff", color = NA),
      legend.background = element_rect(fill = "#ffffff", color = NA),
      legend.key = element_rect(fill = "#ffffff", color = NA)
    ) +
    facet_wrap(as.formula(paste("~", facet_by)), scales = "free_y")
  
  return(p)
}

#' Create volcano plot with sex × fuel type faceting (with background colors)
#'
#' @param data Data frame with Sex, sample_name, log2FC, adj.P.Val, DEG, delabel
#' @param selected_genes Vector of gene names to highlight
#' @param sex_choice "Male" or "Female"
#'
#' @return ggplot object
#'
plot_volcano_deg_sex <- function(data, selected_genes, sex_choice) {
  require(ggplot2)
  require(ggrepel)
  
  # Set factor order (Male on top)
  data$Sex <- factor(data$Sex, levels = c("Male", "Female"))
  
  # Filter data based on the specified sex
  filtered_data <- subset(data, Sex == sex_choice)
  
  # Specify the background color for the sex chosen
  background_color <- ifelse(sex_choice == "Female", "#fff0f3", "#F0FFF0")
  
  # Count UP and DOWN genes for each sample_name
  up_counts <- filtered_data %>%
    dplyr::filter(DEG == "UP") %>%
    dplyr::group_by(Sex, sample_name) %>%
    dplyr::summarise(n = n(), .groups = "drop")
  
  down_counts <- filtered_data %>%
    dplyr::filter(DEG == "DOWN") %>%
    dplyr::group_by(Sex, sample_name) %>%
    dplyr::summarise(n = n(), .groups = "drop")
  
  # Generate the plot
  p <- ggplot(data = filtered_data, aes(x = log2FC, y = -log10(adj.P.Val), color = DEG, label = delabel)) +
    geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
    
    # Plot all genes first
    geom_point(alpha = 0.7, size = 2) +
    
    # Overlay labeled genes (highlight them with black outlines)
    geom_point(data = subset(filtered_data, !is.na(delabel) | X %in% selected_genes),
               aes(x = log2FC, y = -log10(adj.P.Val)),
               color = "black", size = 2, shape = 21, stroke = 1) +
    
    # Label all genes with a delabel, including both labeled and selected genes
    geom_label_repel(data = subset(filtered_data, !is.na(delabel) | X %in% selected_genes),
                     aes(label = ifelse(!is.na(delabel), delabel, X)),
                     size = 3, fontface = "bold", fill = "white", color = "black",
                     box.padding = 0.4, segment.color = "black", max.overlaps = 20) +
    
    # Annotate UP count at the top-right corner of each graph
    geom_text(data = up_counts, aes(x = 3.5, y = 20.5, label = paste0("UP: ", n)),
              color = "#bb0c00", size = 3, fontface = "bold", hjust = 1, inherit.aes = FALSE) +
    
    # Annotate DOWN count at the top-left corner of each graph
    geom_text(data = down_counts, aes(x = -4.5, y = 20.5, label = paste0("DOWN: ", n)),
              color = "#00AFBB", size = 3, fontface = "bold", hjust = 0, inherit.aes = FALSE) +
    
    scale_color_manual(
      values = c("DOWN" = "#00AFBB", "NO" = "grey", "UP" = "#bb0c00"),
      labels = c("Downregulated", "Not significant", "Upregulated")
    ) +
    coord_cartesian(ylim = c(0, 20), xlim = c(-8, 8)) +
    labs(
      x = expression("log"[2] * "FC"),
      y = expression("-log"[10] * "adj.P.Val")
    ) +
    theme_minimal(base_size = 14) +
    theme(
      strip.text = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
      panel.background = element_rect(fill = background_color, color = NA),  # BACKGROUND COLOR APPLIED HERE
      plot.background = element_rect(fill = "#ffffff", color = NA),
      legend.background = element_rect(fill = "#ffffff", color = NA),
      legend.key = element_rect(fill = "#ffffff", color = NA),
      panel.grid.major = element_line(color = "gray80"),
      panel.grid.minor = element_line(color = "gray90")
    ) +
    facet_grid(Sex ~ sample_name, scales = "free_y")
  
  return(p)
}

#' Create UpSet plot showing DEG overlap across fuel types
#'
#' @param data Wide-format data frame with gene presence/absence by fuel type
#' @param title Plot title
#' @param fuel_colors Named vector of fuel type colors
#'
#' @return UpSet plot object
#'
plot_upset_degs <- function(data, title = "DEG Overlap", fuel_colors = NULL) {
  require(ComplexUpset)
  
  if (is.null(fuel_colors)) {
    fuel_colors <- c(
      "Red Oak 25 µg/cm²" = "brown3",
      "Pine 25 µg/cm²" = "deepskyblue3",
      "Peat 25 µg/cm²" = "gold3",
      "Eucalyptus 25 µg/cm²" = "darkorchid"
    )
  }
  
  fuel_sets <- names(fuel_colors)
  
  p <- upset(
    data,
    intersect = fuel_sets,
    set_sizes = upset_set_size(
      geom = geom_bar(aes(fill = group), width = 0.6)
    ),
    queries = list(
      upset_query(set = "Red Oak 25 µg/cm²", color = fuel_colors["Red Oak 25 µg/cm²"], fill = fuel_colors["Red Oak 25 µg/cm²"]),
      upset_query(set = "Pine 25 µg/cm²", color = fuel_colors["Pine 25 µg/cm²"], fill = fuel_colors["Pine 25 µg/cm²"]),
      upset_query(set = "Peat 25 µg/cm²", color = fuel_colors["Peat 25 µg/cm²"], fill = fuel_colors["Peat 25 µg/cm²"]),
      upset_query(set = "Eucalyptus 25 µg/cm²", color = fuel_colors["Eucalyptus 25 µg/cm²"], fill = fuel_colors["Eucalyptus 25 µg/cm²"])
    ),
    base_annotations = list(
      'Intersection size' = intersection_size(
        counts = TRUE,
        bar_number_threshold = 1,
        text = list(size = 3)
      )
    ),
    matrix = (
      intersection_matrix(
        geom = geom_point(aes(color = group), size = 3, shape = 19),
        segment = geom_segment(color = "black", size = 0.5),
        outline_color = list(active = "black", inactive = "grey70")
      )
      + scale_color_manual(values = fuel_colors, breaks = fuel_sets)  # EXPLICIT COLOR MAPPING
    ),
    themes = upset_default_themes(text = element_text(size = 12))
  ) +
    scale_fill_manual(values = fuel_colors, name = "Fuel Type") +
    ggtitle(title) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right"
    )
  
  return(p)
}

#' Create barplot of gene expression by exposure and sex
#'
#' @param vst_data VST-transformed count data with Geneid column
#' @param metadata Metadata with SAMPLEID, EXPOSURE, SEX, PATIENTCODE
#' @param genes Vector of gene names to plot
#'
#' @return ggplot object
#'
plot_gene_expression_bars <- function(vst_data, metadata, genes) {
  require(ggplot2)
  
  # Filter and reshape data
  plot_data <- vst_data %>%
    filter(Geneid %in% genes) %>%
    pivot_longer(cols = -Geneid, names_to = "SAMPLEID", values_to = "expression") %>%
    left_join(metadata, by = "SAMPLEID") %>%
    filter(grepl("25", EXPOSURE) | EXPOSURE == "PBS_Control")
  
  # Calculate fold change relative to PBS
  plot_data <- plot_data %>%
    group_by(PATIENTCODE, Geneid) %>%
    mutate(log2FC = log2(expression / expression[EXPOSURE == "PBS_Control"])) %>%
    ungroup() %>%
    filter(EXPOSURE != "PBS_Control")
  
  # Summary statistics
  summary_data <- plot_data %>%
    group_by(EXPOSURE, Geneid, SEX) %>%
    summarise(
      mean_expression = mean(log2FC, na.rm = TRUE),
      sd_expression = sd(log2FC, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Exposure colors
  exposure_colors <- c(
    "Eucalyptus_25" = "darkorchid",
    "Peat_25" = "gold2",
    "Pine_25" = "deepskyblue3",
    "RedOak_25" = "brown3"
  )
  
  p <- ggplot(summary_data, aes(x = EXPOSURE, y = mean_expression,
                                fill = interaction(SEX, EXPOSURE),
                                color = EXPOSURE)) +
    geom_bar(stat = "identity", position = "dodge", size = 1.2, width = 0.7,
             alpha = ifelse(summary_data$SEX == "M", 0.6, 1)) +
    geom_errorbar(aes(ymin = mean_expression - sd_expression,
                      ymax = mean_expression + sd_expression),
                  position = position_dodge(width = 0.7), width = 0.25, size = 1.2) +
    geom_jitter(data = plot_data,
                aes(x = EXPOSURE, y = log2FC, color = EXPOSURE),
                position = position_dodge(width = 0.8), alpha = 0.8, size = 1.5) +
    scale_fill_manual(values = c(
      "M.Eucalyptus_25" = scales::alpha("darkorchid", 0.6), "F.Eucalyptus_25" = "white",
      "M.Peat_25" = scales::alpha("gold", 0.6), "F.Peat_25" = "white",
      "M.Pine_25" = scales::alpha("deepskyblue3", 0.6), "F.Pine_25" = "white",
      "M.RedOak_25" = scales::alpha("brown3", 0.6), "F.RedOak_25" = "white"
    )) +
    scale_color_manual(values = exposure_colors) +
    facet_wrap(~ Geneid, scales = "free_y", nrow = 1) +
    labs(x = NULL, y = "log2FC") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(size = 12, face = "bold"),
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(labels = function(x) gsub("_25", "", x))
  
  return(p)
}

# ============================================================================
# 4. GSEA PLOTS
# ============================================================================

#' Create barplot of GSEA pathway enrichment scores
#'
#' @param gsea_data Data frame with Pathway, NES, FDR_qval columns
#' @param n_top Number of top pathways to show per group
#' @param facet_by Column for faceting
#'
#' @return ggplot object
#'
plot_gsea_barplot <- function(gsea_data, n_top = 12, facet_by = "sample_name") {
  require(ggplot2)
  
  # Select top pathways
  plot_data <- gsea_data %>%
    group_by(!!sym(facet_by)) %>%
    arrange(FDR_qval) %>%
    slice_head(n = n_top) %>%
    ungroup() %>%
    mutate(
      Path = gsub("^HALLMARK_", "", Pathway),
      Path = gsub("_", " ", Path)
    )
  
  p <- ggplot(plot_data, aes(x = reorder(Path, NES), y = NES, fill = FDR_qval)) +
    geom_col() +
    scale_fill_gradient(low = "darkblue", high = "gray") +
    geom_hline(yintercept = 0, color = "black", size = 0.5) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(face = "bold", angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(face = "bold"),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      panel.grid.major.x = element_blank()
    ) +
    coord_flip() +
    labs(
      x = "Pathway",
      y = "Normalized Enrichment Score",
      title = "Hallmark Pathways Enrichment",
      fill = "FDR qval"
    ) +
    facet_wrap(as.formula(paste("~", facet_by)), scales = "free_y")
  
  return(p)
}

#' Create dot plot of GSEA results (size = |NES|, color = NES)
#'
#' @param gsea_data Data frame with sample_name, Pathway, NES, FDR_qval
#'
#' @return ggplot object
#'
plot_gsea_dotplot <- function(gsea_data) {
  require(ggplot2)
  
  # Clean pathway names
  gsea_data <- gsea_data %>%
    mutate(
      clean_path = gsub("^HALLMARK_", "", Pathway),
      clean_path = gsub("_", " ", clean_path)
    )
  
  p <- ggplot(gsea_data, aes(x = sample_name, y = clean_path, 
                             size = abs(NES), color = NES)) +
    geom_point() +
    scale_size_continuous(range = c(2, 10), limits = c(1, 2.4)) +
    scale_color_gradient2(
      low = "blue", mid = "white", high = "red", midpoint = 0,
      limits = c(-2.3, 2.3),
      name = "NES"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      panel.grid.major = element_line(color = "black"),
      panel.grid.minor = element_line(color = "black")
    ) +
    labs(
      title = "GSEA Results",
      x = "Fuel Type",
      y = "Pathway",
      size = "Absolute NES",
      color = "NES"
    ) +
    guides(color = guide_colorbar(title = "NES", title.position = "top"))
  
  return(p)
}

#' ComplexHeatmap of pathways × samples with GSEA results
#'
#' @param pathway_data Data frame with pathway names and sample enrichment scores
#' @param title Heatmap title
#'
#' @return ComplexHeatmap object
#'
plot_gsea_heatmap <- function(pathway_data, title = "GSEA Pathway Heatmap") {
  require(ComplexHeatmap)
  require(circlize)
  
  matrix_data <- as.matrix(pathway_data[, -1])
  rownames(matrix_data) <- pathway_data[[1]]
  
  ht <- Heatmap(
    matrix_data,
    name = "Enrichment Score",
    column_title = title,
    col = colorRamp2(c(min(matrix_data), 0, max(matrix_data)), 
                     c("blue", "white", "red")),
    clustering_distance_rows = "euclidean"
  )
  
  return(ht)
}

# ============================================================================
# 5. SAVE FUNCTIONS & COLOR SCALES
# ============================================================================

#' Save plot with standardized settings
#'
#' @param filename Output filename
#' @param plot ggplot object
#' @param width Width in inches (default 10)
#' @param height Height in inches (default 8)
#' @param dpi Resolution (default 300)
#' @param path Output path (default "output/figures/")
#' @param bg Background color (default "white")
#'
#' @return Invisibly TRUE
#'
save_plot <- function(filename, plot, width = 10, height = 8, dpi = 300, 
                      path = "output/figures/", bg = "white") {
  require(ggplot2)
  
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  
  ggsave(
    filename = file.path(path, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    units = "in",
    bg = bg
  )
  
  cat("✓ Saved:", file.path(path, filename), "\n")
  invisible(TRUE)
}

#' Save table with standardized settings
#'
#' @param data Data frame to save
#' @param filename Output filename
#' @param path Output path (default "output/tables/")
#'
#' @return Invisibly TRUE
#'
save_table <- function(data, filename, path = "output/tables/") {
  require(readr)
  
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  
  write_csv(data, file = file.path(path, filename))
  
  cat("✓ Saved:", file.path(path, filename), "\n")
  invisible(TRUE)
}

# ============================================================================
# COLOR SCALE FUNCTIONS (FOR GGPLOT)
# ============================================================================

#' Apply exposure colors to ggplot
#' @param use_fill Use fill instead of color (default FALSE)
scale_exposure_colors <- function(use_fill = FALSE) {
  if (use_fill) {
    scale_fill_manual(values = color_exposure, name = "Exposure")
  } else {
    scale_color_manual(values = color_exposure, name = "Exposure")
  }
}

#' Apply sex colors to ggplot
#' @param use_fill Use fill instead of color (default FALSE)
scale_sex_colors <- function(use_fill = FALSE) {
  if (use_fill) {
    scale_fill_manual(values = color_sex, name = "Sex")
  } else {
    scale_color_manual(values = color_sex, name = "Sex")
  }
}

#' Apply DEG status colors to ggplot
#' @param use_fill Use fill instead of color (default FALSE)
scale_deg_colors <- function(use_fill = FALSE) {
  if (use_fill) {
    scale_fill_manual(values = color_deg, name = "Status")
  } else {
    scale_color_manual(values = color_deg, name = "Status")
  }
}

cat("✓ Plotting functions loaded (general + RNA-seq + GSEA)\n")
# Configuration settings for Inflammatory and RNA-seq analysis

# Statistical thresholds
alpha_level <- 0.05  # Significance level
beta_level <- 0.8     # Power level

# Color palettes
color_palette_inflammatory <- c("#FF5733", "#33FF57", "#3357FF") # Example colors
color_palette_rnaseq <- c("#C70039", "#900C3F", "#581845") # Example colors

# LLOD reference data
llod_reference <- list(
  inflammatory = 10,  # Lower limit of detection for inflammatory markers
  rnaseq = 5          # Lower limit of detection for RNA-seq
)

# Plot settings
plot_settings <- list(
  xlim = c(0, 100),  # Set x-axis limits
  ylim = c(0, 1),    # Set y-axis limits
  point_size = 3,     # Size of points
  line_width = 1.5    # Width of lines
)

# File paths
file_paths <- list(
  data_inflammatory = "data/inflammatory_data.csv",
  data_rnaseq = "data/rnaseq_data.csv",
  output_directory = "results/"
)
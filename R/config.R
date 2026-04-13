# ============================================================================
# PROJECT CONFIGURATION & CONSTANTS
# ============================================================================
# All hardcoded values, thresholds, color palettes, and data references

# ============================================================================
# STATISTICAL THRESHOLDS
# ============================================================================
ZERO_CUTOFF <- 0.30                # For imputation (proportion of zeros)
ALPHA_SIGNIFICANCE <- 0.9           # For Grubbs test
LOG2FC_CUTOFF <- 1.0                # For DEG classification
ADJ_P_CUTOFF <- 0.051               # For p-value significance
FDR_THRESHOLD <- 0.25               # For GSEA pathways
FDR_THRESHOLD_HIGH <- 0.05          # For stricter GSEA filtering
SHAPIRO_P_THRESHOLD <- 0.05         # For normality testing
NES_CUTOFF <- 0.6                   # For GSEA NES cutoff

# ============================================================================
# COLOR PALETTES
# ============================================================================

# Exposure/Fuel Type Colors
color_exposure <- c(
  "Eucalyptus_25" = "darkorchid",
  "Eucalyptus_5" = "darkorchid",
  "Peat_25" = "gold2",
  "Peat_5" = "gold2",
  "Pine_25" = "deepskyblue3",
  "Pine_5" = "deepskyblue3",
  "Red_Oak_25" = "brown3",
  "Red_Oak_5" = "brown3",
  "PBS_Control" = "gray40",
  "Eucalyptus 25" = "darkorchid",
  "Peat 25" = "gold2",
  "Pine 25" = "deepskyblue3",
  "Red Oak 25" = "brown3",
  "PBS Control" = "gray40",
  "Untreated Control" = "forestgreen"
)

# Sex Colors
color_sex <- c(
  "M" = "lemonchiffon",
  "F" = "thistle1",
  "Male" = "lemonchiffon",
  "Female" = "thistle1"
)

# DEG Status Colors
color_deg <- c(
  "UP" = "#bb0c00",
  "DOWN" = "#00AFBB",
  "NO" = "grey"
)

# Significance Colors
color_significance <- c(
  "significant" = "red",
  "non-significant" = "grey"
)

# ============================================================================
# LLOD (Lower Limit of Detection) REFERENCE - V-PLEX 29-Plex Kit
# ============================================================================
cytokine_llod <- data.frame(
  Analyte = c(
    "EOTAXIN", "EOTAXIN3", "GMCSF", "IFNY", "IL1A", "IL1B", "IL2", "IL4", "IL5", "IL6",
    "IL7", "IL10", "IL12P40", "IL12P70", "IL13", "IL15", "IL16", "IL17", "IP10", "MCP1",
    "MCP4", "MDC", "MIP1A", "MIP1B", "TARC", "TNFA", "TNFB", "VEGF", "IL8"
  ),
  LLOD = c(
    0.20, 1.44, 0.16, 0.37, 0.09, 0.05, 0.09, 0.02, 0.14, 0.06,
    0.12, 0.07, 0.04, 0.33, 0.11, 0.24, 0.15, 2.83, 0.31, 0.12,
    0.16, 0.05, 1.00, 0.05, 0.09, 0.11, 0.04, 0.08, 1.12
  )
)

# ============================================================================
# PLOT DIMENSIONS & SETTINGS
# ============================================================================
PLOT_DPI <- 300
PLOT_WIDTH_SINGLE <- 6
PLOT_HEIGHT_SINGLE <- 4
PLOT_WIDTH_MULTI <- 14
PLOT_HEIGHT_MULTI <- 8
PLOT_WIDTH_WIDE <- 16
PLOT_HEIGHT_TALL <- 10
PLOT_FONT_BASE <- 14

# ============================================================================
# FILE PATHS (relative to project root)
# ============================================================================
PATH_DATA_RAW <- "data/raw/"
PATH_DATA_PROCESSED <- "data/processed/"
PATH_OUTPUT_FIGURES <- "output/figures/"
PATH_OUTPUT_TABLES <- "output/tables/"
PATH_R_FUNCTIONS <- "R/"

# ============================================================================
# ANALYSIS PARAMETERS
# ============================================================================

# DEG highlighting genes
highlighted_genes <- c("SPDEF", "IL1A", "IL1B", "SRXN1", "CYP1A1", "SOX9", 
                       "PDIA2", "HIC1", "TXNRD1", "PTGER2", "GPR68", "SFRP2", "F2RL3")

# GSEA pathway list (priority pathways)
priority_pathways <- c(
  "REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "ESTROGEN_RESPONSE_EARLY",
  "ESTROGEN_RESPONSE_LATE",
  "XENOBIOTIC_METABOLISM",
  "ANDROGEN_RESPONSE",
  "OXIDATIVE_PHOSPHORYLATION",
  "WNT_BETA_CATENIN_SIGNALING",
  "GLYCOLYSIS",
  "CHOLESTEROL_HOMEOSTASIS",
  "PI3K_AKT_MTOR_SIGNALING",
  "TGF_BETA_SIGNALING",
  "TNFA_SIGNALING_VIA_NFKB",
  "NOTCH_SIGNALING"
)

# ============================================================================
# RNA-SEQ ANALYSIS PARAMETERS
# ============================================================================

# DREAM analysis settings
DREAM_N_CORES <- 4
DREAM_VOOM_PLOT <- TRUE

# Gene filtering thresholds
GENE_BACKGROUND_THRESHOLD <- 0.20  # Proportion of samples for filtering
GENE_ZERO_THRESHOLD <- 0.10        # Max proportion of zero samples

# Variance partition thresholds
VARIANCE_MIN_CONTRIBUTION <- 0.01  # Minimum variance contribution to report

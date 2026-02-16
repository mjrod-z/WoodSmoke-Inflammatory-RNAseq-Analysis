# Inflammatory & RNA-Seq Analysis - First Paper

## 📋 Overview

This project contains refactored analysis code for studying inflammatory response and differential gene expression in human bronchial epithelial cells (hBEC) exposed to biomass smoke.

**Key Features:**
- ✅ Modular function organization (config, data processing, statistics, plotting)
- ✅ Complete LLOD imputation for 29-plex cytokine panel
- ✅ RNA-seq normalization, DEG analysis (DREAM), and GSEA pathway analysis
- ✅ Statistical testing (LMER, ART ANOVA, nonparametric)
- ✅ Professional visualizations (volcano, UpSet, barplots) with consistent theming

---

## 📁 Project Structure

First Paper Code/ ├── analysis/ │ ├── 01_inflammatory_analysis.Rmd │ └── 02_rnaseq_analysis.Rmd ├── data/ │ ├── raw/ │ └── processed/ ├── output/ │ ├── figures/ │ ├── tables/ │ └── reports/ ├── R/ │ ├── _load_all.R │ ├── config.R │ ├── functions_analysis.R │ ├── functions_data_rnaseq.R │ └── functions_plots.R └── README.md

Code

---

## 🚀 Quick Start

### 1. Load Configuration & Functions

```r
source("R/_load_all.R")
```r
This automatically sources:

config.R - All constants, colors, and paths
functions_analysis.R - Data processing & statistical testing
functions_data_rnaseq.R - RNA-seq processing
functions_plots.R - All visualization functions

### 2. Run Analysis
Option A: Inflammatory Analysis

Open and knit analysis/01_inflammatory_analysis.Rmd:

Loads raw cytokine data (MSD, IL-8 ELISA, LDH, TEER)
Performs LLOD imputation
Runs LMER, ART ANOVA, and nonparametric tests
Generates summary tables with significance stars
Outputs to output/tables/ and output/reports/

Option B: RNA-Seq Analysis

Open and knit analysis/02_rnaseq_analysis.Rmd:

Gene filtering, normalization, and VST transformation
PCA and dendrogram QC
DREAM differential expression analysis (exposure + sex interaction)
DEG classification (UP/DOWN/NO)
GSEA pathway enrichment analysis
Volcano plots, UpSet plots, barplots, and dot plots

### 3. View Results
Tables: output/tables/*.csv
Figures: output/figures/*.png
Reports: output/reports/*.html

📊 Analysis Pipeline
Inflammatory Analysis (01_inflammatory_analysis.Rmd)
Data Processing:

LLOD Imputation: LLOD/√2 method for values <30% zeros (sex-stratified)
Outlier Detection: Grubbs test on PBS controls
Transformations: Log10 & sqrt with Shapiro-Wilk normality testing
IL-8 Replacement: MSD panel IL-8 replaced with ELISA data
LDH Cytotoxicity: Calculated relative to positive control
Statistical Tests:

LMER (Linear Mixed-Effects Model)

Random intercept by patient
Treatment-vs-control contrasts (FDR corrected)
Stratified by sex (All, M, F)
ART ANOVA (Aligned Rank Transform)

Non-parametric alternative to ANOVA
Handles non-normal distributions
Nonparametric Tests

Friedman test (paired) with pairwise Wilcoxon
Kruskal-Wallis test (unpaired) with pairwise Wilcoxon
Sex interaction tests

Output Files:

nonparametric_cytokine_results.csv
lmer_contrasts_all.csv
art_contrasts_all.csv
MSD_hbec_lmer_with_significance.csv
MSD_hbec_art_with_significance.csv
LDH_TEER_summary.csv

RNA-Seq Analysis (02_rnaseq_analysis.Rmd)
Data Processing:

Gene filtering (remove low-abundance genes, >20% samples above median)
Median-of-ratios normalization
Variance stabilizing transformation (VST)
Quality Control:

PCA plot (sample clustering by exposure)
Dendrogram analysis (correlation-based)
Differential Expression (DREAM):

Linear mixed-effects model with voomWithDreamWeights
Random effect: (1 | PATIENTCODE)
Fixed effects: EXPOSURE and EXPOSURE * SEX
Treatment-vs-control contrasts (PBS Control as reference)
Log2FC > 1, adj.P.Val < 0.05
GSEA Pathway Analysis:

Hallmark pathway enrichment
NES (Normalized Enrichment Score) calculation
FDR < 0.25 significance threshold
Sex-stratified and fuel-type analyses
Visualizations:

Volcano plots (fuel type + sex-stratified with background colors)
UpSet plots (DEG overlap across fuel types)
GSEA barplots (pathway enrichment by fuel/sex)
GSEA dot plots (size = |NES|, color = NES)
Euler diagrams (DEG overlap by direction)
Gene expression barplots (IL1A, SPDEF, CYP1A1, etc.)
Output Files:

all_deg_results.csv
significant_degs.csv
degs_by_fuel.csv
degs_by_sex.csv
IPA_DEGs_FUEL.csv
IPA_DEGs_SEX.csv
MJR_seq_normalized.csv
MJR_seq_vst.csv

🛠️ Key Functions
Analysis Functions (functions_analysis.R)
Data Processing:

load_cytokine_data() - Load and clean CSV files
rename_for_coding() - Standardize column names
clean_gene_names() - Format gene identifiers
impute_lod() - LLOD imputation with sex-specific zero checks
detect_outliers_grubbs() - Outlier detection
test_normality() - Shapiro-Wilk test
test_transformations() - Test log10/sqrt transforms
log10_transform(), sqrt_transform() - Apply transformations

Statistical Testing:

exposure_lmer_pairwise() - LMER with emmeans contrasts
interaction_lmer_pairwise() - LMER with interaction terms
exposure_art_pairwise() - ART ANOVA with contrasts
nonparametric_pairwise() - Friedman, Kruskal-Wallis, Wilcoxon

Fold Change & Summaries:

calculate_fold_change_PBS() - PBS-normalized fold change
summarize_to_wide() - Create mean ± SD tables
add_significance_stars() - Add * / # annotations

RNA-Seq Functions (functions_data_rnaseq.R)

normalize_counts() - Median-of-ratios normalization
filter_genes() - Remove low-abundance genes
prepare_metadata_seq() - Match metadata to samples
classify_degs() - Classify UP/DOWN/NO based on FC/padj
read_gsea_html() - Parse GSEA HTML output
filter_gsea_results() - Filter by FDR threshold
prepare_deg_export() - Format DEGs for IPA

Plotting Functions (functions_plots.R)
General Plots:

theme_project() - Consistent project theme
plot_exposure_boxplot() - Boxplot with jitter
plot_pca() - PCA scatter with ellipses
plot_mds() - MDS plot

RNA-Seq Plots:

plot_volcano_deg() - Volcano plot by fuel type
plot_volcano_deg_sex() - Volcano with sex faceting (background colors)
plot_upset_degs() - UpSet plot for DEG overlap
plot_gene_expression_bars() - Gene expression barplots

GSEA Plots:

plot_gsea_barplot() - Horizontal bars (NES)
plot_gsea_dotplot() - Dot plot (size = |NES|, color = NES)
plot_gsea_heatmap() - ComplexHeatmap

Utilities:

save_plot() / save_table() - Standardized output functions
scale_exposure_colors(), scale_sex_colors(), scale_deg_colors() - Color scales

⚙️ Configuration (config.R)
Key Parameters
Statistical Thresholds:

ZERO_CUTOFF = 0.30 - Max proportion of zeros for imputation
ADJ_P_CUTOFF = 0.051 - Significance threshold
LOG2FC_CUTOFF = 1.0 - Log2 fold-change cutoff for DEGs
FDR_THRESHOLD = 0.25 - GSEA FDR cutoff
SHAPIRO_P_THRESHOLD = 0.05 - Normality test threshold

RNA-Seq Parameters:

DREAM_N_CORES = 4 - Parallel processing cores
GENE_BACKGROUND_THRESHOLD = 0.20 - Gene filtering threshold
GENE_ZERO_THRESHOLD = 0.10 - Max zero proportion
LLOD Reference

29-plex cytokine panel with LLOD values for:

Pro-inflammatory: IL-1A (0.09), IL-1B (0.05), IL-6 (0.06), TNF-α (0.11)
Chemokines: IL-8 (1.12), IP-10 (0.31), MCP-1 (0.12), MIP-1A (1.00)
Th1/Th2: IFN-γ (0.37), IL-2 (0.09), IL-4 (0.02), IL-5 (0.14), IL-13 (0.11)
Growth factors: VEGF (0.08), GM-CSF (0.16)
(Full table in config.R)

Color Palettes
Exposures: Eucalyptus (darkorchid), Peat (gold2), Pine (deepskyblue3), Red Oak (brown3), PBS (gray40)
Sex: Male (lemonchiffon), Female (thistle1)
DEG Status: UP (red #bb0c00), DOWN (blue #00AFBB), NO (grey)
Highlighted Genes

SPDEF, IL1A, IL1B, SRXN1, CYP1A1, SOX9, PDIA2, HIC1, TXNRD1, PTGER2, GPR68, SFRP2, F2RL3

Priority GSEA Pathways

Reactive oxygen species pathway
Estrogen/androgen response
Xenobiotic metabolism
Oxidative phosphorylation
PI3K/AKT/mTOR signaling
WNT/β-catenin signaling
(Full list in config.R)

📝 Notes
All file paths are relative to project root (uses here::here())
Output files automatically save to output/tables/, output/figures/, and output/reports/
Knit reports generate HTML with table of contents
All results include significance stars (*, **, ***, #, ##, ###, ~) and LLOD filtering
DEG analysis cached in data/processed/dream_*_results.rds (takes 1-2 hours to run)

📊 Data Availability

Raw data files are provided as supplemental tables with the published manuscript. To run the analysis:

Download supplemental data files from the journal website
Place them in data/raw/ directory (will be created on first run)
Run the analysis scripts in analysis/
Required data files:

MSD cytokine panel measurements
IL-8 ELISA measurements
LDH & TEER measurements
Sample metadata
RNA-seq raw count matrix
Processed data and all figures will be generated automatically in the output/ directory.

🐍 Python Configuration (Optional)

If you encounter reticulate Python errors, they can be safely ignored or suppressed by adding to _load_all.R:

```r
Sys.setenv(RETICULATE_PYTHON = "")
options(reticulate.python = NULL)
```r

👤 Author
mjrod-z

📅 Last Updated
January 2025

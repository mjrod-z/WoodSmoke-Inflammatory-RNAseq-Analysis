# Functions for Data Loading, Cleaning, and Transformation

## Function Definitions

### Load Cytokine Data
load_cytokine_data <- function(file_path) {
    data <- read.csv(file_path)
    return(data)
}

### Impute LOD
impute_lod <- function(data, lod_value) {
    data[is.na(data)] <- lod_value
    return(data)
}

### Calculate Log2 Fold Change PBS
calculate_log2fc_pbs <- function(control, treatment) {
    return(log2(treatment / control))
}

### Calculate Fold Change PBS
calculate_fc_pbs <- function(control, treatment) {
    return(treatment / control)
}

### Normalize Counts
normalize_counts <- function(counts) {
    return(counts / rowSums(counts))
}

### Rename for Coding
rename_for_coding <- function(data) {
    colnames(data) <- gsub("_", ".", colnames(data))
    return(data)
}

### Summarize to Wide
summarize_to_wide <- function(data) {
    library(tidyr)
    return(data %>% spread(key = Variable, value = Value))
}

# End of Functions
# ============================================================================
# DATA ANALYSIS FUNCTIONS
# ============================================================================
# Consolidated module for data processing, statistical testing, and analysis
# Covers both general and inflammatory (cytokine) analyses
#
# Sections:
#   1. Data Loading & Cleaning
#   2. Data Processing (LLOD imputation, outlier detection, transformations)
#   3. Statistical Testing (LMER, ART, Nonparametric)
#   4. Fold Change & Normalization
#   5. Summary & Formatting

library(tidyverse)
library(lme4)
library(emmeans)
library(ARTool)
library(rstatix)
library(outliers)

# ============================================================================
# 1. DATA LOADING & CLEANING FUNCTIONS
# ============================================================================

#' Load and Clean Cytokine Data
#'
#' Reads raw CSV and removes blank values
#'
#' @param file_path Path to CSV file
#'
#' @return Data frame with cleaned cytokine data
#'
load_cytokine_data <- function(file_path) {
  data <- read_csv(file_path)
  data[data == ""] <- NA
  data <- na.omit(data)
  return(data)
}

#' Rename Columns for Coding Convenience
#'
#' Standardizes column names by replacing spaces/hyphens with underscores
#'
#' @param data Data frame
#'
#' @return Data frame with cleaned column names
#'
rename_for_coding <- function(data) {
  colnames(data) <- gsub(" ", "_", colnames(data))
  colnames(data) <- gsub("^MJR-", "MJR_", colnames(data))
  colnames(data) <- gsub("Red_Oak_(\\d+)", "RedOak_\\1", colnames(data))
  return(data)
}

#' Clean and Format Gene Names
#'
#' Removes gene ID prefixes, standardizes naming
#'
#' @param gene_name Character vector of gene names
#'
#' @return Cleaned gene names
#'
clean_gene_names <- function(gene_name) {
  gsub("^(HALLMARK_|LOC|LINC)", "", gene_name) %>%
    gsub("_", " ", .)
}

# ============================================================================
# 2. DATA PROCESSING FUNCTIONS
# ============================================================================

#' Impute Values Below LLOD (Lower Limit of Detection)
#'
#' Uses LLOD/√2 method for imputation. Only imputes cytokines with <30% zeros.
#' Accounts for sex-specific zero proportions.
#'
#' @param input_data Data frame with cytokine measurements
#' @param cols Column names to impute (NULL = all numeric)
#' @param cytokine_llod Data frame with Analyte and LLOD columns
#' @param zero_cutoff Max proportion of zeros to allow (default 0.30)
#'
#' @return List with: data (imputed), valid_cytokines (vector), skipped_cytokines (df)
#'
impute_lod <- function(input_data, cols = NULL, cytokine_llod, zero_cutoff = 0.30) {
  
  if (is.null(cols)) {
    cols <- names(input_data)[sapply(input_data, is.numeric)]
  }
  
  data <- input_data
  valid_cytokines <- c()
  skipped <- data.frame(Analyte = character(), Reason = character())
  
  for (col in cols) {
    if (!(col %in% names(data))) next
    
    # Find LLOD value
    llod_value <- cytokine_llod$LLOD[cytokine_llod$Analyte == col]
    if (length(llod_value) == 0) next
    
    # Count zeros/missing per sex
    zero_props <- data %>%
      group_by(SEX) %>%
      summarise(zero_prop = mean(get(col) == 0 | is.na(get(col)), na.rm = TRUE), .groups = "drop")
    
    # Skip if any sex group >30% zeros
    if (!all(zero_props$zero_prop <= zero_cutoff)) {
      skipped <- rbind(skipped, data.frame(
        Analyte = col,
        Reason = paste0("Too many zeros (", round(max(zero_props$zero_prop)*100, 1), "%)")
      ))
      next
    }
    
    # Impute with LLOD/sqrt(2)
    impute_value <- llod_value / sqrt(2)
    data[[col]][data[[col]] == 0 | is.na(data[[col]])] <- impute_value
    valid_cytokines <- c(valid_cytokines, col)
  }
  
  list(data = data, valid_cytokines = valid_cytokines, skipped_cytokines = skipped)
}

#' Detect Outliers Using Grubbs Test
#'
#' @param data Data frame with numeric columns
#' @param alpha Significance level (default 0.9)
#'
#' @return Data frame with outlier detection results
#'
detect_outliers_grubbs <- function(data, alpha = 0.9) {
  
  results <- data.frame()
  
  for (col in names(data)) {
    if (!is.numeric(data[[col]])) next
    
    test <- try(grubbs.test(data[[col]], opposite = FALSE), silent = TRUE)
    if (inherits(test, "try-error")) next
    
    results <- rbind(results, data.frame(
      Variable = col,
      p_value = test$p.value,
      Outlier = ifelse(test$p.value < (1-alpha), "Yes", "No"),
      stringsAsFactors = FALSE
    ))
  }
  
  return(results)
}

#' Test Normality Using Shapiro-Wilk Test
#'
#' @param data Numeric vector or data frame column
#'
#' @return List with p.value and test result
#'
test_normality <- function(data) {
  data <- as.numeric(data)
  data <- data[!is.na(data)]
  
  if (length(data) < 3) {
    cat("Warning: Sample size (n=", length(data), ") too small for Shapiro-Wilk test\n")
    return(list(p.value = NA, statistic = NA, method = "Too few samples"))
  }
  
  if (length(data) > 5000) {
    cat("Warning: Sample size (n=", length(data), ") too large, using first 5000\n")
    data <- data[1:5000]
  }
  
  shapiro.test(data)
}

#' Test Normality of Original, Log10, and Sqrt Transformed Data
#'
#' @param data Numeric vector
#' @param p_threshold P-value threshold for normality (default 0.05)
#'
#' @return Data frame with transformation results
#'
test_transformations <- function(data, p_threshold = 0.05) {
  
  data_clean <- as.numeric(data)
  data_clean <- data_clean[!is.na(data_clean)]
  
  # Test original
  if (length(data_clean) >= 3) {
    test_orig <- shapiro.test(data_clean)
    p_orig <- test_orig$p.value
    normal_orig <- p_orig > p_threshold
  } else {
    p_orig <- NA
    normal_orig <- NA
  }
  
  # Test log10
  log_data <- log10(data_clean + 1)
  if (length(log_data) >= 3) {
    test_log <- shapiro.test(log_data)
    p_log <- test_log$p.value
    normal_log <- p_log > p_threshold
  } else {
    p_log <- NA
    normal_log <- NA
  }
  
  # Test sqrt
  sqrt_data <- sqrt(data_clean + 1)
  if (length(sqrt_data) >= 3) {
    test_sqrt <- shapiro.test(sqrt_data)
    p_sqrt <- test_sqrt$p.value
    normal_sqrt <- p_sqrt > p_threshold
  } else {
    p_sqrt <- NA
    normal_sqrt <- NA
  }
  
  data.frame(
    Original_p = p_orig,
    Original_Normal = normal_orig,
    Log10_p = p_log,
    Log10_Normal = normal_log,
    Sqrt_p = p_sqrt,
    Sqrt_Normal = normal_sqrt
  )
}

#' Apply Log10 Transformation with Shift
#'
#' @param input_data Data frame or vector
#' @param columns_to_transform Column names to transform
#'
#' @return Transformed data
#'
log10_transform <- function(input_data, columns_to_transform) {
  transformed_data <- as.data.frame(input_data)
  for (col in columns_to_transform) {
    min_value <- min(transformed_data[[col]], na.rm = TRUE)
    shift_value <- ifelse(min_value <= 0, abs(min_value) + 1, 0)
    transformed_data[[col]] <- log10(transformed_data[[col]] + shift_value)
  }
  transformed_data
}

#' Apply Square Root Transformation
#'
#' @param input_data Data frame or vector
#' @param columns_to_transform Column names to transform
#'
#' @return Transformed data
#'
sqrt_transform <- function(input_data, columns_to_transform) {
  transformed_data <- as.data.frame(input_data)
  for (col in columns_to_transform) {
    transformed_data[[col]] <- sqrt(transformed_data[[col]] + 1)
  }
  transformed_data
}

# ============================================================================
# 3. STATISTICAL TESTING FUNCTIONS
# ============================================================================

#' Aligned Rank Transform ANOVA with Pairwise Contrasts
#'
#' @param data Data frame with EXPOSURE, SEX, PATIENTCODE and response variables
#' @param group Filter by group ("All", "M", "F")
#' @param adjust_method P-value adjustment method (default "fdr")
#' @param response_columns Columns to test
#' @param ctrl_level Control level name (default "PBS Control")
#'
#' @return Data frame with results and group column
#'
exposure_art_pairwise <- function(data, group = "All", adjust_method = "fdr",
                                  response_columns = NULL, ctrl_level = "PBS Control") {
  
  filtered_data <- switch(group,
                          "M" = data %>% filter(EXPOSURE != "Untreated Control", SEX == "M"),
                          "F" = data %>% filter(EXPOSURE != "Untreated Control", SEX == "F"),
                          "All" = data %>% filter(EXPOSURE != "Untreated Control"),
                          stop("Invalid group")
  )
  
  filtered_data <- filtered_data %>%
    mutate(
      SEX = factor(SEX),
      EXPOSURE = factor(EXPOSURE),
      PATIENTCODE = factor(PATIENTCODE)
    )
  
  if (is.null(response_columns)) {
    response_columns <- intersect(names(filtered_data), valid_cytokines)
  }
  
  results_list <- list()
  
  for (response in response_columns) {
    if (!response %in% names(filtered_data)) next
    
    formula <- if (group == "All") {
      as.formula(paste(response, "~ EXPOSURE * SEX + Error(PATIENTCODE)"))
    } else {
      as.formula(paste(response, "~ EXPOSURE + Error(PATIENTCODE)"))
    }
    
    m.art <- try(art(formula, data = filtered_data), silent = TRUE)
    if (inherits(m.art, "try-error")) {
      warning(paste("ART model failed for", response))
      next
    }
    
    anova_results <- anova(m.art)
    exposure_row <- anova_results[grepl("^EXPOSURE$", anova_results[[1]], ignore.case = TRUE), ]
    exposure_p <- if (nrow(exposure_row) > 0) exposure_row[["Pr(>F)"]][1] else NA
    
    interaction_p <- if (group == "All") {
      inter_row <- anova_results[grepl("EXPOSURE:SEX", anova_results[[1]], ignore.case = TRUE), ]
      if (nrow(inter_row) > 0) inter_row[["Pr(>F)"]][1] else NA
    } else NA
    
    m.art.con <- artlm.con(m.art, "EXPOSURE")
    emm <- emmeans(m.art.con, ~EXPOSURE)
    pairwise <- contrast(emm, method = "trt.vs.ctrl", ref = ctrl_level, adjust = adjust_method)
    
    results <- as.data.frame(pairwise)
    results$response <- response
    results$group <- group
    results$exposure_p <- exposure_p
    results$interaction_p <- interaction_p
    
    results_list[[response]] <- results
  }
  
  bind_rows(results_list)
}

#' Linear Mixed-Effects Model (LMER) with Pairwise Contrasts
#'
#' @param data Data frame with EXPOSURE, SEX, PATIENTCODE and response variables
#' @param group Filter by group ("All", "M", "F")
#' @param adjust_method P-value adjustment method (default "fdr")
#' @param response_columns Columns to test
#' @param ctrl_level Control level name (default "PBS Control")
#'
#' @return Data frame with results and group column
#'
exposure_lmer_pairwise <- function(data, group = "All", adjust_method = "fdr",
                                   response_columns = NULL, ctrl_level = "PBS Control") {
  
  if (group != "All" && "SEX" %in% names(data)) {
    data <- data %>% filter(SEX == group)
  }
  
  stopifnot("EXPOSURE" %in% names(data), "PATIENTCODE" %in% names(data))
  
  data <- data %>%
    mutate(EXPOSURE = factor(EXPOSURE), PATIENTCODE = factor(PATIENTCODE))
  
  if (ctrl_level %in% levels(data$EXPOSURE)) {
    data$EXPOSURE <- relevel(data$EXPOSURE, ref = ctrl_level)
  }
  
  if (is.null(response_columns)) {
    response_columns <- intersect(names(data), valid_cytokines)
  }
  
  results_list <- list()
  
  for (resp in response_columns) {
    if (!resp %in% names(data)) next
    df <- data %>% filter(!is.na(.data[[resp]]))
    if (length(unique(df$EXPOSURE)) < 2) next
    
    model_formula <- as.formula(paste(resp, "~ EXPOSURE + (1 | PATIENTCODE)"))
    model <- try(lmer(model_formula, data = df), silent = TRUE)
    if (inherits(model, "try-error")) next
    
    emm <- emmeans(model, ~EXPOSURE)
    pairwise <- contrast(emm, method = "trt.vs.ctrl", ref = ctrl_level, adjust = adjust_method)
    pairwise_df <- as.data.frame(summary(pairwise))
    pairwise_df$response <- resp
    pairwise_df$group <- group
    
    results_list[[resp]] <- pairwise_df
  }
  
  bind_rows(results_list)
}

#' Interaction LMER (EXPOSURE × SEX)
#'
#' @param data Data frame with EXPOSURE, SEX, PATIENTCODE and response variables
#' @param group Group filter (usually "All" for interaction)
#' @param adjust_method P-value adjustment method (default "fdr")
#' @param response_columns Columns to test
#' @param ctrl_level Control level name (default "PBS Control")
#'
#' @return Data frame with results for main effects and interaction
#'
interaction_lmer_pairwise <- function(data, group = "All", adjust_method = "fdr",
                                      response_columns = NULL, ctrl_level = "PBS Control") {
  
  if (group != "All" && "SEX" %in% names(data)) {
    data <- data %>% filter(SEX == group)
  }
  
  stopifnot("EXPOSURE" %in% names(data), "PATIENTCODE" %in% names(data))
  
  data <- data %>%
    mutate(EXPOSURE = factor(EXPOSURE), PATIENTCODE = factor(PATIENTCODE))
  
  if (ctrl_level %in% levels(data$EXPOSURE)) {
    data$EXPOSURE <- relevel(data$EXPOSURE, ref = ctrl_level)
  }
  
  if (is.null(response_columns)) {
    response_columns <- intersect(names(data), valid_cytokines)
  }
  
  results_list <- list()
  
  for (resp in response_columns) {
    if (!resp %in% names(data)) next
    df <- data %>% filter(!is.na(.data[[resp]]))
    if (length(unique(df$EXPOSURE)) < 2) next
    
    model_formula <- as.formula(paste(resp, "~ EXPOSURE * SEX + (1 | PATIENTCODE)"))
    model <- try(lmer(model_formula, data = df), silent = TRUE)
    if (inherits(model, "try-error")) next
    
    emm_exposure <- emmeans(model, ~EXPOSURE)
    pairwise_exposure <- contrast(emm_exposure, method = "trt.vs.ctrl", ref = ctrl_level, adjust = adjust_method)
    exposure_df <- as.data.frame(summary(pairwise_exposure))
    exposure_df$type <- "Exposure_vs_Control"
    
    emm_interaction <- emmeans(model, ~SEX | EXPOSURE)
    interaction_contrast <- contrast(emm_interaction, "pairwise", simple = "SEX", combine = TRUE)
    interaction_df <- as.data.frame(summary(interaction_contrast))
    interaction_df$type <- "Exposure*SEX_Interaction"
    
    exposure_df$response <- resp
    interaction_df$response <- resp
    
    results_list[[resp]] <- bind_rows(exposure_df, interaction_df)
  }
  
  bind_rows(results_list)
}

#' Nonparametric Tests (Friedman, Kruskal-Wallis, Wilcoxon)
#'
#' @param data Data frame with EXPOSURE, PATIENTCODE, SEX and response variables
#' @param group Filter by group ("All", "M", "F")
#' @param adjust_method P-value adjustment method (default "BH")
#' @param response_columns Columns to test
#' @param ctrl_level Control level name (default "PBS Control")
#' @param run_friedman Run Friedman paired test (default TRUE)
#' @param run_kruskal Run Kruskal-Wallis unpaired test (default TRUE)
#' @param run_sex Run sex interaction test (default TRUE)
#'
#' @return Data frame with test results
#'
nonparametric_pairwise <- function(data, group = "All", adjust_method = "BH",
                                   response_columns = NULL, ctrl_level = "PBS Control",
                                   run_friedman = TRUE, run_kruskal = TRUE, run_sex = TRUE) {
  
  # Filter data by group
  filtered_data <- switch(group,
                          "M" = data %>% filter(EXPOSURE != "Untreated Control", EXPOSURE != "10X Lysis Buffer", SEX == "M"),
                          "F" = data %>% filter(EXPOSURE != "Untreated Control", EXPOSURE != "10X Lysis Buffer", SEX == "F"),
                          "All" = data %>% filter(EXPOSURE != "Untreated Control", EXPOSURE != "10X Lysis Buffer"),
                          stop("Invalid group. Choose 'All', 'M', or 'F'.")
  )
  
  # Ensure factors
  if ("SEX" %in% names(filtered_data)) {
    filtered_data <- filtered_data %>%
      mutate(
        SEX = factor(SEX),
        EXPOSURE = factor(EXPOSURE),
        PATIENTCODE = factor(PATIENTCODE)
      )
  }
  
  # Determine response columns
  if (is.null(response_columns)) {
    response_columns <- setdiff(
      names(filtered_data)[sapply(filtered_data, is.numeric)],
      grep("_PBS$", names(filtered_data), value = TRUE)
    )
  } else if (is.numeric(response_columns)) {
    response_columns <- names(filtered_data)[response_columns]
  }
  
  all_results <- list()
  
  for (response in response_columns) {
    res_row <- data.frame(response = response, stringsAsFactors = FALSE)
    
    # =========================
    # 1. FRIEDMAN TEST (paired)
    # =========================
    if (run_friedman) {
      friedman_formula <- as.formula(paste(response, "~ EXPOSURE | PATIENTCODE"))
      friedman_p <- NA_real_
      tryCatch({
        friedman_result <- friedman.test(friedman_formula, data = filtered_data)
        friedman_p <- friedman_result$p.value
      }, error = function(e) {
        warning(paste("Friedman test failed for", response, ":", e$message))
      })
      res_row$friedman_p <- friedman_p
      
      # Pairwise Wilcoxon (paired) vs PBS Control
      exposure_groups <- setdiff(levels(filtered_data$EXPOSURE), ctrl_level)
      reference <- filtered_data %>% filter(EXPOSURE == ctrl_level) %>% arrange(PATIENTCODE)
      
      for (exp_group in exposure_groups) {
        exp_data <- filtered_data %>% filter(EXPOSURE == exp_group) %>% arrange(PATIENTCODE)
        if (nrow(exp_data) == nrow(reference) &&
            all(exp_data$PATIENTCODE == reference$PATIENTCODE)) {
          # Check for sufficient non-missing observations
          x_vals <- exp_data[[response]]
          ref_vals <- reference[[response]]
          
          if (sum(!is.na(x_vals)) >= 1 && sum(!is.na(ref_vals)) >= 1) {
            wilcox_p <- suppressWarnings(
              wilcox.test(x_vals, ref_vals,
                          paired = TRUE, exact = FALSE)$p.value
            )
            res_row[[paste0("friedman_wilcox_", exp_group, "_p")]] <- wilcox_p
          }
        }
      }
    }
    
    # =========================
    # 2. KRUSKAL-WALLIS TEST (independent)
    # =========================
    if (run_kruskal) {
      kw_formula <- as.formula(paste(response, "~ EXPOSURE"))
      kw_p <- NA_real_
      tryCatch({
        # Check if we have multiple groups
        n_groups <- length(unique(filtered_data$EXPOSURE[!is.na(filtered_data[[response]])]))
        if (n_groups < 2) {
          warning(paste("Kruskal-Wallis test skipped for", response, "- only", n_groups, "group(s)"))
          kw_p <- NA_real_
        } else {
          kw_result <- kruskal.test(kw_formula, data = filtered_data)
          kw_p <- kw_result$p.value
        }
      }, error = function(e) {
        warning(paste("Kruskal-Wallis test failed for", response, ":", e$message))
      })
      res_row$kruskal_p <- kw_p
      
      # Pairwise Wilcoxon (unpaired) vs PBS Control
      exposure_groups <- setdiff(levels(filtered_data$EXPOSURE), ctrl_level)
      reference <- filtered_data %>% filter(EXPOSURE == ctrl_level)
      
      for (exp_group in exposure_groups) {
        exp_data <- filtered_data %>% filter(EXPOSURE == exp_group)
        
        # Check for sufficient non-missing observations
        x_vals <- exp_data[[response]]
        ref_vals <- reference[[response]]
        
        if (sum(!is.na(x_vals)) >= 1 && sum(!is.na(ref_vals)) >= 1) {
          tryCatch({
            wilcox_p <- suppressWarnings(
              wilcox.test(x_vals, ref_vals,
                          paired = FALSE, exact = FALSE)$p.value
            )
            res_row[[paste0("kruskal_wilcox_", exp_group, "_p")]] <- wilcox_p
          }, error = function(e) {
            # Silently skip if Wilcoxon fails
          })
        }
      }
    }
    
    # =========================
    # 3. SEX TEST (only if group == "All")
    # =========================
    if (run_sex && group == "All" && "SEX" %in% names(filtered_data)) {
      sex_p <- NA_real_
      tryCatch({
        sex_data <- filtered_data %>%
          mutate(SEX_EXPOSURE = interaction(SEX, EXPOSURE, sep = "_"))
        kw_sex_formula <- as.formula(paste(response, "~ SEX_EXPOSURE"))
        kw_sex_result <- kruskal.test(kw_sex_formula, data = sex_data)
        sex_p <- kw_sex_result$p.value
      }, error = function(e) {
        warning(paste("Sex comparison failed for", response, ":", e$message))
      })
      res_row$sex_p <- sex_p
    }
    
    all_results[[response]] <- res_row
  }
  
  # Combine results + FDR correction
  combined_results <- bind_rows(all_results)
  p_cols <- grep("_p$", names(combined_results), value = TRUE)
  for (pcol in p_cols) {
    padj_col <- sub("_p$", "_adj", pcol)
    combined_results[[padj_col]] <- p.adjust(combined_results[[pcol]], method = adjust_method)
  }
  
  return(combined_results)
}

# ============================================================================
# 4. FOLD CHANGE & NORMALIZATION FUNCTIONS
# ============================================================================

#' Calculate Fold Change Relative to PBS Control (PBS Normalization)
#'
#' @param data Data frame with PATIENTCODE, GROUP, and value columns
#' @param group_col Column containing group information
#' @param value_cols Column names to calculate fold change for
#'
#' @return Data frame with new columns named `{col}_PBS`
#'
calculate_fold_change_PBS <- function(data, group_col, value_cols) {
  data %>%
    group_by(PATIENTCODE) %>%
    mutate(across({{ value_cols }},
                  ~ (. - .[GROUP == "Control"]) / .[GROUP == "Control"],
                  .names = "{.col}_PBS"
    ))
}

#' Calculate Sample Fold Change (Sample - Control) / Control
#'
#' @param data Data frame with EXPOSURE, PATIENTCODE, and response columns
#' @param response_col Column name with measurements
#' @param control_name Name of control group (default "PBS Control")
#'
#' @return Data frame with FoldChange column added
#'
calculate_sample_fold_change <- function(data, response_col, control_name = "PBS Control") {
  
  data <- data %>%
    group_by(PATIENTCODE) %>%
    mutate(
      control_value = mean(!!sym(response_col)[EXPOSURE == control_name], na.rm = TRUE),
      FoldChange = (!!sym(response_col) - control_value) / control_value
    ) %>%
    ungroup() %>%
    select(-control_value)
  
  return(data)
}

#' Summarize Fold Change by Exposure
#'
#' @param data Data frame with FoldChange column
#' @param exclude_control Exclude control from summary (default TRUE)
#'
#' @return Data frame with mean FC, SD, N per exposure
#'
summarize_fold_change <- function(data, exclude_control = TRUE) {
  
  if (exclude_control) {
    data <- data %>% filter(EXPOSURE != "PBS Control")
  }
  
  summary <- data %>%
    group_by(EXPOSURE) %>%
    summarise(
      N = n(),
      MeanFC = mean(FoldChange, na.rm = TRUE),
      SD = sd(FoldChange, na.rm = TRUE),
      SE = SD / sqrt(N),
      Median = median(FoldChange, na.rm = TRUE),
      Min = min(FoldChange, na.rm = TRUE),
      Max = max(FoldChange, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(MeanFC))
  
  return(summary)
}

# ============================================================================
# 5. SUMMARY & FORMATTING FUNCTIONS
# ============================================================================

#' Summarize to Wide Format (Mean ± SD by Exposure)
#'
#' Creates summary table with mean ± SD by exposure
#'
#' @param data Data frame with EXPOSURE and measurement columns
#' @param measure_vars Column names to summarize
#'
#' @return Wide format data frame with mean ± SD
#'
summarize_to_wide <- function(data, measure_vars) {
  
  summary_list <- list()
  
  for (var in measure_vars) {
    summary <- data %>%
      group_by(EXPOSURE) %>%
      summarise(
        !!paste0(var, "_mean") := mean(!!sym(var), na.rm = TRUE),
        !!paste0(var, "_sd") := sd(!!sym(var), na.rm = TRUE),
        .groups = "drop"
      )
    
    if (is.null(summary_list[[1]])) {
      summary_list[[1]] <- summary
    } else {
      summary_list[[1]] <- left_join(summary_list[[1]], summary, by = "EXPOSURE")
    }
  }
  
  return(summary_list[[1]])
}

#' Add Significance Stars to Results
#'
#' Adds star annotations based on p-value and estimate direction
#'
#' @param data Data frame with p.value and estimate columns
#'
#' @return Data frame with new 'stars' column
#'
add_significance_stars <- function(data) {
  data %>%
    mutate(stars = case_when(
      p.value < 0.001 & estimate > 0  ~ "***",
      p.value < 0.001 & estimate < 0  ~ "###",
      p.value < 0.01  & estimate > 0  ~ "**",
      p.value < 0.01  & estimate < 0  ~ "##",
      p.value < 0.05  & estimate > 0  ~ "*",
      p.value < 0.05  & estimate < 0  ~ "#",
      p.value < 0.1   & estimate != 0 ~ "~",
      TRUE ~ ""
    ))
}

cat("✓ Data analysis functions loaded\n")

### 4 #### Metrics Calculation Script ####

# This script calculates various metrics to evaluate the performance of differential expression analysis methods.
# It includes precision, recall, F1-score, accuracy, False Discovery Rate (FDR), Area Under the Receiver Operating
# Characteristic Curve (AUROC), and Precision-Recall Area Under the Curve (PRAUC).

# The script includes the following sections:

# 1. Metric Calculations:
#    - Calculate precision, recall, F1-score, accuracy, and FDR for each method.

# 2. AUROC Calculations:
#    - Calculate the AUROC for each method.

# 3. PRAUC Calculations:
#    - Calculate the PRAUC for each method.

# This script helps evaluate the effectiveness of different statistical methods for identifying differentially expressed genes (DEGs).

### Section 4: Metrics ######

# Load necessary libraries
library(caret)

# Function to calculate precision, recall, F1-score, accuracy, and FDR
calculate_metrics <- function(expected, found, all_genes) {
  tp <- length(intersect(expected, found)) # True Positives
  fp <- length(setdiff(found, expected))  # False Positives
  fn <- length(setdiff(expected, found)) # False Negatives
  tn <- length(setdiff(setdiff(all_genes, expected), found)) # True Negatives
  
  precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
  f1_score <- ifelse(precision + recall == 0, 0, 2 * (precision * recall) / (precision + recall))
  accuracy <- ifelse(tp + fp + fn + tn == 0, 0, (tp + tn) / (tp + fp + fn + tn))
  fdr <- ifelse(tp + fp == 0, 0, fp / (tp + fp))
  
  return(list(precision = precision, recall = recall, f1_score = f1_score, accuracy = accuracy, fdr = fdr))
}

# Lists to store results
metrics_findallmarkers <- list()

# All genes
all_genes <- rownames(sim.groups)

# Comparison for DEFacGroup > 1.5

# Expected DEFacGroup
expected <- summary_over_df_filtered_1$Gene

# Calculate metrics for different test methods
found_wilcox <- summary_over_df_wilcox$Gene
metrics_findallmarkers[["Wilcoxon"]] <- calculate_metrics(expected, found_wilcox, all_genes)

found_t_test <- summary_over_df_t_test$Gene
metrics_findallmarkers[["T-Test"]] <- calculate_metrics(expected, found_t_test, all_genes)

found_deseq2 <- summary_over_df_deseq2$Gene
metrics_findallmarkers[["DESeq2"]] <- calculate_metrics(expected, found_deseq2, all_genes)

found_mast <- summary_over_df_mast$Gene
metrics_findallmarkers[["MAST"]] <- calculate_metrics(expected, found_mast, all_genes)

# Convert lists to dataframes for easy visualization
metrics_findallmarkers_df <- do.call(rbind, lapply(metrics_findallmarkers, function(x) unlist(x)))
metrics_findallmarkers_df <- data.frame(Test = names(metrics_findallmarkers), metrics_findallmarkers_df)

# Display metrics for FindAllMarkers with different test methods
print("Metrics for FindAllMarkers with different test methods:")
print(metrics_findallmarkers_df)

### Section 5: AUROC ######

library(pROC)
library(PRROC)
library(dplyr)

# Function to calculate AUROC
calculate_auroc <- function(expected, all_genes, scores) {
  # Create a binary vector for expected genes
  labels <- as.numeric(all_genes %in% expected)
  
  # Filter NA and infinite values
  valid_indices <- !is.na(scores) & is.finite(scores)
  labels <- labels[valid_indices]
  scores <- scores[valid_indices]
  
  # Calculate AUROC
  roc_obj <- roc(labels, scores)
  auroc <- auc(roc_obj)
  
  return(auroc)
}

all_genes <- rownames(sim.groups)

# Comparison for DEFacGroup > 1.5
# Expected DEFacGroup
expected <- summary_over_df_filtered_1$Gene

# Calculate AUROC for Wilcoxon
scores_wilcox <- -log10(comparison_result_wilcox$p_val_adj)
scores_wilcox <- scores_wilcox[match(all_genes, rownames(comparison_result_wilcox))]
auroc_wilcox <- calculate_auroc(expected, all_genes, scores_wilcox)

# Calculate AUROC for T-Test
scores_t_test <- -log10(comparison_result_t_test$p_val_adj)
scores_t_test <- scores_t_test[match(all_genes, rownames(comparison_result_t_test))]
auroc_t_test <- calculate_auroc(expected, all_genes, scores_t_test)

# Calculate AUROC for DESeq2
scores_deseq2 <- -log10(comparison_result_deseq2$p_val_adj)
scores_deseq2 <- scores_deseq2[match(all_genes, rownames(comparison_result_deseq2))]
auroc_deseq2 <- calculate_auroc(expected, all_genes, scores_deseq2)

# Calculate AUROC for MAST
scores_mast <- -log10(comparison_result_mast$p_val_adj)
scores_mast <- scores_mast[match(all_genes, rownames(comparison_result_mast))]
auroc_mast <- calculate_auroc(expected, all_genes, scores_mast)

# AUROC for each method
print("AUROC for Wilcoxon:")
print(auroc_wilcox)

print("AUROC for T-Test:")
print(auroc_t_test)

print("AUROC for DESeq2:")
print(auroc_deseq2)

print("AUROC for MAST:")
print(auroc_mast)

### Section 6: PRAUC ######

# Function to calculate PRAUC
calculate_prauc <- function(expected, all_genes, scores) {
  # Create a binary vector for expected genes
  labels <- as.numeric(all_genes %in% expected)
  
  # Filter NA and infinite values
  valid_indices <- !is.na(scores) & is.finite(scores)
  labels <- labels[valid_indices]
  scores <- scores[valid_indices]
  
  # Calculate PRAUC
  pr_obj <- pr.curve(scores.class0 = scores[labels == 0],
                     scores.class1 = scores[labels == 1],
                     curve = TRUE)
  prauc <- pr_obj$auc.integral
  
  return(prauc)
}

# Comparison for DEFacGroup > 1.5
# Expected DEFacGroup
expected <- summary_over_df_filtered_1$Gene

# Calculate PRAUC for Wilcoxon
scores_wilcox <- -log10(comparison_result_wilcox$p_val_adj)
scores_wilcox <- scores_wilcox[match(all_genes, rownames(comparison_result_wilcox))]
prauc_wilcox <- calculate_prauc(expected, all_genes, scores_wilcox)

# Calculate PRAUC for T-Test
scores_t_test <- -log10(comparison_result_t_test$p_val_adj)
scores_t_test <- scores_t_test[match(all_genes, rownames(comparison_result_t_test))]
prauc_t_test <- calculate_prauc(expected, all_genes, scores_t_test)

# Calculate PRAUC for DESeq2
scores_deseq2 <- -log10(comparison_result_deseq2$p_val_adj)
scores_deseq2 <- scores_deseq2[match(all_genes, rownames(comparison_result_deseq2))]
prauc_deseq2 <- calculate_prauc(expected, all_genes, scores_deseq2)

# Calculate PRAUC for MAST
scores_mast <- -log10(comparison_result_mast$p_val_adj)
scores_mast <- scores_mast[match(all_genes, rownames(comparison_result_mast))]
prauc_mast <- calculate_prauc(expected, all_genes, scores_mast)

# PRAUC for each method
print("PRAUC for Wilcoxon:")
print(prauc_wilcox)

print("PRAUC for T-Test:")
print(prauc_t_test)

print("PRAUC for DESeq2:")
print(prauc_deseq2)

print("PRAUC for MAST:")
print(prauc_mast)

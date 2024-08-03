## 3 ###### Differential Expression Analysis Script ####
# This script performs differential expression analysis using various statistical methods
# to identify differentially expressed genes (DEGs) across different clusters of single-cell RNA sequencing data.

# Methods used in this analysis:
# 1. Wilcoxon Rank Sum Test (Wilcoxon):
#    - A non-parametric test that compares two independent groups to determine if they come from the same distribution.
#    - Suitable for small sample sizes and does not assume normality.
#    - Used to identify DEGs by comparing gene expression levels between clusters.

# 2. T-Test:
#    - A parametric test that compares the means of two independent groups.
#    - Assumes that the data is normally distributed.
#    - Used to identify DEGs by comparing the average expression levels of genes between clusters.

# 3. DESeq2:
#    - A method based on the negative binomial distribution, designed for RNA-Seq data.
#    - Accounts for differences in sequencing depth and variability.
#    - Identifies DEGs by comparing gene expression counts between clusters and provides statistical significance.

# 4. MAST (Model-based Analysis of Single-cell Transcriptomics):
#    - A method specifically designed for single-cell RNA-Seq data.
#    - Models gene expression as a hurdle model with two parts: a logistic regression for the detection probability and a Gaussian linear model for the expression levels.
#    - Identifies DEGs by accounting for both the presence/absence and the level of expression of genes in single-cell data.

# The script proceeds by performing these analyses on a Seurat object containing single-cell data,
# and generates summary tables and visualizations for the identified DEGs.

### Section 3: Differential Expression Analysis Methods ######

# Install and update necessary packages
if (!requireNamespace("MAST", quietly = TRUE)) {
  BiocManager::install("MAST")
}

# Update all installed packages
update.packages(ask = FALSE)

# Reinstall MAST and its dependencies
BiocManager::install("MAST", type = "source")

# Function to perform DEG analysis using different test methods
perform_comparison_all <- function(data, test_use = "wilcox", only_pos = TRUE, min_pct = 0, logfc_threshold = 0) {
  comparison_result <- FindAllMarkers(data, 
                                      only.pos = only_pos, 
                                      min.pct = min_pct, 
                                      logfc.threshold = logfc_threshold, 
                                      test.use = test_use)
  comparison_result$p_val_adj <- p.adjust(comparison_result$p_val, method = "fdr")
  return(comparison_result)
}

# Perform comparison for DEFacGroup > 1.5 (only overexpressed) using different test methods
start_time <- Sys.time()
comparison_result_wilcox <- perform_comparison_all(data_1_cl, test_use = "wilcox", only_pos = TRUE, logfc_threshold = log2(1.5))
comparison_result_t_test <- perform_comparison_all(data_1_cl, test_use = "t", only_pos = TRUE, logfc_threshold = log2(1.5))
comparison_result_deseq2 <- perform_comparison_all(data_1_cl, test_use = "DESeq2", only_pos = TRUE, logfc_threshold = log2(1.5))
comparison_result_mast <- perform_comparison_all(data_1_cl, test_use = "MAST", only_pos = TRUE, logfc_threshold = log2(1.5))
end_time <- Sys.time()

# Summary of the number of DEGs per cluster for each case DEFacGroup > 1.5
deg_summary_wilcox <- comparison_result_wilcox %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  summarise(Num_DEGs = n())

deg_summary_t_test <- comparison_result_t_test %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  summarise(Num_DEGs = n())

deg_summary_deseq2 <- comparison_result_deseq2 %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  summarise(Num_DEGs = n())

deg_summary_mast <- comparison_result_mast %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  summarise(Num_DEGs = n())

# Include clusters with no significant DEGs, with a value of 0
all_clusters <- sort(unique(Idents(data_1_cl)))
deg_summary_wilcox <- merge(data.frame(Cluster = all_clusters), deg_summary_wilcox, by.x = "Cluster", by.y = "cluster", all.x = TRUE)
deg_summary_wilcox$Num_DEGs[is.na(deg_summary_wilcox$Num_DEGs)] <- 0

deg_summary_t_test <- merge(data.frame(Cluster = all_clusters), deg_summary_t_test, by.x = "Cluster", by.y = "cluster", all.x = TRUE)
deg_summary_t_test$Num_DEGs[is.na(deg_summary_t_test$Num_DEGs)] <- 0

deg_summary_deseq2 <- merge(data.frame(Cluster = all_clusters), deg_summary_deseq2, by.x = "Cluster", by.y = "cluster", all.x = TRUE)
deg_summary_deseq2$Num_DEGs[is.na(deg_summary_deseq2$Num_DEGs)] <- 0

deg_summary_mast <- merge(data.frame(Cluster = all_clusters), deg_summary_mast, by.x = "Cluster", by.y = "cluster", all.x = TRUE)
deg_summary_mast$Num_DEGs[is.na(deg_summary_mast$Num_DEGs)] <- 0

# Create summary tables for FindAllMarkers
summary_over_df_wilcox <- comparison_result_wilcox %>%
  filter(p_val_adj < 0.05) %>%
  mutate(Gene = rownames(comparison_result_wilcox)) %>%
  select(Gene, cluster, avg_log2FC) %>%
  rename(DEFacGroup = avg_log2FC)

summary_over_df_t_test <- comparison_result_t_test %>%
  filter(p_val_adj < 0.05) %>%
  mutate(Gene = rownames(comparison_result_t_test)) %>%
  select(Gene, cluster, avg_log2FC) %>%
  rename(DEFacGroup = avg_log2FC)

summary_over_df_deseq2 <- comparison_result_deseq2 %>%
  filter(p_val_adj < 0.05) %>%
  mutate(Gene = rownames(comparison_result_deseq2)) %>%
  select(Gene, cluster, avg_log2FC) %>%
  rename(DEFacGroup = avg_log2FC)

summary_over_df_mast <- comparison_result_mast %>%
  filter(p_val_adj < 0.05) %>%
  mutate(Gene = rownames(comparison_result_mast)) %>%
  select(Gene, cluster, avg_log2FC) %>%
  rename(DEFacGroup = avg_log2FC)

# Sort and display the DEG summaries by comparison
deg_summary_wilcox <- deg_summary_wilcox[order(deg_summary_wilcox$Cluster), ]
deg_summary_t_test <- deg_summary_t_test[order(deg_summary_t_test$Cluster), ]
deg_summary_deseq2 <- deg_summary_deseq2[order(deg_summary_deseq2$Cluster), ]
deg_summary_mast <- deg_summary_mast[order(deg_summary_mast$Cluster), ]

print(deg_summary_wilcox) # FindAllMarkers - Wilcoxon
print(deg_summary_t_test) # FindAllMarkers - T-Test
print(deg_summary_deseq2) # FindAllMarkers - DESeq2
print(deg_summary_mast) # FindAllMarkers - MAST

# Display DEFacGroup results by cluster
print(deg_counts_df_filtered_1) # DEFacGroup Foldchange > 1.5

# Display summary tables for FindAllMarkers
print(head(summary_over_df_wilcox))
print(dim(summary_over_df_wilcox))

print(head(summary_over_df_t_test))
print(dim(summary_over_df_t_test))

print(head(summary_over_df_deseq2))
print(dim(summary_over_df_deseq2))

print(head(summary_over_df_mast))
print(dim(summary_over_df_mast))


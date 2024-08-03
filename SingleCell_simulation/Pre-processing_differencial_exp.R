### 2 ####  Pre-processing and Simulation Comparison Script ####
# This script performs pre-processing on single-cell RNA sequencing (scRNA-seq) data and compares simulated data
# with real data to evaluate the quality and characteristics of the simulation.

# The script includes the following sections:

# 1. Pre-processing:
#    - Loading the count matrix from the simulated data.
#    - Creating a Seurat object for data handling and analysis.
#    - Normalizing the data using log normalization.
#    - Identifying highly variable features (genes) to focus on the most informative genes.
#    - Scaling the data to center the gene expression values.
#    - Performing Principal Component Analysis (PCA) to reduce dimensionality and capture the major sources of variation.
#    - Finding neighbors based on PCA results to construct a Shared Nearest Neighbor (SNN) graph.
#    - Clustering the cells into distinct groups based on the SNN graph.
#    - Performing Uniform Manifold Approximation and Projection (UMAP) for visualization of clusters in lower dimensions.

# 2. Comparison with Real Data:
#    - Comparing the simulated data with real scRNA-seq data to assess the simulation's accuracy.
#    - Generating comparison plots to visualize the mean expression levels between simulated and real data.
#    - Providing insights into the similarities and differences between simulated and real datasets.

# This script ensures that the simulated data closely resembles the real data, aiding in the development
# and validation of computational models for scRNA-seq analysis.


### Section 1: DEFacGroup Analysis ######

# Obtain differential expression factor data for each group
de_fac_groups <- lapply(1:6, function(i) {
  sim.groups@rowRanges@elementMetadata[[paste0("DEFacGroup", i)]]
})

# Calculate minimum and maximum values for each group
min_max_values <- data.frame(
  Group = paste0("DEFacGroup", 1:6),
  Min = sapply(de_fac_groups, min),
  Max = sapply(de_fac_groups, max)
)

# Display the results
print(min_max_values)

# Define expected DEGs for each cluster
expected_degs_filtered <- list()
for (i in 1:6) {
  cluster_name <- paste0("Cluster_", i-1)
  expected_degs_filtered[[cluster_name]] <- rownames(sim.groups)[rowData(sim.groups)[[paste0("DEFacGroup", i)]] > 1.5]
}

# Create a table with the number of expected DEGs in each cluster
deg_counts_filtered <- sapply(expected_degs_filtered, length)
deg_counts_df_filtered_1 <- data.frame(
  Cluster = names(deg_counts_filtered),
  Num_DEGs_Overexpressed = deg_counts_filtered
)

# Display the table with the number of expected DEGs
print(deg_counts_df_filtered_1) # DEFacGroup Foldchange > 1.5

# Function to create summary tables
create_summary <- function(expected_degs) {
  gene_names <- c()
  clusters <- c()
  defac_groups <- c()
  
  for (i in 1:6) {
    cluster_name <- paste0("Cluster_", i-1)
    group_name <- paste0("DEFacGroup", i)
    
    genes <- expected_degs[[cluster_name]]
    
    gene_names <- c(gene_names, genes)
    clusters <- c(clusters, rep(cluster_name, length(genes)))
    defac_groups <- c(defac_groups, rowData(sim.groups)[genes, group_name])
  }
  
  data.frame(
    Gene = gene_names,
    Cluster = clusters,
    DEFacGroup = defac_groups
  )
}

# Create summary table for DEFacGroup > 1.5
summary_over_df_filtered_1 <- create_summary(expected_degs_filtered)

# Display the summary tables
print(head(summary_over_df_filtered_1)) # DEFacGroup Foldchange > 1.5
print(dim(summary_over_df_filtered_1)) # DEFacGroup Foldchange > 1.5

# Find all differentially expressed genes (DEGs) between clusters
all_DEGs <- FindAllMarkers(data_0_cl, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1.5)

# Display the first results
head(all_DEGs)
tail(all_DEGs)

# Count the number of DEGs per cluster
DEGs_by_cluster <- all_DEGs %>%
  group_by(cluster) %>%
  summarize(num_DEGs = n())

# Display the number of DEGs per cluster
print(DEGs_by_cluster)

# Ensure all clusters are represented, even if they have no significant DEGs
all_clusters <- sort(unique(Idents(data_0_cl)))
deg_summary_mastREF <- merge(data.frame(Cluster = all_clusters), deg_summary_mastREF, by.x = "Cluster", by.y = "cluster", all.x = TRUE)
deg_summary_mastREF$Num_DEGs[is.na(deg_summary_mastREF$Num_DEGs)] <- 0

# Sort and display the DEG summary by comparison
deg_summary_mastREF <- deg_summary_mastREF[order(deg_summary_mastREF$Cluster), ]
print(deg_summary_mastREF)
head(comparison_result_mastREF)

### Section 2: Simulation ######

# Pre-processing
matrix_counts <- sim.groups@assays@data@listData[["counts"]] # Not normalized
data_1 <- CreateSeuratObject(counts = matrix_counts)
data_1 <- NormalizeData(data_1, normalization.method = "LogNormalize", scale.factor = 10000)
data_1 <- FindVariableFeatures(data_1, selection.method = "vst", nfeatures = 2000)
data_1 <- ScaleData(data_1)
data_1_pca <- RunPCA(data_1, features = VariableFeatures(object = data_1))
ElbowPlot(data_1_pca)
data_1_SNN <- FindNeighbors(data_1_pca, dims = 1:19) # Adjust the number of dimensions
data_1_cl <- FindClusters(data_1_SNN, resolution = 0.8) # Adjust the resolution
data_1_UMAP <- RunUMAP(data_1_cl, dims = 1:10)
DimPlot(data_1_UMAP, label = TRUE, reduction = "umap", pt.size = 0.5) + NoLegend()

# Display the number of cells per cluster
table(Idents(data_1_cl))
# Verify the distribution of simulated groups
table(colData(sim.groups)$Group)

### Section 3: Comparison with Real Data ####

# Comparison with real data
comparison <- compareSCEs(list(Splat = sim.groups, Real = as.SingleCellExperiment(data_0)))
plot_means_comparison <- comparison$Plots$Means
plot_means_comparison
# Save plot (uncomment to save)
# ggsave(filename = paste0(outdir,"/Mean_expression_simulation_vs_realdata_1.png"), plot = plot_means_comparison)
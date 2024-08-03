
### 1 ### Single-cell RNA-seq Simulation Script ####
# This script simulates single-cell RNA sequencing data using the Splatter package.
# It includes loading reference data, pre-processing steps, and performing PCA and clustering
# to analyze and visualize the simulated data. The script also identifies differentially expressed genes (DEGs)
# and generates simulation data to validate the results.

####### Libraries ####

# Clear all variables from the workspace
rm(list=ls(all=TRUE))

# Set working directory and output directory for plots
setwd("/Users/MRo/Desktop/git_RedMra")
outdir <- "/Users/MRo/Desktop/git_RedMra/plots"
par(mfrow = c(1, 1)) 
plot.new()

# Load necessary libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("splatter") 
BiocManager::install("Seurat") 
install.packages("ggplot2")
install.packages("tidyverse")
BiocManager::install("scater")
BiocManager::install("SeuratData")
devtools::install_github('satijalab/seurat-data')

# Load libraries
library("ggplot2")
library("scater")
library("splatter")
library("Seurat")
library("SeuratData")
library("tidyverse")
library(Matrix)
library(SingleCellExperiment)


## Load Reference Data ####
# Article - (2017) Splatter: simulation of single-cell RNA sequencing data. https://doi.org/10.1186/s13059-017-1305-0
# Tutorial - Splatter: https://bioconductor.statistik.tu-dortmund.de/packages/3.13/bioc/vignettes/splatter/inst/doc/splatter.html

#### 1. Reference Data pbmc3k 

# Load pbmc3k data
AvailableData()
InstallData("pbmc3k") 

pbmc3k <- LoadData("pbmc3k.SeuratData")
pbmc3k <- UpdateSeuratObject(pbmc3k)

# Data inspection
data_0 <- pbmc3k
class(data_0) # SeuratObject

# Frequency table of cell types
tab0 <- round(table(data_0$seurat_annotations) / sum(table(data_0$seurat_annotations)), 2)
tab0

# Extract gene expression counts
data_counts <- as.matrix(data_0@assays$RNA$counts)
class(data_counts) # matrix
dim(data_counts) # 13714 x 2700

# Estimate parameters for simulation
params <- splatEstimate(data_counts)
params
getParams(params, names = c("nGenes", "nCells", "mean.rate", "mean.shape", "nGroups"))

# Data preprocessing
data_0 <- ScaleData(data_0)
data_0 <- FindVariableFeatures(data_0)

# Perform PCA
data_0_PCA <- RunPCA(data_0, features = VariableFeatures(object = data_0))
DimPlot(data_0_PCA, reduction = "pca") + ggtitle("Real Data (PCA)")

# Clustering the cells
data_0_kNN <- FindNeighbors(data_0_PCA, dims = 1:8)
data_0_cl <- FindClusters(data_0_kNN, resolution = 0.5)
DimPlot(data_0_cl, group.by = "seurat_clusters")

# UMAP visualization
data_0_UMAP <- RunUMAP(data_0_cl, dims = 1:10)
DimPlot(data_0_UMAP, reduction = "umap", pt.size = 0.5) + NoLegend()

###### DEGs Reference ###### 

# Identify DEGs
all_DEGs <- FindAllMarkers(data_0_cl, only.pos = TRUE, logfc.threshold = log2(1.5))
head(all_DEGs)
tail(all_DEGs)

# Count DEGs per cluster
DEGs_by_cluster <- all_DEGs %>%
  group_by(cluster) %>%
  summarize(num_DEGs = n())
print(DEGs_by_cluster)

# Identify DEGs between clusters 0 and 1
cluster_DEGs <- FindMarkers(data_0_cl, ident.1 = 0, ident.2 = 1)
head(cluster_DEGs)


###### # Simulation of scRNA data  ########

#### Simulation of scRNA data for ONE individual with different types of cells ####

set.seed(42)

# Simulation data 1
# sim.groups <- splatSimulate(
# params, 
# group.prob = c(0.45, 0.23, 0.15, 0.10, 0.04, 0.03), # Group membership probabilities
# de.prob = c(0.10, 0.09, 0.07, 0.06, 0.04, 0.04),    # Differential expression probabilities for each group
# de.facLoc = c(0.4, 1, 0.8, 0.8, 0.5, 1),            # Location of differential expression factor, adjusted for each group
# de.facScale = c(0.9, 0.9, 0.7, 0.8, 0.8, 0.9),      # Scale of differential expression factor, adjusted for each group
# batchCells = 3000,                                  # Number of cells per batch in the simulation (3000)
# method = "groups",                                  # Simulation method (groups)
# verbose = TRUE
# )

# Reference data 2 

data_counts <- as.matrix(data_0@assays$RNA$counts)

params_improved <- splatEstimate(data_counts)
params_improved <- setParams(params_improved, dropout.type = "none", dropout.mid = -0.15, dropout.shape = -1.06)

sim.groups <- splatSimulate(
  params_improved, 
  group.prob = c(0.45, 0.23, 0.15, 0.10, 0.04, 0.03),
  de.prob = c(0.10, 0.09, 0.07, 0.06, 0.04, 0.04),
  de.facLoc = c(0.4, 1, 0.8, 0.8, 0.5, 1),
  de.facScale = c(0.9, 0.9, 0.7, 0.8, 0.8, 0.9),
  batchCells = 4000,
  nGenes = 15000,
  method = "groups",
  verbose = TRUE
)

# Save simulation results
saveRDS(sim.groups, file = "simulation_F1.rds")

# Explore simulation results
str(sim.groups)
print(sim.groups)
table(colData(sim.groups)$Group)








### Workshop Script 5 MR Notes ####
# This script processes gene expression data from RNA-Seq experiments.
# It includes loading data, filtering differentially expressed genes (DEGs), 
# performing annotation, PCA analysis, and GO enrichment analysis.

# Main sections of the script:

# 1. **Loading Data**:
#    - Reads in gene expression data and DEGs from specified files.
#    - Normalizes the count data to make it suitable for further analysis.

# 2. **Filtering DEGs**:
#    - Filters DEGs based on log2 fold change and p-value.
#    - Separates the data into samples and controls.

# 3. **Annotation**:
#    - Uses Bioconductor annotation packages to annotate genes with GO terms, pathways, and other relevant information.
#    - Generates summary tables of annotations.

# 4. **PCA Analysis**:
#    - Performs Principal Component Analysis (PCA) to reduce dimensionality and visualize the data.
#    - Clusters genes based on their expression profiles using k-means clustering.
#    - Creates visualizations of the PCA results and the clusters.

# 5. **GO Enrichment Analysis**:
#    - Conducts Gene Ontology (GO) enrichment analysis to identify overrepresented GO terms in the set of DEGs.
#    - Visualizes the results using dot plots.


##### Path #####
rm(list=ls(all=TRUE))
setwd("/Users/")
par(mfrow = c(1, 1))

#### Libraries #####

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db") # Install human genome annotation database

install.packages("ade4") # PCA, Multivariate analysis, Principal Component Analysis
library(ade4)

# Load Bioconductor annotation libraries
options(connectionObserver = NULL)
library(org.Hs.eg.db)

install.packages("corrplot")
install.packages("pheatmap")
library(corrplot)
library(pheatmap)

##### Load Data #####

## Data for DEGs
Degs1 <- read.delim("Esqz_Result_05.tsv", header = TRUE)
head(Degs1)
dim(Degs1)
rownames(Degs1)
# Subset with FC > 0.585
sub_Degs1 <- subset(Degs1, abs(log2FoldChange) > 0.585)
head(sub_Degs1)
dim(sub_Degs1)
rownames(sub_Degs1)

# Normalized count data
Data_org0 <- read.delim("GSE203082_rnaseq_count_matrix.tsv", header = TRUE)
head(Data_org0)
dim(Data_org0)

library(vsn)
Data_org <- justvsn(as.matrix(Data_org0))
dim(Data_org)
head(Data_org, 30)
tail(Data_org)

###### Adjust Labels #####
colnames(Data_org) <- c("M_1442_KCl_1hr", "M_1442_KCl_6hr", "M_1442_PBS", "M_2513_KCl_1hr", "M_2513_KCl_6hr", "M_2513_PBS",
                        "M_2620_KCl_1hr", "M_2620_KCl_6hr", "M_2620_PBS", "M_2962_KCl_1hr", "M_2962_KCl_6hr", "M_2962_PBS",
                        "CTR_3084_KCl_1hr", "CTR_3084_KCl_6hr", "CTR_3084_PBS", "CTR_3130_KCl_1hr", "CTR_3130_KCl_6hr", 
                        "CTR_3130_PBS","CTR_3234_KCl_1hr", "CTR_3234_KCl_6hr", "CTR_3234_PBS", "M_499_KCl_1hr", 
                        "M_499_KCl_6hr", "M_499_PBS","CTR_553_KCl_1hr", "CTR_553_KCl_6hr", "CTR_553_PBS", "M_581_KCl_1hr", 
                        "M_581_KCl_6hr", "M_581_PBS","CTR_690_KCl_1hr", "CTR_690_KCl_6hr", "CTR_690_PBS")
E <- Data_org
head(E)

# Boxplot ordered data
par(mar = c(7, 2, 2, 2))
boxplot(E, las = 2, cex.axis = 0.7)
dim(E)

# Samples / Controls
samples <- E[, grepl("M", colnames(E))]
head(samples)
controls <- E[, grepl("CTR", colnames(E))]
head(controls)

##### New Data Frame and Ordered Data #####

EM <- cbind(controls, samples)
head(EM)
dim(EM)
boxplot(EM, las = 2, cex.axis = 0.7)

# meanSdPlot

data_matrix <- as.matrix(EM)
dim(data_matrix)
meanSdPlot(data_matrix[1:19031,], ranks = TRUE, main = "After Normalization")

##### Filter DEGs #####

# Extract gene names from the filtered table Sub_Degs1

# Filter DEGs

# With fold change gens_DEGsFC

gens_DEGsFC <- rownames(sub_Degs1)
cont_gens_FC <- EM[rownames(EM) %in% gens_DEGsFC, ]
head(cont_gens_FC, 3) # Count of DEGs
dim(cont_gens_FC)

# Without fold change gens_DEG
gens_DEGs <- rownames(Degs1)
Cont_05 <- EM[rownames(EM) %in% gens_DEGs, ]
head(Cont_05)
dim(Cont_05)

# Plots

corrplot(cor(cont_gens_FC), cl.lim = c(0.95, 1), is.corr = FALSE, tl.col = "black",
         col = colorRampPalette(c("blue", "white", "black"))(200))
pheatmap(cont_gens_FC[1:69, ], 
         cluster_rows = TRUE, cluster_cols = TRUE,
         fontsize = 7)

pheatmap(Cont_05, 
         cluster_rows = TRUE, cluster_cols = TRUE,
         fontsize = 4)

boxplot(cont_gens_FC, las = 2)

##### For Annotation #####

library(AnnotationDbi)

columns(org.Hs.eg.db)
ensids <- gens_DEGsFC
cols <- c("GO", "PATH", "PFAM", "GENENAME")
select(org.Hs.eg.db, keys = ensids, columns = cols, keytype = "ENSEMBL") # Generate sub-table

GenomicVariables <- select(org.Hs.eg.db, keys = ensids, columns = cols, keytype = "ENSEMBL")
write.table(GenomicVariables, "GenomicVariables_Annotation.csv")
GenomicVariables <- read.table("GenomicVariables_Annotation.csv")

dim(GenomicVariables)
head(GenomicVariables, 10)
tail(GenomicVariables)
PATH <- as.factor(GenomicVariables$PATH) 
GO <- as.factor(GenomicVariables$GO)  

par(mar = c(7, 3, 2, 2) + 0.3)
barplot(summary(PATH), las = 2)
barplot(summary(GO), las = 2)

###### PCA Data #####

# Data visualization
# Cont_05: data with p > 0.05
Z_05 <- Cont_05
dim(Z_05)
head(Z_05)

# Cont_gens_FC: data with p > 0.05 and log2foldchange > abs(0.585) 
Z <- cont_gens_FC
dim(Z)
class(Z)
head(Z)

geneSymbols <- Z[, 0]
geneSymbols_05 <- Z_05[, 0]

par(mar = c(9, 3, 3, 2) + 0.3)
boxplot(Z_05, las = 2)
title("Genes with (p < 0.05)")
boxplot(Z, las = 2)
title("DEGs (p > 0.05) log2FoldChange > |0.585|)")
meanSdPlot(Z)
meanSdPlot(Z_05)

####### PCA for Z_05 (Cont_05: data with p > 0.05) #####

library(ade4)
par(mar = c(3, 3, 3, 1) + 0.8)
pcalym <- dudi.pca(Z_05) 
names(pcalym)

# Calculate the % variance
pve <- 100 * pcalym$eig / sum(pcalym$eig)
pve
cumsum(pve)
head(pcalym$li)
plot(pcalym$li) # y-axis fold change, x-axis expression level

# Improve filtering by data -> DEGs at 0.05

s.corcircle(pcalym$co) 
corrplot::corrplot(cor(Z_05))

# K-means clustering
set.seed(42) 
cl <- kmeans(Z_05, 3)
class(cl)
dim(Z_05)
names(cl)
table(cl$cluster)
s.class(pcalym$li, as.factor(cl$cluster), col = c(1, 2, 3, 4, 5))

####### PCA for Z (Cont_05: data with p > 0.05 and FC |0.585|) #####

par(mar = c(3, 3, 3, 1) + 0.8)
pcalymFC <- dudi.pca(Z)
names(pcalymFC)

# Calculate the % variance
pve <- 100 * pcalymFC$eig / sum(pcalymFC$eig)
pve
cumsum(pve)
head(pcalymFC$li)
plot(pcalymFC$li) # y-axis fold change, x-axis expression level

par(mar = c(3, 4, 6, 6) + 0.8)
s.corcircle(pcalymFC$co) 
corrplot::corrplot(cor(Z))

##### K-means PCA for Z (data with p > 0.05 and FC |0.585|) ######

set.seed(42) 
cl_05FC <- kmeans(Z, 3)
class(cl_05FC)
dim(Z)
table(cl_05FC$cluster)
names(cl_05FC)
s.class(pcalymFC$li, as.factor(cl_05FC$cluster), col = c(1, 2, 3)) 

# Cluster analysis
cluster_fc <- cl_05FC$cluster
table(cl$cluster)
genes_clusters <- data.frame(Genes = rownames(Z), Cluster = cluster_fc)
head(genes_clusters)

##### PCA and GO Mapping #####

columns(org.Hs.eg.db)
dst <- gens_DEGsFC
cols_GO <- c("GO")
GO_result <- select(org.Hs.eg.db, keys = dst, columns = cols_GO, keytype = "ENSEMBL") 
dim(GO_result)
head(GO_result)
go_result <- as.data.frame(GO_result)

### Combine DataFrames PCA and GO ###

library(tidyverse)
Tab_GO_PCA <- merge(x = GO_result, y = genes_clusters, by.x = "ENSEMBL", by.y = "Genes", all.x = TRUE)
# Grouping
grouped <- Tab_GO_PCA %>%
  group_by(Cluster, ONTOLOGY) %>%
  summarise(cnt = n())
# Totals
totals <- Tab_GO_PCA %>%
  group_by(Cluster) %>%
  summarise(total = n())
# Combine grouped and totals
Tab_GO_PCA <- merge(x = grouped, y = totals, by = "Cluster", all.x = TRUE)
Tab_GO_PCA$percentage <- round(100 * Tab_GO_PCA$cnt / Tab_GO_PCA$total)
Tab_GO_PCA

##### Clusters by Expression Level PATH #####

library(org.Hs.eg.db)
columns(org.Hs.eg.db)
dst <- gens_DEGs
cols_path <- c("PATH")

PATH_result <- select(org.Hs.eg.db, keys = dst, columns = "PATH", keytype = "ENSEMBL")
dim(PATH_result)
names(PATH_result)
head(PATH_result)
PATH_result <- as.data.frame(PATH_result)
summary(as.factor(PATH_result$PATH))
barplot(summary(na.omit(as.factor(PATH_result$PATH))), las = 2)

indpath_MAPK <- which(PATH_result$PATH == "04010") 
length(indpath_MAPK)
genesMAPK <- PATH_result$ENSEMBL[indpath_MAPK]
genesMAPKT <- as.data.frame(genesMAPK)

dim(Z_05)
factorMapk <- rep(0, 634)
ind_table <- which(geneSymbols %in% genesMAPK)
factorMapk[ind_table] <- 1

set.seed(42) 
cl_05 <- kmeans(Z_05, 3)

s.class(pcalym$li, as.factor(cl_05$cluster), col = c(1, 2))
table(factorMapk, cl_05$cluster)

cluster_05 <- cl_05$cluster
genes_clustersP <- data.frame(Genes = rownames(Z_05), Cluster = cluster_05)
head(genes_clustersP)
dim(genes_clustersP)

Tab_PCA_PATH <- merge(x = genesMAPKT, y = genes_clustersP, by.x = "genesMAPK", by.y = "Genes", all.x = TRUE)

grouped <- Tab_PCA_PATH %>%
  group_by(Cluster) %>%
  summarise(cnt = n())

##### PCA and Grouping on Samples #####
pcalym2 <- dudi.pca(t(Z), scannf = FALSE, nf = 2)
names(pcalym2)
plot(pcalym2$li)
head(Z)

sample_names <- rownames(t(Z))
colors <- ifelse(grepl("M_", sample_names), "red", "black")
plot(pcalym2$li, col = colors)
legend("topright", legend = c("Controls", "Patients"), fill = c("red", "black"))

pve <- 100 * pcalym2$eig / sum(pcalym2$eig)
cumsum(pve)

cl <- kmeans(t(Z), 2)
s.class(pcalym2$li, as.factor(cl$cluster), col = c(1, 2))

# Condition vector
fac <- rep(c(1, 2), c(15, 18))
s.class(pcalym2$li, as.factor(fac), col = c(1, 2))
table(cl$cluster, fac)

#### Enrichment in R ####

#### ClusterProfiler
BiocManager::install("clusterProfiler", force = TRUE)
BiocManager::install("pathview")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

df <- sub_Degs1
df$X <- rownames(df)
df <- df[, c("X", names(df)[-which(names(df) == "X")])]
names(df)
head(df)
dim(df)

original_gene_list <- df$log2FoldChange
names(original_gene_list) <- df$X
gene_list <- na.omit(original_gene_list)
gene_list <- sort(gene_list, decreasing = TRUE)
gene_list

gse <- gseGO(geneList = gene_list, 
             ont = "ALL", 
             keyType = "ENSEMBL", 
             nPerm = 1000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

dotplot(gse, showCategory = 10, font.size = 10, split = ".sign")
dotplot(gse, showCategory = 10, font.size = 10, split = ".sign") + facet_grid(. ~ .sign)
dotplot(gse, showCategory = 5, font.size = 10, split = ".sign")

### DESeq2 Workshop Script M.Rojo ####
# This script performs differential expression analysis using DESeq2 on RNA-Seq data.
# It includes loading data, pre-filtering, adjusting labels, creating new data frames,
# generating initial plots, and performing DESeq2 analysis with and without considering time as a factor.

# Clear current environment
rm(list=ls(all=TRUE))

##### Packages ####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("apeglm")

install.packages("umap")
install.packages("car")

require(DESeq2)
require(ggplot2)
require(corrplot)
library(apeglm)

#### Data ####

# Load data
# Link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203082
setwd("/Users/MRo/Desktop/git_RedMra/workshops")
tbl_gen <- read.delim("GSE203082_rnaseq_count_matrix.tsv", stringsAsFactors = FALSE)
head(tbl_gen)
dim(tbl_gen)
row.names(tbl_gen)

# Pre-filtering 

# Pre-filter low count genes
# Keep genes with at least 2 counts > 10
keep <- rowSums(tbl_gen >= 10) >= 2
tbl_gen <- tbl_gen[keep, ]
head(tbl_gen)
dim(tbl_gen)

##### High Reads Filter #####
E = tbl_gen
modE = tbl_gen[1e5 > apply(E, 1, max),] # Remove genes with very high reads
modE = modE[0 < apply(modE[, 16:33], 1, var),] # Remove genes with no variance

dim(modE)
head(modE)

###### Adjust Labels and New Data Frame #####

# Expression matrix

colnames(modE) <- c("M_1442_KCl_1hr", "M_1442_KCl_6hr", "M_1442_PBS", "M_2513_KCl_1hr", "M_2513_KCl_6hr", "M_2513_PBS",
                    "M_2620_KCl_1hr", "M_2620_KCl_6hr", "M_2620_PBS", "M_2962_KCl_1hr", "M_2962_KCl_6hr", "M_2962_PBS",
                    "CTR_3084_KCl_1hr", "CTR_3084_KCl_6hr", "CTR_3084_PBS", "CTR_3130_KCl_1hr", "CTR_3130_KCl_6hr", 
                    "CTR_3130_PBS", "CTR_3234_KCl_1hr", "CTR_3234_KCl_6hr", "CTR_3234_PBS", "M_499_KCl_1hr", 
                    "M_499_KCl_6hr", "M_499_PBS", "CTR_553_KCl_1hr", "CTR_553_KCl_6hr", "CTR_553_PBS", "M_581_KCl_1hr", 
                    "M_581_KCl_6hr", "M_581_PBS", "CTR_690_KCl_1hr", "CTR_690_KCl_6hr", "CTR_690_PBS")
E <- modE
head(E)
samples <- E[, grepl("M", colnames(E))]
head(samples)
controls <- E[, grepl("CTR", colnames(E))]
head(controls)

# Ordered Data Frame
EM <- cbind(controls, samples)
head(EM)
corrplot(cor(EM))

#### Initial Boxplots ####

# Without ordering
par(mar = c(9, 4, 4, 2) + 0.1)
boxplot(modE, las = 2)
boxplot(modE, boxwex = 0.7, notch = TRUE, main = "Boxplot Controls and Samples", outline = FALSE, 
        las = 2, cex.axis = 0.8)

# Initial box-and-whisker plot
par(mar = c(9, 4, 2, 1))
boxplot(EM, las = 2)
boxplot(EM, boxwex = 0.7, notch = TRUE, main = "GSE203082", ylab = "lg(cnt + 1)", outline = FALSE, las = 2)

boxplot(EM, boxwex = 0.7, notch = TRUE, main = "Boxplot Controls and Patients", outline = FALSE, 
        las = 2, cex.axis = 0.8)

# Without ordering 
corrplot(cor(modE), cl.lim = c(0.95, 1), is.corr = FALSE, tl.col = "black",
         col = colorRampPalette(c("blue", "white", "black"))(200))

# Ordering

# Set margins
par(mar = c(7, 4, 4, 2) + 0.1)

corrplot(cor(EM), cl.lim = c(0.95, 1), is.corr = FALSE, tl.col = "black",
         col = colorRampPalette(c("blue", "white", "black"))(200))

head(EM) 
row.names(EM)

##### DESeq2 ####

# Create condition vector
coldata <- data.frame(
  Time_hours = factor(rep(c("1h", "6h", "PBS"), 11)),  
  condition = factor(rep(c(rep('Control', 15), rep('Chil_Esqz', 18))))
)
head(coldata)

rownames(coldata) = colnames(EM)
coldata$condition <- relevel(coldata$condition, ref = 'Control')
coldata

# Model without considering time

# Model for each gene with condition vs control

dds0 <- DESeqDataSetFromMatrix(
  countData = EM[, 1:33],
  colData = coldata[1:33, ],
  design = ~ condition)

######## Normalization 

dds0 <- estimateSizeFactors(dds0)
sizeFactors(dds0)

# Remove samples with size factors greater than 3 or less than 0.3

head(counts(dds0, normalized = TRUE))
boxplot(counts(dds0, normalized = TRUE), las = 2)
corrplot(cor(counts(dds0, normalized = TRUE)))

corrplot(cor(counts(dds0, normalized = TRUE)), cl.lim = c(0.95, 1), is.corr = FALSE, tl.col = "black",
         col = colorRampPalette(c("blue", "white", "black"))(200))

###### Likelihood Ratio Test ######

dds0 <- DESeq(dds0, test = "LRT", reduced = ~1) # Likelihood ratio test compares current fit vs null model
# dds1 <- DESeq(dds0, test="Wald") # Wald test only for two conditions
res0 <- results(dds0)
class(res0)
res0
resultsNames(dds0) # Coefficients
plotDispEsts(dds0) # Dispersion estimates

plotMA(dds0)
resultsNames(dds0)
rownames(dds0)

#### Shrinkage Method ####

# Once the model is fitted, apply the shrinkage method
resLFC0 <- lfcShrink(dds0, coef = "condition_Chil_Esqz_vs_Control", type = "apeglm")

# Function from the "apeglm" package
resLFC0 # Reduced LFC

plotMA(resLFC0) # Scale change

# P-values without shrinkage

Padj0 = data.frame(p = res0$padj[res0$padj < 0.95]) # Creation of data frames 
# Write table
ggplot(Padj0, aes(x = p)) + geom_histogram(breaks = seq(0, 0.95, 0.05)) +
  labs(y = 'Frequency') + scale_x_continuous(breaks = seq(0, 0.95, 0.05))

plotCounts(dds0, gene = which.min(res0$padj), intgroup = "condition")

# Results with shrinkage

Padj0 = data.frame(p = res0$padj[resLFC0$padj < 0.95])

ggplot(Padj0, aes(x = p)) + geom_histogram(breaks = seq(0, 0.95, 0.05)) +
  labs(y = 'Frequency') + scale_x_continuous(breaks = seq(0, 0.95, 0.05))

# Shrinkage made no difference
plotCounts(dds0, gene = which.min(res0$padj), intgroup = "condition")

# To know how many genes are differentially expressed -->

resSig0 <- subset(resLFC0, padj < 0.05)
resSig0 
resOrdered0 <- resSig0[order(resSig0$padj), ]
resOrdered0
dim(resOrdered0)
rownames(resOrdered0)

#####  9 genes with the lowest p-value in the model without time. #####

padjwt = resOrdered0$padj
head(padjwt)
geneswt <- rownames(resOrdered0)[padjwt < 0.01]
head(geneswt)
par(mfrow = c(3, 3), mar = c(2, 2, 1, 1))

for(i in 1:9) {
  plotCounts(dds0, gene = geneswt[i], intgroup = c("condition"), 
             pch = c(20), col = c("black"), 
             legend = TRUE, main = paste("Gene:", geneswt[i]))
}

write.csv(as.data.frame(resSig0), file = "Mod_cond_controls_vs_schizophrenia_patients.csv")

##### Full Design with Time and Condition Factor #####

dds <- DESeqDataSetFromMatrix(
  countData = EM[, 1:33],
  colData = coldata[1:33, ],
  design = ~ Time_hours * condition) # Here is the difference
# This time, time is included as a covariate and many more differences are seen

dds <- DESeq(dds, test = "LRT", reduced = ~1)
# No Wald
res <- results(dds)
resultsNames(dds) # Coefficients intercepts generated with all conditions and times.
res

padj = data.frame(p = res$padj[res$padj < 0.95])
head(padj)

ggplot(padj, aes(x = p)) + geom_histogram(breaks = seq(0, 0.95, 0.05)) +
  labs(y = 'Frequency') + scale_x_continuous(breaks = seq(0, 0.95, 0.05))

plotMA(res, ylim = c(-4, 4))

plotCounts(dds, gene = which.min(res$padj), intgroup = "condition")

##### 9 genes with the lowest p-value with time #####

padj = res$padj
head(padj)
genes = rownames(modE)[order(padj)]
head(genes)

par(mfrow = c(3, 3), mar = c(2, 2, 1, 1))

for(i in 1:9) {
  plotCounts(dds, gene = genes[i], intgroup = c("Time_hours"), 
             pch = c(20, 2), col = c("black", "red"), 
             legend = TRUE, main = paste("Gene:", genes[i]))
}

par(mfrow = c(1, 1))

###### Download CSV ####

resOrdered <- res[order(res$padj), ]
resSig <- subset(resOrdered, padj < 0.05)
resSig
dim(resSig)
row.names(resSig)

write.csv(as.data.frame(resSig), file = "Schizophrenia_results_05.csv")
write.table(resSig, "Schizophrenia_Result_05.tsv", sep = "\t", row.names = TRUE)

#### Heatmap of DEGs ######

# Load necessary packages
install.packages("pheatmap")
library(pheatmap)

# Obtain the first 100 DEGs
DEG_names <- rownames(resSig)[1:100]

# Filter expressions of the first 100 DEGs in the table EM
DEGs_exp <- DEGs_exp[resSig$padj[1:100], ]

# Create heatmap with ordered genes
pheatmap(DEGs_exp, 
         cluster_rows = TRUE, cluster_cols = TRUE,
         fontsize = 8)

# Install the package if not installed
install.packages("openxlsx")

# Load the package
library(openxlsx)

# Save the dataframe as an Excel file
write.xlsx(as.data.frame(resSig), file = "Schizophrenia_condition_results_0.01.xlsx", rowNames = TRUE)

### Workshop Script 2_MR ####
# This script processes microarray data, including loading data, initial plotting, descriptive statistics,
# summarization from probes to genes, normalization, and correlation analysis.

# Clear current environment
rm(list=ls(all=TRUE))

######## Bioconductor Packages ####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("affy")
install.packages("simpleaffy")

BiocManager::install("CLL")
BiocManager::install("Rcpp")
BiocManager::install("affydata")
BiocManager::install("vsn")
BiocManager::install("gcrma")
BiocManager::install("latticeExtra")
BiocManager::install("RColorBrewer")
BiocManager::install("oligo")
BiocManager::install("pd.e.coli.2")
BiocManager::install("pd.hg.u133.plus.2")
BiocManager::install("hgu133plus2.db")
BiocManager::install("hgu95av2.db")
BiocManager::install("AnnotationDbi")

# Install CRAN packages
install.packages("Matrix")
install.packages("corrplot")
install.packages("pheatmap")

####### 1. Load Data #####

library(GEOquery)
gset <- getGEO("GSE179789", GSEMatrix = TRUE, getGPL = FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
class(gset)
gset
dim(gset)

###### Initial Plots #####

# Boxplot

# Order samples
final_numbers <- as.numeric(sub(".*_([0-9]+)$", "\\1", colnames(exprs(gset))))
order <- order(final_numbers)
ordered_data <- exprs(gset)[, order]
colors <- ordered_data

# Boxplot 

colors <- c(rep("#FF6F61", 8), rep("#6B5B95", 8), rep("#88B04B", 8))  

boxplot(exprs(gset), boxwex = 0.7, cex.axis = 0.7, notch = TRUE, 
        main = "GSE179789", outline = FALSE, las = 2, col = colors)

legend("topright", legend = c("CHD patients without DLT treatment", 
                              "CHD patients with DLT treatment for 8 weeks", 
                              "Healthy controls"), 
       fill = c("#FF6F61", "#6B5B95", "#88B04B"))
dim(exprs(gset))
class(exprs(gset))

# Correlation 

library(corrplot)
corrplot(cor(exprs(gset)))
corrplot.mixed(cor(exprs(gset)), lower="number", upper="circle", number.cex = 0.6, tl.cex = 0.6)

hist(exprs(gset)[, 20])

head(exprs(gset))

# Calculate the correlation matrix between samples
sample_correlation <- cor(exprs(gset))

# Heatmap
library(pheatmap)
pheatmap(exprs(gset)[1:50,])

# Visualize the heatmap of sample correlation
pheatmap(sample_correlation)

##############

# Descriptive Statistics

rowSums(exprs(gset)) -> row_sums
quantile(row_sums)

# Boxplot for outliers
boxplot.stats(row_sums) # Outlier data
boxplot(row_sums, 
        main = "Boxplot of Gene Expression Sums per Row",
        xlab = "Samples",
        ylab = "Sum of Gene Expressions",
        col = "lightblue",
        border = "darkblue",
        notch = TRUE)  

# Heatmap of outlier values
pheatmap(exprs(gset)[row_sums > 300,])
pheatmap(exprs(gset)[row_sums > 300,][1:50,])

# Probes
my_first_exset <- exprs(gset)
head(my_first_exset)
probes <- row.names(my_first_exset)
probes[1:10]

######### 1.3 Summarization (Probes -> Genes) ########

library("hgu133plus2.db") 
PROBES <- as.character(probes)
keytypes(hgu133plus2.db)
OUT <- select(hgu133plus2.db, keys = PROBES, columns = c("SYMBOL", "ENTREZID", "GENENAME"), keytype = "PROBEID")
dim(OUT)
OUT

# Summarization

first <- rep("Affymetrix", 57156)
probe <- cbind(as.character(first), OUT[1:2])  
head(probe)
dim(probe)
summary(probe)
names(probe) <- c("platform", "probeid", "MSU6")
head(probe)

# Function to map probes to genes
do_map <- function(a, platform, probe) {
  probeB <- probe[probe$platform == platform, names(probe) %in% c("probeid", "MSU6")]
  probeC <- as.data.frame(probeB, row.names = probeB$probeid)
  b <- new_cbind(probeC, as.data.frame(a))
  c <- b[, -c(1, 2)]
  c <- as.matrix(c)
  c <- apply(c, 2, as.numeric)
  rownames(c) <- b[, 2]
  c <- aggregate(c, list(rownames(c)), max)
  rownames(c) <- c[, 1]
  c <- c[, -1]
  return(c)
}

# Function to combine tables by row names
new_cbind <- function(...){
  input <- eval(substitute(list(...), env = parent.frame()))  
  names.orig <- NULL
  nrows <- numeric()
  for (i in 1:length(input)) {
    nrows[i] <- nrow(input[[i]])
    names.orig <- c(names.orig, colnames(input[[i]])) 
  }  
  idx <- (1:length(input))[order(nrows, decreasing = TRUE)]
  x <- NULL
  for (i in 1:length(input)) {
    x <- c(x, rownames(input[[idx[i]]]))
  }  
  r <- data.frame(row.names = unique(x))
  for (i in 1:length(input)) {
    r <- cbind(r, data.frame(input[[i]][match(rownames(r), rownames(input[[i]])), ]))
  }  
  colnames(r) <- names.orig  
  return(r)
}

# Replace probes with genes -> Table
E <- do_map(my_first_exset, "Affymetrix", probe)
# Check dimensions and first genes stored in E:
dim(my_first_exset)
dim(E) # Observe the reduction in rows, now genes
head(E)
pheatmap(E[1:50,])

######## Other Plots. 

boxplot(E, boxwex = 0.7, cex.axis = 0.7, notch = TRUE, 
        main = "GSE179789", outline = FALSE, las = 2, col = colors)

legend("topright", legend = c("CHD patients without DLT treatment", 
                              "CHD patients with DLT treatment for 8 weeks", 
                              "Healthy controls"), 
       fill = c("#FF6F61", "#6B5B95", "#88B04B"))
dim(E)
class(E)

# Correlation 
library(corrplot)
corrplot(cor(E))
corrplot.mixed(cor(E), lower = "number", upper = "circle", number.cex = 0.6, tl.cex = 0.6)
hist((E)[, 20])

####### NORMALIZATION #######

##### Normalization ######

library(vsn) # VSN method by Huber
citation("vsn")
my_first_exset_norm <- justvsn(my_first_exset)
write.csv(my_first_exset_norm, "my_first_exset_norm.csv")
class(my_first_exset_norm)

# Standard deviation plot
meanSdPlot(my_first_exset[1:20000,], ranks = TRUE, main = "Before Normalization")
meanSdPlot(my_first_exset_norm[1:20000,], ranks = TRUE, main = "After Normalization")

# Boxplot of normalized data
boxplot(my_first_exset_norm, boxwex = 0.8, cex.axis = 0.7, notch = TRUE, main = "Boxplot of Normalized Data",  outline = FALSE, las = 2, col = colors)

legend("topright", legend = c("CHD patients without DLT treatment", 
                              "CHD patients with DLT treatment for 8 weeks", 
                              "Healthy controls"), 
       fill = c("#FF6F61", "#6B5B95", "#88B04B")) 

# Heatmap of normalized data
library(pheatmap)
pheatmap(my_first_exset_norm[1:50,])
pheatmap(exprs(gset)[1:50,])

# Correlation 
library(corrplot)
corrplot(cor(my_first_exset_norm))
corrplot.mixed(cor(my_first_exset_norm), lower = "number", upper = "circle", number.cex = 0.6, tl.cex = 0.6)

hist(exprs(gset)[, 20])
head(my_first_exset_norm)
head(gset)

######## Sum Probes / Outliers #####
rowSums(my_first_exset_norm) -> row_sums_norm
quantile(row_sums_norm)

# Boxplot for outliers
boxplot.stats(row_sums_norm) # Outlier data 
boxplot(row_sums_norm, 
        main = "Boxplot of Gene Expression Sums per Row",
        xlab = "Samples",
        ylab = "Sum of Gene Expressions",
        col = "lightblue",
        border = "darkblue",
        notch = TRUE)  

# Heatmap of outlier values 
pheatmap(my_first_exset_norm[row_sums > 108,])
pheatmap(my_first_exset_norm[row_sums > 108,][1:100,])

# Probes
probes_norm <- row.names(my_first_exset_norm)
probes_norm[1:10]

###### Summarization ####

# Replace probes with genes -> Table in normalized data
EN <- do_map(my_first_exset_norm, "Affymetrix", probe)
# Check dimensions and first genes stored in E:
dim(my_first_exset_norm)
dim(EN) # Observe the reduction in rows, now genes
head(EN)
pheatmap(EN[1:50,])

# Visualize the heatmap of sample correlation
# Order samples
final_numbers <- as.numeric(sub(".*_([0-9]+)$", "\\1", colnames(exprs(gset))))
order <- order(final_numbers)
ordered_data_norm <- my_first_exset_norm[, order]

# Order the correlation matrix according to the same order as the ordered data
sample_correlation_norm <- sample_correlation[order, order]

# Create the heatmap
pheatmap(sample_correlation_norm)
 
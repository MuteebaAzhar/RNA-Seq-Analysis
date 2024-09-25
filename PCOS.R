
# RNA-Seq Differential Expression Analysis using DESeq2

# Install required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("GEOquery", "edgeR", "DESeq2", "org.Hs.eg.db", "clusterProfiler", "pheatmap", "EnhancedVolcano"), force = TRUE)

# Load necessary libraries
library(tidyverse)
library(GEOquery)
library(edgeR)
library(DESeq2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(pheatmap)
library(EnhancedVolcano)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(RColorBrewer)

# Load RNA-seq dataset
GSE138518_RNA <- read_csv("path/to/your/GSE138518_RNA.csv")
head(GSE138518_RNA)  # View first few rows

# Prepare count matrix
count_columns <- c("N20", "N21", "N25", "P14", "P15", "P16")
count_matrix <- GSE138518_RNA[, count_columns]
rownames(count_matrix) <- GSE138518_RNA$ENSEMBL

# Create sample metadata
sample_info <- data.frame(
  Sample = c("N20", "N21", "N25", "P14", "P15", "P16"),
  Condition = factor(c("Normal", "Normal", "Normal", "PCOS_patient", "PCOS_patient", "PCOS_patient"))
)
rownames(sample_info) <- sample_info$Sample

# Check that row names match column names
stopifnot(all(rownames(sample_info) == colnames(count_matrix)))

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = sample_info, design = ~Condition)

# Preprocessing: Filter low count reads
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Run differential expression analysis
dds$Condition <- relevel(dds$Condition, ref = "Normal")
deg <- DESeq(dds)
res <- results(deg)
write.csv(res, "DESeq2_results.csv")

# Filter significant results (adjusted p-value < 0.05)
res_0.05 <- results(deg, alpha = 0.05)
write.csv(as.data.frame(res_0.05), "DESeq2_significant_results.csv")

# Gene annotation
res_0.05$Symbol <- mapIds(org.Hs.eg.db, rownames(res_0.05), keytype = "ENSEMBL", column = "SYMBOL")
write.csv(as.data.frame(res_0.05), "DESeq2_annotated_results.csv")

# PCA Plot
vsd <- vst(deg, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "Condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, label = name)) + 
  geom_point(size = 3) + 
  geom_text_repel() + 
  labs(title = "PCA Plot", x = paste0("PC1: ", percentVar[1], "% variance"), y = paste0("PC2: ", percentVar[2], "% variance"))

# Volcano Plot
res_df <- as.data.frame(res)
res_df$Expression <- ifelse(res_df$log2FoldChange > 1 & res_df$padj < 0.05, "Upregulated", 
                            ifelse(res_df$log2FoldChange < -1 & res_df$padj < 0.05, "Downregulated", "Not Significant"))

ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = Expression)) +
  geom_point() +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  theme_classic() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 p-value") +
  theme(legend.title = element_blank())


# Generate MA plot of differential expression results
plotMA(res)

########## Cook's Distance Outlier Detection ##########
# Calculate Cook's distance to detect outlier genes
cooks_dist <- assays(dds)[["cooks"]]
n <- ncol(dds)
cooks_cutoff <- 4 / n

# Identify and filter out outlier genes
outlier_genes <- apply(cooks_dist, 1, function(x) any(x > cooks_cutoff))
outlier_genes

filtered_dds <- dds[!outlier_genes, ]

# Re-run DESeq analysis on filtered dataset
filtered_dds <- DESeq(filtered_dds)
res_filtered <- results(filtered_dds)
summary(res_filtered)

# Generate filtered MA plot
plotMA(res_filtered, main="MA Plot", ylim=c(-6, 6))

########## Volcano Plot ##########
# Convert results to a dataframe
res_filtered_df <- as.data.frame(res_filtered)

# Classify genes as Up, Down, or Not Significant
res_filtered_df <- res_filtered_df %>%
  mutate(Expression = case_when(
    log2FoldChange > 1 & pvalue < 0.05 ~ "Up Regulated",
    log2FoldChange < -1 & pvalue < 0.05 ~ "Down Regulated",
    TRUE ~ "Not Significant"
  ))

# Check column names in res_filtered_df
colnames(res_filtered_df)

# Add gene symbols
res_filtered_df$Symbol <- rownames(res_filtered_df)

# Define custom colors for volcano plot
custom_colors <- c("Up Regulated" = "red", "Down Regulated" = "blue", "Not Significant" = "grey")

# Plot the volcano plot using ggplot2
ggplot(data = res_filtered_df, aes(x = log2FoldChange, y = -log10(pvalue), color = Expression)) +
  geom_point() +
  scale_color_manual(values = custom_colors) +
  theme_classic() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 p-value") +
  geom_text_repel(aes(label = Symbol), max.overlaps = 10) +
  scale_x_continuous(limits = c(-10, 10), breaks = seq(-10, 10, by = 2)) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 20, by = 2))


# Save top DEGs
top_genes <- res_df %>% arrange(padj) %>% head(30)
write.csv(top_genes, "top_DEGs.csv")
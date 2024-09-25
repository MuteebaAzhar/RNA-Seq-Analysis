# Differential Gene Expression Analysis of Polycystic Ovary Syndrome (PCOS) Using RNA-Seq Data and DESeq2

## Project Overview

This project aims to analyze the differential gene expression in human ovarian granulosa cells in the context of Polycystic Ovary Syndrome (PCOS). Using RNA-seq data from the GEO database (GSE138518), we identified key genes and pathways involved in PCOS pathogenesis. By leveraging the DESeq2 package, we conducted a robust differential expression analysis to uncover potential therapeutic targets for treating PCOS.

The analysis focuses on filtering outliers using Cook's distance and visualizing results through MA plots, volcano plots, and identifying top differentially expressed genes. This project is significant in understanding molecular mechanisms underlying PCOS and exploring novel therapeutic avenues.

## Key Features
- **RNA-Seq Data Analysis**: Raw sequencing data processed through DESeq2 for normalization, differential expression analysis, and visualization.
- **Cook's Distance for Outlier Detection**: Outlier genes identified and filtered using Cook's distance, ensuring accurate results.
- **MA and Volcano Plots**: Clear visual representation of the differential expression results using MA plots and volcano plots.
- **Top Gene Selection**: Identification and saving of the top 30 differentially expressed genes based on adjusted p-values.
  
## Tools and Technologies
- **R programming language**
- **DESeq2**: For differential gene expression analysis.
- **ggplot2 and ggrepel**: For advanced data visualization.
- **dplyr**: For data manipulation.
- **EnhancedVolcano**: For volcano plot creation.
  
## Analysis Workflow

1. **Data Preprocessing**:
    - RNA-seq data from GSE138518 was downloaded and prepared for analysis.
    - DESeq2 was used for normalization and differential gene expression analysis.

2. **Outlier Detection**:
    - Cook's distance was calculated to detect outlier genes, which were removed to refine the analysis.

3. **Differential Expression**:
    - Differential expression analysis was performed using DESeq2, with significance determined based on adjusted p-values and log2 fold changes.
    - Key genes with significant differential expression were identified.

4. **Visualization**:
    - MA plots and volcano plots were generated to visualize gene expression changes.
    - The top 30 most significant genes were saved for further investigation.

## Visualizations

### MA Plot
The MA plot provides a quick overview of the distribution of log2 fold changes across the mean expression values, highlighting up- and down-regulated genes.

### Volcano Plot
The volcano plot visualizes the significance (p-value) against log2 fold change of gene expression, allowing for the clear identification of the most differentially expressed genes.

## Results

The analysis identified several genes that are differentially expressed between PCOS and control samples. These genes are involved in key biological pathways such as hormonal regulation, inflammatory responses, and metabolic processes. The filtered dataset, after removing outlier genes, improved the reliability of the results. The top 30 differentially expressed genes were saved in the `top_DEGs.csv` file.

## Acknowledgments

The RNA-seq data for this project was sourced from the publicly available GEO dataset [GSE138518](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138518). We gratefully acknowledge the original contributors of the dataset.



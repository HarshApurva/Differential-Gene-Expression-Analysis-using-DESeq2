🧬 Differential Gene Expression Analysis using DESeq2

This repository provides a complete R pipeline for performing Differential Gene Expression Analysis (DGEA) using the DESeq2 package. It supports raw non-decimal count data from RNA-seq (e.g., from HTSeq, featureCounts, etc.) and compares control vs disease conditions.


🚀 Features
Input: CSV file with raw integer counts (genes as rows, samples as columns)
Automatically installs and loads required packages

Generates:
✅ MA Plot
✅ Volcano Plot
✅ PCA Plot
✅ Heatmap of top differentially expressed genes

Outputs: differential_expression_results.csv with log2FC, p-values, padj, etc.


📂 How to Use
Step 1: Modify These Two Things in the Script
Path to your raw count file (in CSV format):
count_file <- "C:/path/to/your_file.csv"
Sample conditions (e.g., "control" or "treatment"), in the order your columns appear:

condition_groups <- c("control", "control", "treatment", "treatment")
⚠ Make sure your count data only contains non-decimal integer values, and column order matches the condition labels!


📦 Dependencies
The script uses the following R packages:
DESeq2
ggplot2
pheatmap
RColorBrewer
EnhancedVolcano
These will be automatically installed if not present.

📊 Sample Output
The script will create a results/ folder containing:
MA_plot.png
volcano_plot.png
heatmap.png
PCA_plot.png
differential_expression_results.csv

📎 Sample outputs are attached to this repo for reference.



💡 Notes
Input must be raw (unnormalized) integer counts — not TPM, FPKM, or log-transformed data.
Genes with low counts are filtered before DE analysis.
If your data is merged from multiple files, make sure all sample columns are numeric.

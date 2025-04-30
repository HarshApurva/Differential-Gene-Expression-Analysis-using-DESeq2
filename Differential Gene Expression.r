# Differential Gene Expression (Control vs Disease) using DESeq2

# Required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in c("DESeq2", "pheatmap", "ggplot2", "RColorBrewer")) {
  if (!require(pkg, character.only = TRUE)) BiocManager::install(pkg)
  library(pkg, character.only = TRUE)
}

# === USER INPUT ===
count_file <- "your_file.csv"  # <-- Replace with your CSV path
groups <- c("control", "control", "disease", "disease", "disease")  # match sample order
output_dir <- "results"; dir.create(output_dir, showWarnings = FALSE)

# === Load data ===
counts <- read.csv(count_file, row.names = 1, check.names = FALSE)
counts[] <- lapply(counts, as.integer)
coldata <- data.frame(condition = factor(groups), row.names = colnames(counts))

# === DESeq2 analysis ===
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)
res_df <- as.data.frame(res[order(res$padj), ])
res_df$gene <- rownames(res_df)

# === Save results ===
write.csv(res_df, file.path(output_dir, "differential_expression_results.csv"), row.names = FALSE)

# === MA Plot ===
png(file.path(output_dir, "MA_plot.png")); plotMA(res); dev.off()

# === Volcano Plot ===
vol <- na.omit(res_df)
vol$log10padj <- -log10(vol$padj)
vol$sig <- with(vol, padj < 0.05 & abs(log2FoldChange) > 1)
png(file.path(output_dir, "volcano_plot.png"))
plot(vol$log2FoldChange, vol$log10padj, col=ifelse(vol$sig, "red", "gray"),
     pch=20, main="Volcano Plot", xlab="log2 Fold Change", ylab="-log10(padj)")
abline(h=-log10(0.05), col="blue", lty=2); abline(v=c(-1, 1), col="blue", lty=2)
dev.off()

# === Heatmap (Top DE genes) ===
top_genes <- head(order(res$padj), 30)
mat <- log2(counts(dds, normalized=TRUE)[top_genes, ] + 1)
pheatmap(mat, annotation_col = coldata,
         scale = "row", show_rownames = TRUE, 
         color = colorRampPalette(rev(brewer.pal(7, "RdBu")))(100),
         filename = file.path(output_dir, "heatmap.png"))

# === PCA ===
pca_data <- prcomp(t(log2(counts(dds, normalized=TRUE) + 1)))
png(file.path(output_dir, "PCA_plot.png"))
plot(pca_data$x[,1:2], col=as.numeric(coldata$condition), pch=16,
     xlab="PC1", ylab="PC2", main="PCA Plot")
legend("topright", legend=levels(coldata$condition), col=1:2, pch=16)
dev.off()

cat("✅ Analysis complete! Check the 'results' folder.\n")
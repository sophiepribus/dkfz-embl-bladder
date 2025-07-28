# This code was used for the GSEA analyses set-ups (deseq) (L1s - figure 5E, LOY - supplementary fig. 5C). 

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")

library(DESeq2)
library(dplyr)

# for figure 5E: define L1-high and L1-low samples
coldata <- data.frame(
  sample = c("B42tumor", "B123tumor", "B134tumor", "B39tumor", "B156tumor", "B60tumor"),
  condition = c("L1_high", "L1_high", "L1_high", "L1_low", "L1_low", "L1_low"),
  row.names = c("B42tumor", "B123tumor", "B134tumor", "B39tumor", "B156tumor", "B60tumor")
)

# for supplementary fig. 5C: define LOY and normal samples
coldata <- data.frame(
  sample = c("B134tumor", "B123tumor", "B154tumor", "B74tumor", "B157tumor", "B60tumor", "B24tumor", "B39tumor", "B94tumor", "B125tumor", "B175tumor"),
  condition = c("LOY", "LOY", "LOY", "normal_Y", "normal_Y", "normal_Y", "normal_Y", "normal_Y", "normal_Y", "normal_Y", "normal_Y"),
  row.names = c("B134tumor", "B123tumor", "B154tumor", "B74tumor", "B157tumor", "B60tumor", "B24tumor", "B39tumor", "B94tumor", "B125tumor", "B175tumor")
)

# load merged salmon gene counts file for all samples
counts_df <- read.delim("salmon.merged.gene_counts.tsv", row.names = 1)

counts_df <- counts_df[, rownames(coldata)]
counts_df <- round(counts_df)

# construct DESeq run
dds <- DESeqDataSetFromMatrix(countData = counts_df,
                              colData = coldata,
                              design = ~ condition)

# ref = L1_low for fig. 5E, normal_Y for supp. fig. 5C
dds$condition <- relevel(dds$condition, ref = "")

dds <- DESeq(dds)
res <- results(dds)

# construct ranked genes list by log2FoldChange
ranked_list <- data.frame(
  gene_id = rownames(res),
  log2FoldChange = res$log2FoldChange
)

# sort ranked_list
ranked_list <- ranked_list[!is.na(ranked_list$log2FoldChange), ]
ranked_list <- ranked_list[order(-ranked_list$log2FoldChange), ]
ranked_list_df <- as.data.frame(ranked_list)

# save ranked list table to be loaded into GSEA application (GSEAPreranked)
write.table(ranked_list, "ranked_genes.rnk", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
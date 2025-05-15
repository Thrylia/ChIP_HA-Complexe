# This is an example file
# As each dataset retrieved was in a different format, it was necessary to adapt each time.  


#################################################################################### DESeq2
library(DESeq2)

# Load the count table
counts <- read.csv("GSE108388_baf_complex.mdf.csv", sep = "\t", row.names = 1)

# Select replicate columns (R indexing is 1-based)
selected_counts <- counts[, c(
  2, 3, 4,     # arid2
  5, 6, 7,     # brd7
  8, 9, 10,    # pbrm1
  11, 12, 13,  # smarca2
  14, 16, 18,  # smarca4_2878_4
  15, 17, 19,  # smarca4_2878_6
  20, 22, 24   # wt
)]

# Rename columns for clarity
colnames(selected_counts) <- c(
  paste0("arid2_rep", 1:3),
  paste0("brd7_rep", 1:3),
  paste0("pbrm1_rep", 1:3),
  paste0("smarca2_rep", 1:3),
  paste0("smarca4_4_rep", 1:3),
  paste0("smarca4_6_rep", 1:3),
  paste0("wt_rep", 1:3)
)

# Define experimental conditions
condition <- factor(c(
  rep("arid2", 3),
  rep("brd7", 3),
  rep("pbrm1", 3),
  rep("smarca2", 3),
  rep("smarca4_4", 3),
  rep("smarca4_6", 3),
  rep("wt", 3)
))

# Create metadata DataFrame
coldata <- data.frame(row.names = colnames(selected_counts), condition = condition)

# Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = selected_counts,
                              colData = coldata,
                              design = ~ condition)

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)

# Extract results
res_arid2 <- results(dds, contrast = c("condition", "arid2", "wt"))
res_brd7 <- results(dds, contrast = c("condition", "brd7", "wt"))
res_pbrm1 <- results(dds, contrast = c("condition", "pbrm1", "wt"))
res_smarca2 <- results(dds, contrast = c("condition", "smarca2", "wt"))
res_smarca4_4 <- results(dds, contrast = c("condition", "smarca4_4", "wt"))
res_smarca4_6 <- results(dds, contrast = c("condition", "smarca4_6", "wt"))

# Save the DESeq2 results to tab-separated files without quotes
write.table(as.data.frame(res_arid2), file = "res_arid2_vs_wt.tsv", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(as.data.frame(res_brd7), file = "res_brd7_vs_wt.tsv", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(as.data.frame(res_pbrm1), file = "res_pbrm1_vs_wt.tsv", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(as.data.frame(res_smarca2), file = "res_smarca2_vs_wt.tsv", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(as.data.frame(res_smarca4_4), file = "res_smarca4_4_vs_wt.tsv", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(as.data.frame(res_smarca4_6), file = "res_smarca4_6_vs_wt.tsv", sep = "\t", quote = FALSE, row.names = TRUE)






#################################################################################### Volcano
library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)
library(stringr)

# Thresholds
pvalue_threshold <- 0.05
log2FC_threshold <- 0.58

# Read DESeq2 results
data <- read_tsv("DESeq2_results_PerKO_vs_WT.csv", col_names = TRUE)

# Rename first column if empty or contains "1"
if (colnames(data)[1] == "" || grepl("1", colnames(data)[1])) {
  colnames(data)[1] <- "gene"
}

# Read genes to annotate
annotations <- read_tsv("KEGG_ensembl-annot.csv", col_names = FALSE)
kept_genes <- tolower(str_trim(annotations[[2]]))

# Prepare data
data <- data %>%
  mutate(
    gene_lower = tolower(str_trim(gene)),
    is_in_list = gene_lower %in% kept_genes,
    is_significant = pvalue < pvalue_threshold & abs(log2FoldChange) > log2FC_threshold,
    color = case_when(
      is_significant & log2FoldChange > 0 ~ "red",
      is_significant & log2FoldChange < 0 ~ "blue",
      is_in_list ~ "grey",
      TRUE ~ NA_character_
    ),
    label = ifelse(is_in_list, gene, NA_character_)  # Show all gene names from the list
  ) %>%
  filter(is_in_list)  # Keep only genes in the list

# Volcano plot
p <- ggplot(data, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = color), size = 1.5, alpha = 0.7) +
  scale_color_identity() +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = Inf) +
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold),
             linetype = "dashed", color = c("blue", "red"), linewidth = 0.5) +
  geom_hline(yintercept = -log10(pvalue_threshold),
             linetype = "dashed", color = "grey50", linewidth = 0.5) +
  theme_classic() +
  labs(x = "log2 Fold Change", y = "-log10(p-value)")

ggsave("Per1-2_DESeq2_results.filter.pvalue.png", p, width = 7, height = 6, dpi = 300)

# setwd("~/Desktop/depression/E-GEO-24095")

library(limma)
library(clusterProfiler)
library(org.Hs.eg.db) 
library(AnnotationDbi)
library(dplyr)


# pass a character vector of file names
files <- list.files("~/Desktop/depression/E-GEO-24095/data",
                    pattern="\\.gpr$",
                    full.names=TRUE)

## Read .gpr files 
# R->Red channel intensity 
# G->green channel intensity
# Rb,Gb -> background intensity 
# genes -> prob annotation 
RG <- read.maimages(files, source="genepix")
# Background correction use normexp model
RG <- backgroundCorrect(RG, method="normexp", offset=50)
# Normalization M = log2(R/G), A=1/2log2(R*G)
MA <- normalizeWithinArrays(RG, method="loess")
MA <- normalizeBetweenArrays(MA, method="Aquantile")
# extract M and A matrices
# M: log2 fold changes for each array
# A: average log intensity
M_matrix <- as.data.frame(MA$M)
A_matrix <- as.data.frame(MA$A)
# add prefix to column names
colnames(M_matrix) <- paste0("M_", colnames(M_matrix))
colnames(A_matrix) <- paste0("A_", colnames(A_matrix))
# Extract gene annotation 
gene_info <- as.data.frame(MA$genes)
# Combine everything
MA_df <- cbind(gene_info, M_matrix, A_matrix)
head(MA_df)
# Extract only expression matrix
expr_matrix <- MA_df[, grep("^M_", colnames(MA_df))]
expr_matrix <- as.matrix(expr_matrix) # convert to matrix 
rownames(expr_matrix) <- MA_df$ID # set gene IDs as rownames
# Clean column name
colnames(expr_matrix) <- gsub(".*GSM", "GSM", colnames(expr_matrix))



## Prepare sample annotation
# Load sample information to R
targets <- read.delim("E-GEOD-24095.sdrf.txt", 
                      stringsAsFactors = FALSE)

sp_info <- targets[, c(
  "Source.Name",
  "Characteristics..disease.state.",
  "Characteristics..sample.id.",
  "Characteristics..tissue.",
  "Characteristics..age.",
  "Label",
  "Factor.Value..TISSUE.",
  "Factor.Value..DISEASE.STATE."
)]

sp_info$GSM <- sub(" [12]$", "", sp_info$Source.Name)
sp_unique <- sp_info[!duplicated(sp_info$GSM), ]
annotation <- data.frame(Tissue = sp_unique$Factor.Value..TISSUE.)
rownames(annotation) <- sp_unique$GSM
# Reorder to match expression matrix
annotation <- annotation[colnames(expr_matrix), , drop = FALSE]





# Calculate LogFC use Limma
design <- matrix(1, ncol = 1, nrow = ncol(MA$M))
colnames(design) <- "MDDvsCTRL"
fit <- lmFit(MA, design)
fit <- eBayes(fit)
results <- topTable(fit, number = Inf, adjust.method = "BH")
head(results)
## ---- 4) Thresholding & labeling ----
thr <- log2(1.3)     
results$negLogFDR <- -log10(results$adj.P.Val)

results$regulation <- "Not significant"
results$regulation[results$adj.P.Val < 0.05 & results$logFC > thr] <- "Upregulated"
results$regulation[results$adj.P.Val < 0.05 & results$logFC < -thr] <- "Downregulated"
# Remove non-biological Probes
keep <- !(results$ID %in% c("-", "Null", "Dye Marker"))
results <- results[keep, ]

results <- results[results$ID != "" &
    !is.na(results$ID),
]

# Keep probe with lowest FDR per gene
results <- results[order(results$adj.P.Val), ]
results <- results[!duplicated(results$ID), ]
results <- results[, !colnames(results) %in% 
                               c("Block", "Row", "Column", "Name")] # remove unwanted column
results <- results[!is.na(results$ID) &
                               results$ID != "", ] # Remove rows with missing ID
rownames(results) <- results$ID  


# Convert gene name to ENTREZID
results$ENTREZID <- mapIds(
  org.Hs.eg.db,
  keys = results$ID,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

res_mapped <- results[, c(
  names(results)[1],     
  "ENTREZID",                 
  setdiff(names(results), c(names(results)[1], "ENTREZID"))
)]
res_mapped <- res_mapped[!is.na(res_mapped$ENTREZID), ]

write.csv(res_mapped, file = "EGEO24095_full_results_mapped.csv", row.names = TRUE, quote = FALSE)



res_sig <- results[results$regulation != "Not significant", ]


# plot use ggplot2
gene_to_label <- c("GPX4", "TFR2", "ALDH9A1", "HSPA8", "PAWR", 
                   "PTGS2", "MAPK3", "SIRT6", "NOS1")
label_data <- results[results$ID %in% gene_to_label, ]



library(ggplot2)
p <- ggplot(data = results, aes(x = logFC, y = negLogFDR)) +
  geom_point(alpha = 0.6, size = 1.2,aes(color = regulation)) +
  scale_color_manual(values = c(
    "Upregulated" = "#F6C0CC",
    "Downregulated" = "#ABDAEC",
    "Not significant" = "grey"
  )) +
  geom_vline(xintercept = c(-thr, thr), linetype = "dashed", color = "grey70") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey70") +
  geom_text_repel(data = subset(results, ID %in% gene_to_label),
                  aes(label = ID), size = 3, max.overlaps = Inf) + 
  geom_point(data = subset(results, ID %in% gene_to_label),
             aes(x = logFC, y = negLogFDR), color = "black", fill = NA, shape = 21, size = 2) + 
  labs(
    title = "Micrroarray: MDD vs Control",
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)",
    color = ""
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.text = element_text(size = 11)
  )
print(p)

# Save 
ggsave("Fig2F_volcano.jpeg", plot = p, width = 4.76, height = 4.88, dpi = 300)





# Extract significant DEGs

res_sig <- results_mapped[results_mapped$regulation != "Not significant", ]







# Heatmap 
# Keep only significant genes
sig_genes <- results[
  results$adj.P.Val < 0.05 &
    abs(results$logFC) > thr,   
]
# Extract Expression Matrix for significant genes
expr_sig <- expr_matrix[rownames(expr_matrix) %in% sig_genes$ID, ]


# heatmap for genes of interest 

gene_list <- c("GPX4", "TFR2", "ALDH9A1", "HSPA8", "PAWR", 
                   "PTGS2", "MAPK3", "SIRT6", "NOS1", "TXNIP",
               "FTL1", "CD63")


gene_list <- intersect(gene_list, rownames(expr_sig))
expr_subset <- expr_sig[gene_list, ]

library(RColorBrewer)

breaks <- seq(-2, 2, length.out = 101)

colors <- colorRampPalette(c("#67A5cc", "white", "#E282A7"))(100)
p <- pheatmap(expr_subset,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colors,
         breaks = breaks,
         show_colnames = FALSE,
         main = "Microarray: MDD vs Control",
         legend = TRUE,
         legend_breaks = c(-4, -2, 0, 2, 4),
         legend_labels = c("-4", "-2", "0", "2", "4"),
         annotation_legend = TRUE,
         name = "LogFC")


pheatmap::pheatmap(expr_subset,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   color = colors,
                   breaks = breaks,
                   show_colnames = FALSE,
                   main = "Microarray: MDD vs Control",
                   legend = TRUE,
                   legend_breaks = c(-4, -2, 0, 2, 4),
                   legend_labels = c("-4", "-2", "0", "2", "4"),
                   annotation_legend = TRUE,
                   name = "LogFC",
                   filename = "heatmap_array.jpeg",
                   width = 9.4,
                   height = 3.01,
                   res = 300)



jpeg("heatmap_array.jpeg", width = 9.4, height = 3.01, units = "in", res = 300)
grid::grid.draw(p$gtable)
dev.off()

# Heatmap
scaled_data <- t(scale(t(heatmap_data)))  # Scale across rows
library(pheatmap)




# GESA KEGG use .rnk file
gene_rank <- res_mapped[, c("ID", "t")]
gene_rank <- gene_rank[!duplicated(gene_rank$ID), ]
gene_rank <- gene_rank[order(gene_rank$t, decreasing = TRUE), ]
write.table(gene_rank,
            "MDD_microarray_GSEA.rnk",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)


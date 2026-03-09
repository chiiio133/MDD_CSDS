##setwd("~/Desktop/depression/GSE109315/sus_vs_ctrl")

df1 <- read.csv("vHPC_log2CPM_CvsS_full_results_mapped.csv", row.names = 1)


## ----add -log10(adjusted p-value) ----
thr <- log2(1.3)     
df1$negLogFDR <- -log10(df1$adj.P.Val)

df1$regulation <- "Not significant"

df1$regulation[
  df1$adj.P.Val < 0.05 & df1$logFC > thr
] <- "Upregulated"

df1$regulation[
  df1$adj.P.Val < 0.05 & df1$logFC < -thr
] <- "Downregulated"





# Select genes to label

genes_to_label <- c("Hspa8", "Pawr", "Aldh1a1", "Txnip", "Ptgs2","Kdm6b",
                     "Mapk3", "Nfe2l1", "Cd36", "Prdx6b", "Edn1", "Slc7a11",
                    "Cyp1b1", "Park7", "Prdx4", "Sirt6")


# Volcano plot 

p <- ggplot(data = df1, aes(x = logFC, y = negLogFDR)) + 
  geom_point(aes(color = regulation), 
             stroke = 0, alpha = 0.7, size = 2) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-thr, thr), linetype = "dashed", color = "black") +
  # add frame to genes_to_label
  geom_point(data = subset(df1, SYMBOL %in% genes_to_label),
             aes(x = logFC, y = negLogFDR), color = "black", fill = NA, shape = 21, size = 2) + 
  # add text
  geom_text_repel(data = subset(df1, SYMBOL %in% genes_to_label),
                  aes(label = SYMBOL), size = 3, max.overlaps = Inf) + 
  
  # add annotation

  annotate("text", x = -1.5, y = 8, label = "Downregulated", 
           size = 3.4, color = "black") +
  annotate("text", x = 1.5, y = 8, label = "Upregulated", 
           size = 3.4, color = "black") + 
  scale_color_manual(values = c(
    
  "Downregulated" = "#35978f",
  "Upregulated" = "#bf812d",
  "Not significat" = "gray")) + 
  scale_x_continuous(limits = c(-3, 3)) + 
  labs(title = "vHPC: resilient vs control", x = "Log2 Fold Change", 
       y = "-Log10 (ajusted p-value)")
theme_classic() + 
  theme(panel.grid = element_blank(), 
        plot.title = element_text(
          vjust = 0.5, hjust = 0.5, size = 10, color = "black"), 
        legend.position = "none")
        

ggsave("vHPC_sus_v.jpeg", plot = p, width = 5.4, height = 4.2, dpi = 300)


  
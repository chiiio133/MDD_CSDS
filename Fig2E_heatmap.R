library(gridtext)
library(tibble)

##setwd("~/Desktop/depression/E-GEO-24095")

#for expr_matrix and annotation see EGEO24095_norm_vocalno


#---------prepare the data frame-------
expr_matrix <- expr_matrix[!rownames(expr_matrix) %in% c("Dye Marker", "-"), ]


expr_matrix <- expr_matrix %>%
  as.data.frame() %>%
  rownames_to_column("ID")

annotation <- annotation %>%
  as.data.frame() %>%
  rownames_to_column("sample_id")


# Genes to highlight in the heatmap

htmap_highlight_genes <- c(marker_gene)

# Heatmap visualization of DEGs of human hippocampus microarray ----------


column_title <- c("DG", "CA1")

DEG <- res_mapped %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > thr) %>%
  dplyr::select(ID, regulation)

heatmap_meta <- annotation %>%
  mutate(Tissue = factor(Tissue, c("Dentate Gyrus (DG)", "CA1 subregions of the hippocampus")))


htmap_df <- expr_matrix %>% 
  dplyr::select(ID, heatmap_meta$sample_id) %>%
  inner_join(DEG, by = "ID") %>%
  mutate(fontface = ifelse(ID %in% htmap_highlight_genes, 4, 3)) %>%
  mutate(regulation = factor(regulation, c("Upregulated", "Downregulated"))) %>%
  mutate(col = ifelse(regulation == "Upregulated", "red", "royalblue3"))

htmap_df_scaled <- htmap_df %>%
  select(heatmap_meta$sample_id) %>%
  t() %>% scale() %>% t()

htmap_col <- circlize::colorRamp2(
  breaks = c(-5, 0, 5), colors = c("navy", "white", "firebrick2"))


label_df <- htmap_df %>%
  mutate(index = 1:n()) %>%
  filter(ID %in% htmap_highlight_genes)


annot_row <- rowAnnotation(
  link = anno_mark(
    at = label_df$index,
    labels = label_df$ID,
    labels_gp = gpar(fontsize = 10, fontface = 3)),
  annotation_legend_param = list(title_gp = gpar(fontsize = 10.5, fontface = "bold")))

annot_col <- columnAnnotation(
  Tissue = heatmap_meta$Tissue,
  col = list(
    Tissue = setNames(
      c("#3b58a7", "#ed2e23"),
      c("Dentate Gyrus (DG)", "CA1 subregions of the hippocampus")
    )
  ),
  annotation_name_gp = gpar(fontsize = 0, fontface = "bold"),
  annotation_name_side = "left",
  simple_anno_size = unit(3, "mm"),
  annotation_legend_param = list(
    Tissue = list(labels = gt_render(column_title))
  )
)


p <- Heatmap(
  htmap_df_scaled,
  use_raster = FALSE,
  col = htmap_col,
  row_title = NULL,
  row_split = htmap_df$regulation,
  show_row_names = FALSE,
  show_row_dend = FALSE,
  column_title = gt_render(column_title),
  column_names_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  column_split = heatmap_meta$Tissue,
  show_column_names = TRUE,
  show_column_dend = FALSE,
  top_annotation = annot_col,
  right_annotation = annot_row,
  name = "logFC",
  cluster_columns = FALSE,
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  row_gap = unit(1.5, "mm"),
  column_gap = unit(1.5, "mm")
)

filename <- "htmap_array"


pdf(paste0(filename, ".pdf"), width = 6.22, height = 8.71)
draw(p)
dev.off()

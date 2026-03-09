library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(cowplot)
library(dplyr)
library(RColorBrewer)

## palette
pal_text <- as.list(
  setNames(c("darkgoldenrod2", "orangered", "#FB9A99", "darkred", ggpubfigs::friendly_pals$nickel_five[5],
             brewer.pal(9, "Greens")[5], "forestgreen", "steelblue1", "steelblue4", "grey15"),
           c(paste0(c("Kdm6b", "Mapk3", "Nfe2l1", "Sirt6"), ".up"),
             paste0(c("Hspa8", "Pawr", "Aldh1a1", "Txnip", "Ptgs2"), ".down"))))

pal_text <- as.list(
  setNames(
    c(
      "darkgoldenrod2", "orangered", "#FB9A99", "darkred",
      "#56B4E9",                       # replacement for ggpubfigs color
      brewer.pal(9, "Greens")[5], "forestgreen",
      "steelblue1", "steelblue4", "grey15"
    ),
    c(
      paste0(c("Kdm6b", "Mapk3", "Nfe2l1", "Sirt6"), ".up"),
      paste0(c("Hspa8", "Pawr", "Aldh1a1", "Txnip", "Ptgs2"), ".down")
    )
  )
)




## Function for permutation test comparing 2 sets of DEGs ----------------------
perm_test <- function(set1, set2, all_genes1, all_genes2, name1, name2, 
                      prefix1 = "", prefix2 = "", n_perm = 10000) {
  n_overlap <- inner_join(set1, set2, by = "SYMBOL",
                          suffix = c(".set1", ".set2")) %>%
    mutate(direction.set1 = paste(direction.set1, "in", name1),
           direction.set2 = paste(direction.set2, "in", name2)) %>%
    mutate(direction.set1 = forcats::fct_rev(direction.set1),
           direction.set2 = forcats::fct_rev(direction.set2)) %>%
    group_by(direction.set1, direction.set2) %>%
    summarise(n_overlap = n())
  
  set1_ls <- split(set1, set1$direction)
  set2_ls <- split(set2, set2$direction)
  
  set.seed(42)
  n_overlap_perm <- replicate(n_perm, {
    set1_perm <- lapply(set1_ls, function(set1_sub) {
      sample(all_genes1, size = nrow(set1_sub), replace = FALSE)
    })
    set2_perm <- lapply(set2_ls, function(set2_sub) {
      sample(all_genes2, size = nrow(set2_sub), replace = FALSE)
    })
    
    mapply(
      FUN = function(x, y) {length(intersect(x, y))},
      x = list(up.up = set1_perm$up, down.up = set1_perm$down, 
               up.down = set1_perm$up, down.down = set1_perm$down),
      y = list(up.up = set2_perm$up, down.up = set2_perm$up,
               up.down = set2_perm$down, down.down = set2_perm$down),
      USE.NAMES = TRUE) }) %>% 
    as.data.frame() %>% rownames_to_column("dirs") %>%
    separate(dirs, c("direction.set1", "direction.set2"), sep = "\\.") %>%
    mutate(direction.set1 = paste(direction.set1, "in", name1),
           direction.set2 = paste(direction.set2, "in", name2)) %>%
    mutate(direction.set1 = forcats::fct_rev(direction.set1),
           direction.set2 = forcats::fct_rev(direction.set2)) %>%
    pivot_longer(cols = 3:ncol(.), names_to = "index", values_to = "n_overlap_perm")
  
  p_overlap_perm <- n_overlap_perm %>% 
    left_join(n_overlap, by = c("direction.set1", "direction.set2")) %>%
    group_by(direction.set1, direction.set2) %>%
    summarise(p.perm = sum(n_overlap_perm > n_overlap)/n_perm) %>%
    # mutate(p.perm = ifelse(p.perm > 0.5, 1 - p.perm, p.perm)) %>%
    as.data.frame() %>%
    pivot_wider(names_from = "direction.set1", values_from = "p.perm")
  
  n_bins <- min(max(n_overlap_perm$n_overlap_perm, na.rm = TRUE), 15)
  hist_overlap_perm <- n_overlap_perm %>%
    ggplot(aes(x = n_overlap_perm)) +
    geom_bar(stat = "count", size = 0.4) +
    geom_vline(data = n_overlap, aes(xintercept = n_overlap),
               color = "red", lty = 2) + 
    ggh4x::facet_grid2(direction.set2 ~ direction.set1, 
                       scales = "free", independent = "all") +
    theme_cowplot(font_size = 11)
  
  ggsave(paste0("figures/perm_test/perm_test_", prefix1, name1, "_", prefix2, name2, ".png"), 
         hist_overlap_perm, width = 5, height = 4, dpi = 400)
  return(list(n_overlap = n_overlap, p_overlap_perm = p_overlap_perm,
              hist_overlap_perm = hist_overlap_perm))
}



## Functions for drawing circos plots -------------------------------------
library(dplyr)
library(tidyr)
library(stringr)
library(circlize)
library(grid)

circo_overlap <- function(res1, res2, name1, name2, niceFacing = TRUE, 
                          circo_pal = NULL, gene_width_short = 2, nudge_x = rep(0, 4),
                          suffix = "", Fig_num = "", size, gene_width = 10, gene_cex = 0.6, num_cex = 0.75, p_cex = 0.85,
                          degree = -250, big_gap = 5, small_gap = 1, gene_list = "", step = 150,
                          phago_sig = NULL, show_selected_genes = FALSE, break_p_label = FALSE,
                          overlap_perm_test = NULL, conflict_facing = "outside") {
  
  # use dplyr::select
  res1 <- res1 %>% 
    as.data.frame() %>% 
    dplyr::select(SYMBOL, direction)  
  
  res2 <- res2 %>% 
    as.data.frame() %>% 
    dplyr::select(SYMBOL, direction)  
  
  names <- c(paste0(name2, c(".up", ".down")), paste0(name1, c(".up", ".down")))
  names <- factor(names, names)
  
  overlap_genes <- res1 %>%
    inner_join(res2, by = "SYMBOL", suffix = c(".1", ".2")) %>% 
    mutate(direction.1 = ifelse(is.na(direction.1), NA, paste0(name1, ".", direction.1)),
           direction.2 = ifelse(is.na(direction.2), NA, paste0(name2, ".", direction.2))) %>%
    pivot_longer(cols = 2:3, names_to = "group", values_to = "from") %>%
    dplyr::rename("to" = "SYMBOL") %>%  # 重命名为to，保持逻辑一致
    dplyr::select(-group) %>%  # 明确使用dplyr::select
    mutate(from = str_replace(from, "2$", name2)) %>%
    mutate(value = 1, value2 = gene_width) %>% 
    dplyr::select(from, to, value, value2)  # 明确使用dplyr::select
  
  if(show_selected_genes) {
    overlap_genes <- overlap_genes %>%
      mutate(value2 = ifelse(to %in% c(gene_list, phago_sig), value2, gene_width_short))
  }
  
  count_ds <- rbind(
    dplyr::count(res1, direction) %>% mutate(direction = paste0(name1, ".", direction)),
    dplyr::count(res2, direction) %>% mutate(direction = paste0(name2, ".", direction))) %>%
    dplyr::rename("from" = "direction") %>% 
    full_join(dplyr::count(overlap_genes, from), by = "from", suffix = c(".total", ".overlap")) %>% 
    mutate(n.overlap = replace_na(n.overlap, 0)) %>%
    mutate(to = from, value = (n.total - n.overlap)/2, value2 = value) %>% 
    dplyr::select(from, to, value, value2)  # use dplyr::select
  
  group_genes <- overlap_genes %>%
    group_by(to) %>%
    summarise(group = paste(from, collapse = "_"), .groups = "drop") %>%  # 添加.groups参数避免警告
    mutate(group = str_remove_all(group, paste0(name1, "\\.|", name2, "\\."))) %>%
    arrange(group)
  
  group <- structure(c(rep("dirs", 4), as.character(group_genes$group)), 
                     names = c(as.character(names), group_genes$to)) %>% 
    factor(c(unique(as.character(group_genes$group)), "dirs"))
  
  group_mid <- split(group, group) %>%
    lapply(function(x) {x[ceiling(length(x)/2)]}) %>%
    lapply(names) 
  
  if(is.null(circo_pal)) {
    grid.col.ds <- setNames(c("#f7970c", "#4a4e4d", "brown1", "#0275d8"), names)
  } else {
    grid.col.ds <- unlist(circo_pal[as.character(names)])
  }
  grid.col <- c(grid.col.ds, setNames(rep("grey", nrow(group_genes)), group_genes$to))
  
  circo_df <- rbind(overlap_genes, count_ds) 
  if(show_selected_genes) {
    genes_display <- intersect(group_genes$to, c(gene_list, phago_sig))
  } else { 
    genes_display <- group_genes$to 
  }
  
  par(bg="white")
  
  # create the out put directory if not exist
  if(!dir.exists("figures/circos")) {
    dir.create("figures/circos", recursive = TRUE)
  }
  if(!dir.exists("Source_Data")) {
    dir.create("Source_Data", recursive = TRUE)
  }
  
  # Correct the file path generation logic
  pdf_path <- paste0("figures/circos/Fig.", Fig_num, "circos_", name1, "_", name2, "_", suffix, ".pdf")
  pdf_path <- str_replace_all(pdf_path, "\\\n| ", "_") %>% 
    str_replace_all("_+", "_") %>% 
    str_remove_all(",")
  
  dev.copy(pdf, pdf_path, height = size, width = size)
  
  circos.clear()
  circos.par(clock.wise = FALSE, start.degree = degree, message = FALSE)
  chordDiagram(circo_df, grid.col = grid.col,
               link.visible = circo_df[[1]] != circo_df[[2]],
               self.link = 2, 
               big.gap = big_gap, small.gap = small_gap,
               order = names(group),
               group = group,
               link.sort = TRUE, 
               transparency = 0.1,
               annotationTrack = c("grid"), 
               preAllocateTracks = list(track.height = max(strwidth(names(group)))))
  
  for(si in genes_display) {
    xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
    ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    if(si %in% phago_sig) {col = "#00bf46"; face = 4}
    else if(si %in% intersect(subset(group_genes, group == "up_up")$to, gene_list)) {col = "brown1"; face = 4}
    else if(si %in% intersect(subset(group_genes, group == "down_down")$to, gene_list)) {col = "dodgerblue3"; face = 4}
    else if (si %in% gene_list) {col = "black"; face = 4}
    else {col = "grey50"; face = 3}
    
    circos.text(mean(xlim), ylim[1], si, sector.index = si, track.index = 1, 
                facing = "clockwise", adj = c(0, 0.5), col = col,
                cex = gene_cex, niceFacing = TRUE, font = face)
  }
  
  # Check is overlap_perm_test exist
  if(!is.null(overlap_perm_test)) {
    perm_test_res <- overlap_perm_test$p_overlap_perm %>%
      pivot_longer(cols = 2:3, names_to = "direction.set1", values_to = "p") %>%
      mutate(dir = paste(str_extract(direction.set1, "up|down"),
                         str_extract(direction.set2, "up|down"), sep = "_"))
    
    p_facing <- list(up_up = "outside", up_down = conflict_facing,
                     down_up = conflict_facing, down_down = "inside")
    nudge_x <- setNames(nudge_x, names)
    
    for(x in names(group_mid)[-length(group_mid)]) {
      si <- group_mid[[x]]
      xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
      ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
      
      p <- subset(perm_test_res, dir == x)$p
      label <- paste0("p = ", round(ifelse(p > 0.5, 1 - p, p), 4), 
                      ifelse(p > 0.5, " *", ""))
      if(p == 0 | p == 1) {
        label <- substitute(
          paste(bold(x)^bold(y), bold(z)), list(x = "p < 10", y = "-5", z = ifelse(p > 0.5, " *", ""))) %>% as.expression()
      }
      if(sum(group == x) < 4 & break_p_label) {
        label <- str_replace(label, "= ",  "=\n")
      }
      
      circos.text(mean(xlim), 0.55, label, sector.index = si, track.index = 1,
                  facing = p_facing[[x]], cex = p_cex, niceFacing = FALSE, font = 2)
    }
  }
  
  for(si in unique(circo_df$from)) {
    si_dir <- paste0(str_to_title(str_extract(si, "up|down")), " in")
    si_name <- ifelse(str_detect(si, name1), name1, name2)
    
    xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
    ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    
    if(str_detect(si_name, "cKO|\\+|\\-")) {
      si_sup <- str_extract(si_name, "cKO|\\+|\\-")
      si_base <- str_remove(si_name, "cKO|\\+|\\-")
      y = c(0.37, 0.26)
      label_dir <- substitute(bold(dir), list(dir = si_dir)) 
      label <- c(label_dir, substitute(bolditalic(x)^bold(y), list(x = si_base, y = si_sup))) %>%
        as.expression()
    } else { 
      y = 0.3; label <- paste0(si_dir, "\n", si_name) 
    }
    circos.text(
      mean(xlim) + nudge_x[si], y, label,
      sector.index = si, track.index = 1, col = grid.col.ds[si],
      facing = "inside", font = 2,
      cex=num_cex, niceFacing = niceFacing
    )
    
    circos.axis(h = 0, major.at = seq(0, 3000, step),
                labels.cex = 0.5, labels.facing = "outside",
                sector.index = si, track.index = 1)
  }
  
  dev.off()
  
  overlap_genes %>% select(-value2) %>%
    pivot_wider(names_from = "from", values_from = "value") %>%
    write.csv(paste0("Source_Data/Fig.", Fig_num, "_circos_", name1, "_", name2, "_", suffix, ".csv"), row.names = FALSE)
}



# genes to highlight

gene_to_label <- c("Hspa8", "Pawr", "Aldh1a1", "Txnip", "Ptgs2","Kdm6b",
                    "Mapk3", "Nfe2l1", "Cd36", "Prdx6b", "Edn1", "Slc7a11",
                    "Cyp1b1", "Park7", "Prdx4", "Sirt6")


# Top DEGs in susceptible vs control

sus <- read.csv("vHPC_log2CPM_CvsS_full_results_mapped.csv", row.names = 1)

thr <- log2(1.3)     

sus$negLogFDR <- -log10(sus$adj.P.Val)
sus$regulation <- "Not significant"
sus$regulation[sus$adj.P.Val < 0.05 & sus$logFC > thr] <- "Upregulated"
sus$regulation[sus$adj.P.Val < 0.05 & sus$logFC < -thr] <- "Downregulated"
sus <- sus %>% arrange(desc(negLogFDR)) %>% distinct(SYMBOL, .keep_all = TRUE)

sus_sig <- sus %>%arrange(desc(negLogFDR)) %>%slice(1:1000) %>%
  mutate(direction = ifelse(logFC > 0, "up", "down"))




# Top 100 DEGs in resilient vs control

res <- read.csv("vHPC_log2CPM_CvsR_full_results_mapped.csv", row.names = 1)

thr <- log2(1.3)     

res$negLogFDR <- -log10(res$adj.P.Val)
res$regulation <- "Not significant"
res$regulation[res$adj.P.Val < 0.05 & res$logFC > thr] <- "Upregulated"
res$regulation[res$adj.P.Val < 0.05 & res$logFC < -thr] <- "Downregulated"


res <- res %>% arrange(desc(negLogFDR)) %>% distinct(SYMBOL, .keep_all = TRUE)

res_sig <- res %>%arrange(desc(negLogFDR)) %>%slice(1:1000) %>%
  mutate(direction = ifelse(logFC > 0, "up", "down"))

# Premutation test on the overlapping gens------

overlap_sus_res_sig <- perm_test(
  set1 = sus_sig, set2 = res_sig,
  all_genes1 = sus$SYMBOL,
  all_genes2 = res$SYMBOL,
  name1 = "sus", name2 = "res",
  n_perm = 10000
)

write.csv(overlap_sus_res_sig$p_overlap_perm,
          "Source_Data/Fig.2D_perm_test.csv")


# Circos plot comparison of up- and down-regulated DEGs ---------------


circo_overlap(res1 = sus_sig, res2 = res_sig, Fig_num = "2D",
              name1 = "sus", name2 = "res", 
              suffix = "sus_res_sub", size = 7.5, gene_width = 25,
              num_cex = 0.9, gene_cex = 0.75, gene_width_short = 3, 
              big_gap = 5, small_gap = 0.4, 
              gene_list = c(gene_to_label, "Gpx4"),
              degree = 100, step = 300, show_selected_genes = TRUE, circo_pal = pal_text,
              overlap_perm_test = overlap_sus_res_sig, 
              nudge_x = c(0, 0, 50, -50), niceFacing = FALSE)


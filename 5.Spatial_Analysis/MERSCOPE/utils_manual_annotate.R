#===============================================================================|
# ANNOTATION FUNCTION ---------
#===============================================================================|

annotate_clusters_manual <- function(
    deseq_results,
    marker_panels = standard_markers(),
    score_method = c("weighted", "count"),
    min_log2fc = 0,
    max_padj = 0.05,
    top_n = 3,
    output_prefix = "annotation",
    ambiguity_threshold = 0.7  # ratio of top score to second score
) {
  
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(patchwork)
  
  score_method <- match.arg(score_method)
  
  cat("Starting cluster annotation...\n")
  
  #===============================================================================|
  # 1. PREPARE DATA --------
  #===============================================================================|
  
  df <- deseq_results %>%
    mutate(
      gene_upper = toupper(gene),
      padj_safe = ifelse(is.na(padj) | padj <= 0, 1e-300, padj),
      neglog10p = -log10(padj_safe),
      is_up = log2FoldChange > min_log2fc & padj < max_padj
    )
  
  up_genes <- df %>% filter(is_up)
  clusters <- sort(unique(up_genes$cluster))
  
  cat(sprintf("Analyzing %d clusters with %d upregulated genes\n", 
              length(clusters), nrow(up_genes)))
  
  #===============================================================================|
  # 2. SCORE EACH CLUSTER --------
  #===============================================================================|
  
  score_matrix <- matrix(
    0, 
    nrow = length(clusters), 
    ncol = length(marker_panels),
    dimnames = list(as.character(clusters), names(marker_panels))
  )
  
  summary_list <- list()
  contributors_list <- list()
  
  for(clust in clusters) {
    
    cluster_genes <- up_genes %>% filter(cluster == clust)
    
    # Calculate scores for each cell type
    scores <- numeric(length(marker_panels))
    names(scores) <- names(marker_panels)
    
    for(celltype in names(marker_panels)) {
      markers_upper <- toupper(marker_panels[[celltype]])
      matching <- cluster_genes %>% filter(gene_upper %in% markers_upper)
      
      if(nrow(matching) == 0) {
        scores[celltype] <- 0
      } else if(score_method == "weighted") {
        #For a given cluster × cell type pair: 
        #find which of that cell type's known marker genes are 
        #upregulated in that cluster (from DESeq2), 
        #then sum up log2FC × -log10(padj) across those matching markers. 
        scores[celltype] <- sum(matching$log2FoldChange * matching$neglog10p, na.rm = TRUE)
      } else {
        # Or simply count the number of matching markers
        scores[celltype] <- nrow(matching)
      }
      
      # Store contributors for this celltype
      if(nrow(matching) > 0) {
        contrib <- matching %>%
          arrange(desc(neglog10p)) %>%
          mutate(
            cluster = clust,
            celltype = celltype,
            score = log2FoldChange * neglog10p
          ) %>%
          dplyr::select(cluster, celltype, gene, log2FoldChange, padj, score) %>%
          head(5)
        
        contributors_list[[paste(clust, celltype, sep = "_")]] <- contrib
      }
    }
    
    score_matrix[as.character(clust), ] <- scores
    
    # Top predictions
    top_idx <- order(scores, decreasing = TRUE)[1:min(top_n, length(scores))]
    
    #===============================================================================|
    # Check for ambiguity ---------
    #===============================================================================|
    
    top1_score <- scores[top_idx[1]]
    top2_score <- scores[top_idx[2]]
    
    # Calculate ratio of top score to second score (Label as ambiguous if ratio < threshold)
    score_ratio <- if(top2_score > 0) top1_score / top2_score else Inf
    
    is_ambiguous <- score_ratio < (1 / ambiguity_threshold) | 
      (top1_score > 0 & top2_score > 0 & score_ratio < 1.5)
    
    celltype1 <- names(marker_panels)[top_idx[1]]
    
    summary_row <- tibble(
      cluster = clust,
      n_up_genes = nrow(cluster_genes),
      n_marker_hits = sum(cluster_genes$gene_upper %in% 
                            unique(toupper(unlist(marker_panels)))),
      is_ambiguous = is_ambiguous,
      score_ratio = score_ratio
    )
    
    for(i in 1:top_n) {
      if(i <= length(top_idx)) {
        idx <- top_idx[i]
        summary_row[[paste0("rank", i, "_celltype")]] <- names(marker_panels)[idx]
        summary_row[[paste0("rank", i, "_score")]] <- scores[idx]
      } else {
        summary_row[[paste0("rank", i, "_celltype")]] <- NA_character_
        summary_row[[paste0("rank", i, "_score")]] <- NA_real_
      }
    }
    
    summary_list[[as.character(clust)]] <- summary_row
  }
  
  summary_df <- bind_rows(summary_list) %>% arrange(cluster)
  contributors_df <- bind_rows(contributors_list)
  
  #===============================================================================|
  # 3. NORMALIZE SCORES -------
  #===============================================================================|
  
  score_matrix_norm <- t(scale(t(score_matrix)))
  score_matrix_norm[is.na(score_matrix_norm)] <- 0
  
  #===============================================================================|
  # 3B. REFINED AMBIGUITY DETECTION using normalized scores
  #===============================================================================|
  
  for(i in 1:nrow(summary_df)) {
    clust <- summary_df$cluster[i]
    
    # Get z-scores for this cluster
    cluster_zscores <- score_matrix_norm[as.character(clust), ]
    
    # Count how many cell types have high z-scores (> 0.5)
    n_high_scores <- sum(cluster_zscores > 0.5)
    
    # Check if top 2 z-scores are similar
    top_zscores <- sort(cluster_zscores, decreasing = TRUE)[1:2]
    zscore_diff <- top_zscores[1] - top_zscores[2]
    
    # Mark as ambiguous if:
    # 1. Multiple cell types have z-score > 0.5, OR
    # 2. Top 2 z-scores are very close (diff < 0.5)
    summary_df$is_ambiguous[i] <- n_high_scores >= 2 | zscore_diff < 0.5
    summary_df$n_high_scores[i] <- n_high_scores
    summary_df$zscore_diff[i] <- zscore_diff
  }
  
  #===============================================================================|
  # 4. ADD "-like" SUFFIX FOR AMBIGUOUS ANNOTATIONS --------
  #===============================================================================|
  
  summary_df <- summary_df %>%
    mutate(
      rank1_celltype_final = ifelse(
        is_ambiguous & !is.na(rank1_celltype) & rank1_score > 0,
        paste0(rank1_celltype, "-like"),
        rank1_celltype
      ),
      confidence = case_when(
        is.na(rank1_celltype) | rank1_score == 0 ~ "None", # no marker hit
        !is_ambiguous & zscore_diff > 1.0 ~ "High",
        !is_ambiguous ~ "Medium",
        TRUE ~ "Low (ambiguous)"
      )
    )
  
  #===============================================================================|
  # 5. ORDER CLUSTERS BY CELL TYPE -------------
  #===============================================================================|
  
  cluster_order <- summary_df %>%
    arrange(rank1_celltype_final, desc(rank1_score)) %>%
    pull(cluster)
  
  #===============================================================================|
  # 6. CREATE FIGURE -------------
  #===============================================================================|
  
  cat("Creating visualizations...\n")
  
  # Panel A: Barplot with confidence indicated by alpha
  p_bar <- ggplot(summary_df, 
                  aes(x = factor(cluster, levels = cluster_order), 
                      y = rank1_score, 
                      fill = rank1_celltype_final,
                      alpha = confidence)) +
    geom_col(width = 0.7) +
    scale_alpha_manual(
      values = c("High" = 1.0, "Medium" = 0.8, "Low (ambiguous)" = 0.6, "None" = 0.3),
      name = "Confidence"
    ) +
    coord_flip() +
    scale_fill_brewer(palette = "Set3", name = "Cell Type") +
    labs(
      x = "Cluster",
      y = ifelse(score_method == "weighted", 
                 "Score (Sum of log2FC x -log10 padj)",
                 "# Markers"),
      title = "A. Top Predicted Cell Type"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      plot.title = element_text(face = "bold", size = 12)
    )
  
  # Panel B: Heatmap
  score_long <- as.data.frame(score_matrix_norm) %>%
    tibble::rownames_to_column("cluster") %>%
    pivot_longer(-cluster, names_to = "celltype", values_to = "zscore") %>%
    mutate(
      cluster = factor(cluster, levels = cluster_order),
      celltype = factor(celltype, levels = colnames(score_matrix_norm))
    )
  
  p_heat <- ggplot(score_long, aes(x = celltype, y = cluster, fill = zscore)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0, 
      name = "z-score",
      limits = c(-3, 3),
      oob = scales::squish
    ) +
    labs(
      x = "Cell Type",
      y = "Cluster",
      title = "B. Cell Type Enrichment (Row z-score)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "right",
      panel.grid = element_blank()
    )
  
  # Panel C: Dotplot
  p_dot <- ggplot(score_long %>% filter(abs(zscore) > 0.5), 
                  aes(x = cluster, y = celltype)) +
    geom_point(aes(size = abs(zscore), color = zscore)) +
    scale_size_continuous(range = c(2, 10), name = "abs(z-score)") +
    scale_color_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0, 
      name = "z-score"
    ) +
    labs(
      x = "Cluster (ordered by annotation)",
      y = "Cell Type",
      title = "C. Significant Enrichments (|z| > 0.5)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(face = "bold", size = 12),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
  
  # Panel D: Summary with confidence
  summary_stats <- summary_df %>%
    group_by(rank1_celltype_final, confidence) %>%
    summarise(n_clusters = n(), .groups = "drop") %>%
    arrange(desc(n_clusters)) %>%
    head(10)
  
  p_stats <- ggplot(summary_stats, 
                    aes(x = reorder(rank1_celltype_final, n_clusters), 
                        y = n_clusters,
                        fill = confidence)) +
    geom_col(width = 0.7) +
    scale_fill_manual(
      values = c("High" = "#2ecc71", "Medium" = "#f39c12", 
                 "Low (ambiguous)" = "#e74c3c", "None" = "#95a5a6"),
      name = "Confidence"
    ) +
    geom_text(aes(label = n_clusters), hjust = -0.2, size = 3.5) +
    coord_flip() +
    labs(
      x = "Cell Type",
      y = "Number of Clusters",
      title = "D. Cluster Distribution by Cell Type"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    ) +
    ylim(0, max(summary_stats$n_clusters) * 1.1)
  
  # Combine panels
  combined_plot <- (p_bar | p_heat) / 
    (p_dot + plot_layout(widths = c(2)) | p_stats + plot_layout(widths = c(1))) +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(
      title = "Cluster Cell Type Annotation Summary",
      subtitle = sprintf("Method: %s | Thresholds: log2FC > %.2f, padj < %.2f ",
                         score_method, min_log2fc, max_padj),
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5)
      )
    )
  
  #===============================================================================|
  # 7. CREATE MAPPINGS
  #===============================================================================|
  
  celltype_to_clusters <- summary_df %>%
    group_by(rank1_celltype_final, confidence) %>%
    summarise(
      n_clusters = n(),
      clusters = paste(cluster, collapse = ", "),
      avg_score = mean(rank1_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_clusters), desc(avg_score))
  
  cluster_to_celltype <- summary_df %>%
    select(cluster, 
           celltype = rank1_celltype_final, 
           celltype_base = rank1_celltype,
           score = rank1_score,
           confidence,
           is_ambiguous,
           n_high_scores,
           zscore_diff) %>%
    arrange(celltype, desc(score))
  # if confidence is None (no marker hit)
  cluster_to_celltype = mutate(cluster_to_celltype, celltype = case_when(confidence=="None" ~ "Unassigned", TRUE ~ celltype))
  
  #===============================================================================|
  # 8. SAVE OUTPUTS ----------
  #===============================================================================|
  
  ggsave(
    paste0(output_prefix, "_combined_figure.pdf"),
    combined_plot,
    width = 18,
    height = 12,
    units = "in"
  )
  cat(sprintf("Saved: %s_combined_figure.pdf\n", output_prefix))
  
  result <- list(
    summary = summary_df,
    contributors = contributors_df,
    score_matrix_raw = score_matrix,
    score_matrix_norm = score_matrix_norm,
    celltype_to_clusters = celltype_to_clusters,
    cluster_to_celltype = cluster_to_celltype,
    cluster_order = cluster_order,
    combined_plot = combined_plot,
    parameters = list(
      marker_panels = marker_panels,
      score_method = score_method,
      min_log2fc = min_log2fc,
      max_padj = max_padj,
      ambiguity_threshold = ambiguity_threshold,
      n_clusters = length(clusters),
      n_celltypes = length(marker_panels)
    )
  )
  
  class(result) <- c("cluster_annotation", "list")
  
  saveRDS(result, paste0(output_prefix, "_results.rds"))
  cat(sprintf("Saved: %s_results.rds\n", output_prefix))
  
  write.csv(
    summary_df %>% dplyr::select(cluster, rank1_celltype_final, rank1_score, 
                                 confidence, is_ambiguous, n_high_scores,
                                 rank2_celltype, rank2_score, 
                                 n_up_genes, n_marker_hits),
    paste0(output_prefix, "_summary.csv"),
    row.names = FALSE
  )
  
  write.csv(celltype_to_clusters,
            paste0(output_prefix, "_celltype_to_clusters.csv"),
            row.names = FALSE)
  
  write.csv(cluster_to_celltype,
            paste0(output_prefix, "_cluster_to_celltype.csv"),
            row.names = FALSE)
  
  cat("\n=== ANNOTATION SUMMARY ===\n")
  print(summary_df %>% 
          dplyr::select(cluster, rank1_celltype_final, confidence, rank1_score))
  
  cat("\n=== AMBIGUOUS CLUSTERS ===\n")
  ambiguous <- summary_df %>% filter(is_ambiguous)
  if(nrow(ambiguous) > 0) {
    print(ambiguous %>% select(cluster, rank1_celltype_final, rank2_celltype, 
                               zscore_diff, n_high_scores))
  } else {
    cat("No ambiguous clusters detected.\n")
  }
  
  return(invisible(result))
}
#===============================================================================|==|
# HELPER FUNCTIONS ----------
#===============================================================================|==|

get_celltype_mapping <- function(annotation_result) {
  annotation_result$celltype_to_clusters
}

get_cluster_mapping <- function(annotation_result) {
  annotation_result$cluster_to_celltype
}

add_annotations_to_seurat <- function(seurat_obj, annotation_result, 
                                      cluster_col = "seurat_clusters") {
  
  mapping <- annotation_result$cluster_to_celltype
  
  celltype_map <- setNames(mapping$celltype, mapping$cluster)
  
  seurat_obj$de.novo.celltype <- plyr::mapvalues(
    seurat_obj@meta.data[[cluster_col]],
    from = str_remove(names(celltype_map),"g"),
    to = as.character(celltype_map),
    warn_missing = FALSE
  )
  
  score_map <- setNames(mapping$score, mapping$cluster)
  seurat_obj$de.novo.celltype_score <- plyr::mapvalues(
    seurat_obj@meta.data[[cluster_col]],
    from = str_remove(names(score_map),"g"),
    to = score_map,
    warn_missing = FALSE
  )
  
  return(seurat_obj)
}

get_markers <- function(annotation_result, cluster_id, celltype = NULL) {
  if(is.null(celltype)) {
    celltype <- annotation_result$summary %>%
      filter(cluster == cluster_id) %>%
      pull(rank1_celltype)
  }
  
  annotation_result$contributors %>%
    filter(cluster == cluster_id, celltype == celltype) %>%
    arrange(desc(score))
}

compare_celltype <- function(annotation_result, celltype) {
  scores <- annotation_result$score_matrix_norm[, celltype, drop = FALSE]
  
  data.frame(
    cluster = rownames(scores),
    zscore = scores[, 1]
  ) %>%
    arrange(desc(zscore))
}

load_annotation <- function(rds_file) {
  result <- readRDS(rds_file)
  class(result) <- c("cluster_annotation", "list")
  return(result)
}

# Standard markers 
gold_standard_markers <- function() {
  list(
    Astrocyte         = c("AQP4", "SLC1A3", "ALDH1L1","GJA1","FGFR3"),
    Oligodendrocyte   = c("MOG", "PLP1", "MBP", "MOBP","OPALIN","OLIG1","OLIG2"),
    OPC               = c("PDGFRA", "CSPG4"),
    Microglia         = c("P2RY12", "TMEM119","HEXB","CX3CR1"),
    Macrophage        = c("CD163", "MRC1", "LYVE1"),
    Monocyte          = c("CD14","LYZ"),
    Endothelial       = c("CLDN5", "PECAM1", "CDH5","VWF"),
    Pericyte          = c("PDGFRB", "RGS5"),
    Smooth_Muscle     = c("ACTA2", "MYH11", "TAGLN"),
    Fibroblast        = c("COL1A1", "DCN", "LUM"),
    Excitatory_Neuron = c("SLC17A7", "CAMK2A","RBFOX3"),
    Inhibitory_Neuron= c("GAD1", "GAD2","PVALB", "SST", "VIP"),
    T_Cell            = c("CD3D", "CD3E","CCR7","CD4")
  )
}

create_batch_colors <- function() {
  c("Phase2-1-region_R1" = "#ff6d6f", "Phase2-1-region_R2" = "#ed6668", "Phase2-1-region_R3" = "#b06060", "Phase2-2-region_R1" = "#88caff", "Phase2-2-region_R2" = "#7aa9d0", "Phase2-2-region_R3" = "#6b889e", "Phase2-3-region_R1" = "#9df89a", "Phase2-3-region_R2" = "#8fd98d", "Phase2-3-region_R3" = "#82ba80", "Phase2-3-region_R4" = "#749b72", "Phase2-4-region_R1" = "#e39eed", "Phase2-4-region_R2" = "#c890d0", "Phase2-4-region_R3" = "#ad82b3", "Phase2-4-region_R4" = "#927496", "Phase2-5-region_R1" = "#ffcc55", "Phase2-5-region_R3" = "#ffb555", "Phase2-5-region_R4" = "#e89e55", "Phase2-5-region_R5" = "#bb8855")
  
}
create_celltype_colors <- function() {

  colors <- c(
    "Ambiguous" = "gray50",
    "Oligodendrocytes" = "#789cff",
    "Oligodendrocyte" = "#789cff",
    "Oligodendrocyte-like" = "#B5C6F3",
    "Astrocyte" = "#00CED1",
    "Astrocyte-like" = "#7FFFD4",
    "Microglia" = "#8B00FF",
    "Microglia-like" = "#DA70D6",
    "OPCs" = "#FF1493",
    "OPC" = "#FF1493",
    "OPC-like" = "#FF69B4",
    "Macrophages" = "blue",
    "Macrophage" = "blue",
    "Macrophage-like" = "#79CCFF",
    "Monocytes" = "black",
    "Monocyte" = "black",
    "Monocyte-like" = "#C0C0C0",
    "Endothelial" = "#fa003f",
    "Endothelial-like" = "#F56C7A",
    "Pericytes" = "#6F05E9",
    "Pericyte" = "#6F05E9",
    "Pericyte-like" = "#C5A8E6",
    "SMC" = "#FFD4C1",
    "Smooth_Muscle" = "#FFD4C1",
    "Smooth_Muscle-like" = "#EE6123",
    "Fibroblast" = "#006400",
    "Fibroblast-like" = "#32CD32",
    "Excitatory Neurons" =  "#FFD000",
    "Excitatory_Neuron" = "#FFD000",
    "Excitatory_Neuron-like" = "#fff1b2",
    "Inhibitory Neurons"= "#9ACD32",
    "Inhibitory_Neuron" = "#9ACD32",
    "Inhibitory_Neuron-like" = "#ADFF2F",
    "CD8+ T Cells" = "#3E0085",
    "T_Cell" = "#3E0085",
    "T_Cell-like" = "#736089"
  )
  
  return(colors)
}


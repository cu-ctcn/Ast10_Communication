### Pipeline for White Matter/Grey Matter Annotations ------ #####################################
library(Seurat)
library(tidyverse)
library(scCustomize)
library(patchwork)

celltype_pal = scCustomize_Palette(num_groups = 36) #c("#0072B2", "#E69F00", "#009E73", "#F0E442", "#CC79A7", "#56B4E9", "#D55E00", "#000000", "#006064", "#B36C8D", "#A8A8A8", "#C85122", "#36A96E", "#7F3C8D", "#FFC20A","#7F3C8D", "#FFC20A")

annotate_wm_gm_pipeline <- function(seurat_object, region, output_dir,figure_dir,  annotation_column = "predicted.cell.type",
                                    oligodendrocyte_celltype = "Oligodendrocyte", neuronal_celltype =  "Excitatory_Neuron",
                                    wm_markers, gm_markers,
                                    fov_id = "fov1", niches_k = 2, neighbors_k = 50, assay = "Vizgen",
                                    mypal = c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B")) {
  

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  sobj <- seurat_object
  sobj@meta.data =  sobj@meta.data%>%select(-contains("niche"))
  
  assay_list = sobj@assays%>%names
  DefaultAssay(sobj) = "Vizgen"
  if ("niche"%in%assay_list){
    sobj[["niche"]]= NULL
    }
  
  sobj <- BuildNicheAssay(object = sobj, fov = fov_id, group.by = annotation_column,
                          niches.k = niches_k, neighbors.k = neighbors_k)
  
  
  # Automatically determine WM/GM based on oligodendrocyte fraction
  message("Determining WM/GM niches...")
  niche_fraction <- sobj@meta.data %>%
    group_by(niches) %>%
    summarise(
      total_cells_in_niche = n(),
      oligodendrocyte_count = sum(.data[[annotation_column]] %in% oligodendrocyte_celltype, na.rm = TRUE),
      oligodendrocyte_fraction = oligodendrocyte_count / total_cells_in_niche,
      exc_neuron_count = sum(.data[[annotation_column]] %in% neuronal_celltype, na.rm = TRUE),
      exc_neuron_fraction = exc_neuron_count / total_cells_in_niche,
      .groups = 'drop'
    ) %>%
    arrange(desc(oligodendrocyte_fraction))
    
 
  max_oli_frac = max(niche_fraction$oligodendrocyte_fraction)
  max_neu_frac = max(niche_fraction$exc_neuron_fraction)
  wm_niche_id = filter(niche_fraction)%>%slice_max(oligodendrocyte_fraction)%>%pull(niches)
  # assume grey matter if: got the highest fraction of exc neurons, OR (not GM AND Exc fraction > 5%)
  gm_niche_id = c(filter(niche_fraction)%>%slice_max(exc_neuron_fraction)%>%pull(niches), filter(niche_fraction,(!niches==wm_niche_id) & exc_neuron_fraction >=0.05)%>%pull(niches))%>%unique
  
  if (is.na(wm_niche_id) || sum(is.na(gm_niche_id)) ){
    sobj$predicted.niche <- "Unassigned"
  } else {
    message(paste0("WM niche: ", wm_niche_id, " (Oli%: ", round(niche_fraction$oligodendrocyte_fraction[1]*100, 1), ")"))
    message(paste0("GM niche: ", gm_niche_id, " (Oli%: ", round(niche_fraction$oligodendrocyte_fraction[2]*100, 1), ")"))
    
    sobj$predicted.niche <- NA
    sobj$predicted.niche[sobj$niches == wm_niche_id] <- "White Matter"
    sobj$predicted.niche[sobj$niches %in% gm_niche_id] <- "Grey Matter"
    sobj$predicted.niche <- factor(sobj$predicted.niche, levels = c("Grey Matter","White Matter"))
  }
  
  # 1. Output visualization of the tissue
  message("Generating tissue visualization...")
  celltype.plot <- ImageDimPlot(sobj,fov = fov_id, group.by = annotation_column, size = 0.75, cols = "polychrome",
                                dark.background = FALSE) + ggtitle("Cell Type Distribution")
  
  niche.plot <- ImageDimPlot(sobj,fov = fov_id, group.by = "predicted.niche", size = 0.75, dark.background = FALSE) +
    ggtitle("Niches Identified") +
    scale_fill_manual(values = mypal) 
  
  combined_plot <- celltype.plot | niche.plot
  ggsave(file.path(figure_dir, sprintf("%s_tissue_and_niche_visualization.pdf", region)), combined_plot, width = 11.5, height = 5.5, dpi = 300)
  
  
  # 2. Barplot: Fraction of cell types in WM/GM
  message("Generating cell type fraction plot...")
  # Calculate cell counts and fractions
  cell_counts <- sobj@meta.data %>%
    filter(!is.na(predicted.niche)) %>%
    group_by(predicted.niche, .data[[annotation_column]]) %>%
    summarise(count = n(), .groups = 'drop')
  
  total_wm_cells <- sum(cell_counts$count[cell_counts$predicted.niche == "White Matter"], na.rm = TRUE)
  total_gm_cells <- sum(cell_counts$count[cell_counts$predicted.niche == "Grey Matter"], na.rm = TRUE)
  
  fraction_data <- cell_counts %>%
    group_by(predicted.niche) %>%
    mutate(fraction = count / sum(count)) %>%
    ungroup()
  
  bar_plot <- ggplot(fraction_data, aes(x = predicted.niche, y = fraction, fill = .data[[annotation_column]])) +
    geom_bar(stat = "identity", position = "fill") +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "Fraction of Cell Types in White Matter and Grey Matter",
         x = "Tissue Annotation", y = "Fraction of Cells", fill = "Cell Type") +
    theme_minimal() +
    scale_fill_manual(values = celltype_pal)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5)) +
    annotate("text", x = "White Matter", y = 1.05, label = paste("WM Cells:", total_wm_cells), vjust = 0, size = 3) +
    annotate("text", x = "Grey Matter", y = 1.05, label = paste("GM Cells:", total_gm_cells), vjust = 0, size = 3)
  
  ggsave(file.path(figure_dir,  sprintf("%s_cell_type_fraction_bar_plot.pdf", region)), bar_plot, width = 8, height = 6, dpi = 300)
  
  # Tissue plot with WM/GM markers
  message("Generating marker plots...")
 
  all_genes_in_object <- rownames(GetAssayData(object = sobj, assay = assay, slot = "counts"))
  

  wm_markers_to_plot <- wm_markers[wm_markers %in% all_genes_in_object]
  gm_markers_to_plot <- gm_markers[gm_markers %in% all_genes_in_object]
  
 
    DefaultAssay(sobj) = assay
    wm_marker_plots <- lapply(wm_markers_to_plot, plot_spatial_gene_expression, sobj, fov_id, assay)

    gm_marker_plots <- lapply(gm_markers_to_plot,  plot_spatial_gene_expression, sobj, fov_id, assay)
    
    # Combine marker plots
    all_marker_plots <- c(wm_marker_plots, gm_marker_plots)
    if (length(all_marker_plots) > 0) {
      combined_marker_plot <- ggpubr::ggarrange(plotlist = all_marker_plots, ncol = length(all_marker_plots))
      ggsave(file.path(figure_dir, sprintf("%s_wm_gm_marker_tissue_plots.pdf", region)), combined_marker_plot,
             width = 4 * min(length(all_marker_plots), 3.5),
             height = 4 * ceiling(length(all_marker_plots) / min(length(all_marker_plots), 5)),
             dpi = 300)
    } else {
      message("No valid markers found. Skipping.")
    }
  
 
  # 3. Save the annotations 
  annotations_df <- sobj@meta.data %>%
    select(one_of(c(annotation_column, "predicted.niche"))) %>%
    mutate(cell_barcode = rownames(sobj@meta.data)) %>%
    filter(!is.na(predicted.niche))
  
  
  annotations_path <- file.path(output_dir, "wm_gm_annotations.csv")
  write_csv(annotations_df, annotations_path)
  return(sobj)
}

#==========================================================================|
## Reads in annotation CSV and inserts into Seurat object
insert_annotations_into_seurat <- function(seurat_object, annotations_csv_path,
                                           barcode_column_csv, annotation_column_csv,
                                           new_metadata_column_name) {
  
  annotations_df <- read_csv(annotations_csv_path)
  
  # Ensure the barcode column exists in the CSV
  if (!barcode_column_csv %in% colnames(annotations_df)) {
    stop(paste0("Barcode column '", barcode_column_csv, "' not found."))
  }
  # Ensure the annotation column exists in the CSV
  if (!annotation_column_csv %in% colnames(annotations_df)) {
    stop(paste0("Annotation column '", annotation_column_csv, "' not found."))
  }
  
  annotations_df <- annotations_df %>%
    column_to_rownames(var = barcode_column_csv)
  
  # Get cell names
  seurat_cell_names <- colnames(seurat_object)
  
  new_annotations <- rep(NA, length(seurat_cell_names))
  names(new_annotations) <- seurat_cell_names
  
  # Merge annotations
  common_cells <- intersect(seurat_cell_names, rownames(annotations_df))
  
  if (length(common_cells) == 0) {
    warning("No matching barcodes found. No annotations added.")
  } else {
    new_annotations[common_cells] <- annotations_df[common_cells, annotation_column_csv]
  }

  # Add the new metadata column to the Seurat object
  seurat_object <- AddMetaData(object = seurat_object, metadata = new_annotations, col.name = new_metadata_column_name)
  
  return(seurat_object)
}




#' Custom function to plot gene expression on spatial coordinates
plot_spatial_gene_expression <- function( gene_name, seurat_object,fov_id, assay , plot_slot = "data", alpha_range = c(0.1, 1), point_size = 0.25) {
  DefaultAssay(seurat_object) = assay 
  if (!("data" %in% Layers(seurat_object@assays[[DefaultAssay(seurat_object)]]))){
    warning("Data not normalized, normalizing...")
    seurat_object <- NormalizeData(seurat_object)
  }
  # Check if the gene is present in the specified slot
  assay_data <- FetchData(object = seurat_object, assay = assay,vars = gene_name, slot = plot_slot)
  
  if (is.null(assay_data)) {
    warning(paste0("Could not retrieve expression for ", gene_name, ". Skipping."))
    return(NULL)
  }
  
  if (!(gene_name %in% rownames(seurat_object))) {
    warning(paste0("Gene '", gene_name, "' not found. Skipping."))
    return(NULL)
  }
  
  # Extract spatial coordinates
  spatial_coords <- tryCatch({
    GetTissueCoordinates(seurat_object@images[[fov_id]])
  }, error = function(e) {
    message(paste0("Error getting coordinates for FOV '", fov_id, "': ", e$message))
    return(NULL)
  })
  
  if (is.null(spatial_coords) || nrow(spatial_coords) == 0) {
    warning(paste0("No coordinates for FOV '", fov_id, "'. Skipping ", gene_name, "."))
    return(NULL)
  }
  
  # Ensure cell names match between expression data and spatial coordinates
  common_cells <- intersect(rownames(assay_data), spatial_coords$cell)
  if (length(common_cells) == 0) {
    warning(paste0("No common cells for FOV '", fov_id, "'. Skipping ", gene_name, "."))
    return(NULL)
  }
  
  # Subset data to common cells
  
  plot_data <- assay_data%>%rownames_to_column("cell")%>%inner_join(spatial_coords, by = "cell")

  # Create the plot
  p <- ggplot(plot_data, aes(x = x, y  = y , color = .data[[gene_name]])) +
    geom_point(size = point_size, alpha = alpha_range[2]) + 
    scale_color_viridis_c(option = "plasma", name = "Expression") + 
    labs(title = paste0(gene_name, " Expression")) +
    theme_void() + 
    coord_fixed() 
  
  return(p)
}

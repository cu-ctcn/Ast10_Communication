rm(list = ls())
#===============================================================================|
# NICHE ANNOTATION FOR INTEGRATED OBJECT ----------
#===============================================================================|

library(Seurat)
library(tidyverse)
library(scCustomize)
library(patchwork)

setwd("Code/MERSCOPE_Analysis/WM_GM_Annotation")
source("utils_niche_annotation_step.R")

niche_colors <- c(
  "White Matter" = "#442288",
  "Grey Matter" = "#6CA2EA",
  "Unassigned" = "#CCCCCC"
)
#===============================================================================|
# SETUP ----------
#===============================================================================|

data_sub_path = "f_results_npc=100"
data_path = "/path/to/integration_dir"
data <- readRDS(file.path(data_path, data_sub_path, "integrated.rds"))

data@meta.data = mutate(data@meta.data, wm.celltype = case_when(de.novo.celltype %in% c("Oligoendrocyte-like") ~ "Oligodendrocyte", TRUE ~ de.novo.celltype ))

assay_list = data@assays%>%names

DefaultAssay(data) = "Vizgen"
if ("niche"%in%assay_list){data[["niche"]]= NULL}
output_dir <- file.path("output", data_sub_path,"niche_annotation")
figure_dir <- output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)


annotation_column <- "de.novo.celltype" 

# Use oligodendrocyte for niche determination, since it is most abundant in WM 
reference_celltype <- c("Oligodendrocyte")

# Markers for visualization
wm_markers <- c("MOG", "MBP","MOBP", "PLP1","OPALIN")
gm_markers <- c("SLC17A7", "RBFOX3","RASGRF2","RORB")

#===============================================================================|
# PROCESS EACH SAMPLE SEPARATELY ----------
#===============================================================================|

samples <- unique(data$sample_id)

annotated_objects <- list()

for(sample in samples) {

  message(sprintf("\n========== Processing Sample: %s ==========\n", sample))
  niches_k = 2     # 2 niches: WM and GM
  neighbors_k = 200
  # # These samples need fine-tune parameters
  # if (sample %in% c(samples[1],samples[6], samples[11])){
  #   niches_k = 6
  #   neighbors_k = 400
  # }
  # if (sample %in% c(samples[10])){
  #   niches_k = 2
  #   neighbors_k = 100
  # }

  sample_obj <- subset(data, subset = sample_id == sample)
  # The function expects a single FOV - check what FOV ID exists
  fov_ids <- names(sample_obj@images)
  
  if(length(fov_ids) == 0) {
    warning(sprintf("No FOV found for sample %s. Skipping...", sample))
    next
  }
  
  # Use first FOV
  fov_id <- fov_ids[1]
  
  # Run annotation pipeline
  sample_annotated <- annotate_wm_gm_pipeline(
    seurat_object = sample_obj,
    region = sample,  # Use sample name as region identifier
    output_dir = output_dir,
    figure_dir = figure_dir,
    annotation_column = annotation_column,
    oligodendrocyte_celltype = reference_celltype,
    wm_markers = wm_markers,
    gm_markers = gm_markers,
    fov_id = fov_id,
    niches_k = niches_k,      # 2 niches: WM and GM
    neighbors_k = neighbors_k,  
    assay = "Vizgen"  
  )
  
  ImageDimPlot(
    sample_annotated,
    fov = fov_id,
    group.by = "predicted.niche",
    size = 1,
    cols = niche_colors,
    dark.background = FALSE,coord.fixed = T,
  ) +
    ggtitle(sprintf("Sample: %s", sample), subtitle = paste0("niches_k =", niches_k, ", neighbors_k =", neighbors_k)) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right"
    )
  annotated_objects[[sample]] <- sample_annotated
}

#===============================================================================|
# MERGE NICHE ANNOTATIONS BACK TO INTEGRATED OBJECT ----------
#===============================================================================|

message("\n========== Merging niche annotations back to integrated object ==========\n")

# Extract niche annotations from all samples
niche_annotations <- do.call(rbind, lapply(names(annotated_objects), function(sample) {
  obj <- annotated_objects[[sample]]
  data.frame(
    cell_barcode = colnames(obj),
    sample_id = sample,
    predicted.niche = obj$predicted.niche,
    niches = obj$niches,
    stringsAsFactors = FALSE
  )
}))

# Add to integrated object
data$predicted.niche <- niche_annotations$predicted.niche[
  match(colnames(data), niche_annotations$cell_barcode)
]

data$niches <- niche_annotations$niches[
  match(colnames(data), niche_annotations$cell_barcode)
]

# Save annotated integrated object
saveRDS(data, file.path(output_dir,"integrated_object_with_niches.rds"))

sample_id_col <- "sample_id" 
sobj_list <- SplitObject(data, split.by = sample_id_col)
saveRDS(sobj_list, file.path(output_dir,"integrated_list_object_with_niches.rds"))
#===============================================================================|
# PLOT ALL SAMPLES WITH NICHE ANNOTATIONS ----------
#===============================================================================|

message("\n========== Creating spatial plots for all samples ==========\n")



# Plot each sample spatially
sample_plots <- list()

for(sample in samples) {
  
  # Subset to sample
  sample_obj <- subset(data, subset = sample_id == sample)
  
  # Skip if no spatial coordinates
  if(length(sample_obj@images) == 0) next
  
  fov_id <- names(sample_obj@images)[1]
  
  # Create spatial plot colored by niche
  p1 <- ImageDimPlot(
    sample_obj,
    fov = fov_id,
    group.by = "predicted.niche",
    size = 1,
    cols = niche_colors,
    dark.background = FALSE,coord.fixed = T,
  ) +
    ggtitle(sprintf("Sample: %s", sample), subtitle = paste0("niches_k =", niches_k, ", neighbors_k =", neighbors_k)) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right"
    )
  p2 =  plot_spatial_gene_expression( gene_name = "SLC17A7", sample_obj, fov_id, "Vizgen")
  sample_plots[[sample]] <- p1+p2
}

if(length(sample_plots) > 0) {
  combined_spatial <- wrap_plots(sample_plots, ncol = 4)
  
  ggsave(
    file.path(figure_dir, "all_samples_niche_spatial.pdf"),
    combined_spatial,
    width = 25,
    height = 30,
    limitsize = FALSE
  )
  ggsave(
    file.path(figure_dir, "all_samples_niche_spatial.png"),
    combined_spatial,
    width = 25,
    height = 20,
    limitsize = FALSE
  )
  
  message(sprintf("Saved combined spatial plot: %s", 
                  file.path(figure_dir, "all_samples_niche_spatial.pdf")))
}

#===============================================================================|
# CREATE SUMMARY PLOTS ----------
#===============================================================================|

# 1. Bar plot: Niche distribution across samples
niche_summary <- data@meta.data %>%
  filter(!is.na(predicted.niche)) %>%
  group_by(sample_id, predicted.niche) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample_id) %>%
  mutate(fraction = count / sum(count))

p_niche_by_sample <- ggplot(niche_summary, 
                            aes(x = sample_id, y = fraction, fill = predicted.niche)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = niche_colors) +
  labs(
    title = "Niche Distribution Across Samples",
    x = "Sample",
    y = "Fraction",
    fill = "Niche"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(figure_dir, "niche_distribution_by_sample.pdf"),
  p_niche_by_sample,
  width = 8,
  height = 6
)

# 2. Cell type composition in each niche
celltype_by_niche <- data@meta.data %>%
  filter(!is.na(predicted.niche), !is.na(.data[[annotation_column]])) %>%
  group_by(predicted.niche, .data[[annotation_column]]) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(predicted.niche) %>%
  mutate(fraction = count / sum(count))

p_celltype_by_niche <- ggplot(celltype_by_niche,
                              aes(x = predicted.niche, y = fraction, 
                                  fill = .data[[annotation_column]])) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Cell Type Composition by Niche",
    x = "Niche",
    y = "Fraction",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

ggsave(
  file.path(figure_dir, "celltype_composition_by_niche.pdf"),
  p_celltype_by_niche,
  width = 10,
  height = 6
)

# 3. UMAP colored by niche 
if("umap" %in% names(data@reductions)) {
  p_umap_niche <- DimPlot_scCustom(
    data,
    reduction = "umap",
    group.by = "predicted.niche",
    colors_use = niche_colors,
    pt.size = 0.5
  ) +
    ggtitle("UMAP colored by Niche") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(
    file.path(figure_dir, "umap_by_niche.pdf"),
    p_umap_niche,
    width = 8,
    height = 6
  )
}

#===============================================================================|
# SAVE SUMMARY TABLE ----------
#===============================================================================|

# Create summary table
summary_table <- data@meta.data %>%
  filter(!is.na(predicted.niche)) %>%
  group_by(sample_id, predicted.niche, .data[[annotation_column]]) %>%
  summarise(cell_count = n(), .groups = "drop")

write.csv(
  summary_table,
  file.path(output_dir, "niche_celltype_summary.csv"),
  row.names = FALSE
)


rm(list = ls())
setwd("Code/MERSCOPE_Analysis/Integration")

library(Seurat)
library(tidyverse)
library(scCustomize)
source("funcs_lisi.R")
source("utils_pseudobulk_deg.R")
source("utils_manual_annotate_v2.R")
options(future.globals.maxSize = 3e+09)

npc = 100
output_path = sprintf("f_results_npc=%d", npc)
dir.create(output_path, recursive = T)
#===============================================================================|
# MERGE -----------------
#===============================================================================|
#obj_list <- readRDS("/path/to/obj_list.rds")
#data = merge(obj_list[[1]], obj_list[-1])
#saveRDS(data,"/path/to/merged.obj.rds")

data = readRDS("/path/to/merged.obj.rds")

#===============================================================================|
# NORMALIZATION ---------
#===============================================================================|
data <- NormalizeData(data, assay = "Vizgen") %>%
  ScaleData(assay = "Vizgen", features = rownames(data))

data <- SCTransform(data,
                    assay = "Vizgen",
                    method = "glmGamPoi",
                    seed.use = 1234,
                    variable.features.n = nrow(data),
                    vars.to.regress = "volume") # regress out cell size

# Force all SCT models to reference Vizgen
for(i in seq_along(data@assays$SCT@SCTModel.list)) {
  if(!is.null(data@assays$SCT@SCTModel.list[[i]])) {
    data@assays$SCT@SCTModel.list[[i]]@umi.assay <- "Vizgen"
  }
}


# Now PrepSCTFindMarkers should work
data <- PrepSCTFindMarkers(data, assay = "SCT", verbose = TRUE)

#===============================================================================|
# INTEGRATION -----------
#===============================================================================|

## PCA --------------------------

data <- RunPCA(data, assay = "SCT", npcs = 100, seed.use = 1234, features = rownames(data))

pdf(file.path( "pca_elbow_plot.pdf"), height = 3, width = 6)
ElbowPlot(data, ndims = 100) + geom_hline(yintercept = 1, lty = 2, col = "blue")
dev.off()

pdf(file.path("pca_dim_loadings.pdf"), height = 16, width = 12)
VizDimLoadings(data, dims = 1:10, reduction = "pca")
dev.off()


## INTEGRATION w/ HARMONY --------
set.seed(1234)
pdf(file.path(output_path, "harmony_convergence.pdf"), height = 3, width = 5)
data <- IntegrateLayers(
  object = data, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",npcs = npc,
  verbose = FALSE
)
dev.off()

## UMAP w/o Integration -----------
set.seed(1234)
tmp <- FindNeighbors(data,
                     assay = "SCT",
                     reduction = "pca",
                     dims = 1:npc,
                     k.param = 30)
set.seed(1234)
tmp <- RunUMAP(object = tmp,
               assay = "SCT",
               reduction = "pca",
               dims = 1:npc,
               umap.method = "uwot",
               reduction.key = "UMAP_",
               seed.use = 1234,
               n.neighbors = 30,
               return.model = TRUE,
               reduction.name = "umap.unintegrated")

tmp$batch <- tmp$sample_id
tmp <- calculate_lisi(data = tmp, split_by = "batch", reduction = "umap.unintegrated") #ilisi
tmp <- calculate_lisi(data = tmp, split_by = "predicted.cell.type", reduction = "umap.unintegrated") #clisi
saveRDS(DietSeurat(tmp, assays = "SCT", dimreducs = "umap.unintegrated"),  file.path(output_path,"unintegrated.rds"))

## UMAP w/ Harmony ---------------
set.seed(1234)
data <- FindNeighbors(data,
                      assay = "SCT",
                      reduction = "harmony",
                      dims = 1:npc,
                      k.param = 30)
set.seed(1234)
data <- RunUMAP(data,
                assay = "SCT",
                reduction = "harmony",
                dims = 1:npc,
                umap.method = "uwot",
                reduction.key = "UMAP_",
                seed.use = 1234,
                n.neighbors = 30,
                return.model = TRUE,
                reduction.name = "umap")

data$batch <- data$sample_id
data <- calculate_lisi(data = data, split_by = "batch", reduction = "umap")
data <- calculate_lisi(data = data, split_by = "predicted.cell.type", reduction = "umap")
#===============================================================================
# DE NOVO CLUSTERING -----------
#===============================================================================

data <- FindClusters(data,
                     resolution = 0.7,
                     algorithm = 1,
                     random.seed = 1234)

#===============================================================================
# FIND MARKERS -----------
#===============================================================================
## DEG analysis -------------
cluster_resoln = "SCT_snn_res.0.7"
data@meta.data$de_novo_clus = data@meta.data[[cluster_resoln ]]
Idents(data) = "de_novo_clus"

results <- run_pseudobulk_deseq2(
  seurat_obj = data,
  cluster_var = "de_novo_clus",
  sample_var = "sample_id"
)

#===============================================================================
# ANNOTATE -----------
#===============================================================================

annotation_results <- annotate_clusters_manual(
  deseq_results = results$markers_df,
  marker_panels = gold_standard_markers(),
  min_log2fc = 0,
  ambiguity_threshold =  0.7,
  output_prefix = file.path(output_path,"annotation",paste0("de_novo_annotation",cluster_resoln))
)

data <- add_annotations_to_seurat(
  seurat_obj = data,
  annotation_result = annotation_results,
  cluster_col = "de_novo_clus"  
)


#group the ambiguous (for npc = 100, "SCT_snn_res.0.7")
# The Unassigned cluster expresses key excitatory markers so we manually assign it to Excitatory_Neuron
data@meta.data = data@meta.data%>%mutate(de.novo.celltype =case_when(de.novo.celltype%in%c("Unassigned") ~ "Excitatory_Neuron",
                                                                     de.novo.celltype%in%c("Astrocyte-like") ~ "Astrocyte",
                                                                TRUE ~ de.novo.celltype) )

#===============================================================================
# EVALUATE INTEGRATION & ANNOTATIONS -----------
#===============================================================================

p1 = DimPlot_scCustom(data, reduction = "umap", group.by = c("sample_id"), colors_use = create_batch_colors())
p2 = DimPlot_scCustom(data, reduction = "umap", group.by = cluster_resoln)
p3 = DimPlot_scCustom(data, reduction = "umap", group.by = c("de.novo.celltype","predicted.cell.type"), colors_use = create_celltype_colors())
ggsave((p1+p2)/(p3), filename= file.path(output_path,paste0("umap.harmony_",cluster_resoln, ".pdf")), width = 16, height = 10)

pdf(file.path( output_path,paste0("celltype_core_marker_dp_de.novo.cell.type_",cluster_resoln,".pdf")), height = 7, width = 6)
p = Clustered_DotPlot(data,features = gold_standard_markers()%>%unlist,group.by = "de.novo.celltype",colors_use_exp = colorspace::diverge_hsv(100)  )
print(p)
dev.off()


saveRDS(data, file.path(output_path, "integrated.rds"))


## Calculate percentages ---------------------------
celltype_comp <- data@meta.data %>%
  group_by(sample_id, de.novo.celltype) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(sample_id) %>%
  mutate(percentage = (count / sum(count)) * 100) %>%
  ungroup()
median_data <- celltype_comp %>%
  group_by(de.novo.celltype) %>%
  summarise(median_pct = median(percentage)) %>%
  ungroup()

ggplot(celltype_comp, aes(x = reorder(de.novo.celltype, percentage, median),
                          y = percentage,
                          fill = de.novo.celltype)) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1.5) +
  geom_boxplot(outliers = F) +
  geom_text(data = median_data,
            aes(x = reorder(de.novo.celltype, median_pct, median),
                y = 0.05,
                label = sprintf("%.2f%%", median_pct)),
            hjust = -0.2, vjust = 0.5,size = 3.5,
            color = "black") +
  scale_y_log10(
    breaks = c(0.1, 1, 5, 10, 20, 50, 100),
    labels = c( "0.1", "1", "5", "10", "20", "50", "100")
  ) +
  labs(
    title = "Cell Type Percentage Distribution Across Tissues",
    x = "Cell Type",y = "Percentage (%)"
  ) +
  scale_fill_manual(values = create_celltype_colors())+
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed")
  ) +
  coord_flip(clip = "off")
ggsave(file.path(output_path, "celltype_boxplot_log.pdf"), width = 6, height = 5)



saveRDS(data, file.path(output_path, "integrated.rds"))

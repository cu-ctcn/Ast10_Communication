#----------------------------------------------------------#|
#     Identify astrocytes from mouse scRNAseq data 
#               (data from Huang et al)
#----------------------------------------------------------#|
#                Natacha Comandante-Lou 
#                 (nc3018@columbia.edu)
#----------------------------------------------------------#|
rm(list = ls())
setwd("6.PLXNB1_Validation/mouse")
library(tidyverse)
library(Seurat)
library(SeuratWrappers)
library(scCustomize)
library(nichenetr)
library(SeuratDisk)
sources('utils/utils_mouse_celltype_mapping_pipeline.R')
################################################################################|
# Read Reference
################################################################################|
reference.sub = LoadH5Seurat("/home/mf3362/20220728-Frauke-cytokine/20220922-downsample-rosmap/merge-propotionally-downsampled-celltypes.out/prop_0.1/step5-umap.h5Seurat")%>%
  subset(x=.,downsample = 8000)

state.levels = reference.sub@meta.data$cell.type%>%unique
my_pal = scCustomize_Palette(length(state.levels)+1,ggplot_default_colors = FALSE)[-2]
names(my_pal) = state.levels
################################################################################|
# Read Query
################################################################################|

query.list <- readRDS("../seurat.obj/mouse_plxnb1_ko_sobj_list.rds")

# Prepare query
## convert mouse to human genes
query.list = lapply(query.list, function(x){
  m_genes = x@assays[["RNA"]]@data%>%rownames()
  h_genes = convert_mouse_to_human_symbols(m_genes)%>%na.omit()
  counts_h = x@assays[["RNA"]]@counts[names(h_genes),]

  rownames(counts_h) = h_genes
  x[["RNA_human"]] = CreateAssayObject(counts =counts_h)
  DefaultAssay(x) = "RNA_human"
  return(x)
})

################################################################################|
# Run mapping
################################################################################|
source("/mnt/mfs/ctcn/team/natacha/PLXNB1-KO-mouse_Huang_et_al/Analysis/20240826/script/plxnb1_ko_mouse_ref_mappping_pipeline_v2.R")
h.markers = c("RBFOX3", "SLC17A7", "SLC17A6", "GAD1", "GAD2", "CUX2",
              "RORB", "THEMIS", "CTGF", "PVALB", "SST", "VIP", "LAMP5", "AQP4", "ALDH1L1",
              "FGFR3", "GFAP", "MOG", "MAG", "MBP", "PDGFRA", "AIF1", "TMEM119", "PTPRC", "CLDN5", "PDGFRB", "RGS5")

for (n_features in c(4000)){
  for (npcs in c(30)){
    for (k in c(10)){
  reference.sub <- FindVariableFeatures(reference.sub, selection.method = "vst", nfeatures = n_features)%>%
    RunPCA(npcs = npcs)%>%
    FindNeighbors(features = VariableFeatures(object = reference.sub),k.param = k)%>%
    RunUMAP( dims = 1:npcs,return.model = TRUE)
  

  output.folder = file.path("ref_mapping_res")

  output = plxnb1_ko_mouse_label_transfer_pipeline(query.list,query.assay = "RNA_human",
                                                                 reference.sub,
                                                                 n_features = n_features, # for SCTransform
                                                                 normalization.method ="SCT",
                                                                 k = k,
                                                                 npcs = npcs,
                                                                 state.levels = state.levels,
                                                                 output.folder = output.folder,
                                                                 label_oi = "cell.type",
                                                                 save.query = T,
                                                                 save.ref = FALSE,
                                                                 plot.markers = h.markers,
                                                                 pal = my_pal
  )
      
    }
  }
}


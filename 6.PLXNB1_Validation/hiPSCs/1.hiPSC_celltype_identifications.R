#----------------------------------------------------------#|
#      Identify most astrocyte-like cells from iPSCs
#----------------------------------------------------------#|
#                Natacha Comandante-Lou 
#                 (nc3018@columbia.edu)
#----------------------------------------------------------#|
rm(list = ls())
setwd("6.PLXNB1_Validation/hiPSCs")
library(tidyverse)
library(Seurat)
library(SeuratWrappers)
library(scCustomize)
library(nichenetr)
library(SeuratDisk)
sources('utils/utils_celltype_mapping_pipeline.R')
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

query.list <- readRDS(file.path("../data/seurat_obj", "sobj_list_postqc.rds")) # post-qc object

query.list = lapply(query.list, function(s){
  s@meta.data = s@meta.data%>%mutate(condition = case_when(orig.ident=="BIZH02_5533C4PB0WT_0" ~ "Control" ,
                                                         orig.ident=="BIZH02_5533C4PB7KO_0" ~ "PLXNB1-KO"))
  return(s)
})

################################################################################|
# Run mapping
################################################################################|
source("/mnt/mfs/ctcn/team/natacha/PLXNB1-KO-hiPSCs/script/1._celltype_mapping_pipeline_v2.R")
# celltype markers 
h.markers = c("RBFOX3", "SLC17A7", "SLC17A6", "GAD1", "GAD2", "CUX2",
              "RORB", "THEMIS", "CTGF", "PVALB", "SST", "VIP", "LAMP5", "AQP4", "ALDH1L1",
              "FGFR3", "GFAP", "MOG", "MAG", "MBP", "PDGFRA", "AIF1", "TMEM119", "PTPRC", "CLDN5", "PDGFRB", "RGS5")

n_features = 1000
npcs = 30
k = 10

reference.sub <- FindVariableFeatures(reference.sub, selection.method = "vst", nfeatures = n_features)%>%
  RunPCA(npcs = npcs)%>%
  FindNeighbors(features = VariableFeatures(object = reference.sub),k.param = k)%>%
  RunUMAP( dims = 1:npcs,return.model = TRUE)


output.folder = file.path("../cell-type_ref_mapping_res")

output = PLXNB1_ko_hipsc_label_transfer_pipeline(query.list,query.assay = "RNA",
                                                               reference.sub,
                                                               n_features = n_features, # for SCTransform
                                                               normalization.method ="SCT",
                                                               k = k,
                                                               npcs = npcs,
                                                               state.levels = state.levels,
                                                               output.folder = output.folder,
                                                               label_oi = "cell.type",
                                                               harmony_query = T,
                                                               save.query = T,
                                                               save.ref = FALSE,
                                                               plot.markers = h.markers,
                                                               pal = my_pal
)
      


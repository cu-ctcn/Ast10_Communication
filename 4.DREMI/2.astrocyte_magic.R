#------------------------------------------------#|
#         MAGIC imputation of astrocyte 
#           gene expression from CUIMC1
#------------------------------------------------#|
#           Natacha Comandante-Lou 
#             (nc3018@columbia.edu)
#------------------------------------------------#|
setwd("4.DREMI")
rm(list=ls())
gc()
library(Seurat)
library(SeuratDisk)
library(Rmagic)
library(ggplot2)
library(readr)
library(viridis)
library(phateR)
library(logger)
library(tidyverse)
library(scCustomize)
###############################################################################|
# DATA PREP ---------
###############################################################################|
celltype_oi = 'astrocytes'
sobj <- LoadH5Seurat(sprintf("/home/mf3362/snrnaseq/annotation.2024-02-26/Seurat/%s.h5Seurat",celltype_oi))

set.seed(10023)
ncells = 25000 # randomly sample ncells for each cell-state to ensure equal representation
sobj.sm <- sobj[, sample(colnames(sobj), size =ncells, replace=F)]
sobj.sm = NormalizeData(sobj.sm, normalization.method = "LogNormalize", scale.factor = 10000)
###############################################################################|
# MAGIC IMPUTATION ----------
###############################################################################|
reticulate::use_condaenv("/home/nc3018/conda/my-envs/magic")
source("/mnt/mfs/ctcn/team/natacha/RScript/run_magic.R")
DefaultAssay(sobj.sm) = "RNA"
data = GetAssayData(sobj.sm,slot = "count",assay = "RNA")%>%t()

# Keep genes expressed in at least 10 cells
keep_cols <- colSums(data > 0) > 10
data <- data[,keep_cols]

# Normalize and transform
data <- library.size.normalize(data)
data <- sqrt(data)

# run MAGIC
all_genes = colnames(data)


## selected params
t = 5; knn=10
#auto-selected t=11
magic_output_dir  = "magic_output_n=25000"
final_res = lapply(t,run_magic,data,all_genes,
                   knn=knn,
                   plot=TRUE,
                   plot_only = FALSE,
                   init = NULL,
                   gene_x_viz="PLXNB1",gene_y_viz="SLC38A2",
                   save.as = c("magic"),
                   output.dir = magic_output_dir )


data_MAGIC = final_res[[1]][["data_MAGIC"]]
magic_assay <- CreateAssayObject(data = data_MAGIC[["result"]]%>%as.matrix()%>%t())

# add this assay to the previously created Seurat object
sobj.sm[["MAGIC"]] <- magic_assay
Misc(sobj.sm, slot = "magic") = data_MAGIC$operator

saveRDS(sobj.sm, file.path(magic_output_dir , sprintf("knn=%d_t=%d/sobj_after_magic.rds", knn, t)))



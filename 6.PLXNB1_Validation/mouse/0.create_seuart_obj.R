#----------------------------------------------------------#|
#       Preprocess mouse scRNAseq data (from Huang et al)
#               with or without PLXNB1 KO
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

data_dir <- '~/PLXNB1-KO-mouse_Huang_et_al/Single cell RNA-seq/Cell Ranger outs filtered_feature_bc_matrix'
meta = read_tsv("~/PLXNB1-KO-mouse_Huang_et_al/Single cell RNA-seq/meta2.tsv")
file_list = list.files(data_dir)%>%
  .[!grepl(pattern = ".tsv",x = .)]
  # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz

############################################################|
# QC:
# - removing cells with ≤200 or ≥5,000 genes
# - removing cells with  ≥50,000 UMI counts
# - keep genes express more than 2 cells
# - removing cells with ≥ 10% mitochondrial reads
############################################################|

sobj_list = lapply(file_list, function(f){
  data <- Read10X(data.dir = file.path(data_dir,f))
  sobj = CreateSeuratObject(counts = data, project = meta%>%filter(Sample==str_replace(f,"P",""))%>%paste(., collapse  = "-"))
  
  # QC removing cells with ≤200 or ≥5,000 genes,umi count ≥ 50000, keep genes express more than 2 cells
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
  sobj <- subset(sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt <= 10 & nCount_RNA <50000)
  all_genes = sobj@assays[["RNA"]]@counts%>%rownames()
  keep_genes = all_genes[(sobj@assays[["RNA"]]@counts%>%rowSums()>2)]
  sobj  <- subset(sobj , features = keep_genes)
  
  sobj <- SCTransform(sobj, verbose = FALSE)

  return(sobj)
  })
dir.create("../seurat.obj", recursive = T)

saveRDS(sobj_list, file = "../seurat.obj/mouse_plxnb1_ko_sobj_list.rds")

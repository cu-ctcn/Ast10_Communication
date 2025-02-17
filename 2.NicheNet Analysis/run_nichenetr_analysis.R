#-------------------------------------------------------------------#|
#                 NicheNet Analysis on Ast10 signature
#-------------------------------------------------------------------#|
#                   Natacha Comandante-Lou 
#                   (nc3018@columbia.edu)
#-------------------------------------------------------------------#|
# Use Ast10 signature genes to infer upstream ligands microglia and/or astrocyte subsets

rm(list = ls())
gc()
renv::load("2.NicheNet Analysis")
setwd("2.NicheNet Analysis")
source("nichenetr_pipeline.R")

library(nichenetr) 
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(SeuratDisk)
##############################################################################|
# READ DATA & NICHENET RESOURCES -------------
##############################################################################|
# A combined seurat objects with the microglia and astrocytes from Green et al
seuratObj = LoadH5Seurat('/mnt/mfs/ctcn/team/natacha/ROSMAP_snRNAseq/MultiNicheNet Analysis/resources/mic_ast_2024-02-26.h5Seurat')

DefaultAssay(seuratObj) = "RNA"

seuratObj.437 = subset(seuratObj,subset = individual.in.landscape ==T) #remove macrophages and monocytes before the analysis
log_info(seuratObj.437@meta.data$projid%>%unique()%>%length)

log_info("LogNormalizing..")
print(seuratObj.437@assays[["RNA"]]@data[1:10,1:10])
seuratObj.437 <- NormalizeData(seuratObj.437, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(seuratObj.437)
log_info("Scaling Data..")
seuratObj.437 <- ScaleData(seuratObj.437, features = all.genes)
print(seuratObj.437@assays[["RNA"]]@data[1:10,1:10])
print(seuratObj.437@assays[["RNA"]]@scale.data[1:10,1:10])
log_info("Done preprocessing.")
##############################################################################|
# PARAMETERS -------------
##############################################################################|
# Prioritization critera:
# (See nichenetr for more detail)
# * Upregulation of the ligand in a sender cell type compared to other cell types: `de_ligand`
# * Upregulation of the receptor in a receiver cell type: `de_receptor`
# * Average expression of the ligand in the sender cell type: `exprs_ligand`
# * Average expression of the receptor in the receiver cell type: `exprs_receptor`
# * Condition-specificity of the ligand across all cell types: `ligand_condition_specificity` (none in this case)
# * Condition-specificity of the receptor across all cell types: `receptor_condition_specificity` (none in this case)
args = commandArgs(trailingOnly = T)
pct_thres = as.numeric(args[1])
prioritize = as.logical(args[2])
w.de_ligand = as.numeric(args[3])
w.de_receptor = as.numeric(args[4])
w.activity_scaled = as.numeric(args[5])
w.exprs_ligand = as.numeric(args[6])
w.exprs_receptor = as.numeric(args[7])
w.ligand_condition_specificity = as.numeric(args[8])
w.receptor_condition_specificity = as.numeric(args[9])
output_folder = args[10]

## DE output path ----
de_output_folder = "ROSMAP_snRNAseq/Differential_Expression/de_analysis_output"
subtype.cutoff = "adjp=0.05_logfc=1.00"


## NicheNet Parameters----
receiver = "Ast.10"
sender = c("Ast.1","Ast.2","Ast.3","Ast.4","Ast.5","Ast.6","Ast.7","Ast.8","Ast.9","Ast.10",
           "Mic.1","Mic.2","Mic.3","Mic.4","Mic.5","Mic.6","Mic.7","Mic.8","Mic.9","Mic.10",
           "Mic.11","Mic.12","Mic.13","Mic.14","Mic.15","Mic.16")

geneset_oi = read.csv(sprintf("%s/%s_de_genes_pairwise_%s.csv",de_output_folder,receiver,subtype.cutoff))%>%pull(gene) %>%unique()

### nichenet de info

DE_table = NULL# no pre-calculated DE table
expression_info = NULL # no pre-calculated expression_info
##############################################################################|
# RUN NICHENET PIPELINE -------------
##############################################################################|
nn_output = run_nichenet_pipeline(seuratObj.437,
                                  receiver = receiver,
                                  sender = sender,
                                  
                                  run_nichenet = TRUE,
                                  expressed_genes_receiver = NULL,
                                  expressed_genes_sender = NULL,
                                  gene_set_oi = genset_oi,
                                  pct_thres = pct_thres, # threshold for defining expressed genes
                                  
                                  # Prioritization weights
                                  prioritize = TRUE,
                                  w.de_ligand = w.de_ligand,
                                  w.de_receptor = w.de_receptor,
                                  w.activity_scaled = w.activity_scaled,
                                  w.exprs_ligand = w.exprs_ligand,
                                  w.exprs_receptor = w.exprs_receptor,
                                  w.ligand_condition_specificity =w.ligand_condition_specificity,
                                  w.receptor_condition_specificity =w.receptor_condition_specificity,
                                  celltype_colname = "state",
                                  DE_table = DE_table, 
                                  expression_info = expression_info,
                                  
                                  # Plotting
                                  plot.prior.table = NULL,
                                  plot.sender = sender,
                                  top_n = 100,
                                  
                                  figure_folder = file.path(output_folder, receiver),
                                  output_folder = file.path(output_folder, receiver)
                                  
                                  
)


saveRDS(nn_output,file = file.path(output_folder, receiver,"nn_output.rds"))

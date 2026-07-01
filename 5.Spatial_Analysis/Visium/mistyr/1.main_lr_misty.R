#------------------------------------------------#|
#     MISTy analysis to predict Ast10 markers
#         based on local ligand expression
#------------------------------------------------#|
#           Natacha Comandante-Lou 
#             (nc3018@columbia.edu)
#------------------------------------------------#|
rm(list=ls())
setwd("5.Spatial_Analysis/mistyr")
options(future.globals.maxSize= 10*1024^3)
renv::load("5.Spatial_Analysis/mistyr")

library(tidyverse)
library(Seurat)
future::plan(future::multisession) 
library(mistyR)
source("lr_misty_pipeline.R")
source("../utils/utils_load_genes_oi.R")
## Read DEGs -------------------------------------------------
receiver = "Ast.10"

## Read NicheNet-inferred ligands and receptors---------------
nn_output_folder = "../../2.NicheNet Analysis/nn_output"
plsr_output_folder = "../../3.PLSR/plsr_analysis"

resource = load_genes_oi(nn_output_folder,plsr_output_folder, receiver)
nn_ligands = resource[["nn_ligands"]];nn_receptors = resource[["nn_receptors"]]
all_ligands = resource[["all_ligands"]];all_receptors = resource[["all_receptors"]]
pls_ligands = resource[["pls_ligands"]];vip_ligands = resource[["vip_ligands"]]
markers_oi = resource[["markers_oi"]]

################################################################################|
# Misty Params-------
################################################################################|
## Views ---------------------------------
receiver.view = markers_oi
receptor.view = nn_receptors
ligand.view = all_ligands

##-----------------------------------------

view_types <- list("main" = "intra",
                   "juxta" = "juxta")
view_assays = list("main" = "SCT",
                   "juxta" =  "SCT")
view_params <- list("main" = NULL,
                    "juxta" = 160*3)
view_features = list("main" = c(receiver.view),
                     "juxta" =  c(ligand.view))
results_folder = sprintf("mistyr_results", view_params[["juxta"]])

###############################################################################|
# Set paths -------
###############################################################################|

res.dir = "output"
## Getting sample annotations --------------------------------------------------
slide_files_folder <- "objects_rm_wm/"
slide_files <- list.files(slide_files_folder)
slide_ids <- gsub("[.]rds", "", slide_files)


###############################################################################|
# Run MistyR--------
###############################################################################|
misty_outs <- map(slide_files, run_misty_pipeline, view_types = view_types,
                  view_assays=view_assays,view_params = view_params, view_features = view_features, 
                  score_module = F,
                  out.dir = file.path(res.dir,results_folder))

suffix = view_assays[["main"]]
slide_files <- list.dirs(file.path(getwd(),res.dir, results_folder),recursive = F)%>%.[grepl(suffix,.)]
pdf(file.path(res.dir[1],results_folder,"collected_improvement_stats.pdf"), width = 30,height = 20)
print(collect_results(slide_files) %>% plot_improvement_stats())
dev.off()

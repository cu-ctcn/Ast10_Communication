
rm(list=ls())
gc()
setwd("3.PLSR")

library(logger)
library(Seurat) 
library(tidyverse)
library(SeuratDisk)
library(ggpubr)
library(reshape2)
library(abind)
library(reshape2)
library(colorspace)
library(abind)
library(logger)

source('../utils/plsr_helper_func.R')
source('utils/run_plsr_pipeline.R')
###############################################################################|
# LOAD DATA ------
###############################################################################|
args = commandArgs(trailingOnly = T)

nichenet_output_folder = args[1]

# Gene Expression
mic.ast = LoadH5Seurat('/mnt/mfs/ctcn/team/natacha/ROSMAP_snRNAseq/MultiNicheNet Analysis/resources/mic_ast_2024-02-26.h5Seurat')
DefaultAssay(mic.ast) = "RNA"
Idents(mic.ast)
mic.ast = subset(mic.ast,idents = c("Macrophages", "Monocytes"), invert = TRUE) #remove macrophages and monocytes before the analysis
mic.ast.437 = subset(mic.ast,subset = individual.in.landscape ==T) 

log_info("LogNormalizing..")
mic.ast.437 <- NormalizeData(mic.ast.437, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(mic.ast.437)
log_info("Scaling Data..")
mic.ast.437 <- ScaleData(mic.ast.437, features = all.genes)
log_info("Done preprocessing.")

# Get cell state frequency
freq.w = readRDS("1.snucRNAseq_Cell-state_Meta-analysis/cell-state_freq/freq_rds/Discovery/state_freq.rds")
# Frequency correlation
freq.cor = freq.w%>%select(-c("projid"))%>%cor(use="complete.obs")  #discard the entire row if an NA is present


plsr_output_folder = "receptor_plsr_output"
dir.create(plsr_output_folder,recursive = TRUE)
###############################################################################|
# PSEUDOBULK EXPRESSION -------
###############################################################################|
receiver = "Ast.10"
top_n_signals = 100

## Get ligands of interests

prior_table = read.csv(file = file.path(nichenet_output_folder,receiver,"prioritized_table.csv"))
prior_table.top = prior_table[c(1:top_n_signals),]

ligand_oi = prior_table.top%>%pull(ligand)%>%unique()
receptor_oi = prior_table.top%>%pull(receptor)%>%unique()

## Get population average ----------------------------------------------
### average data per projid, cell type ---------------------

# * Assumes that the data has been log normalized (assays = "RNA", slot = "data)
# * AverageExpression  = log(Mean(exp(log_norm_gene_exp)))

data.id.state = AverageExpression(mic.ast.437,features = unique(c(receptor_oi)),group.by = c("cell.type","projid"),
                                  assays = "RNA",slot = "data",
                                  return.seurat = TRUE)
snames <- str_split(rownames(data.id.state@meta.data), "_",simplify = TRUE)
data.id.state  = AddMetaData(data.id.state , metadata = snames[,1], col.name = "state")
data.id.state  = AddMetaData(data.id.state , metadata = snames[,2], col.name = "projid")

df.id.state = as.data.frame(t(data.id.state@assays[["RNA"]]@data ))
df.id.state$state = snames[,1]
df.id.state$projid = snames[,2]


gc()


## 
# Count receiver cells
min_cells = 0
meta.data = mic.ast.437@meta.data
rec_count = group_by(meta.data, projid,state)%>%tally()%>%filter(state%in%"Ast.10")
keep_donor = rec_count%>%filter(n >=min_cells)%>%pull(projid)
###############################################################################|
# CONTRUCT MODEL IN/OUTPUT -----
###############################################################################|

df.plsr = filter(df.id.state,state%in%"Astrocyte")%>%filter(projid %in% keep_donor)
df.plsr = merge(df.plsr, select(freq.w, c("projid",receiver)),by="projid")


### input -----
X = select(df.plsr, -c("projid",receiver,"state"))%>%`rownames<-`(df.plsr$projid)

### output -------
Y =log(select(df.plsr,all_of(receiver))*100+1)%>%`rownames<-`(df.plsr$projid)



###############################################################################|
# RUN PLSR PIPELINE -----
###############################################################################|

source("/utils/run_plsr_pipeline.R")
pls_res = run_plsr_pipeline (X,Y,
                             tts.seed = 132,
                             test.size = 0.1,
                             npls = NULL,
                             n_splits = 5,
                             n_repeats = 100,
                             output.dir =plsr_output_folder)

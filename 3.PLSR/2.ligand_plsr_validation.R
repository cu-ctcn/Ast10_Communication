#-----------------------------------------#|
# Predict Ast.10 frequency in new data
#       using trained ligand model
#-----------------------------------------#|
#       Natacha Comandante-Lou 
#       (nc3018@columbia.edu)
#-----------------------------------------#|
rm(list = ls())
setwd("3.PLSR")
library(Seurat)
library(logger)
library(tidyverse)

source('../utils/plsr_helper_func.R')
source('utils/run_plsr_pipeline.R')

args = commandArgs(trailingOnly = T)

data.res.dir = args[1] #folder with cell-state annnotated seurat objects
study_name = args[2] # dataset name
freq_path = args[3] # path to frequency data
plsr_res_path = args[4] # path to trained plsr model
donor_oi_path = args[5] # selected list of donor projids
output.dir = args[6]

receiver = "Ast.10"
celltype_oi = c("microglia","astrocytes")
donor_oi = readRDS(donor_oi_path)%>%filter(reference==study_name)%>%pull(projid)


###############################################################################|
# Get Single-cell Data ------
###############################################################################|
df_params <- read_csv(file.path(data.res.dir, "optimized_params.csv"))%>%
  mutate(file_path = file.path(data.res.dir,cell.type,"ref_mapping_res",sprintf("n_features=%d/npcs=%d_k=%d", n_features, npcs, k)))%>%
  mutate(ref_file_path = file.path(data.res.dir,cell.type,"freq_dj.rds"))%>%
  filter(cell.type %in% celltype_oi)

sobj_list = lapply(df_params$file_path, function(x){readRDS(file.path(x,"annotated_query.rds" ))})
seuratObj <- merge(sobj_list[[1]], y = sobj_list[[2]], add.cell.ids = celltype_oi, project = study_name)
seuratObj@meta.data[["projid"]] = as.character(seuratObj@meta.data[["projid"]])%>%stringr::str_pad(.,8,pad = "0")
###############################################################################|
# Get ligand Pseudobulk --------
###############################################################################|
## Updated NicheNet Database --------
# (Source: https://zenodo.org/record/7074291/files/)

# NicheNet’s ligand-receptor data sources
lr_network = readRDS("/mnt/mfs/ctcn/team/natacha/NicheNet_Database/7074291/lr_network_human_21122021.rds")

# Ligand-target model: 
ligand_target_matrix = readRDS("/mnt/mfs/ctcn/team/natacha/NicheNet_Database/7074291/ligand_target_matrix_nsga2r_final.rds") # target genes in rows, ligands in columns

# Weighted integrated network
weighted_networks = readRDS("/mnt/mfs/ctcn/team/natacha/NicheNet_Database/7074291/weighted_networks_nsga2r_final.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))

all_ligands = lr_network$from%>%unique()
all_receptors = lr_network$to%>%unique()

genes_oi = intersect(c(all_ligands, all_receptors), rownames(seuratObj@assays[["RNA"]]@data))


expr = GetPseudobulk_seurat_v5(seuratObj, assay_oi = "RNA",normalization.method = "LogNormalize",
                     genes_oi,group.by.vars = c("predicted.id","projid"),
                     save = TRUE, output.dir = output.dir, output.type = c("seurat","data.frame"))

###############################################################################|
# Assemble pls input output--------
###############################################################################|

## Cell-state frequencies----
freq.w = readRDS(freq_path)


## Correlation from Discovery dataset---
cor.res <- readRDS(file.path(plsr_res_path, "ligand-oi-in-sender_pearson_corr_Ast.10.rds"))

df = expr[["df.expr"]]%>%
  merge(select(freq.w, c("projid", receiver)), by = "projid")%>%
  filter(.data[[receiver]]> 0)
df.plsr = AssembleInput(cor.res, df, method = "mean.scaled",receiver = receiver,
                        r.thres = 0.1, p.thres = 0.05)

## Load prior model to get x_var -----
jb = reticulate::import("joblib")
sk = reticulate::import("sklearn")
## Load Model------
plsr.res <- readRDS(file.path(plsr_res_path ,"plsr_res_response=Ast.10.rds"))
pls_final= jb$load(file.path(plsr_res_path ,"pls_final_response=Ast.10.joblib"))

# note: temporary fix - get the original x order while pulling x_var from vip columns
genes_oi = tibble(gene = plsr.res$loadings%>%rownames()%>% # same order as model
  setdiff(.,receiver))
df.vip = plsr.res$vip_scores
df.vip$gene = lapply(df.vip$x_var, function(x){ str_split(x,"(\ |/|\\(|\\))",simplify = TRUE)[1]})%>%do.call("c",.)
df.vip = df.vip%>%left_join(genes_oi,.)
x_var_oi = df.vip%>%pull(x_var)%>%na.omit

max_nan = max(colSums(is.na(df.plsr%>%dplyr::select(contains(x_var_oi)))))
df.plsr = df.plsr[,colSums(is.na(df.plsr))<=max_nan] #filter out columns with more than 4 NaN

#merge freq column
df.plsr = merge(df.plsr, freq.w%>%select(c("projid",receiver)), by = "projid")%>%
  mutate(logfreq = log(.data[[receiver]]*100+1))
#filter donors with at least 1 Nan
df.plsr = df.plsr[!(rowSums(is.na(df.plsr))>=1),]

###############################################################################|
# Validation with trained model --------
###############################################################################|


## Set up input and output ----

prior_pls_pred = function(opt,df.plsr){
  if(opt=="non-overlap"){
    
    # Non overlapping donors with CUIMC1
    X = df.plsr%>%`rownames<-`(df.plsr$projid)%>%filter(projid%in% donor_oi)%>%select(contains(x_var_oi))
    Y = df.plsr%>%`rownames<-`(df.plsr$projid)%>%filter(projid%in% donor_oi)%>%select(logfreq)
  }
  if(opt=="overlap"){
    # Overlapping donors with CUIMC1
    X = df.plsr%>%`rownames<-`(df.plsr$projid)%>%filter(!projid%in% donor_oi)%>%select(contains(x_var_oi))
    Y = df.plsr%>%`rownames<-`(df.plsr$projid)%>%filter(!projid%in% donor_oi)%>%select(logfreq)
  }
  if(opt=="all"){
    # All donors
    X = df.plsr%>%`rownames<-`(df.plsr$projid)%>%select(contains(x_var_oi))
    Y = df.plsr%>%`rownames<-`(df.plsr$projid)%>%select(logfreq)
  }
  
  # scale on all donors who have receiver cells
  if (dim(X)[1]>0){
  scaler_x = sk$preprocessing$StandardScaler()$fit(df.plsr%>%select(colnames(X)))
  scaler_y = sk$preprocessing$StandardScaler()$fit(df.plsr%>%select(colnames(Y)))
  Xz = scaler_x$transform(X)%>%`colnames<-`(colnames(X))%>%`rownames<-`(rownames(X))
  Yz = scaler_y$transform(Y)%>%`colnames<-`(colnames(Y))%>%`rownames<-`(rownames(Y))
  

  # Model prediction  
  Y_pred = pls_final$predict(Xz)
  
  # Calculate Q2 to evaluate model predictions
  exp_var = sk$metrics$r2_score(Yz, Y_pred)
  
  
  df.Y = tibble(data = Yz, model = Y_pred)
  p = ggpubr::ggscatter(df.Y, x = "data", y = "model",
                    add = "reg.line", conf.int = FALSE,
                    add.params = list(color = "black"),
                    cor.coef = TRUE, cor.method = "pearson", color = "black", shape = 21, fill = "purple",alpha = 0.8, size = 2.5,repel = TRUE) + # pearson coefficient (R)
    geom_abline(slope = 1, linetype = "dashed")+
    
    annotate("text",label = sprintf("variance predicted = %0.2f%%",exp_var*100), x = 0.3, y =0.6)+
    ggtitle(label = sprintf("%s (n = %d)", opt, dim(df.Y)[1]))+
    ylab('Model Prediction')+
    xlab(study_name)+
    theme_classic()
  
  # Save data
  if(opt=="all"){
    res = list(df.plsr = df.plsr,
               Xz = Xz,
               Yz = Yz,
               Y_pred = Y_pred,
               explained_variance = exp_var,
               data_path = df_params%>%pull(file_path)
    )
    saveRDS(res, file = file.path(output.dir,"results.rds"))
  }
  return(p)
  }
}


fig= lapply(c("non-overlap","overlap","all"),prior_pls_pred, df.plsr)%>%
  ggpubr::ggarrange(plotlist = .,nrow = 1, ncol = 3)
ggsave(plot = fig, filename = file.path(output.dir,sprintf("ligand_prior_pls_model_vs_%s.pdf",study_name)), width = 9.6 , height = 3.09)


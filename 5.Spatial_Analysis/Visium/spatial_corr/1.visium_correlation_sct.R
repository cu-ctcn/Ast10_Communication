#------------------------------------------------#|
#         Spatial correlation analysis 
#
#         Target spot Ast10 signature
#                     vs.
#         Neighborhood ligand signature 
#                     
#------------------------------------------------#|
#           Natacha Comandante-Lou 
#             (nc3018@columbia.edu)
#------------------------------------------------#|

rm(list = ls())
gc()

setwd("5.Spatial_Analysis/spatial_corr")
library(tidyverse)
library(dplyr)
library(scCustomize)
library(reticulate)
library(colorspace)
library(cowplot)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)
library(ggpubr)

source("../../utils/seurat_get_summary_scores.R")
source("../utils/utils_visium_corr.R")

###############################################################################|
# MODULE SCORES ----------
###############################################################################|
receiver_list = c("Ast.1","Ast.2","Ast.3","Ast.4","Ast.5","Ast.6","Ast.7","Ast.8","Ast.9","Ast.10")

# Load markers for each ast states
markers = read_excel("Supplementary Tables/Table.S2_Differentially_Expressed_Genes.xlsx")%>%
  filter(visium.analysis == T)%>%split(., f = .[["state"]])
module.features.list = lapply(markers, function(x){unique(x$gene)})%>%.[receiver_list]
names(module.features.list) = sprintf("%s_signature", receiver_list)
saveRDS(module.features.list,'ast.modules.rds')

###############################################################################|
# Load Data ----------
###############################################################################|
sobj_path = "objects_rm_wm"
slide_files = list.files(sobj_path)

slide_id_list = str_extract(slide_files,"\\d+_\\d")

## Analysis parameters ----- 
summary_method = "seurat" #method to calculate summary scores
module_assay = "SCT_reg_nuclei" #assay to calculate signature scores (dependent variables)
ligand_assay = "SCT_reg_nuclei" #assay to calculate ligand summary scores (independent variables)
include_center = F #include center spot itself when calculate correlation?
r = 3 #radius for calculating neighborhood
pxl_scale = 156
slot_oi = "data"
ligands = c("JAM2","NCAN","SEMA4D","GNAS","LRRC4B","RGMA") #only include positive VIP ligands predictive of Ast.10 (from plsr)
## paths ----- 
output.dir = "output"
dir.create(output.dir, recursive = T)
dir.create(file.path(output.dir,"example_neighbors"), recursive = T)


run_corr_pipeline = function(s, pos_only=T){

  if (pos_only){ligands = c("JAM2","NCAN","SEMA4D","GNAS","LRRC4B","RGMA")}
  
  sid = str_replace(s, ".rds", "")
  slide = readRDS(file.path(sobj_path,s))
  
  #Normalization --------------
  DefaultAssay(slide) = "Spatial"
  slide =  SCTransform(slide, assay = "Spatial",vst.flavor = "v2", vars.to.regress = c("NumNuclei","percent.mt"),
                       verbose = T,new.assay.name ="SCT_reg_nuclei", return.only.var.genes = F)


  slide= seurat_get_summary_scores(slide,
                                   module.features.list,
                                   module.key =paste0("MODULE_",module_assay),
                                   method =summary_method,
                                   assay = module_assay,
                                   out.fig.dir = file.path(output.dir, "module.scores"),
                                   slot = "data",
                                   plot = F,
                                   plot.col.limits = NULL,
                                   plot.reduction =  "umap",
                                   
                                   )


  ###############################################################################|
  # CORRELATION WITHIN NEIGHBORHOOD----------
  ###############################################################################|
  

  ## Gather Neighbor Features --------
  logger::log_info("**ligands: ", paste(ligands,collapse= ", "))
  write_lines(ligands, file = file.path(output.dir,"ligands_used.txt"))
  ## Neighbor ligand scores --

  df.nb = gatherNeighborGenes(geneset.list = list("ligands"=ligands), slide, 
                           obj_type = "Seurat",assay_oi = ligand_assay, slot = slot_oi, method = summary_method, neighbor_method = "mean",
                       r, pxl_scale, shape = "disk",include_center = include_center)
  # Neighbor ligand scores --
  df.rnd = gatherNeighborGenes(geneset.list = list("ligands"=ligands), slide, 
                           obj_type = "Seurat",assay_oi = ligand_assay, slot = slot_oi,method = summary_method, neighbor_method = "mean",
                           r, pxl_scale, shape = "random",include_center = include_center)
  
  df = merge(df.nb ,df.rnd%>%dplyr::select(-c('ligands',"x","y")), by = "SpotID")
  ## Merge with Module scores ---
  module_assay.list = paste0("MODULE_",module_assay,".", summary_method)
  module.scores = lapply(module_assay.list, function(x){
      tmp = DefaultAssay(slide)
      DefaultAssay(slide) = x
      FetchData(slide, assay = x, vars = rownames(slide@assays[[x]]@data))%>%`colnames<-`(paste0(x,"_", colnames(.)))
      })%>%
    do.call("cbind",.)%>%rownames_to_column("SpotID")
  

  meta.data = slide@meta.data
  
  df = df%>%merge(.,module.scores, by = "SpotID")%>%merge(.,meta.data, by = "SpotID")
  df$SpatialCluster = factor(df$SpatialCluster,levels = c("L1","L2","L3","L4","L5","L6"))
  
  # Correlation within neighborhood ---------
  
  x_var_list = paste0("MODULE_", module_assay,".",summary_method,"_",receiver_list,"-signature")
  y_var_list = c(sprintf("ligands_mean-disk=%0.1f",r),sprintf("ligands_mean-random=%0.1f",r))
  
  
  ### Run across y_var_list:
  cor.res = lapply(y_var_list, function(y_var){
    # Plot individual samples
    
    lapply(x_var_list, function(x){ggplot(df, aes(x = .data[[x]], y = .data[[y_var]]))+
        geom_point(color = "#615EFC")+
        stat_cor(method = "pearson")+
        geom_smooth(method = "lm", color = "black")+
        scale_color_manual(values = lay_pal)+
        theme_classic()})%>%ggarrange(plotlist = .)%>%ggsave(plot = ., filename = file.path(output.dir, paste0(sid,"_", y_var,"_corr.pdf")),width = 9, height = 8)
    
    # Output Correlation Scores
    cor.res.per.layer = lapply(x_var_list, function(x){
      lapply(unique(df$SpatialCluster), function(layer){
        .df = filter(df, SpatialCluster==layer)
        res = cor.test(.df[[x]], .df[[y_var]])
        tibble(r = res[["estimate"]][["cor"]], p.value = res[["p.value"]], x = x, y = y_var, SpatialCluster = layer)
        
      })%>%do.call("rbind",.)
      
      })%>%do.call("rbind",.)%>%mutate(sample_id = sid)
    
    cor.res.all = lapply(x_var_list, function(x){
        res = cor.test(df[[x]], df[[y_var]])
        tibble(r = res[["estimate"]][["cor"]], p.value = res[["p.value"]], x = x, y = y_var, SpatialCluster = "All")
        
      })%>%do.call("rbind",.)%>%mutate(sample_id = sid)
    
    cor.res = rbind(cor.res.per.layer,cor.res.all)
    return(cor.res)
  })%>%do.call("rbind",.)
    
  ### Plot sample neighborhood ------
  coord = dplyr::select(df,c("x","y"))
  center.x = range(coord$x)[1]+(range(coord$x)%>%diff)/2
  center.y = range(coord$y)[1]+(range(coord$y)%>%diff)/2
  center = mutate(coord,dist = (x-center.x)^2 + (y-center.y)^2)%>%
    arrange(dist)%>%dplyr::slice(1)
  xc = center%>%pull(x)
  yc = center%>%pull(y)
  example_neighbor_plot(coord, xc, yc,r,save.as = file.path(output.dir,"example_neighbors",paste0(sid,"_example_neighbors.pdf")))

  return(cor.res)

}

res = lapply(slide_files,run_corr_pipeline, pos_only  )
res = res%>%do.call("rbind",.)%>%group_by(sample_id)%>%mutate(p.adj = p.adjust(p.value,method = "BH"))
res$x = str_replace(res$x,paste0("MODULE_",module_assay,".",summary_method,"_"),"")
res$x = factor(res$x, levels = paste0(receiver_list, "-signature"))
saveRDS(res,file.path(output.dir,"cor.res.rds"))

#flip the sign of Ast.1 correlation because their signature genes marked their downregulation
flip.sig = unique(res$x)%>%.[str_detect(.,"Ast.1-")]
res=  res%>%mutate(r = case_when(x%in%flip.sig ~ -r, TRUE ~ r))
saveRDS(res,file.path(output.dir,"cor.res.flip.rds"))

#significant corr
res.sig = res%>%filter(p.adj<0.05 & abs(r) >0)

###############################################################################|
# Plot all samples -----------
###############################################################################|

sig.lab = str_replace(res$x%>%unique,paste0("MODULE_",module_assay,".",summary_method,"_"),"")

# Plot correlation (separate panels for each signature)
var_oi = sprintf("ligands_mean-disk=%0.1f",r)
p3 = ggplot()+
  geom_point(data =  res%>%filter(SpatialCluster=="All" & y==var_oi ), aes(x = r, y = -log10(p.adj)), color= "grey")+
  geom_point(data = res.sig%>%filter(SpatialCluster=="All"& y==var_oi ) , aes(x = r, y = -log10(p.adj),color = x))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5)+
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)+
  scale_color_manual(values = my_pal)+
  facet_wrap(facets ="x", nrow = 2, axis.labels = "all",axes = "all")+
  guides(color=guide_legend(title="Signature",override.aes = list(size=5)))+
  ggtitle("Correlation")+
  theme_classic()

ggsave(plot = p3, filename = file.path(output.dir,"cor.volcano.facet.pdf"), width =8, height = 3)




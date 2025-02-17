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
library(dplyr)
library(ggpubr)
library(cowplot)
load("../utils/pathway_resources.rda")

reticulate::use_condaenv("/home/nc3018/conda/my-envs/magic")

nn.pathway.dir = "pathways/nn_pathway_Ast.10"
magic_param = "knn=10_t=5"
dremi.output.dir = file.path("dremi/dremi_pathway",magic_param)


# ###############################################################################|
# # LOAD MAGIC output ---------
# ###############################################################################|

magic_output_path = file.path("magic_output_n=25000",magic_param)
sobj.sm =readRDS(file.path(magic_output_path,'sobj_after_magic.rds'))


sobj.sm = NormalizeData(sobj.sm, normalization.method = "LogNormalize", scale.factor = 10000)


###############################################################################|
# MODULE SCORES ----------
###############################################################################|
# Calculate Ast.10 signature
source("../utils/seurat_get_summary_scores.R")
receiver = "Ast.10"
#remove ligand and receptors
full_network <- readRDS(file.path(nn.pathway.dir,"SEMA4B.JAM2.NRG3.EFEMP1.LRIG1.FN1.NCAN.LRRC4B.SEMA4D.RGMA.LGALS3.GNAS.NRXN1/network.rds"))
full_attri= full_network$attributes
rm_pathway_genes = full_attri%>%filter(!Role%in%c("target"))%>%pull(Gene)
targets = setdiff(pos_targets, rm_pathway_genes) #calculate ast.10 module score using positive target genes there are not in the inferred pathway

module.features.list = list(targets)
names(module.features.list) = sprintf("%s_signatures", receiver)

output.dir = "module.score"
seed = 1865
sobj.sm= seurat_get_summary_scores(sobj.sm,
                                 module.features.list,
                                 module.key = "magic.module",
                                 method = "seurat",
                                 assay = "MAGIC",
                                 out.fig.dir = file.path(output.dir,"modules.score"),
                                 plot = TRUE,
                                 plot.col.limits = NULL,
                                 plot.reduction =  "umap",
                                 seed = seed
                                 )

saveRDS(sobj.sm, file.path(magic_output_path, "sobj_after_magic.rds"))

###############################################################################|
# DREMI FUNCTIONS-------------
###############################################################################|
source("utils/run_dremi.R")

# run dremi on the list of genes
run_dremi_wrapper = function(xx, data.magic, k = 10, n_bins = 36, shuffle=T, plot_dremi,role,output.dir){
  
  # gene x
  X = data.magic[["result"]]%>%as.data.frame%>%select(xx)
  dremi = run_dremi(X = X, Y = Y,output.dir = file.path(output.dir,"dremi.heatmap",role),
                    k = k, n_bins = n_bins, n_mesh =2,dremi = NULL,
                    plot = T, plot_only = FALSE, xlab = xx,ylab, limits = c(0, 0.15), scale.axis = TRUE)

  if (shuffle){
    # random shuffled gene x
    set.seed(235)
    Xr = sample(X[,1],length(X[,1]),replace = FALSE)
    dremi.rand = run_dremi(X = Xr, Y = Y,output.dir = file.path(output.dir,"dremi.heatmap",role),
                           k = k, n_bins = n_bins, n_mesh = 2,dremi = NULL,
                           plot = FALSE, plot_only = FALSE)
    
    res = data.frame(gene = xx, dremi = dremi[["dremi"]], role = role)%>%
      rbind(data.frame(gene = sprintf("shuffled %s", xx), dremi = dremi.rand[["dremi"]], role = "random"))
  }else{
    
    res = data.frame(gene = xx, dremi = dremi[["dremi"]], role = role)
  }
  return(res)
  
}

run_dremi_on_pathways = function(pathway_oi,Y,ylab,output.dir,data.magic,sobj,rm_pathway_genes,
                                 plot_hist = TRUE, plot_only = FALSE, plot_dremi = TRUE, output = NULL){
  
  # run dremi on all genes within a pathway of interests vs Y
  if (plot_dremi | !plot_only){
    log_info("Getting genes from pathway: ", pathway_oi)
    network = readRDS(file.path(nn.pathway.dir,pathway_oi,"network.rds"))
    attri= network$attributes%>%filter(!Role %in% c("target"))
    
    dremi_res_list = list()
    
    for (role in attri$Role%>%unique()){
    #  for (role in "receptor"){
      if (!dir.exists(file.path(output.dir,role))){dir.create(file.path(output.dir,"dremi.heatmap",role))}    
      
      log_info(" * Calculating DREMI on : ", role,"s")
      
      # get genes 
      gene_x_list = attri%>%filter(Role %in% role)%>%pull(Gene)%>%intersect(colnames(data.magic[["result"]]))
      
      
      dremi_res = lapply(gene_x_list, run_dremi_wrapper, data.magic = data.magic, plot_dremi = plot_dremi, role = role,output.dir = output.dir)%>%do.call("rbind",.)
      dremi_res$pathway =  pathway_oi 
      dremi_res_list[[role]] = dremi_res
    }
    
    # randomly select 100 genes that are not in the pathway
    # n_rnd_genes = 1000#length(attri$Gene)#1000
    # log_info(" * Calculating DREMI on : ", n_rnd_genes," random genes")
    # random_gene_list = sample(colnames(data.magic[["result"]]),n_rnd_genes)%>%setdiff(rm_pathway_genes)
    # dremi_res_r = lapply(random_gene_list, run_dremi_wrapper, shuffle=F,plot_dremi = F, role = "rndGene")%>%do.call("rbind",.)
    # tmp = dremi_res_r %>% filter(!role%in% "random")
    # tmp$role = "rndGene"
    # dremi_res_list[["rndGene"]] = tmp

    

    output = do.call("rbind",dremi_res_list)
  
  }
  
  saveRDS(output, file = file.path(output.dir, "dremi_output.rds"))

  
  return(output)
}

plot_dremi = function(output,output.dir){
  
  roles = c(unique(output$role))
  my_pal = c(`random` = "grey",`ligand` = "#F0E442",`associated_ligand` = "#F0E442", `receptor` = "#009E73" , `signalingMediator` = "#B016C1", `transcriptionFactor` =  "#FFA500",
             `rnd_ligand` = "#B2A84D", `rnd_receptor` = "#1F755A" , `rnd_signalingMediator` = "#60116B", `rnd_transcriptionFactor` =  "#B57002"
  )
  names(my_pal) =roles 
  if (is.null(output)){log_info("Please provide output to plot.")}
  p1 = ggplot(output, aes(x=dremi, fill=role, color = role)) +
    geom_histogram(alpha=0.5, position="identity",bins = 50)+
    scale_fill_manual(values = my_pal)+
    scale_color_manual(values = my_pal)+
    theme_classic()
  
  roles_to_plot = c("ligand","receptor","signalingMediator","transcriptionFactor")
  p2 = lapply(roles_to_plot, function(r){
    ggplot(output%>%filter(role%in%c("random",paste0("rnd_",r),r)), aes(x=dremi)) +
      #geom_histogram(aes(y=..density..), colour="black", fill="white",bins = 10)+
      geom_density(alpha=.5, aes(fill = role), color = "black")+ 
      scale_fill_manual(values = my_pal)+
      theme_classic()})%>%ggarrange(plotlist = ., nrow = 1)
  
  p3 =ggplot(output%>%mutate(role = case_when(role %in% roles_to_plot ~ "pathway", TRUE ~ role))%>%filter(!role%in%"random"), aes(x=dremi)) +
    #geom_histogram(aes(y=..density..), colour="black", fill="white",bins = 10)+
    geom_density(alpha=.5, aes(fill = role), color = "black")+ 
    theme_classic()
  
  p4=lapply(roles_to_plot, function(r){
    ggplot(output%>%filter(role%in%c("random",r)), aes(x=dremi)) +
      geom_density(alpha=.5, aes(fill = role), color = "black")+ 
      scale_fill_manual(values = my_pal)+
      geom_segment(data = output%>%filter(role%in%c(r)), aes(x = dremi, y= 0, yend = -2.5, xend = dremi), color = 5, lwd = 1)+
      ggtitle(r)+
      theme_classic()})%>%ggarrange(plotlist = ., nrow = 1)
  
  
  ggsave(ggarrange(plot_grid(p1,p3, rel_widths = c(0.3, 0.7),nrow = 1), p2,p4, nrow = 3, ncol = 1), 
         filename = file.path(output.dir,"dremi_distr_histogram.pdf"),width = 4*length(roles_to_plot), height = 9) 
  
}
################################################################################|
# RUN ---------------
################################################################################|
#sobj.sm = readRDS(file.path(magic_output_path, "sobj_after_magic_add_module.rds"))
## Set up -------------------------------------------------------------------
final_res <- readRDS(file.path(magic_output_path, "magic_output.rds"))

pathways_list = list.files(nn.pathway.dir)
pathway_oi = pathways_list[1]

## Run DREMI on MAGIC-imputed module scores--------
for (assay_oi in c("magic")){

    
  output.dir = file.path(dremi.output.dir,sprintf("dremi_output_n=25000.%s",assay_oi))
  data.magic = final_res[["data_MAGIC"]]
  Y = sobj.sm@assays[[paste0(assay_oi,".module.seurat")]]@data%>%t()
  ylab =  colnames(Y)
  
  ### random genes from each categories --------
  
  rnd_genelist = annotations%>%filter(expressed ==T & !gene %in% c(rm_pathway_genes,targets)) # all genes expressed in astrocytes (>1%) but not in the pathways or targets

  dremi_res_r = lapply(unique(rnd_genelist$role), function(r){
    random_gene_list = rnd_genelist%>%filter(role ==r )%>%pull(gene)
    n_rnd_genes = length(random_gene_list )
    log_info(" * Calculating DREMI on : ", n_rnd_genes," random ", r, " genes ...")
    dremi_res_r = lapply(random_gene_list, run_dremi_wrapper, data.magic = data.magic, shuffle=F,plot_dremi = F, role = paste0("rnd_",r),output.dir = output.dir)%>%do.call("rbind",.)%>%
      filter(!role%in% "random")
    dremi_res_r$role =paste0("rnd_",r)
    return(dremi_res_r)
  })%>%do.call("rbind",.)
  dremi_res_r$pathway = "random"
  saveRDS(dremi_res_r, file = file.path(output.dir,"dremi_random_genes.rds"))

  ### pathway genes --------
  res.list = lapply( pathway_oi, function(p,output.dir){
    print(p)
    res = run_dremi_on_pathways(p,Y,ylab, file.path(output.dir, p),data.magic,sobj.sm,rm_pathway_genes = rm_pathway_genes,
                                plot_hist = TRUE, plot_only = F, plot_dremi = F)
    
    res = rbind(res, dremi_res_r)
    dir.create(file.path(output.dir, p), recursive = T)
    saveRDS(res, file = file.path(output.dir, p,"dremi_output.rds"))
    
    plot_dremi(res,file.path(output.dir, p))
    return(res)
  },output.dir)

}
###############################################################################|
full_pathway = "SEMA4B.JAM2.NRG3.EFEMP1.LRIG1.FN1.NCAN.LRRC4B.SEMA4D.RGMA.LGALS3.GNAS.NRXN1"
full_network <- readRDS(file.path(nn.pathway.dir,full_pathway,"network.rds"))
full_attri= full_network$attributes

ligands_oi = str_split(full_pathway, "\\.")%>%unlist

attri_assignment = lapply(ligands_oi, function(x){
  attri = readRDS(file.path(nn.pathway.dir,x,"network.rds"))[["attributes"]]
  out = tibble(full_attri$Gene%in% attri$Gene)%>%`colnames<-`(x)
  #out$gene = full_attri$Gene
  return(out)
})%>%do.call("cbind",.)
attri_assignment$sum = rowSums(attri_assignment)
attri_assignment$Gene = full_attri$Gene
attri_assignment = merge(attri_assignment, full_attri, by = "Gene")

saveRDS(attri_assignment, file.path(dremi.output.dir,"pathway_assignment.rds"))

rm(list = ls())
gc()
setwd("5.Spatial_Analysis/mistyr")

library(tidyverse)
library(dplyr)
library(colorspace)
library(gplots)
library(ggpubr)
library(cowplot)
library(Seurat)
library(ComplexHeatmap)
################################################################################|
# Input -------
################################################################################|
results_folder <- "output"
view_oi = "juxta_480"
suffix = "_SCT"
deg_thres = 0 #-Inf
rm_targets = c("AC016745.1", "AC002069.2", "AL512662.2", "AC006148.1") # genes not found in original snRNAseq data


################################################################################|
# Define genes of interest -------
################################################################################|

rm_projid = NULL

## Read DEGs -------------------------------------------------
receiver = "Ast.10"
### Markers from DE analysis -----
de_output_folder = "ast_mic_subtype_signature_genes_summary/"
subtype.cutoff = "adjp=0.05_logfc=1.00"
de_exp <- read.csv(sprintf("%s/%s_de_genes_summary_%s.csv",de_output_folder,receiver,subtype.cutoff))
markers_oi = de_exp%>%filter(avg.logFC>deg_thres)%>%pull(gene)%>%setdiff(rm_targets)


## Read NicheNet-inferred ligands and receptors---------------
nn_output_folder = "../../2.NicheNet Analysis/nn_output"
prioritized_table = read_csv(file.path(nn_output_folder, receiver, "prioritized_table.csv"))

all_ligands = prioritized_table%>%pull(ligand)%>%unique()
all_receptors = prioritized_table%>%pull(receptor)%>%unique()

top_n = 100
nn_ligands = unique(prioritized_table%>%slice(1:top_n)%>%pull(ligand))
nn_receptors = unique(prioritized_table%>%slice(1:top_n)%>%pull(receptor))

## Use top vip ligands ----------------------------------------

plsr_output_folder = "../../3.PLSR/plsr_analysis"
plsr_res <- readRDS(file.path(plsr_output_folder,"plsr_res_response=Ast.10.rds"))
df_vip = plsr_res$vip_scores
df_vip$ligand = stringr::str_split(df_vip$x_var,pattern = "\\ \\(", simplify = TRUE)[,1]
df_vip = filter(df_vip, !(ligand=="NRG3" & vip < 1)) #NRG3 has both positive and negative ligand summary score, here we keep the one with significant vip

pls_ligands = pull(df_vip, ligand)
vip_ligands = filter(df_vip, vip >=1)%>%pull(ligand)
top_vip_ligands = filter(df_vip, vip >=1.1)%>%pull(ligand)

## Create Predictor annotations ------
pred_annot = tibble(Predictor = all_ligands)%>%
  mutate(ligand = Predictor)%>%
  mutate(nn_ligand = case_when(Predictor%in%nn_ligands ~ 1, T ~ 0))%>%
  mutate(pls_ligand = case_when(Predictor%in%pls_ligands ~ 1, T ~ 0))%>%
  mutate(vip_ligand = case_when(Predictor%in%vip_ligands ~ 1, T ~ 0))%>%
  mutate(top_vip_ligand = case_when(Predictor%in%top_vip_ligands ~ 1, T ~ 0))%>%
  merge(df_vip, by = "ligand",all = T)
################################################################################|
# Get Results -------
################################################################################|

## Results ---

slide_files <- list.files(results_folder)%>%.[grepl(suffix,.)]
slide_ids <- gsub(suffix, "", slide_files)
rm_slides = lapply(rm_projid, function(id){slide_ids[str_starts(slide_ids, id)]})%>%unlist%>%
  paste0(., suffix)



dir.create(file.path(results_folder,view_oi), recursive = T)
################################################################################|
# Collect importance score --------
################################################################################|

collect_importance = function(results_folder, view_oi, .rm_slides=NULL){
  .rm_slides = file.path(results_folder, .rm_slides)
  slide_files = list.dirs(results_folder,recursive = FALSE)%>%.[grepl(suffix,.)]%>%setdiff(.rm_slides)
  logger::log_info("number of slides: ",length(slide_files))
  lapply(slide_files, function(slide_file){
    misty_res_slide <- readRDS(file.path(slide_file,"misty_res_slide.rds"))
    res = misty_res_slide[["importances"]]
    df = res%>%filter(view%in%view_oi )
    df$sample = slide_file
    return(df)
  })%>%do.call("rbind",.)
  
  
}

# All importance scores from all samples
imp.df = collect_importance(results_folder, view_oi,.rm_slides=rm_slides)
# Summarised importance scores across samples
summary.imp.df = imp.df%>%group_by(Predictor, Target)%>%
  summarise(mean_imp = mean(Importance, na.rm = T))%>% # average across samples
  merge( pred_annot, by = "Predictor", all.x = T)%>%
  mutate(color_bar = case_when(vip >=1 ~ "blue",vip < 1 ~ "lightblue",nn_ligand ==1 ~ "orange",TRUE ~ "white"))%>%
  mutate(predictor_prior= case_when(vip >=1 ~ "vip ligand",
                                    vip < 1 ~ "pls ligand",
                                    nn_ligand ==1 ~ "nn ligand",
                                    TRUE ~ "n.s."
                                    
  ))

# Summarised importance scores on the targets of interests
summary.imp.df.f = filter(summary.imp.df)%>%filter(Target %in% markers_oi)

# Calculate P-val
summary.imp.df.f  = summary.imp.df.f%>%mutate(rank = rank(-mean_imp))%>%
  mutate(max_rank = length(rank))
################################################################################|
# Figure 5D: Averaged Importance Heatmap of Predictors_oi --------
################################################################################|

plot_hm_importance = function(.df, var_oi = "mean_imp", results_folder){
  p = ggplot(.df, aes(y = Predictor, x = Target, fill = .data[[var_oi]]))+
    geom_tile(color = "black")+
    scale_fill_gradient2(high = "blue" , low = "white", mid = "white", limits = c(-1,1),oob = scales::squish)+
    theme(axis.text.x = element_text(angle = 45, hjust=1),text = element_text(size = 12))
  ggsave(filename = file.path(results_folder,sprintf("hm_%s.pdf",var_oi)), width= length(unique(.df$Target))*0.5, height= length(predictor_oi)*0.36, limitsize = F)
  return(p)
}

predictor_oi = pls_ligands
summary.imp.df.f.selected_ligands = filter(summary.imp.df.f, Predictor %in% predictor_oi)
hm_mean = plot_hm_importance(summary.imp.df.f.selected_ligands, var_oi = "mean_imp", results_folder = file.path(results_folder, view_oi))
hm_mean

summary.imp.df.f.selected_ligands%>%saveRDS(., file = file.path(results_folder,view_oi, "summary_importance_pls_ligands.rds"))

################################################################################|
# Figure 5E: Permutation test--------
################################################################################|

var_oi = "mean_imp"

# summarize the mean importance scores for each gene by summing up 
# the mean importance across ast10 targets
var_oi = "mean_imp"
gene_summary = summary.imp.df.f%>%reshape2::dcast(Predictor ~ Target ,value.var = var_oi)%>%
  column_to_rownames("Predictor")%>%
  mutate(across(colnames(.), function(df) (df-mean(df, na.rm = T))/sd(df, na.rm = T)))%>%
  rowSums(na.rm = T)%>%as.data.frame%>%`colnames<-`(paste0("sum_",var_oi))%>%
  rownames_to_column("Predictor")%>%
  merge( pred_annot, by = "Predictor", all.x = T)%>%
  filter(!is.na(.data[[paste0("sum_",var_oi)]]))
gene_summary[is.na(gene_summary)]= 0

gene_summary = gene_summary%>%mutate(rank = rank(-sum_mean_imp))

# randomly select ligand set of the same size as the set of vip ligands (abs(vip)>=1)
set.seed(5653)

n_iter = 1e5
n = sum(gene_summary$vip_ligand)
perm_imp = function(prior, gene_summary, n_iter){
  n = sum(gene_summary[[prior]]) # number of ligands in the set
  #permute
  rnd_mean = lapply(c(1:n_iter), function(s){
    mean(sample(gene_summary$sum_mean_imp, n))
  })%>%do.call("c",.)
  
  
  rnd_mean.df = tibble(iter = c(1:n_iter), sum_mean_imp = rnd_mean )
  
  prior.mean = filter(gene_summary, .data[[prior]] == 1)%>%pull(sum_mean_imp)%>%mean
  prior.p = sum(rnd_mean.df$sum_mean_imp>prior.mean)/n_iter
  
  p =  ggplot(rnd_mean.df, aes(x = sum_mean_imp)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 1, fill = "white") +
    geom_density(lwd = 1, colour = "#90AD1C",
                 fill = "#E0E75C", alpha = 0.25)+
    #geom_text(x = prior.mean, y = 0.03, label = sprintf("P = %0.1e",prior.p ), color = "#D81C5C")+
    geom_vline(xintercept = prior.mean, color = "#747BFF", linewidth = 1.5)+
    xlab("Mean total importance")+
    ylab("Density")+
    ggtitle(prior)+
    theme_classic()
  
  return(list("results" =  rnd_mean.df, "plot" = p ))
}

# ligands with different prioritization levels
p_vip = perm_imp("vip_ligand", gene_summary, n_iter)
p_pls = perm_imp("pls_ligand", gene_summary, n_iter)
p_nn = perm_imp("nn_ligand", gene_summary, n_iter)

perm_fig = plot_grid(p_vip[["plot"]],p_pls[["plot"]], p_nn[["plot"]], nrow = 1, ncol = 3)
ggsave(plot = perm_fig,filename =  file.path(results_folder,view_oi,"permutation_emp_p.png"),dpi = 350, width = 6, height = 3)
ggsave(plot = p_vip[["plot"]],filename =  file.path(results_folder,view_oi,"permutation_emp_p_vip.png"), dpi = 350,width = 2, height = 1.8)

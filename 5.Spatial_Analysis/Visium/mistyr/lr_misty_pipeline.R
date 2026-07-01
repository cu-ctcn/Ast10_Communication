# ** Adapted from Ricardo O. Ramirez Flores **
#------------------------------------------------#|
#     MISTy analysis to predict Ast10 markers
#         based on local ligand expression
#------------------------------------------------#|
#           Natacha Comandante-Lou 
#             (nc3018@columbia.edu)
#------------------------------------------------#|


library(tidyverse)
library(Seurat)
future::plan(future::multisession) 
library(mistyR)
library(logger)
source("../utils/utils_mistyr.R")
source("../../utils/seurat_get_summary_scores.R")
future::plan(future::multisession)
devtools::session_info()
library(OmnipathR)

################################################################################|
# MISTy Definition -------
################################################################################|

# Pipeline definition:
run_colocalization <- function(slide,
                               view_types, 
                               view_assays, 
                               view_features, 
                               view_params,
                               out_label,
                               misty_out_alias = "results/lr_") {
  
 
  # Define spatial context of each view -----

  misty_out <- paste0(misty_out_alias, 
                      out_label, "_", view_assays[["main"]])
  
  run_misty_seurat(visium.slide = slide,
                   view.assays = view_assays,
                   view.features = view_features,
                   view.types = view_types,
                   view.params = view_params,
                   spot.ids = NULL,
                   out.alias = misty_out)
  
  return(misty_out)
}




run_misty_pipeline = function(slide_file,view_types, view_assays, view_features,view_params,
                              score_module = F,
                              out.dir = "mistyr_results/"){
  
  if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}
  dir.create(paste0(out.dir,"logs"), recursive = T)
  #############################################################################|
  # Log -------
  #############################################################################|
  log_layout(layout_simple)
  log_file = file.path(paste0(out.dir,"logs"), paste0(gsub("[.]rds", "", slide_file),"_mistyr_pipeline.log"))
  if(file.exists(log_file)){
    log_warn("Overwriting previous log: ", log_file)
    file.remove(log_file)
  }
  
  log_appender(appender_tee(file = log_file))
  
  logger::log_info("###########################################################")
  logger::log_info("#               Running MistyR on ",slide_file,"          #")
  logger::log_info("###########################################################")
  
  logger::log_info("* View Assays: ")
  lapply(names(view_types), function(v){
    logger::log_info("    ","-", v,": ",view_assays[[v]])
    
  })

  # Read spatial transcriptomics data
  slide_id <- gsub("[.]rds", "", slide_file)
  slide <- readRDS(paste0(slide_files_folder, slide_file))
  
  
  
  ## Normalization --------------
  DefaultAssay(slide) = "Spatial"

  ## sctransform
  slide =  SCTransform(slide, assay = "Spatial", vst.flavor = "v2",
                       verbose = T,vars.to.regress = c("percent.mt","NumNuclei"))
  #remove genes with zero expression:
  keep_genes = which(rowSums(slide@assays[["SCT"]]@data>0)>0)%>%names
  slide <- subset(slide, features = keep_genes)

  logger::log_info("Assays available: ", paste(slide@assays%>%names(),collapse = ", "))

  # Summarize score for recevier
  if (score_module){
    slide= seurat_get_summary_scores(slide,
                                     view_features,
                                     module.key = "module",
                                     method = "seurat",
                                     assay = "SCT",
                                     out.fig.dir = file.path("modules.score"),
                                     plot = F,
                                     plot.col.limits = NULL,
                                     plot.reduction =  "umap" )
    view_features[["main"]] = "main"
  }
  
  
  
  logger::log_info("* Saving to ...",file.path(out.dir,"objects"))
  dir.create(file.path(out.dir,"objects"), recursive = T)
  saveRDS(slide, file.path(out.dir,"objects",sprintf("%s_slide.rds", slide@meta.data$orig.ident%>%unique) ))
  
  logger::log_info("* Running MistyR: ")
  DefaultAssay(slide) <- "SCT"
  view_features = lapply(view_features, function(x){rownames(slide)%>%intersect(x)})
  print(view_features)
  
  
  mout <- run_colocalization(slide = slide,
                             view_types = view_types,
                             view_assays = view_assays,
                             view_features = view_features,
                             view_params = view_params,
                             out_label = slide_id,
                             misty_out_alias = out.dir)
  
  misty_res_slide <- collect_results(mout)
  

  
  
  saveRDS(misty_res_slide,  paste0(mout,"/misty_res_slide.rds"))
  
  logger::log_info("* Plotting: ")
  
  plot_folder <- paste0(mout, "/plots")
  
  dir.create(plot_folder, recursive = T)
  
  pdf(file = paste0(plot_folder, "/", slide_id, "_", "summary_plots.pdf"), width = 15, height = 20)
  
  mistyR::plot_improvement_stats(misty_res_slide)
  mistyR::plot_view_contributions(misty_res_slide)
  
  # plot interaction heatmap
  lapply(view_types, function(v){
    if (v=="intra"){mistyR::plot_interaction_heatmap(misty_res_slide, "intra", cutoff = 0)}else{
      mistyR::plot_interaction_heatmap(misty_res_slide, paste0(v,"_",view_params[[v]]), cutoff = 0)
    }
    
  })
  
  dev.off()
  
  return(mout)
  
}
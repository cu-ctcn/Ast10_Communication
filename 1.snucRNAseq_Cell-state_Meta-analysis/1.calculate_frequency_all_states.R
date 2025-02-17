#-----------------------------------------#|
# Calculate Cell-State Frequency
#-----------------------------------------#|
#       Natacha Comandante-Lou 
#       (nc3018@columbia.edu)
#-----------------------------------------#|
rm(list = ls())
setwd("/Code/1.snucRNAseq_Cell-state_Meta-analysis")
source('../utils/plsr_helper_func.R')


library(plyr) 
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(scCustomize)
library(ggpubr)

output_folder = "output"
dir.create(output_folder)

cell_metadata_file = "cell-state_freq/cell_annotations" #contains data.frames of cell-state annotations for each cell in the four snRNAseq data
cell_metadata_list = lapply(list.files(cell_metadata_file),function(m){readRDS(m)})
#################################################################################|
# Get freq -----------
#################################################################################|
get_freq_meta= function(cell_meta, freq_df= NULL, save_file_name, meta_file, qc_meta_file, study_name, group.by = "predicted.state"){
  # Calculate cell-state frequency based on cell-state annotations
  freq_df = GetStateFreq_v4(meta =cell_meta, group.by = group.by, obs.key = "projid",h5ad = F,
                            logfreq = FALSE, output.dir = file.path(output_folder,"freq_rds",study_name), save=TRUE)
  # Standardize projid format to 8 digit
  freq_df$projid = as.character(freq_df$projid)%>%stringr::str_pad(.,8,pad = "0")
  rownames(freq_df) = freq_df$projid
  freq_df$reference = study_name

  meta_file$projid = as.character(meta_file$projid)%>%stringr::str_pad(.,8,pad = "0")
  qc_meta_file$projid = as.character(qc_meta_file$projid)%>%stringr::str_pad(.,8,pad = "0")
  

  # Assemble metadata
  meta.data = filter(meta_file, projid%in%rownames(freq_df))%>%
    arrange(match(projid,  rownames(freq_df)))%>%
    select(var_oi)%>%
    mutate(sqrt.tangles=sqrt(tangles))%>%
    mutate(sqrt.amyloid=sqrt(amyloid))%>%
    mutate(reference = study_name)%>%
    merge(qc_meta_file, all.x = T, by = "projid")%>%
    mutate(nhw.group = case_when(race7==1 & spanish==2 ~ "NHW", TRUE ~ "nonNHW"))
  
  rownames(meta.data) = meta.data$projid
  
  # Square-root transform of frequency
  freq_df.s = sqrt(freq_df%>%select(-c("projid","reference")))
  freq_df.s$projid = rownames(freq_df.s)
  freq_df.s$reference = study_name
  rownames(freq_df.s) = rownames(freq_df)
  
  # Log10 transform of frequency
  freq_df.l = log(freq_df%>%select(-c("projid","reference"))*100+1)
  freq_df.l$projid = rownames(freq_df.l)
  freq_df.l$reference = study_name
  rownames(freq_df.l) = rownames(freq_df)
  
  output = list(meta.data = meta.data, freq = freq_df, sqrt_freq = freq_df.s, log_freq = freq_df.l)
  saveRDS(output,save_file_name)
  return(output)
}

dir.create(file.path(output_folder,"output_rds"), recursive = T)

out_mit = get_freq_meta(cell_metadata_list[["MIT"]] , save_file_name = file.path(output_folder ,"output_rds","mit_state_freq.rds"), meta_file=all_meta, study_name = "MIT", qc_meta_file = qc.mit.meta )
out_cuimc2= get_freq_meta(cell_metadata_list[["CUIMC2"]], save_file_name = file.path(output_folder ,"output_rds","cuimc2_state_freq.rds"), meta_file=all_meta, study_name = "CUIMC2", qc_meta_file = qc.cuimc2.meta)
out_diversity= get_freq_meta(cell_metadata_list[["Diversity-All"]], save_file_name = file.path(output_folder ,"output_rds","diversity_state_freq.rds"), meta_file=all_meta, study_name = "Diversity", qc_meta_file = qc.diversity.meta)
out_discovery = get_freq_meta(cell_metadata , save_file_name = file.path(output_folder ,"output_rds","discovery_state_freq.rds"), meta_file=all_meta, study_name = "Discovery", qc_meta_file = qc.discovery.meta, group.by = "state")


####################################
# Split Diversity --------
####################################
# split Diversity set by NHW (removed from analysis) or nonNHW
split_metadata = function(split.by, .df, group_oi,save_file_name, rename.ref = T){
  projid_oi = filter(.df[["meta.data"]], .data[[split.by]]==group_oi)%>%pull(projid)
  out = lapply(.df, function(x){ filter(x, projid%in%projid_oi)})
  if (rename.ref){
    out = lapply(out, function(x){ x = mutate(x, reference = paste(reference, group_oi, sep = "-"));return(x)})
  }
  saveRDS(out,save_file_name)
  return(out)
}

out_diversity_nhw = split_metadata(split.by = "nhw.group", .df = out_diversity, group_oi = "NHW",
                                   save_file_name = file.path(output_folder,"output_rds", "diversity_state_freq_nhw.rds"))
out_diversity_other= split_metadata(split.by = "nhw.group", .df = out_diversity, group_oi = "nonNHW", 
                                    save_file_name = file.path(output_folder, "output_rds", "diversity_state_freq_nonNHW.rds"))

####################################
# Merge -------------
####################################

# all data
out_merge = lapply(names(out_mit), function(x){ plyr::rbind.fill(out_discovery[[x]],
                                                                 out_mit[[x]],
                                                                 out_cuimc2[[x]],
                                                                 out_diversity_nhw[[x]],
                                                                 out_diversity_other[[x]])})

names(out_merge) = c("meta.data" ,"freq","sqrt_freq","log_freq" )
saveRDS(out_merge,file.path(output_folder,"snuc_state_freq_all_data.rds"))

# data analyzed (removed overlap)
out_merge_analyzed = lapply(names(out_mit), function(x){ 
  df = plyr::rbind.fill(out_discovery[[x]],out_mit[[x]],out_cuimc2[[x]],out_diversity_other[[x]])%>%
    rename_study%>%analyzed_donors
  return(df)
})
names(out_merge_analyzed) = c("meta.data" ,"freq","sqrt_freq","log_freq" )
saveRDS(out_merge_analyzed,file.path(output_folder,"snuc_state_freq.rds"))

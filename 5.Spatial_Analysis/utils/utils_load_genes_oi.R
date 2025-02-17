#------------------------------------------------#|
#           Get genes of interests
#------------------------------------------------#|
#           Natacha Comandante-Lou 
#             (nc3018@columbia.edu)
#------------------------------------------------#|
# Genes of interests: 
# - all ligands, all receptors (from existin nichenet and omnipath resources)
# - nichenet ligands, receptors (from nichenet prioritization, figure 2)
# - plsr ligands, receptors (from ligands in plsr model, figure 3)
# - vip ligands, receptors (from significant ligands in plsr model, figure 3)
# - receiver marker genes (from differential expression analysis)
library(tidyverse)
library(OmnipathR)


load_genes_oi = function(nn_output_folder, plsr_output_folder, receiver){
  prioritized_table = read_csv(file.path(nn_output_folder, "prioritized_table.csv"))
  
  # Get All ligands
  
  ## ligands from NicheNet
  all_ligands_nn = prioritized_table%>%pull(ligand)%>%unique()
  
  ## ligands from Omnipath
  lig_rec <- OmnipathR::import_intercell_network(interactions_param = list(datasets = c('ligrecextra', 'omnipath', 'pathwayextra')),
                                                 transmitter_param = list(parent = 'ligand'),
                                                 receiver_param = list(parent = 'receptor'))
  all_ligands <- unique(lig_rec$source_genesymbol)%>%c(.,all_ligands_nn)%>%unique
  
  # Get All receptors
  all_receptors = prioritized_table%>%pull(receptor)%>%unique()
  
  # Get Prioritized ligands receptors
  top_n = 100
  nn_ligands = unique(prioritized_table%>%slice(1:top_n)%>%pull(ligand))
  nn_receptors = unique(prioritized_table%>%slice(1:top_n)%>%pull(receptor))
  
  ## Get top vip ligands ----------------------------------------
  plsr_res <- readRDS(file.path( plsr_output_folder,"ligand_model","plsr_res_response=Ast.10.rds"))
  df_vip = plsr_res$vip_scores
  df_vip$ligand = stringr::str_split(df_vip$x_var,pattern = "\\ \\(", simplify = TRUE)[,1]
  
  pls_ligands = pull(df_vip, ligand)
  vip_ligands = filter(df_vip, vip >=1)%>%pull(ligand) # significant ligands
  
  ### Markers from DE analysis -----
  de_output_folder = "ast_mic_subtype_signature_genes_summary/"
  subtype.cutoff = "adjp=0.05_logfc=1.00"
  de_exp <- read.csv(sprintf("%s/%s_de_genes_summary_%s.csv",de_output_folder,receiver,subtype.cutoff))
  markers_oi = de_exp%>%pull(gene)
  
  return(list("all_ligands" = all_ligands, "all_receptors" = all_receptors, 
               "nn_ligands" =  nn_ligands,  "nn_receptors" = nn_receptors,
              "pls_ligands" = pls_ligands, "vip_ligands" = vip_ligands, "markers_oi" = markers_oi))
}

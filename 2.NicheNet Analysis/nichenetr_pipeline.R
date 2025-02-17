#' NicheNet Pipeline using nichenetr v2
require(nichenetr)
require(SeuratDisk)
require(Seurat) 
require(tidyverse)
require(dplyr)
require(logger)
require(cowplot)

source('2.NicheNet Analysis/utils/nichenet_utils.R') # plot functions
run_nichenet_pipeline = function(seuratObj,
                                 receiver = "Ast.10",
                                 sender,
                                 
                                run_nichenet = TRUE,
                                 expressed_genes_receiver = NULL,
                                 expressed_genes_sender = NULL,
                                 gene_set_oi,
                                 pct_thres = 0.1, # threshold for defining expressed genes
                                 
                                 # Prioritization weights
                                prioritize = TRUE,
                                  w.de_ligand = 1,
                                  w.de_receptor = 1,
                                  w.activity_scaled = 2,
                                  w.exprs_ligand = 1,
                                  w.exprs_receptor = 1,
                                  w.ligand_condition_specificity = 0,
                                  w.receptor_condition_specificity = 0,
                                  celltype_colname = "state",
                                  DE_table,
                                  expression_info,
                  
                                # Plotting
                                  plot.prior.table = NULL,
                                  plot.sender,
                                  top_n = 100,
                                
                                 figure_folder = "nichenet_figures",
                                 output_folder = "nichenet_output"

  
  ){
  if(!dir.exists(output_folder)){dir.create(output_folder,recursive = TRUE)}
  if(!dir.exists(figure_folder)){dir.create(figure_folder,recursive = TRUE)}
  #############################################################################|
  # LOG -------
  #############################################################################|
  log_layout(layout_simple)
  log_file = file.path(output_folder , sprintf("nichenetr_pipeline_%s.log", receiver))
  if(file.exists(log_file)){
    log_warn("Overwriting previous log: ", log_file," \n")
    file.remove(log_file)
  }
  
  log_appender(appender_tee(file = log_file))
  
  if (run_nichenet){
    #############################################################################|
    # LOAD NICHENET DATABASE (V2) -----------
    #############################################################################|
    log_info("Default Assay: ",DefaultAssay(seuratObj))
    log_info("Number of participants: ", seuratObj@meta.data$projid%>%unique%>%length())
    # NicheNet’s ligand-receptor data sources
    lr_network_path = "/mnt/mfs/ctcn/team/natacha/NicheNet_Database/7074291/lr_network_human_21122021.rds"
    ligand_target_matrix_path = "/mnt/mfs/ctcn/team/natacha/NicheNet_Database/7074291/ligand_target_matrix_nsga2r_final.rds"
    weighted_networks_path = "/mnt/mfs/ctcn/team/natacha/NicheNet_Database/7074291/weighted_networks_nsga2r_final.rds"
    
    log_info("Reading NicheNet Data Sources from: ")
    log_info("* lr_network: ", lr_network_path )
    log_info("* ligand_target_matrix: ", ligand_target_matrix_path)
    log_info("* weighted_networks: ", weighted_networks_path)
    log_info("\n")
    lr_network = readRDS(lr_network_path)
    lr_network = lr_network %>% dplyr::distinct(from, to)
    
    # Ligand-target model: This model denotes the prior potential that a particular 
    # ligand might regulate the expression of a specific target gene.
    ligand_target_matrix = readRDS(ligand_target_matrix_path) # target genes in rows, ligands in columns
    
    # Weighted integrated network
    weighted_networks = readRDS(weighted_networks_path)
    weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))
    
    DefaultAssay(seuratObj) = "RNA"
    ##############################################################################|
    # 1. DEFINE SENDERS & RECEIVER -------------
    ##############################################################################|
    log_info("1. Define receiver and senders")
    log_info("* Receiver: ",receiver)
    log_info("* Sender: ", paste(sender, collapse=", "),"\n")
    log_info("\t")
    ## Receiver ----------------
    # Receiver expressed genes
    if (is.null(expressed_genes_receiver)){
      log_info("Getting genes expressed in >=", pct_thres*100, "% of receiver cells...")
      expressed_genes_receiver = get_expressed_genes(receiver, seuratObj, pct = pct_thres,assay_oi = "RNA")
      
      if(!dir.exists(file.path( output_folder, "expressed_genes"))){dir.create(file.path( output_folder, "expressed_genes"), recursive = TRUE)}
      saveRDS(expressed_genes_receiver, file = file.path( output_folder, "expressed_genes", sprintf("expressed_genes_%s_pct=%0.1f.rds",receiver,pct_thres)))
  
    }
    
    # Receiver background genes
    background_expressed_genes = rownames(ligand_target_matrix)
    
    ## Sender --------------------
    Idents(seuratObj) = "state"
  
    # Sender expressed genes
    if (is.null(expressed_genes_sender)){
      log_info("Getting genes expressed in >=", pct_thres*100, "% of sender cells...")
      expressed_genes_sender = sender %>% unique() %>% lapply(get_expressed_genes, seuratObj, pct_thres,assay_oi = "RNA") %>% unlist() %>% unique()
      
      if(!dir.exists(file.path( output_folder, "expressed_genes"))){dir.create(file.path( output_folder, "expressed_genes"), recursive = TRUE)}
      saveRDS(expressed_genes_sender, file = file.path( output_folder, "expressed_genes", sprintf("expressed_genes_senders_pct=%0.1f.rds",pct_thres)))
      
    }
    ##############################################################################|
    # 2. DEFINE GENESET OF INTEREST-------------
    ##############################################################################|
    log_info("2. Define geneset of interests")
  
    # Gene set of interests: cell-state specific signature genes (genes significantly 
    #  de-expressed only in receiver) & exist in NicheNet database
    geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
    
    log_info("* Number of genes in the geneset of interests: ", length(geneset_oi),"\n")
    log_info("\t")
    ##############################################################################|
    # 3. DEFINE POTENTIAL LIGANDS-------------
    ##############################################################################|
    log_info("3. Define potential ligands")
    
    # Potential ligands: these are ligands that are expressed by the “sender/niche” 
    # cell population and bind a (putative) receptor expressed by the “receiver/target” population
    
    ligands = lr_network %>% pull(from) %>% unique()
    receptors = lr_network %>% pull(to) %>% unique()
    
    expressed_ligands = intersect(ligands,expressed_genes_sender)
    expressed_receptors = intersect(receptors,expressed_genes_receiver)
    
    potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
    
    log_info("* Number of potential ligands: ", length(potential_ligands))
    log_info("\t")
    ##############################################################################|
    # 4. PERFORM NICHENET LIGAND ACTIVITY ANLAYSIS------------
    ##############################################################################|
    log_info("4. Perform NicheNet ligand activity analysis\n")
    ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
                                                  background_expressed_genes = background_expressed_genes, 
                                                  ligand_target_matrix = ligand_target_matrix, 
                                                  potential_ligands = potential_ligands)
    
    ligand_activities = ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected))) # AUPR was the most informative measure in NN v2
    log_info("\t")
    ##############################################################################|
    # 5. PRIORITIZATION------------
    ##############################################################################|
    if (prioritize){
      log_info("5. PRIORITIZATION")
      #Default weights
      prioritizing_weights = c("de_ligand" = w.de_ligand,
                               "de_receptor" = w.de_receptor,
                               "activity_scaled" = w.activity_scaled,
                               "exprs_ligand" = w.exprs_ligand,
                               "exprs_receptor" = w.exprs_receptor,
                               "ligand_condition_specificity" = w.ligand_condition_specificity,
                               "receptor_condition_specificity" = w.receptor_condition_specificity)
      
      log_info("Prioritization weights: ")
      log_info(" * Scaled ligand activity: activity_scaled = ",w.activity_scaled)
      log_info(" * Upregulation of the ligand in a sender cell type compared to other cell types: de_ligand = ", w.de_ligand)
      log_info(" * Upregulation of the receptor in a receiver cell type: de_receptor = ", w.de_receptor)
      log_info(" * Average expression of the ligand in the sender cell type: exprs_ligand = ",w.exprs_ligand)
      log_info(" * Average expression of the receptor in the receiver cell type: exprs_receptor = ", w.exprs_receptor)
      log_info(" * Condition-specificity of the ligand across all cell types: ligand_condition_specificity = ",w.ligand_condition_specificity)
      log_info(" * Condition-specificity of the receptor across all cell types: receptor_condition_specificity = ",w.receptor_condition_specificity)
      
      states <- unique(seuratObj$state)
      lr_network_renamed <- lr_network %>% rename(ligand=from, receptor=to)
      
      
      # ! Takes a long time to run, save for future reload ! --
      # Only calculate DE with genes that are in the ligand-receptor network
      if (is.null(DE_table)){
        DE_table <- calculate_de(seuratObj, celltype_colname = celltype_colname,
                                 features = union(expressed_ligands, expressed_receptors))
        saveRDS(DE_table,file = file.path(output_folder,"lig_rec_de_table.rds"))
      }
      
      if (is.null(expression_info)){
        # Average expression information
        expression_info <- get_exprs_avg(seuratObj, celltype_colname)
        saveRDS(expression_info,file = file.path(output_folder,"expression_info.rds"))
      }
    
      processed_DE_table <- process_table_to_ic(DE_table, table_type = "celltype_DE", lr_network_renamed,
                                                senders_oi = sender, receivers_oi = receiver)
      
      processed_expr_table <- process_table_to_ic(expression_info, table_type = "expression", lr_network_renamed)
      
      prior_table <- generate_prioritization_tables(processed_expr_table,
                                                    processed_DE_table,
                                                    ligand_activities,
                                                    prioritizing_weights = prioritizing_weights)
      
      saveRDS(prior_table, file = file.path(output_folder,"prioritized_table.rds"))
      write_csv(prior_table,file = file.path(output_folder,"prioritized_table.csv"))
      
    }
    ##############################################################################|
    # 6. OUTPUT -----
    ##############################################################################|
    
    nn_output = list(
      "receiver" = receiver,
      "sender" = sender,
      "geneset_oi" = geneset_oi,
      "background_expressed_genes" = background_expressed_genes,
      "expressed_ligands" = expressed_ligands,
      "expressed_receptors" = expressed_receptors,
      "ligand_activities" = ligand_activities,
      "prioritization_weights" = prioritizing_weights,
      "prior_table" = prior_table)
    return(nn_output)
  }
  ###############################################################################|
  # 7. PLOT------------
  ###############################################################################|
  if (!run_nichenet){
    prior_table = plot.prior.table
  }
  # plot only the top_n interactions ( ligand x sender x receptor)
  prior_table_top = prior_table%>%slice(1:top_n)
  prior_table_top$ligand= factor(prior_table_top$ligand , levels = rev(prior_table_top$ligand%>%unique))
  prior_table_top$receptor = factor(prior_table_top$receptor , levels = rev(prior_table_top$receptor%>%unique))
  
  
  
  dir.create(figure_folder, recursive = T)
  
  

  ## Ligand activity-----
  p1 = plot_ligand_activity(prior_table_top,var_oi = "activity_zscore")
  leg1 = get_legend(p1)
  p1 = p1+ guides(fill=FALSE)
  p1
  ## Ligand dotplot ----
  p2 = dotplot_expression(prior_table,ligands_oi =  prior_table_top$ligand%>%unique%>%as.character, fill.var = "lfc")
  leg2 = get_legend(p2)
  p2 = p2+ guides(fill=FALSE, size = "none")
  p2
  ## Receptor heatmap ----
  p3 = plot_ligand_receptor(prior_table_top,var_oi = "lfc_receptor")
  leg3 = get_legend(p3)
  p3 = p3 + guides(fill = F)
  p3
  
  # Assemble plots -------
  # Assemble heatmap
  top_fig = plot_grid(p1,p2,p3, nrow = 1 , rel_widths = c(0.25, 1,1))
  bottom_fig = plot_grid(leg1, leg2, leg3,  nrow = 1)
  fig = plot_grid(top_fig, bottom_fig, nrow = 2, rel_heights = c(5, 1))
  
  pdf(file.path(figure_folder,sprintf("%s_heatmap_topn = %d _output.pdf",receiver,top_n)), width = 12.5, height =10 )
  print(fig)
  dev.off()


}

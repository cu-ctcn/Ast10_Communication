#' Calculate signature z-scores
#' 
#' 
#' 
require(logger)
require(Seurat)
require(dplyr)
require(scCustomize)
require(ggrastr)
require(ggplot2)

seurat_get_summary_scores = function(data, feature.list, 
                                     module.key = "MODULE",
                                     method = "mean.zscore",
                                     assay = "RNA",
                                     slot = "data",
                                     out.fig.dir = "figures/modules_score",
                                     plot = TRUE,
                                     plot.col.limits = NULL,
                                     plot.reduction = "umap",
                                     nbin = 24,
                                     seed = 1
                                     ){
  log_info("-----------------------------------------------------------------")
  log_info("               GET GENE MODULE SUMMARY SCORES                    ")
  log_info("-----------------------------------------------------------------")
  set.seed(seed)
  start = Sys.time()
  
  
  
  # Check whether specified features are expressed, remove them otherwise:
  
  DefaultAssay(data) = assay
  all.genes = rownames(data)
  
  log_info("Summary score for the following gene list(s):")
  gene.list.names = names(feature.list)
  log_info(" - ",paste0(gene.list.names,collapse = ", "))
  
  if ( method == "mean.zscore" ){
    log_info("")
    log_info("Calculating summary score using ",method," method....")
    
    
    avgerage_zscore = function(gs){
      
      # Detect whether genes in gene sets are found
      log_info(" - Module ",gs,":")
      genes_found = intersect(feature.list[[gs]],all.genes)
      genes_not_found = setdiff(feature.list[[gs]], genes_found )
      
      if(length(genes_not_found)>0){
        log_info("   * The following genes are not found:")
        log_info("     ",paste0(genes_not_found,collapse = ", "))}
      
      log_info("   * Number of features found: ", length(genes_found)," (",round(length(genes_found)*100/length(feature.list[[gs]])),"%)")
      
      # Calculate Summary Score
      if (length(genes_found)>1){
        gene_oi_exp= GetAssayData(data, slot = slot)[genes_found,]
        gene_oi_zexp = (gene_oi_exp-rowMeans(gene_oi_exp))/apply(gene_oi_exp, 1, sd)
        na_genes = rownames(gene_oi_zexp)[which(rowSums(gene_oi_zexp)%in%NaN)]
        
        if (length(na_genes)>0){
          log_info("Following genes contain NaN, removed them in calculation: ", paste0(na_genes,collapse = ", ")) 
          gene_oi_zexp = gene_oi_zexp[!rownames(gene_oi_zexp)%in%na_genes,]%>%as.matrix()
        }
        #gene_oi_zexp = gene_oi_zexp%>%as.matrix()%>%na.omit()
        summary_score = colMeans(gene_oi_zexp)%>%t()
        rownames(summary_score) = gs
      }else{
        log_warn("   * Less than 2 genes are found in ",gs,". Skipping calculation when using mean.zscore method")
        summary_score = NA
      }

      return(summary_score)
    }
    # Calculate summary scores  
    summary_scores = do.call("rbind",lapply(gene.list.names,avgerage_zscore))%>%na.omit()
    modules = CreateAssayObject(data = summary_scores)
    
    data[[paste0(module.key,".",method)]]  = modules
    
  }
  
  if (method == "seurat"){
    log_info("")
    log_info("Calculating summary score using ",method," method....")
    
    gene.list.names = names(feature.list)
    
    seurat_score_module = function(gs){
      # Detect whether genes in gene sets are found
      log_info(" - Module ",gs,":")
      genes_found = intersect(feature.list[[gs]],all.genes)
      genes_not_found = setdiff(feature.list[[gs]], genes_found )
      
      if(length(genes_not_found)>0){
        log_info("   * The following genes are not found:")
        log_info("     ",paste0(genes_not_found,collapse = ", "))}
      
      log_info("   * Number of features found: ", length(genes_found)," (",round(length(genes_found)*100/length(feature.list[[gs]])),"%)")
      
      # Calculate Summary Score
      data = AddModuleScore(data,features = list(feature.list[[gs]]), assay = assay, name = gs,nbin = nbin, slot = slot, seed = seed)
      
      gs_col = colnames(data@meta.data)[stringr::str_detect(colnames(data@meta.data),gs)]
      summary_score = data@meta.data%>%dplyr::select(all_of(gs_col))%>%as.matrix()%>%t()
      rownames(summary_score) = gs
      return(summary_score)
    }
    
    summary_scores = do.call("rbind",lapply(gene.list.names,seurat_score_module))
    
    modules = CreateAssayObject(data = summary_scores)
    
    data[[paste0(module.key,".",method)]]  = modules
  }
  end = Sys.time()
  log_info("Elapsed Time (Add Module Scores): ", round(end-start, 2) ," ", units(end-start))
  
  if (plot){
    if(!dir.exists(out.fig.dir)){dir.create(out.fig.dir,recursive = TRUE)}
    curr_assay = DefaultAssay(data)
    DefaultAssay(data) = paste0(module.key,".",method)
    
    #plot umap
    log_info("Plot module scores on ",plot.reduction," ....")
    
    
    
    plot_score_dim_reduc = function(gs){
      if (is.null(plot.col.limits)){
        plot.col.limits = quantile(GetAssayData(data, slot = slot)[gs,],c(0.15,0.85))
      }
      
        FeaturePlot_scCustom(data,slot = "data",reduction = plot.reduction, features = gs,na_cutoff=NA)+
        scale_color_viridis_c(limits = plot.col.limits,oob = scales::squish,direction = 1)

      #FeaturePlot_scCustom(data,slot = "data",reduction = plot.reduction, features = gs,na_cutoff=NA)
    }
    
    p_list = lapply(rownames(data[[paste0(module.key,".",method)]]),plot_score_dim_reduc)
    fig =  ggpubr::ggarrange(plotlist = p_list, nrow = 1,ncol = length(p_list),common.legend =FALSE)
    ggsave(plot = fig, filename = file.path(out.fig.dir,sprintf("%s_%s.png",plot.reduction,paste0(module.key,".",method))),width = length(p_list)*5.5, height = 5, limitsize = FALSE)
    
    
    #plot violin
    log_info("Plot module scores as violin plot....")
    p_vln = Stacked_VlnPlot(seurat_object = data, features = rownames(data[[paste0(module.key,".",method)]]),
                    x_lab_rotate = TRUE)
    p_vln = ggpubr::annotate_figure(p_vln,fig.lab =paste0(module.key,".",method), fig.lab.pos = "top.left",fig.lab.face = "bold",fig.lab.size = 14)
    ggsave(plot = p_vln, filename = file.path(out.fig.dir,sprintf("%s_%s.png","vln",paste0(module.key,".",method))),width = 7, height = 9)
    
    DefaultAssay(data) =curr_assay
  }
  
  return(data)

 }
 

###TEST####
# state_oi = "Ast.10"
# state.deg = read.csv(sprintf("/mnt/mfs/ctcn/team/natacha/ROSMAP_snRNAseq/Differential_Expression/20230124/de_analysis_output/%s_de_genes_pairwise_adjp=0.05_logfc=1.00.csv",state_oi))
# state.markers = list(filter(state.deg,( (ident.2==state_oi&logFC<0)| (ident.1==state_oi&logFC>0)) )%>%pull(gene)%>%unique())
# names(state.markers) =  sprintf("%s.module",state_oi)
# state.markers[["gene.set.2"]] = state.markers[[1]]
# 
# data_test = seurat_get_summary_scores(data, feature.list = state.markers, 
#                                      module.key = "MODULE",
#                                      method = "mean.zscore",
#                                      assay = "SCT",
#                                      out.fig.dir = "figures/modules_score",
#                                      plot = TRUE,
#                                      plot.col.limits = c(-0.5,0.5),
#                                      plot.reduction = "umap")

#----------------------------------------------------------#|
#      Pipeline for label transfer with seurat
#----------------------------------------------------------#|
#                Natacha Comandante-Lou 
#                 (nc3018@columbia.edu)
#----------------------------------------------------------#|
require(Seurat)
require(SeuratDisk)
require(SeuratWrappers)
require(tidyverse)
require(scCustomize)
require(ggpubr)
require(logger)
require(dplyr)

PLXNB1_ko_hipsc_label_transfer_pipeline = function(query.sub.list,query.assay = "Spatial",
                                   reference.sub,
                                   n_features = 4000, # for SCTransform
                                   normalization.method ="SCT",
                                   k = 10,
                                   npcs = 30,
                                   state.levels = NULL,
                                   output.folder = "ref_mapping_res",
                                   harmony_query = T,
                                   label_oi = "state",
                                   save.query = FALSE,
                                   save.ref = FALSE,
                                   plot.markers,
                                   pal
){
  
  output.dir = file.path(output.folder,sprintf("n_features=%d/npcs=%d_k=%d", n_features, npcs,k))
  if (!dir.exists(output.dir)){dir.create(output.dir, recursive = TRUE)}
  #############################################################################|
  # Log -------
  #############################################################################|
  log_layout(layout_simple)
  log_file = file.path(output.dir , "PLXNB1_ko_annotations_pipeline.log")
  if(file.exists(log_file)){
    log_warn("Overwriting previous log: ", log_file)
    file.remove(log_file)
  }
  
  log_appender(appender_tee(file = log_file))
  
  log_info("-----------------------------------------------------------------")
  log_info("      PLXNB1 KO DATA LABEL TRANSFER               ")
  log_info("-----------------------------------------------------------------")
  
  #############################################################################|
  # Transfer data -------
  #############################################################################|
  annotations = list()
  ref.umap_list = list()
  pca_list = list()
  
  log_info("Parameters:")
  log_info(" - n_features = ",n_features)
  log_info(" - nPCs = ",npcs)
  log_info(" - k = ",k)
  log_info("   ")

for (s in c(1:length(query.sub.list))){
  set.seed(1245)
  log_info("# Run ",s,":")
  start = Sys.time()
  query.sub = query.sub.list[[s]]
  #query.sub <- SCTransform(query.sub,variable.features.n = n_features,assay = query.assay)
  query.sub <- SCTransform(query.sub,vars.to.regress = "percent.mt",variable.features.n = n_features,assay = query.assay)
  log_info("* Transfering data ...")
  
  # anchors <- FindTransferAnchors(reference = reference.sub , k.anchor = k,
  #                                query = query.sub ,normalization.method = normalization.method,
  #                                dims = 1:npcs, reference.reduction = "pca")
  anchors <- FindTransferAnchors(reference = reference.sub , 
                                 query = query.sub ,normalization.method = normalization.method,
                                 dims = 1:npcs, reference.reduction = "pca")
  if (is.null(state.levels)){state.levels = levels(reference.sub$state)}
  
  #############################################################################|
  # Map Query -------
  #############################################################################|
  
  query.sub <- TransferData(anchorset = anchors, reference = reference.sub, query = query.sub,
                            refdata = list(state = label_oi),dims = 1:npcs,
                            k.weight = k
  )
  
  log_info("* Mappig query data ...")
  query.sub <- IntegrateEmbeddings(anchorset = anchors, reference = reference.sub,
                                   query = query.sub, new.reduction.name = "ref.pca")
  
  query.sub <- RunPCA(query.sub, npcs = npcs)%>%
               FindNeighbors(k.param = k,
                  features = VariableFeatures(object = query.sub)%>%intersect(VariableFeatures(object = reference.sub))) # use hvgs from reference that exists in query
  
  
  query.sub <- ProjectUMAP(query = query.sub, query.reduction = "ref.pca", reference = reference.sub,
                           k.param = k,
                           reference.reduction = "pca", reduction.model = "umap")
  
  query.sub@meta.data$predicted.state = factor(query.sub@meta.data$predicted.state, levels =state.levels )
  
 
  #############################################################################|
  # Output annotation table -------
  #############################################################################|
  log_info("* Outputting Annotations and Reference UMAP coordinates...")
  annotations[[s]] = query.sub@meta.data%>%select(contains("predict"))
  
  
  end = Sys.time()
  log_success("Elapsed Time"," (Run ",s,"): ", round(end-start, 2) ," ", units(end-start))


  
  query.sub.list[[s]] = query.sub
  ref.umap_list[[s]] = query.sub@reductions[["ref.umap"]]@cell.embeddings
  pca_list[[s]] = query.sub@reductions[["pca"]]@cell.embeddings
  
  rownames(ref.umap_list[[s]]) = paste(rownames(ref.umap_list[[s]]),s,sep="_")
  rownames(pca_list[[s]]) = paste(rownames(pca_list[[s]]),s,sep="_")
  gc()
  }
  # Merge query --------
  if (length(query.sub.list)>1){
    query = merge(query.sub.list[[1]],query.sub.list[-1],
                  project = "hiPSC PLXNB1 ko")
    log_info("Merged query...")
  }else{
    query = query.sub.list[[1]]
    
  }
  if (harmony_query & length(query.sub.list)>1){
    library(harmony)
    log_info("Runing Harmony....")
    query_harmony = query
  
    query_harmony <- SCTransform(query_harmony)
    query_harmony <- RunPCA(query_harmony,assay = "SCT")
    query_harmony <- RunHarmony(query_harmony, group.by.var = "condition",assay.use = "SCT",plot_convergence = TRUE)
    query_harmony <- RunUMAP(query_harmony, reduction = "harmony", dims = 1:30)
    
    p1 = DimPlot_scCustom(query_harmony, reduction = "umap", group.by = "condition", label = F, label.size = 3.5,
                     repel = TRUE) + ggtitle("Harmony Integration")
    p2 = DimPlot_scCustom(query_harmony, reduction = "umap", group.by = "predicted.state", label = T, label.size = 3,
                          repel = TRUE,colors_use = pal) + ggtitle("Harmony Integration")
    ggsave(p1+p2, filename = file.path(output.dir,sprintf("harmony_umap_condition_%s.pdf",DefaultAssay(query))), width = 13, height = 4.73)
    
    query = query_harmony
  }
  ref.umap.embeddings = do.call("rbind",ref.umap_list)
  pca.embeddings = do.call("rbind",pca_list)
  dim(ref.umap.embeddings)
  log_info("Merged embeddings")
  log_info(DefaultAssay(query))
  query[["ref.umap"]] <- CreateDimReducObject(embeddings = ref.umap.embeddings, key = "refUMAP_", assay = DefaultAssay( query))
  query[["pca"]] <- CreateDimReducObject(embeddings =pca.embeddings, key = "PCA_", assay = DefaultAssay( query))
  
  gc()

  #############################################################################|
  # Plot -------
  #############################################################################|

  log_info("* Generating plots....")
  # Compare Marker Distr -------
  p1 = Stacked_VlnPlot(seurat_object = reference.sub, features = plot.markers,group.by = label_oi,colors_use = pal,x_lab_rotate = T)+ ggtitle("Reference")

  p2 = Stacked_VlnPlot(seurat_object = query, features = plot.markers,colors_use = pal,group.by = "predicted.state",x_lab_rotate = T)+ ggtitle("Query (Seurat Label Transfer)")

  ggsave(ggarrange(p1, p2, nrow = 1),filename = file.path(output.dir, "marker_vln.pdf"),width = 12, height = 20)


  # Compare UMAP-----
  p1 = DimPlot_scCustom(reference.sub, reduction = "umap", group.by = label_oi, label = TRUE, label.size = 3.5,
                        repel = TRUE,colors_use = pal) + ggtitle("Reference annotations")+coord_cartesian(xlim = c(-15,15), ylim =  c(-15,15))

  p2 = DimPlot_scCustom(query, reduction = "ref.umap", group.by = "predicted.state", label = FALSE,pt.size = 0.8,
                        label.size = 4, repel = TRUE,colors_use = pal) + ggtitle("Query transferred labels") +coord_cartesian(xlim = c(-15,15), ylim =  c(-15,15))

  ggsave(p1+p2,filename = file.path(output.dir, "ref.umaps.pdf"),width = 15, height = 5)

  # Plot Markers expression on UMAP -----
  p = FeaturePlot_scCustom(query,features =plot.markers,reduction = "ref.umap",num_columns = 3)
  ggsave(p, filename = file.path(output.dir,sprintf("markers_expr_ref.umap_%s.pdf",DefaultAssay(query))), width = 12, height = 30)
  
  
  # Plot PLXNB1 expression on UMAP -----
  p = FeaturePlot_scCustom(query,features = "PLXNB1",reduction = "ref.umap",num_columns = s,split.by = "condition")
  ggsave(p, filename = file.path(output.dir,sprintf("PLXNB1_expr_ref.umap_%s.pdf",DefaultAssay(query))), width = s*5, height = 5)
  if (harmony_query & length(query.sub.list)>1){
  p = FeaturePlot_scCustom(query,features = "PLXNB1",reduction = "umap",num_columns = s,split.by = "condition")
  ggsave(p, filename = file.path(output.dir,sprintf("PLXNB1_expr_umap_%s.pdf",DefaultAssay(query))), width = s*5, height = 5)
  }
  
  # Plot PLXNB1 expression on Vln -----
  p = Stacked_VlnPlot(seurat_object = query, features = c("PLXNB1","GFAP"),colors_use = pal,group.by = "condition",x_lab_rotate = T)
  ggsave(p,filename = file.path(output.dir, "PLXNB1_vln.pdf"),width = 5.6, height = 3.5)
#############################################################################|
# Save-------
#############################################################################|
#Save annotations and reduction embeddings
saveRDS(annotations,file = file.path(output.dir,"annotations.rds"))
reductions = cbind(ref.umap_list,pca_list)
saveRDS(reductions,file = file.path(output.dir,"reductions.rds"))

# Save query seurat obj
if (save.query){
  log_info("*Saving query seurat object...")
  
  saveRDS(query,file = file.path(output.dir,"annotated_query.rds"))
  
}
# Save reference seurat obj
if (save.ref){
  log_info("*Saving refer seurat object...")
  saveRDS(reference.sub,file = file.path(output.dir,"reference.rds"))
}

return(list(query = query,
            annotations = annotations,
            reductions = reductions
))
  
}
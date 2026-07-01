#------------------------------------------------#|
#         Spatial correlation analysis 
#           Supporting functions
#------------------------------------------------#|
#           Natacha Comandante-Lou 
#             (nc3018@columbia.edu)
#------------------------------------------------#|
require(scCustomize)
require(Seurat)
lay_pal = c("#F0027F", "#377EB8", "#4DAF4A", "#984EA3", "#FFD700", "#FF7F00", "#1A1A1A", "#666666")
my_pal = scCustomize_Palette(14,ggplot_default_colors = F)[-c(1:2)]
###############################################################################|
# Get Neighbor --------------
###############################################################################|
gatherGenes = function(genes_oi, .st, obj_type = "Seurat", assay_oi = "SCT", method = "mean.zscores",geneset_name){
  # SpatialExperient as Input Object--------------------
  if (obj_type == "SpatialExperiment"){
    # Merge spot info with selected gene expression
    
    coord = data.frame(x = .st@int_colData@listData[["spatialCoords"]][["pxl_col_in_fullres"]],
                       y = .st@int_colData@listData[["spatialCoords"]][["pxl_row_in_fullres"]])
    coord$SpotID = rownames(coord)
    
    # Get expressions of genes of interests
    if (length(genes_oi) > 1){
      exp = data.frame(t(as.matrix(assay(.st,assay_oi)[genes_oi,])))
    }else{
      exp = data.frame(as.matrix(assay(.st,assay_oi)[genes_oi,]))
      colnames(exp) = genes_oi
    }
    
    
    
    umap = data.frame(.st@int_colData@listData[["reducedDims"]]@listData[["UMAP"]])
    umap$SpotID = rownames(umap)
    
    
    if (method == "none"){
      exp$SpotID = rownames(exp)
      
      output = merge(coord,exp, by = "SpotID")
    }
    
    if (method == "zscores"){
      exp.z = as.data.frame(sapply(exp, function(exp) (exp-mean(exp))/sd(exp)))
      
      rm_col = colSums(is.na(exp.z))>1
      exp.z = exp.z[,!rm_col]
      
      if (sum(rm_col) > 1){warning(sprintf("removing NaN columns: %s", paste0(colnames(exp.z)[rm_col],collapse = ", ")))}
      
      rownames(exp.z) = rownames(exp)
      exp.z$SpotID = rownames(exp.z)
      output = merge(coord,exp.z, by = "SpotID")
    }
    
    if (method == "mean.zscores"){
      exp.z = as.data.frame(sapply(exp, function(exp) (exp-mean(exp))/sd(exp)))
      
      rm_col = colSums(is.na(exp.z))>1
      exp.z = exp.z[,!rm_col]
      
      if (sum(rm_col) > 1){warning(sprintf("removing NaN columns: %s", paste0(colnames(exp.z)[rm_col],collapse = ", ")))}
      
      rownames(exp.z) = rownames(exp)
      
      exp.z.avg = data.frame(SpotID = rownames(exp.z), mean.zscores = rowMeans(exp.z))
      colnames(exp.z.avg)[colnames(exp.z.avg)%in%"mean.zscores"] = paste0(geneset_name,"-mean.zscores")
      
      output = merge(coord,exp.z.avg, by = "SpotID")
    }
  }
  
  # Seurat as Input Object--------------------
  if(obj_type == "Seurat"){
    # Merge spot info with selected gene expression
    #coord = data.frame(x = .st@images$slice1@coordinates$col, y = .st@images$slice1@coordinates$row)
    coord = data.frame(x = .st@images[["slice1"]]@coordinates[["imagecol"]], y = .st@images[["slice1"]]@coordinates[["imagerow"]])
    coord$SpotID = rownames(.st@meta.data)
    
    # Get expressions of genes of interests
    
    exp = data.frame(t(as.matrix(.st@assays[[assay_oi]]@data[genes_oi,])))
    colnames(exp) = genes_oi
    
    umap = data.frame(.st@reductions[["umap"]]@cell.embeddings)
    umap$SpotID = rownames(umap)
    
    
    if (method == "none"){
      exp$SpotID = rownames(exp)
      
      output = merge(coord,exp, by = "SpotID")
    }
    
    if (method == "zscores"){
      exp.z = as.data.frame(sapply(exp, function(exp) (exp-mean(exp))/sd(exp)))
      
      rm_col = colSums(is.na(exp.z))>1
      exp.z = exp.z[,!rm_col]
      
      if (sum(rm_col) > 1){warning(sprintf("removing NaN columns: %s", paste0(colnames(exp.z)[rm_col],collapse = ", ")))}
      
      rownames(exp.z) = rownames(exp)
      exp.z$SpotID = rownames(exp.z)
      output = merge(coord,exp.z, by = "SpotID")
    }
    
    if (method == "mean.zscores"){
      exp.z = as.data.frame(sapply(exp, function(exp) (exp-mean(exp))/sd(exp)))
      
      rm_col = colSums(is.na(exp.z))>1
      exp.z = exp.z[,!rm_col]
      
      if (sum(rm_col) > 1){warning(sprintf("removing NaN columns: %s", paste0(colnames(exp.z)[rm_col],collapse = ", ")))}
      
      rownames(exp.z) = rownames(exp)
      
      exp.z.avg = data.frame(SpotID = rownames(exp.z), mean.zscores = rowMeans(exp.z))
      colnames(exp.z.avg)[colnames(exp.z.avg)%in%"mean.zscores"] = paste0(geneset_name,"-mean.zscores")
      
      output = merge(coord,exp.z.avg, by = "SpotID")
    }
    
    if(method == "seurat"){
      
      
    }
  }  
  
  output = merge(output,umap,by="SpotID")
  return(output)
  
} 

source('/mnt/mfs/ctcn/team/natacha/RScript/seurat_get_summary_scores.R')
prepareSummaryScores = function(.st,feature.list,assay,method, slot = "data"){
  #slot = slot to fetch data from
  .st = seurat_get_summary_scores(.st, feature.list, 
                                  module.key = paste0("MODULE_",assay),
                                  method = method,
                                  assay = assay,
                                  slot = slot,
                                  out.fig.dir = NULL,
                                  plot = F,
                                  plot.col.limits = NULL,
                                  plot.reduction = "umap"
  )
  var_oi = rownames(.st@assays[[paste0("MODULE_",assay,".",method)]]@data)
  print(var_oi)
  DefaultAssay(.st) = paste0("MODULE_",assay,".",method)
  output = FetchData(.st, assay = paste0("MODULE_",assay), vars = var_oi)%>%
    rownames_to_column("SpotID")
  
  
  coord = data.frame(x = .st@images[["slice1"]]@coordinates[["imagecol"]], y = .st@images[["slice1"]]@coordinates[["imagerow"]])
  coord$SpotID = rownames(.st@meta.data)
  
  output = merge(output,coord, by = "SpotID")
  
  return(output)
}


getSpotsInDisk = function(.coord,xc,yc,r,pxl_scale=156){
  # r: radius (in pxl)
  # (x-xc)^2+(y-yc)^2 <= r^2
  #pxl_scale = 156
  x = .coord[["x"]]
  y = .coord[["y"]]
  
  d = (x-xc)^2+(y-yc)^2
  
  if (d <=(r*pxl_scale)^2){in_disk = TRUE}else{in_disk=FALSE}
  if(x==xc&y==yc){in_disk="CENTER"}
  return(in_disk)
}
getSpotsInRing = function(.coord,xc,yc,r,pxl_scale=156){
  # r: radius (in pxl)
  # (x-xc)^2+(y-yc)^2 = r^2
  #pxl_scale = 156
  x = .coord[["x"]]
  y = .coord[["y"]]
  
  d = (x-xc)^2+(y-yc)^2
  
  if (d > ((r-0.5)*pxl_scale)^2 & d < ((r+0.5)*pxl_scale)^2){in_ring = TRUE}else{in_ring=FALSE}
  if(x==xc&y==yc){in_ring="CENTER"}
  return(in_ring)
}


sumNeighhors = function(id,.df, r, feature,shape = "disk",pxl_scale,include_center = TRUE){
  coord = dplyr::select(.df,c("x","y"))
  coord_c = dplyr::filter(.df, SpotID%in%id)%>%dplyr::select(c("x","y"))
  xc = coord_c[["x"]]; yc = coord_c[["y"]]
  if (shape=="disk"){
    .df$in_disk = apply(coord,MARGIN = 1, FUN = getSpotsInDisk,xc,yc, r,pxl_scale)
    if (include_center){
      df.disk = filter(.df, !in_disk%in%FALSE)}
    else{
      df.disk = filter(.df, in_disk%in%TRUE)
    }
    output = sum(df.disk[[feature]])
  }
  if (shape=="ring"){
    .df$in_ring = apply(coord,MARGIN = 1, FUN = getSpotsInRing,xc,yc, r,pxl_scale)
    if (include_center){
      df.ring = filter(.df, !in_ring%in%FALSE)
    }else{
      df.ring = filter(.df, in_ring%in%TRUE)
    }
    output = sum(df.ring[[feature]])
  }
  if (shape=="random"){
    #shuffle
    .df.shuffle = .df
    .df.shuffle$in_disk = .df$in_disk[sample(c(1:nspots),nspots)]
    .df.shuffle$in_disk[.df.shuffle$in_disk=="CENTER"]=FALSE 
    .df.shuffle$in_disk[.df$in_disk=="CENTER"]="CENTER" #replace center as in coord
    if (include_center){
      df.disk = filter(.df.shuffle, !in_disk%in%FALSE)}
    else{
      df.disk = filter(.df.shuffle, in_disk%in%TRUE)
    }
    output = sum(df.disk[[feature]])
  }
  return(output)
} 

averageNeighhors = function(id,.df, r, feature,shape = "disk",pxl_scale, include_center = TRUE){
  coord = dplyr::select(.df,c("x","y"))
  coord_c = dplyr::filter(.df, SpotID%in%id)%>%dplyr::select(c("x","y"))
  xc = coord_c[["x"]]; yc = coord_c[["y"]]
  if (shape=="disk"){
    .df$in_disk = apply(coord,MARGIN = 1, FUN = getSpotsInDisk,xc,yc, r,pxl_scale)
    
    if (include_center){
      df.disk = filter(.df, !in_disk%in%FALSE)}
    else{
      df.disk = filter(.df, in_disk%in%TRUE)
    }
    output = mean(df.disk[[feature]])
  }
  if (shape=="ring"){
    .df$in_ring = apply(coord,MARGIN = 1, FUN = getSpotsInRing,xc,yc, r,pxl_scale)
    
    if (include_center){
      df.ring = filter(.df, !in_ring%in%FALSE)
    }else{
      df.ring = filter(.df, in_ring%in%TRUE)
    }
    
    output = mean(df.ring[[feature]])
  }
  if (shape=="random"){
    .df$in_disk = apply(coord,MARGIN = 1, FUN = getSpotsInDisk,xc,yc, r,pxl_scale)
    nspots = length(.df$in_disk)
    
    
    #shuffle
    .df.shuffle = .df
    .df.shuffle$in_disk = .df$in_disk[sample(c(1:nspots),nspots)]
    .df.shuffle$in_disk[.df.shuffle$in_disk=="CENTER"]=FALSE 
    .df.shuffle$in_disk[.df$in_disk=="CENTER"]="CENTER" #replace center as in coord
    
    if (include_center){
      df.disk = filter(.df.shuffle, !in_disk%in%FALSE)}
    else{
      df.disk = filter(.df.shuffle, in_disk%in%TRUE)
    }
    
    output = mean(df.disk[[feature]])
  }
  return(output)
} 

gatherNeighborGenes = function(geneset.list, .st, obj_type = "Seurat",assay_oi = "SCT", slot = "data",method = "mean.zscores", neighbor_method = "mean", r, pxl_scale, shape = "disk",include_center = TRUE){
  geneset_name = names(geneset.list)
  #.df = gatherGenes(genes_oi = genes_oi ,obj_type = obj_type,.st = .st ,assay_oi = assay_oi, method = method ,geneset_name = geneset_name)
  
  .df = prepareSummaryScores(.st,geneset.list,assay_oi,method, slot)
  # for naming new columns
  
  feature_list = colnames(.df)%>%setdiff(c("SpotID","x","y"))
  
  
  if (neighbor_method %in% "sum"){
    
    feature.neighbor = sprintf("%s_sum-%s=%0.1f",feature,shape,r)
    .df[[feature.neighbor]] = apply(.df%>%dplyr::select("SpotID"), MARGIN=1, FUN = sumNeighhors, .df, r, feature, shape,pxl_scale,include_center)
  }
  
  if (neighbor_method %in% "mean"){
    for(feature in feature_list){
      feature.neighbor = sprintf("%s_mean-%s=%0.1f",feature,shape,r)
      .df[[feature.neighbor]] = apply(.df%>%dplyr::select("SpotID"), MARGIN=1, FUN = averageNeighhors, .df, r, feature,shape,pxl_scale,include_center)
    }
  }
  # test = lapply(.df$SpotID, averageNeighhors,.df,r,feature)
  
  return(.df)
}


example_neighbor_plot = function(coord, xc, yc,r, save.as = "example_neighborhood.pdf"){
  ## Plot neighbors of a given dot as an example --------

  coord$in_disk = apply(coord[,c(1:2)],MARGIN = 1, FUN = getSpotsInDisk,xc,yc, r,pxl_scale=156)
  coord$in_disk = factor(coord$in_disk, levels= c("CENTER","TRUE","FALSE"))
  coord = arrange(coord, -as.numeric(in_disk))
  nspots = length(coord$in_disk)
  
  #shuffle
  coord.shuffle = coord
  coord.shuffle$in_disk = coord$in_disk[sample(c(1:nspots),nspots)]
  coord.shuffle$in_disk[coord.shuffle$in_disk=="CENTER"]=FALSE 
  coord.shuffle$in_disk[coord$in_disk=="CENTER"]="CENTER" #replace center as in coord
  coord.shuffle$in_disk = factor(coord.shuffle$in_disk, levels= c("CENTER","TRUE","FALSE"))
  coord.shuffle = arrange(coord.shuffle, -as.numeric(in_disk))
  p1 = ggplot(coord, aes(x = x, y = y, color =  in_disk)) +
    geom_point(size = 1.2) +
    scale_color_manual(values = c(lay_pal[5],lay_pal[4],"#B5C0D0"))+
    coord_fixed() +
    ggtitle("Neighbors") +
    theme_void()
  p2 = ggplot(coord.shuffle, aes(x = x, y = y, color =  in_disk)) +
    geom_point(size = 1.2) +
    scale_color_manual(values = c(lay_pal[5],lay_pal[4],"#B5C0D0"))+
    coord_fixed() +
    ggtitle("Shuffled Neighbors") +
    theme_void()
  ggarrange(p1,p2,common.legend = TRUE)
  ggsave(save.as,device = "pdf",width = 7.69, height = 2.69)

}
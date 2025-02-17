require(ggrepel)

reticulate::use_condaenv('~/conda/my-envs/magic' )
my_pal = c("#8338ec","#ffbe0b","#3a86ff","#fb5607")
##############################################################################|
# MAGIC IMPUTATION ----------
###############################################################################|
source("../../utils/run_magic.R")


run_magic_pipeline = function(x,magic_output_dir="magic", t_list = c(3,5,10), knn=5,min_cells = 10,
                              gene_x_viz="PLXNB1", gene_y_viz="SMTN",
                              default_assay = "RNA", default_slot = "counts",plot = T, normalize = T, sqrt = T){
  logger::log_info("Running MAGIC on ", default_assay, " ...")
  DefaultAssay(x) = default_assay

  data = GetAssayData(x,slot = default_slot,assay = default_assay)%>%t()

  # Keep genes expressed in at least 10 cells
  keep_cols <- colSums(data > 0) > min_cells
  data <- data[,keep_cols]
  
  # Normalize and transform
  if (normalize){  data <- library.size.normalize(data)}
  if (sqrt){   data <- sqrt(data)}
  
  # run MAGIC
  all_genes = colnames(data)#%>%setdiff(mt_genes)
  
  output.dir = magic_output_dir
  dir.create(output.dir, recursive = T)
  final_res = lapply(t_list,run_magic,data,all_genes,
                     knn=knn,
                     plot=plot,
                     plot_only = FALSE,
                     init = NULL,
                     solver='approximate',
                     gene_x_viz=gene_x_viz,gene_y_viz=gene_y_viz,
                     save.as = "magic",
                     output.dir = output.dir )
  # Add Results
  for (i in 1:length(t_list)){
    t = t_list[i]
    magic_assay <- CreateAssayObject(data = final_res[[i]]$data_MAGIC[["result"]]%>%as.matrix()%>%t())
    # add this assay to the previously created Seurat object
    x[[sprintf("MAGIC_knn=%d_t=%d",knn,t)]] <- magic_assay
    Misc(x, slot = "magic") = final_res[[i]]$data_MAGIC[["result"]]$operator
  }
  
  # Plot
  assay_list = Assays(x)[Assays(x)%>%grep(pattern="MAGIC_")]
  lapply(assay_list, function(y){
    tmp_assay = DefaultAssay(x) 
    DefaultAssay(x) = y
    x = x%>%ScaleData()%>% RunPCA( .,assay = y,features = rownames(x))%>%RunUMAP(dims = 1:30,reduction.name = paste0(y,"_umap"))
    p = DimPlot_scCustom(x, reduction = paste0(y,"_umap"),label = F, group.by = "condition",colors_use = my_pal) 
    DefaultAssay(x) =  tmp_assay
    return(p)
  })%>%ggarrange(plotlist = ., nrow = 1)%>%
    #annotate_figure(p = ., fig.lab = unique(x$orig.ident), fig.lab.size = 14, fig.lab.face = "bold")%>%
    ggsave(file.path(output.dir, "umap_magic.pdf"),plot = ., width = 4*length(assay_list)+1, height = 4,create.dir = T)
  
  return(x)
  
  
}


##############################################################################|
# Calculate Module Score----------
###############################################################################|
source("../../../utils/seurat_get_summary_scores.R")

## Get list of ast10 markers -------
de_output_dir = "ast_mic_subtype_signature_genes_summary/"
adjp_thres.n = 0.05
fc_thres.n = 1.00
state_oi = "Ast.10"
de_file_name = sprintf("%s%s_de_genes_summary_adjp=%0.2f_logfc=%0.2f.csv",de_output_dir,state_oi, adjp_thres.n,fc_thres.n)
degs<- read.csv(de_file_name,check.names = FALSE)%>%
  arrange(-avg.logFC) # DEGs from my pseudobulk differential expression analysis


##############################################################################|
# Density Plot -------------
##############################################################################|

plot_density = function(s,.df = NULL, y_assay_oi=NULL, x_assay_oi=NULL,   
                        density_oi = "ndensity",  
                        yvar = "Ast.10.pos.n", xvar = "PLXNB1",
                        yinter = NULL,
                        xinter = NULL,
                        gate.method = "gaussian", # or median,
                        gate.ref, # reference condition for gating
                        out_file
                        ){

  if (!is.null(s)){
    # balance condition to fit gmm
    cellcount = group_by(s@meta.data, condition)%>%tally()
    Idents(s) = "condition"
    n_sample = min(cellcount$n)
    set.seed(66362)
    s.sampled <- subset(x = s, downsample = n_sample)
    
    ## Fit gaussian mixed model for gating
    if(is.null(yinter) & gate.method == "gaussian"){
    y_gate = get_gauss_mix_gate(s.sampled,assay_oi=y_assay_oi, var_oi = yvar, n = 2)%>%
      group_by( grp)%>%summarise(mean = mean(x))
    yinter = y_gate$mean%>%mean
    }
    
    if(is.null(xinter)& gate.method == "gaussian"){
    x_gate = get_gauss_mix_gate(s.sampled,assay_oi=x_assay_oi, var_oi =xvar, n = 2)%>%
      group_by( grp)%>%summarise(mean = mean(x))
    xinter = x_gate$mean%>%mean
    }
    

    # extract data to plot
    df.plot = s@assays[[y_assay_oi]]@data%>%t()%>%as.data.frame
    df.plot= s@assays[[x_assay_oi]]@data%>%t()%>%as.data.frame%>%select(xvar)%>%cbind(.,df.plot)
    df.plot$cell_id = rownames(df.plot)
    s@meta.data$cell_id = rownames(s@meta.data)
    df.plot = merge(df.plot, s@meta.data, by = "cell_id")
  
    ## Use median for gating
    if(is.null(yinter) & gate.method == "median"){
      if (!is.null(gate.ref)){
        # use cells from condition == gate.ref to calculate median
        s.ref = s%>%subset(subset=condition==gate.ref)
      }else{
        s.ref = s
      }
      
      df.ref = s.ref@assays[[y_assay_oi]]@data%>%t()%>%as.data.frame
      yinter = median(df.ref[[yvar]])
    }
    if(is.null(xinter) & gate.method == "median"){
      if (!is.null(gate.ref)){
        s.ref = s%>%subset(subset=condition==gate.ref)
      }else{
        s.ref = s
      }
      
      df.ref = s.ref@assays[[x_assay_oi]]@data%>%t()%>%as.data.frame
      xinter = median(df.ref[[xvar]])
    }
    
  }else{
    
    df.plot = .df
    ## Fit gaussian mixed model for gating
    if(is.null(yinter) & gate.method == "gaussian"){
      y_gate = get_gauss_mix_gate(s.sampled,assay_oi=y_assay_oi, var_oi = yvar, n = 2)%>%
        group_by( grp)%>%summarise(mean = mean(x))
      yinter = y_gate$mean%>%mean
    }
    
    if(is.null(xinter)& gate.method == "gaussian"){
      x_gate = get_gauss_mix_gate(s.sampled,assay_oi=x_assay_oi, var_oi =xvar, n = 2)%>%
        group_by( grp)%>%summarise(mean = mean(x))
      xinter = x_gate$mean%>%mean
    }
    
    
    ## Use median for gating
    if(is.null(yinter) & gate.method == "median"){
      yinter = median(df.plot[[yvar]])
    }
    if(is.null(xinter) & gate.method == "median"){
      xinter = median(df.plot[[xvar]])
    }
    
  }

  
  
  log_info("x-intercept: ",xinter)
  log_info("y-intercept: ",yinter)
  nudgey = abs(yinter)*0.6#.1
  nudgex = abs(xinter)*0.6#.1
  #-----------------------------------
  
  df.plot = df.plot%>%mutate(gated.y = case_when(.data[[yvar]]>yinter ~ "y+" , TRUE ~ "y-"))%>%
    mutate(gated.x = case_when(.data[[xvar]]>xinter ~ "x+" , TRUE ~ "x-"))%>%
    mutate(gated = paste0(gated.x, gated.y))
  # calculate percent of cells in each quadrant
  pct = group_by(df.plot, condition, gated)%>%tally()%>%mutate(pct = sprintf("%0.1f %%",n*100/sum(n)))%>%
    mutate(x = case_when( str_detect(gated, "x\\+") ~ xinter + nudgex, TRUE ~ xinter - nudgex))%>%
    mutate(y = case_when( str_detect(gated, "y\\+") ~ yinter+ nudgey, TRUE ~ yinter - nudgey))
  
  p2 = ggplot(df.plot, aes(x = .data[[xvar]], y = .data[[yvar]], group = condition))+
    geom_point(alpha = 0.5, size = 0.3, color = "#36454F")+
    stat_density_2d(geom = "polygon",contour_var =  density_oi,
                    aes(alpha = after_stat(level), fill = ..level..),
                    bins = 20) +
    geom_hline(yintercept =  yinter, linetype = "dashed", color = "gray", linewidth = 0.5)+
    geom_vline(xintercept =  xinter, linetype = "dashed", color = "gray", linewidth = 0.5)+
    geom_text_repel(data = pct, aes(x, y , label = pct), color = "blue")+
    scale_alpha_binned(range = c(0.1,0.7))+
    #coord_cartesian(ylim = c(-0.03, 0.06), xlim = c(0, 0.42))+
    facet_wrap(facets = vars(condition))+
    # scale_fill_continuous_sequential(palette = "ag_Sunset", rev = FALSE) +
    scale_fill_viridis_c(option = "turbo", oob=scales::squish)+
    #scale_fill_viridis_c(option = "turbo",limits = c(50,600), oob=scales::squish)+
    theme_classic()
  
  p2
  #-----------------------------------
  p3 = ggplot(df.plot, aes(x= .data[[xvar]], y = .data[[yvar]], color = condition))+
    geom_point(size = 0.2)+
    geom_smooth(method = "lm")+
    scale_color_manual(values = my_pal)+
    theme_classic()
  p3
  ggsave(plot = cowplot::plot_grid(p2,p3, rel_widths = c(0.5,0.3)), filename = paste0(out_file,"_", density_oi,".pdf"), width = 12, height = 3.5,create.dir = T)
  
  return(p2+p3)
}


###############################################################################|
# GATE AVERAGED DENSITY ------
###############################################################################|
reticulate::source_python("/mnt/mfs/ctcn/team/natacha/RScript/py_gauss_mix.py")
np = reticulate::import("numpy")
sk = reticulate::import("sklearn")
scprep = reticulate::import("scprep")

get_gauss_mix_gate = function(.sobj,.df = NULL, assay_oi=NULL, var_oi, n = 2){
  
  if(!is.null(.sobj)){
    DefaultAssay(.sobj) = assay_oi
    X = GetAssayData(.sobj,assay =assay_oi)%>%t()%>%as.data.frame()
  }else{ 
    X = .df
  }
  x_oi = np$array(select(X, var_oi)%>%as.matrix())
  # Fit mixture
  mixture_model = sk$mixture$GaussianMixture(n_components=as.integer(n))
  classes = mixture_model$fit_predict(x_oi)
  

  if(!is.null(.sobj)){
  # Reorder the class labels so that 0 has the lowest average likelihood and 2 has the highest.
  classes = tibble(grp = scprep$utils$sort_clusters_by_values(classes, x_oi),
                   x = pull(X, var_oi)
                   )%>%
    `rownames<-`(rownames(.sobj@meta.data))
  }else{
    # Reorder the class labels so that 0 has the lowest average likelihood and 2 has the highest.
    classes = tibble(grp = scprep$utils$sort_clusters_by_values(classes, x_oi),
                     x = pull(X, var_oi)
    )%>%
      `rownames<-`(.df$cell_id)
    
  }
  return(classes)
}



##############################################################################|
# Violin plot --------
##############################################################################|


require(colorspace)
require(tidyverse)
require(ggplot2)
require(rstatix)
require(ggbeeswarm)
jitter_boxplot_w_stats = function( class_oi,vars_oi,assay_oi = "MODULE.mean.zscore",.sobj,  alternative = "two.sided", 
                                   save_file="figures/boxplot.pdf",method = "t.test",
                                   my_comparisons, class_pal=NULL){
  
  DefaultAssay(.sobj)= assay_oi
  #if(!dir.exists(save_dir)){dir.create(save_dir, recursive = TRUE)}
  .df = FetchData(.sobj,assay = assay_oi, vars = vars_oi)
  .df$cell.id = rownames(.df)
  
  meta = .sobj@meta.data%>%select(all_of(class_oi))
  meta$cell.id = rownames(meta)
  
  .df = merge(.df,meta, by = "cell.id")
  
  if (length(my_comparisons)>1){p_label = "p.adj"; p_signif = "p.adj.signif"}else{p_label = "p"; p_signif = "p.signif"}
  
  if(is.null(class_pal)){  class_pal = c("#8338ec","#3a86ff","#fb5607","#ffbe0b")}
  
  plist = lapply(class_oi,
                 function(class_oi){
                   if (method == "t.test"){
                     logger::log_info("t.test")
                   stat.test <- .df%>% ungroup()%>%t_test(comparisons = my_comparisons,alternative = alternative,
                                                          formula = as.formula(sprintf("%s ~ %s",vars_oi, class_oi)))%>%
                     add_xy_position(x = class_oi, dodge = 0.9)%>%
                     add_significance(
                       p.col = p_label,
                       output.col = NULL,
                       cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),
                       symbols = c("****", "***", "**", "*", "ns")
                     )
                   logger::log_info("calculated stats")
                   }
                   if (method == "wilcox.test"){
                     logger::log_info("wilcox.test")
                     stat.test <- .df%>% ungroup()%>%wilcox_test(comparisons = my_comparisons,alternative = alternative,
                                                            formula = as.formula(sprintf("%s ~ %s",vars_oi, class_oi)))%>%
                       add_xy_position(x = class_oi, dodge = 0.9)%>%
                       add_significance(
                         p.col = p_label,
                         output.col = NULL,
                         cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),
                         symbols = c("****", "***", "**", "*", "ns")
                       )
                     logger::log_info("calculated stats")
                     }
                   
                   ggplot(data = .df, aes(x = .data[[class_oi]], y = .data[[vars_oi]],group = .data[[class_oi]])) +
                     geom_quasirandom(groupOnX=TRUE,size = 1,aes(color=.data[[class_oi]]))+
                     geom_boxplot(alpha = 0.3,outlier.shape = NA)+
                     scale_color_manual(values = class_pal)+
                     guides(color=guide_legend(class_oi,override.aes = list(size= 3)) )+
                     stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01)+
                     #stat_pvalue_manual(stat.test, label = p_signif, tip.length = 0.01,step.increase = 0.03)+
                     ggtitle(vars_oi)+
                     #ylim(c(0.5,2))+
                     #xlab(paste0(stringr::str_split(class_oi,"\\.",simplify = TRUE)[,1:2], collapse = "-"))+
                     theme(panel.background = element_blank(),
                           legend.text = element_text(size = 14),
                           legend.title = element_text(size = 14),
                           axis.line = element_line(colour = "black"),
                           axis.text.x = element_text(angle = 45, hjust = 1,size = 14),
                           axis.text.y = element_text(size = 14),
                           axis.title.x = element_text(size = 16),
                           axis.title.y = element_text(size = 16))})
  fig = ggarrange(plotlist = plist, nrow = 1, common.legend = TRUE)
  
  if (!is.null(save_file)){
    ggsave(fig, file = save_file,width = length(plist)*3+3,height = 6, dpi = 300)
  }
  return(fig)
}


##############################################################################|
# Correlation matrix --------
##############################################################################|
# see https://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization

plot_gene_corr_mat = function(sobj.sm, vars, assay_oi, save.as = "corr.mat.pdf",width = 6, height = 6){
  M = FetchData(sobj.sm, assays=assay_oi,vars =vars%>%intersect(rownames(sobj.sm)))
  cormat <-cor(M)
 
  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }
  # Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  # Reorder the correlation matrix
  cormat <- reorder_cormat(cormat)
  
  full_cormat = cormat
  upper_tri <- get_upper_tri(cormat)
  
  # Melt the correlation matrix
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  # Create a ggheatmap
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed()
  
  # Print the heatmap
  pdf(save.as, width = width, height = height)
  print(ggheatmap)
  dev.off()
  return( tibble(hm = list(ggheatmap), 
                 dd = list(as.dist((1-full_cormat)/2)),
                 hc = list(hclust(as.dist((1-full_cormat)/2))),
                 cormat = list(melted_cormat), vars = list(vars), assay = assay_oi))
}

plot_gene_corr_mat_v2 = function(sobj.sm, vars,  assay_oi, gene_clusters = NULL, k = 2, dist = "correlation", save.as = "corr.mat.pdf", main = "Correlation Matrix"){
  DefaultAssay(sobj.sm) = assay_oi
  M = FetchData(sobj.sm, assays=assay_oi,vars =vars%>%intersect(rownames(sobj.sm)))
  cormat <-cor(M)
  
  if(is.null(gene_clusters)){
    logger::log_info("Clustering...")
    # hierarchical clustering of gene-gene relationship
    res = pheatmap(cormat, cutree_rows = k, scale = "none",cluster_rows = T ,cluster_cols = T,
                   clustering_distance_rows = dist, clusteringer_distance_cols = dist,main = main,
                   treeheight_col=0,breaks = seq(-1,1, by = 0.02), 
                   color = colorRampPalette(c("blue", "white", "red"))(100),border_color = "black")
  
    #output cluster
    gene_clusters = tibble(gene = rownames(cormat))
    gene_clusters <- cbind(gene_clusters, 
                           cluster = cutree(res$tree_row, k = k))
    
    out = tibble(cormat = list(cormat), 
                 hc = list(res$tree_row),
                 order = list(res$tree_row$order),
                 gene_clusters = list(gene_clusters[res$tree_row$order,]),
                 hm = list(res)
    )
  }else{
    logger::log_info("Not clustering...")
    #recorder cormat based on gene_cluster
    
    cormat = cormat[gene_clusters$gene, gene_clusters$gene]

    gaps_row = which(diff(as.numeric(factor(gene_clusters$cluster)))!=0)
    res = pheatmap(cormat,cluster_rows = F ,cluster_cols = F,  scale = "none",
                   gaps_row = gaps_row, gaps_col = gaps_row, main = main,breaks = seq(-1,1, by = 0.02), 
                   color = colorRampPalette(c("blue", "white", "red"))(100),border_color = "black")
    
    
    out = tibble(cormat = list(cormat), 
                 hc = list(NULL),
                 order = list(gene_clusters$gene),
                 gene_clusters = list(gene_clusters),
                 hm = list(res)
    )
  }
  
  #plot
  pdf(save.as);print(res);dev.off()
  
  return(out)
}


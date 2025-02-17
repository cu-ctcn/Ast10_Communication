#----------------------------------------------------------#|
#             Preprocess hiPSC scRNAseq data 
#               with or without PLXNB1 KO
#----------------------------------------------------------#|
#                Natacha Comandante-Lou 
#                 (nc3018@columbia.edu)
#----------------------------------------------------------#|
rm(list = ls())
gc()

library(Seurat)
library(tidyverse)
library(scCustomize)
library(cowplot)
#Load data ---------
data.dir = "/mnt/mfs/ctcn/team/natacha/PLXNB1-KO-hiPSCs/iAstro_CRISPR_line_553/cellranger"
file_list = list.files(data.dir,recursive = F)%>%.[grepl("BIZ",.)]
out.dir = "../data/seurat_obj"
fig.dir = "../figures"
# dir.create(out.dir, recursive = T)
# dir.create(fig.dir, recursive = T)
###############################################################################|
# Load Data ---------------------
###############################################################################|
# QC:
# - removing cells with ≤800 or ≥7,000 genes
# - removing cells with ≤1800 or ≥50,000 UMI counts
# - keep genes express more than 2 cells
# - removing cells with >= 10 % mitochondrial reads

umi_thres = c(1800,50000)
feat_thres = c(800,7000) 
sobj_list = lapply(file_list, function(f){
  data <- Read10X_h5(file.path(data.dir,f,"filtered_feature_bc_matrix.h5"))
  sobj = CreateSeuratObject(counts = data, project = "iAstro_CRISPR_553")
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
  sobj$orig.ident = f
  
  # Visualize QC metrics as a violin plot
  f1 = VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = "orig.ident")
  f2a <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "percent.mt")+theme(legend.position = "none")
  f2b <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+theme(legend.position = "none")
  
  meta = sobj@meta.data


  f3 = ggplot(meta, aes(x=.data[[ "nCount_RNA"]])) +
    geom_histogram(binwidth=50) +
    geom_vline(xintercept = umi_thres, linetype = "dashed",color = "blue", linewidth = 0.5)+
    ggtitle(sprintf("%0.1f - %0.1f", umi_thres[1], umi_thres[2]))+
    xlab("UMI count")+
   # coord_cartesian(xlim = c(0,1550), ylim = c(0,20500))+
    theme_classic()
  
  f4 = ggplot(meta, aes(x=.data[[ "nFeature_RNA"]])) +
    geom_histogram(binwidth=50) +
    geom_vline(xintercept =feat_thres, linetype = "dashed",color = "blue", linewidth = 0.5)+
    xlab("Gene count")+
    ggtitle(sprintf("%0.1f - %0.1f", feat_thres[1], feat_thres[2]))+
    #coord_cartesian(xlim = c(0,10010),ylim = c(0,150))+
    theme_classic()
  
  
  fig = plot_grid(f1, plot_grid(f2a, f2b),plot_grid(f3, f4), nrow = 3)
  ggsave(fig, filename = file.path(fig.dir, paste0(f,"_qc.png")), width = 11, height = 15)
  
  return(sobj)
})

###############################################################################|

sobj_list.f = lapply(sobj_list, function(sobj){
  sobj <- subset(sobj, subset = nFeature_RNA > feat_thres[1] & nFeature_RNA < feat_thres[2] & percent.mt <= 10 & nCount_RNA < umi_thres[2]&  nCount_RNA > umi_thres[1])
  all_genes = sobj@assays[["RNA"]]@counts%>%rownames()
  keep_genes = all_genes[(sobj@assays[["RNA"]]@counts%>%rowSums()>2)]
  sobj  <- subset(sobj , features = keep_genes)

  return(sobj)

})

saveRDS(sobj_list, file.path(out.dir ,"sobj_list.rds"))
saveRDS(sobj_list.f, file.path(out.dir, "sobj_list_postqc.rds"))

###############################################################################|

hipsc = merge(sobj_list.f[[1]],sobj_list.f[[2]] )%>%SCTransform(., verbose = FALSE)
condition = hipsc@meta.data%>%select("orig.ident")%>%mutate(condition = case_when(orig.ident=="BIZH02_5533C4PB0WT_0" ~ "Control" ,
                                                                                  orig.ident=="BIZH02_5533C4PB7KO_0" ~ "PLXNB1-KO"))
hipsc = AddMetaData(hipsc, condition)

saveRDS(hipsc , file.path(out.dir ,"sobj_merged.rds"))

# UMAP and clustering
DefaultAssay(hipsc) = "SCT"
hipsc <- RunPCA(hipsc, verbose = FALSE)
hipsc <- RunUMAP(hipsc, dims = 1:30, verbose = FALSE)
hipsc <- FindNeighbors(hipsc, dims = 1:30, verbose = FALSE)
hipsc <- FindClusters(hipsc, verbose = FALSE,resolution = 0.5)

# Visualize marker expression
DefaultAssay(hipsc) = "RNA"
fig1 = DimPlot_scCustom(hipsc, label = TRUE)+DimPlot_scCustom(hipsc, label = TRUE, group.by = "orig.ident")
fig2 = FeaturePlot_scCustom(hipsc, features = c(ast.markers, "SOX9","GFAP","S100"))
fig3 = DotPlot_scCustom(hipsc, flip_axes = T,
                        features =  c(ast.markers,"SOX9","GFAP","S100","NFIB","PLXNB1","SMTN","SLC38A2",
                                      pull(Cell_type_Markers, Gene))%>%intersect(rownames(hipsc)))
pdf(file = file.path(fig.dir,"mssm_ipscs_control.pdf"), width = 12, height = 12)
print(cowplot::plot_grid(fig1, fig2, nrow = 2, rel_heights = c(1,3)))

print(cowplot::plot_grid(fig3, "", ncol = 2, nrow = 1, rel_width = c(1,0.5)))

dev.off()

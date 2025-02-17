#----------------------------------------------------------#|
#     Analysis of PLXNB1-KO response in mouse astrocytes 
#----------------------------------------------------------#|
#                Natacha Comandante-Lou 
#                 (nc3018@columbia.edu)
#----------------------------------------------------------#|

rm(list = ls())
setwd("6.PLXNB1_Validation/mouse")
reticulate::use_condaenv('~/conda/my-envs/magic' )
library(tidyverse)
library(Seurat)
library(scCustomize)
library(SeuratDisk)
library(ggpubr)
library(phateR)
library(ggrepel)

source('../hiPSCs/utils/utils_plxnb1_ko.R')
t=3
knn=5

conf_thres = 0.9

condition_oi = c("Control-WT","Control-PLXNB1-KO")
my_comparisons = list(c("Control-WT","Control-PLXNB1-KO"))

# A tibble: 2 × 2
# condition             n
# <fct>             <int>
#   1 Control-WT         1231
# 2 Control-PLXNB1-KO  2444
################################################################################|
# Extract Astrocytes -------
################################################################################|

query <- readRDS("/mnt/mfs/ctcn/team/natacha/PLXNB1-KO-mouse_Huang_et_al/Analysis/20241008/ref_mapping_res/n_features=4000/npcs=30_k=10/annotated_query.rds")
#ast = query%>%subset(., subset=predicted.state=="Astrocyte")%>%subset(., subset =predicted.state.score>0.6)
ast = query%>%subset(., subset=predicted.state=="Astrocyte")%>%subset(., subset =predicted.state.score>=conf_thres)

ast@meta.data = mutate(ast@meta.data, KO = case_when(grepl(x =condition,pattern="KO")==T ~ "PLXNB1-KO", TRUE ~ "WT"))
ast@meta.data  = mutate(ast@meta.data , AD = case_when(grepl(x =condition,pattern="APP")==T ~ "APP", TRUE ~ "Control"))
ast@meta.data$condition = paste0(ast@meta.data$AD , "-",ast@meta.data$KO)
ast@meta.data$condition = factor(ast@meta.data$condition, levels = c("Control-WT","Control-PLXNB1-KO","APP-WT","APP-PLXNB1-KO"))

# Subset ast -focusing on Ctrl only
ast = subset(ast, subset=condition%in%condition_oi)

ast@meta.data%>%group_by(condition)%>%tally()

min_cells = 10 #keep genes with at least 10cells
output_dir = sprintf("output/score_%0.1f_output_t=%d_knn_%d_min-cells=%d",conf_thres,t,knn, min_cells)
dir.create(output_dir,recursive = T)

DefaultAssay(ast) = "RNA"
# noramlize
ast = NormalizeData(ast, normalization.method = "LogNormalize", scale.factor = 10000)
# noramlize
ast = SCTransform(ast, assay = "RNA")

saveRDS(ast, "../seurat.obj/mouse_plxnb1_ko_ast.rds")


##############################################################################|
# MAGIC IMPUTATION ----------
###############################################################################|

default_assay = "RNA"


# magic to impute gene expression on mouse astrocytes
sobj.magic = run_magic_pipeline(ast,default_assay =default_assay, magic_output_dir= file.path(output_dir, "magic", default_assay),
                                gene_x_viz="Plxnb1", gene_y_viz="Smtn",
                                default_slot = "count",t_list = t,knn=knn,normalize = T, sqrt = T, min_cells = min_cells)



saveRDS(sobj.magic, file.path(output_dir, "magic", "sobj.magic.rds"))

##############################################################################|
# Calculate Module Score----------
###############################################################################|
# Add features based on gene-gene correlation analysis
output <- readRDS("/mnt/mfs/ctcn/team/natacha/Analysis/Compare_Ast_Cross-Systems/output_v2/output.rds")%>%filter(grouping=="mice-ast")
feature.list = output[["gene_list"]][[1]]
names(feature.list) = paste0("ast.10.",names(feature.list))



#Calculate summary score
for (default_module_assay in c(sprintf("MAGIC_knn.%d_t.%d",knn,t),"SCT","RNA")){
  sobj.magic= seurat_get_summary_scores(sobj.magic , feature.list,
                                   module.key = paste0("MODULE.",default_module_assay),
                                   method = "seurat",
                                   assay = default_module_assay,
                                   out.fig.dir = file.path(output_dir, "module_score",default_module_assay ),
                                   plot = F,
                                   plot.col.limits = NULL,
                                   plot.reduction = "umap")

}


##############################################################################|
# Save Dataframe --------
##############################################################################|
genes_oi = c("Plxnb1",feature.list[["ast.10.all"]]) # all astrocyte 10 markers plus plxnb1 for plotting
x_assay_oi = sprintf("MAGIC_knn.%d_t.%d",knn,t)
y_assay_oi = sprintf("MODULE.MAGIC_knn.%d_t.%d.seurat",knn,t)

output_df = sobj.magic@meta.data%>%rownames_to_column("cell_id")%>%
  full_join(.,FetchData(sobj.magic,vars = genes_oi, assay = x_assay_oi)%>%rownames_to_column("cell_id"), by = "cell_id")

tmp = DefaultAssay(sobj.magic)
DefaultAssay(object = sobj.magic) <- y_assay_oi
module_scores = FetchData(sobj.magic,vars = names(feature.list))%>%rownames_to_column("cell_id")
output_df = full_join(output_df,module_scores, by = "cell_id")
DefaultAssay(object = sobj.magic)<- tmp

saveRDS(output_df, file = file.path(output_dir,"output_magic_df.rds" ))
##############################################################################|
# Density Plot -------------
##############################################################################|
for( yvar in names(feature.list)){
  xvar = "Plxnb1"
  fig_dir = file.path(output_dir, "density",yvar)
  dir.create(fig_dir , recursive = T)
  x_assay_oi = sprintf("MAGIC_knn.%d_t.%d",knn,t)
  
  
  plot_density(sobj.magic,y_assay_oi = sprintf("MODULE.MAGIC_knn.%d_t.%d.seurat",knn,t), 
               x_assay_oi = x_assay_oi,yvar = yvar, xvar = xvar,
               yinter = NULL, xinter = NULL,gate.method ="median", gate.ref = "Control-WT",
               out_file = file.path(fig_dir,"magic_magic" ))

}

##############################################################################|
# Violin plot (Overall Signature)--------
##############################################################################|
DefaultAssay(sobj.magic)= sprintf("MODULE.MAGIC_knn.%d_t.%d.seurat",knn,t)

class_oi = "condition"
alternative = "two.sided"
assay_oi= sprintf("MODULE.MAGIC_knn.%d_t.%d.seurat",knn,t)

# Wilcox test
p_list1 = lapply(names(feature.list), function(v){
  jitter_boxplot_w_stats( class_oi,v,assay_oi = assay_oi,.sobj = sobj.magic, alternative =alternative,
                          save_file=NULL,method = "wilcox.test", my_comparisons, my_pal)
})

ggarrange(plotlist = p_list1, nrow =1, ncol = 3)%>%
  ggsave(plot = ., filename = file.path(output_dir,  sprintf("wilcox.test.boxplot_magic_%s.pdf", alternative)), width = 12, height = 6)

##############################################################################|
# Violin plot (Expression) --------
##############################################################################|

vars_oi_list = genes_oi

DefaultAssay(sobj.magic)= assay_oi
vars_oi_list = vars_oi_list%>%intersect(.,rownames(sobj.magic))
p5 = lapply(vars_oi_list, function(x){jitter_boxplot_w_stats( class_oi,x,assay_oi = assay_oi,.sobj = sobj.magic,alternative =alternative,
                                                              save_file=NULL,method = "wilcox.test",my_comparisons, my_pal)})%>%ggarrange(plotlist = ., nrow = 5, ncol = 7)

ggsave(plot = p5, filename = file.path(output_dir, paste0("wilcox.expr_boxplot_all_",assay_oi,"_",alternative,".pdf")), width = 20, height = 25)


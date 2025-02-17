#-----------------------------------------#|
# Train ligand PLSR on CUIMC1 data
#-----------------------------------------#|
#       Natacha Comandante-Lou 
#       (nc3018@columbia.edu)
#-----------------------------------------#|

rm(list=ls())
gc()
renv::load("2.NicheNet Analysis")
setwd("3.PLSR")
library(Seurat) 
library(tidyverse)
library(SeuratDisk)
library(ggpubr)
library(reshape2)
library(abind)
library(reshape2)
library(colorspace)
library(abind)
library(logger)

###############################################################################|
# LOAD DATA ------
###############################################################################|
args = commandArgs(trailingOnly = T)

nichenet_output_folder = args[1]

# Gene Expression
# A combined seurat objects with the microglia and astrocytes from Green et al

mic.ast = LoadH5Seurat('/mnt/mfs/ctcn/team/natacha/ROSMAP_snRNAseq/MultiNicheNet Analysis/resources/mic_ast_2024-02-26.h5Seurat')
DefaultAssay(mic.ast) = "RNA"
Idents(mic.ast)
mic.ast = subset(mic.ast,idents = c("Macrophages", "Monocytes"), invert = TRUE) #remove macrophages and monocytes before the analysis
mic.ast.437 = subset(mic.ast,subset = individual.in.landscape ==T) 

log_info("LogNormalizing..")
mic.ast.437 <- NormalizeData(mic.ast.437, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(mic.ast.437)
log_info("Scaling Data..")
mic.ast.437 <- ScaleData(mic.ast.437, features = all.genes)
log_info("Done preprocessing.")

# Get cell state frequency
freq.w = readRDS("1.snucRNAseq_Cell-state_Meta-analysis/cell-state_freq/freq_rds/Discovery/state_freq.rds")
# Frequency correlation
freq.cor = freq.w%>%select(-c("projid"))%>%cor(use="complete.obs")  #discard the entire row if an NA is present


plsr_output_folder = "ligand_plsr_output"
dir.create(plsr_output_folder,recursive = TRUE)
###############################################################################|
# PSEUDOBULK EXPRESSION -------
###############################################################################|

receiver = "Ast.10"
top_n_signals = 100

## Get ligands of interests

prior_table = read.csv(file = file.path(nichenet_output_folder,receiver,"prioritized_table.csv"))
prior_table.top = prior_table[c(1:top_n_signals),]

ligand_oi = prior_table.top%>%pull(ligand)%>%unique()
receptor_oi = prior_table.top%>%pull(receptor)%>%unique()

## Get population average ----------------------------------------------
### average data per projid, cell state ---------------------

# * Assumes that the data has been log normalized (assays = "RNA", slot = "data)
# * AverageExpression  = log(Mean(exp(log_norm_gene_exp)))

data.id.state = AverageExpression(mic.ast.437 ,features = unique(c(ligand_oi)),group.by = c("state","projid"),
                                  assays = "RNA",slot = "data",
                                  return.seurat = TRUE)
snames <- str_split(rownames(data.id.state @meta.data), "_",simplify = TRUE)
data.id.state  = AddMetaData(data.id.state , metadata = snames[,1], col.name = "state")
data.id.state  = AddMetaData(data.id.state , metadata = snames[,2], col.name = "projid")

df.id.state = as.data.frame(t(data.id.state@assays[["RNA"]]@data ))
df.id.state$state = snames[,1]
df.id.state$projid = snames[,2]%>%stringr::str_pad(.,8,pad = "0")


gc()
###############################################################################|
# CELL-STATE SPECIFIC LIGAND CORRELATION -----
###############################################################################|
## Cell-state Specific Ligands Correlation with Frequency --------
meta.data = mic.ast.437@meta.data
ncell_min = 0
# Select participants with non-zero Ast.10 cells
cell.count = meta.data%>%group_by(projid,state)%>%tally()
projid.use = filter(cell.count, state%in% receiver & n > ncell_min )%>%pull(projid)%>%as.character%>%stringr::str_pad(.,8,pad = "0")
df.id.state.w =  dcast(melt(df.id.state, id.vars=c("projid", "state")), projid~variable+state)%>%
  merge(select(freq.w, c("projid", receiver)), by = "projid")%>%
  filter(projid %in% projid.use)%>%
  mutate(log_freq = log(.data[[receiver]]*100+1))

# Correlation function
cor_func = function(x_var,y_var,.df,method = "pearson"){
  res = cor.test(.df[[x_var]],.df[[y_var]],method = method)
  
  output = data.frame(x_var = x_var, y_var = y_var, r = res[["estimate"]], p.value = res[["p.value"]])
  colnames(output)[colnames(output)=="r"] = sprintf("%s.r",method)
  return(output)
}

## Pearson Corr ------------
method = "pearson"
variable_list = colnames(select(df.id.state.w,-c("projid",receiver,"log_freq")))
y_var = "log_freq"
cor.res = lapply(variable_list, cor_func, y_var, df.id.state.w,method)

df.cor = as.data.frame(abind(lapply(cor.res, function(X) select(X, c("x_var","y_var","pearson.r","p.value"))), along = 1))
df.cor$pearson.r = as.numeric(df.cor$pearson.r)
df.cor$p.value = as.numeric(df.cor$p.value)
df.cor$ligand = as.factor(str_split_fixed(df.cor$x_var,"_",n = 2)[,1])
df.cor$sender = as.factor(str_split_fixed(df.cor$x_var,"_",n = 2)[,2])
df.cor = mutate(df.cor, p_val_category = case_when(p.value<=0.005 ~ 8,
                                                   p.value>0.005 & p.value<=0.05 ~ 4,
                                                   p.value>0.05 ~ 1))
saveRDS(df.cor,sprintf("%s/ligand-oi-in-sender_%s_corr_%s.rds",plsr_output_folder,method,receiver))


#order ligands and senders
df.cor$ligand = factor(df.cor$ligand, levels = ligand_oi)
state.list = df.cor$sender%>%unique()
state.split = stringr::str_split(state.list , "\\.",simplify = TRUE)
sid = gtools::mixedorder(state.split[,2])
state.levels = state.list[sid]
df.cor$sender = factor(df.cor$sender, levels = c(state.levels[state.levels%>%grep("Ast.",.)],state.levels[state.levels%>%grep("Mic.",.)]))
df.cor =  df.cor%>% mutate(pearson.thres = abs(pearson.r) >= 0.1)%>%
  mutate(alpha_scale = case_when(pearson.thres ~ 1, TRUE ~ 0))

# Bubble plot in Extended Data Figure S9A
p = ggplot(data = df.cor,aes(x=sender, y=ligand, size=p_val_category, fill=pearson.r)) +
  geom_point(shape=21, color="black",aes(alpha = alpha_scale)) +
  scale_alpha(range = c(0.1, 1), name="|R| >= 0.1",breaks = c(0.1,1),labels = c("|R| < 0.1","|R| >= 0.1")) +
  scale_size(range = c(1, 8), name="P-value",breaks = c(1,4,8),labels = c("> 0.05","0.005 - 0.05","< 0.005")) +
  scale_fill_continuous_diverging(palette = "Blue-Red 2", cmax = 200, l1 = 30, l2 = 100,
                                  limits = c(-1,1),
                                  name = "R" )+
  scale_y_discrete(limits=rev)+
  ggtitle(sprintf("%s Log-frequency Correlation",receiver))+
  theme_classic() +
  theme(title = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1,size = 12),axis.text.y = element_text(size = 12))
ggsave(p, filename = sprintf("%s/bubble_cell-state_specifc_ligand_pearsons_corr_with_%s.pdf",plsr_output_folder,receiver),device = "pdf", width = 10, height = 13.5 )

# Ligands expressed by senders that pass the following threshold will be considered when calculating ligand summary scores
df.cor.f = filter(df.cor, abs(pearson.r) >= 0.1 & p.value <= 0.05 )%>%filter(!sender %in% receiver)



###############################################################################|
# CONTRUCT MODEL IN/OUTPUT -----
###############################################################################|

## Calculate ligand summary scores --------
aggregate_expression = function(selected_ligand, .df){
  #average significant senders for a given ligand, separated by directionality
  var_list.pos = df.cor.f%>%filter(pearson.r > 0) %>% filter(ligand == selected_ligand)%>%pull(x_var) #pull all significantly correlated senders for a given ligand
  var_list.neg = df.cor.f%>%filter(pearson.r < 0) %>%filter(ligand == selected_ligand)%>%pull(x_var)

  #scale to zero mean across participants before averaging across senders
  .dfs = sapply(.df, function(x)(x-mean(x%>%na.omit())))%>%as.data.frame()
  xx.pos = select(.dfs, var_list.pos)
  xx.neg = select(.dfs, var_list.neg)
  
  xx.pos.mean = data.frame(rowMeans(xx.pos,na.rm=TRUE))
  xx.neg.mean = data.frame(rowMeans(xx.neg,na.rm=TRUE))
  
  sname.pos = str_split_fixed(colnames(xx.pos),n=2,pattern = "_")
  colnames(xx.pos.mean) = sprintf("%s (%s)",unique(sname.pos[,1]),paste0(sname.pos[,2],collapse  = "/"))
  
  sname.neg = str_split_fixed(colnames(xx.neg),n=2,pattern = "_")
  colnames(xx.neg.mean) = sprintf("%s (%s)",unique(sname.neg[,1]),paste0(sname.neg[,2],collapse  = "/"))
  
  output = cbind(xx.pos.mean,xx.neg.mean)
  
  return(output)
}


ligand_list = unique(df.cor.f$ligand)

df.plsr = as.data.frame(abind(lapply(ligand_list, aggregate_expression,df.id.state.w%>%select(-"projid"))))
df.plsr = df.plsr[,colSums(is.na(df.plsr))<1] #remove ligand columns that have at least on NaN
df.plsr$projid = df.id.state.w$projid 
df.plsr[[receiver]] = df.id.state.w[[receiver]]

### input -----
# ligand summary score
X = select(df.plsr, -c("projid",receiver))%>%`rownames<-`(df.plsr$projid)

### output -------
# log receiver frequency
Y =log(select(df.plsr,all_of(receiver))*100+1)%>%`rownames<-`(df.plsr$projid)


###############################################################################|
# RUN PLSR PIPELINE -----
###############################################################################|
source("/utils/run_plsr_pipeline.R")
pls_res = run_plsr_pipeline (X,Y,
                   tts.seed = 132,
                   test.size = 0.1,
                   npls = NULL,
                   n_splits = 5,
                   n_repeats = 100,
                   output.dir =plsr_output_folder)

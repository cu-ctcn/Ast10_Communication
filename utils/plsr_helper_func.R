require(Seurat)
require(tidyverse)
require(dplyr)
require(logger)
#############################################################################|
# Pseudobulk -------
#############################################################################|
GetPseudobulk = function(seuratObj, assay_oi = "RNA",normalization.method = "LogNormalize",
                         genes_oi,group.by.vars = c("state","projid"),
                         save = TRUE, output.dir, output.type = c("seurat","data.frame")){
  if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}
  #############################################################################|
  ## LOG -------
  #############################################################################|
  log_layout(layout_simple)
  log_file = file.path(output.dir , "get_pseudobulk_expr.log")
  if(file.exists(log_file)){
    log_warn("Overwriting previous log: ", log_file," \n")
    file.remove(log_file)
  }
  
  log_appender(appender_tee(file = log_file))
  
  
  #############################################################################|
  ## Average Expression-------
  #############################################################################|
  log_info("-----------------------------------------------")
  log_info("            PSEUDOBULK EXPRESSION             ")
  log_info("-----------------------------------------------")
  log_info("Getting Pseudobulk Expression...")
  log_info("* current default assay: ",  DefaultAssay(seuratObj))
  log_info("* set assay to: ", assay_oi)
  log_info("* normalization: ",normalization.method )
  
  ### average data per projid, cell state 
  start = Sys.time()
  seuratObj = NormalizeData(object = seuratObj, normalization.method = normalization.method,assay = assay_oi)
  #seuratObj = NormalizeData(seuratObj,scale.factor = 10000, do.log = TRUE,normalization.method = normalization.method,assay = assay_oi)
  DefaultAssay(seuratObj) = assay_oi

  # * Assumes that the data has been log normalized (assays = "RNA", slot = "data)
  # * AverageExpression  = log(Mean(exp(log_norm_gene_exp)))
  
  
  data = AverageExpression(seuratObj,features = unique(genes_oi),group.by = group.by.vars,
                                    assays = assay_oi,slot = "data",
                                    return.seurat = TRUE)
  snames <- stringr::str_split(rownames(data@meta.data), "_",simplify = TRUE)
  data  = AddMetaData(data , metadata = snames[,1], col.name = group.by.vars[1])
  data  = AddMetaData(data , metadata = snames[,2], col.name = group.by.vars[2])
  
  gc()
  
  df = data@assays[[assay_oi]]@data%>%t()%>%as.data.frame()
  df[[group.by.vars[1]]] = stringr::str_split(rownames(df),'_',simplify = TRUE)[,1]
  df[[group.by.vars[2]]]= stringr::str_split(rownames(df),'_',simplify = TRUE)[,2]
  
  df.w =  reshape2::dcast(reshape2::melt(df, id.vars=c(group.by.vars[2], group.by.vars[1])), as.formula(paste(group.by.vars[2], "~ variable + ",group.by.vars[1])))
  # some entries are NA because not all projid have all states
  gc()
  
  if (save & "seurat" %in% output.type){
    log_info("Saving seurat object ...")
    saveRDS(data, file = file.path(output.dir,"pseudobulk_expr_sobj.rds"))
    
  }
  
  if (save & "data.frame" %in% output.type){
    log_info("Saving data.frame ...")

    saveRDS(df.w, file = file.path(output.dir,"pseudobulk_expr_df.rds"))
    
  }
  end = Sys.time()
  log_success("Elapsed Time: ", round(end-start, 2) ," ", units(end-start))
  return(list(expr = data, df.expr = df.w))
}


GetPseudobulk_seurat_v5 = function(seuratObj, assay_oi = "RNA",normalization.method = "LogNormalize",
                         genes_oi,group.by.vars = c("state","projid"),
                         save = TRUE, output.dir, output.type = c("seurat","data.frame")){
  if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}
  #############################################################################|
  ## LOG -------
  #############################################################################|
  log_layout(layout_simple)
  log_file = file.path(output.dir , "get_pseudobulk_expr.log")
  if(file.exists(log_file)){
    log_warn("Overwriting previous log: ", log_file," \n")
    file.remove(log_file)
  }
  
  log_appender(appender_tee(file = log_file))
  
  
  #############################################################################|
  ## Average Expression-------
  #############################################################################|
  log_info("-----------------------------------------------")
  log_info("            PSEUDOBULK EXPRESSION             ")
  log_info("-----------------------------------------------")
  log_info("Getting Pseudobulk Expression...")
  log_info("* current default assay: ",  DefaultAssay(seuratObj))
  log_info("* set assay to: ", assay_oi)
  log_info("* normalization: ",normalization.method )
  
  ### average data per projid, cell state 
  start = Sys.time()
  seuratObj = NormalizeData(object = seuratObj, normalization.method = normalization.method,assay = assay_oi)
  #seuratObj = NormalizeData(seuratObj,scale.factor = 10000, do.log = TRUE,normalization.method = normalization.method,assay = assay_oi)
  DefaultAssay(seuratObj) = assay_oi
  
  # * Assumes that the data has been log normalized (assays = "RNA", slot = "data)
  # * AverageExpression  = log(Mean(exp(log_norm_gene_exp)))
  
  
  data = AverageExpression(seuratObj,features = unique(genes_oi),group.by = group.by.vars,
                           assays = assay_oi,slot = "data",
                           return.seurat = TRUE)
  snames <- stringr::str_split(rownames(data@meta.data), "_",simplify = TRUE)
  data  = AddMetaData(data , metadata = snames[,1], col.name = group.by.vars[1])
  data  = AddMetaData(data , metadata = snames[,2], col.name = group.by.vars[2])
  
  gc()
  
  
  df = data@assays[[assay_oi]]@layers[["data"]]%>%t()%>%as.data.frame()# Seurat v5
  rownames(df) = colnames(data)
  colnames(df) = rownames(data)
  df[[group.by.vars[1]]] = stringr::str_split(rownames(df),'_',simplify = TRUE)[,1]
  df[[group.by.vars[2]]]= stringr::str_split(rownames(df),'_',simplify = TRUE)[,2]
  
  df.w =  reshape2::dcast(reshape2::melt(df, id.vars=c(group.by.vars[2], group.by.vars[1])), as.formula(paste(group.by.vars[2], "~ variable + ",group.by.vars[1])))
  # some entries are NA because not all projid have all states
  gc()
  
  if (save & "seurat" %in% output.type){
    log_info("Saving seurat object ...")
    saveRDS(data, file = file.path(output.dir,"pseudobulk_expr_sobj.rds"))
    
  }
  
  if (save & "data.frame" %in% output.type){
    log_info("Saving data.frame ...")
    
    saveRDS(df.w, file = file.path(output.dir,"pseudobulk_expr_df.rds"))
    
  }
  end = Sys.time()
  log_success("Elapsed Time: ", round(end-start, 2) ," ", units(end-start))
  return(list(expr = data, df.expr = df.w))
}

#############################################################################|
# Get cell-state frequency -------
#############################################################################|
GetStateFreq_v4 = function(seuratObj=NULL,meta = NULL, group.by = "predicted.id", obs.key = "projid", rm_immune = F,
                           logfreq = TRUE, output.dir, save=T, h5ad = T){
  
  # Calculate relative cell-state proportions:
  # 
  if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}
  #############################################################################|
  ## LOG -------
  #############################################################################|
  library(reticulate)
  
  log_layout(layout_simple)
  log_file = file.path(output.dir , "get_freq.log")
  if(file.exists(log_file)){
    log_warn("Overwriting previous log: ", log_file," \n")
    file.remove(log_file)
  }
  
  log_appender(appender_tee(file = log_file))
  
  
  
  log_info("-----------------------------------------------")
  log_info("      CALCULATE CELL-STATE FREQUENCY           ")
  log_info("-----------------------------------------------")
  log_info("* Calculating cell-state frequency...")
  log_info("* Log transform frequency = ", logfreq)
  start = Sys.time()
  if (is.null(meta)){
    meta = seuratObj@meta.data
  }
  if (is.integer(meta$projid)){
    meta$projid = as.character(meta$projid)%>%stringr::str_pad(.,8,pad = "0")
  }
  # Regroup cell states as in Green et al
  abbv = str_split(meta[[group.by]], pattern = '\\.', simplify = T)
  meta$abbv = abbv[,1]
  meta = mutate(meta, cell.type = case_when( abbv == "Ast"~ "Astrocyte",
                                             abbv == "CD8+ T Cells" ~ "CD8+ T Cells",
                                             abbv %in% c("Arteriole", "End","Venule") ~ "Endothelial",
                                             abbv == "Erythrocytes" ~ "Erythrocytes",
                                             abbv == "Exc" ~"Excitatory Neurons",
                                             abbv == "Fib" ~"Fibroblast",
                                             abbv == "Inh" ~ "Inhibitory Neurons",
                                             abbv == "Macrophages" ~ "Macrophages",
                                             abbv == "Mic" ~ "Microglia",
                                             abbv == "Monocytes" ~ "Monocytes",
                                             abbv == "NK Cells" ~ "NK Cells",
                                             abbv == "Neutrophils" ~ "Neutrophils",
                                             abbv %in% c("COP", "MFOL","OPC") ~ "OPCs",
                                             abbv == "Oli" ~"Oligodendrocytes",
                                             abbv == "Peri" ~ "Pericytes",
                                             abbv == "SMC" ~ "SMC"
  ))
  
  meta = mutate(meta, grouping.by = case_when(cell.type =="Astrocyte" ~ "Astrocyte",
                                              cell.type =="Excitatory Neurons" ~ "Excitatory Neurons",
                                              cell.type %in% c("CD8+ T Cells","Erythrocytes","NK Cells","Neutrophils","Macrophages","Monocytes") ~ "Immune",
                                              cell.type =="Inhibitory Neurons" ~ "Inhibitory Neurons",
                                              cell.type %in% c("Microglia") ~ "Microglia",
                                              cell.type =="OPCs" ~ "OPCs",
                                              cell.type =="Oligodendrocytes" ~ "Oligodendrocytes",
                                              cell.type == "Endothelial" ~ "Endothelial",
                                              cell.type == "Pericytes" ~ "Pericytes",
                                              cell.type == "SMC" ~ "SMC",
                                              cell.type == "Fibroblast" ~ "Fibroblast"
  ))
  # Calculate participant-wise subpopulation-proportion (within cell-type)
  if (rm_immune){
    df <- meta %>% 
      filter(grouping.by != "Immune") %>% # remove NK cells and T-cells for which we do not have sufficient cells
      dplyr::count(grouping.by, cell.type, projid, .data[[group.by]]) %>%
      dplyr::group_by(grouping.by, projid) %>%
      mutate(prevalence=n/sum(n)) %>%
      mutate(log_freq =  log(prevalence*100+1))%>%
      ungroup()
    gc()
  }else{
    
    df <- meta %>% 
      dplyr::count(grouping.by, cell.type, projid, .data[[group.by]]) %>%
      dplyr::group_by(grouping.by, projid) %>%
      mutate(prevalence=n/sum(n)) %>%
      mutate(log_freq =  log(prevalence*100+1))%>%
      ungroup()
    gc()
  }
  
  
  # Creating proportion- and count matrices of participants' subpopulations
  freq <- reshape2::dcast(df, as.formula(sprintf("projid~%s",group.by)), value.var = "prevalence", fill = 0, fun.aggregate = sum) %>% tibble::column_to_rownames(.,"projid")
  log.freq <- reshape2::dcast(df, as.formula(sprintf("projid~%s",group.by)), value.var = "log_freq", fill = 0, fun.aggregate = sum) %>% tibble::column_to_rownames(.,"projid")
  counts <- reshape2::dcast(df, as.formula(sprintf("projid~%s",group.by)), value.var = "n", fill = 0, fun.aggregate = sum) %>% tibble::column_to_rownames("projid")
  
  
  
  if (h5ad){
    use_condaenv("BEYOND.ViaEnv")
    sc <- import("scanpy")
    np <- import("numpy")
    pd <- import("pandas")
    anndata = import ("anndata")
    ids <- rownames(freq)
    data <- anndata::AnnData(
      # Subpopulation proportions: participants (rows) over subpopulations (columns)
      X = freq ,
      
      # Counts and sqrt(prop) of subpopulations
      layers = list(counts = counts[ids, colnames(freq )],
                    sqrt.freq = sqrt(freq ),
                    log.freq = log.freq
      ),
      
      # General information about the different subpopulations
      var = df %>% dplyr::select(grouping.by, cell.type, group.by) %>% unique %>% column_to_rownames(group.by) %>% `[`(colnames(freq),),
      
      # # General information about the participants
      # obs = list(batches = donor.batches[ids,],
      #            main.batch = main.batch[ids,]),
      
      # Participant-wise QCs
      # obsm = list(QCs = qcs[ids,] %>% `rownames<-`(rownames(.) %>% as.character()))
    )
    
    anndata::write_h5ad(data, file.path(output.dir,"freq.h5ad"))
  }
  
  
  
  if(save){
    saveRDS(df, file.path(output.dir,"cell.annotations.rds"))
    freq$projid = rownames(freq)
    saveRDS(freq, file.path(output.dir,"state_freq.rds"))
    log.freq$projid = rownames(log.freq)
    saveRDS(log.freq, file.path(output.dir,"state_logfreq.rds"))
    counts$projid = rownames(counts)
    saveRDS(counts, file.path(output.dir,"counts.rds"))
    
  }
  
  end = Sys.time()
  log_success("Elapsed Time: ", round(end-start, 2) ," ", units(end-start))
  
  return(freq)
}

GetCellTypeFreq = function(seuratObj=NULL,meta = NULL, group.by = "predicted.id", obs.key = "projid",
                           logfreq = TRUE, output.dir, save){
  
  if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}
  #############################################################################|
  ## LOG -------
  #############################################################################|
  log_layout(layout_simple)
  log_file = file.path(output.dir , "get_freq.log")
  if(file.exists(log_file)){
    log_warn("Overwriting previous log: ", log_file," \n")
    file.remove(log_file)
  }
  
  log_appender(appender_tee(file = log_file))
  
  
  
  log_info("-----------------------------------------------")
  log_info("      CALCULATE CELL-TYPE Frequency           ")
  log_info("-----------------------------------------------")
  log_info("* Calculating cell-type frequency...")
  log_info("* Log transform frequency = ", logfreq)
  start = Sys.time()
  if (is.null(meta)){
    meta = seuratObj@meta.data
  }
  
  count = group_by(meta, .data[[obs.key]], .data[[group.by ]])%>%
    count(.drop = FALSE)%>%spread( .data[[obs.key]], n, fill = 0 )%>%
    gather(!!obs.key, n, -!!group.by)
  count$subset = str_split(count[[group.by]], pattern = '\\.', simplify = T)[,1]#substr(count[[group.by]],1,3)
  freq = group_by(count, .data[[obs.key]])%>%mutate(freq = n/sum(n)) #frequencies normalized to projid
  

  if(logfreq){freq = mutate(freq, log.freq = log(freq*100+1))}
  
  if(!logfreq){freq.w = reshape2::dcast(freq, as.formula(sprintf("%s ~ %s ",obs.key, group.by)),value.var="freq")}
  if(logfreq){freq.w = reshape2::dcast(freq, as.formula(sprintf("%s ~ %s ",obs.key, group.by)),value.var="log.freq")}
  
  if(save & logfreq){
    saveRDS(freq.w, file.path(output.dir,"state_logfreq.rds"))
  }
  if(save & !logfreq){
    saveRDS(freq.w, file.path(output.dir,"state_freq.rds"))
  }
  
  end = Sys.time()
  log_success("Elapsed Time: ", round(end-start, 2) ," ", units(end-start))
  return(freq.w)
}

#############################################################################|
# Prepare PLSR input -------
#############################################################################|
## Calculate ligand summary scores --------
aggregate_expression = function(selected_ligand, df.cor.f,.df){
  var_list.pos = df.cor.f%>%filter(pearson.r > 0) %>% filter(ligand == selected_ligand)%>%pull(x_var) #pull all significantly correlated senders for a given ligand
  var_list.neg = df.cor.f%>%filter(pearson.r < 0) %>%filter(ligand == selected_ligand)%>%pull(x_var)
  
  #scale to zero mean for each variable
  .dfs = sapply(.df, function(x)(x-mean(x%>%na.omit())))%>%as.data.frame()
  xx.pos = select(.dfs, var_list.pos)
  xx.neg = select(.dfs, var_list.neg)
  
  #Average scaled mean across variables
  xx.pos.mean = data.frame(rowMeans(xx.pos,na.rm=TRUE))
  xx.neg.mean = data.frame(rowMeans(xx.neg,na.rm=TRUE))
  
  sname.pos = str_split_fixed(colnames(xx.pos),n=2,pattern = "_")
  colnames(xx.pos.mean) = sprintf("%s (%s)",unique(sname.pos[,1]),paste0(sname.pos[,2],collapse  = "/"))
  
  sname.neg = str_split_fixed(colnames(xx.neg),n=2,pattern = "_")
  colnames(xx.neg.mean) = sprintf("%s (%s)",unique(sname.neg[,1]),paste0(sname.neg[,2],collapse  = "/"))
  
  output = cbind(xx.pos.mean,xx.neg.mean)
  
  return(output)
}


AssembleInput = function(cor.res, expr.df, method = "mean.scaled",receiver = "Ast.10",
                         r.thres = 0.1, p.thres = 0.05){
  
  #############################################################################|
  ## LOG -------
  #############################################################################|
  log_layout(layout_simple)
  log_file = file.path(output.dir , "assemble_input.log")
  if(file.exists(log_file)){
    log_warn("Overwriting previous log: ", log_file," \n")
    file.remove(log_file)
  }
  
  log_appender(appender_tee(file = log_file))
  
  log_info("-----------------------------------------------")
  log_info("            ASSEMBLE PLSR LIGAND INPUT             ")
  log_info("-----------------------------------------------")
 
  log_info("* method: ",  method)
  log_info("* correlation threshold: ", r.thres )
  log_info("* p-value threshold: ",p.thres )
  log_info("* receiver: ", receiver )
  
  start = Sys.time()
  # significant threshold
  df.cor.f = filter(cor.res, abs(pearson.r) >= r.thres & p.value <= p.thres )%>%filter(!sender %in% receiver)
  ligand_list = unique(df.cor.f$ligand)
  
  df.plsr = as.data.frame(abind::abind(lapply(ligand_list, aggregate_expression,df.cor.f,expr.df%>%select(-"projid"))))
  #df.plsr = df.plsr[,colSums(is.na(df.plsr))<1] #remove ligand columns that have at least on NaN
  df.plsr$projid = expr.df$projid 
  
  end = Sys.time()
  log_success("Elapsed Time: ", round(end-start, 2) ," ", units(end-start))
  return(df.plsr)
}

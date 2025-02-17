#-------------------------------------------------------------------#|
# Meta-analysis of associations between cell-states and AD traits
#-------------------------------------------------------------------#|
#                   Natacha Comandante-Lou 
#                   (nc3018@columbia.edu)
#-------------------------------------------------------------------#|
rm(list = ls())
setwd("/Code/1.snucRNAseq_Cell-state_Meta-analysis")
library(tidyverse)
library(dplyr)
library(scCustomize)
library(ggpubr)
library(colorspace)

output_folder = "trait_associations_final"
dir.create(output_folder)
purple2orange <- colorRampPalette(c("#3a0ca3", "#5d77d5", "white", "#ffb627", "#e36414"))
pal_ref = c(`CUIMC 1 (Discovery)`="#A79CEB", `Diversity` ="#DCE25A",`MIT`="#FF3B53",`CUIMC 2` ="#13DDC6",`Replication`="#26A9E0",`Overall`="black")

freq.list = readRDS("Code/1.snucRNAseq_Cell-state_Meta-analysis/cell-state_freq/snuc_state_freq.rds")
# [[1]] : metadata
# [[2]] : freq
# [[3]] : sqrt freq
# [[4]] : log freq
dataset_levels = c("CUIMC 1 (Discovery)","CUIMC 2","MIT","Diversity")
full_meta = freq.list[[1]]
freq_df = freq.list[[3]] # use sqrt freq

full_meta$reference = factor(full_meta$reference, levels = dataset_levels )
freq_df$reference = factor(freq_df$reference, levels = dataset_levels )


donor_oi_discov = full_meta%>%filter(reference=="CUIMC 1 (Discovery)")%>%pull(projid)
donor_oi_multiome = full_meta%>%filter(reference=="CUIMC 2")%>%pull(projid)
donor_oi_mit = full_meta%>%filter(reference=="MIT")%>%pull(projid)
donor_oi_diversity_other = full_meta%>%filter(reference=="Diversity")%>%pull(projid)


run_trait_association = function(ref_oi, full_meta, freq_df){
  print(ref_oi)
  if(ref_oi == "CUIMC 1 (Discovery)"){donor_oi = donor_oi_discov}
  if(ref_oi == "CUIMC 2"){donor_oi = donor_oi_multiome}
  if(ref_oi == "MIT"){donor_oi = donor_oi_mit}
  if(ref_oi == "Diversity"){donor_oi = donor_oi_diversity_other}
  metadata = full_meta%>%filter(reference%in%ref_oi & projid%in%donor_oi)
  metadata$age_death = as.numeric(metadata$age_death)
  query.freq.log = freq_df%>%filter(reference%in%ref_oi & projid%in%donor_oi)
  query.freq.log[is.na(query.freq.log)] = 0
  
  ############################################################################|
  
  
  cov_oi = c("age_death","msex","pmi","study","educ","Cell_Count","Gene_Median_Genes_Detected")
  var_oi = c("projid",cov_oi,"age_bl","amyloid","tangles","cogng_demog_slope")
  
  state_oi_list = colnames(query.freq.log)%>%setdiff(.,c("projid","reference","CD8+ T Cells","Erythrocytes","NK Cells","Neutrophils"))
  
  metadata <- metadata%>%dplyr::select(var_oi)%>%
    merge(.,query.freq.log%>%dplyr::select(c("projid",state_oi_list)), by = "projid")%>%
    mutate(sqrt_amyloid = sqrt(.data[["amyloid"]]))%>%
    mutate(sqrt_tangles= sqrt(.data[["tangles"]]))
  metadata = metadata[rowSums(is.na(metadata%>%dplyr::select(var_oi)))==0,] #remove donors without patho information
  
  # count number of indv per state
  counts = lapply(state_oi_list, function(x){
    .df = metadata%>%dplyr::select(x)>0
    return(data.frame(state = x, n_indiv_per_state = sum(.df)))
  })%>%do.call("rbind",.)
  ######################################|
  # Linear Regression -----
  ######################################|
  y_var_list  = c("sqrt_amyloid","sqrt_tangles","cogng_demog_slope")
  
  lm_func = function(y_var, state_oi, cov_oi, .df ){
    
    x_var = paste0(c(cov_oi,state_oi),collapse = " + ")
    print(x_var)
    fit = lm(as.formula(sprintf("%s ~ %s", y_var,x_var)), .df)
    res = summary(fit)[["coefficients"]]%>%as.data.frame()
    res$terms = rownames(res)
    
    output = filter(res, terms == state_oi)
    if(dim(output)[1]==0){
      tmp = colnames(output)
      output[1,] = rep(NA,length(tmp))
      output$terms = state_oi
      output$y_var = y_var
    }else{
      output$y_var = y_var
    }
    return(output)
    
  }
  
  lm.res = lapply(state_oi_list, function(x){lapply(y_var_list,lm_func, x, cov_oi, metadata)%>%do.call("rbind",.)})%>%
    do.call("rbind",.)%>%group_by(y_var)%>%
    mutate(p.signif = case_when(`Pr(>|t|)` <= 0.0001 ~ "****",
                                    `Pr(>|t|)` <= 0.001 & `Pr(>|t|)` > 0.0001 ~ "***",
                                    `Pr(>|t|)` <= 0.01 & `Pr(>|t|)` > 0.001 ~ "**",
                                    `Pr(>|t|)` <= 0.05 & `Pr(>|t|)` > 0.01 ~ "*",
                                    TRUE ~ ""
    ))%>%
    mutate(p.adj = p.adjust(.data[["Pr(>|t|)"]], method = "BH"))%>%ungroup()%>%
    mutate(p.adj.signif = case_when(p.adj <= 0.0001 ~ "****",
                                p.adj <= 0.001 & p.adj > 0.0001 ~ "***",
                                p.adj <= 0.01 & p.adj > 0.001 ~ "**",
                                p.adj <= 0.05 & p.adj > 0.01 ~ "*",
                                TRUE ~ ""
    ))
  lm.res$y_var = factor(lm.res$y_var, levels = y_var_list )
  lm.res$n = dim(metadata)[1]
  lm.res$state = lm.res$terms
  lm.res = merge(lm.res, counts, by = "state")
  lm.res$reference = ref_oi
  return(lm.res)
}
       

lm.res = lapply(dataset_levels, run_trait_association,  full_meta, freq_df)%>%do.call("rbind",.)

saveRDS(lm.res, file.path(output_folder,"summary_stats_ad_patho-state_asso.rds"))


#################################################|
# Meta-analysis (inverse variance based) -------
#################################################|
library(meta)
metagen_wrap = function(df, effect_var = "Estimate",  se_var = "Std. Error", study_label = "reference", method.tau = 'REML'){
  # wrapper for the metagen package
  meta_analysis <- metagen(
    TE = df[[effect_var]],    # Effect
    seTE = df[[se_var]], # Standard errors
    data = df,
    studlab = df[[study_label]], # Study labels
    method.tau =  method.tau  # Random-effects model
  )

    meta_analysis_df = tibble(terms = unique(df$terms),
                              y_var = unique(df$y_var),
                               SEi = meta_analysis$seTE.random,
                               Bi = meta_analysis$TE.random,
                               wi = NA,
                               Zi = meta_analysis$zval.random,
                               P = meta_analysis$pval.random,
                               model = "random.effect",
                               P_hetero = meta_analysis$pval.Q,
                               I2 = meta_analysis$I2
                               )%>%rbind(
                       tibble(terms = unique(df$terms),
                              y_var = unique(df$y_var), 
                              SEi = meta_analysis$seTE.fixed,
                              Bi = meta_analysis$TE.fixed,
                              wi = NA,
                              Zi = meta_analysis$zval.fixed,
                              P = meta_analysis$pval.fixed,
                              model = "fixed.effect",
                              P_hetero = meta_analysis$pval.Q,
                              I2 = meta_analysis$I2
                               
                       ))

  return(meta_analysis_df)

}



#################################################|
# Meta-analysis (inverse variance based) -------
#################################################|

run_meta_analysis = function(.lm.res, min_n_indiv_per_state = 50, levels_oi =  c("Overall","CUIMC 1 (Discovery)","MIT","CUIMC 2","Diversity")){

  # remove states that are less than n=50
  lm.res.f = filter(.lm.res, n_indiv_per_state> min_n_indiv_per_state)
  
  dff = lm.res.f%>% mutate(Bi = Estimate)%>%mutate(SEi = `Std. Error`)%>%
    mutate(wi = 1/(SEi^2))%>%
    mutate(Zi= Bi/SEi) %>% mutate(P = `Pr(>|t|)`)%>%dplyr::select(terms, y_var, SEi, Bi,wi,  Zi, P,reference,n_indiv_per_state)
  

  # Meta analysis
  dff.overall = group_by(dff, terms, y_var)%>%group_split()%>%
    lapply(.,metagen_wrap, effect_var = "Bi",  se_var = "SEi", study_label = "reference", method.tau = 'REML')%>%do.call("rbind",.,)

  
  # Benjimini-Hochberg correction
  dff.overall = group_by(dff.overall, y_var, model)%>%mutate(p.adj = p.adjust(.data[["P"]], method = "BH"))%>%
    mutate(p.signif = case_when(p.adj <= 0.0001 ~ "****",
                                    p.adj <= 0.001 & p.adj > 0.0001 ~ "***",
                                    p.adj <= 0.01 & p.adj > 0.001 ~ "**",
                                    p.adj <= 0.05 & p.adj > 0.01 ~ "*",
                                    TRUE ~ ""
    ))%>%
    mutate(reference = "Overall")
  
  
  dff = group_by(dff, reference, y_var)%>% mutate(p.adj = p.adjust(.data[["P"]], method = "BH"))%>%
    mutate(p.signif= case_when(p.adj <= 0.0001 ~ "****",
                               p.adj <= 0.001 & p.adj > 0.0001 ~ "***",
                               p.adj <= 0.01 & p.adj > 0.001 ~ "**",
                               p.adj <= 0.05 & p.adj > 0.01 ~ "*",
                               TRUE ~ ""))%>%ungroup
  
  df.meta = plyr::rbind.fill(dff, dff.overall)
  df.meta$reference = factor(df.meta$reference, levels = levels_oi)

  df.meta = mutate(df.meta, cell.type  = case_when(grepl(x=.data[["terms"]], pattern = "Ast.") ~ "Astrocytes",
                                                   grepl(x=.data[["terms"]], pattern = "Mic.") ~ "Microglia",
                                                   grepl(x=.data[["terms"]], pattern = "Inh.") ~ "InhibitoryNeurons",
                                                   grepl(x=.data[["terms"]], pattern = "Exc.") ~ "ExcitatoryNeurons",
                                                   grepl(x=.data[["terms"]], pattern = "O.") ~ "Oligodendrocytes",
                                                   grepl(x=.data[["terms"]], pattern = "End.")|.data[["terms"]]%in%c("Arteriole","Venule") ~ "Endothelial",
                                                   TRUE ~ "Other"
  ))
  
  return(df.meta)
}


## Run meta analysis including discovery --------
df.meta = run_meta_analysis(.lm.res = lm.res, min_n_indiv_per_state = 50)
saveRDS(df.meta, file.path(output_folder,"meta-analysis_ad_patho-state_asso_v2.rds"))
## Run meta analysis without discovery (replication only)--------
df.meta.rep = run_meta_analysis(.lm.res = lm.res%>%filter(!reference%in%"CUIMC 1 (Discovery)"), min_n_indiv_per_state = 50, 
                                levels_oi = c("Overall","MIT","CUIMC 2","Diversity"))
saveRDS(df.meta.rep, file.path(output_folder,"meta-analysis_ad_patho-state_asso_v2_rep.rds"))



#################################################|
# Heatmap -------
#################################################|

model_oi = "fixed.effect"

meta_res =df.meta%>%
  filter((!reference=="Overall")|(model==model_oi))
meta_res_rep = df.meta.rep%>%
  filter((!reference=="Overall")|(model==model_oi))

# pull list of states that are significant in at least one pathology
state_oi.1 = filter(meta_res,  reference%in%"Overall" & (!p.signif==""))%>%pull(terms)
state_oi.2 = filter(meta_res_rep,  reference%in%"Overall" & (!p.signif==""))%>%pull(terms)
state_oi = unique(c(state_oi.1, state_oi.2))

plot_cluster_map = function(.meta_res, state_oi ,ref_oi = "Overall", row_order = NULL ,savefile = "clustermap.pdf", fig_title = "Meta-analysis (with Discovery)"){
  save.folder = str_split(savefile,"/", simplify = T)%>%.[-length(.)]%>%paste(.,collapse = "/")
  dir.create(save.folder, recursive = T)
  df.plot = filter(.meta_res, reference%in%ref_oi & terms%in% state_oi)
  df.plot.w = reshape2::dcast(df.plot, terms ~ y_var,value.var = "Zi")
  rownames(df.plot.w) = df.plot.w$terms
  
  p.lab = reshape2::dcast(df.plot, terms ~ y_var,value.var = "p.signif")
  rownames(p.lab) = p.lab$terms
  
  if (is.null(row_order)){
    df.plot.w = df.plot.w[,-1]
    p.lab = p.lab[,-1]
  }else{
    df.plot.w = df.plot.w[match(row_order, df.plot.w$terms),-1]
    p.lab = p.lab[match(row_order, p.lab$terms),-1]
  }
  library("gplots")
  #my_pal = bluered(100)
  #my_pal = rev(diverging_hcl( n = 100,palette = "Purple-Green"))
  my_pal = rev(purple2orange(100))
  Rowv = ifelse(is.null(row_order),T,F)
  mat = df.plot.w%>%as.matrix
  
  pdf(savefile,width = 6, height = 8)
  hm = heatmap.2(mat, scale = "none", Colv = "none", Rowv=Rowv, col = my_pal, cellnote = p.lab,notecol="black",dendrogram = "row",
                 key.xlab = "t-stat",lhei=c(0.6,4), lwid=c(1.5,3.5), keysize=0.7, key.par = list(cex=0.6),breaks=seq(-7,7,14/100),
                 sepwidth =c(0.0005,0.0005),sepcolor = "black",colsep=1:ncol(mat),
                 rowsep=1:nrow(mat),
                 trace = "none", density.info = "none",margins = c(8,15),cexCol = 1.2)
  print(hm)
  title(fig_title)
  dev.off()
  
  
  return(hm)
}



# Heatmap (ranked from top)-----

row_order = reshape2::dcast(meta_res%>%filter(reference%in%"Overall" & terms%in% state_oi)%>%mutate(signed.P = sign(Bi)*(-log10(P))), terms ~ y_var,value.var = "signed.P")%>%
  arrange(-cogng_demog_slope,sqrt_tangles,sqrt_amyloid)%>%
  distinct(terms)%>%pull(terms)
plot_cluster_map(meta_res, state_oi ,row_order = rev(row_order), savefile = file.path(output_folder,"ranked_by_pval",paste0(model_oi,"_p-ranked_hm_purple-orange.pdf")), fig_title = paste0("Meta-analysis (with Discovery) - ", model_oi))
plot_cluster_map(meta_res_rep, state_oi ,row_order = rev(row_order), savefile = file.path(output_folder,"ranked_by_pval",paste0(model_oi,"_p-ranked_hm_rep_purple-orange.pdf")), fig_title = paste0("Meta-analysis (Replication) - ", model_oi))

# Heatmap (Discovery)---
plot_cluster_map(meta_res, state_oi ,row_order = rev(row_order), ref_oi = "CUIMC 1 (Discovery)", savefile = file.path(output_folder,"ranked_by_pval",paste0(model_oi,"_p-ranked_hm_purple-orange_Discovery.pdf")), fig_title = "Discovery")
plot_cluster_map(meta_res, state_oi ,row_order = rev(row_order), ref_oi = "MIT", savefile = file.path(output_folder,"ranked_by_pval",paste0(model_oi,"_p-ranked_hm_purple-orange_MIT.pdf")), fig_title = "MIT")


##############################################################################|
# Forest plot ------
##############################################################################|
state_oi = "Ast.10"
df.plot = filter(meta_res_rep, reference=="Overall")%>%mutate(reference = "Replication")%>% #pull out metanalysis results on replication data only
  rbind(., meta_res)
df.plot$reference = factor(df.plot$reference,levels = c("Overall","Replication","Diversity","MIT","CUIMC 2","CUIMC 1 (Discovery)"))

y_var_list = unique(meta_res$y_var)


custom_forest = function(y, .meta_res){
  # add 0.95 CI:
  .meta_res = mutate(.meta_res, lower = Bi-1.96*SEi, upper = Bi+1.96*SEi)%>%filter(y_var==y)
  ggplot(.meta_res, aes(x=Bi, xmin=lower, xmax=upper, y=reference, shape = is.meta, fill=reference)) +
    #geom_pointrange(show.legend=FALSE, color = "black") + 
    geom_linerange( size=1, aes(color=reference),show.legend=FALSE) +
    geom_point(size = 3, color = "black",show.legend=FALSE)+
    geom_vline(xintercept=0, lty=2) +
    scale_fill_manual(values = pal_ref) +
    scale_color_manual(values = pal_ref) +
    xlab(paste0(" Effect size (",model_oi,")")) +
    scale_shape_manual(values = c( 21L,23L))+
    ggtitle(y) +
    theme_light()
}

p2 = lapply(y_var_list,
            custom_forest,
            df.plot%>%filter(terms%in%state_oi)%>%
              mutate(is.meta = case_when(reference %in% c("Overall","Replication") ~ "Y", TRUE ~ "N")))
pdf(file.path(output_folder,paste0(model_oi,"_forest_plot_ast.10.pdf")), width = 10, height = 2)
ggpubr::ggarrange(plotlist = p2, common.legend = TRUE, nrow = 1, ncol = length(y_var_list))%>%print()
dev.off()

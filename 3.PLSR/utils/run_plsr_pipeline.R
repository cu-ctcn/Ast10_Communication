#' PLSR pipeline
#' - Plot model performance (R2 and Q2) with varying # of PLS components
#' - Perform cross-validation
#' - Compare model predictions vs data in held out test set
require(reticulate)
require(logger)
require(tidyverse)
require(ggpubr)

use_condaenv("/home/nc3018/conda/my-envs/sklearn-env")
sk = import("sklearn")
np = import("numpy")
jb = import("joblib")

run_plsr_pipeline = function(X,Y,
                             tts.seed = 132,
                             test.size = 0.2,
                             npls = 1,
                             n_splits = 5,
                             n_repeats = 10,
                             
                             plot.scores.pc = c(1,2),
                             
                             plot.loadings.pc = c(1,2),
                             output.dir = "plsr_output"
){
  #############################################################################|
  # LOG --------
  #############################################################################|
  if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}
  log_file = file.path(output.dir, "run_plsr_pipeline.log")
  
  if(file.exists(log_file)){
    log_warn("Overwriting previous log: ", log_file)
    file.remove(log_file)
  }
  
  log_appender(appender_tee(file = log_file))
  
  log_info("|##################################################################|")
  log_info("          Partial Least-Squares Regression Modeling               ")
  log_info("|##################################################################|")
  
  
  log_info(date())
  
  ### Splitting data into training and test sets ------------------------
  log_info("Splitting data into training and test sets: ", test.size*100, "% test data")
  tts = sk$model_selection$train_test_split(X, Y, test_size=test.size, random_state=as.integer(tts.seed),shuffle =TRUE )
  X_train = tts[[1]]
  X_test = tts[[2]]
  Y_train = tts[[3]]
  Y_test = tts[[4]]
  
  
  
  ### Model Performace -----------------------------------------------------------
  log_info("Calculating model performance across pls components... ")
  start = Sys.time()
  
  log_info(" * Calculating R2... ")
  #Full Model R2
  r2_scores = list()
  Xz_train = sk$preprocessing$StandardScaler()$fit_transform(X_train)%>%`colnames<-`(colnames(X_train))%>%`rownames<-`(rownames(X_train))
  Yz_train = sk$preprocessing$StandardScaler()$fit_transform(Y_train)%>%`colnames<-`(colnames(Y_train))%>%`rownames<-`(rownames(Y_train))
  
  
  for (k in c(1:dim(X_train)[2])){
    pls = sk$cross_decomposition$PLSRegression(n_components=k)
    # fit model
    pls$fit(Xz_train, Yz_train)
    # evaluate model
    Yhat = pls$predict(Xz_train)
    # store
    #r2_score = sk$metrics$explained_variance_score(Yz_train, Yhat)
    r2_score = sk$metrics$r2_score(Yz_train, Yhat)
    #print(r2_score)
    r2_scores[[k]] = r2_score
    
  }
  
  # K-fold Q2
  log_info(" * Calculating Q2 using repeated ", n_splits, " fold cross-validation, repeat ", n_repeats, " times...")
  q2_scores = list()
  
  for (k in c(1:dim(Xz_train)[2])){
    pls = sk$cross_decomposition$PLSRegression(n_components=k)
    
    #kf = sk$model_selection$KFold(n_splits = as.integer(10),shuffle = TRUE)
    rkf = sk$model_selection$RepeatedKFold(n_splits=as.integer(n_splits), n_repeats=as.integer(n_repeats), random_state=as.integer(2652124))
    
    kf_list <- iterate(rkf$split(Xz_train,Yz_train))
    
    
    Y_pred = list()
    
    for (idx in  c(1:length(kf_list))){
      train_ix = kf_list[[idx]][[1]]+as.integer(1)
      test_ix = kf_list[[idx]][[2]]+as.integer(1)
      # fit model
      pls$fit(Xz_train[train_ix, ],Yz_train[train_ix])
      # evaluate model
      Yhat = pls$predict((Xz_train[test_ix,]))
      
      # store
      Y_pred[[idx]] = data.frame("Yhat" = Yhat, "Replicate" = ceiling(idx/n_splits), "ID" = test_ix, "Y" =Yz_train[test_ix,] )
      
      
      
    }
    
    Y_pred = do.call("rbind",Y_pred)
    
    
    q2_score = lapply(c(1:n_repeats), function(r){sk$metrics$r2_score(Y_pred%>%filter(Replicate%in%r)%>%pull(Y), Y_pred%>%filter(Replicate%in%r)%>%pull(Yhat))})%>%
      do.call("rbind",.)
    
    q2_scores[[k]] = data.frame(mean = mean(q2_score), sd = sd(q2_score))
  }
  
  Q2 = do.call("rbind",q2_scores)
  
  # Plot model performance
  log_info("Plotting...")
  model_performance = data.frame(R2 = c(0,do.call(rbind,r2_scores)*100),Q2_mean = c(0,Q2$mean*100), Q2_sd = c(0, Q2$sd*100), k = c(0:dim(Xz_train)[2]))
  p1 = ggplot(data = model_performance)+
    geom_line(aes(x = k, y = R2), color = 'blue')+
    geom_point(aes(x = k, y = R2), color = 'blue', size = 2.5)+
    geom_line(aes(x = k, y = Q2_mean), color = 'red')+
    geom_errorbar(aes(x = k, ymin = Q2_mean-Q2_sd, ymax = Q2_mean+Q2_sd), color = 'red')+
    geom_point(aes(x = k, y = Q2_mean), color = 'red', size = 2.5)+
    ylab('Percent variance')+
    xlab('Number of PLS components')+
    ylim(c(0,100))+
    theme_classic()
  ggsave(p1, filename = file.path(output.dir,"model_performance.pdf"),device = 'pdf',width = 3.5, height = 3.5,dpi = 300)
  end = Sys.time()
  log_info("[Model Performance] Elapsed: ", round(end-start, 2) ," ", units(end-start))
  
  
  ### Score plots-------
  pls = sk$cross_decomposition$PLSRegression(n_components=min(dim(X_train)))
  # fit model
  pls$fit(Xz_train, Yz_train)
  x_scores = as.data.frame(pls$x_scores_)%>%`rownames<-`(rownames(X_train))%>%`colnames<-`(paste0("PLS ", c(1:min(dim(X_train)))))
  y_scores = as.data.frame(pls$y_scores_)%>%`rownames<-`(rownames(Y_train))%>%`colnames<-`(paste0("PLS ", c(1:min(dim(X_train)))))
  x_scores_scaled =as.data.frame(sapply(x_scores, function(x){x/(max(x)-min(x))}))
  y_scores_scaled =as.data.frame(sapply(y_scores, function(x){x/(max(x)-min(x))}))
  
  if (!is.null(plot.scores.pc)){
    xi = colnames(x_scores_scaled)[plot.scores.pc[1]]
    yi = colnames(x_scores_scaled)[plot.scores.pc[2]]
    df.plot = x_scores_scaled%>%select(all_of(c(xi,yi)))%>%cbind(Y_train)
    
    response = colnames(Y_train)[1]
    p2 = ggplot(df.plot,aes(x =.data[[xi]], y = .data[[yi]], fill = .data[[response]]))+
      geom_point(shape = 21, size = 3,color = 'black' )+
      scale_fill_viridis_c()+
      ggtitle(label = 'PLSR Scores')+
      ylab(yi)+
      xlab(xi)+
      theme_classic()
    ggsave(p2,filename = file.path(output.dir, sprintf("x-scores_scaled_%s_vs_%s.pdf",xi,yi)),device = 'pdf',width = 4.12, height = 3.47,dpi = 300)
    
  }
  
  
  
  
  
  ### Pick optimal pls -------------------------------------------------------
  if (is.null(npls)){
    # Determine optimal PLS component
    df.plot = model_performance[-1,]
    df.plot$r2.diff = diff(model_performance$R2)
    df.plot$q2.diff = diff(model_performance$Q2_mean)
    df.plot$ratio = df.plot$Q2_mean/df.plot$R2                
    
    # metric 1: difference between r2-q2
    co1 = min(filter(df.plot, ratio == max(df.plot$ratio))%>%pull(k))
    log_info("Minimum PC where Q2 captures most of the R2: PC ", co1)
    
    ## metric 2: max Q2
    
    co2 <- min(filter(df.plot, Q2_mean == max(df.plot$Q2_mean))%>%pull(k))
    
    # last point where change of % of variation is more than 0.1%.
    log_info("Minimum PC where Q2 is the highest: PC ", co2)
    
    npls = min(co1, co2)
    
    log_info("Optimal PC selected: PC ", npls)
    
    
  }
  ### Final eval-------------------------------------------------------------------
  log_info("Final Evaluation on training and test data...")
  k_final = as.integer(npls)
  pls_final = sk$cross_decomposition$PLSRegression(n_components=k_final)
  response = colnames(Y)
  
  
  ## Training data (k-fold cross-val prediction vs data) -------
  pls_crossval = sk$cross_decomposition$PLSRegression(n_components=k_final)
  
  rkf = sk$model_selection$RepeatedKFold(n_splits=as.integer(n_splits), n_repeats=as.integer(n_repeats), random_state=as.integer(2652124))
  
  kf_list <- iterate(rkf$split(Xz_train,Yz_train))
  
  Y_pred = list()
  Y_data = list()
  
  for (idx in  c(1:length(kf_list))){
    train_ix = kf_list[[idx]][[1]]+as.integer(1)
    test_ix = kf_list[[idx]][[2]]+as.integer(1)
    # fit model
    pls_crossval$fit(Xz_train[train_ix, ],Yz_train[train_ix])
    # evaluate model
    Yhat = pls_crossval$predict((Xz_train[test_ix,]))
    
    # store
    Y_pred[[idx]] = data.frame("Yhat" = Yhat, "Replicate" = ceiling(idx/n_splits), "ID" = test_ix)
    Y_data[[idx]] = Yz_train[test_ix,]
    
  }
  
  #print(q2_score)
  Y_pred = do.call("rbind",Y_pred)
  Y_data = do.call("rbind",lapply(Y_data,as.data.frame))%>%`colnames<-`("Y")
  exp_var_train = sk$metrics$r2_score(Y_data, Y_pred$Yhat)
  train_prediction = data.frame(data = Y_data$Y, model.prediction = Y_pred$Yhat, replicate = Y_pred$Replicate, id = Y_pred$ID)
  
  # Summarize model prediction across n_repeats
  data_summary <- function(data, varname, groupnames){
    
    summary_func <- function(x, col){
      c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE))
    }
    data_sum<-plyr::ddply(data, groupnames, .fun=summary_func,
                          varname)
    data_sum <- plyr::rename(data_sum, c("mean" = varname))
    return(data_sum)
  }
  
  train_prediction_summary = data_summary(train_prediction, "model.prediction", "id")%>%
    merge(data_summary(train_prediction, "data", "id"), by="id")
  
  p3 = ggscatter(train_prediction_summary, x = "data", y = "model.prediction",
                 add = "reg.line", conf.int = FALSE,
                 add.params = list(color = "black"),
                 cor.coef = TRUE, cor.method = "pearson", color = "purple", alpha = 0.5, size = 2.5) + # pearson coefficient (R)
    geom_abline(slope = 1, linetype = "dashed")+
    geom_errorbar(data = train_prediction_summary,aes(ymin=model.prediction-sd.x, ymax=model.prediction+sd.x, x = data),color = "purple", alpha = 0.5)+
    annotate("text",label = sprintf("variance \npredicted = %0.2f%%",exp_var_train*100), x = -0.1, y = -1)+
    ggtitle(label = response)+
    ylab('Model Prediction')+
    xlab('Training Data')+
    theme_classic()
  ggsave(p3,filename = file.path(output.dir, sprintf("training_performance_npls=%d.pdf",k_final)),device = 'pdf',width = 3.34, height = 3.18,dpi = 300)
  
  
  
  ## Test data ----
  Xz_test = sk$preprocessing$StandardScaler()$fit(X_train)$transform(X_test)
  Yz_test = sk$preprocessing$StandardScaler()$fit(Y_train)$transform(Y_test)
  
  pls_final$fit(Xz_train,Yz_train)
  
  
  
  Yhat_test = pls_final$predict(Xz_test)
  exp_var_test = sk$metrics$r2_score(Yz_test, Yhat_test)
  test_prediction = data.frame(data = as.numeric(Yz_test), model.prediction = Yhat_test)
  
  
  p4 = ggscatter(test_prediction, x = "data", y = "model.prediction",
                 add = "reg.line", conf.int = FALSE, 
                 add.params = list(color = "black"),
                 cor.coef = TRUE, cor.method = "pearson", color = "blue", alpha = 0.5, size = 2.5) + # pearson coefficient (R)
    geom_abline(slope = 1, linetype = "dashed")+
    #stat_cor(method = "pearson",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),color = "blue")+ # R^2
    annotate("text",label = sprintf("variance \npredicted = %0.2f%%",exp_var_test*100), x = -0.1, y = -0.6)+
    #coord_cartesian(ylim = c(-2, 2), xlim = c(-2,2))+
    ggtitle(label = response)+
    ylab('Model Prediction')+
    xlab('Test Data')+
    theme_classic()
  
  ggsave(p4, filename = file.path(output.dir, sprintf("test_performance_npls=%d.pdf",k_final)),device = 'pdf',width = 3.34, height = 3.18,dpi = 300)
  
  
  ### Variable Importance (VIP Scores)--------------------------------------------
  log_info("Calculating VIP scores....")
  source_python('3.PLSR/utils/vip.py')
  pearson_func = function(x_var,y_var,.df){
    res = cor.test(.df[[x_var]],.df[[y_var]])
    return(c(x_var, y_var, res[["estimate"]][["cor"]], res[["p.value"]]))
  }
  variable_list = colnames(X)
  df.plsr = cbind(X,Y)
  
  y_var = colnames(Y)
  
  df.cor = as.data.frame(t(as.matrix(sapply(variable_list,pearson_func,y_var,df.plsr))))
  colnames(df.cor) = c("x_var","response","Pearson's r","P-value")
  df.cor$`Pearson's r` = as.numeric(df.cor$`Pearson's r`)
  df.cor$`P-value` = as.numeric(df.cor$`P-value`)
  df.cor$sign = sign(df.cor$`Pearson's r`)
  
  df_vip = data.frame(x_var = colnames(Xz_train), vip = vip(pls_final))
  
  df_vip = dplyr::left_join( df.cor,df_vip)
  df_vip$sign = sign(df_vip$`Pearson's r`)
  
  plot_vip = function(.df){
    df2 = .df%>%mutate(abs_vip = abs(`vip`))%>%mutate(vip.cat = case_when(.data[["vip"]]*.data[["sign"]]>=1 ~ "VIP>=1",
                                                                          .data[["vip"]]*.data[["sign"]]<=-1 ~ "VIP<-1",
                                                                          TRUE ~ "-1 < VIP < 1"))
    df2$vip.cat = factor( df2$vip.cat,levels = c("VIP>=1","-1 < VIP < 1","VIP<-1"))
    
    
    mypal = c("#F94144","#d3d3d3","#277da1")
    ggplot(df2, aes(x = reorder(x_var, abs_vip),y = .data[["vip"]]*.data[["sign"]],fill=as.factor(.data[["vip.cat"]])))+
      geom_bar(data = df2,stat="identity",color = "black")+
      geom_hline(yintercept = 0)+
      geom_hline(yintercept = 1,color = "gray",linetype = "dashed")+
      geom_hline(yintercept = -1,color = "gray",linetype = "dashed")+
      scale_fill_manual(values = mypal[sort(as.numeric(unique(df2$vip.cat)))],name = "")+
      ylab("VIP scores")+
      xlab(sprintf("Predictors of %s ",response))+
      coord_cartesian(ylim = c(-2, 2))+
      #ylim(c(-3,3))+
      theme(axis.text.y = element_text(hjust = 1,size = 14),axis.title.y = element_text(hjust = 0.5,size = 16),
            axis.text.x = element_text(angle = 45, hjust = 1,size = 9),axis.title.x = element_text(hjust = 0.5,size = 16),
            legend.text = element_text(size = 10),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            strip.text = element_text(size = 16))
  }
  
  p5 = plot_vip(df_vip)
  
  
  ggsave(p5, filename = sprintf("%s/vip_response=%s_npls=%d.pdf",output.dir, response,k_final),device = 'pdf',width = length(unique(df_vip$x_var))/2, height = 10,dpi = 300)
  
  
  ### Loading plots ---------
  #Note loadings are not weights (weights are used to calculate VIP, measures predictability of Y; loadings are used to project to PCs)
  x_var_labels = stringr::str_split(colnames(X_train),pattern = "\ ", simplify = TRUE)[,1]%>%make.unique(sep='*')
  x_loadings = as.data.frame(pls$x_loadings_)%>%`rownames<-`(x_var_labels)%>%`colnames<-`(paste0("PLS ", c(1:min(dim(X_train)))))
  x_loadings$type = "x_loadings"
  y_loadings = as.data.frame(pls$y_loadings_)%>%`rownames<-`(colnames(Y_train))%>%`colnames<-`(paste0("PLS ", c(1:min(dim(X_train)))))
  y_loadings$type = "y_loadings"
  loadings = rbind(x_loadings, y_loadings)
  if (!is.null(plot.loadings.pc)){
    require(ggrepel)
    xi = colnames(loadings)[plot.loadings.pc[1]]
    yi = colnames(loadings)[plot.loadings.pc[2]]
    df.plot = loadings%>%select(all_of(c(xi,yi,"type")))
    df.plot$var = rownames(df.plot)
    
    
    x_lim = max(abs(df.plot[[xi]]))+0.1
    y_lim = max(abs(df.plot[[yi]]))+0.1
    
    sig_var = stringr::str_split(filter(df_vip, vip >=1)%>%pull(x_var), pattern = "\ ", simplify = TRUE)[,1]
    p6 = ggplot()+
      geom_hline(yintercept =  0, linetype= "dashed", color = "gray")+
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray")+
      geom_point(data = df.plot,aes(x =.data[[xi]], y = .data[[yi]], fill = type),shape = 21, size = 2,color = 'gray' , fill = 'gray')+
      geom_point(data = filter(df.plot, var%in% c(sig_var)),aes(x =.data[[xi]], y = .data[[yi]]),fill = "blue",shape = 21, size = 2,color = 'black' )+
      geom_text_repel(data = filter(df.plot, var%in% c(sig_var, response)), aes(x =.data[[xi]], y = .data[[yi]],label = .data[["var"]]),nudge_x = 0.01, force= 0.5 )+
      scale_fill_viridis_d()+
      ggtitle(label = 'PLSR Loadings')+
      ylab(yi)+
      xlab(xi)+
      # xlim(-x_lim, x_lim)+
      # ylim(-y_lim, y_lim)+
      coord_cartesian(ylim = c(-y_lim, y_lim),xlim = c(-x_lim, x_lim))+
      theme_classic()
    ggsave(p6,filename = file.path(output.dir, sprintf("loadings_%s_vs_%s.pdf",xi,yi)),device = 'pdf',width = 6.3, height = 4.8,dpi = 300)
  }
  
  
  ### Save output-------
  jb$dump(pls_final, sprintf("%s/pls_final_response=%s.joblib",output.dir, response)) 
  scaler_x = sk$preprocessing$StandardScaler()$fit(X_train)
  scaler_y = sk$preprocessing$StandardScaler()$fit(Y_train)
  jb$dump(scaler_x, sprintf("%s/scaler_x.joblib",output.dir)) 
  jb$dump(scaler_y, sprintf("%s/scaler_y.joblib",output.dir)) 
  
  
  jb$dump(pls_final, sprintf("%s/pls_final_response=%s.joblib",output.dir, response)) 
  output = list("model" = pls_final,
                "model_performance" = model_performance,
                "data" = df.plsr,
                "vip_scores" = df_vip,
                "final_npls" = k_final,
                "loadings" = loadings,
                "x_scores" = x_scores,
                "y_scores" = y_scores,
                "x_scores_scaled" = x_scores_scaled,
                "y_scores_scaled" = y_scores_scaled
  )
  log_info("Saving....")
  saveRDS(output,file = sprintf("%s/plsr_res_response=%s.rds",output.dir, response))
  return(output)
}
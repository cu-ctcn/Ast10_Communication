#' Wrapper function for MAGIC
#' @import Rmagic
#' @import logger
#' @export
#' 
require(Rmagic)
require(logger)

run_magic = function(t, knn,data, genes, #genes to run MAGIC on
                     plot=TRUE,
                     plot_only = TRUE,
                     init,
                     gene_x_viz, gene_y_viz, #genes to visualize
                     save.as = c("seurat","magic"),
                     solver = 'exact',
                     sobj, #seuart object to save
                     output.dir = "magic_output",
                     new.assay.name = "MAGIC",
                     save.obj.path = NULL
){
  
  if(!dir.exists(output.dir)){dir.create(output.dir, recursive = T)}
  log_info("#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#")
  log_info("                MAGIC              ")
  log_info("#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#")
  
  log_info("Parameters:")
  log_info(" * knn = ",knn)
  log_info(" * t = ",t)
  
  if (!length(init)>1){
    log_info("No input MAGIC object, calculate MAGIC without initialization...")
    init = NULL
  }
  if (!plot_only){
    start = Sys.time()
    data_MAGIC <- magic(data, genes=genes,init = init,knn = knn, t = t, solver = solver, seed = 3751)
    end = Sys.time() 
    log_info("Time Elapse: ",round(end-start,2)," ",units(end-start))
    
  }

  if (plot|plot_only){
    
    log_info("Plotting...")
    
    p = ggplot(data_MAGIC%>%as.data.frame()) +
      geom_point(aes(.data[[gene_x_viz]],.data[[gene_y_viz]], color=.data[[gene_y_viz]])) +
      viridis::scale_color_viridis(option="B")+
      ggtitle(sprintf("t = %d, knn = %d",t,knn))+
      theme_classic()
    ggsave(p, filename = file.path(output.dir,sprintf("knn=%d_t=%d",knn,t),sprintf("magic_scatter__%s_vs_%s_t=%d_knn=%d.png",gene_x_viz, gene_y_viz, t, knn)),width = 6, height = 5,create.dir = T)
    

  }
  
  output = list()
  output$t = t
  output$knn = knn
  output$data_MAGIC = data_MAGIC
  
  if (!is.null(save.as)){
   
    if ("magic" %in% save.as){
      log_info("Saving magic object...")
      if(!dir.exists(file.path(output.dir,sprintf("knn=%d_t=%d",knn,t)))){dir.create(file.path(output.dir,sprintf("knn=%d_t=%d",knn,t)), recursive = TRUE)}
      saveRDS(output,file.path(output.dir,sprintf("knn=%d_t=%d",knn,t),"magic_output.rds"))
      
    }
    if ("seurat" %in% save.as){
      log_info("Add magic to seurat object, saving...")
      
      magic_assay <- CreateAssayObject(data = data_MAGIC[["result"]]%>%as.matrix()%>%t())
      
      # add this assay to the previously created Seurat object
      sobj[[new.assay.name]] <- magic_assay
      Misc(sobj, slot = "magic") = data_MAGIC$operator
      
      if (!is.null(save.obj.path)){
        if(!dir.exists(file.path(output.dir,sprintf("knn=%d_t=%d",knn,t)))){dir.create(file.path(output.dir,sprintf("knn=%d_t=%d",knn,t)), recursive = TRUE)}
        save.obj.path = file.path(output.dir,sprintf("knn=%d_t=%d",knn,t),"sobj_after_magic.rds")
        saveRDS(sobj,save.obj.path)
        }

      
      output$seurat = sobj
    }
    
  }
  return(output)
}

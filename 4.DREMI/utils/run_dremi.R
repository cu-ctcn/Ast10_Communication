#' Plot DREMI output
#' @import ggplot2
#' 

#' k : int, range=[0:n_samples), optional (default: 10)
#'     Number of neighbors
#' n_bins : int, range=[0:inf), optional (default: 20)
#'          Number of bins for density resampling
#' n_mesh : int, range=[0:inf), optional (default: 3)
#'          In each bin, density will be calculcated around (mesh ** 2) points
#reticulate::use_condaenv("/home/nc3018/conda/my-envs/sctools-3.8")
reticulate::source_python('/mnt/mfs/ctcn/team/natacha/Packages/scprep-custom/scprep/knnDREMI_ncl.py')

run_dremi = function(X,Y, k = 10, n_bins = 32, n_mesh = 2,dremi = NULL,
                     output.dir = "dremi.output",
                     plot = TRUE, plot_only = FALSE, xlab,ylab, limits = c(0, 0.35), scale.axis = TRUE){
  
  set.seed(5654)
  if(!dir.exists(output.dir)){dir.create(output.dir, recursive = TRUE)}
  if(!plot_only|is.null(dremi)){
    k = as.integer(k)
    n_bins = as.integer(n_bins)
    n_mesh = as.integer(n_mesh)
    
    dremi = knnDREMI(x = X, y = Y, k = k, n_bins = n_bins, n_mesh = n_mesh,
                     plot = TRUE, return_drevi = TRUE)
    
    names(dremi) = c("dremi","density","x_bins_z","y_bins_z","x_bins","y_bins")
  }

  
  if (plot|plot_only){
   p =  plotDREMI(dremi, xlab, ylab, limits = limits,scale.axis = scale.axis)
   dir.create(file.path(output.dir,"dremi.heatmap"),recursive = TRUE)
   ggsave(plot = p,filename =file.path(output.dir,"dremi.heatmap",sprintf("dremi_x=%s_y=%s.pdf",xlab,ylab)), width = 4.5, height = 3.9 )
   dremi[["plot"]] = p
   }
  return(dremi)
}



plotDREMI = function(dremi, xlab, ylab, limits = c(0,0.3),scale.axis = TRUE){
  if (scale.axis){  xi = 3; yi = 4} else {xi = 5;  yi = 6}

  
  x_bins = dremi[[xi]]
  y_bins = dremi[[yi]]
  bin_width_x = x_bins[2]-x_bins[1]
  bin_width_y = y_bins[2]-y_bins[1]
  c_density = dremi[[2]]
  
  # Center point of each bin
  x_bins_center = x_bins[1:length(x_bins)-1]+bin_width_x/2
  y_bins_center = y_bins[1:length(y_bins)-1]+bin_width_x/2
  
  
  mesh_pts= plot3D::mesh(x_bins_center, y_bins_center)
  
  df.dremi = data.frame(x = c(t(mesh_pts[["x"]])), y = c(t(mesh_pts[["y"]])), c_density = c(c_density))
  
  
  ggplot() +
    geom_tile(data = df.dremi, aes(x = x, y = y, fill = c_density)) +
    xlab(xlab)+
    ylab(ylab)+
    ggtitle(sprintf("DREMI score: %0.2f",dremi[[1]]))+
    scale_fill_viridis_c(limits = limits, oob = scales::squish) +
    theme_classic()
  

}

#####################
# Test

# dremi1 = sc$stats$knnDREMI(x=X, y=Y,k = as.integer(10),n_bins = as.integer(20),plot = TRUE, return_drevi = TRUE,
#                            filename='test.pdf')

# plotDREMI = function(dremi, xlab, ylab, limits = c(0,0.3)){
#   xi = 5; yi = 6
#   
#   x_bins = dremi[[xi]]
#   y_bins = dremi[[yi]]
#   bin_width_x = x_bins[2]-x_bins[1]
#   bin_width_y = y_bins[2]-y_bins[1]
#   c_density = dremi[[2]]
#   
#   # Center point of each bin
#   x_bins_center = x_bins[1:length(x_bins)-1]+bin_width_x/2
#   y_bins_center = y_bins[1:length(y_bins)-1]+bin_width_x/2
#   
#   
#   mesh_pts= plot3D::mesh(x_bins_center, y_bins_center)
#   
#   df.dremi = data.frame(x = c(t(mesh_pts[["x"]])), y = c(t(mesh_pts[["y"]])), c_density = c(c_density))
#   
#   #Fit
#   # 
#   # max_idx = apply(c_density,2,which.max)
#   # yy = max_idx
#   # scale_factor = (max(y_bins)-min(y_bins))/length(y_bins)
#   # xx = x_bins_center
#   # fit <- nls(yy ~ SSlogis(xx, Asym, xmid, scal), data = data.frame(xx,yy))
#   # df_pred = data.frame(x_pred = seq(min(xx), max(xx), length.out = 100),
#   #                      y_pred = min(y_bins)+scale_factor*predict(fit, data.frame(xx = seq(min(xx), max(xx), length.out = 100))))
#   # 
#   # plot(xx, yy)
#   # lines(seq(min(xx), max(xx), length.out = 100), 
#   #       predict(fit, data.frame(xx = seq(min(xx), max(xx), length.out = 100))))
#   
#   ggplot() +
#     geom_tile(data = df.dremi, aes(x = x, y = y, fill = c_density)) +
#     #geom_point(data = df_pred, aes(x = x_pred, y=y_pred),color = "white", size = 0.3)+
#     xlab(xlab)+
#     ylab(ylab)+
#     # xlim(c(-0.5,0.7))+
#     # ylim(c(-0.5, 2.5))+
#     # xlim(c(-3,4))+
#     # ylim(c(-0.55, 10))+
#     ggtitle(sprintf("DREMI score: %0.2f",dremi[[1]]))+
#     scale_fill_viridis_c(limits = limits, oob = scales::squish) +
#     #scale_fill_continuous_sequential("plasma",rev = FALSE)+
#     theme_classic()
#   
#   # 
#   # df.test = data.frame(x = seq(min(xx), max(xx), length.out = 100), y1 = test1, y2 = test2) 
#   # ggplot()+
#   #   geom_point(data = df.test, aes(x = x, y = y1),color = "blue")+
#   #   geom_point(data = df.test, aes(x = x, y = y2),color = "red")+
#   #   theme_classic()
# }

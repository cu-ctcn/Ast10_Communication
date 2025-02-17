require(nichenetr)
require(ggforce)
require(ggnewscale)
require(shadowtext)
require(cowplot)
require(gplots)



purple2orange <- colorRampPalette(c("#3a0ca3", "#9a93e7", "white", "#ffb627", "#fb8f3f"))
blue <- colorRampPalette(c("white","#dbf0ff","blue","#3a0ca3"))
orange <- colorRampPalette(c("white","#fcf6bd","#FFb500","#e36414"))
green <-colorRampPalette(c("white","#AECC9A","#195500"))


###############################################################################|
# Functions ------------
###############################################################################|
plot_ligand_activity = function(prior.tbl, var_oi = "activity_zscore"){

   ggplot(prior.tbl%>%select(c("ligand", var_oi) )%>%distinct(ligand,.keep_all = T),
          aes(y = ligand, x = 1, fill = .data[[var_oi]]))+
          geom_tile(color = "black",na.rm = F)+
          theme_classic()+
          scale_fill_gradientn(colors = green(100))+
          theme( axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
                 axis.ticks.x  = element_blank()
                 )
  
  
}

plot_ligand_receptor = function(prior.tbl, var_oi = "lfc_receptor"){

  ggplot(prior.tbl%>%select(ligand, receptor,var_oi ), aes(y = ligand, x = receptor, fill = .data[[var_oi]]))+
    geom_tile(color = "black",na.rm = F)+
    theme_classic()+
    scale_fill_gradientn(colours =  orange(100))+
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 7))

}

dotplot_expression = function(prior.tbl, ligands_oi, fill.var = "lfc"){
  .df.l = filter(prior.tbl, ligand%in%ligands_oi)%>%
    select(ligand, sender,scaled_avg_exprs_ligand, lfc_ligand, avg_ligand, pct_expressed_sender)
  .df.l$ligand = factor(.df.l$ligand, levels = rev(ligands_oi))
  var_ligand = colnames(.df.l)[.df.l%>%colnames()%>%str_detect(fill.var)]
  p1 = ggplot(.df.l, aes(x = sender, y = ligand, fill = .data[[var_ligand]], size = pct_expressed_sender))+
    geom_point(shape = 21)+
    scale_fill_gradientn(colors = purple2orange(100), limits = c(-1,1),oob = scales::oob_squish)+
    scale_size_binned_area(max_size = 6)+
    #scale_y_reverse()+
    theme_classic()+
    theme(axis.text.x = element_text(angle=45, hjust=1),legend.position = "top")

  return(p1)
}

###############################################################################|


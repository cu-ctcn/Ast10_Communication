#----------------------------------------------------------#|
#      Analysis of intensity data from CellProfiler 
#   Single-cell model to predict SLC39A11 (ZIP11) level
#----------------------------------------------------------#|
#                Natacha Comandante-Lou 
#                 (nc3018@columbia.edu)
#----------------------------------------------------------#|
rm(list = ls())
setwd("6.PLXNB1_Validation/immunofluorescence")

library(tidyverse)
library(ggpubr)
library(afex)
library(readxl)
data = read_excel("/mnt/mfs/ctcn/team/natacha/Manuscript/Code/SupplementaryTables/Supplementary Tables/Table.S8_CellProfiler_Analysis.xlsx",
                  sheet = "Raw CellProfiler Output"
                  )
data$ImageNumber = as.factor(data$ImageNumber)

scale_factor = 2^16-1 # for 16 bit images
fig_dir = "figures"
dir.create(fig_dir, recursive = T)


###############################################################################|
# Plot 2d density ---------------
###############################################################################|

plot_2d_density = function(df, xvar, yvar, var_type = "Intensity" ,scale_factor = 2^16-1, split.by = NULL){
  xvar_plot = paste0(var_type, "_", xvar)
  yvar_plot = paste0(var_type, "_", yvar)
  
  df[[xvar]] = df[[xvar_plot]]*scale_factor
  df[[yvar]] = df[[yvar_plot]]*scale_factor
  p = ggplot(df, aes(x = .data[[xvar]], y = .data[[yvar]]))+
    geom_point(alpha = 0.5, size = 0.3, color = "#36454F")+
    stat_density_2d(geom = "polygon",
                    aes(alpha = after_stat(level), fill = ..level..),
                    bins = 20) +
    scale_alpha_binned(range = c(0.1,0.7))+
    scale_fill_viridis_c(option = "turbo", oob=scales::squish)+
    scale_x_continuous(transform = "log10")+
    scale_y_continuous(transform = "log10")+
    theme_classic(base_size = 14)
  
  if (!is.null(split.by)){
    p = p+facet_wrap(split.by)
  }
  return(p)
}

p1 = plot_2d_density(data, xvar = "MeanIntensity_GFAP", yvar = "MeanIntensity_PLXNB1")
p2 = plot_2d_density(data, xvar = "MeanIntensity_GFAP", yvar = "MeanIntensity_ZIP11")
p3 = plot_2d_density(data, xvar = "MeanIntensity_PLXNB1", yvar = "MeanIntensity_ZIP11")



fig = ggarrange(p1,p2,p3, nrow = 1, ncol = 3)
ggsave(file.path(fig_dir, "intensity_2d_density.pdf"), fig, width = 12, height = 7.5)



###############################################################################|
# Linear Mixed Effects model ---------------
###############################################################################|
library(lme4)
library(stats)
library(afex)
# generate models
# log10-transform intensity then z-score
data.scaled <- data%>%mutate(across(c(Intensity_MeanIntensity_GFAP, Intensity_MeanIntensity_PLXNB1, Intensity_MeanIntensity_ZIP11),
                                    function(df) log10(df*scale_factor)))%>%
  mutate(across(c(Intensity_MeanIntensity_GFAP, Intensity_MeanIntensity_PLXNB1, Intensity_MeanIntensity_ZIP11),
                function(df) (df-mean(df))/sd(df)))


m0.lm <- lm(Intensity_MeanIntensity_ZIP11 ~ 1, data = data.scaled )
m0.lmer = lmer(Intensity_MeanIntensity_ZIP11 ~ 1 + (1|ImageNumber), REML = T, data = data.scaled )

AIC(logLik(m0.lm))
AIC(logLik(m0.lmer))

# m0.lmer has the lowest AIC - include random effect is justified
# > AIC(logLik(m0.lm))
# [1] 8675.552
# > AIC(logLik(m0.lmer))
# [1] 6274.816

# Final model:
# - Fixed effect - GFAP, PLXNB1
# - Random effect - image, slopes and intercepts of GFAP and PLXNB1

# Use afex to get p-value
m.lmer <- mixed(Intensity_MeanIntensity_ZIP11 ~  (1+ Intensity_MeanIntensity_GFAP +Intensity_MeanIntensity_PLXNB1|ImageNumber) + Intensity_MeanIntensity_GFAP +Intensity_MeanIntensity_PLXNB1, 
                data=data.scaled, check_contrasts = FALSE, test_intercept = TRUE, method="KR")
summary(m.lmer)

saveRDS(m.lmer, file = file.path(fig_dir, "final_model.rds"))

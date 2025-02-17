#------------------------------------------------#|
# Calculate statistics for dataset demographics
#------------------------------------------------#|
#       Natacha Comandante-Lou 
#       (nc3018@columbia.edu)
#------------------------------------------------#|
rm(list=ls())
gc()

library(gt)
library(gtExtras)
library(readr)
library(tidyverse)
library(ggpubr)
library(colorspace)
library(cowplot)
library(ggbeeswarm)
library(gtsummary)
library(ggpubr)
library(rstatix)
source('utils_plt.R')

fig_dir = "data_description"
dir.create(fig_dir)

meta.data = readRDS("/Code/1.snucRNAseq_Cell-state_Meta-analysis/cell-state_freq/all_meta.rds")%>%
  filter( !reference=="Diversity-NHW") #focus on participants other than non-hispanic white in the diversity data
meta.data$reference = droplevels(meta.data$reference)

#-----------------------------------------------------------------------|
# Table: Demographics Statistical Comparison
#-----------------------------------------------------------------------|
summary.tbl2 <-
  meta.data |>
  select( sex,age_death, pathoAD, sqrt_amyloid, sqrt_tangles, cogng_demog_slope, reference) |>
  mutate(pathoAD = case_when(pathoAD==1 ~"AD", TRUE ~ "Non-AD"))%>%
  rename("Sex"="sex","Pathological AD" = "pathoAD","Age Death"="age_death", "Amyloid"="sqrt_amyloid","Tangles" ="sqrt_tangles", "Cognitive Slope"="cogng_demog_slope")|>
  tbl_summary(
    by = reference,
    missing = "no",
    statistic = all_continuous() ~ "{mean} ± {sd}",
    type = all_dichotomous() ~ "categorical"
  )|>
  add_p(all_continuous()~"oneway.test")%>%
  modify_header(label = "**Variable**") |> # update the column header
  italicize_levels() |>
  as_gt()%>%
  fmt_scientific(columns = "p.value",decimals = 2)%>%
  tab_header(title = "Demographics")%>%
  gt_theme_538_v2()
summary.tbl2 %>%gtsave(file.path(fig_dir,"summary_tbl.html"))


##############################################################################|
# Figures ------
##############################################################################|

meta.data = mutate(meta.data, sqrt_amyloid = sqrt(amyloid))%>%mutate(sqrt_tangles = sqrt(tangles))    
# Figure 1 Age Distribution
figure1 = boxplot_quasirnd(y_var = "age_death", x_var = "reference", yrange = NULL,meta.data,dataset_pal)
ggsave(plot = figure1,filename = file.path(fig_dir, "age_death.pdf"), width = 6, height = 5)

# Figure 2-4 Pathologies
figure2 = boxplot_quasirnd(y_var = "sqrt_amyloid", x_var = "reference", yrange = NULL,meta.data,dataset_pal)
ggsave(plot = figure2,filename = file.path(fig_dir, "sqrt_amyloid.pdf"), width = 6, height = 5)

figure3 = boxplot_quasirnd(y_var = "sqrt_tangles", x_var = "reference", yrange = NULL,meta.data,dataset_pal)
ggsave(plot = figure3,filename = file.path(fig_dir, "sqrt_tangles.pdf"), width = 6, height = 5)

figure4 = boxplot_quasirnd(y_var = "cogng_demog_slope", x_var = "reference", yrange = NULL,meta.data,dataset_pal)
ggsave(plot = figure4,filename = file.path(fig_dir, "cogng_demog_slope.pdf"), width = 6, height = 5)

# Figure 5 Sex
# Stacked + percent
summary.count = meta.data %>%filter(!is.na(sex))%>%group_by(reference, sex)%>%tally()
figure5 = bar_pct(fill_var = "sex", y_var = "n",x_var = "reference", summary.count,sex_pal)
ggsave(plot = figure5,filename = file.path(fig_dir, "pct_sex.pdf"), width = 3, height = 5)

# Figure 6 PathoAD
# Stacked + percent
summary.pathoAD = meta.data %>% filter(!is.na(pathoAD))%>%
  group_by(reference,  pathoAD) %>%tally()%>%mutate(pct.pathoAD = n*100/sum(n))%>%
  mutate(pathoAD = case_when(pathoAD==0 ~ "non-AD", TRUE ~ "AD"))
figure6 = bar_pct(fill_var = "pathoAD", y_var = "n",x_var = "reference", summary.pathoAD,ad_pal,method="chisq",lgd.pos = "top")
ggsave(plot = figure6,filename = file.path(fig_dir, "pct_pathoAD.pdf"), width = 5, height = 5)

# Figure 7 Demographics
summary.pop = meta.data %>% filter(!is.na(population))%>%
  group_by(reference,  population) %>%tally()%>%mutate(pct = n*100/sum(n))
figure7 = bar_pct(fill_var = "population", y_var = "n",x_var = "reference", summary.pop,population_pal,method="chisq",lgd.pos = "right")
ggsave(plot = figure7,filename = file.path(fig_dir, "pct_population.pdf"), width = 5, height = 5)


panel1 = ggarrange(plotlist = list(figure1, figure2, figure3, figure4), nrow = 1, ncol = 4, common.legend = T)
panel2 = cowplot::plot_grid(figure5, figure6, figure7, NULL, nrow = 1, ncol = 4, rel_widths= c(1,1,2,0.5))

ggsave(plot = ggarrange(panel1, panel2, nrow = 2), filename = file.path(fig_dir, 'data_description_plots.pdf'), width = 10, height = 10)



###############################################################################|

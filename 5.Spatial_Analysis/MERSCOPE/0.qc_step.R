#!/usr/bin/env Rscript

rm(list = ls(all.names = TRUE))
gc()

library(argparse)

parser <- ArgumentParser(description = 'MERSCOPE Quality Control Step')
parser$add_argument('--region_dir', required = TRUE, help = 'Path to region directory')
parser$add_argument('--region', required = TRUE, help = 'Region name')
parser$add_argument('--sample_name', required = TRUE, help = 'Sample name')
parser$add_argument('--output_dir', required = TRUE, help = 'Output directory')
parser$add_argument('--renv_path', required = TRUE, help = 'Path to renv library')
parser$add_argument('--utils_path', required = TRUE, help = 'Path to utils_merscope.R')
parser$add_argument('--metadata_file', required = TRUE, help = 'Path to sample metadata CSV')

args <- parser$parse_args()


library(renv)
renv::load(args$renv_path)
library(tidyverse)
library(Seurat)
library(progressr)
library(logger)
library(scCustomize)

source(args$utils_path)


data_raw_dir <- file.path(args$output_dir, "data_raw", args$sample_name)
data_qc_dir <- file.path(args$output_dir, "data_qc", args$sample_name)
figure_dir <- file.path(args$output_dir, "figures", args$sample_name)
output_dir <- file.path(args$output_dir, "output", args$sample_name)

dir.create(data_raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

log_info("Starting QC for sample: {args$sample_name}, region: {args$region}")

log_info("Loading MERSCOPE data from: {args$region_dir}")
sobj <- load_merscope(
  folder_path = args$region_dir,
  fov_name = "fov1",
  filter = "^Blank-"  # exclude blank control probes
)

# Add sample metadata ====
sample.meta <- read_csv(args$metadata_file)
selected_sample.meta = filter(sample.meta,`Experiment name`== args$sample_name & str_detect( args$region,Region))
sobj@meta.data$orig.ident = args$sample_name
sobj@meta.data$region = args$region
sobj@meta.data = cbind(sobj@meta.data, selected_sample.meta%>%select(-c("Region","Combo Name")))
sobj@meta.data$brain_region = "Human-Brain-MFC"
sobj@meta.data$sample_id = paste0(sobj@meta.data$`Short Sample ID`%>%unique, "-", args$region )

raw_file <- file.path(data_raw_dir, sprintf("raw_%s.rds", args$region))
saveRDS(sobj, raw_file)
log_info("Saved raw object to: {raw_file}")

log_info("Performing quality control analysis")
output <- qc_plot(
  sobj@meta.data, 
  save.as = file.path(figure_dir, paste0(args$region, "_qc_hist.pdf")),
  q_umi = c(0.05, 0.99), 
  min_umi = 0, 
  max_umi = 4000,
  q_area = c(0.05, 1), 
  min_area = 0, 
  max_area = 5000,
  q_feat = c(0.05, 0.99), 
  min_feat = 30, 
  max_feat = Inf
)

save.var <- c("area_thres", "umi_thres", "feat_thres", "stats", "meta.data", "fig")
qc_results_file <- file.path(output_dir, sprintf("qc_results_%s.rds", args$region))
saveRDS(output[save.var], qc_results_file)
log_info("Saved QC results to: {qc_results_file}")

sub <- subset(output$sobj, subset = pass_qc == TRUE)

qc_file <- file.path(data_qc_dir, sprintf("qc_%s.rds", args$region))
saveRDS(sub, qc_file)
log_info("Saved QC-filtered object to: {qc_file}")

log_info("Generating FOV plots")
coord <- sobj@images[["fov1"]]@boundaries[["centroids"]]@coords
xlim.fov <- range(coord[, "x"])
ylim.fov <- range(coord[, "y"])

p1 <- qc_fov_feature_contin(
  data = output$sobj, 
  feature = "nCount_Vizgen", 
  fov = "fov1", 
  limit = c(0, 3000), 
  xlim = xlim.fov, 
  ylim = ylim.fov
) + ggtitle("Before QC")

p2 <- qc_fov_feature_contin(
  data = sub, 
  feature = "nCount_Vizgen", 
  fov = "fov1", 
  limit = c(0, 3000), 
  xlim = xlim.fov, 
  ylim = ylim.fov
) + ggtitle("After QC")

cols <- c("#0474BA", "orange")
p3 <- qc_fov_feature_discrete(
  data = output$sobj, 
  feature = "pass_qc", 
  fov = "fov1", 
  cols = cols, 
  xlim = xlim.fov, 
  ylim = ylim.fov
) + ggtitle("Pass QC")

fov_plot_file <- file.path(figure_dir, paste0(args$region, "_qc_fov.pdf"))
ggsave(p1 + p2 + p3, filename = fov_plot_file, width = 11.5, height = 3.3)
log_info("Saved FOV plots to: {fov_plot_file}")

log_info("QC step completed successfully for sample: {args$sample_name}, region: {args$region}")


rm(list =ls())
gc()
setwd("Code/MERSCOPE_Analysis/SpatialLagAnalysis")
source('utils_spatial_lag.R')
source('utils_results_summary.R')
library(tidyverse)
##################################################################|
# SENSITIVITY ANALYSIS CONFIGURATION
##################################################################|


# Define spatial neighborhood size (k-nearest neighbors) to test
spatial_configs <- data.frame(
  config_id = c("k3", "k6", "k9", "k12"),
  k = c(3, 6, 9, 12), #this range capture both contact and paracrine signaling
  stringsAsFactors = FALSE
)



##################################################################|
# STATIC PARAMETERS
##################################################################|
data_path <- "/path/to/wm_gm_annotation_dir"
data_subpath <-"f_results_npc=100/niche_annotation"
layer_oi <- "Grey Matter"


# Predictor (ligand-receptors) parameters
normalize <- "log1p_scale"    # log1p ligand/receptor raw counts then z-score per gene
input_cells <- "ast_mic"      # model neighborhood includes astrocytes + microglia
cell_types <- c("Astrocyte", "Microglia") #model neighborhood includes astrocytes + microglia only

# Response variable (Ast10 score) parameters
score_method <- "mean"        # aggregate Ast10 marker expression as mean
response_normalize <- "NULL"       # No normalization of Ast10 marker expression between averging
response_type <- "continuous" # Use continuousAst10 score as response
ast_markers <- c("SLC38A2", "SMTN", "SLC39A11", "CACNA1D", "MT1E") # Ast10 markers

# Mixed Effect Model parameters
model_per_ligand <- TRUE      # separate model per ligand (not pooled)
use_lr_pairs <- TRUE          # test ligand-receptor pairs, not ligands alone
family_type <- "gaussian"     # linear mixed model for continuous response
random_effect <- "sample_id"  # random intercept per donor
zero_inflated <- FALSE        # whether use zero-inflated distribution in the glmmTMB
lag_covariates <- "nCount_Vizgen"  # control for cell density in spatial lag
covariates <- NULL                 # no additional fixed-effect covariates
add_lr_interaction <- FALSE   # no ligand × receptor interaction term

# Output and Logging
main_out_dir <- file.path("..",data_subpath, "final_results_seurat_v5.1.0_sel_>=0.3")
dir.create(main_out_dir, recursive = TRUE, showWarnings = FALSE)


log_file <- file.path(main_out_dir, "sensitivity_analysis.log")
log_appender(appender_tee(log_file))

log_info(paste(rep("=", 70), collapse = ""))
log_info("SPATIAL SENSITIVITY ANALYSIS")
log_info(paste(rep("=", 70), collapse = ""))


write.csv(spatial_configs,
          file.path(main_out_dir, "spatial_configs.csv"),
          row.names = FALSE)

##################################################################|
# LOAD DATA 
##################################################################|

log_info("\n### Loading multi-sample MERSCOPE data ###")

sobj_list <- LoadMerscope_MultiSample(
  obj_path = file.path(data_path,data_subpath, "integrated_list_object_with_niches.rds"),
  celltype = NULL,
  samples = NULL,
  id_name = "sample_id",
  layer_oi = layer_oi,
  niche_var = "predicted.niche",
  min_cells_per_sample = 
)

log_info("Successfully loaded {length(sobj_list)} samples")

### Load Ligand Information ###
log_info("Loading ligand VIP data...")
ligand_vip <- load_ligand_vip()%>%filter(selection_rate >=0.30)

# Get all available PLSR ligands
all_plsr_ligands <- ligand_vip %>% 
  pull(ligand) %>% 
  intersect(rownames(sobj_list[[1]]))

log_info("Found {length(all_plsr_ligands)} PLSR ligands in data")



##################################################################|
# SCORE CELLS 
##################################################################|

log_info("\n### Calculating astrocyte subset scores ###")

response_col <- "ast10_score"
response_threshold <- 0

sobj_list_scored <- GetCellSubset_MultiSample(
  sobj_list = sobj_list,
  markers = ast_markers,
  slot = "count",
  method = score_method,
  threshold = response_threshold,
  response_type = response_type,
  validate_markers = TRUE,
  response_name = response_col,
  normalize = response_normalize,
  output_format = "seurat"
)

# Filter to astrocytes and microglia
log_info("Filtering to astrocytes and microglia...")

input_sobj_list <- list()

for (i in seq_along(sobj_list_scored)) {
  cells_to_keep <- colnames(sobj_list_scored[[i]])[
    sobj_list_scored[[i]]$de.novo.celltype %in% cell_types
  ]
  if (length(cells_to_keep) > 0) {
    input_sobj_list[[i]] <- sobj_list_scored[[i]][, cells_to_keep]
  }
}
names(input_sobj_list) <- names(sobj_list_scored)
input_sobj_list <- input_sobj_list[lengths(input_sobj_list) > 0]

log_info("Retained {length(input_sobj_list)} samples")

# Get astrocyte cells for response
get_celltype <- function(sobj_list, celltype = "Astrocyte") {
  unlist(lapply(sobj_list, function(x) {
    colnames(x)[x$de.novo.celltype == celltype]
  }))
}

astrocyte_cells <- get_celltype(input_sobj_list, "Astrocyte")
log_info("Found {length(astrocyte_cells)} astrocyte cells for response")

##################################################################|
# GET ALL L-R PAIRS TO TEST
##################################################################|

log_info("\n### Determining L-R pairs to test ###")

# Get all L-R pairs for all PLSR ligands
all_lr_pairs <- get_ligand_receptors(all_plsr_ligands)
all_lr_pairs <- all_lr_pairs %>% 
  filter(to %in% rownames(input_sobj_list[[1]])) %>%
  rename(ligand = from, receptor = to)

log_info("Testing {nrow(all_lr_pairs)} unique L-R pairs")

##################################################################|
# SENSITIVITY ANALYSIS LOOP - ALL L-R PAIRS AT ALL CONFIGS
##################################################################|

log_info("\n### Starting sensitivity analysis ###")
log_info("Testing all {nrow(all_lr_pairs)} L-R pairs at all {nrow(spatial_configs)} configurations\n")

all_results <- list()
all_summaries <- list()

# Get unique ligands and receptors
ligands_for_model <- unique(all_lr_pairs$ligand)
receptors_for_model <- unique(all_lr_pairs$receptor)

for (i in 1:nrow(spatial_configs)){
  
  config <- spatial_configs[i, ]
  config_id <- config$config_id
  
  log_info(paste(rep("=", 70), collapse = ""))
  log_info("CONFIG {i}/{nrow(spatial_configs)}: k={config$k}")
  log_info("  Testing all {length(ligands_for_model)} ligands × {length(receptors_for_model)} receptors")
  log_info(paste(rep("=", 70), collapse = ""))
  
  # Create config directory
  config_dir <- file.path(main_out_dir, config_id)
  dir.create(config_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Fit models for ALL L-R pairs at this configuration
  tryCatch({
    
    log_info("  Fitting spatial lag models...")
    
    results <- SpatialLagX_MultiSample_Model(
      sobj_list = input_sobj_list,
      ligands_oi = ligands_for_model,
      receptors_oi = receptors_for_model,
      assay = "Vizgen",
      slot = "count",
      k = config$k,
      normalize = normalize,
      model_per_ligand = model_per_ligand,
      add_lr_interaction = add_lr_interaction,
      selected_cells = astrocyte_cells,
      response_col = response_col,
      response_type = response_type,
      response_threshold = response_threshold,
      validate_inputs = TRUE,
      random_effect_structure = random_effect,
      covariates = covariates,
      lag_covariates = lag_covariates,
      model_engine = "glmmTMB",
      family_type = family_type,
      zero_inflated = zero_inflated
    )
    

    log_success("  Fitted {length(results)} L-R pair models")
    
    # Extract summaries - p_adj calculated wihtinthis config only
    log_info("  Extracting summaries (p_adj within config)...")
    summary_result <- extract_multisample_summary(
      results, 
      combine_effects = FALSE,
      include_random_effects = FALSE,
      include_glance = FALSE,
      covariates = c(covariates, lag_covariates),
      rm.covariates = TRUE,
      adjust_group.by = "predictor_type",  # p.adj within predictor type
      adjust_method = "BH"
    )
    
    # Get fixed effects dataframe
    summary_df <- summary_result[["fixed_effects"]]
    
    if (is.null(summary_df) || nrow(summary_df) == 0) {
      log_warn("  No summary data for {config_id}")
      next
    }
    
    # Add configuration metadata
    summary_df <- summary_df %>%
      mutate(
        config_id = config_id,
        k = config$k,
        n_samples = length(input_sobj_list),
        n_astrocytes = length(astrocyte_cells),
        model_type = "separate_lr_pairs",
        input_cells = input_cells,
        family_type = family_type
      )
    
    # Extract ligand/receptor from model_name if not present
    if (!"ligand" %in% colnames(summary_df)) {
      summary_df <- summary_df %>%
        mutate(
          ligand = str_extract(model_name, "^[^_]+"),
          receptor = str_extract(model_name, "(?<=_)[^_]+$")
        )
    }
    
    # Add expected configuration and mechanism info
    summary_df <- summary_df %>%
      mutate(
        # Create L-R pair identifier
        lr_pair = paste0(ligand, "_", receptor)
      )
    
    # Store results
    all_results[[config_id]] <- results
    all_summaries[[config_id]] <- summary_df
    
    # Save config outputs
    saveRDS(results, file.path(config_dir, "results.rds"))
    write.csv(summary_df, file.path(config_dir, "summary.csv"), row.names = FALSE)
    
    # Generate forest plot
    if (nrow(summary_df) > 0) {
      log_info("  Generating forest plot...")
      tryCatch({
        create_comparison_forests(
          successful_summaries = summary_df %>% filter(effect_type == "fixed"),
          main_output_dir = config_dir,
          vip_scores = ligand_vip 
        )
      }, error = function(e) {
        log_warn("  Forest plot failed: {e$message}")
      })
    }
    
    # Summary statistics
    n_sig_ligand <- sum(summary_df$p_adj < 0.05 & summary_df$predictor_type == "ligand", na.rm = TRUE)
    n_total_ligand <- sum(summary_df$predictor_type == "ligand", na.rm = TRUE)
    log_info("  Significant ligand terms: {n_sig_ligand}/{n_total_ligand} ({round(100*n_sig_ligand/n_total_ligand,1)}%)")
    

    
  }, error = function(e) {
    log_error("  Error in {config_id}: {e$message}")
    print(e)
  })
  
  log_info("")
}


# Combine all summaries
combined_summary <- bind_rows(all_summaries)
write.csv(combined_summary, 
          file.path(main_out_dir, "all_configs_combined.csv"), 
          row.names = FALSE)

##################################################################|
# SUMMARIZE L-R PAIR PERFORMANCE ACROSS CONFIGS
##################################################################|

log_info("\n### Analyzing L-R pair significance at expected vs other scales ###")


lr_pair_performance <- combined_summary %>%
  filter(effect_type == "fixed", predictor_type == "ligand") %>%
  filter(!is.na(ligand), !is.na(receptor)) %>%
  # Keep only main lag terms (not covariates)
  filter(str_detect(term, paste0("^`?lag_", ligand))) %>%
  filter(!str_detect(term, ":")) %>%
  group_by(lr_pair, ligand, receptor, config_id, k) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(is_sig = p_adj < 0.05) %>%
  group_by(lr_pair, ligand, receptor) %>%
  summarise(
    n_configs_tested = n_distinct(config_id),
    n_configs_sig = sum(is_sig, na.rm = TRUE),
    configs_sig = paste(sort(unique(config_id[is_sig])), collapse = ", "),
    mean_effect_size = mean(estimate, na.rm = TRUE),
    sd_effect_size = sd(estimate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_configs_sig), desc(abs(mean_effect_size)))

write.csv(lr_pair_performance,
          file.path(main_out_dir, "lr_pair_performance.csv"),
          row.names = FALSE)

log_info("\n=== L-R Pair Performance Summary ===")
log_info("Total L-R pairs tested: {nrow(lr_pair_performance)}")

##################################################################|
# FINAL SUMMARY
##################################################################|

log_info("\n" %+% paste(rep("=", 70), collapse = ""))
log_info("SENSITIVITY ANALYSIS COMPLETE")
log_info(paste(rep("=", 70), collapse = ""))

# Save complete session
saveRDS(
  list(
    all_results = all_results,
    all_lr_pairs = all_lr_pairs,
    all_summaries = all_summaries,
    combined_summary = combined_summary,
    lr_pair_performance = lr_pair_performance,
    spatial_configs = spatial_configs,
    session_info = sessionInfo()
  ),
  file.path(main_out_dir, "sensitivity_analysis_complete.rds")
)

log_success("\nAll analyses complete!")
log_info("Results: {main_out_dir}")
saveRDS(sobj_list_scored, file.path(main_out_dir,"sobj_list_scored.rds"))

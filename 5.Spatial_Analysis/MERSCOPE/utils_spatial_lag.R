library(Seurat)
library(ggpubr)
library(spdep)
library(readr)
library(readxl)
library(stringr)
library(logger)
library(purrr)
library(dplyr)
library(tibble)
library(lme4)      
library(glmmTMB)   
library(broom.mixed) 


##################################################################|
# Data Loading Functions
##################################################################|

# Load ligand VIP data with error handling

load_ligand_vip <- function(vip_path = "/path/to/ligand_vip_detailed_results.rds"){
  tryCatch({
    ligand_vip <- readRDS(vip_path) %>% .[["ligand_stability"]]
    ligand_vip = mutate(ligand_vip, vip = median_vip_signed, ligand = ligand_clean) # define column that match format
    log_info("Ligand VIP data loaded successfully: {nrow(ligand_vip)} entries")
    return(ligand_vip)
  }, error = function(e) {
    log_error("Failed to load ligand VIP data: {e$message}")
    return(NULL)
  })
}

# Enhanced ligand-receptor lookup with caching
get_ligand_receptors <- function(ligands_oi = NULL, 
                                 receptors_oi = NULL,
                                 lr_network_path = "/path/to/lr_network_human.rds",
                                 cache_network = TRUE) {
  
  # Load or retrieve cached network
  if (cache_network && exists(".lr_network_cache", envir = .GlobalEnv)) {
    lr_network <- get(".lr_network_cache", envir = .GlobalEnv)
  } else {
    tryCatch({
      lr_network <- readRDS(lr_network_path)
      if (cache_network) {
        assign(".lr_network_cache", lr_network, envir = .GlobalEnv)
      }
      log_info("Loaded L-R network: {nrow(lr_network)} interactions")
    }, error = function(e) {
      log_error("Failed to load L-R network: {e$message}")
      return(NULL)
    })
  }
  
  # Filter based on input
  if (!is.null(ligands_oi) && is.null(receptors_oi)) {
    return(filter(lr_network, from %in% ligands_oi))
  } else if (is.null(ligands_oi) && !is.null(receptors_oi)) {
    return(filter(lr_network, to %in% receptors_oi))
  } else if (!is.null(ligands_oi) && !is.null(receptors_oi)) {
    return(filter(lr_network, from %in% ligands_oi & to %in% receptors_oi))
  } else {
    return(lr_network)
  }
}

##################################################################|
# MERSCOPE Multi-Sample Loading Function
##################################################################|

LoadMerscope_MultiSample <- function(obj_path = "/path/to/annot_obj_list.rds",
                                     celltype = "Astrocyte", 
                                     samples = NULL, #List of sample ID
                                     id_name = "sample_id", #Meta data column to name sample by
                                     layer_oi = NULL,
                                     niche_var = 'predicted.niche',
                                     min_cells_per_sample = 50) {
  
  tryCatch({
    annot_obj_list <- readRDS(obj_path)
    all_samples <- lapply(annot_obj_list, function(x){x@meta.data[[id_name]]%>%unique})%>%do.call("c",.)
    names(annot_obj_list) = all_samples
    
    if (is.null(samples)) {
      samples = all_samples
    }
    
    
    
    # Process each sample
    processed_samples <- list()
    
    for (i in seq_along(samples)) {
      sample_idx <- samples[i]
      sobj <- annot_obj_list[[sample_idx]]
      log_info("Processing sample {sample_idx}...")
      
      
      # Cell type subsetting
      if (!is.null(celltype)) {
        Idents(sobj) <- "de.novo.celltype"
        available_types <- unique(Idents(sobj))
        
        if (!celltype %in% available_types) {
          log_warn("Cell type '{celltype}' not found in sample {sample_idx}. Skipping...")
          next
        }
        
        sobj <- subset(sobj, idents = celltype)
        log_info("Sample {sample_idx}: Subset to {celltype}: {ncol(sobj)} cells")
      }
      
      # Layer/niche filtering
      if (!is.null(layer_oi)) {
        Idents(sobj) <- niche_var
        available_layers <- unique(Idents(sobj))
        
        if (!layer_oi %in% available_layers) {
          log_warn("Layer '{layer_oi}' not found in sample {sample_idx}. Skipping...")
          next
        }
        
        sobj <- subset(sobj, idents = layer_oi)
        log_info("Sample {sample_idx}: Subset to {layer_oi}: {ncol(sobj)} cells")
      }
      
      # Check minimum cell count
      if (ncol(sobj) < min_cells_per_sample) {
        log_warn("Sample {sample_idx} has only {ncol(sobj)} cells (< {min_cells_per_sample}). Skipping...")
        next
      }
      
      # Add sample identifier to metadata
      
      processed_samples[[sample_idx]] <- sobj
    }
    
    if (length(processed_samples) == 0) {
      stop("No valid samples found after filtering")
    }
    
    log_info("Successfully processed {length(processed_samples)} samples")
    return(processed_samples)
    
  }, error = function(e) {
    log_error("Error loading MERSCOPE multi-sample data: {e$message}")
    return(NULL)
  })
}

##################################################################|
# Cell Subset Scoring Function
##################################################################|

GetCellSubset_MultiSample <- function(sobj_list, 
                                      markers = c("SLC38A2","SMTN","SLC39A11","CACNA1D"),
                                      slot = "count",
                                      celltype = "Astrocyte",
                                      method = c("multiplicative", "additive", "mean"),
                                      threshold = 0,
                                      response_type = c("binary", "continuous"),
                                      validate_markers = TRUE,
                                      response_name = "ast10_score",
                                      normalize_method = NULL,
                                      output_format = c("seurat","data.frame")
) {
  
  method <- match.arg(method)
  response_type <- match.arg(response_type)
  if (is.null(normalize_method)|normalize_method=="NULL"){
    log_info("Scoring markers (no normalization)")
    }else{
    log_info("Scoring markers (normalize: {normalize_method})")
    }
  tryCatch({
    if (output_format == "data.frame") combined_results <- list()
    
    
    for (sample_id in names(sobj_list)) {
      sobj <- sobj_list[[sample_id]]
      
      # Cell type subsetting
      if (!is.null(celltype)) {
        Idents(sobj) <- "de.novo.celltype"
        available_types <- unique(Idents(sobj))
        
        if (!celltype %in% available_types) {
          log_warn("Cell type '{celltype}' not found in {sample_id}. Skipping...")
          next
        }
        
        sobj <- subset(sobj, idents = celltype)
        log_info("{sample_id}: Subset to {celltype}: {ncol(sobj)} cells")
      }
      
      # Validate markers exist in the object
      if (validate_markers) {
        available_features <- rownames(sobj)
        missing_markers <- setdiff(markers, available_features)
        if (length(missing_markers) > 0) {
          log_warn("{sample_id}: Missing markers: {paste(missing_markers, collapse=', ')}")
          current_markers <- intersect(markers, available_features)
          if (length(current_markers) == 0) {
            log_warn("{sample_id}: No valid markers found. Skipping...")
            next
          }
        } else {
          current_markers <- markers
        }
      } else {
        current_markers <- markers
      }
      
      # Calculate scores
      df <- FetchData(sobj, assay = "Vizgen", current_markers, layer = slot)
      if (!is.null(normalize_method) & !normalize_method=="NULL"){
        df  = .normalize_expression(df, normalize_method = normalize_method)
      }
      
      
      # Calculate composite score based on method
      if (method == "multiplicative") {
        df$score <- apply(df[, current_markers, drop = FALSE], 1, prod)
      } else if (method == "additive") {
        df$score <- rowSums(df[, current_markers, drop = FALSE])
      } else if (method == "mean") {
        df$score <- rowMeans(df[, current_markers, drop = FALSE])
      }
      
      # Add region information
      df$sample_id <- sobj@meta.data$sample_id
      df$experiment <- sobj@meta.data$`Experiment name`
      df$region <- sobj@meta.data$region
      df$chemistry <- sobj@meta.data$Chemistry
      df$cell_names <- rownames(df)
      
      # Determine response variable
      if (response_type == "binary") {
        if (is.null(threshold)){threshold=median(df$score)}
        df$response <- as.numeric(df$score > threshold)
        cells_oi <- df %>% filter(response == 1) %>% pull("cell_names")
        log_info("{sample_id}: {length(cells_oi)} cells above threshold ({threshold})")
      } else {
        df$response <- df$score
        cells_oi <- rownames(df)
      }
      
      if (output_format == "data.frame") {
        combined_results[[sample_id]] <- list(
          df = df,
          cells_oi = cells_oi,
          method = method,
          threshold = threshold,
          response_type = response_type,
          markers_used = current_markers
        )
      }
      
      if (output_format == "seurat") {
        
        tmp = df%>%select( response)%>%`colnames<-`(response_name)
        sobj_list[[sample_id]] <- AddMetaData(sobj_list[[sample_id]], tmp)
        
        # check number of cells
        if (sum(!is.na(sobj_list[[sample_id]][[response_name]]))==dim(sobj)[2]){
          log_success("Added {response_name} to seurat object for {sample_id}.")
        }else{
          log_error("Some cells within {celltype} did not have a response score.")
        }
      }
    }
    
    
    
    if (output_format == "data.frame"){
      # Combine all data frames
      combined_df <- do.call(rbind, lapply(combined_results, function(x) x$df))
      all_cells_oi <- unlist(lapply(combined_results, function(x) x$cells_oi))
      
      log_info("Multi-sample cell subset scoring completed: {nrow(combined_df)} total cells across {length(combined_results)} samples")
      
      return(list(
        combined_df = combined_df,
        sample_results = combined_results,
        all_cells_oi = all_cells_oi,
        method = method,
        threshold = threshold,
        response_type = response_type,
        markers_used = markers
      ))
    }
    if (output_format == "seurat"){
      total_cells = lapply(sobj_list, function(m){
        if (!is.null(celltype)) {
          m@meta.data%>%filter(de.novo.celltype == celltype)%>%nrow
        }else{
          m@meta.data%>%nrow
        }
      })%>%unlist%>%sum
      log_info("Multi-sample cell subset scoring completed: {total_cells} total cells across {length(sobj_list)} samples")
      
      return(sobj_list)
    }
    
    
  }, error = function(e) {
    log_error("Error in GetCellSubset_MultiSample: {e$message}")
    return(NULL)
  })
}

##################################################################|
# Spatial Lag Mixed Effects Model
##################################################################|

SpatialLagX_MultiSample_Model <- function(
    sobj_list,               # List of Seurat objects by sample
    ligands_oi,              # Vector of ligand genes of interest
    receptors_oi = NULL,     # Optional vector of receptor genes of interest
    assay = "Vizgen",        # Assay slot name in Seurat object
    slot = "count",          # Which data slot to use
    k = 6,                   # Number of neighbors for KNN
    normalize = c("raw", "log1p", "log1p_scale"), # Normalization method
    model_per_ligand = TRUE, # Whether to run separate model per ligand-receptor pair
    add_lr_interaction = TRUE, # Whether to add ligand-receptor interactions
    selected_cells = NULL,   # Subset of cells to use for modeling
    response_col = "response", # Column name for response variable
    response_type = c("binary", "continuous"), # Type of response variable
    response_threshold = 0,  # Threshold for binary response (optional, for binar)
    validate_inputs = TRUE,  # Whether to validate gene names
    random_effect_structure = c("sample_id"), # Random effect structure
    covariates = NULL, # other covariates to include in mixed model, column names from sobj@meta.data
    lag_covariates = c("nCount_Vizgen"),
    model_engine = c("glmmTMB", "lme4"), # Mixed model engine
    family_type = NULL,      # Will be auto-determined based on response_type
    zero_inflated = T, #whether use zero-inflated distribution in the glmmTMB
    control_params = NULL    # Additional control parameters
) { 
  
  normalize <- match.arg(normalize)
  response_type <- match.arg(response_type)
  random_effect_structure <- match.arg(random_effect_structure)
  model_engine <- match.arg(model_engine)
  
  # Auto-determine family if not specified
  if (is.null(family_type)) {
    family_type <- if (response_type == "binary") "binomial" else "gaussian"
  }
  
  tryCatch({
    log_info("Fitting spatial lag models across {length(sobj_list)} samples...")
    log_info("  Family: {family_type} | Random: (1|{random_effect_structure}) | Covariates: {paste(c(covariates, paste0('lag_', lag_covariates)), collapse=', ')}")

    # Input validation across all samples
    if (validate_inputs) {
      all_features <- unique(unlist(lapply(sobj_list, rownames)))
      
      # Check ligands
      missing_ligands <- setdiff(ligands_oi, all_features)
      if (length(missing_ligands) > 0) {
        log_warn("Missing ligands across all samples: {paste(missing_ligands, collapse=', ')}")
        ligands_oi <- intersect(ligands_oi, all_features)
      }
      
      # Check receptors
      if (!is.null(receptors_oi)) {
        missing_receptors <- setdiff(receptors_oi, all_features)
        if (length(missing_receptors) > 0) {
          log_warn("Missing receptors across all samples: {paste(missing_receptors, collapse=', ')}")
          receptors_oi <- intersect(receptors_oi, all_features)
        }
      }
      
      if (length(ligands_oi) == 0) {
        stop("No valid ligands found across all samples")
      }
    }
    
    # Process each sample and combine data
    combined_data <- .process_multisample_spatial_data(
      sobj_list = sobj_list,
      ligands_oi = ligands_oi,
      receptors_oi = receptors_oi,
      assay = assay,
      slot = slot,
      k = k,
      normalize = normalize,
      response_col = response_col,
      selected_cells = selected_cells,
      covariates = covariates
    )

    if (is.null(combined_data)) {
      stop("Failed to process multi-sample spatial data")
    }
    
    # Run models
    if (model_per_ligand) {
      
      
      # Create input dataframe for modeling
      if (!is.null(receptors_oi)) {
        log_info("Running individual mixed-effects models for {length(ligands_oi)} ligand-receptor pairs")
        input_df <-  get_ligand_receptors(ligands_oi)%>%select(from, to)%>%filter(to %in% receptors_oi)
        colnames(input_df) = c("ligand", "receptor")
        
        results = purrr::pmap(input_df, function(ligand, receptor){
          
          .fit_mixed_effects_spatial_model(
            combined_data = combined_data,
            .ligand = ligand,
            .receptor = receptor,
            add_lr_interaction = add_lr_interaction,
            random_effect_structure = random_effect_structure,
            covariates = covariates,
            lag_covariates = lag_covariates,
            model_engine = model_engine,
            zero_inflated = zero_inflated,
            family_type = family_type,
            control_params = control_params
          )
          
        })%>%setNames(paste0(input_df$ligand, ifelse(is.na(input_df$receptor), "", paste0("_", input_df$receptor))))
      } else {
        log_info("Running individual mixed-effects models for {length(ligands_oi)} ligands")
        input_df <- data.frame(ligand = ligands_oi)
        
        results = pmap(input_df, function(ligand){
          
          .fit_mixed_effects_spatial_model(
            combined_data = combined_data,
            .ligand = ligand,
            .receptor = NULL,
            add_lr_interaction = add_lr_interaction,
            random_effect_structure = random_effect_structure,
            covariates = covariates,
            lag_covariates = lag_covariates,
            model_engine = model_engine,
            zero_inflated = zero_inflated,
            family_type = family_type,
            control_params = control_params
          )
          
        })%>%setNames(paste0(input_df$ligand, ifelse(is.na(input_df$receptor), "", paste0("_", input_df$receptor))))
      }
      
    } else {
      log_info("Running joint mixed-effects model across all ligands")
      
      results <- .fit_mixed_effects_spatial_model(
        combined_data = combined_data,
        .ligand = ligands_oi,  # Pass all ligands
        .receptor = receptors_oi,
        add_lr_interaction = add_lr_interaction,
        random_effect_structure = random_effect_structure,
        covariates = covariates,
        lag_covariates = lag_covariates,
        model_engine = model_engine,
        family_type = family_type,
        zero_inflated = zero_inflated,
        control_params = control_params
      )
    }
    
    return(results)
    
  }, error = function(e) {
    log_error("Error in SpatialLagX_MultiSample_Model: {e$message}")
    return(NULL)
  })
}

##################################################################|
# Core Multi-Sample Data Processing Function
##################################################################|

# Helper function to apply chosen normalization

.normalize_expression <- function(expression_df, normalize_method, total_counts = NULL) {
  if (is.null(expression_df)) {
    return(NULL)
  }
  
  # Ensure expression_df is a data frame for dplyr operations
  expression_df <- as.data.frame(expression_df)
  
  # Apply normalization based on the method
  if (normalize_method == "log1p") {
    normalized_df <- expression_df %>%
      mutate(across(everything(), ~ log1p(.x)))

  } else if (normalize_method == "log1p_scale") {
    normalized_df <- expression_df %>%
      mutate(across(everything(), ~ log1p(.x))) %>%
      mutate(across(everything(), ~ as.numeric(scale(.x))))

  } else if (normalize_method == "raw") {
    # No normalization
    normalized_df <- expression_df
    
  } else {
    warning("Unknown normalization method: ", normalize_method, ". No normalization applied.")
    normalized_df <- expression_df # Return original if method is unknown
  }
  
  rownames(normalized_df) <- rownames(expression_df)
  return(normalized_df)
}


.process_multisample_spatial_data <- function(
    sobj_list,
    ligands_oi,
    receptors_oi = NULL,
    assay = "Vizgen",
    slot = "count",
    k = 6,
    normalize = "log1p",
    response_col = "response",
    selected_cells = NULL,
    covariates = NULL
) {
  
  tryCatch({
    
    combined_df_list = lapply(names(sobj_list), 
                              function(sample_id){
                                # Process each sample
                                sobj <- sobj_list[[sample_id]]
                                
                                # Extract coordinates
                                coords <- GetTissueCoordinates(sobj, scale = "tissue") %>% 
                                  column_to_rownames("cell")
                                
                                # Validate response variable
                                if (!response_col %in% colnames(sobj@meta.data)) {
                                  log_warn("Response column '{response_col}' not found in {sample_id} metadata. Skipping...")
                                  next
                                }
                                
                                # Extract response variable
                                cell_state <- sobj@meta.data[[response_col]]
                                names(cell_state) <- rownames(sobj@meta.data)
                                
                                # Filter ligands available in this sample
                                available_features <- rownames(sobj)
                                current_ligands <- intersect(ligands_oi, available_features)
                                
                                if (length(current_ligands) == 0) {
                                  log_warn("No ligands available in {sample_id}. Skipping...")
                                  next
                                }
                                
                                # Fetch expression data
                                lig_expr <- FetchData(sobj, vars = current_ligands, assay = assay, layer = slot)
                                
                                # Fetch receptor data if specified
                                rec_expr <- NULL
                                if (!is.null(receptors_oi)) {
                                  current_receptors <- intersect(receptors_oi, available_features)
                                  if (length(current_receptors) > 0) {
                                    rec_expr <- FetchData(sobj, vars = current_receptors, assay = assay, layer = slot)
                                  }
                                }
                                
                                # Normalize expression data
                                total_counts = sobj$nCount_Vizgen
                                
                                # Normalize ligand expression
                                lig_expr_norm <- .normalize_expression(
                                  expression_df = lig_expr,
                                  normalize_method = normalize,
                                  total_counts = total_counts
                                )%>%as.matrix()
                                # Normalize receptor data if available
                                if (!is.null(rec_expr)) {
                                  rec_expr_norm <- .normalize_expression(
                                    expression_df = rec_expr,
                                    normalize_method = normalize,
                                    total_counts = total_counts
                                  )%>%as.matrix
                                } else {
                                  rec_expr_norm <- NULL
                                }
                                
                                # Create spatial weights
                                spatial_weights <- .create_spatial_weights(coords, k)
                                
                                if (is.null(spatial_weights)) {
                                  log_warn("Failed to create spatial weights for {sample_id}. Skipping...")
                                  next
                                }
                                
                                # Compute spatial lags for ligands ------
                                lag_matrix <- matrix(NA, nrow = nrow(lig_expr_norm), ncol = ncol(lig_expr_norm))
                                colnames(lag_matrix) <- paste0("lag_", colnames(lig_expr_norm))
                                rownames(lag_matrix) <- rownames(lig_expr_norm)
                                for (i in 1:ncol(lig_expr_norm)) {
                                  lag_matrix[, i] <- lag.listw(spatial_weights, lig_expr_norm[, i], zero.policy = TRUE)
                                }
                                # Prepare dataframe for this sample
                                df_sample <- as.data.frame(lag_matrix)
                                # Compute spatial lags for covariates if available--------
                                # Add covariates
                                if (!is.null(lag_covariates)){
                                  df_lag_cov = sobj@meta.data%>%select(lag_covariates)%>%as.matrix
                                  
                                  # Compute spatial lags for ligands ------
                                  lag_cov_matrix <- matrix(NA, nrow = nrow(df_lag_cov), ncol = ncol(df_lag_cov))
                                  colnames(lag_cov_matrix) <- paste0("lag_", colnames(df_lag_cov))
                                  rownames(lag_cov_matrix) <- rownames(df_lag_cov)
                                  for (i in 1:ncol(df_lag_cov)) {
                                    lag_cov_matrix[, i] <- lag.listw(spatial_weights, df_lag_cov[, i], zero.policy = TRUE)
                                  }
                                  df_sample = cbind(df_sample, lag_cov_matrix)
                                }
                                
                                df_sample$response <- cell_state
                                df_sample$sample_id <- sobj@meta.data$sample_id
                                df_sample$region <- sobj@meta.data$region
                                df_sample$cell_name <- rownames(df_sample)
                                
                                if (!is.null(rec_expr_norm)) {
                                  colnames(rec_expr_norm) <- paste0("expr_", colnames(rec_expr_norm))
                                  df_sample <- cbind(df_sample, rec_expr_norm)
                                }
                                
                                # Add covariates
                                if (!is.null(covariates)){
                                  df_cov = sobj@meta.data%>%select(covariates)
                                  df_sample = cbind(df_sample , df_cov)
                                }
                                return(df_sample)
                                
                              })%>%setNames(.,names(sobj_list))
    
    if (length(combined_df_list) == 0) {
      log_error("No samples processed successfully")
      return(NULL)
    }
    
    # Combine all samples
    combined_df <- do.call(rbind, combined_df_list)
    
    # Subset cells if specified
    if (!is.null(selected_cells)) {
      combined_df <- combined_df %>% filter(cell_name %in% selected_cells)
      log_info("Using {nrow(combined_df)} selected cells for modeling")
    }
    
    log_info("Multi-sample spatial data processing completed: {nrow(combined_df)} total cells")
    
    return(combined_df)
    
  }, error = function(e) {
    log_error("Error in .process_multisample_spatial_data: {e$message}")
    return(NULL)
  })
}

##################################################################|
# Spatial Weights Creation Helper Function
##################################################################|

.create_spatial_weights <- function(coords, k) {

  tryCatch({
    knn_obj <- knearneigh(coords, k = k)
    nb <- knn2nb(knn_obj)
    listw <- nb2listw(nb, style = "W", zero.policy = TRUE)
    return(listw)
    
  }, error = function(e) {
    log_error("Error creating spatial weights: {e$message}")
    return(NULL)
  })
}

##################################################################|
# Mixed Effects Model Fitting Function
##################################################################|

.fit_mixed_effects_spatial_model <- function(
    combined_data,
    .ligand,
    .receptor = NULL,
    add_lr_interaction = TRUE,
    random_effect_structure = "sample_id",
    covariates = NULL,
    lag_covariates = NULL,
    model_engine = "glmmTMB",
    family_type = "binomial",
    zero_inflated,
    control_params = NULL
) {
  
  tryCatch({
    # Prepare formula components
    lag_terms <- paste0("lag_", .ligand)

    # Define fixed effects with  covariates
    if (!is.null(lag_covariates)){
      fixed_effects <- c(lag_terms, covariates, paste0("lag_", lag_covariates))
    }else{
      fixed_effects <- c(lag_terms, covariates)
    }
    
    
    # Add receptor terms if specified
    if (!is.null(.receptor)) {
      expr_rec_terms <- paste0("expr_", .receptor)
      fixed_effects <- c(fixed_effects, expr_rec_terms)
      
      # Add interactions if requested
      if (add_lr_interaction) {
        interaction_terms <- paste(paste0("`",lag_terms,"`"), paste0("`",expr_rec_terms,"`"), sep = ":")
        
      }
    }
    
    # Construct random effects
    if (random_effect_structure == "region") {
      random_effects <- "(1|region_id)"
    } else if (random_effect_structure == "sample_id") {
      random_effects <- "(1|sample_id)"
    }
    
    
    # Create formula
    formula_str <- paste("response ~", paste(paste0("`", fixed_effects, "`"), collapse = " + "), "+", random_effects)
    if (add_lr_interaction){
      formula_str <- paste0(formula_str, " + ", interaction_terms)
    }
    formula <- as.formula(formula_str)
    
    # Check for sufficient variation in response
    if (length(unique(combined_data$response)) < 2) {
      warning("Insufficient variation in response variable")
      return(list(
        model = NULL,
        formula = formula_str,
        error = "Insufficient response variation"
      ))
    }
    
    # Fit model based on engine
    if (model_engine == "glmmTMB") {
      # Set up family
      if (family_type == "binomial") {
        family_obj <- binomial()
      } else if (family_type == "gaussian") {
        family_obj <- gaussian()
      } else {
        family_obj <- get(family_type)()
      }
      
      # Set up control parameters
      if (is.null(control_params)) {
        control_params <- glmmTMBControl(optimizer = optim,
                                         optArgs = list(method = "BFGS",rel.tol = 1e-18,iter.max = 50000000, eval.max = 50000000))
      }
      
      if (zero_inflated){
        
        fit <- glmmTMB(formula = formula, 
                       data = combined_data, 
                       family = family_obj,
                       ziformula=~1,
                       control = control_params)
      }else{
        fit <- glmmTMB(formula = formula, 
                       data = combined_data, 
                       family = family_obj,
                       control = control_params)
      }
      
    } else if (model_engine == "lme4") {
      if (family_type == "binomial") {
        fit <- glmer(formula = formula, 
                     data = combined_data, 
                     family = binomial,
                     control = glmerControl(optimizer = "bobyqa"))
      } else {
        fit <- lmer(formula = formula, 
                    data = combined_data,
                    control = lmerControl(optimizer = "bobyqa"))
      }
    }
    
    return(list(
      model = fit,
      formula = formula_str,
      data = combined_data,
      model_engine = model_engine,
      family = family_type,
      random_effect_structure = random_effect_structure
    ))
    
  }, error = function(e) {
    log_error("Error fitting mixed-effects model: {e$message}")
    return(list(
      model = NULL,
      formula = formula_str,
      error = e$message
    ))
  })
}



library(broom)
library(broom.mixed)  # For mixed models
library(dplyr)

### Streamlined Extract Multi-Sample Model Summary Function ######################

extract_multisample_summary <- function(results, 
                                        adjust_method = "BH", 
                                        adjust_group.by = "predictor_type",
                                        include_random_effects = TRUE,
                                        include_glance = TRUE,
                                        covariates = "lag_nCount_Vizgen",
                                        rm.covariates = T,
                                        combine_effects = FALSE) {
  if (is.null(results)) {
    warning("Results is NULL")
    return(list(fixed_effects = data.frame()))
  }
  
  if (!is.list(results) || !is.null(results$model)) {
    results <- list(combined_model = results)
  }

  fixed_list <- list()
  random_list <- list()
  glance_list <- list()
  
  for (model_name in names(results)) {
    model_result <- results[[model_name]]
    curr_ligand <- str_split(model_name, "\\_",simplify = T)[1]
    curr_receptor <- str_split(model_name, "\\_",simplify = T)[1,2]
    if (is.null(model_result$model)) {
      warning(paste("Model", model_name, "is NULL, skipping"))
      next
    }
    
    model_obj <- model_result$model
    random_effect_var = model_obj[["modelInfo"]][["grpVar"]]
    fix_effect_var = paste0(model_obj[["frame"]]%>%colnames()%>%setdiff(c("response",random_effect_var)),collapse = " + ")
    formula_label = sprintf("%s + (1|%s)",fix_effect_var , random_effect_var)
    tryCatch({
      fixed_df <- tidy(model_obj, effects = "fixed",conf.int = T) %>%
        filter(term != "(Intercept)") %>%
        mutate(
          model_name = model_name,
          formula = formula_label,
          fix_effect_var =  fix_effect_var,
          effect_type = "fixed",
          predictor = paste0(term, " [", fix_effect_var, "]"),
          predictor_type = ifelse(term%in%covariates, "covariate",
                                  ifelse(str_detect(term,"lag_"),"ligand",
                                         ifelse(str_detect(term,"expr_"),"receptor", NA))),
          is_spatial_lag = grepl("lag_", term),
          is_interaction = grepl(":", term),
          term_type = case_when(
            is_spatial_lag ~ "spatial_lag",
            is_interaction ~ "interaction",
            grepl("_expr$", term) ~ "expression",
            TRUE ~ "other"
          ),
          group = NA_character_,
          effect = NA_character_
        )
      
      if(rm.covariates){
        fixed_df = filter(fixed_df, !term%in%covariates)
      }
      fixed_list[[model_name]] <- fixed_df
      
      if (include_random_effects) {
        random_pars <- tidy(model_obj, effects = "ran_pars",conf.int = T) %>%
          mutate(
            model_name = model_name,
            formula = formula_label,
            fix_effect_var =  fix_effect_var,
            effect_type = "random",
            predictor = paste0(term, " [", fix_effect_var, "]"),
            predictor_type = "random",
            is_spatial_lag = FALSE,
            is_interaction = FALSE,
            term_type = case_when(
              grepl("sd__", term) ~ "random_sd",
              grepl("var__", term) ~ "random_var", 
              grepl("cor__", term) ~ "random_cor",
              group == "Residual" ~ "residual",
              TRUE ~ "random_other"
            )
          )
        
        random_list[[model_name]] <- random_pars
      }
      
      if (include_glance) {
        glance_df <- glance(model_obj) %>%
          mutate(model_name = model_name)
        glance_list[[model_name]] <- glance_df
      }
      
    }, error = function(e) {
      warning(paste("Error processing model", model_name, ":", e$message))
    })
  }
  
  output <- list()
  
  if (length(fixed_list) > 0) {
    combined_fixed <- bind_rows(fixed_list) %>%
      group_by(.data[[adjust_group.by]])%>%
      # Apply multiple testing correction only to fixed effects
      mutate(
        p_adj = if_else(effect_type == "fixed", 
                        p.adjust(p.value, method = adjust_method), 
                        NA_real_),
        significant = if_else(effect_type == "fixed", 
                              p_adj < 0.05, 
                              NA),
        effect_size = if_else(effect_type == "fixed",
                              categorize_effect_size(estimate),
                              NA_character_)
      ) %>%
      arrange(effect_type, p_adj, desc(abs(estimate)))
    
    output$fixed_effects <- combined_fixed
  }
  
  if (length(random_list) > 0) {
    output$random_effects <- bind_rows(random_list)
  }
  
  if (length(glance_list) > 0) {
    output$model_fit <- bind_rows(glance_list)
  }
  
  if (combine_effects && length(fixed_list) > 0 && length(random_list) > 0) {
    all_effects <- bind_rows(
      combined_fixed, #
      bind_rows(random_list)
    ) %>%arrange(model_name, effect_type, p_adj, desc(abs(estimate)))

    return(all_effects)
  }
  
  return(output)
}

categorize_effect_size <- function(estimates) {
  abs_est <- abs(estimates)
  case_when(
    is.na(abs_est) ~ "unknown",
    abs_est < 0.1 ~ "negligible", 
    abs_est < 0.3 ~ "small",
    abs_est < 0.5 ~ "medium",
    TRUE ~ "large"
  )
}



### Create Comparison Forest Plots with VIP Ordering and Directional Coloring ###############
create_comparison_forests = function(successful_summaries, main_output_dir, 
                                     effect.size.var = "estimate",
                                     p.var = "p_adj", 
                                     ci.low.var = "conf.low", 
                                     ci.up.var = "conf.high",
                                     vip_scores = NULL) {
  
  if(nrow(successful_summaries) == 0) {
    log_warn("No data provided for forest plots")
    return(NULL)
  }
  
  # Directional colors (red = positive, blue = negative)
  direction_colors_fill <- c(
    "Positive" = "#d73027",
    "Negative" = "#4575b4",
    "Non-significant" = "gray85"
  )

  direction_colors_border <- c("Positive" = "#d73027",
                               "Negative" = "#4575b4",
                               "Non-significant" = "gray85")
  
  successful_summaries <- successful_summaries %>%
    filter(effect_type == "fixed") %>%
    filter(term != "lag_nCount_Vizgen") %>%
    filter(!is.na(.data[[ci.low.var]]), !is.na(.data[[ci.up.var]])) %>%
    mutate(
      ligand = str_split(model_name, "_", simplify = T)[,1],
      receptor = str_split(model_name, "_", simplify = T)[,2]
    )
  
  if (!is.null(vip_scores)) {
    # Get unique ligand-VIP pairs to avoid many-to-many
    vip_unique <- vip_scores %>% 
      select(ligand, vip) %>% 
      distinct(ligand, .keep_all = TRUE)
    
    successful_summaries <- successful_summaries %>%
      left_join(vip_unique, by = "ligand", relationship = "many-to-one") %>%
      mutate(
        vip = ifelse(is.na(vip), 0, vip),
        vip_group = factor(
          ifelse(vip >= 1, "Positive VIP (>=1)", "Negative VIP (<1)"),
          levels = c("Positive VIP (>=1)", "Negative VIP (<1)")
        )
      )
  } else {
    successful_summaries <- successful_summaries %>%
      mutate(
        vip = 0,
        vip_group = factor("All Ligands", levels = "All Ligands")
      )
    log_warn("No VIP scores provided - plotting without VIP grouping")
  }
  
  successful_summaries <- successful_summaries %>%
    mutate(
      direction = factor(
        case_when(
          is.na(.data[[p.var]]) | .data[[p.var]] >= 0.05 ~ "Non-significant",
          .data[[effect.size.var]] > 0 ~ "Positive",
          .data[[effect.size.var]] < 0 ~ "Negative",
          TRUE ~ "Non-significant"
        ),
        levels = c("Positive", "Negative", "Non-significant")
      )
    )
  
  successful_summaries <- successful_summaries %>%
    arrange(vip_group, desc(abs(vip)), desc(abs(.data[[effect.size.var]]))) %>%
    mutate(
      model_name = factor(model_name, levels = unique(model_name))
    )
  
  
  n_pred = length(unique(successful_summaries$term))
  n_modeltype = length(unique(successful_summaries$model_type))
  n_family_type = length(unique(successful_summaries$family_type))

  separate_configs <- successful_summaries %>%
    filter(model_type == "separate_lr_pairs") %>% 
    pull(config_id) %>% 
    unique()
  
  if(length(separate_configs) > 0) {
    fig1_plots <- lapply(separate_configs, function(c){
      .df = successful_summaries %>% filter(config_id == c)
      
      plot_lab = sprintf("%s\n%s-%s", unique(.df$config_id), unique(.df$input_cells),
                         unique(.df$family_type))

      p <- ggplot(.df,
                  aes(x = .data[[effect.size.var]], y = model_name,
                      xmin = .data[[ci.low.var]], xmax = .data[[ci.up.var]])) + 
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.6) +
        geom_linerange(aes(color = direction), linewidth = 4, 
                       position = position_dodge(width = 0.5)) +
        geom_point(aes(fill = direction), size = 3.5, shape = 21, color = "white", 
                   stroke = 0.5, position = position_dodge(width = 0.5)) +
        scale_fill_manual(values = direction_colors_fill, name = "Association") +
        scale_color_manual(values = direction_colors_border, name = "Association") +
        scale_y_discrete(name = "Ligand-Receptor Pair") +
        scale_x_continuous(name = "Effect Size") +
        facet_grid(vip_group ~ predictor_type, scales = "free", space = "free_y") +
        theme_minimal(base_size = 12) +
        ggtitle(plot_lab) +
        theme(
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major.x = element_line(color = "grey90"),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 10),
          legend.position = "top",
          legend.title = element_text(face = "bold", size = 11),
          strip.text.y = element_text(angle = 0, hjust = 0, face = "bold", size = 10),
          strip.text.x = element_text(face = "bold", size = 10),
          strip.background = element_rect(fill = "grey95", color = "black"),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
        )
      return(p)
    })
    
    fig1 <- tryCatch({
      ggpubr::ggarrange(plotlist = fig1_plots,
                        nrow = max(1, n_modeltype),
                        ncol = max(1, n_family_type),
                        common.legend = TRUE, 
                        legend = "top")
    }, error = function(e) {
      log_warn("Could not create combined plot for separate_lr_pairs: {e$message}")
      if(length(fig1_plots) == 1) return(fig1_plots[[1]]) else return(NULL)
    })
  } else {
    fig1 <- NULL
    log_info("No separate_lr_pairs models to plot")
  }
  
  combined_configs <- successful_summaries %>%
    filter(model_type == "combined_lr_pairs") %>% 
    pull(config_id) %>% 
    unique()
  
  if(length(combined_configs) > 0) {
    fig2_plots <- lapply(combined_configs, function(c){
      .df = successful_summaries %>% filter(config_id == c)
      
      plot_lab = sprintf("%s\n%s-%s", unique(.df$config_id), unique(.df$input_cells),
                         unique(.df$family_type))

      .df <- .df %>%
        arrange(vip_group, desc(vip), desc(abs(.data[[effect.size.var]]))) %>%
        mutate(term = factor(term, levels = unique(term)))

      p <- ggplot(.df,
                  aes(x = .data[[effect.size.var]], y = term,
                      xmin = .data[[ci.low.var]], xmax = .data[[ci.up.var]])) + 
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.6) +
        geom_linerange(aes(color = direction), linewidth = 2.5, 
                       position = position_dodge(width = 0.5)) +
        geom_point(aes(fill = direction), size = 4, shape = 21, color = "black", 
                   stroke = 0.5, position = position_dodge(width = 0.5)) +
        scale_fill_manual(values = direction_colors_fill, name = "Association") +
        scale_color_manual(values = direction_colors_border, name = "Association") +
        scale_y_discrete(name = "Term") +
        scale_x_continuous(name = "Effect Size") +
        facet_grid(vip_group ~ ., scales = "free", space = "free_y") +
        theme_minimal(base_size = 12) +
        ggtitle(plot_lab) +
        theme(
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major.x = element_line(color = "grey90"),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 10),
          legend.position = "top",
          legend.title = element_text(face = "bold", size = 11),
          strip.text.y = element_text(angle = 0, hjust = 0, face = "bold", size = 10),
          strip.background = element_rect(fill = "grey95", color = "black"),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
        )
      return(p)
    })
    
    fig2 <- tryCatch({
      ggpubr::ggarrange(plotlist = fig2_plots,
                        nrow = max(1, n_modeltype),
                        ncol = max(1, n_family_type),
                        common.legend = TRUE, 
                        legend = "top")
    }, error = function(e) {
      log_warn("Could not create combined plot for combined_lr_pairs: {e$message}")
      if(length(fig2_plots) == 1) return(fig2_plots[[1]]) else return(NULL)
    })
  } else {
    fig2 <- NULL
    log_info("No combined_lr_pairs models to plot")
  }
  
  if(!is.null(fig1) && !is.null(fig2)) {
    p1 <- tryCatch({
      ggpubr::ggarrange(fig1, fig2, nrow = 1, common.legend = TRUE, legend = "top")
    }, error = function(e) {
      log_warn("Could not combine final plots: {e$message}. Saving separately.")
      return(list(fig1 = fig1, fig2 = fig2))
    })
  } else if(!is.null(fig1)) {
    p1 <- fig1
  } else if(!is.null(fig2)) {
    p1 <- fig2
  } else {
    log_warn("No plots generated")
    return(NULL)
  }
  
  tryCatch({
    if(is.list(p1) && !inherits(p1, "ggplot")) {
      if(!is.null(p1$fig1)) {
        ggsave(file.path(main_output_dir, "forest_plot_separate_lr_pairs_vip.pdf"), p1$fig1,
               width = 22, height = max(8, 1 + n_pred * 0.4), dpi = 300, bg = "white", limitsize = FALSE)
      }
      if(!is.null(p1$fig2)) {
        ggsave(file.path(main_output_dir, "forest_plot_combined_lr_pairs_vip.pdf"), p1$fig2,
               width = 22, height = max(8, 1 + n_pred * 0.4), dpi = 300, bg = "white", limitsize = FALSE)
      }
    } else {
      ggsave(file.path(main_output_dir, "forest_plot_by_model_type_vip.pdf"), p1,
             width = 24, height = max(10, 1 + n_pred * 0.4), dpi = 300, bg = "white", limitsize = FALSE)
    }
    log_info("VIP-ordered forest plots saved successfully")
  }, error = function(e) {
    log_error("Failed to save forest plots: {e$message}")
  })
  
  return(p1)
}

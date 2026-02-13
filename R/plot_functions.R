#' @import tidyverse
plot_DeltaCTSE <- function(ctse_res_mod, ctse_res_ori, title = NULL, subtitle = NULL,
                           group_colors = NULL, mod_name = "Modified", ori_name = "Original"){
  
  cells <- names(ctse_res_ori$ExpMAE_Summary)
  
  if (!is.character(title)){
    title = paste0(mod_name, " Performance Gain over ", ori_name)
  }
  
  change_df <- data.frame(
    cell_group = cells, 
    ExpSCorr = (ctse_res_mod$ExpSCorr_Summary[cells] - ctse_res_ori$ExpSCorr_Summary[cells]), 
    CvSpe = (ctse_res_mod$CvSpe_Summary[cells] - ctse_res_ori$CvSpe_Summary[cells]),
    ExpMAE = (ctse_res_mod$ExpMAE_Summary[cells] - ctse_res_ori$ExpMAE_Summary[cells]),
    ExpRMSE = (ctse_res_mod$ExpRMSE_Summary[cells] - ctse_res_ori$ExpRMSE_Summary[cells])
  )
  
  change_long <- change_df %>%
    pivot_longer(cols = -cell_group, names_to = "Metric", values_to = "Delta") %>%
    group_by(Metric) %>%
    mutate(
      # Direction of improvement
      Is_Improved = case_when(
        Metric == "ExpSCorr" ~ Delta > 0,
        Metric %in% c("CvSpe", "ExpMAE", "ExpRMSE") ~ Delta < 0
      ),
      # STANDARDIZATION: Scale the deltas within each metric to a common intensity
      # We use abs() because we want to see the magnitude of change for the color
      Scaled_Intensity = abs(as.numeric(scale(Delta)))
    ) %>%
    ungroup() %>%
    # Create a score where positive is 'Good Change' and negative is 'Bad Change'
    mutate(Improvement_Score = ifelse(Is_Improved, Scaled_Intensity, -Scaled_Intensity))
  
  
  if (!is.vector(group_colors)){
    ggplot(change_long, aes(x = Delta, y = cell_group, fill = Improvement_Score)) +
      geom_col() +
      facet_wrap(~Metric, scales = "free_x", nrow = 1) +
      scale_fill_gradient2(
        low = "#B2182B", mid = "#F7F7F7", high = "#1B7837",
        midpoint = 0, name = "Intensity of\nImprovement"
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4) +
      theme_minimal() +
      labs(
        title = title,
        subtitle = subtitle,
        x = paste0("Difference (", mod_name, " - ", ori_name,")"),
        y = NULL
      ) +
      theme(
        strip.text = element_text(face = "bold", size = 12),
        # Keep text black, but use the color vector for the axis ticks or margins
        axis.text.y = element_text(size = 9), 
        # This creates a colored 'border' effect on the left of the text
        panel.spacing = unit(1.5, "lines")
      )
  } else {
    label_colors <- group_colors[levels(factor(change_long$cell_group))]
    
    ggplot(change_long, aes(x = Delta, y = cell_group, fill = Improvement_Score)) +
      geom_col() +
      facet_wrap(~Metric, scales = "free_x", nrow = 1) +
      scale_fill_gradient2(
        low = "#B2182B", mid = "#F7F7F7", high = "#1B7837",
        midpoint = 0, name = "Intensity of\nImprovement"
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4) +
      theme_minimal() +
      labs(
        title = title,
        subtitle = subtitle,
        x = paste0("Difference (", mod_name, " - ", ori_name,")"),
        y = NULL
      ) +
      theme(
        strip.text = element_text(face = "bold", size = 12),
        # Keep text black, but use the color vector for the axis ticks or margins
        axis.text.y = element_text(size = 9), 
        # This creates a colored 'border' effect on the left of the text
        axis.ticks.y = element_line(color = label_colors, size = 2),
        axis.ticks.length.y = unit(0.3, "cm"),
        panel.spacing = unit(1.5, "lines")
      )
  }
}

#' @import tidyverse
#' @import scales
plot_dumbbellCTSE <- function(ctse_res_list, modified, original, 
                              truth = "Ground Truth", 
                              baseline = "Plain Bulk", 
                              title = NULL, subtitle = NULL, group_colors = NULL) {
  
  if (!is.character(title)){
    title = paste0(modified, " Improvement over ", original)
  }
  if (!is.character(subtitle)){
    subtitle = paste0("Track starts at ", baseline, " and ends at ", truth)
  }
  
  all_possible_metrics <- c("CvSpe", "ExpMAE", "ExpRMSE", "ExpSCorr")
  available_summaries <- names(ctse_res_list[[original]])
  
  metrics_to_process <- all_possible_metrics[paste0(all_possible_metrics, "_Summary") %in% available_summaries]
  
  ctse_res_df <- data.frame()
  
  methods_to_plot <- c(truth, baseline, original, modified)
  methods_to_plot <- methods_to_plot[methods_to_plot %in% names(ctse_res_list)]
  
  for (m in methods_to_plot) {
    first_metric <- paste0(metrics_to_process[1], "_Summary")
    cell_labels <- names(ctse_res_list[[m]][[first_metric]])
    
    tmp_df <- data.frame(
      method = m,
      cell_label = cell_labels,
      stringsAsFactors = FALSE
    )
    
    for (met in metrics_to_process) {
      summary_name <- paste0(met, "_Summary")
      vals <- as.numeric(ctse_res_list[[m]][[summary_name]])
      
      if (met == "CvSpe" && m %in% c(baseline, "Fraction-Regressed Bulk")) {
        vals[is.na(vals)] <- 1
      }
      
      tmp_df[[met]] <- vals
    }
    
    ctse_res_df <- rbind(ctse_res_df, tmp_df)
  }
  
  plot_data <- ctse_res_df %>%
    pivot_longer(cols = all_of(metrics_to_process), 
                 names_to = "Metric", 
                 values_to = "Value") %>%
    pivot_wider(names_from = method, values_from = Value)
  
  orig_sym <- sym(original)
  mod_sym <- sym(modified)
  
  plot_data <- plot_data %>%
    mutate(
      Improvement = case_when(
        Metric %in% c("CvSpe", "ExpMAE", "ExpRMSE") ~ (!!orig_sym - !!mod_sym) / abs(!!orig_sym),
        Metric == "ExpSCorr" ~ (!!mod_sym - !!orig_sym) / abs(!!orig_sym),
        TRUE ~ 0
      )
    ) %>%
    # Cap improvement at 1.5 (150%) to keep legend readable and prevent outlier stretching
    mutate(Improvement_Capped = pmax(pmin(Improvement, 1.5, na.rm = TRUE), -1.5, na.rm = TRUE))
  
  limit_val <- 1.5
  val_points <- c(-limit_val, -1, 0, 1, limit_val)
  col_points <- c("#67000d", "#fb6a4a", "darkgrey", "#addd8e", "#00441b")
  rescaled_points <- rescale(val_points, from = c(-limit_val, limit_val))
  
  base_sym <- sym(baseline)
  truth_sym <- sym(truth)
  
  
  if (!is_vector(group_colors)){
    plot_theme <- theme(
      strip.text = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 9),
      panel.spacing = unit(1.5, "lines")
      )
  } else {
    label_colors <- group_colors[levels(factor(plot_data$cell_label))]
    plot_theme <- theme(
      strip.text = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 9),
      panel.spacing = unit(1.5, "lines"),
      axis.ticks.y = element_line(color = label_colors, size = 2),
      axis.ticks.length.y = unit(0.3, "cm"),
      )
  }
  
  p <- ggplot(plot_data, aes(y = cell_label)) +
    # Background "Track" (Baseline to Truth)
    geom_segment(aes(x = !!base_sym, xend = !!truth_sym, yend = cell_label), 
                 color = "grey", linewidth = 1.2) +
    
    # The Improvement Arrow (Original to Modified)
    geom_segment(aes(x = !!orig_sym, xend = !!mod_sym, yend = cell_label, color = Improvement_Capped), 
                 linewidth = 1.8, 
                 arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
    
    # Markers
    geom_point(aes(x = !!base_sym), color = "grey", size = 2.5) +
    geom_point(aes(x = !!orig_sym, color = Improvement_Capped), size = 2) +
    geom_point(aes(x = !!truth_sym), shape = 18, color = "grey", size = 3.5) +
    
    facet_wrap(~Metric, scales = "free_x", nrow = 1) +
    
    scale_color_gradientn(
      colors = col_points,
      values = rescaled_points,
      limits = c(-limit_val, limit_val),
      name = "Relative\nImprovement",
      breaks = c(-1.5, -1, 0, 1, 1.5),
      labels = c("<<<<", "-100%", "0", "+100%", ">>>>")
    ) +
    theme_minimal() +
    labs(x = "Metric Value", y = NULL,
         title = title,
         subtitle = subtitle) +
    plot_theme
  
  return(p)
}

#' @import viridis
plot_ExpSpe <- function(ctse_res_mod, ctse_res_ori, gene_pressure, title = NULL, subtitle = NULL, mod_name = "Modified", ori_name = "Original"){
  if (!is.character(title)){
    title = "Expression Specificity by Lin’s Concordance Correlation Coefficient (CCC)"
  }
  
  gene_intensity <- apply(abs(gene_pressure), 2, max)
  
  spe_df_wide <- data.frame(
    Gene = names(ctse_res_mod$ExpSpe_Vector),
    BP_CCC = ctse_res_ori$ExpSpe_Vector,
    DP_CCC = ctse_res_mod$ExpSpe_Vector
  )
  
  spe_df_wide$Interaction_Intensity <- gene_intensity[match(spe_df_wide$Gene, names(gene_intensity))]
  
  ggplot(spe_df_wide, aes(x = BP_CCC, y = DP_CCC, color = Interaction_Intensity)) +
    geom_point(alpha = 0.5, size = 1.2) +
    # Add the y=x line (identity line)
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", alpha = 0.6) +
    # Use a perceptually uniform color scale (magma or viridis)
    scale_color_viridis_c(option = "magma", direction = -1, name = "Max Abs\nGene Pressure") +
    theme_minimal() +
    labs(
      title = title,
      subtitle = subtitle,
      x = paste0(ori_name, " CCC"),
      y = paste0(mod_name, " CCC")
    ) +
    theme(legend.position = "right")
}
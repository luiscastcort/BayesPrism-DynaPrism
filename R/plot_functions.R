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
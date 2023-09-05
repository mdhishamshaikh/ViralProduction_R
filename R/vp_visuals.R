#' Visualizations of viral production data
#' 
#' @description
#' A major step in data analyses is `Data Visualization`. Some suggestions to visualize the flow cytometry, viral
#' production or analyzed viral production data are given. 
#' 
#' `plot_overview_counts_over_time`: plots the bacterial and viral counts, retrieved from the flow cytometry, in
#' function of the different time ranges of the assay. The bacterial endpoint is also highlighted in pink. 
#' 
#' `plot_collision_rates`: plots the difference in collision rates over time between VP and VPC samples. The
#' collision rate is defined as the frequency of physical encounter between a bacteriophage and a bacterial cell in the sample.
#' 
#' `plot_comparison_methods`: compares the different methods to calculate viral production. The calculated viral
#' production is compared between the linear regression variants and the VIPCAL variants, also an overview of
#' all the methods is made. Next to that, a plot that compares linear regression with VIPCAL and VIPCAL-SE is also produced.
#'
#' `plot_percentage_cells`: plots the percentage of lytically infected and lysogenic cells in function of the 
#' different burst sizes. Important to note, that this is the percentage of cells at the bacterial endpoint of the assay,
#' so the point where ideally the assay is stopped to retrieve less biased results between samples. 
#' 
#' `plot_nutrient_release`: plots the nutrient release for organic carbon, nitrogen and phosphor in function of the burst
#' size at the end of the assay (T0_T24).
#'
#' @param data Data frame with the output of the flow cytometry.
#' @param original_abundances Data frame with the abundances of bacterial and virus population in the original sample. 
#' @param vp_results Data frame with the viral production calculation results, available in global environment 
#' after running \code{calculate_viral_production}.
#' @param analyzed_vp_results_bacterial_endpoint Data frame with the analyzed viral production results, available in global environment 
#' after running \code{analyze_viral_production}. Important is that the analyzing step is ran on the \code{vp_results_ouput_BP_df}
#' data frame from \code{calculate_viral_production} instead of the \code{vp_results_ouput_BP_df} so that the bacterial
#' endpoint is taken into account. 
#' @param analyzed_vp_results_T0_T24 Data frame with the analyzed viral production results, available in global environment 
#' after running \code{analyze_viral_production}. Important is that the analyzing step is ran on the \code{vp_results_ouput_T24_df}
#' data frame from \code{calculate_viral_production} instead of the \code{vp_results_ouput_BP_df} so that the nutrient release
#' at the end of the assay is calculated. 
#' 
#' @return Plot objects will be stored in variable `plot_list` in the global environment.
#' 
#' @name vp_visuals
#' @rdname vp_visuals
#'
#' @examples \dontrun{
#' # Setup
#' data_NJ2020_all <- read.csv(system.file('extdata', 
#' 'NJ2020_Station_2_and_6_all_populations.csv', package = "viralprod"))
#' 
#' original_abundances_NJ2020 <- read.csv(system.file('extdata',
#' 'NJ2020_original_abundances.csv', package = "viralprod"))
#' 
#' calculate_viral_production(data_NJ2020_all, write_csv = F)
#' 
#' .GlobalEnv$plot_list <- list()
#' 
#' # Visuals
#' plot_overview_counts_over_time(data_NJ2020_all)
#' plot_collision_rates(data_NJ2020_all, original_abundances_NJ2020)
#' plot_comparison_methods(vp_results_output_df) 
#' 
#' analyze_viral_production(vp_results_output_BP_df, data_NJ2020_all, 
#' original_abundances_NJ2020, write_csv = F)
#' plot_percentage_cells(analyzed_vp_results_df)
#' 
#' analyze_viral_production(vp_results_output_T24_df, data_NJ2020_all, 
#' original_abundances_NJ2020, write_csv = F)
#' plot_nutrient_release(analyzed_vp_results_df)
#' }
plot_overview_counts_over_time <- function(data){
  ## 1. Setup
  plot_dataframe <- data %>%
    vp_average_replicate_dataframe() %>%
    dplyr::mutate(Sample_Type = factor(.data$Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                  Plot_Group = ifelse(.data$Population %in% c('c_Bacteria', 'c_Viruses'), 'Total',
                                      ifelse(.data$Microbe == 'Bacteria', 'Bacteria', 'Viruses')),
                  Plot_Group = factor(.data$Plot_Group, levels = c('Total', 'Bacteria', 'Viruses')))
  
  for (combi_tag in unique(plot_dataframe$tag)){
    plot_dataframe_counts <- plot_dataframe %>%
      dplyr::filter(.data$tag == combi_tag) %>%
      dplyr::mutate(Time_Time = factor(.data$Time_Time, levels = unique(.data$Time_Time)))
    
    plot_highlight_bacterial_endpoint <- data %>%
      tidyr::unite(dplyr::all_of(c('Location', 'Station_Number', 'Depth')), col = 'tag', remove = F) %>%
      dplyr::filter(.data$tag == combi_tag)
    
    bacterial_endpoint_index <- vp_bacterial_endpoint(plot_highlight_bacterial_endpoint, visual = T)
    
    ## 2. Make plot
    n <- ggplot2::ggplot(plot_dataframe_counts, ggplot2::aes(x = .data$Timepoint, y = .data$Mean, 
                                                             color = .data$Population, shape = .data$Population)) + 
      ggplot2::geom_point(size = 1.5) + 
      ggplot2::geom_line() + 
      ggplot2::geom_smooth(linewidth = 1.0, method = 'lm', se = F) + 
      ggplot2::geom_hline(yintercept = 0, color = '#000000', size = 0.3, linetype = 'dashed') + 
      ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$Mean - .data$SE, ymax = .data$Mean + .data$SE), 
                             width = 0.5, size = 0.5) + 
      
      ggplot2::facet_grid(.data$Plot_Group + .data$Sample_Type ~ .data$Time_Time) + 
      
      ggplot2::scale_y_continuous(labels = scales::unit_format(unit = NULL, scale = 1e-6)) + 
      ggplot2::scale_x_continuous(breaks = unique(plot_dataframe_counts$Timepoint)) + 
      ggplot2::scale_shape_manual(name = 'Populations', 
                                  labels = unique(plot_dataframe_counts$Population),
                                  values = c(16,17,16,17,15,16,17)) +
      ggplot2::scale_color_manual(name = 'Populations',
                                  labels = unique(plot_dataframe_counts$Population),
                                  values = c("#FF0000", "#0000FF", "#0099CC", "#666666", 
                                             "#666633", "#669900", "#FF0066")) + 
      
      ggplot2::labs(title = paste0(plot_dataframe_counts$tag, " - Overview"),
                    subtitle = 'Bacterial and Viral counts for Lytic and Lysogenic inductions',
                    x = 'Sampling Timepoints\n (in hours)',
                    y = 'FCM Counts\n (in millions)') + 
      
      ggplot2::theme_bw() + 
      ggplot2::theme(strip.background = ggplot2::element_rect(color = 'black', fill = 'white'),
                     axis.title = ggplot2::element_text(face = 'bold'),
                     title = ggplot2::element_text(face = 'bold'),
                     plot.subtitle = ggplot2::element_text(face = 'plain'),
                     legend.position = "bottom") +
      
      ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = FALSE))
    
    ## 3. Highlight bacterial endpoint
    n_gtable <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(n)) 
    strip_index <- which(grepl('strip-', n_gtable$layout$name))
    
    fills <- c(rep("white", length(strip_index)))
    fills[bacterial_endpoint_index] <- "#FF6699"
    
    color <- 1
    for (index in strip_index) {
      current_facet_strip <- which(grepl('rect', n_gtable$grobs[[index]]$grobs[[1]]$childrenOrder)) 
      n_gtable$grobs[[index]]$grobs[[1]]$children[[current_facet_strip]]$gp$fill <- fills[color]
      color <- color + 1 
    }
    
    plot_name <- paste0(unique(plot_dataframe_counts$Location), '_Station_', 
                        unique(plot_dataframe_counts$Station_Number), '_Depth_', 
                        unique(plot_dataframe_counts$Depth), '_Overview')
    .GlobalEnv$plot_list[[plot_name]] <- n_gtable
  }
}


# Helper function of plot_collision_rates for adding second color scale.
#' @noRd
draw_key_cust <- function(data, params, size){
  if (data$colour == "#996633") {
    ggplot2::draw_key_vpath(data, params, size)
  } else {
    ggplot2::draw_key_path(data, params, size)
  }
}


#' @rdname vp_visuals
plot_collision_rates <- function(data, 
                                 original_abundances){
  ## 1. Setup
  # Get correct data frame
  plot_dataframe <- data %>%
    dplyr::select(dplyr::all_of(c("Location", "Station_Number", "Depth", "Sample_Type", 
                                  "Timepoint", "Replicate", "c_Bacteria", "c_Viruses")))
  
  original_abundances_DF <- original_abundances %>%
    dplyr::rename(c_Bacteria = 'Total_Bacteria',
                  c_Viruses = 'Total_Viruses') %>%
    dplyr::select(dplyr::all_of(c('Station_Number', 'c_Bacteria', 'c_Viruses')))
  
  abundance_dimensions <- data.frame(expand.grid('Location' = unique(plot_dataframe$Location),
                                                 'Station_Number' = unique(plot_dataframe$Station_Number),
                                                 'Depth' = unique(plot_dataframe$Depth),
                                                 'Sample_Type' = unique(plot_dataframe$Sample_Type),
                                                 'Timepoint' = -3,
                                                 'Replicate' = unique(plot_dataframe$Replicate)))
  
  abundances_dataframe <- dplyr::left_join(abundance_dimensions, original_abundances_DF, by = 'Station_Number')
  plot_dataframe <- dplyr::full_join(abundances_dataframe, plot_dataframe) %>%
    tidyr::unite(dplyr::all_of(c('Location', 'Station_Number', 'Depth')), col = 'tag', remove = F)
  
  # Collision rates
  plot_dataframe$BV <- plot_dataframe$c_Bacteria * plot_dataframe$c_Viruses
  collision_rate_results <- list()
  
  for (combi_tag in unique(plot_dataframe$tag)){
    for (sample in unique(plot_dataframe$Sample_Type)){
      for (rep in unique(plot_dataframe$Replicate)){
        DF <- plot_dataframe %>%
          dplyr::filter(.data$tag == combi_tag, .data$Sample_Type == sample, .data$Replicate == rep)
        
        for (time in unique(DF$Timepoint)){
          collision_rate <- (DF %>% dplyr::filter(.data$Timepoint == time))$BV / (DF %>% dplyr::filter(.data$Timepoint == 0))$BV
          result <- c(combi_tag, sample, rep, time, collision_rate)
          collision_rate_results[[length(collision_rate_results) + 1]] <- result
        }
      }
    }
  }
  collision_rate_dataframe <- data.table::data.table(t(data.table::as.data.table(collision_rate_results)))
  colnames(collision_rate_dataframe) <- c("tag", "Sample_Type", "Replicate", "Timepoint", "Collision_Rate")
  
  collision_rate_dataframe <- collision_rate_dataframe %>%
    dplyr::mutate_at(c('Collision_Rate', 'Timepoint'), as.numeric) %>%
    dplyr::group_by(.data$tag, .data$Sample_Type, .data$Timepoint) %>%
    dplyr::filter(.data$Sample_Type != '0.22') %>%
    dplyr::summarise(CR_Mean = mean(.data$Collision_Rate)) %>%
    dplyr::left_join(plot_dataframe %>%
                       dplyr::select(dplyr::all_of(c('tag', 'Location', 'Station_Number', 'Depth'))), 
                     by = 'tag', relationship = 'many-to-many') %>%
    dplyr::distinct()
  
  collision_rates_plot_df <- collision_rate_dataframe %>%
    dplyr::filter(.data$Timepoint != -3)
  
  collision_rates_plot_df_diff <- collision_rates_plot_df %>%
    tidyr::pivot_wider(names_from = 'Sample_Type',
                       values_from = 'CR_Mean') %>%
    dplyr::mutate(CR_Diff = .data$VP - .data$VPC)
  
  # Percentage of decrease of collision rates
  decrease_results <- list()
  
  for (station in unique(collision_rate_dataframe$Station_Number)){
    for (sample in unique(collision_rate_dataframe$Sample_Type)){
      DF <- collision_rate_dataframe %>%
        dplyr::filter(.data$Station_Number == station, .data$Sample_Type == sample, .data$Timepoint %in% c(-3,0))
      
      decrease <- (DF[DF$Timepoint == -3,]$CR_Mean - DF[DF$Timepoint == 0,]$CR_Mean) / DF[DF$Timepoint == -3,]$CR_Mean * 100
      result <- c(station, sample, decrease)
      decrease_results[[length(decrease_results)+1]]<- result
    }
  }
  decrease_dataframe <- data.frame(t(sapply(decrease_results, c)))
  colnames(decrease_dataframe) <- c('Station_Number', 'Sample_Type', 'Perc_Decrease')
  decrease_dataframe <- decrease_dataframe %>%
    tidyr::pivot_wider(names_from = 'Sample_Type',
                       values_from = 'Perc_Decrease') %>%
    dplyr::mutate_at(c('VP', 'VPC'), as.numeric)
  
  # Bacterial endpoint
  bp_dataframe <- data %>% 
    tidyr::unite(dplyr::all_of(c('Location', 'Station_Number', 'Depth')), col = 'tag', remove = F)
  bp_results <- list()
  
  for (combi_tag in unique(bp_dataframe$tag)){
    DF <- bp_dataframe %>%
      dplyr::filter(.data$tag == combi_tag)
    
    unique_timepoints <- unique(DF$Timepoint)
    
    bacterial_endpoint_index <- vp_bacterial_endpoint(DF, visual = T) + 1 # +1 since we aren't looking at Time_Ranges but at unique Timepoints
    bacterial_endpoint_timepoint <- unique_timepoints[bacterial_endpoint_index]
    bp_results[[length(bp_results)+1]] <- bacterial_endpoint_timepoint
  }
  bacterial_endpoint_dataframe <- data.frame(t(as.data.frame(bp_results))) %>%
    dplyr::mutate(Station_Number = unique(bp_dataframe$Station_Number))
  colnames(bacterial_endpoint_dataframe) <- c('Bacterial_Endpoint', 'Station_Number')
  
  ## 2. Make plot
  n <- ggplot2::ggplot() + 
    ggplot2::geom_point(data = collision_rates_plot_df, 
                        ggplot2::aes(x = .data$Timepoint, y = .data$CR_Mean, color = .data$Sample_Type, shape = .data$Sample_Type),
                        show.legend = T) + 
    ggplot2::geom_text(data = decrease_dataframe, ggplot2::aes(x = 11, y = 12),
                       label = paste0('VP Decrease: ', round(decrease_dataframe$VP,2), ' %\nVPC Decrease: ', round(decrease_dataframe$VPC,2), ' %'),
                       size = 3, color = 'black', show.legend = F) +
    
    ggplot2::scale_color_manual(values = c('#FF0000', '#0000FF')) +
    ggplot2::scale_fill_manual(values = c('#FF0000', '#0000FF')) +
    ggplot2::scale_shape_manual(values = c(15, 17)) +
    
    ggplot2::guides(color = ggplot2::guide_legend(title = 'Treatment', order = 1),
                    shape = ggplot2::guide_legend(title = 'Treatment', order = 1),
                    fill = ggplot2::guide_legend(title = 'Treatment', order = 1)) +
    
    ggnewscale::new_scale_color() +
    ggplot2::geom_line(data = collision_rates_plot_df_diff, 
                       ggplot2::aes(x = .data$Timepoint, y = .data$CR_Diff, color = "lineCR"), 
                       size = 1, alpha = 0.5, key_glyph = 'cust') +
    ggplot2::geom_vline(data = bacterial_endpoint_dataframe, 
                        ggplot2::aes(xintercept = .data$Bacterial_Endpoint, color = "lineBP"), 
                        linewidth = 1.5, alpha = 0.5, key_glyph = 'cust') +
    
    ggplot2::scale_color_manual(labels = c(lineCR = "Difference in collision rates between VP and VPC treatment", 
                                           lineBP = "Bacterial Endpoint"),
                                values = c(lineCR = "#333333", lineBP = "#996633")) +
    
    ggplot2::guides(color = ggplot2::guide_legend(title = '', order = 2)) +
    
    ggplot2::facet_grid(~ .data$Station_Number) + 
    
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = colorspace::lighten("#323D5E", amount = 0.0) , color = NA),
                   strip.text = ggplot2::element_text(face = 'bold', color = 'white', size = 10),
                   panel.border = ggplot2::element_rect(linewidth = 2),
                   panel.background = ggplot2::element_rect(fill = NA),
                   title = ggplot2::element_text(face = 'bold'),
                   plot.subtitle = ggplot2::element_text(face = 'plain'),
                   legend.title = ggplot2::element_text(face = 'bold', size = 10),
                   legend.text = ggplot2::element_text(size = 9),
                   axis.title = ggplot2::element_text(face = 'bold',size = 10),
                   axis.text = ggplot2::element_text(size = 10)) + 
    
    ggplot2::labs(title = paste0(unique(collision_rates_plot_df$Location), '_Collision_Rates'),
                  subtitle = 'Difference in collision rates between VP and VPC samples over time\nPercentages represents decrease of collision rates in T0_sample compared to original seawater sample (WW_sample)',
                  x = 'Sampling Timepoints\n (in hours)',
                  y = 'Mean Relative Collision Rate') 
  
  plot_name <- paste0(unique(collision_rates_plot_df$Location), '_Collision_Rates')
  .GlobalEnv$plot_list[[plot_name]] <- n
}


#' @rdname vp_visuals
plot_comparison_methods <- function(vp_results){
  ## 1. Setup
  plot_data_linear <- vp_results %>%
    dplyr::filter(stringr::str_starts(.data$VP_Method, 'LM'), .data$Population == 'c_Viruses')
  
  plot_data_VIPCAL <- vp_results %>%
    dplyr::filter(stringr::str_starts(.data$VP_Method, "VPCL"), .data$Population == 'c_Viruses')
  
  plot_data_LM_vs_VPCL <- vp_results %>%
    dplyr::filter(.data$VP_Method %in% c('LM_AR_DIFF','VPCL_AR_DIFF', 'VPCL_AR_DIFF_LMER_SE'), 
                  .data$Population == 'c_Viruses')
  
  plot_data_VPCL_vs_VPCL_SE <- vp_results %>%
    dplyr::filter(.data$VP_Method %in% c('VPCL_AR_DIFF', 'VPCL_AR_DIFF_LMER_SE'), 
                  .data$Population == 'c_Viruses')
  
  plot_data_all_methods <- vp_results %>%
    dplyr::filter(.data$Population == 'c_Viruses')
  
  ## 2. Make plot
  # Linear methods
  n_linear_mean <- ggstatsplot::ggbetweenstats(data = plot_data_linear,
                                               x = 'VP_Method',
                                               y = 'VP', 
                                               type = 'nonparametric',
                                               plot.type = "violin",
                                               violin.args = list(fill = NA),
                                               boxplot.args = list(width = 0),
                                               pairswise.display = 's',
                                               pairwise.comparisons = TRUE,
                                               centrality.plotting = FALSE,
                                               ggsignif.args = list(textsize = 3),
                                               title = 'Kruskal-Wallis Test: Linear Methods Mean')
  
  n_linear_SE <- ggstatsplot::ggbetweenstats(data = plot_data_linear,
                                             x = 'VP_Method',
                                             y = 'VP_SE', 
                                             type = 'nonparametric',
                                             plot.type = "violin",
                                             violin.args = list(fill = NA),
                                             boxplot.args = list(width = 0),
                                             pairswise.display = 's',
                                             pairwise.comparisons = TRUE,
                                             centrality.plotting = FALSE,
                                             ggsignif.args = list(textsize = 3),
                                             title = 'Kruskal-Wallis Test: Linear Methods SE')
  
  n_linear <- ggstatsplot::combine_plots(list(n_linear_mean, n_linear_SE),
                                         plotgrid.args = list(nrow = 2),
                                         annotation.args = list(
                                           title = "Comparison of viral production calculation",
                                           subtitle = "Population: c_Viruses; Calculation method: all linear regression variants")) +
    
    ggplot2::theme(axis.title = ggplot2::element_text(face = 'bold'),
                   plot.title = ggplot2::element_text(face = 'bold'),
                   plot.subtitle = ggplot2::element_text(face = 'plain'))
  
  plot_name <- paste0(unique(plot_data_linear$Location), '_Comparison_Linear_Methods')
  .GlobalEnv$plot_list[[plot_name]] <- n_linear
  
  # VIPCAL methods
  n_VIPCAL_mean <- ggstatsplot::ggbetweenstats(data = plot_data_VIPCAL,
                                               x = 'VP_Method',
                                               y = 'VP', 
                                               type = 'nonparametric',
                                               plot.type = "violin",
                                               violin.args = list(fill = NA),
                                               boxplot.args = list(width = 0),
                                               pairswise.display = 's',
                                               pairwise.comparisons = TRUE,
                                               centrality.plotting = FALSE,
                                               ggsignif.args = list(textsize = 3),
                                               title = 'Kruskal-Wallis Test: VIPCAL Methods Mean')
  
  n_VIPCAL_SE <- ggstatsplot::ggbetweenstats(data = plot_data_VIPCAL,
                                             x = 'VP_Method',
                                             y = 'VP_SE', 
                                             type = 'nonparametric',
                                             plot.type = "violin",
                                             violin.args = list(fill = NA),
                                             boxplot.args = list(width = 0),
                                             pairswise.display = 's',
                                             pairwise.comparisons = TRUE,
                                             centrality.plotting = FALSE,
                                             ggsignif.args = list(textsize = 3),
                                             title = 'Kruskal-Wallis Test: VIPCAL Methods SE')
  
  n_VIPCAL <- ggstatsplot::combine_plots(list(n_VIPCAL_mean, n_VIPCAL_SE),
                                         plotgrid.args = list(nrow = 2),
                                         annotation.args = list(
                                           title = "Comparison of viral production calculation",
                                           subtitle = "Population: c_Viruses; Calculation method: all VIPCAL variants")) +
    
    ggplot2::theme(axis.title = ggplot2::element_text(face = 'bold'),
                   plot.title = ggplot2::element_text(face = 'bold'),
                   plot.subtitle = ggplot2::element_text(face = 'plain'))
  
  plot_name <- paste0(unique(plot_data_VIPCAL$Location), '_Comparison_VIPCAL_Methods')
  .GlobalEnv$plot_list[[plot_name]] <- n_VIPCAL
  
  # LM vs VIPCAL vs VIPCAL_SE
  n_LM_vs_VIPCAL <- ggstatsplot::ggbetweenstats(data = plot_data_LM_vs_VPCL,
                                                x = 'VP_Method',
                                                y = 'VP', 
                                                type = 'parametric',
                                                plot.type = "violin",
                                                violin.args = list(fill = NA),
                                                pairswise.display = 's',
                                                pairwise.comparisons = TRUE,
                                                centrality.plotting = TRUE,
                                                ggsignif.args = list(textsize = 3),
                                                title = 'Welch test: LM vs VIPCAL vs VIPCAL_SE')
  
  n_VIPCAL_vs_VIPCAL_SE <- ggstatsplot::ggbetweenstats(data = plot_data_VPCL_vs_VPCL_SE,
                                                       x = 'VP_Method',
                                                       y = 'VP', 
                                                       type = 'parametric',
                                                       plot.type = "violin",
                                                       violin.args = list(fill = NA),
                                                       pairswise.display = 's',
                                                       pairwise.comparisons = TRUE,
                                                       centrality.plotting = TRUE,
                                                       ggsignif.args = list(textsize = 3),
                                                       title = 'Welch test: VIPCAL vs VIPCAL_SE')
  
  n_compared <- ggstatsplot::combine_plots(list(n_LM_vs_VIPCAL, n_VIPCAL_vs_VIPCAL_SE),
                                           plotgrid.args = list(nrow = 2),
                                           annotation.args = list(
                                             title = "Comparison of LM, VIPCAL and VIPCAL_SE",
                                             subtitle = "Population: c_Viruses")) +
    
    ggplot2::theme(axis.title = ggplot2::element_text(face = 'bold'),
                   plot.title = ggplot2::element_text(face = 'bold'),
                   plot.subtitle = ggplot2::element_text(face = 'plain'))
  
  plot_name <- paste0(unique(plot_data_LM_vs_VPCL$Location), '_Comparison_LM_VIPCAL_VIPCAL_SE')
  .GlobalEnv$plot_list[[plot_name]] <- n_compared
  
  # All methods
  n_all <- ggplot2::ggplot(data = plot_data_all_methods, 
                           ggplot2::aes(x = .data$VP_Method, y = .data$VP, fill = .data$VP_Method)) + 
    ggplot2::geom_violin() + 
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = 0.2), 
                        size = 1.5, shape = 16, alpha = 0.6) + 
    ggplot2::geom_hline(yintercept = 0) + 
    
    ggplot2::scale_fill_brewer(palette = 'Spectral') + 
    
    ggplot2::labs(x = 'VP calculation method',
                  y = 'Viral Production',
                  title = 'Comparison of viral production calculation methods',
                  subtitle = 'Population: c_Viruses') + 
    
    ggplot2::theme_classic() + 
    ggplot2::theme(axis.title = ggplot2::element_text(face = 'bold'),
                   title = ggplot2::element_text(face = 'bold'),
                   plot.subtitle = ggplot2::element_text(face = 'plain'),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  plot_name <- paste0(unique(plot_data_all_methods$Location), '_Comparison_ALL')
  .GlobalEnv$plot_list[[plot_name]] <- n_all
}


#' @rdname vp_visuals
plot_percentage_cells <- function(analyzed_vp_results_bacterial_endpoint){
  ## 1. Setup
  plot_data <- analyzed_vp_results_bacterial_endpoint %>%
    dplyr::filter(.data$Sample_Type != 'VPC', .data$VP_Method == 'VPCL_AR_DIFF_LMER_SE', 
                  .data$Population == 'c_Viruses') %>%
    dplyr::select('Location', 'Station_Number', 'Time_Range', 'Population', 'Sample_Type', 
                  tidyr::starts_with('P_Cells_')) %>%
    dplyr::group_by(.data$Station_Number, .data$Sample_Type) %>%
    dplyr::mutate(Timepoint = as.numeric(gsub("[^0-9.]+", "", .data$Time_Range))) %>%
    tidyr::pivot_longer(cols = tidyr::starts_with('P_Cells_'), 
                        names_to = 'Burst_Size', 
                        values_to = 'P_Cells') %>%
    dplyr::mutate(Burst_Size = substring(.data$Burst_Size, nchar(.data$Burst_Size) - 4))
  
  ## 2. Make plot
  n <- ggplot2::ggplot(data = plot_data) + 
    ggplot2::geom_col(mapping = ggplot2::aes(x = .data$Burst_Size, y = .data$P_Cells, fill = .data$Sample_Type), 
                      position = 'dodge') + 
    ggplot2::geom_text(data = plot_data, ggplot2::aes(x = 'BS_25', y = 100),
                       label = paste0('Timepoint of the assay\n(Bacterial Endpoint): T0_T', plot_data$Timepoint),
                       size = 3, color = 'black', show.legend = F) +
    
    ggplot2::scale_fill_manual(name = 'Percentage cells',
                               labels = c(Diff = 'Lysogenic cells', VP = 'Lytically infected cells'),
                               values = c(Diff = "#66CC00", VP = "#CCCC66")) + 
    
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE, order = 1)) +
    
    ggplot2::facet_grid(~ .data$Station_Number) + 
    
    ggplot2::labs(x = 'Burst_Size', 
                  y = 'Percentage of cells', 
                  title = 'Percentage of lytically infected and lysogenic cells for different burst sizes',
                  subtitle = 'Population: c_Viruses; Calculation method: VPCL_AR_DIFF_LMER_SE; Bacterial endpoint taken into account') +
    
    ggplot2::theme_bw() + 
    ggplot2::theme(strip.background = ggplot2::element_rect(color = 'black', fill = '#999999'),
                   axis.title = ggplot2::element_text(face = 'bold'),
                   title = ggplot2::element_text(face = 'bold'),
                   plot.subtitle = ggplot2::element_text(face = 'plain'),
                   legend.position = "right")
  
  plot_name <- paste0(unique(plot_data$Location), '_Percentage_Cells')
  .GlobalEnv$plot_list[[plot_name]] <- n
}


#' @rdname vp_visuals
plot_nutrient_release <- function(analyzed_vp_results_T0_T24){
  ## 1. Setup
  plot_data <- analyzed_vp_results_T0_T24 %>%
    dplyr::filter(.data$Population == 'c_Viruses', .data$Sample_Type == 'VP', 
                  .data$VP_Method == 'VPCL_AR_DIFF_LMER_SE') %>%
    dplyr::select('Location', 'Station_Number', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'VP_Method', tidyr::matches('Total_DO')) %>%
    tidyr::pivot_longer(cols = tidyr::matches('Total_DO'), 
                        names_to = 'Nutrient_per_BS', 
                        values_to = 'Nutrient_release') %>%
    dplyr::mutate(Burst_Size = substring(.data$Nutrient_per_BS, nchar(.data$Nutrient_per_BS) - 1),
           Nutrient = gsub(".*DO(.).*", "\\1", .data$Nutrient_per_BS))
  
  ## 2. Make plot
  n <- ggplot2::ggplot(data = plot_data) + 
    ggplot2::geom_col(mapping = ggplot2::aes(x = .data$Nutrient_release, y = .data$Burst_Size, fill = .data$Nutrient), 
                      position = 'dodge', orientation = 'y') + 
    
    ggplot2::scale_fill_manual(name = 'Type of nutrient',
                               values = c('#CC6666', '#339900', '#3399CC')) +
    ggplot2::scale_x_continuous() + 
    
    ggplot2::facet_grid(.data$Station_Number ~ .) + 
    
    ggplot2::labs(title = 'Total nutrient release per burst size',
                  subtitle = 'Population: c_Viruses; Sample_Type: VP; Calculation method: VPCL_AR_DIFF_LMER_SE; Time of assay: T0_T24',
                  x = 'Total nutrient release',
                  y = 'Burst size') + 
    
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title = ggplot2::element_text(face = 'bold'),
                   title = ggplot2::element_text(face = 'bold'),
                   plot.subtitle = ggplot2::element_text(face = 'plain'))
  
  plot_name <- paste0(unique(plot_data$Location), '_Total_Nutrient_Release')
  .GlobalEnv$plot_list[[plot_name]] <- n
}

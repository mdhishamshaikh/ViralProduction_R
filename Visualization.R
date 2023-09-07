###---Visualization---###

## 1. Source necessary functions and create folder to store all figures
source('viral_production_step2_source.R') # Some of the functions like dataframe functions and bacterial endpoint are necessary

folder_name <- 'Figures'

if (!file.exists(folder_name)){ # Check if folder not exists already
  dir.create(folder_name)
}

# Importing data from step 1, step 2 and the analyzing step
data <- read.csv('NJ2020.csv') %>%
  rename(Station_Number = Expt_No)

df_abundance <- read.csv('NJ2020_abundance.csv') %>% # Consist of the abundances of all populations in original seawater sample for each experiment
  rename(Station_Number = Expt_No)

vp_calc_NJ2020 <- read.csv(paste0(output_dir, '/vp_calc_ALL.csv'))
vp_calc_NJ2020_analyzed <- read.csv(paste0(output_dir, '/vp_calc_ANALYZED.csv'))

## 2. Visualization
# 2.1 An overview of the bacterial and viral counts for lytic and lysogenic inductions over time
# Data consist of the output csv.file of step 1, which is the input file for step 2
overview_plot_counts_over_time <- function(data){
  
  # Add tag column to data that represents the unique ID (Location_Station_Depth)
  data_tag <- data %>%
    unite(all_of(c('Location', 'Station_Number', 'Depth')), col = 'tag', remove = F)
  
  
  ## !!! FIX subgroup thing => not in dataframe anymore !!!
  
  
  # Create dataframe for plot: averaging the sample replicates, calculate the differences by subtraction and add timepoints
  df_plot_cot <- df_AVG(data) %>%
    mutate(Sample_Type = factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
           Subcrobe = ifelse(Subgroup == 'Parent', 'Total', as.character(Microbe)),
           Subcrobe = factor(Subcrobe, levels= c('Total', 'Bacteria', 'Viruses'))) %>% # Added a column for facet_grid later on => Subgroup and Microbe to one
    unite(all_of(c('Location', 'Station_Number', 'Depth')), col = 'tag', remove = F) # Add unique tag
  
  # Filter data based on Location, Station_Number and Depth and make plot
  for (combi_tag in unique(df_plot_cot$tag)){
    
    # Dataframe for ggplot
    df_plot_cot_sub <- df_plot_cot %>%
      filter(tag == combi_tag) %>%
      mutate(Time_Time = factor(Time_Time, levels = unique(Time_Time))) # Different Time_Ranges over experiments, so factor on subset and not on df_plot_cot
    
    # Original data for bacterial endpoint
    data_bp <- data_tag %>%
      filter(tag == combi_tag)
    
    # Create ggplot object
    n <- ggplot(df_plot_cot_sub, aes(x = Timepoint, y = Mean, color = Population, shape = Population)) +
      
      # Adding layers to represent the data
      geom_point(size = 1.5) +
      geom_line() +
      geom_smooth(linewidth = 1.0, method = 'lm', se = F) +
      geom_hline(yintercept = 0, color = "#000000", size = 0.3, linetype = 'dashed') +
      geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.5, size = 0.5) +
      
      # Divide plot in subplots for each population: facet_grid(rows ~ columns)
      facet_grid(Subcrobe + Sample_Type ~ Time_Time) +
      
      # Scaling of axes, colors and shapes
      scale_y_continuous(labels = unit_format(unit = NULL, scale = 1e-6)) + 
      scale_x_continuous(breaks = unique(df_plot_cot_sub$Timepoint)) + 
      scale_shape_manual(name = 'Populations', 
                         labels = c("Total Bacteria", "HNA Bacteria", "LNA Bacteria",
                                    "V1 Viruses", "V2 Viruses", "V3 Viruses", "Total Viruses"),
                         values = c(16,17,16,17,15,16,17)) +
      scale_color_manual(name = 'Populations',
                         labels = c("Total Bacteria", "HNA Bacteria", "LNA Bacteria",
                                    "V1 Viruses", "V2 Viruses", "V3 Viruses", "Total Viruses"),
                         values = c(c_Bacteria = "#FF0000", c_HNA = "#0000FF", c_LNA = "#0099CC",
                                    c_Viruses = "#FF0066", c_V1 = "#666666", c_V2 = "#666633",
                                    c_V3 = "#669900")) + 
      
      # Labels
      labs(title = paste0(paste0(unique(df_plot_cot_sub$Location), '_Station_', unique(df_plot_cot_sub$Station_Number), '_Depth_', unique(df_plot_cot_sub$Depth))," - Overview"),
           subtitle = 'Bacterial and Viral counts for Lytic and Lysogenic inductions',
           x = 'Sampling Timepoints\n (in hours)',
           y = 'FCM Counts\n (in millions)') + 
      
      # Theme
      theme_bw() + 
      theme(strip.background = element_rect(color = 'black', fill = 'white'),
            axis.title = element_text(face = 'bold'),
            title = element_text(face = 'bold'),
            plot.subtitle = element_text(face = 'plain'),
            legend.position = "bottom") +
      
      # Guide
      guides(color = guide_legend(nrow = 2, byrow = FALSE))
    
    # Determine bacterial endpoint
    bp_cot <- bacterial_endpoint(data_bp, visual = T) - 1 # We want to highlight the timepoint before the one were GT of bacterie is lower then 24
    
    # Highlight the bacterial endpoint facet grid (= timepoint)
    n_gtable <- ggplot_gtable(ggplot_build(n)) # ggplot_build makes object that can be rendered from plot_object, ggplot_gtable creates structured grid representation of plot
    strip_index <- which(grepl('strip-', n_gtable$layout$name)) # Find indices of the facet strip elements
    
    fills <- c(rep("white", length(strip_index))) # Vector of fill colors
    fills[bp_cot] <- "#FF6699" # Change color of bacterial endpoint
    
    # Iterate through facet strip element in gtable and modify the fill color of the background rectangle (rect grob)
    color <- 1
    for (i in strip_index) { # iterate thorugh facet strip elements
      j <- which(grepl('rect', n_gtable$grobs[[i]]$grobs[[1]]$childrenOrder)) # index that corresponds to rect grob within current facet strip
      n_gtable$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[color] # Modify the color
      color <- color + 1 # Go to next color
    }
    
    # Save plot as svg file
    ggsave(paste0('Figures/', paste0(unique(df_plot_cot_sub$Location), '_Station_', unique(df_plot_cot_sub$Station_Number), '_Depth_', unique(df_plot_cot_sub$Depth)), '_Overview.svg'), plot = n_gtable, width = 10, height = 10)
  }
}

overview_plot_counts_over_time(data)

# 2.2 Difference in collision rates between VP and VPC samples over time
# To determine the collision rates: output csv.file from Step 1 + raw abundances in seawater sample (WW)
# Custom Key Glyph
draw_key_cust <- function(data, params, size) {
  if (data$colour == "#996633") {
    draw_key_vpath(data, params, size)
  } else {
    draw_key_path(data, params, size)
  }
}

collision_rates_plot <- function(data, abundance){
  
  # Select correct data
  df_contact_rates <- data %>%
    select(all_of(c("Location", "Station_Number", "Sample_Type", "Timepoint", "Replicate", "c_Bacteria", "c_Viruses")))
  
  df_abundance <- abundance %>% # Consist of the abundances of all populations in original seawater sample for each experiment
    rename(c_Bacteria = Total_Bacteria, c_Viruses = Total_Viruses) %>%
    select(all_of(c('Location', 'Station_Number', 'c_Bacteria', 'c_Viruses'))) 
  
  # Abundance file consists of just one abundance per experiment => in flow cytometry, we analysed three sample types with each three replicates
  # Expand abundance file so that dimensions match
  abundance_dims <- data.frame(expand.grid('Location' = 'NJ2020',
                                           'Station_Number' = 1:7,
                                           'Sample_Type' = c('VP', 'VPC', '0.22'),
                                           'Timepoint' = -3,
                                           'Replicate' = 1:3))
  
  # Join to one dataframe
  df_abundance <- full_join(abundance_dims, df_abundance)
  df_contact_rates <- full_join(df_abundance, df_contact_rates)
  
  # Calculate the collision rates: Bacterial_abundance * Viral_Abundance at a certain timepoint
  # The change in collision rate: CR_T / CR_T0
  df_contact_rates$BV <- df_contact_rates$c_Bacteria * df_contact_rates$c_Viruses
  cr_res <- list()
  
  for (loc in unique(df_contact_rates$Location)){
    for (station in unique(df_contact_rates$Station_Number)){
      for(sample in unique(df_contact_rates$Sample_Type)){
        for (rep in unique(df_contact_rates$Replicate)){
          
          df <- df_contact_rates %>%
            filter(Location == loc, Station_Number == station, Sample_Type == sample, Replicate == rep)
          
          for(time in unique(df$Timepoint)){
            
            CR <- (df %>% filter(Timepoint == time))$BV / (df %>% filter(Timepoint == 0))$BV
            res <- c(loc, station, sample, rep, time, CR)
            cr_res[[length(cr_res)+1]]<- res
          }
        }
      }
    }
  }
  # All results from nested list to table
  df_CR <- data.table(t(as.data.table(cr_res))) 
  colnames(df_CR) <- c("Location", "Station_Number", "Sample_Type", "Replicate", "Timepoint", "Collision_Rate")
  df_CR <- df_CR %>%
    mutate_at('Collision_Rate', as.numeric) %>%
    filter(Sample_Type != '0.22') %>%
    group_by(Location, Station_Number, Sample_Type, Timepoint) %>%
    summarise(CR_Mean = mean(Collision_Rate), CR_SE = plotrix::std.error(Collision_Rate)) %>% # Average over the replicates
    arrange('Location', 'Station_Number') 

  # Calculate the percentage of decrease of CR between seawater sample (WW) and T0
  decrease_res <- list()
  
  for (station in unique(df_CR$Station_Number)){
    for (sample in unique(df_CR$Sample_Type)){
      df <- df_CR %>%
        filter(Station_Number == station, Sample_Type == sample, Timepoint == c(-3,0))
      
      decr <- (df[df$Timepoint == '-3',]$CR_Mean - df[df$Timepoint == '0',]$CR_Mean) / df[df$Timepoint == '-3',]$CR_Mean * 100
      res <- c(station, sample, decr)
      decrease_res[[length(decrease_res)+1]]<- res
    }
  }
  decrease_df <- data.frame(t(sapply(decrease_res, c)))
  colnames(decrease_df) <- c('Station_Number', 'Sample_Type', 'Perc_Decrease')
  decrease_df <- decrease_df %>%
    pivot_wider(names_from = Sample_Type,
                values_from = Perc_Decrease) %>%
    mutate_at(c('VP', 'VPC'), as.numeric)
  
  # For plotting: original samples can be removed, start from T0
  df_CR_plot <- df_CR %>%
    filter(Timepoint != '-3')
  
  # Calculate the difference in collision rate between VP and VPC samples
  df_CR_plot_diff <- df_CR_plot %>%
    select(-CR_SE) %>%
    pivot_wider(names_from = Sample_Type,
                values_from = CR_Mean) %>%
    mutate(CR_Diff = VP - VPC)
  
  # Calculate the bacterial endpoint to show on plot
  data_bp <- data %>% 
    unite(all_of(c('Location', 'Station_Number', 'Depth')), col = 'tag', remove = F)
  
  bp_res <- list()
  
  for (combi_tag in unique(data_bp$tag)){
    data_bp_tag <- data_bp %>%
      filter(tag == combi_tag)
    
    timepoints <- unique(data_bp_tag$Timepoint)
    
    bp <- bacterial_endpoint(data_bp_tag, visual = T)
    timepoint <- timepoints[bp]
    bp_res[[length(bp_res)+1]] <- timepoint
  }
  
  # Result dataframe consists of the bacterial endpoint for every unique tag (experiment)
  bp_df <- data.frame(t(as.data.frame(bp_res)))
  colnames(bp_df) <- c('Bacterial_Endpoint')
  rownames(bp_df) <- c(1:7)
  bp_df <- bp_df %>%
    mutate(Station_Number = rownames(bp_df))
  
  # Create the main plot
  n <- ggplot() + 
    
    geom_point(data = df_CR_plot, aes(x = as.numeric(Timepoint), y = CR_Mean,
                                      color = Sample_Type, shape = Sample_Type), show.legend = T) + 
    geom_text(data = decrease_df, aes(x = 11, y = 12),
              label = paste0('VP Decrease: ', round(decrease_df$VP,2), ' %\nVPC Decrease: ', round(decrease_df$VPC,2), ' %'),
              size = 3, color = 'black', show.legend = F) +
    
    scale_color_manual(values = c('#FF0000', '#0000FF')) +
    scale_fill_manual(values = c('#FF0000', '#0000FF')) +
    scale_shape_manual(values = c(15, 17)) +
    
    guides(color = guide_legend(title = 'Treatment', order = 1),
           shape = guide_legend(title = 'Treatment', order = 1),
           fill = guide_legend(title = 'Treatment', order = 1)) +
    
    ggnewscale::new_scale_color() +
    geom_line(data = df_CR_plot_diff, aes(x = as.numeric(Timepoint), y = CR_Diff,
              color = "lineCR"), size = 1, alpha = 0.5, key_glyph = 'cust') +
    geom_vline(data = bp_df, aes(xintercept = Bacterial_Endpoint,
               color = "lineBP"), linewidth = 1.5, alpha = 0.5, key_glyph = 'cust') +
    
    scale_color_manual(labels = c(lineCR = "Difference in collision rates between VP and VPC treatment", lineBP = "Bacterial Endpoint"),
                       values = c(lineCR = "#333333", lineBP = "#996633")) +
    
    guides(color = guide_legend(title = '', order = 2)) +
    
    facet_wrap(~ Station_Number, ncol = 3) + 
    
    theme_bw() +
    theme(strip.background = element_rect(fill = lighten("#323D5E", amount = 0.0) , color = NA),
          strip.text = element_text(face = 'bold', color = 'white', size = 10),
          panel.border = element_rect(linewidth = 2),
          panel.background = element_rect(fill = NA),
          title = element_text(face = 'bold'),
          plot.subtitle = element_text(face = 'plain'),
          legend.title = element_text(face = 'bold', size = 10),
          legend.text = element_text(size = 9),
          axis.title = element_text(face = 'bold',size = 10),
          axis.text = element_text(size = 10)) + 
    
    labs(title = paste0(unique(df_CR_plot$Location), '_Collision_Rates'),
         subtitle = 'Difference in collision rates between VP and VPC samples over time\nPercentage represents decrease of collision rates in T0_sample compared to original seawater sample (WW_sample)',
         x = 'Sampling Timepoints\n (in hours)',
         y = 'Mean Relative Collision Rate') 
  
  # Save plot as svg file
  ggsave(paste0('Figures/', unique(df_CR_plot$Location), '_Collision_Rates.svg'), plot = n, width = 15, height = 8)
}

collision_rates_plot(data, df_abundance)

# 2.3 Comparison of viral production calculation by different methods
compare_variants_vp <- function(data){
  
  ## 1. Linear regression variants
  # Plot data
  plot_data_LM <- data %>%
    filter(str_starts(VP_Type, "LM"), Population == 'c_Viruses')
  
  # Mean
  n1 <- ggstatsplot::ggbetweenstats(data = plot_data_LM,
                                    x = VP_Type,
                                    y = VP, 
                                    type = 'nonparametric',
                                    plot.type = "violin",
                                    violin.args = list(fill = NA),
                                    boxplot.args = list(width = 0),
                                    pairswise.display = 's',
                                    pairwise.comparisons = TRUE,
                                    centrality.plotting = FALSE,
                                    ggsignif.args = list(textsize = 3),
                                    title = 'Kruskal-Wallis Test: LM Mean')
  
  # SE
  n2 <- ggstatsplot::ggbetweenstats(data = plot_data_LM,
                                    x = VP_Type,
                                    y = VP_SE, 
                                    type = 'nonparametric',
                                    plot.type = "violin",
                                    violin.args = list(fill = NA),
                                    boxplot.args = list(width = 0),
                                    pairswise.display = 's',
                                    pairwise.comparisons = TRUE,
                                    centrality.plotting = FALSE,
                                    ggsignif.args = list(textsize = 3),
                                    title = 'Kruskal-Wallis Test: LM SE')
  
  # Combine and save plot
  n <- ggstatsplot::combine_plots(
    list(n1, n2),
    plotgrid.args = list(nrow = 2),
    annotation.args = list(
      title = "Comparison of viral production calculation",
      subtitle = "Population: c_Viruses; Calculation method: all linear regression variants"
    )
  ) + 
    theme(axis.title = element_text(face = 'bold'),
          plot.title = element_text(face = 'bold'),
          plot.subtitle = element_text(face = 'plain'))
  
  ggsave(paste0('Figures/', unique(data$Location), '_Comparison_Linear_Methods.svg'), plot = n, width = 8, height = 10)
  
  ## 2. VIPCAL variants
  # Plot data
  plot_data_VPCL <- data %>%
    filter(str_starts(VP_Type, "VPCL"), Population == 'c_Viruses')
  
  # Mean
  m1 <- ggstatsplot::ggbetweenstats(data = plot_data_VPCL,
                                    x = VP_Type,
                                    y = VP, 
                                    type = 'nonparametric',
                                    plot.type = "violin",
                                    violin.args = list(fill = NA),
                                    boxplot.args = list(width = 0),
                                    pairswise.display = 's',
                                    pairwise.comparisons = TRUE,
                                    centrality.plotting = FALSE,
                                    ggsignif.args = list(textsize = 3),
                                    title = 'Kruskal-Wallis Test: VIPCAL Mean')
  
  # SE
  m2 <- ggstatsplot::ggbetweenstats(data = plot_data_VPCL,
                                    x = VP_Type,
                                    y = VP_SE, 
                                    type = 'nonparametric',
                                    plot.type = "violin",
                                    violin.args = list(fill = NA),
                                    boxplot.args = list(width = 0),
                                    pairswise.display = 's',
                                    pairwise.comparisons = TRUE,
                                    centrality.plotting = FALSE,
                                    ggsignif.args = list(textsize = 3),
                                    title = 'Kruskal-Wallis Test: VIPCAL SE')
  
  # Combine and save plot
  m <- ggstatsplot::combine_plots(
    list(m1, m2),
    plotgrid.args = list(nrow = 2),
    annotation.args = list(
      title = "Comparison of viral production calculation",
      subtitle = "Population: c_Viruses; Calculation method: all VIPCAL variants"
    )
  ) +
    theme(axis.title = element_text(face = 'bold'),
          plot.title = element_text(face = 'bold'),
          plot.subtitle = element_text(face = 'plain'))
  
  ggsave(paste0('Figures/', unique(data$Location), '_Comparison_VIPCAL_Methods.svg'), plot = m, width = 10, height = 10)
  
  ## 3. VIPCAL/VIPCAL_SE/LM
  # Plot data
  plot_data1 <- data %>%
    filter(VP_Type %in% c('VPCL_AR_DIFF', 'VPCL_AR_DIFF_LMER_SE'), Population == 'c_Viruses') %>%
    mutate(VP = VP / 1e6)
  
  plot_data2 <- data %>%
    filter(VP_Type %in% c('VPCL_AR_DIFF', 'VPCL_AR_DIFF_LMER_SE', 'LM_AR_DIFF'), Population == 'c_Viruses') %>%
    mutate(VP = VP / 1e6)
  
  # ggplot object
  x1 <- ggstatsplot::ggbetweenstats(data = plot_data1,
                                    x = VP_Type,
                                    y = VP, 
                                    type = 'parametric',
                                    plot.type = "violin",
                                    violin.args = list(fill = NA),
                                    pairswise.display = 's',
                                    pairwise.comparisons = TRUE,
                                    centrality.plotting = TRUE,
                                    ggsignif.args = list(textsize = 3),
                                    title = 'Welch test: VIPCAL vs VIPCAL_SE')
  
  x2 <- ggstatsplot::ggbetweenstats(data = plot_data2,
                                    x = VP_Type,
                                    y = VP, 
                                    type = 'parametric',
                                    plot.type = "violin",
                                    violin.args = list(fill = NA),
                                    pairswise.display = 's',
                                    pairwise.comparisons = TRUE,
                                    centrality.plotting = TRUE,
                                    ggsignif.args = list(textsize = 3),
                                    title = 'Welch test: LM vs VIPCAL vs VIPCAL_SE')
  
  # Combine and save plot
  x <- ggstatsplot::combine_plots(
    list(x1, x2),
    plotgrid.args = list(nrow = 2),
    annotation.args = list(
      title = "Comparison of LM, VIPCAL and VIPCAL_SE",
      subtitle = "Population: c_Viruses"
    )
  ) +
    theme(axis.title = element_text(face = 'bold'),
          plot.title = element_text(face = 'bold'),
          plot.subtitle = element_text(face = 'plain'))
  
  ggsave(paste0('Figures/', unique(data$Location), '_Comparison_LM_VIPCAL_VIPCAL_SE.svg'), plot = x, width = 8, height = 10)
  
  ## 4. Overview all methods
  # Plot data
  plot_data_all <- data %>%
    filter(Population == 'c_Viruses')
  
  # ggplot object
  all <- ggplot(data = plot_data_all, aes(x = VP_Type, y = VP, fill = VP_Type)) + 
    geom_violin() + 
    geom_point(position = position_jitterdodge(jitter.width = 0.2), size = 1.5, shape = 16, alpha = 0.6) + 
    geom_hline(yintercept = 0) + 
    scale_fill_brewer(palette = 'Spectral') + 
    labs(x = 'VP calculation method',
         y = 'Viral Production',
         title = 'Comparison of viral production calculation methods',
         subtitle = 'Population: c_Viruses') + 
    theme_classic() + 
    theme(axis.title = element_text(face = 'bold'),
          title = element_text(face = 'bold'),
          plot.subtitle = element_text(face = 'plain'),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  # Save plot as svg
  ggsave(paste0('Figures/', unique(data$Location), '_Comparison_ALL.svg'), plot = all, width = 15, height = 10)
}

compare_variants_vp(vp_calc_NJ2020)

# 2.4 Percentage lytically infected and lysogenic cells over time
percentage_cells <- function(data){
  
  # Select data
  plot_data <- data %>%
    filter(Sample_Type != 'VPC', VP_Type == 'VPCL_AR_DIFF_LMER_SE', Population == 'c_Viruses') %>%
    select('Location', 'Station_Number', 'Time_Range', 'Population', 'Sample_Type', starts_with('P_Cells_')) %>%
    group_by(Station_Number, Sample_Type) %>%
    mutate(Timepoint = as.numeric(gsub("[^0-9.]+", "", Time_Range))) %>%
    pivot_longer(cols = starts_with('P_Cells_'), names_to = 'Burst_Size', values_to = 'P_Cells') %>%
    mutate(Burst_Size = substring(Burst_Size, nchar(Burst_Size) - 4))
  
  
  # ggplot object
  n <- ggplot(data = plot_data) + 
    geom_col(mapping = aes(x = Burst_Size, y = P_Cells, fill = Sample_Type), position = 'dodge') + 
    geom_text(data = plot_data, aes(x = 'BS_25', y = 100),
              label = paste0('Timepoint of the assay\n(Bacterial Endpoint): T0_T', plot_data$Timepoint),
              size = 3, color = 'black', show.legend = F) +
    scale_fill_manual(name = 'Percentage cells',
                      labels = c(Diff = 'Lysogenic cells', VP = 'Lytically infected cells'),
                      values = c(Diff = "#66CC00", VP = "#CCCC66")) + 
    
    guides(fill = guide_legend(nrow = 2, byrow = TRUE, order = 1)) +
    
    facet_wrap(~Station_Number, ncol = 3) + 
    labs(x = 'Burst_Size', 
         y = 'Percentage of cells', 
         title = 'Percentage of lytically infected and lysogenic cells for different burst sizes',
         subtitle = 'Population: c_Viruses; Calculation method: VPCL_AR_DIFF_LMER_SE; Bacterial endpoint taken into account') +
    theme_bw() + 
    theme(strip.background = element_rect(color = 'black', fill = '#999999'),
          axis.title = element_text(face = 'bold'),
          title = element_text(face = 'bold'),
          plot.subtitle = element_text(face = 'plain'),
          legend.position = "right")
  
  # Save plot
  ggsave(paste0('Figures/', unique(plot_data$Location), '_Percentage_Cells.svg'), plot = n, width = 15, height = 8)
}

# Use the bacterial endpoint csv instead of all the analyzed results
vp_calc_NJ2020_bp <- read.csv(paste0(output_dir, '/vp_calc_BP.csv'))
vp_calc_NJ2020_bp_analyzed <- analyze_vpres(vp_calc_NJ2020_bp, data, df_abundance, write_output = F)

percentage_cells(vp_calc_NJ2020_bp_analyzed)

# 2.5 Total nutrient release
nutrient_release <- function(data){
  
  # Select data
  plot_data <- data %>%
    filter(Population == 'c_Viruses', Sample_Type == 'VP', VP_Type == 'VPCL_AR_DIFF_LMER_SE', Time_Range == 'T0_T24') %>%
    select('Location', 'Station_Number', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'VP_Type', matches('Total_DO')) %>%
    pivot_longer(cols = matches('Total_DO'), names_to = 'Nutrient_per_BS', values_to = 'Nutrient_release') %>%
    mutate(Burst_Size = substring(Nutrient_per_BS, nchar(Nutrient_per_BS) - 1),
           Nutrient = gsub(".*DO(.).*", "\\1", Nutrient_per_BS),
           Nutrient_release = Nutrient_release * 1e11)
  
  # ggplot object
  n <- ggplot(data = plot_data) + 
    geom_col(mapping = aes(x = Nutrient_release, y = Burst_Size, fill = Nutrient), 
             position = 'dodge', orientation = 'y') + 
    scale_fill_manual(name = 'Type of nutrient',
                      values = c('#CC6666', '#339900', '#3399CC'))+
    scale_x_continuous() + 
    facet_grid(Station_Number ~ .) + 
    labs(title = 'Total nutrient release per burst size',
         subtitle = 'Population: c_Viruses; Sample_Type: VP; Calculation method: VPCL_AR_DIFF_LMER_SE; Time of assay: T0_T24',
         x = 'Total nutrient release (x10e-11 g nutrient/mLh)',
         y = 'Burst size') + 
    theme_bw() +
    theme(axis.title = element_text(face = 'bold'),
          title = element_text(face = 'bold'),
          plot.subtitle = element_text(face = 'plain'))
  
  # Save plot
  ggsave(paste0('Figures/', unique(plot_data$Location), '_Total_Nutrient_Release.svg'), plot = n, width = 10, height = 10)
}
  
nutrient_release(vp_calc_NJ2020_analyzed)

# Test
## Compare VPCL_AR_DIFF vs VPCL_AR_DIFF_LMER_SE viral production
test1 <- function(vp_results){
  plot_DF <- vp_results %>%
    filter(VP_Method %in% c('VPCL_AR_DIFF', 'VPCL_AR_DIFF_LMER_SE')) %>%
    select(-all_of(c('abs_VP', 'VP_SE', 'VP_R_Squared'))) %>%
    mutate(VP = VP / 1e6) %>%
    pivot_wider(names_from = 'VP_Method', values_from = 'VP')
  
  n <- ggplot(data = plot_DF, aes(x = VPCL_AR_DIFF, y = VPCL_AR_DIFF_LMER_SE,
                                  color = Sample_Type, shape = Sample_Type,
                                  fill = Sample_Type)) + 
    geom_point() + 
    geom_abline(intercept = 0, slope = 1) + 
    
    ggplot2::labs(title = 'Comparison VIPCAL, VIPCAL-SE vp_values') + 
    
    ggplot2::theme_bw() + 
    ggplot2::theme(strip.background = ggplot2::element_rect(color = 'black', fill = 'white'),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.title = ggplot2::element_text(face = 'plain'),
                   title = ggplot2::element_text(face = 'bold'))
  
  return(n)
}

test2 <- function(vp_results){
  plot_DF <- vp_results %>%
    filter(VP_Method %in% c('VPCL_AR_DIFF', 'VPCL_AR_DIFF_LMER_SE')) %>%
    select(-all_of(c('abs_VP', 'VP_SE', 'VP_R_Squared'))) 
  
  n <- ggplot(data = plot_DF, aes(x = VP, y = VP_Method)) + 
    geom_point(position = 'jitter') +
    
    ggplot2::labs(title = 'Comparison VIPCAL, VIPCAL-SE ROGME distribution') + 
    
    ggplot2::theme_bw() + 
    ggplot2::theme(strip.background = ggplot2::element_rect(color = 'black', fill = 'white'),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.title = ggplot2::element_text(face = 'plain'),
                   title = ggplot2::element_text(face = 'bold'))
  
  return(n)
}

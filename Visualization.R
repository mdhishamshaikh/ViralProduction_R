###---Visualization---###

## 1. Source necessary functions and create folder to store all figures
source('viral_production_step2_source.R') # Some of the functions like dataframe functions and bacterial endpoint are necessary

folder_name <- 'Figures'

if (!file.exists(folder_name)){ # Check if folder not exists already
  dir.create(folder_name)
}

## 2. Visualization
# 2.1 An overview of the bacterial and viral counts for lytic and lysogenic inductions over time
# Data consist of the output csv.file of step 1, which is the input file for step 2
overview_plot_counts_over_time <- function(data){
  
  # Add tag column to data that represents the unique ID (Location_Station_Depth)
  data_tag <- read.csv(data) %>%
    rename(Station_Number = Expt_No) %>%
    unite(c('Location', 'Station_Number', 'Depth'), col = 'tag', remove = F)
  
  # Create dataframe for plot: averaging the sample replicates, calculate the differences by subtraction and add timepoints
  df_plot_cot <- df_AVG(data_tag) %>%
    mutate(Sample_Type = factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
           Time_Time = factor(Time_Time, levels = unique(Time_Time)),
           Subcrobe = ifelse(Subgroup == 'Parent', 'Total', as.character(Microbe)),
           Subcrobe = factor(Subcrobe, levels= c('Total', 'Bacteria', 'Viruses'))) %>% # Added a column for facet_grid later on => Subgroup and Microbe to one
    unite(c('Location', 'Station_Number', 'Depth'), col = 'tag', remove = F) # df_AVG() doesn't take the tag column with it
  
  # Filter data based on Location, Station_Number and Depth and make plot
  for (combi_tag in unique(data_tag$tag)){
    
    # Dataframe for ggplot
    df_plot_cot_sub <- df_plot_cot %>%
      filter(tag == combi_tag)
    
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
      scale_y_continuous(labels = unit_format(unit = NULL, scale = 1e-6),
                         breaks = seq(-4e+6,16e+6, 4e+6),
                         limits = c(-4e+6, 16e+6)) + 
      scale_x_continuous(breaks = unique(df_plot_cot$Timepoint)) + 
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
      theme(strip.text = element_text(face = "bold"),
            strip.background = element_rect(color = 'black', fill = 'white'),
            axis.title = element_text(face = 'bold'),
            title = element_text(face = 'bold'),
            legend.position = "bottom") +
      
      # Guide
      guides(color = guide_legend(nrow = 2, byrow = TRUE))
    
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
    
    # Draw and return plot
    grid::grid.draw(n_gtable)
    
    # Save plot as svg file
    ggsave(paste0('Figures/', paste0(unique(df_plot_cot_sub$Location), '_Station_', unique(df_plot_cot_sub$Station_Number), '_Depth_', unique(df_plot_cot_sub$Depth)), '_Overview.svg'), plot = n_gtable, width = 10, height = 10)
  }
}

#overview_plot_counts_over_time('NJ1.csv')
overview_plot_counts_over_time('NJ2020.csv')

# 3.2 Difference in collision rates between VP and VPC samples over time
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
  
  # Read in both files
  df_contact_rates <- read.csv(data) %>%
    rename(Station_Number = Expt_No) %>%
    select(c("Location", "Station_Number", "Sample_Type", "Timepoint", "Replicate", "c_Bacteria", "c_Viruses", "VBR"))
  
  df_abundance <- read.csv(abundance) %>% # Consist of the abundances of all populations in original seawater sample for each experiment
    rename(Station_Number = Expt_No, c_Bacteria = Total_Bacteria, c_Viruses = Total_Viruses) %>%
    select(c('Location', 'Station_Number', 'c_Bacteria', 'c_Viruses', 'VBR')) 
  
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
  
  # Calculate the collision rates
  # Determine the collision rate: Bacterial_abundance * Viral_Abundance at a certain timepoint
  # The change in collision rate: CR_T / CR_T0
  df_contact_rates$BV <- df_contact_rates$c_Bacteria * df_contact_rates$c_Viruses
  cr_res <- list()
  
  for (loc in unique(df_contact_rates$Location)){
    for (station in unique(df_contact_rates$Station_Number)){
      for(sample in unique(df_contact_rates$Sample_Type)){
        for (rep in unique(df_contact_rates$Replicate)){
          
          df <- df_contact_rates %>%
            filter(Location == loc, Station_Number == station, Sample_Type == sample, Replicate == rep)
          
          for(time in unique(df_contact_rates$Timepoint)){
            
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
    arrange('Location', 'Station_Number',
            factor(Timepoint, levels = c("-3", "0", "3", "6", "9", "12", "24"))) 
  
  
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
  
  # For plotting, start from T0
  df_CR_plot <- df_CR %>%
    filter(Timepoint != '-3')
  
  # Calculate the differnce in collision rate between VP and VPC samples
  df_CR_plot_diff <- df_CR_plot %>%
    select(-CR_SE) %>%
    pivot_wider(names_from = Sample_Type,
                values_from = CR_Mean) %>%
    mutate(CR_Diff = VP - VPC)
  
  # Calculate the bacterial endpoint to show on plot
  data_bp <- read.csv(data) %>%
    rename(Station_Number = Expt_No) %>% 
    unite(c('Location', 'Station_Number', 'Depth'), col = 'tag', remove = F)
  bp_res <- list()
  timepoints <- unique(data_bp$Timepoint)
  
  for (combi_tag in unique(data_bp$tag)){
    data_bp_tag <- data_bp %>%
      filter(tag == combi_tag)
    
    bp <- bacterial_endpoint(data_bp_tag, visual = T)
    timepoint <- timepoints[bp]
    bp_res[[length(bp_res)+1]] <- timepoint
  }
  bp_df <- data.frame(t(as.data.frame(bp_res)))
  colnames(bp_df) <- c('Bacterial_Endpoint')
  rownames(bp_df) <- c(1:7)
  bp_df <- bp_df %>%
    mutate(Station_Number = rownames(bp_df))
  
  # Create the main plot
  n <- ggplot() + 
    
    geom_point(data = df_CR_plot, aes(x = as.numeric(Timepoint), y = CR_Mean,
                                      color = Sample_Type, shape = Sample_Type), show.legend = T) + 
    geom_text(data = decrease_df, aes(x = 12, y = 12),
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
    
    scale_color_manual(name = 'Lines',
                       labels = c(lineCR = "Difference in collision rates between VP and VPC treatment", lineBP = "Bacterial Endpoint"),
                       values = c(lineCR = "#333333", lineBP = "#996633"),
                       guide = guide_legend(order = 2)) +
    
    facet_wrap(~ Station_Number, ncol = 3) + 
    
    theme_bw() +
    theme(strip.background = element_rect(fill = lighten("#323D5E", amount = 0.0) , color = NA),
          strip.text = element_text(face = 'bold', color = 'white', size = 10),
          panel.border = element_rect(linewidth = 2),
          panel.background = element_rect(fill = NA),
          title = element_text(face = 'bold'),
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

collision_rates_plot('NJ2020.csv', 'NJ2020_abundance.csv')











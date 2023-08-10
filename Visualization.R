###---Visualization---###

## 1. Source necessary functions and create folder to store all figures
source('viral_production_step2_source.R') # Some of the functions like dataframe functions and bacterial endpoint are necessary

folder_name <- 'Figures'

if (!file.exists(folder_name)){ # Check if folder not exists already
  dir.create(folder_name)
}

## 2. Importing data from flow cytometer
data <- read.csv('NJ1.csv')
names(data)[names(data) == 'Expt_No'] <- 'Station_Number'

data_all <- read.csv('NJ2020.csv')
names(data_all)[names(data_all) == 'Expt_No'] <- 'Station_Number'

## 3. Visualization
# 3.1.  An overview of the bacterial and viral counts for lytic and lysogenic inductions over time
# This is the output of step 1 and the input for step 2
overview_plot_counts_over_time <- function(data){
  
  # Create dataframe for plot: averaging the sample replicates, calculate the differences by subtraction and add timepoints
  df_plot_cot <- df_AVG(data) %>%
    mutate(Sample_Type = factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
           Time_Time = factor(Time_Time, levels = unique(Time_Time)),
           Subcrobe = ifelse(Subgroup == 'Parent', 'Total', as.character(Microbe)),
           Subcrobe = factor(Subcrobe, levels= c('Total', 'Bacteria', 'Viruses'))) %>% # Added a column for facet_grid later on => Subgroup and Microbe to one
    unite(c('Location', 'Station_Number', 'Depth'), col = 'tag', remove = F)
  
  # Filter data based on Location, Station_Number and Depth and make plot
  for (combi_tag in unique(df_plot_cot$tag)){
    df_plot_cot_sub <- df_plot_cot %>%
      filter(tag == combi_tag)
    
    # Create ggplot object
    n <- ggplot(df_plot_cot_sub, aes(x = Timepoint, y = Mean, color = Population, shape = Population)) +
      
      # Adding layers to represent the data
      geom_point(size = 1.5) +
      geom_line() +
      geom_smooth(size = 1.0, method = 'lm', se = F) +
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
                                    c_Viruses = "#FF0066", c_V1 = "#666666", c_V2 = "#666600",
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
    bp_cot <- bacterial_endpoint(df_plot_cot_sub, visual = T) - 1 # We want to highlight the timepoint before the one were GT of bacterie is lower then 24
    
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

overview_plot_counts_over_time(data)
overview_plot_counts_over_time(data_all)










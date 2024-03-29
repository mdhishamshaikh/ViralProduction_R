#A script to install and load required packages, and to source global functions

####1. Installing and loading libraries####

#Installing BiocManager
{
  if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
}

#List of packages to install
packages_to_load<- c("tidyverse", 
                     "flowWorkspace",
                     "scales",
                     "readxl")


#Checking if packages are already present. If absent, then installing packages from BiocManager 
for (pack in packages_to_load){
  if(!requireNamespace(pack))
    BiocManager::install(pack, force = T)
  }

#Loading libraries  
for (pack in packages_to_load){
   library(pack, character.only = T)
}



####2. Overview plot functions ####

#To create an overview table with time ranges. Needed for overview plotting
overview_df_tp_avg<- function(df, keep_0.22 = F) {
  #only works if we remove 0.22
  DF<- df[df$Sample_Type != '0.22',]
  
  DF<- DF %>%
    gather('c_Bacteria', 'c_HNA', 'c_LNA', 'c_Viruses', 'c_V1', 'c_V2', 'c_V3', key="count", value="value") %>%
    #gathering all the count information in two columns
    group_by(Location, Expt_No, Depth, Timepoint, Sample_Type, count ) %>% #grouping the dataframe in this order
    summarise(n =n(), mean=mean(value), sd=sd(value)) #calculating means and sd
  
  if ('VPC' %in% DF$Sample_Type){ #Here, if we have VPC samples, it'll create necessary columns for it
    colnames_mean<- c("VP", "VPC")
    colnames_sd<- c("VP", "VPC")
    
  } else {
    colnames_mean<- c("VP")
    colnames_sd<- c("VP")
    
  }
  
  DF_mean<- select(DF, -c('sd')) %>% #splitting the dataframe cause I haven't figure out how to spread the table without adding NAs
    spread('Sample_Type', 'mean')
  if (length(colnames_mean)==2){
    colnames(DF_mean)[7:8]<- colnames_mean
  } else if (length(colnames_mean)==1){
    colnames(DF_mean)[7]<- colnames_mean
  } 
  
  DF_sd<- select(DF, -c('mean')) %>%
    spread('Sample_Type', 'sd')
  if (length(colnames_sd)==2){
    colnames(DF_sd)[7:8]<- colnames_sd
  } else if (length(colnames_sd)==1){
    colnames(DF_sd)[7]<- colnames_sd
  } 
  
  if ('VPC' %in% DF$Sample_Type){ #If VPC exists, then we need to calculate Diff curve. Here we use simple subtarction between mean VPC and mean VP, and add their SDs
    DF_mean$Diff <- with(DF_mean, VPC-VP) #calculating Diff mean
    DF_mean<- pivot_longer(DF_mean, cols = c("VP", "VPC", "Diff"), names_to= 'Sample_Type', values_to='mean_value')
    DF_sd$Diff <- with(DF_sd, VPC+VP) #Calculating Diff sd, which is addition of the other sds
    DF_sd<- pivot_longer(DF_sd, cols = c("VP", "VPC", "Diff"), names_to='Sample_Type', values_to= 'sd_value')
  }
  
  #Merging both mean and sd dataframes
  DF<- merge(DF_mean, DF_sd, by= c('Location', 'Expt_No', 'Depth', 
                                   'Timepoint', 'count', 'n', 'Sample_Type'))
  
  DF<-  mutate(DF, Subgroup = case_when(DF$count == 'c_Bacteria' ~ "Total",
                                        DF$count == 'c_Viruses' ~ "Total",
                                        DF$count == 'c_HNA' ~ "Bacteria",
                                        DF$count == 'c_LNA' ~ "Bacteria",
                                        DF$count == 'c_V1' ~ "Viruses",
                                        DF$count == 'c_V2' ~ "Viruses",
                                        DF$count == 'c_V3' ~ "Viruses")) %>%
    arrange(Timepoint)
  
  rm('DF_mean', 'DF_sd', 'colnames_mean', 'colnames_sd')
  
  #For the time ranges
  TP<- unique(DF$Timepoint) #identifying the timepoints present
  colnames<- c()
  for (col in 2: length(TP)){
    a<- paste("T", TP[1], "_T", TP[col], sep = "")
    colnames[length(colnames)+1]<- a
  }
  colvalues<- c()
  for (col in 2: length(TP)){
    a<- paste("T", TP[1], ":T", TP[col], sep = "")
    colvalues[length(colvalues)+1]<- a
  }
  
  ncol<- ncol(DF)
  DF[colnames]<- NA
  
  #This works for 6 timepoints
  DF[,ncol+1]<- case_when(DF$Timepoint == TP[1] ~ colvalues[1],
                          DF$Timepoint == TP[2] ~ colvalues[1])
  
  
  DF[,ncol+2]<- case_when(DF$Timepoint == TP[1] ~ colvalues[2],
                          DF$Timepoint == TP[2] ~ colvalues[2],
                          DF$Timepoint == TP[3] ~ colvalues[2])
  
  DF[,ncol+3]<- case_when(DF$Timepoint == TP[1] ~ colvalues[3],
                          DF$Timepoint == TP[2] ~ colvalues[3],
                          DF$Timepoint == TP[3] ~ colvalues[3],
                          DF$Timepoint == TP[4] ~ colvalues[3])
  
  DF[,ncol+4]<- case_when(DF$Timepoint == TP[1] ~ colvalues[4],
                          DF$Timepoint == TP[2] ~ colvalues[4],
                          DF$Timepoint == TP[3] ~ colvalues[4],
                          DF$Timepoint == TP[4] ~ colvalues[4],
                          DF$Timepoint == TP[5] ~ colvalues[4])
  
  DF[,ncol+5]<- case_when(DF$Timepoint == TP[1] ~ colvalues[5],
                          DF$Timepoint == TP[2] ~ colvalues[5],
                          DF$Timepoint == TP[3] ~ colvalues[5],
                          DF$Timepoint == TP[4] ~ colvalues[5],
                          DF$Timepoint == TP[5] ~ colvalues[5],
                          DF$Timepoint == TP[6] ~ colvalues[5])
  
  
  DF<- DF %>%
    pivot_longer(cols = colnames, names_to = "Time_Range", values_to = "Time_Time")%>%
    drop_na()
  
  rm('colnames', 'colvalues', 'TP', 'a', 'ncol')
  return(DF)
}

#To create an overview plot with time ranges
overview_plots_tp_avg<- function(df, ...){
  n<- ggplot(df, aes(x= Timepoint, y= mean_value, color= count , shape=count))+
    geom_point(size= 1.5)+
    geom_smooth(size= 1.0, method = 'lm', se = F)+
    geom_line()+
    geom_hline(yintercept = 0, color= '#636363', size= 0.3, linetype = "dashed")+
    facet_grid(factor(Subgroup, levels = c("Total", "Bacteria", "Viruses")) +
                 factor(Sample_Type, levels = c("VP", "VPC", "Diff"))~
                 factor(Time_Time, levels = unique(df$Time_Time))  
    )+
    scale_color_manual(name= 'Populations',
                       labels=c("Total Bacteria", "HNA Bacteria", "LNA Bacteria",
                                "Total Viruses", "V1 Viruses", "V2 Viruses", "V3 Viruses"),
                       values= c(c_Bacteria = "#AD002A99", 
                                 c_HNA = "#00468B99",
                                 c_LNA = "#0099B499",
                                 c_Viruses = "#ED000099",
                                 c_V1 = "#1B191999",
                                 c_V2 = "#7E6148B2",
                                 c_V3 = "#9C9C9C"))+
    scale_shape_manual(name = 'Populations', 
                       values = c(16,16,16,17,17,17,17),
                       labels=c("Total Bacteria", "HNA Bacteria", "LNA Bacteria",
                                "Total Viruses", "V1 Viruses", "V2 Viruses", "V3 Viruses"))+
    scale_y_continuous(labels = unit_format(unit = NULL, scale = 1e-6),
                       breaks= seq(-4e+6,13e+6, 4e+6),
                       limits = c(-4e+6, 14e+6))+
    theme_bw()+
    geom_errorbar(aes(ymin=mean_value + sd_value, ymax= mean_value - sd_value), width = 0.5, size = 0.5)+
    scale_x_continuous(breaks = unique(NJ1$Timepoint))+
    labs(title = 'NJ1 - Viral Production - Timepoint Sloughing', subtitle = 'Overview - Bacterial and Viral counts for Lytic and Lysogenic inductions',
         x= 'Sampling Timepoints\n (in hours)', y='FCM Counts\n (in millions)')+
    theme(strip.text = element_text(face = "bold"),
          strip.background = element_rect( color = 'black', fill = 'white'),
          axis.title = element_text(face = 'bold'),
          title = element_text(face = 'bold'),
          legend.position = "bottom")+
    guides(color = guide_legend(nrow = 2, byrow = TRUE))
  
  return(n)
  
}


#### Bacterial Generation Time####
Ba_gt<- function(x){
  GT<- (log10(2)*(Ba$Timepoint[x]-Ba$Timepoint[1]))/(log10(Ba$VP_mean[x])-log10(Ba$VP_mean[1]))
  print(GT) 
}


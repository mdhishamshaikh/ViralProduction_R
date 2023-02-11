
####Dataframe functions####


#Creating four dataframe functions.
#1. separate replicates
#2. averaged replicate
#3. SR with time ranges
#4. AVG with time ranges

#5. Time ranges
tp<- function(DF){
  TP<- unique(DF$Timepoint)
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

#separate replicate dataframe
df_sr<- function(df, keep_0.22 = F){
  DF<- df%>%
    select(c('Location', 'Expt_No', 'Depth', 'Sample_Type', 'Timepoint', 'Replicate',
             'c_Bacteria', 'c_HNA', 'c_LNA', 'c_Viruses', 'c_V1', 'c_V2', 'c_V3'))%>%
    gather('c_Bacteria', 'c_HNA', 'c_LNA', 'c_Viruses', 'c_V1', 'c_V2', 'c_V3', key="Population", value="count")%>%
    mutate(Microbe = if_else(Population == 'c_Bacteria' | Population == 'c_HNA' | Population == 'c_LNA', "Bacteria", "Viruses"))%>%
    arrange('Location', 'Expt_No', 'Depth', 'Sample_Type', 'Timepoint', 'Replicate','Population')
  if (keep_0.22 == F){
    DF<- DF[DF$Sample_Type != '0.22',]
  }
  return(DF)
}


df_avg<- function(df, keep_0.22 = F) {
  #only works if we remove 0.22
  DF<- df[df$Sample_Type != '0.22',]
  
  DF<- DF %>%
    gather('c_Bacteria', 'c_HNA', 'c_LNA', 'c_Viruses', 'c_V1', 'c_V2', 'c_V3', key="Population", value="count") %>%
    group_by(Location, Expt_No, Depth, Sample_Type, Timepoint, Population ) %>%
    summarise(n =n(), mean=mean(count), sd=sd(count)) #calculating means and sd
  
  if ('VPC' %in% DF$Sample_Type){
    colnames_mean<- c("VP", "VPC")
    colnames_sd<- c("VP", "VPC")
    
  } else {
    colnames_mean<- c("VP")
    colnames_sd<- c("VP")
    
  }
  
  DF_mean<- select(DF, -c('sd')) %>% #splitting the dataframe cause I haven't figure out how to spread teh table without adding NAs
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
  
  if ('VPC' %in% DF$Sample_Type){
    DF_mean$Diff <- with(DF_mean, VPC-VP) #calculating Diff mean
    DF_mean<- pivot_longer(DF_mean, cols = c("VP", "VPC", "Diff"), names_to= 'Sample_Type', values_to='mean')
    DF_sd$Diff <- with(DF_sd, VPC+VP) #Calculating Diff sd, which is addition of the other sds
    DF_sd<- pivot_longer(DF_sd, cols = c("VP", "VPC", "Diff"), names_to='Sample_Type', values_to= 'sd')
  }
  
  DF<- merge(DF_mean, DF_sd, by= c('Location', 'Expt_No', 'Depth',
                                   'Timepoint', 'Population', 'n', 'Sample_Type')) %>%
    mutate(Microbe = if_else(Population == 'c_Bacteria' | Population == 'c_HNA' | Population == 'c_LNA', "Bacteria", "Viruses"))%>%
    mutate(Subgroup = if_else(Population == 'c_Bacteria' | Population == 'c_Viruses', "Parent", "Subgroup"))
  rm('DF_mean', 'DF_sd', 'colnames_mean', 'colnames_sd')
  
  DF<- DF%>%
    arrange('Location',
            'Expt_No',
            'Depth',
            'Sample_Type',
            'Population',
            Timepoint)
  
  return(DF)
}


#df average replicates timepoints

df_avg_tp<- function(df, keep_0.22 = F){
  DF<- df_avg(df)
  
  TP<- unique(DF$Timepoint)
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


df_sr_tp<- function(df, keep_0.22 = F){
  DF<- df_sr(df)
  
  TP<- unique(DF$Timepoint)
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
  
  rm('colnames', 'colvalues', 'TP',  'a', 'ncol')
  return(DF)
}



####Peaks and Valleys####

#findpeaks and findvalleys are derived from findPeaks and findValleys, respectively,
#from quantmod R package

#without standard deviation
# findpeaks<- function(x){
#   which(diff(sign(diff(x, na.pad = FALSE)), na.pad = FALSE) < 
#           0) + 1
# }
# 
# findvalleys<-function(x){
#   which(diff(sign(diff(x, na.pad = FALSE)), na.pad = FALSE) > 
#           0) + 1
# }
{
  library(tidyverse)
  library(ggsci)
  library(scales)
}

#with standard deviation

peaks<- function(values){
  list<- c()
  
  for (len in 1:(length(values)-1)){
    d<- sign(values[len+1] - values[len])
    #print(d)
    list[[length(list)+1]]<- d
  }
  which(diff(as.numeric(list))<0) +2 -1
  
}


valleys<- function(values){
  list<- c()
  for (len in 1:(length(values)-1)){
    d<- sign(values[len+1] - values[len])
    #print(d)
    list[[length(list)+1]]<- d
  }
  which(diff(as.numeric(list))>0) +2 -1
}


peaks_sd<- function(values, sd){
  list<- c()
  
  for (len in 1:(length(values)-1)){
    d<- sign((values[len+1] - sd[len+1]) - (values[len] + sd[len]))
    print(d)
    list[[length(list)+1]]<- d
    #print(list)
  }
  which(diff(as.numeric(list))<0) +2 -1
  
}


valleys_sd<- function(values, sd){
  list<- c()
  
  for (len in 1:(length(values)-1)){
    d<- sign((values[len+1] - sd[len+1]) - (values[len] + sd[len]))
    print(d)
    list[[length(list)+1]]<- d
  }
  which(diff(as.numeric(list))>0) +2 -1
  
}

plots_vipcal<- function(df){

ggplot(df, aes(x= Timepoint, y= mean_value, color= count , shape=count))+
  geom_point(size= 1.5)+
  geom_line(size= 0.6)+
  geom_hline(yintercept = 0, color= '#636363', size= 0.3, linetype = "dashed")+
  facet_grid(factor(Subgroup, 
                    levels = c("Total", "Bacteria", "Viruses")) ~ 
               factor(Sample_Type, levels = c("VP", "VPC", "Diff")))+
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
                     breaks= seq(-4e+6,13e+6, 2e+6),
                     limits = c(-4e+6, 14e+6))+
  theme_bw()+
  geom_errorbar(aes(ymin=mean_value + sd_value, ymax= mean_value - sd_value), width = 0.5, size = 0.5)+
  scale_x_continuous(breaks = unique(NJ1$Timepoint))+
  labs(title = 'NJ1 - Viral Production - VIPCAL', subtitle = 'Overview - Bacterial and Viral counts for Lytic and Lysogenic inductions',
       x= 'Sampling Timepoints\n (in hours)', y='FCM Counts\n (in millions)')+
  theme(strip.text = element_text(face = "bold"),
        strip.background = element_rect( color = 'black', fill = 'white'),
        axis.title = element_text(face = 'bold'),
        title = element_text(face = 'bold'))
}


plots_lm<- function(df){
  
  ggplot(df, aes(x= Timepoint, y= mean_value, color= count , shape=count))+
    geom_point(size= 1.5)+
    geom_smooth(size= 1.0, method = 'lm', se = T)+
    geom_hline(yintercept = 0, color= '#636363', size= 0.3, linetype = "dashed")+
    facet_grid(factor(Subgroup, 
                      levels = c("Total", "Bacteria", "Viruses")) ~ 
                 factor(Sample_Type, levels = c("VP", "VPC", "Diff")))+
    scale_color_manual(name= 'Populations',
                       labels=c("Total Bacteria", "HNA Bacteria", "LNA Bacteria",
                                "Total Viruses", "V1 Viruses", "V2 Viruses", "V3 Viruses"),
                       breaks=c("Total Bacteria", "HNA Bacteria", "LNA Bacteria",
                                "Total Viruses", "V1 Viruses", "V2 Viruses", "V3 Viruses"),
                       values= c(c_Bacteria = "#AD002A99", 
                                 c_HNA = "#00468B99",
                                 c_LNA = "#0099B499",
                                 c_Viruses = "#ff8a8a",
                                 c_V1 = "#1B191999",
                                 c_V2 = "#7E6148B2",
                                 c_V3 = "#9C9C9C"))+
    scale_shape_manual(name = 'Populations', 
                       values = c(16,16,16,17,17,17,17),
                       labels=c("Total Bacteria", "HNA Bacteria", "LNA Bacteria",
                                "Total Viruses", "V1 Viruses", "V2 Viruses", "V3 Viruses"))+
    scale_y_continuous(labels = unit_format(unit = NULL, scale = 1e-6),
                       breaks= seq(-4e+6,13e+6, 2e+6),
                       limits = c(-4e+6, 14e+6))+
    theme_bw()+
    geom_errorbar(aes(ymin=mean_value + sd_value, ymax= mean_value - sd_value), width = 0.5, size = 0.5)+
    scale_x_continuous(breaks = unique(NJ1$Timepoint))+
    labs(title = 'NJ1 - Viral Production - LM', subtitle = 'Overview - Bacterial and Viral counts for Lytic and Lysogenic inductions',
         x= 'Sampling Timepoints\n (in hours)', y='FCM Counts\n (in millions)')+
    theme(strip.text = element_text(face = "bold"),
          strip.background = element_rect( color = 'black', fill = 'white'),
          axis.title = element_text(face = 'bold'),
          title = element_text(face = 'bold'))
}



####Overview Plots####

df_lm_tp<- function(df, ...){
  df<- df[df$Sample_Type != '0.22',] #Lose 0.22 values
  
  df<- gather(df, 'c_Bacteria', 'c_HNA', 'c_LNA', 'c_Viruses', 'c_V1', 'c_V2', 'c_V3', key="count", value="value") %>%
    group_by(Location, Expt_No, Depth, Sample_Type, Timepoint, count ) %>%
    summarise(n =n(), mean=mean(value), sd=sd(value)) #calculating means and sd
  
  df_mean<- df[,1:8] %>% #splitting the dataframe cause I haven't figure out how to spread the table without adding NAs
    spread('Sample_Type', 'mean')
  colnames(df_mean)[7:8]<- c("VP_mean", "VPC_mean")
  df_mean$Diff_mean <- with(df_mean, VPC_mean-VP_mean) #calcualting Diff mean
  
  df_sd<- df[,c(1:7,9)] %>%
    spread('Sample_Type', 'sd')
  colnames(df_sd)[7:8]<- c("VP_sd", "VPC_sd")
  df_sd$Diff_sd <- with(df_sd, VPC_sd+VP_sd) #Calculating Diff sd, which is addition of the other sds
  
  
  df<- merge(df_mean, df_sd, by= c('Location', 'Expt_No', 'Depth',
                                      'Timepoint', 'count', 'n')) %>%
    mutate(Microbe = if_else(count == 'c_Bacteria' | count == 'c_HNA' | count == 'c_LNA', "Bacteria", "Viruses"))%>%
    mutate(Subgroup = if_else(count == 'c_Bacteria' | count == 'c_Viruses', "Parent", "Subgroup"))
  
  rm(df_mean)
  rm(df_sd)
  TP<- unique(NJ1$Timepoint)
  colnames<- c()
  for (col in 2: length(TP)){
    a<- paste("T", TP[1], "_T", TP[col], sep = "")
    colnames[length(colnames)+1]<- a
  }
  df<- df%>%
    mutate("T0_T3" = case_when(Timepoint == '0' ~ "T0:T3",
                               Timepoint == '3' ~ "T0:T3"))%>%
    
    mutate("T0_T6" = case_when(Timepoint == '0' ~ "T0:T6",
                               Timepoint == '3' ~ "T0:T6", 
                               Timepoint == '6' ~ "T0:T6"))%>%
    
    mutate("T0_T17" = case_when(Timepoint == '0' ~ "T0:T17",
                                Timepoint == '3' ~ "T0:T17", 
                                Timepoint == '6' ~ "T0:T17",
                                Timepoint == '17' ~ "T0:T17"))%>%
    
    mutate("T0_T20" = case_when(Timepoint == '0' ~ "T0:T20",
                                Timepoint == '3' ~ "T0:T20",
                                Timepoint == '6' ~ "T0:T20",
                                Timepoint == '17' ~ "T0:T20",
                                Timepoint == '20' ~ "T0:T20"))%>%
    
    mutate("T0_T24" = case_when(Timepoint == '0' ~ "T0:T24",
                                Timepoint == '3' ~ "T0:T24",
                                Timepoint == '6' ~ "T0:T24",
                                Timepoint == '17' ~ "T0:T24",
                                Timepoint == '20' ~ "T0:T24",
                                Timepoint == '24' ~ "T0:T24")) %>%
    pivot_longer(cols = colnames, names_to = "Time_Range", values_to = "Time_Time")%>%
    drop_na()
  
  
  return(df)
}

plots_lm_tp<- function(df, ...){
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
  
  
  o<- ggplot_gtable(ggplot_build(n))
  strip_both<- which(grepl('strip-', o$layout$name))
  fills<- c(rep("white", 14))
  fills[bp_endpoint]<- "#E489A6" #Change the index to that of the time range after bacterial generation time
  k <- 1
  
  for (i in strip_both) {
    j<- which(grepl('rect', o$grobs[[i]]$grobs[[1]]$childrenOrder))
    o$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k<- k+1
  }
  #https://ojkftyk.blogspot.com/2019/01/r-ggplot2-change-colour-of-font-and.html
  grid::grid.draw(o)
}




#####Slope Functions####
slope_lm_sr<- function(df_sr){ #takes SR dataframe as an input
  
  lm_vp<- list()
  slope_lm_sr_df<- data.frame()
  for (location in unique(df_sr$Location)){
    for (expt_no in unique(df_sr$Expt_No)){
      for (depth in unique(df_sr$Depth)){
        for (time in unique(df_sr$Time_Range)){
          for (viruses in unique(df_sr[df_sr$Microbe == "Viruses",]$Population)){
            for (prod in unique(df_sr$Sample_Type)){
              for (rep in unique(df_sr$Replicate)){
                
                df2<- df_sr[df_sr$Location == location,]
                df2<- df2[df2$Expt_No == expt_no,]
                df2<- df2[df2$Depth == depth,]
                df2<- df2[df2$Time_Range == time,]
                df2<- df2[df2$Population == viruses,]
                df2<- df2[df2$Sample_Type == prod,]
                df2<- df2[df2$Replicate == rep,]
                
                lm<- summary(lm(data = df2, count ~ as.numeric(Timepoint)))
                print(lm)
                slope<- c(location, expt_no, depth, time, viruses, prod, rep, lm$coefficients[c(2,4)], lm$r.squared)
                print(lm$coefficients[,3])
                print(lm$coefficients[[2]])
                lm_vp[[length(lm_vp)+1]] <- slope  
              }
            }
          }
        }
      }
    }
  }
  
  slope_lm_sr_df<- data.frame(t(sapply(lm_vp, c)))
  colnames(slope_lm_sr_df)<- c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'Replicate', 'LM_SR_Slope', 'LM_SR_SE', 'LM_SR_R_squared')
  slope_lm_sr_df[, c('LM_SR_Slope', 'LM_SR_SE', 'LM_SR_R_squared')]<- lapply(slope_lm_sr_df[, c('LM_SR_Slope', 'LM_SR_SE', 'LM_SR_R_squared')], as.numeric)
  return(slope_lm_sr_df)
}

slope_lm_avg<- function(df_avg){ #takes AR dataframe as an input
  
  lm_vp<- list()
  slope_lm_avg_df<- data.frame()
  for (location in unique(df_avg$Location)){
    for (expt_no in unique(df_avg$Expt_No)){
      for (depth in unique(df_avg$Depth)){
        for (time in unique(df_avg$Time_Range)){
          for (viruses in unique(df_avg[df_avg$Microbe == "Viruses",]$Population)){
            for (prod in unique(df_avg$Sample_Type)){
            
                df2<- df_avg[df_avg$Location == location,]
                df2<- df2[df2$Expt_No == expt_no,]
                df2<- df2[df2$Depth == depth,]
                df2<- df2[df2$Time_Range == time,]
                df2<- df2[df2$Population == viruses,]
                df2<- df2[df2$Sample_Type == prod,]
                
                lm<- summary(lm(data = df2, mean ~ Timepoint))
                print(lm)
                slope<- c(location, expt_no, depth, time, viruses, prod, lm$coefficients[c(2,4)], lm$r.squared)
                lm_vp[[length(lm_vp)+1]] <- slope  
              }
            }
          }
        }
      }
    }

  
  slope_lm_avg_df<- data.frame(t(sapply(lm_vp, c)))
  colnames(slope_lm_avg_df)<- c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'LM_AVG_Slope', 'LM_AVG_SE', 'LM_AVG_R_squared')
  slope_lm_avg_df[, c('LM_AVG_Slope', 'LM_AVG_SE', 'LM_AVG_R_squared')]<- lapply(slope_lm_avg_df[, c('LM_AVG_Slope', 'LM_AVG_SE', 'LM_AVG_R_squared')], as.numeric)
  return(slope_lm_avg_df)
}

slope_all_points<- function(df_sr){ #takes df_sr as the input
  lm_vp<- list()
  slope_lm_all_points_df<- data.frame()
  for (location in unique(df_sr$Location)){
    for (expt_no in unique(df_sr$Expt_No)){
      for (depth in unique(df_sr$Depth)){
        for (time in unique(df_sr$Time_Range)){
          for (viruses in unique(df_sr[df_sr$Microbe == "Viruses",]$Population)){
            for (prod in unique(df_sr$Sample_Type)){
              
                df2<- df_sr[df_sr$Location == location,]
                df2<- df2[df2$Expt_No == expt_no,]
                df2<- df2[df2$Depth == depth,]
                df2<- df2[df2$Time_Range == time,]
                df2<- df2[df2$Population == viruses,]
                df2<- df2[df2$Sample_Type == prod,]
                
                lm<- summary(lm(data = df2, count ~ as.numeric(Timepoint)))
                print(lm)
                slope<- c(location, expt_no, depth, time, viruses, prod, lm$coefficients[c(2,4)], lm$r.squared)
                lm_vp[[length(lm_vp)+1]] <- slope  
            }
          }
        }
      }
    }
  }
  
  slope_lm_all_points_df<- data.frame(t(sapply(lm_vp, c)))
  colnames(slope_lm_all_points_df)<- c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'LM_AP_Slope', 'LM_AP_SE', 'LM_AP_R_squared')
  slope_lm_all_points_df[, c('LM_AP_Slope', 'LM_AP_SE', 'LM_AP_R_squared')]<- lapply(slope_lm_all_points_df[, c('LM_AP_Slope', 'LM_AP_SE', 'LM_AP_R_squared')], as.numeric)
  return(slope_lm_all_points_df)
}


####VIPCAL Functions####
vipcal_sr<- function(df_sr){ #takes SR dataframe as an input
  
  vipcal_vp<- list()
  vpcl_vp_sr_df<- data.frame()
  for (location in unique(df_sr$Location)){
    for (expt_no in unique(df_sr$Expt_No)){
      for (depth in unique(df_sr$Depth)){
        for (time in unique(df_sr$Time_Range)){
          for (viruses in unique(df_sr[df_sr$Microbe == "Viruses",]$Population)){
            for (prod in unique(df_sr$Sample_Type)){
              for (rep in unique(df_sr$Replicate)){
                
                df2<- df_sr[df_sr$Location == location,]
                df2<- df2[df2$Expt_No == expt_no,]
                df2<- df2[df2$Depth == depth,]
                df2<- df2[df2$Time_Range == time,]
                df2<- df2[df2$Population == viruses,]
                df2<- df2[df2$Sample_Type == prod,]
                df2<- df2[df2$Replicate == rep,]
                
                
                p<- peaks(c(+10e+10, df2$count, -10e+10))-1
                v<- valleys(c(+10e+10, df2$count, -10e+10))-1
                
                if(identical(length(p), length(v))){
                  print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
                }else{
                  print("Number peaks and valleys are not the same. This will lead to erroneous viral production calculations")
                }
                
                if(length(p)==0){
                  vp<- NA
                  
                } else if (length(p)==1) {
                  vp<- (df2$count[p[1]] - df2$count[v[1]])/(df2$Timepoint[p[1]] - df2$Timepoint[v[1]])
                  
                } else if (length(p)==2) {
                  vp<- ((df2$count[p[1]] - df2$count[v[1]])/(df2$Timepoint[p[1]] - df2$Timepoint[v[1]]) + 
                          (df2$count[p[2]] - df2$count[v[2]])/(df2$Timepoint[p[2]] - df2$Timepoint[v[2]]))/2
                  
                } else if (length(p)==3) {
                  vp<-  ((df2$count[p[1]] - df2$count[v[1]])/(df2$Timepoint[p[1]] - df2$Timepoint[v[1]]) + 
                           (df2$count[p[2]] - df2$count[v[2]])/(df2$Timepoint[p[2]] - df2$Timepoint[v[2]]) +
                           (df2$count[p[3]] - df2$count[v[3]])/(df2$Timepoint[p[3]] - df2$Timepoint[v[3]]))/3
                }
                
                
                vipcal<- c(location, expt_no, depth, time, viruses, prod, rep, vp)
                
                vipcal_vp[[length(vipcal_vp)+1]] <- vipcal  
              }
            }
          }
        }
      }
    }
  }
  
  vpcl_vp_sr_df<- data.frame(t(sapply(vipcal_vp, c)))
  colnames(vpcl_vp_sr_df)<- c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'Replicate', 'VPCL_SR')
  vpcl_vp_sr_df$VPCL_SR<- as.numeric(vpcl_vp_sr_df$VPCL_SR)
  return(vpcl_vp_sr_df)
}

vipcal_avg_sd<- function(df_avg){ #takes avg_tp dataframe as an input
  
  vipcal_vp<- list()
  vpcl_vp_avg_df<- data.frame()
  for (location in unique(df_avg$Location)){
    for (expt_no in unique(df_avg$Expt_No)){
      for (depth in unique(df_avg$Depth)){
        for (time in unique(df_avg$Time_Range)){
          for (viruses in unique(df_avg[df_avg$Microbe == "Viruses",]$Population)){
            for (prod in unique(df_avg$Sample_Type)){
              
              
              df2<- df_avg[df_avg$Location == location,]
              df2<- df2[df2$Expt_No == expt_no,]
              df2<- df2[df2$Depth == depth,]
              df2<- df2[df2$Time_Range == time,]
              df2<- df2[df2$Population == viruses,]
              df2<- df2[df2$Sample_Type == prod,]
              
              
              p<- peaks_sd(c(+10e+10, df2$mean, -10e+10),
                           c(+10e+10, df2$sd, -10e+10))-1
              v<- valleys_sd(c(+10e+10, df2$mean, -10e+10),
                             c(+10e+10, df2$sd, -10e+10))-1
              
              if(identical(length(p), length(v))){
                print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
              }else{
                print("Number peaks and valleys are not the same. This will lead to erroneous viral production calculations")
              }
              
              if(length(p)==0){
                vp<- NA
                
              } else if (length(p)==1) {
                vp<- (df2$mean[p[1]] - df2$mean[v[1]])/(df2$Timepoint[p[1]] - df2$Timepoint[v[1]])
                
              } else if (length(p)==2) {
                vp<- ((df2$mean[p[1]] - df2$mean[v[1]])/(df2$Timepoint[p[1]] - df2$Timepoint[v[1]]) + 
                        (df2$mean[p[2]] - df2$mean[v[2]])/(df2$Timepoint[p[2]] - df2$Timepoint[v[2]]))/2
                
              } else if (length(p)==3) {
                vp<-  ((df2$mean[p[1]] - df2$mean[v[1]])/(df2$Timepoint[p[1]] - df2$Timepoint[v[1]]) + 
                         (df2$mean[p[2]] - df2$mean[v[2]])/(df2$Timepoint[p[2]] - df2$Timepoint[v[2]]) +
                         (df2$mean[p[3]] - df2$mean[v[3]])/(df2$Timepoint[p[3]] - df2$Timepoint[v[3]]))/3
              }
              
              
              vipcal<- c(location, expt_no, depth, time, viruses, prod, vp)
              
              vipcal_vp[[length(vipcal_vp)+1]] <- vipcal  
              
            }
          }
        }
      }
    }
  }
  
  vpcl_vp_avg_df<- data.frame(t(sapply(vipcal_vp, c)))
  colnames(vpcl_vp_avg_df)<- c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population', 'Sample_Type',  'VPCL_AVG_SE')
  vpcl_vp_avg_df$VPCL_AVG_SE<- as.numeric(vpcl_vp_avg_df$VPCL_AVG_SE)
  return(vpcl_vp_avg_df)
}

vipcal_avg<- function(df_avg){ #takes avg_tp dataframe as an input
  
  vipcal_vp<- list()
  vpcl_vp_avg_df<- data.frame()
  for (location in unique(df_avg$Location)){
    for (expt_no in unique(df_avg$Expt_No)){
      for (depth in unique(df_avg$Depth)){
        for (time in unique(df_avg$Time_Range)){
          for (viruses in unique(df_avg[df_avg$Microbe == "Viruses",]$Population)){
            for (prod in unique(df_avg$Sample_Type)){
              
              
              df2<- df_avg[df_avg$Location == location,]
              df2<- df2[df2$Expt_No == expt_no,]
              df2<- df2[df2$Depth == depth,]
              df2<- df2[df2$Time_Range == time,]
              df2<- df2[df2$Population == viruses,]
              df2<- df2[df2$Sample_Type == prod,]
              
              
              p<- peaks(c(+10e+10, df2$mean, -10e+10))-1
              v<- valleys(c(+10e+10, df2$mean, -10e+10))-1
              
              if(identical(length(p), length(v))){
                print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
              }else{
                print("Number peaks and valleys are not the same. This will lead to erroneous viral production calculations")
              }
              
              if(length(p)==0){
                vp<- NA
                
              } else if (length(p)==1) {
                vp<- (df2$mean[p[1]] - df2$mean[v[1]])/(df2$Timepoint[p[1]] - df2$Timepoint[v[1]])
                
              } else if (length(p)==2) {
                vp<- ((df2$mean[p[1]] - df2$mean[v[1]])/(df2$Timepoint[p[1]] - df2$Timepoint[v[1]]) + 
                        (df2$mean[p[2]] - df2$mean[v[2]])/(df2$Timepoint[p[2]] - df2$Timepoint[v[2]]))/2
                
              } else if (length(p)==3) {
                vp<-  ((df2$mean[p[1]] - df2$mean[v[1]])/(df2$Timepoint[p[1]] - df2$Timepoint[v[1]]) + 
                         (df2$mean[p[2]] - df2$mean[v[2]])/(df2$Timepoint[p[2]] - df2$Timepoint[v[2]]) +
                         (df2$mean[p[3]] - df2$mean[v[3]])/(df2$Timepoint[p[3]] - df2$Timepoint[v[3]]))/3
              }
              
              
              vipcal<- c(location, expt_no, depth, time, viruses, prod, vp)
              
              vipcal_vp[[length(vipcal_vp)+1]] <- vipcal  
              
            }
          }
        }
      }
    }
  }
  
  vpcl_vp_avg_df<- data.frame(t(sapply(vipcal_vp, c)))
  colnames(vpcl_vp_avg_df)<- c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population', 'Sample_Type',  'VPCL_AVG')
  vpcl_vp_avg_df$VPCL_AVG<- as.numeric(vpcl_vp_avg_df$VPCL_AVG)
  return(vpcl_vp_avg_df)
}



####To calculate lysogeny from all points data####
calc_diff_lm_AP<- function(df){
  
  VP<- df[df$Sample_Type == "VP",]
  VPC<- df[df$Sample_Type == "VPC",]
  Diff<- as.data.frame(cbind(VP[,1:6], VPC[,7] -VP[,7], VPC[,8] + VPC[,8]))
  Diff[,6]<- "Diff"
  colnames(Diff)[7]<- 'LM_AP_Slope'
  colnames(Diff)[8]<- 'LM_AP_SE'
  Diff$LM_AP_R_squared <- NA
  
  df<- full_join(df, Diff, by = NULL)
  return(df)
}

calc_diff_lm_SR<- function(df){
  
  VP<- df[df$Sample_Type == "VP",]
  VPC<- df[df$Sample_Type == "VPC",]
  Diff<- as.data.frame(cbind(VP[,1:6], VPC[,7] -VP[,7], VPC[,8] + VPC[,8]))
  Diff[,4]<- "Diff"
  colnames(Diff)[7]<- 'LM_SR_Slope_Mean'
  colnames(Diff)[8]<- 'LM_SR_Slope_SE'
  
  df<- full_join(df, Diff, by = NULL)
  return(df)
}

calc_diff_lm_AR<- function(df){
  
  VP<- df[df$Sample_Type == "VP",]
  VPC<- df[df$Sample_Type == "VPC",]
  Diff<- as.data.frame(cbind(VP[,1:6], VPC[,7] -VP[,7], VPC[,8] + VPC[,8]))
  Diff[,6]<- "Diff"
  colnames(Diff)[7]<- 'LM_AVG_Slope'
  colnames(Diff)[8]<- 'LM_AVG_SE'
  
  df<- full_join(df, Diff, by = NULL)
  return(df)
}


####LMER Diff Curve Extraction####
lmer_model<- function(df, value = count){
  
  lmer_data<- data.frame()
  model_plots<- list()
  n<-0
  for (rep in c(1,2,3)){
    df$Replicate[df$Replicate == rep & df$Sample_Type == 'VPC'] <- rep+3
  }
  
  for (virus in unique(df$Population)){
    
    model<- lme4::lmer(data = df, value ~ Sample_Type*as.factor(Timepoint) + (1+ Sample_Type | Replicate))
    plot<- model_plot(model, df = df)
    emmeans<- emmeans::emmeans(model, ~ Sample_Type|as.factor(Timepoint))
    contrast<- pairs(emmeans) 
    df1<- data.frame(rep(virus, 6),
                     rep("Diff", 6),
                     summary(contrast)$Timepoint, 
                     -(summary(contrast)$estimate),
                     summary(contrast)$SE
    )
    colnames(df1)<- c("Population", "Sample_Type", "Timepoint", "Mean", "SE")
    st<- summary(emmeans)[1]
    tp<- summary(emmeans)[2]
    mean<- summary(emmeans)[3]
    se<- summary(emmeans)[4]
    pop<- rep(virus, 6)
    df2<- data.frame(pop,st,tp,mean,se)
    colnames(df2)<- c("Population", "Sample_Type", "Timepoint", "Mean", "SE")
    
    
    lmer_data<- rbind(lmer_data, df2, df1)
    # DF<- rbind(df2, df1) %>%
    #  arrange(Timepoint)
    #n<- n +1
    #lmer_data[[n]]<- DF
    #model_plots[[n]]<- plot
    
  }
  return(lmer_data)
  return(model_plots)
  
  
}

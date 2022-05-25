
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
peaks<- function(x,y){
  list<- c()
for (i in 1:7){
  g<- (sign((x[i+1]-y[i+1]) -  (x[i]+y[i])))
  #print(g)
  list[[length(list)+1]]<- g
}

which(diff(as.numeric(list))<0) +2 -1 #peaks
}

valleys<- function(x,y){
  list<- c()
  for (i in 1:7){
    g<- (sign((x[i+1]-y[i+1]) -  (x[i]+y[i])))
    #print(g)
    list[[length(list)+1]]<- g
  }
 
  which(diff(as.numeric(list))>0) +2 -1 #peaks
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
    labs(title = 'NJ1 - Viral Production - LM', subtitle = 'Overview - Bacterial and Viral counts for Lytic and Lysogenic inductions',
         x= 'Sampling Timepoints\n (in hours)', y='FCM Counts\n (in millions)')+
    theme(strip.text = element_text(face = "bold"),
          strip.background = element_rect( color = 'black', fill = 'white'),
          axis.title = element_text(face = 'bold'),
          title = element_text(face = 'bold'))
}

df_lm_tp<- function(df){
  df<- df[df$Sample_Type != '0.22',] #Lose 0.22 values
  
  df<- gather(df, 'c_Bacteria', 'c_HNA', 'c_LNA', 'c_Viruses', 'c_V1', 'c_V2', 'c_V3', key="count", value="value") %>%
    group_by(Location, Expt_No, Depth, Sample_Type, Timepoint, count ) %>%
    summarise(n =n(), mean=mean(value), sd=sd(value)) #calculating means and sd
  
  df_mean<- df[,1:8] %>% #splitting the datfram cause I haven't figure out how to spread teh table without adding NAs
    spread('Sample_Type', 'mean')
  colnames(df_mean)[7:8]<- c("VP_mean", "VPC_mean")
  df_mean$Diff_mean <- with(df_mean, VPC_mean-VP_mean) #calcualting Diff mean
  
  df_sd<- df[,c(1:7,9)] %>%
    spread('Sample_Type', 'sd')
  colnames(df_sd)[7:8]<- c("VP_sd", "VPC_sd")
  df_sd$Diff_sd <- with(df_sd, VPC_sd+VP_sd) #Calculating Diff sd, whcih si addition of the other sds
  
  
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

plots_lm_tp<- function(df){
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

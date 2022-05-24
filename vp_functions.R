
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


plots_lm_tp<- function(df){
ggplot(df, aes(x= Timepoint, y= mean_value, color= count , shape=count))+
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
                     breaks= seq(-4e+6,13e+6, 2e+6),
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
}


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

plots_vp<- function(df){

ggplot(df, aes(x= Timepoint, y= mean_value, color= count , shape=count))+
  geom_point(size= 1.5)+
  geom_line(size= 0.8)+
  facet_grid(factor(Subgroup, levels = c("Total", "Bacteria", "Viruses"))~factor(Sample_Type, levels = c("VP", "VPC", "Diff")))+
  scale_color_manual(name= 'Populations',
                     labels=c("Total Bacteria", "HNA Bacteria", "LNA Bacteria",
                              "Total Viruses", "V1 Viruses", "V2 Viruses", "V3 Viruses"),
                     values= c(c_Bacteria = "#00468BFF", 
                               c_HNA = "#ED0000FF",
                               c_LNA = "#42B540FF",
                               c_Viruses = "#0099B4FF",
                               c_V1 = "#925E9FFF",
                               c_V2 = "#FDAF91FF",
                               c_V3 = "#AD002AFF"))+
  scale_shape_manual(name = 'Populations', 
                     values = c(16,16,16,17,17,17,17),
                     labels=c("Total Bacteria", "HNA Bacteria", "LNA Bacteria",
                              "Total Viruses", "V1 Viruses", "V2 Viruses", "V3 Viruses"))+
  scale_y_continuous(labels = unit_format(unit = NULL, scale = 1e-6),
                     breaks= seq(-5e+6,13e+6, 2e+6),
                     limits = c(-5e+6, 14e+6))+
  theme_bw()+
  geom_errorbar(aes(ymin=mean_value + sd_value, ymax= mean_value - sd_value), width = 0.5, size = 1.0)+
  scale_x_continuous(breaks = c(0,3,6,17,20,24))+
  labs(title = 'NJ1 - Viral Production', subtitle = 'Overview - Bacterial and Viral counts for Lytic and Lysogenic inductions',
       x= 'Sampling Timepoints\n (in hours)', y='FCM Counts\n (in millions)')+
  theme(strip.text = element_text(face = "bold"),
        strip.background = element_rect( color = 'black', fill = 'white'),
        axis.title = element_text(face = 'bold'),
        title = element_text(face = 'bold'))
}

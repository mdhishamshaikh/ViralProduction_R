source("vp_functions.R")


#creating different time point plots against VP, VPC and Diff. Then make them separate for total, bacteria, and viruses
{

for (colnames in 2:length(TP)) {
  
  NJ1a$paste("'T",TP[1],"_T",TP[colnames], "'", sep = "") <- NA
}

NJ1a<- NJ1
NJ1a[colnames]<- NA #ADDING columns
}


TP<- unique(NJ1$Timepoint)
tp_char<- c()
for (tp in 1:length(TP)){
  a<- paste("'", TP[tp], "'", sep = "")
  tp_char[length(tp_char)+1]<- a
}

colnames<- c()
for (col in 2: length(TP)){
  a<- paste("T", TP[1], "_T", TP[col], sep = "")
  colnames[length(colnames)+1]<- a
}
timepoint<- c()
for (time in 2: length(TP)){
  a<- paste("T", TP[1], ":T", TP[time], sep = "")
  timepoint[length(timepoint)+1]<- a
}

paste("'", timepoint[2], "'", sep = "")

NJ1b<- NJ1%>%
  mutate("T0_T3" =  case_when(Timepoint ==tp_char[1] ~ timepoint[1],
                                 Timepoint == tp_char[2] ~ timepoint[1]))



NJ1a<- NJ1%>%
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

  








n<- ggplot(NJ1a, aes(x= Timepoint, y= mean_value, color= count , shape=count))+
  geom_point(size= 1.5)+
  geom_smooth(size= 1.0, method = 'lm', se = F)+
  geom_line()+
  geom_hline(yintercept = 0, color= '#636363', size= 0.3, linetype = "dashed")+
  facet_grid(factor(Subgroup, levels = c("Total", "Bacteria", "Viruses")) +
               factor(Sample_Type, levels = c("VP", "VPC", "Diff"))~
                factor(Time_Time, levels = unique(NJ1a$Time_Time))  
              )+
  scale_color_manual(name= 'Populations',
                     labels=c("Total Bacteria", "HNA Bacteria", "LNA Bacteria",
                              "Total Viruses", "V1 Viruses", "V2 Viruses", "V3 Viruses"),
                     values= c(c_Bacteria = "#00468B99", #AD002A99
                               c_HNA = "#00468B99",
                               c_LNA = "#0099B499",
                               c_Viruses = "#AD002A99",  #ED000099
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


o<- ggplot_gtable(ggplot_build(n))
strip_both<- which(grepl('strip-', o$layout$name))
fills<- c(rep("white", 14))
fills[3]<- "#E489A6" #Change the index to that of the time range after bacterial generation time
k <- 1

for (i in strip_both) {
  j<- which(grepl('rect', o$grobs[[i]]$grobs[[1]]$childrenOrder))
  o$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k<- k+1
}
#https://ojkftyk.blogspot.com/2019/01/r-ggplot2-change-colour-of-font-and.html
grid::grid.draw(o)


plots_lm_tp(NJ1a)



{
   library(tidyverse)
  library(ggsci)
  library(scales)
}


NJ1 <- as.data.frame(read_csv("NJ1.csv"))


#For VP, plot total bacteria, viruses, HNA, LNA, V1, V2, V3

VP<- NJ1[NJ1$Sample_Type=='VP',] %>%
  gather('c_Bacteria', 'c_HNA', 'c_LNA', 'c_Viruses', 'c_V1', 'c_V2', 'c_V3', key="count", value="value")%>%
  mutate(Microbe = if_else(VP$count == 'c_Bacteria' | VP$count == 'c_HNA' | VP$count == 'c_LNA', "Bacteria", "Viruses"))%>%
  mutate(Subgroup = if_else(VP$count == 'c_Bacteria' | VP$count == 'c_Viruses', "Parent", "Subgroup"))

#show_col(pal_lancet("lanonc")(9)) #color palette from ggsci::lancet
{
  count_factor_a<- factor(VP$count,
                          levels = c("c_Bacteria", "c_HNA", "c_LNA",
                                     "c_Viruses", "c_V1", "c_V2", "c_V3"))
  labels_a<- c("Total Bacteria", "HNA Bacteria", "LNA Bacteria",
               "Total Viruses", "V1 Viruses", "V2 Viruses", "V3 Viruses")
  a<- ggplot(VP, aes(x=Timepoint, y=value, color=count_factor_a, shape= count_factor_a)) +
    theme_minimal()+
    theme(plot.title = element_text(face = 'bold'),
          legend.title = element_text(face = 'bold'),
          axis.title = element_text(face = 'bold'))+
    geom_point()+
    scale_color_manual(name= 'Populations',
                       labels=labels_a,
                       values= c(c_Bacteria = "#00468BFF", 
                                 c_HNA = "#ED0000FF",
                                 c_LNA = "#42B540FF",
                                 c_Viruses = "#0099B4FF",
                                 c_V1 = "#925E9FFF",
                                 c_V2 = "#FDAF91FF",
                                 c_V3 = "#AD002AFF")) +
    scale_shape_manual(name = 'Populations', 
                       values = c(16,16,16,17,17,17,17),
                       labels=labels_a)+
    labs(title = 'NJ1 - Lytic Viral Production', subtitle = 'Overview - Bacterial and Viral counts',
         x= 'Sampling Timepoints\n (in hours)', y='FCM Counts \n (in millions)')+
    scale_y_continuous(labels = unit_format(unit = NULL, scale = 1e-6),
                       breaks= seq(0,13e+6, 4e+6),
                       limits = c(0, 14e+6))+
    scale_x_continuous(breaks = c(0,3,6,17,20,24))+
    stat_summary(fun= mean, geom = 'line', aes(group= count_factor_a),size= 1.0)+
    stat_summary(fun.min = function(x) mean(x) - sd(x), 
                 fun.max = function(x) mean(x) + sd(x), width = 0.5, size= 0.7,
                 geom = 'errorbar', aes(group = count_factor_a))
  print(a)
  }
{
  VP_sub_b<- VP[VP$Subgroup == 'Parent',]
  count_factor_b<- factor(VP_sub_b$count,
                          levels = c("c_Bacteria", "c_Viruses"))
  labels_b<- c("Total Bacteria", "Total Viruses")
  b<-  ggplot(VP_sub_b, aes(x=Timepoint, y=value, color=count_factor_b, shape= count_factor_b)) +
    theme_minimal()+
    theme(plot.title = element_text(face = 'bold'),
          legend.title = element_text(face = 'bold'),
          axis.title = element_text(face = 'bold'))+
    geom_point()+
    scale_color_manual(name= 'Populations',
                       labels=labels_b,
                       values= c(c_Bacteria = "#00468BFF",
                                 c_Viruses = "#0099B4FF"
                       )) +
    scale_shape_manual(name = 'Populations', 
                       values = c(16,16,17,17,17),
                       labels=labels_b)+
    labs(title = 'NJ1 - Lytic Viral Production - Parents', subtitle = 'Overview - Bacterial and Viral counts',
         x= 'Sampling Timepoints\n (in hours)', y='FCM Counts \n (in millions)')+
    scale_y_continuous(labels = unit_format(unit = NULL, scale = 1e-6),
                       breaks= seq(0,13e+6, 4e+6),
                       limits = c(0, 14e+6))+
    scale_x_continuous(breaks = c(0,3,6,17,20,24))+
    stat_summary(fun= mean, geom = 'line', aes(group= count_factor_b), size= 1.0)+
    stat_summary(fun.min = function(x) mean(x) - sd(x), 
                 fun.max = function(x) mean(x) + sd(x), width = 0.5, size= 0.7,
                 geom = 'errorbar', aes(group = count_factor_b))
  print(b)
}
{
  VP_sub_c<- VP[VP$Subgroup == 'Subgroup',]
  count_factor_c<- factor(VP_sub_c$count,
                          levels = c("c_HNA", "c_LNA","c_V1", "c_V2", "c_V3"))
  labels_c<- c("HNA Bacteria", "LNA Bacteria", "V1 Viruses", "V2 Viruses", "V3 Viruses")
  c<- ggplot(VP_sub_c, aes(x=Timepoint, y=value, color=count_factor_c, shape= count_factor_c)) +
    theme_minimal()+
    theme(plot.title = element_text(face = 'bold'),
          legend.title = element_text(face = 'bold'),
          axis.title = element_text(face = 'bold'))+
    geom_point()+
    scale_color_manual(name= 'Populations',
                       labels=labels_c,
                       values= c(c_HNA = "#ED0000FF",
                                 c_LNA = "#42B540FF",
                                 c_V1 = "#925E9FFF",
                                 c_V2 = "#FDAF91FF",
                                 c_V3 = "#AD002AFF")) +
    scale_shape_manual(name = 'Populations', 
                       values = c(16,16,17,17,17),
                       labels=labels_c)+
    labs(title = 'NJ1 - Lytic Viral Production - Subgroups', subtitle = 'Overview - Bacterial and Viral counts',
         x= 'Sampling Timepoints\n (in hours)', y='FCM Counts \n (in millions)')+
    scale_y_continuous(labels = unit_format(unit = NULL, scale = 1e-6),
                       breaks= seq(0,13e+6, 4e+6),
                       limits = c(0, 14e+6))+
    scale_x_continuous(breaks = c(0,3,6,17,20,24))+
    stat_summary(fun= mean, geom = 'line', aes(group= count_factor_c), size= 1.0)+
    stat_summary(fun.min = function(x) mean(x) - sd(x), 
                 fun.max = function(x) mean(x) + sd(x), width = 0.5, size= 0.7,
                 geom = 'errorbar', aes(group = count_factor_c))
  print(c)
}

plot_list<- list(a,b,c)
gridExtra::grid.arrange(a,b,c, ncol=3)(breaks = c(0,3,6,17,20,24))+
    stat_summary(fun= mean, geom = 'line', aes(group= count_factor_c), size= 1.0)+
    stat_summary(fun.min = function(x) mean(x) - sd(x), 
                 fun.max = function(x) mean(x) + sd(x), width = 0.5, size= 0.7,
                 geom = 'errorbar', aes(group = count_factor_c))
print(c)
}

plot_list<- list(a,b,c)
gridExtra::grid.arrange(a,b,c, ncol=3)


#Making a plot for VP and Diff with All, Parents and Subgroups as separate plots

#For VPC, plot total bacteria, viruses, HNA, LNA, V1, V2, V3

VPC<- NJ1[NJ1$Sample_Type=='VPC',]
VPC<-  gather(VPC,'c_Bacteria', 'c_HNA', 'c_LNA', 'c_Viruses', 'c_V1', 'c_V2', 'c_V3', key="count", value="value")
VPC<- mutate(VPC,Microbe = if_else(VPC$count == 'c_Bacteria' | VPC$count == 'c_HNA' | VPC$count == 'c_LNA', "Bacteria", "Viruses"))%>%
  mutate(Subgroup = if_else(VPC$count == 'c_Bacteria' | VPC$count == 'c_Viruses', "Parent", "Subgroup"))

#show_col(pal_lancet("lanonc")(9)) #color palette from ggsci::lancet
{
  count_factor_d<- factor(VPC$count,
                          levels = c("c_Bacteria", "c_HNA", "c_LNA",
                                     "c_Viruses", "c_V1", "c_V2", "c_V3"))
  labels_d<- c("Total Bacteria", "HNA Bacteria", "LNA Bacteria",
               "Total Viruses", "V1 Viruses", "V2 Viruses", "V3 Viruses")
  d<- ggplot(VPC, aes(x=Timepoint, y=value, color=count_factor_d, shape= count_factor_d)) +
    theme_minimal()+
    theme(plot.title = element_text(face = 'bold'),
          legend.title = element_text(face = 'bold'),
          axis.title = element_text(face = 'bold'))+
    geom_point()+
    scale_color_manual(name= 'Populations',
                       labels=labels_d,
                       values= c(c_Bacteria = "#00468BFF", 
                                 c_HNA = "#ED0000FF",
                                 c_LNA = "#42B540FF",
                                 c_Viruses = "#0099B4FF",
                                 c_V1 = "#925E9FFF",
                                 c_V2 = "#FDAF91FF",
                                 c_V3 = "#AD002AFF")) +
    scale_shape_manual(name = 'Populations', 
                       values = c(16,16,16,17,17,17,17),
                       labels=labels_d)+
    labs(title = 'NJ1 - Lytic + Lysogenic Viral Production', subtitle = 'Overview - Bacterial and Viral counts',
         x= 'Sampling Timepoints\n (in hours)', y='FCM Counts \n (in millions)')+
    scale_y_continuous(labels = unit_format(unit = NULL, scale = 1e-6),
                       breaks= seq(0,13e+6, 4e+6),
                       limits = c(0, 14e+6))+
    scale_x_continuous(breaks = c(0,3,6,17,20,24))+
    stat_summary(fun= mean, geom = 'line', aes(group= count_factor_d),size= 1.0)+
    stat_summary(fun.min = function(x) mean(x) - sd(x), 
                 fun.max = function(x) mean(x) + sd(x), width = 0.5, size= 0.7,
                 geom = 'errorbar', aes(group = count_factor_d))
  print(d)
  }
{
  VPC_sub_e<- VPC[VPC$Subgroup == 'Parent',]
  count_factor_e<- factor(VPC_sub_e$count,
                          levels = c("c_Bacteria", "c_Viruses"))
  labels_e<- c("Total Bacteria", "Total Viruses")
  e<-  ggplot(VPC_sub_e, aes(x=Timepoint, y=value, color=count_factor_e, shape= count_factor_e)) +
    theme_minimal()+
    theme(plot.title = element_text(face = 'bold'),
          legend.title = element_text(face = 'bold'),
          axis.title = element_text(face = 'bold'))+
    geom_point()+
    scale_color_manual(name= 'Populations',
                       labels=labels_e,
                       values= c(c_Bacteria = "#00468BFF",
                                 c_Viruses = "#0099B4FF"
                       )) +
    scale_shape_manual(name = 'Populations', 
                       values = c(16,16,17,17,17),
                       labels=labels_e)+
    labs(title = 'NJ1 - Lytic + Lysogenic Viral Production - Parents', subtitle = 'Overview - Bacterial and Viral counts',
         x= 'Sampling Timepoints\n (in hours)', y='FCM Counts \n (in millions)')+
    scale_y_continuous(labels = unit_format(unit = NULL, scale = 1e-6),
                       breaks= seq(0,13e+6, 4e+6),
                       limits = c(0, 14e+6))+
    scale_x_continuous(breaks = c(0,3,6,17,20,24))+
    stat_summary(fun= mean, geom = 'line', aes(group= count_factor_e), size= 1.0)+
    stat_summary(fun.min = function(x) mean(x) - sd(x), 
                 fun.max = function(x) mean(x) + sd(x), width = 0.5, size= 0.7,
                 geom = 'errorbar', aes(group = count_factor_e))
  print(e)
}
{
  VPC_sub_f<- VPC[VPC$Subgroup == 'Subgroup',]
  count_factor_f<- factor(VPC_sub_f$count,
                          levels = c("c_HNA", "c_LNA","c_V1", "c_V2", "c_V3"))
  labels_f<- c("HNA Bacteria", "LNA Bacteria", "V1 Viruses", "V2 Viruses", "V3 Viruses")
  f<- ggplot(VPC_sub_f, aes(x=Timepoint, y=value, color=count_factor_f, shape= count_factor_f)) +
    theme_minimal()+
    theme(plot.title = element_text(face = 'bold'),
          legend.title = element_text(face = 'bold'),
          axis.title = element_text(face = 'bold'))+
    geom_point()+
    scale_color_manual(name= 'Populations',
                       labels=labels_f,
                       values= c(c_HNA = "#ED0000FF",
                                 c_LNA = "#42B540FF",
                                 c_V1 = "#925E9FFF",
                                 c_V2 = "#FDAF91FF",
                                 c_V3 = "#AD002AFF")) +
    scale_shape_manual(name = 'Populations', 
                       values = c(16,16,17,17,17),
                       labels=labels_f)+
    labs(title = 'NJ1 - Lytic + Lysogenic Viral Production - Subgroups', subtitle = 'Overview - Bacterial and Viral counts',
         x= 'Sampling Timepoints\n (in hours)', y='FCM Counts \n (in millions)')+
    scale_y_continuous(labels = unit_format(unit = NULL, scale = 1e-6),
                       breaks= seq(0,13e+6, 4e+6),
                       limits = c(0, 14e+6))+
    scale_x_continuous(breaks = c(0,3,6,17,20,24))+
    stat_summary(fun= mean, geom = 'line', aes(group= count_factor_f), size= 1.0)+
    stat_summary(fun.min = function(x) mean(x) - sd(x), 
                 fun.max = function(x) mean(x) + sd(x), width = 0.5, size= 0.7,
                 geom = 'errorbar', aes(group = count_factor_f))
  print(f)
}


gridExtra::grid.arrange(a,b,c,d, e,f, ncol=3)



source("vp_functions.R", echo =T)

#Importing Data, calculating means,sd and Differnce values between VP and VPC
{#Works for with and without VPC
  NJ1<- read.csv("NJ1.csv")
NJ1<- NJ1[NJ1$Sample_Type != '0.22',] #Lose 0.22 values

NJ1<- gather(NJ1, 'c_Bacteria', 'c_HNA', 'c_LNA', 'c_Viruses', 'c_V1', 'c_V2', 'c_V3', key="count", value="value") %>%
  group_by(Location, Expt_No, Depth, Sample_Type, Timepoint, count ) %>%
   summarise(n =n(), mean=mean(value), sd=sd(value)) #calculating means and sd
if ('VPC' %in% NJ1$Sample_Type){
  colnames_mean<- c("VP_mean", "VPC_mean")
  colnames_sd<- c("VP_sd", "VPC_sd")
  
} else {
  colnames_mean<- c("VP_mean")
  colnames_sd<- c("VP_sd")
  
}

NJ1_mean<- NJ1[,1:8] %>% #splitting the dataframe cause I haven't figure out how to spread teh table without adding NAs
  spread('Sample_Type', 'mean')
if (length(colnames_mean)==2){
colnames(NJ1_mean)[7:8]<- colnames_mean
} else if (length(colnames_mean)==1){
  colnames(NJ1_mean)[7]<- colnames_mean
} 

NJ1_sd<- NJ1[,c(1:7,9)] %>%
  spread('Sample_Type', 'sd')
if (length(colnames_sd)==2){
  colnames(NJ1_sd)[7:8]<- colnames_sd
} else if (length(colnames_sd)==1){
  colnames(NJ1_sd)[7]<- colnames_sd
} 

if ('VPC' %in% NJ1$Sample_Type){
NJ1_mean$Diff_mean <- with(NJ1_mean, VPC_mean-VP_mean) #calcualting Diff mean
NJ1_mean<- pivot_longer(NJ1_mean, cols = c("VP_mean", "VPC_mean", "Diff_mean"), names_to= 'mean', values_to='mean_value')
NJ1_sd$Diff_sd <- with(NJ1_sd, VPC_sd+VP_sd) #Calculating Diff sd, whcih si addition of the other sds
NJ1_sd<- pivot_longer(NJ1_sd, cols = c("VP_sd", "VPC_sd", "Diff_sd"), names_to='sd', values_to= 'sd_value')
}

NJ1<- merge(NJ1_mean, NJ1_sd, by= c('Location', 'Expt_No', 'Depth',
                                    'Timepoint', 'count', 'n')) 


rm('NJ1_mean', 'NJ1_sd', 'colnames_mean', 'colnames_sd')

print("Viral Production means and standard deviations were calculated")






#Plotting function
#I want to give it a datafram ewith means and sd. This could be just for VP or VPC (including Diff)

NJ1<- NJ1%>% #Adding 'Microbe' and 'Subgroup' to the dataframe for ease. 
  mutate(Microbe = if_else(NJ1$count == 'c_Bacteria' | NJ1$count == 'c_HNA' | NJ1$count == 'c_LNA', "Bacteria", "Viruses"))%>%
  mutate(Subgroup = if_else(NJ1$count == 'c_Bacteria' | NJ1$count == 'c_Viruses', "Parent", "Child")) %>%
  arrange(Timepoint)


}


{ 
levels_total<- c("c_Bacteria", "c_HNA", "c_LNA",
                   "c_Viruses", "c_V1", "c_V2", "c_V3")
levels_parents<- c("c_Bacteria", "c_Viruses")
levels_child<- c("c_HNA", "c_LNA", "c_V1", "c_V2", "c_V3")

labels_total<- c("Total Bacteria", "HNA Bacteria", "LNA Bacteria",
                 "Total Viruses", "V1 Viruses", "V2 Viruses", "V3 Viruses")
labels_parents<- c("Total Bacteria", "Total Viruses")
labels_child<- c("HNA Bacteria", "LNA Bacteria", "V1 Viruses", "V2 Viruses", "V3 Viruses")

color_total<- c(c_Bacteria = "#00468BFF", 
               c_HNA = "#ED0000FF",
               c_LNA = "#42B540FF",
               c_Viruses = "#0099B4FF",
               c_V1 = "#925E9FFF",
               c_V2 = "#FDAF91FF",
               c_V3 = "#AD002AFF")
color_parents<- color_total[c(1,4)]
color_child<- color_total[c(2,3,5:7)]

shape_total<- c(16,16,16,17,17,17,17)
shape_parents<- c(16,17)
shape_child<- c(16,16,17,17,17)
}

{
  plot_list<- list()
  
for (group in c('Total', 'Parent', 'Child')){
  if (group == 'Total') {
    df<- NJ1
    count_factor <- factor(df$count, levels = levels_total)
    labels <- labels_total
    color <- color_total
    shape<- shape_total
    } else if (group == 'Parent') {
      df<- NJ1[NJ1$Subgroup== 'Parent',]
      count_factor <- factor(df$count, levels = levels_parents)
      labels <- labels_parents
      color <- color_parents
      shape<- shape_parents
      } else if (group == 'Child') {
        df<- NJ1[NJ1$Subgroup== 'Child',]
        count_factor <- factor(df$count, levels = levels_child)
        labels <- labels_child
        color <- color_child
        shape<- shape_child
      }
ggplot(df, aes(x=Timepoint, y=VP_mean, color=count_factor, shape= count_factor)) +
    theme_minimal()+
    theme(plot.title = element_text(face = 'bold'),
          legend.title = element_text(face = 'bold'),
          axis.title = element_text(face = 'bold'))+
    geom_point(size=2.0)+
    scale_color_manual(name= 'Populations',
                       labels=labels,
                       values= color) +
    scale_shape_manual(name = 'Populations', 
                       values = shape,
                       labels=labels)+
    labs(title = 'NJ1 - Lytic Viral Production', subtitle = 'Overview - Bacterial and Viral counts',
         x= 'Sampling Timepoints\n (in hours)', y='FCM Counts \n (in millions)')+
    scale_y_continuous(labels = unit_format(unit = NULL, scale = 1e-6),
                       breaks= seq(0,13e+6, 4e+6),
                       limits = c(0, 14e+6))+
    scale_x_continuous(breaks = c(0,3,6,17,20,24)) +
    geom_line(size= 1.0) +
    geom_errorbar(aes(ymin=VP_mean + VP_sd, ymax= VP_mean - VP_sd), width = 0.5, size = 0.7)
 
ggsave(filename = paste("myplot",group,".png",sep=""))
}
  #gridExtra::grid.arrange(grobs= plot_list, ncol=3) 
  
  }#choosing the correct variable. done.

filenames<-dir(pattern='myplot')
foo<-list()
for(j in 1:3) foo[[j]]<-png::readPNG(filenames[j])
                                 
layout(matrix(1:3,nr=1,byr=T))
for (j in 1:3) plot(foo[[j]])                                
                                 
count_factor <- factor(df$count, levels = levels_total)

for  (group in c('Total', 'Parent', 'Child')){
color= c(factor(NJ1$count, levels = if (group == 'Total') {
  print(levels_total)
} else if (group == 'Parent') {
  print(levels_parents)
} else if (group == 'Child') {
  print(levels_child)
}
  ))
}




#Defining within vriable. ew.

{plot_list<- c()
for (group in c('Total', 'Parent', 'Child')){
plot<- ggplot(data = if (group == 'Total') {
  print(NJ1)
} else if (group == 'Parent') {
  NJ1[NJ1$Subgroup== 'Parent',]
} else if (group == 'Child') {
  NJ1[NJ1$Subgroup== 'Child',]
}, 
       aes(x = Timepoint, 
           y = VP_mean,
           color = factor(count, levels = if (group == 'Total') {
                 print(levels_total)
               } else if (group == 'Parent') {
                 print(levels_parents)
               } else if (group == 'Child') {
                 print(levels_child)
               }
               ),
           shape= factor(count, levels = if (group == 'Total') {
                 print(levels_total)
               } else if (group == 'Parent') {
                 print(levels_parents)
               } else if (group == 'Child') {
                 print(levels_child)
               }
               ))) +
  theme_minimal() +
  theme(plot.title = element_text(face = 'bold'),
        legend.title = element_text(face = 'bold'),
        axis.title = element_text(face = 'bold'))+
  geom_point(size=2.0)+
  scale_color_manual(name= 'Populations',
                     labels = if (group == 'Total') {
                       print(labels_total)
                     } else if (group == 'Parent') {
                       print(labels_parents)
                     } else if (group == 'Child') {
                       print(labels_child)
                     },
                     values = if (group == 'Total') {
                       print(color_total)
                     } else if (group == 'Parent') {
                       print(color_parents)
                     } else if (group == 'Child') {
                       print(color_child)
                     }) +
  scale_shape_manual(name = 'Populations', 
                     labels = if (group == 'Total') {
                       print(labels_total)
                     } else if (group == 'Parent') {
                       print(labels_parents)
                     } else if (group == 'Child') {
                       print(labels_child)
                     },
                     values = if (group == 'Total') {
                       print(shape_total)
                     } else if (group == 'Parent') {
                       print(shape_parents)
                     } else if (group == 'Child') {
                       print(shape_child)
                     })+
  labs(title = 'NJ1 - Lytic Viral Production', subtitle = 'Overview - Bacterial and Viral counts',
       x= 'Sampling Timepoints\n (in hours)', y='FCM Counts \n (in millions)')+
  scale_y_continuous(labels = unit_format(unit = NULL, 
                                          scale = 1e-6),
                     breaks= seq(0,13e+6, 4e+6),
                     limits = c(0, 14e+6))+
  scale_x_continuous(breaks = c(0,3,6,17,20,24)) +
  geom_line(size= 1.0) +
  geom_errorbar(aes(ymin=VP_mean + VP_sd, 
                    ymax= VP_mean - VP_sd), 
                width = 0.5, 
                size = 0.7)
print(plot)
plot_list[[length(plot_list)+1]]<- plot
}
}

gridExtra::grid.arrange(grobs= plot_list, ncol=3)

















                             

plots_vp<- function(VP, VPC=NULL, Diff=NULL) {
  
count_factor<- factor(NJ1$count,
                      levels = c("c_Bacteria", "c_HNA", "c_LNA",
                                 "c_Viruses", "c_V1", "c_V2", "c_V3"))
labels<- c("Total Bacteria", "HNA Bacteria", "LNA Bacteria",
           "Total Viruses", "V1 Viruses", "V2 Viruses", "V3 Viruses")
ggplot(NJ1, aes(x=Timepoint, y=VP_mean, color=count_factor, shape= count_factor)) +
  theme_minimal()+
  theme(plot.title = element_text(face = 'bold'),
        legend.title = element_text(face = 'bold'),
        axis.title = element_text(face = 'bold'))+
  geom_point(size=2.0)+
  scale_color_manual(name= 'Populations',
                     labels=labels,
                     values= c(c_Bacteria = "#00468BFF", 
                               c_HNA = "#ED0000FF",
                               c_LNA = "#42B540FF",
                               c_Viruses = "#0099B4FF",
                               c_V1 = "#925E9FFF",
                               c_V2 = "#FDAF91FF",
                               c_V3 = "#AD002AFF")) +
  scale_shape_manual(name = 'Populations', 
                     values = c(16,16,16,17,17,17,17),
                     labels=labels)+
  labs(title = 'NJ1 - Lytic Viral Production', subtitle = 'Overview - Bacterial and Viral counts',
       x= 'Sampling Timepoints\n (in hours)', y='FCM Counts \n (in millions)')+
  scale_y_continuous(labels = unit_format(unit = NULL, scale = 1e-6),
                     breaks= seq(0,13e+6, 4e+6),
                     limits = c(0, 14e+6))+
  scale_x_continuous(breaks = c(0,3,6,17,20,24)) +
  geom_line(size= 1.0) +
  geom_errorbar(aes(ymin=VP_mean + VP_sd, ymax= VP_mean - VP_sd), width = 0.5, size = 0.7)
}


  # stat_summary(fun= mean, geom = 'line', aes(group= count_factor),size= 1.0)+
  # stat_summary(fun.min = function(x) mean(x) - sd(x), 
  #              fun.max = function(x) mean(x) + sd(x), width = 0.5, size= 0.7,
  #              geom = 'errorbar', aes(group = count_factor))




ggplot(NJ1, aes(x= Timepoint, y= mean_value, color= count , shape=count))+
  geom_point(size= 2.0)+
  geom_line(size= 1.0)+
  facet_grid(Subgroup~mean, space= "free_y")+
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
  theme_minimal()+
  geom_errorbar(aes(ymin=mean_value + sd_value, ymax= mean_value - sd_value), width = 0.5, size = 1.0)+
  scale_x_continuous(breaks = c(0,3,6,17,20,24))+
  labs(title = 'NJ1 - Viral Production', subtitle = 'Overview - Bacterial and Viral counts for Lytic and Lysogenic inductions',
       x= 'Sampling Timepoints\n (in hours)', y='FCM Counts')




#Importing Data, calculating means,sd and Differnce values between VP and VPC
{
  NJ1<- read.csv("NJ1.csv")
NJ1<- NJ1[NJ1$Sample_Type != '0.22',] #Lose 0.22 values

NJ1<- gather(NJ1, 'c_Bacteria', 'c_HNA', 'c_LNA', 'c_Viruses', 'c_V1', 'c_V2', 'c_V3', key="count", value="value") %>%
  group_by(Location, Expt_No, Depth, Sample_Type, Timepoint, count ) %>%
   summarise(n =n(), mean=mean(value), sd=sd(value)) #calculating means and sd

NJ1_mean<- NJ1[,1:8] %>% #splitting the datfram cause I haven't figure out how to spread teh table without adding NAs
  spread('Sample_Type', 'mean')
colnames(NJ1_mean)[7:8]<- c("VP_mean", "VPC_mean")
NJ1_mean$Diff_mean <- with(NJ1_mean, VPC_mean-VP_mean) #calcualting Diff mean

NJ1_sd<- NJ1[,c(1:7,9)] %>%
  spread('Sample_Type', 'sd')
colnames(NJ1_sd)[7:8]<- c("VP_sd", "VPC_sd")
NJ1_sd$Diff_sd <- with(NJ1_sd, VPC_sd+VP_sd) #Calculating Diff sd, whcih si addition of the other sds


NJ1<- merge(NJ1_mean, NJ1_sd, by= c('Location', 'Expt_No', 'Depth',
                                    'Timepoint', 'count', 'n')) 
}

#Plotting function

count_factor<- factor("c_Bacteria", "c_HNA", "c_LNA", "c_V1", "c_V2", "c_V3", )
labels<- c("Total Bacteria", "HNA Bacteria", "LNA Bacteria",
           "Total Viruses", "V1 Viruses", "V2 Viruses", "V3 Viruses")
ggplot(NJ1, aes(x=Timepoint, y=VP_mean, color=count_factor, shape= count_factor)) +
  theme_minimal()+
  theme(plot.title = element_text(face = 'bold'),
        legend.title = element_text(face = 'bold'),
        axis.title = element_text(face = 'bold'))+
  geom_point()+
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
  scale_x_continuous(breaks = c(0,3,6,17,20,24))

  # stat_summary(fun= mean, geom = 'line', aes(group= count_factor),size= 1.0)+
  # stat_summary(fun.min = function(x) mean(x) - sd(x), 
  #              fun.max = function(x) mean(x) + sd(x), width = 0.5, size= 0.7,
  #              geom = 'errorbar', aes(group = count_factor))

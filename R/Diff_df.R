source("vp_functions.R", echo =T)

#Importing Data, calculating means,sd and Differnce values between VP and VPC
{#Works for with and without VPC
  NJ1<- read.csv("NJ1.csv")
NJ1<- NJ1[NJ1$Sample_Type != '0.22',] #Lose 0.22 values

NJ1<- gather(NJ1, 'c_Bacteria', 'c_HNA', 'c_LNA', 'c_Viruses', 'c_V1', 'c_V2', 'c_V3', key="count", value="value") %>%
  group_by(Location, Expt_No, Depth, Sample_Type, Timepoint, count ) %>%
   summarise(n =n(), mean=mean(value), sd=sd(value)) #calculating means and sd
if ('VPC' %in% NJ1$Sample_Type){
  colnames_mean<- c("VP", "VPC")
  colnames_sd<- c("VP", "VPC")
  
} else {
  colnames_mean<- c("VP")
  colnames_sd<- c("VP")
  
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
NJ1_mean$Diff <- with(NJ1_mean, VPC-VP) #calcualting Diff mean
NJ1_mean<- pivot_longer(NJ1_mean, cols = c("VP", "VPC", "Diff"), names_to= 'Sample_Type', values_to='mean_value')
NJ1_sd$Diff <- with(NJ1_sd, VPC+VP) #Calculating Diff sd, whcih si addition of the other sds
NJ1_sd<- pivot_longer(NJ1_sd, cols = c("VP", "VPC", "Diff"), names_to='Sample_Type', values_to= 'sd_value')
}

NJ1<- merge(NJ1_mean, NJ1_sd, by= c('Location', 'Expt_No', 'Depth',
                                    'Timepoint', 'count', 'n', 'Sample_Type')) 


rm('NJ1_mean', 'NJ1_sd', 'colnames_mean', 'colnames_sd')

print("Viral Production means and standard deviations were calculated")






#Plotting function
#I want to give it a dataframe with means and sd. This could be just for VP or VPC (including Diff)

NJ1<-  mutate(NJ1, Subgroup = case_when(NJ1$count == 'c_Bacteria' ~ "Total",
                              NJ1$count == 'c_Viruses' ~ "Total",
                              NJ1$count == 'c_HNA' ~ "Bacteria",
                              NJ1$count == 'c_LNA' ~ "Bacteria",
                              NJ1$count == 'c_V1' ~ "Viruses",
                              NJ1$count == 'c_V2' ~ "Viruses",
                              NJ1$count == 'c_V3' ~ "Viruses")) %>%
  arrange(Timepoint)



}


plots_vipcal(NJ1)
plots_lm(NJ1)
  
ggsave("NJ1.png")



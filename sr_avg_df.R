source("vp_functions.R")


NJ1<- read.csv("NJ1.csv")

#separate replicate dataframe
NJ1_sr<- NJ1%>%
  gather('c_Bacteria', 'c_HNA', 'c_LNA', 'c_Viruses', 'c_V1', 'c_V2', 'c_V3', key="count", value="value")%>%
  mutate(Microbe = if_else(count == 'c_Bacteria' | count == 'c_HNA' | count == 'c_LNA', "Bacteria", "Viruses"))%>%
  mutate
NJ1_sr<- NJ1_sr[NJ1_sr$Sample_Type != '0.22',] #removing 0.22 for now

#average replicate dataframe
NJ1_av<- gather(NJ1, 'c_Bacteria', 'c_HNA', 'c_LNA', 'c_Viruses', 'c_V1', 'c_V2', 'c_V3', key="count", value="value") %>%
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

NJ1<- NJ1[NJ1$Sample_Type != '0.22',] #Lose 0.22 values
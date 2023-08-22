###IMPORTANT
####Calculates bacterial production rates and generation time


{
  library(tidyverse)
  
}


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
                                      'Timepoint', 'count', 'n')) %>%
    mutate(Microbe = if_else(count == 'c_Bacteria' | count == 'c_HNA' | count == 'c_LNA', "Bacteria", "Viruses"))%>%
    mutate(Subgroup = if_else(count == 'c_Bacteria' | count == 'c_Viruses', "Parent", "Subgroup"))
  
  rm(NJ1_mean)
  rm(NJ1_sd)
}



Ba_gt<- function(x){
  GT<- (log10(2)*(Ba$Timepoint[x]-Ba$Timepoint[1]))/(log10(Ba$VP_mean[x])-log10(Ba$VP_mean[1]))
  print(GT) 
}

plot<- c()

for (x in c('c_Bacteria', 'c_HNA', 'c_LNA')){
  Ba<- NJ1[NJ1$count== x,] %>%
    arrange(Timepoint)
  
  for (i in 2:6){
    y<-Ba_gt(i)
   plot[length(plot)+1]<- y
   if (y > 24){
     print("Low bacterial production") }
   if (y < 0){
     print("Low bacterial production") }
#   plot(plot, x=c(3,6,17,20,24))
# abline(h=24)
# abline(h=48, col=2)

}
}

#Interesting! LNA bacterial growth isn't that significant

Bacterial_GT <- plot[1:5]
HNA_GT <- plot[6:10]
LNA_GT <- plot[11:15]
bp_df<- data.frame(Bacterial_GT, HNA_GT, LNA_GT) #in hours
bp_endpoint<- intersect(which(Bacterial_GT>0), which(Bacterial_GT<24))[1] #use this as the index for highlighting on plot
bp<- Bacterial_GT[bp_endpoint]
intersect(which(HNA_GT>0), which(HNA_GT<24))[1]
intersect(which(LNA_GT>0), which(LNA_GT<24))[1]



#Add the bacterial generation time and production rate along with time in a tibble or dataframe

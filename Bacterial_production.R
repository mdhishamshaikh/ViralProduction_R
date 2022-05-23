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
                                      'Timepoint', 'count', 'n')) 
  
  rm(NJ1_mean)
  rm(NJ1_sd)
}

#only selecting total bacteria here.
Ba<- NJ1[NJ1$count== 'c_Bacteria',]
plot(x = Ba$Timepoint, y = log10(Ba$VP_mean))

Ba<-arrange(Ba,Timepoint)
#Generation time

#Generation time = (log10(2) *(time2-time1))/(log10(bacteria2)-log10(bacteria1))

Ba_gt<- function(x){
  GT<- (log10(2)*(Ba$Timepoint[x]-Ba$Timepoint[1]))/(log10(Ba$VP_mean[x])-log10(Ba$VP_mean[1]))
  print(GT) 
}

for (i in 1:6){
  y<-Ba_gt(i)
  plot[[length(plot)+1]]<- y
}



for (x in c('c_Bacteria', 'c_HNA', 'c_LNA')){
  Ba<- NJ1[NJ1$count== x,] %>%
    arrange(Timepoint)
  plot<- list()
  for (i in 2:6){
    y<-Ba_gt(i)
   plot[[length(plot)+1]]<- y
   if (0< y > 24){
     print("High bacterial production") }
#   plot(plot, x=c(3,6,17,20,24))
# abline(h=24)
# abline(h=48, col=2)

}
}

#Interesting! LNA bacterial growth isn't that significant




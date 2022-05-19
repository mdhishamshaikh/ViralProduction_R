{
  library(readxl)
  library(quantmod)
  library(magrittr)
  library(ggplot2)
  library(ggpubr)
  library(tidyverse)
  library(gridExtra)
}

{
data<- as.data.frame(read_excel("viralproduction_R.xlsx")) 

  #VIPCAL works on mean values only. So let's calculate the mean.

data$VP_mean <- apply(data[,2:4], 1, mean)
data$VPC_mean <- apply(data[,5:7], 1, mean)
data$Diff_mean <- apply(data[,8:10], 1, mean)
data<- rbind (c(-1e+2, rep(1e+10, times=12)), data)
data<- rbind (data,c(1e+2, rep(-1e+10, times=12)))
data$VP_sd <- apply(data[,2:4], 1, sd)
data$VPC_sd <- apply(data[,5:7], 1, sd)
data$Diff_sd <- apply(data[,8:10], 1, sd)
}

#VIPCAL Viral Production without standard deviation

viralp<- c()
for ( i in 11:12){
  p<- quantmod::findPeaks(data[,i])-1
  print(p)
  v<- quantmod::findValleys(data[,i])-1
  print(v)
  
  print (i)
  
  if (identical(length(p),length(q))) {
    print("great")
  }
  
  if (length(p)==1) {
    vp<- (data[p[1],i] - data[v[1],i])/(data[p[1],1] - data[v[1],1])
    
  } else if (length(p)==2) {
    vp<- ((data[p[1],i] - data[v[1],i])/(data[p[1],1] - data[v[1],1]) + 
            (data[p[2],i] - data[v[2],2])/(data[p[2],i] - data[v[2],1]))/2
    
  } else if (length(p)==3) {
    vp<-  ((data[p[1],i] - data[v[1],i])/(data[p[1],1] - data[v[1],1]) + 
             (data[p[2],i] - data[v[2],i])/(data[p[2],1] - data[v[2],1]) +
             (data[p[3],i] - data[v[1],i])/(data[p[3],1] - data[v[3],1]))/3
  }
  
  print(vp)
  viralp = c(viralp, vp)
}
viralproduction_VP<- viralp[1]
viralproduction_Diff<- viralp[2]

#Here you can see that there are time points overlapping in Diff VIPCAL plot
#Therefore it is important to consider the standard deviation of the points to 
#select for peaks and valleys

ggplot(data[2:7,], aes(Timepoints, VP_mean))+
  geom_point()+
  #geom_smooth()+
geom_line()+
  theme_minimal()+
  geom_errorbar(aes(ymin= VP_mean-VP_sd, ymax=VP_mean+VP_sd), width=.2,
                position=position_dodge(.9))
ggplot(data[2:7,], aes(Timepoints, VPC_mean))+
  geom_point()+
  #geom_smooth()+
  geom_line()+
  theme_minimal()+
  geom_errorbar(aes(ymin= VPC_mean-VPC_sd, ymax=VPC_mean+VPC_sd), width=.2,
                position=position_dodge(.9))
ggplot(data[2:7,], aes(Timepoints, Diff_mean))+
  geom_point()+
  geom_line()+
  theme_minimal()+
  geom_errorbar(aes(ymin= Diff_mean-Diff_sd, ymax=Diff_mean+Diff_sd), width=.2,
                position=position_dodge(.9))




#Trying to see if i could find peaks and valleys

peaks(data$Diff_mean, data$Diff_sd)
valleys(data$Diff_mean, data$Diff_sd)


# works!!!!!!! lol
peaks(data$VP_mean, data$VP_sd)
valleys(data$VP_mean, data$VP_sd)


#Need to work on this. figure out why there are more valleys

viralp<- c()
for ( i in c(11,13)){
  p<- peaks(data[,i])-1
  print(p)
  v<- quantmod::findValleys(data[,i])-1
  print(v)
  
  print (i)
  
  if (identical(length(p),length(q))) {
    print("great")
  }
  
  if (length(p)==1) {
    vp<- (data[p[1],i] - data[v[1],i])/(data[p[1],1] - data[v[1],1])
    
  } else if (length(p)==2) {
    vp<- ((data[p[1],i] - data[v[1],i])/(data[p[1],1] - data[v[1],1]) + 
            (data[p[2],i] - data[v[2],2])/(data[p[2],i] - data[v[2],1]))/2
    
  } else if (length(p)==3) {
    vp<-  ((data[p[1],i] - data[v[1],i])/(data[p[1],1] - data[v[1],1]) + 
             (data[p[2],i] - data[v[2],i])/(data[p[2],1] - data[v[2],1]) +
             (data[p[3],i] - data[v[1],i])/(data[p[3],1] - data[v[3],1]))/3
  }
  
  print(vp)
  viralp = c(viralp, vp)
}
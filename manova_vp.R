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


#We are trying to attempt manova with bonferroni's correction to see if there 
#is a difference between the time points or not


#let's try a few t-test first

#see if the difference between t0 and t3 in VP is significant or not. 
VP<- data[,1:4]
VP<- VP[-c(1,8),]

for (i in 1:5){
  print(i)
  print(t.test(VP[i,2:4], VP[i+1,2:4]))
}


Diff<- data[,c(1, 8:10)]
Diff<- Diff[-c(1,8),]

for (i in 1:5){
  print(i)
  print(t.test(Diff[i,2:4], Diff[i+1,2:4], paired = T))
}


source("vp_functions.R")
NJ1<- read.csv("NJ1.csv")


#1Linear Regression

#1.1 Separate Replicates
#With the help of separate replicates, I can calculate viral production rates
#for all time point ranges, with three slopes, and then average them. 

NJ1_sr_slope<- df_sr_tp(NJ1) %>%
  slope_lm_sr()%>%
  df_avg_slope()


#1.2 Averaged Replicates #JUST WIRTING THIS HERE CAUSE i dont know wherelese to write it. 
#I can't do separate replicates for VIPCAL-model, cause we need to avaergae the replicates

NJ1_avg_slope<- df_avg_tp(NJ1)%>%
  slope_lm_avg()
#don't know what to do about the pre-calculated standard deviations 

#2 VIPCAL

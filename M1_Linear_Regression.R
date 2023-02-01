####0.0 Set Up####
library("tidyverse")
setwd(".")
source("vp_functions.R")
source("0_vp_source.R")

####1.0 Importing Data####
NJ1<- read.csv("NJ1.csv")
unique(NJ1$Sample_Type)
#for now we'll not use the 0.22 samples.

NJ1<- NJ1[NJ1$Sample_Type != "0.22",]
unique(NJ1$Sample_Type)

####2.0 Creating Time Range dataframe####

NJ1<- overview_df_tp_avg(NJ1) #this also calculates difference per timepoint. Therefore the diff curve, which we don't need. We'll exclude it later.
#This diff curve is calculated using averages per timepoint, where we could also apply LMER. We'll do this as a separate method.
#Also, and average replicate dataframe

####3.0 Calculating Slopes####
slopes_NJ1<- slope_lm_avg(NJ1)
slopes_NJ1$Time_Range<- factor(slopes_NJ1$Time_Range, levels = c("T0_T3", "T0_T6", "T0_T17", "T0_T20", "T0_T24"))
slopes_NJ1$LM_AVG_Slope<- as.numeric(slopes_NJ1$LM_AVG_Slope)
#We also need to calculate lysogeny slope now.
#This is done by subtracting VP and VPC per time range per population

Diff_Slopes_NJ1<- slopes_NJ1[slopes_NJ1$Sample_Type == "Diff",]
diff<- spread(Diff_Slopes_NJ1, key = Sample_Type, value = LM_AVG_Slope)

ggplot(data = Diff_Slopes_NJ1, aes(x = Time_Range, y = LM_AVG_Slope)) +
  geom_point()+
  facet_grid(factor(Diff_Slopes_NJ1$Population, levels = c("c_Viruses", "c_V1", "c_V2", "c_V3")))

#calcualting Lysogeny by subtarcting average slopes

lysogeny_df<- slopes_NJ1[slopes_NJ1$Sample_Type != "Diff",] %>%
  spread(key = Sample_Type, value = LM_AVG_Slope)%>%
  mutate(lysogeny = VPC - VP )

ggplot(data = lysogeny_df, aes(x = Time_Range, y = lysogeny)) +
  geom_point()+
  facet_grid(factor(lysogeny_df$Population, levels = c("c_Viruses", "c_V1", "c_V2", "c_V3")))

NJ1_slopes_merge<- left_join(lysogeny_df, diff, byy = c("Time_Range", "Population"))
  
plot(Diff~ lysogeny, data = NJ1_slopes_merge)
abline(0,1)
#redundant!!!!
NJ1_mdf<- read.csv("NJ1.csv") %>%
  df_sr()%>%
  subset( Microbe == "Viruses" & Sample_Type != '0.22')
model_Data<- lmer_model(NJ1_mdf)
model_slope<- slope_lm_avg(model_Data)



#Difference Curve Calculation

tp_df<- tp(df_avg(NJ1))



####Starting again####
####1.0 Linear Regression - Lytic####

library(tidyverse)

data<- read.csv("NJ1.csv")
dim(data)

#1.1 Slope through all the points
data_sr<- df_sr_tp(data)
slopes_LM_AP_lytic<- slope_all_points(data_sr)


#1.2 Slope through every replicate
#we use data_sr for this
slopes_LM_SR_lytic<- slope_lm_sr(data_sr)
slopes_LM_SR_lytic_avg<- slopes_LM_SR_lytic %>%
  group_by(Location, Expt_No, Depth, Sample_Type, Time_Range, Population ) %>%
  summarise(LM_SR_Slope_mean=mean(as.numeric(LM_SR_Slope)), LM_SR_Slope_sd=sd(as.numeric((LM_SR_Slope))))


#we'll have to avergae these slopes then.



#1.3 slope through averaged replicates
data_avg<- df_avg_tp(data)
slopes_LM_AR_lytic<- slope_lm_avg(data_avg)
#this function also calculates the Differences


####2.O Linear Regression - Lysogeny

#2.1 Slope through all point
slopes_LM_AP_lysogenic<- slopes_LM_AP_lytic

  diff<- inner_join(filter(slopes_LM_AP_lytic, Sample_Type == "VPC"), filter(slopes_LM_AP_lytic, Sample_Type == "VP"), by = c("Location", "Expt_No",
                                                                                                                   "Depth", "Time_Range", 
                                                                                                                   "Population")) %>%
  mutate(Sample_Type = "Diff",
          LM_AP_Slope = as.numeric(LM_AP_Slope.x) - as.numeric(LM_AP_Slope.y),
         LM_AP_SE = as.numeric(LM_AP_SE.x) - as.numeric(LM_AP_SE.y),
         LM_AP_R_squared = NA)%>%
  select("Location", "Expt_No",
         "Depth", "Time_Range", 
         "Population", "Sample_Type",
         "LM_AP_Slope", "LM_AP_SE", "LM_AP_R_squared")%>%
    rbind(slopes_LM_AP_lytic)%>%
    arrange(Location, Expt_No,
          Depth, Time_Range, 
          Population, factor(Sample_Type, levels = c("VP", "VPC", "Diff")))
    
ggplot(data = diff, aes(y = as.numeric(LM_AP_R_squared)))+
  geom_boxplot()+
  facet_grid(Sample_Type ~ Population, scales = "fixed")
  
#2.2 Slope trhough every replciate

slopes_LM_SR_lysogenic<- slopes_LM_SR_lytic
slopes_LM_SR_lysogenic_avg<- slopes_LM_SR_lysogenic %>%
  group_by(Location, Expt_No, Depth, Sample_Type, Time_Range, Population ) %>%
  summarise(LM_SR_Slope_mean=mean(as.numeric(LM_SR_Slope)), LM_SR_Slope_sd=sd(as.numeric((LM_SR_Slope))))

diff_sr<- inner_join(filter(slopes_LM_SR_lysogenic_avg, Sample_Type == "VPC"), filter(slopes_LM_SR_lysogenic_avg, Sample_Type == "VP"), by = c("Location", "Expt_No",
                                                                                                                     "Depth", "Time_Range", 
                                                                                                                     "Population")) %>%
  mutate(Sample_Type = "Diff",
         LM_SR_Slope_ = as.numeric(LM_AVG_Slope) - as.numeric(LM_SR_Slope_mean.y),
         LM_SR_Slope_sd = as.numeric(LM_AVG_SE.x) - as.numeric(LM_SR_Slope_sd.y),
         LM_SR_R_squared = NA)%>%
  select("Location", "Expt_No",
         "Depth", "Time_Range", 
         "Population", "Sample_Type",
         "LM_SR_Slope_mean", "LM_SR_Slope_sd", "LM_SR_R_squared")%>%
  rbind(slopes_LM_SR_lysogenic_avg)%>%
  arrange(Location, Expt_No,
          Depth, Time_Range, 
          Population, factor(Sample_Type, levels = c("VP", "VPC", "Diff")))


ggplot(data = diff_sr, aes(y = as.numeric(LM_SR_Slope_mean)))+
  geom_boxplot()+
  facet_grid(Sample_Type ~ Population, scales = "free")


#2.3 Slope through averaged replciates

slopes_LM_AR_lysogenic<- slopes_LM_AR_lytic

diff_ar<- inner_join(filter(slopes_LM_AR_lysogenic, Sample_Type == "VPC"), filter(slopes_LM_AR_lysogenic, Sample_Type == "VP"), by = c("Location", "Expt_No",
                                                                                                                                               "Depth", "Time_Range", 
                                                                                                                                               "Population")) %>%
  mutate(Sample_Type = "Diff",
         LM_AVG_Slope = as.numeric(LM_AVG_Slope.x) - as.numeric(LM_AVG_Slope.y),
         LM_AVG_SE = as.numeric(LM_AVG_SE.x) - as.numeric(LM_AVG_SE.y),
         LM_Avg_R_squared = NA)%>%
  select("Location", "Expt_No",
         "Depth", "Time_Range", 
         "Population", "Sample_Type",
         "LM_AVG_Slope", "LM_AVG_SE", "LM_Avg_R_squared")%>%
  rbind(slopes_LM_AR_lysogenic)%>%
  arrange(Location, Expt_No,
          Depth, Time_Range, 
          Population, factor(Sample_Type, levels = c("VP", "VPC", "Diff")))











#combine the dataframes 

merge_df<- merge(slopes_LM_AP_lytic, slopes_LM_SR_lytic_avg, by.x = c("Location", "Expt_No", "Depth", "Time_Range", "Popualtion"),
                 by.y = c("Location", "Expt_No", "Depth", "Time_Range", "Popualtion"))

join<- left_join(slopes_LM_AP_lytic, slopes_LM_SR_lytic_avg, by = c("Location", "Expt_No", "Depth", "Time_Range", "Population", "Sample_Type"))
#1.0 vs 2.0

par(mfrow = c(4,3))


for(tim in unique(slopes_LM_AR_lytic$Time_Range)){
for(st in c("VP", "VPC")){
for (vir in unique(slopes_LM_AR_lytic$Population)) {
boxplot(as.numeric(slopes_LM_AP_lytic[slopes_LM_AP_lytic$Population == vir & slopes_LM_AP_lytic$Sample_Type == st , ]$LM_AP_R_squared), xlab = paste0(tim, vir, st) )
boxplot(as.numeric(slopes_LM_SR_lytic[slopes_LM_SR_lytic$Population == vir & slopes_LM_SR_lytic$Sample_Type == st,]$LM_SR_R_squared), xlab = paste0(tim, vir, st))
boxplot(as.numeric(slopes_LM_AR_lytic[slopes_LM_AR_lytic$Population == vir & slopes_LM_AR_lytic$Sample_Type == st,]$LM_Avg_R_squared), xlab = paste0(tim, vir, st))
}
}
}

par(mfrow = c(4,3))
for(tim in unique(slopes_LM_AR_lytic$Time_Range)){
  for(st in c("VP", "VPC")){
    for (vir in unique(slopes_LM_AR_lytic$Population)) {
      boxplot(as.numeric(slopes_LM_AP_lytic[slopes_LM_AP_lytic$Population == vir & slopes_LM_AP_lytic$Sample_Type == st , ]$LM_AP_Slope), xlab = paste0(tim, vir, st) )
      boxplot(as.numeric(slopes_LM_SR_lytic[slopes_LM_SR_lytic$Population == vir & slopes_LM_SR_lytic$Sample_Type == st,]$LM_SR_Slope), xlab = paste0(tim, vir, st))
      boxplot(as.numeric(slopes_LM_AR_lytic[slopes_LM_AR_lytic$Population == vir & slopes_LM_AR_lytic$Sample_Type == st,]$LM_AVG_Slope), xlab = paste0(tim, vir, st))
    }
  }
}

par(mfrow = c(4,3))
for(tim in unique(slopes_LM_AR_lytic$Time_Range)){
  for(st in c("VP", "VPC")){
    for (vir in unique(slopes_LM_AR_lytic$Population)) {
      boxplot(as.numeric(slopes_LM_AP_lytic[slopes_LM_AP_lytic$Population == vir & slopes_LM_AP_lytic$Sample_Type == st , ]$LM_AP_SE), xlab = paste0(tim, vir, st) )
      boxplot(as.numeric(slopes_LM_SR_lytic[slopes_LM_SR_lytic$Population == vir & slopes_LM_SR_lytic$Sample_Type == st,]$LM_SR_SE), xlab = paste0(tim, vir, st))
      boxplot(as.numeric(slopes_LM_AR_lytic[slopes_LM_AR_lytic$Population == vir & slopes_LM_AR_lytic$Sample_Type == st,]$LM_AVG_SE), xlab = paste0(tim, vir, st))
    }
  }
}


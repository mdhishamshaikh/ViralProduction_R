source("0_vp_source.R")
source("vp_functions.R")
source("vipcal_functions.R")
library("tidyverse")


####2.1 VIPCAL - with separate replicate.
#We will average all the VP and VPC replicates, followed by suvtracting them to get Diff


vp_sr<- df_sr_tp(data)


VPCL_SR<- vipcal_sr(vp_sr)  %>% arrange( 'Location',
                                         'Expt_No',
                                         'Depth',
                                         factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                                         factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
                                         factor(Time_Range, levels = c('T0_T3', 'T0_T6', 'T0_T9', 'T0_17', 'T0_T24')))

#We shoudl compare this dataframe to LM_SR dataframe to see if how every replicate was calculated. For that I iwll have to adjust the NAs
#here to 0 and the negatives in LM_SR to zeros too.
#calc using vipcal per replicate for VP and VPC.
#I iwll average these replicates and then subtract them to get Diff




VPCL_SR_avg<-VPCL_SR%>% group_by(Location, Expt_No, Depth, Time_Range, Population, Sample_Type)%>%
  summarise(n = n(), VP_Mean = mean(VP, na.rm =T), 
            VP_SE = plotrix::std.error(VP, na.rm = T))%>% 
  calc_diff_vpcl_SR() %>% arrange( 'Location',
                                   'Expt_No',
                                   'Depth',
                                   factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                                   factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
                                   factor(Time_Range, levels = c('T0_T3', 'T0_T6', 'T0_T9', 'T0_17', 'T0_T24')))

calc_diff_vpcl_SR<- function(df){
  #converting all the NAs and NaNs to zeros
  df[is.na(df)] <-0
  
  
  VP<- df[df$Sample_Type == "VP",]
  VPC<- df[df$Sample_Type == "VPC",]
  Diff<- as.data.frame(cbind(VP[,1:6], VPC[,8] - VP[,8], VPC[,9] + VPC[,9]))
  Diff[,6]<- "Diff"
  colnames(Diff)[7]<- 'VP_Mean'
  colnames(Diff)[8]<- 'VP_SE'
  
  df<- full_join(df%>%select(-n), Diff, by = NULL)
  return(df)
}

####2.2 AVERAGE replicates VP and VPC. Followed by subtraction for Diff (SE)

vp_avg<- df_avg_tp(data)
vp_avg<- vp_avg[vp_avg$Sample_Type != 'Diff',]

VPCL_AVG<- vipcal_avg_se(vp_avg)

VPCL_AVG<- VPCL_AVG%>% calc_diff_vpcl_AR()  %>% arrange( 'Location',
                                                         'Expt_No',
                                                         'Depth',
                                                         factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                                                         factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
                                                         factor(Time_Range, levels = c('T0_T3', 'T0_T6', 'T0_T9', 'T0_17', 'T0_T24')))

VPCL_AVG[VPCL_AVG$VPCL_AVG == 0, ]$VPCL_AVG<- NA

calc_diff_vpcl_AR<- function(df){
  df[is.na(df)]<- 0
  
  VP<- df[df$Sample_Type == "VP",]
  VPC<- df[df$Sample_Type == "VPC",]
  Diff<- as.data.frame(cbind(VP[,1:6], VPC[,7] -VP[,7], VPC[,8] + VPC[,8]))
  Diff[,6]<- "Diff"
  colnames(Diff)[7]<- 'VP_Mean'
  colnames(Diff)[8]<- 'VP_SE'
  
  df<- full_join(df, Diff, by = NULL)
  return(df)
}

####2.3 Calculate the differenc curve first and the calc avg vipcal (SE)
vp_avg_diff<- df_avg_tp(data)
VPCL_AVG_Diff<- vipcal_avg_se(vp_avg_diff)  %>% arrange( 'Location',
                                                         'Expt_No',
                                                         'Depth',
                                                         factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                                                         factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
                                                         factor(Time_Range, levels = c('T0_T3', 'T0_T6', 'T0_T9', 'T0_17', 'T0_T24')))
#diff has positive values


####2.4 aVERAGED REPLICATE . NO STANDARD ERROR. 

###2.4.1 FIRST VP, VPC AND THEN SUBTREACT FOR DIFF CURVE
vp_avg_NO_SE<- df_avg_tp(data)
vp_avg_NO_SE<- vp_avg_NO_SE[vp_avg_NO_SE$Sample_Type != 'Diff',]

VPCL_AVG_NO_SE<- vipcal_avg(vp_avg_NO_SE) %>% arrange( 'Location',
                                                       'Expt_No',
                                                       'Depth',
                                                       factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                                                       factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
                                                       factor(Time_Range, levels = c('T0_T3', 'T0_T6', 'T0_T9', 'T0_17', 'T0_T24')))
sum(is.na(VPCL_AVG_NO_SE$VPCL_AVG))

calc_diff_vpcl_AR_No_SE<- function(df){
  df[is.na(df)]<- 0
  
  VP<- df[df$Sample_Type == "VP",]
  VPC<- df[df$Sample_Type == "VPC",]
  Diff<- as.data.frame(cbind(VP[,1:6], VPC[,7] -VP[,7]))
  Diff[,6]<- "Diff"
  colnames(Diff)[7]<- 'VP_Mean'
  
  df<- full_join(df, Diff, by = NULL)
  return(df)
}

VPCL_AVG_NO_SE<- VPCL_AVG_NO_SE%>%calc_diff_vpcl_AR_No_SE()  %>% arrange( 'Location',
                                                                          'Expt_No',
                                                                          'Depth',
                                                                          factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                                                                          factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
                                                                          factor(Time_Range, levels = c('T0_T3', 'T0_T6', 'T0_T9', 'T0_17', 'T0_T24')))

###2.4.2 DIFF VURVE AND THEN VIPCAL no se
vp_avg_NO_SE<- df_avg_tp(data)
VPCL_AVG_NO_SE_Diff<- vipcal_avg(vp_avg_NO_SE)  %>% arrange( 'Location',
                                                             'Expt_No',
                                                             'Depth',
                                                             factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                                                             factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
                                                             factor(Time_Range, levels = c('T0_T3', 'T0_T6', 'T0_T9', 'T0_17', 'T0_T24')))
sum(is.na(VPCL_AVG_NO_SE_Diff$VPCL_AVG))

VPCL_AVG_NO_SE_Diff[is.na(VPCL_AVG_NO_SE_Diff)]<- 0


#####2.5 VIPCAL_LMER
VIPCAL_SR_DIFF_LMER_SE<- vipcal_sr_diff_SE(vp_sr)  %>% arrange( 'Location',
                                                                'Expt_No',
                                                                'Depth',
                                                                factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                                                                factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
                                                                factor(Time_Range, levels = c('T0_T3', 'T0_T6', 'T0_T9', 'T0_17', 'T0_T24')))
sum(is.na(VIPCAL_SR_DIFF_LMER_SE$VPCL_LMER_Diff_Slope))
VIPCAL_SR_DIFF_LMER_SE[is.na(VIPCAL_SR_DIFF_LMER_SE)]<- 0

VIPCAL_SR_DIFF_LMER_No_SE<- vipcal_sr_diff_no_SE(vp_sr)  %>% arrange( 'Location',
                                                                      'Expt_No',
                                                                      'Depth',
                                                                      factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                                                                      factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
                                                                      factor(Time_Range, levels = c('T0_T3', 'T0_T6', 'T0_T9', 'T0_17', 'T0_T24')))
sum(is.na(VIPCAL_SR_DIFF_LMER_No_SE$VPCL_LMER_Diff_Slope))
VIPCAL_SR_DIFF_LMER_No_SE[is.na(VIPCAL_SR_DIFF_LMER_No_SE)]<- 0


#different no of peaks and valleys
"NJ2020 1 1 T0_T24 c_V1 Diff"
time<- 'T0_T24'
viruses<- 'c_V3'
location<- 'NJ2020'
expt_no<- 1
depth<- 1
prod<- 'Diff'

df_sr<- vp_sr
plot(x =c(10e+10, -1305197.6, -1048986.9, 2338073.3, 1132680.2, -691991.6, -1439607.6, -10e+10))

peaks(c(10e+10, -1305197.6, -1048986.9, 2338073.3, 1132680.2, -691991.6, -1439607.6, -10e+10))
valleys(c(10e+10, -1305197.6, -1048986.9, 2338073.3, 1132680.2, -691991.6, -1439607.6, -10e+10)) -1

f<- VIPCAL_SR_DIFF_LMER_No_SE
colnames(f)[colnames(f) =='VPCL_LMER_Diff_Slope' ] <- 'VPCL_LMER_no_SE'
f<- pivot_longer(f, cols = 7, names_to = 'VP_Type', values_to = 'VP')

sum(is.na(f))

f[is.na(f$VP),]$VP<- 0

sum(is.na(f$VP))
f$VP


g<- VIPCAL_SR_DIFF_LMER_SE
colnames(g)[colnames(g) =='VPCL_LMER_Diff_Slope' ] <- 'VPCL_LMER_SE'
g<- pivot_longer(g, cols = 7, names_to = 'VP_Type', values_to = 'VP')
sum(is.na(g))

g[is.na(g$VP),]$VP<- 0

sum(is.na(g$VP))
g$VP

h<- full_join(f,g)

#ADD STD ERROR TO THE VIPCAL CALCULATIONS

ggplot(h, aes(x= VP_Type, y = VP, col = Sample_Type))+  geom_point()
 sum(is.na(f$VP))
 sum(is.na(g$VP))
 
 ggplot()+
   geom_point(aes(x = g$VP, y = f$VP, col = g$Sample_Type, shape = g$Population))+
   geom_abline(intercept = 0)+
   scale_x_continuous(limits = c(-100, 10e+4))+
   scale_y_continuous(limits = c(-100, 10e+4))
 
 ggplot()+
   geom_point(aes(x = g[g$Sample_Type == 'Diff',]$VP, y = f[f$Sample_Type == 'Diff',]$VP))+
   geom_abline(intercept = 0)
 
 summary(lm(g$VP ~ f$VP))
 plot(lm(g$VP ~ f$VP))
 
##Setting Up the environment
##AIM: The purpose of this script is to go through all the possible linear 
#regression methods, and the select the few to compare.

source("0_vp_source.R")
source("vp_functions.R")

library("tidyverse")

data<- read.csv("NJ1.csv")

####1.0 Linear Regression####

####1.1 Slope through all the points####
data_sr<- df_sr_tp(data)
slopes_LM_AP<- slope_all_points(data_sr) %>%
  calc_diff_lm_AP() 

slopes_LM_AP<- slopes_LM_AP%>%
  arrange('Location',
          'Expt_No',
          'Depth',
          factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
          factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
          factor(Time_Range, levels = c('T0_T3', 'T0_T6', 'T0_T9', 'T0_17', 'T0_T24'))) 


####1.2 Slope through every replicate####
#we use data_sr for this
slopes_LM_SR<- slope_lm_sr(data_sr)%>%
  arrange('Location',
          'Expt_No',
          'Depth',
          factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
          factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
          factor(Time_Range, levels = c('T0_T3', 'T0_T6', 'T0_T9', 'T0_17', 'T0_T24')))#3 slopes per sample type
#i need to average the replicates here, and then i can calculate Diff

slopes_LM_SR_avg<- slopes_LM_SR %>%
  group_by(Location, Expt_No, Depth, Sample_Type, Time_Range, Population ) %>%
  summarise(LM_SR_Slope_Mean=mean(LM_SR_Slope), LM_SR_Slope_SE=plotrix::std.error((LM_SR_Slope))) %>%
  arrange('Location',
          'Expt_No',
          'Depth',
          factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
          factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
          factor(Time_Range, levels = c('T0_T3', 'T0_T6', 'T0_T9', 'T0_17', 'T0_T24')))



#plot(slopes_LM_AP$LM_AP_Slope ~ slopes_LM_SR_avg$LM_SR_Slope_Mean)
#summary(lm(slopes_LM_AP$LM_AP_Slope ~ slopes_LM_SR_avg$LM_SR_Slope_Mean))



####1.3 slope through averaged replicates####

data_avg<- df_sr_tp(data)
data_avg<- data_avg[data_avg$Microbe == 'Viruses',]
data_avg<- data_avg %>%
  group_by(Location, Expt_No, Depth, Sample_Type, Timepoint, Time_Range, Population ) %>% 
  summarise(n =n(), mean=mean(count), sd=sd(count))

data_avg$Microbe <- 'Viruses'

slopes_LM_AR<- slope_lm_avg(data_avg)

slopes_LM_AR<- calc_diff_lm_AR(slopes_LM_AR) %>%
  arrange('Location',
          'Expt_No',
          'Depth',
          factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
          factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
          factor(Time_Range, levels = c('T0_T3', 'T0_T6', 'T0_T9', 'T0_17', 'T0_T24')))


#plot(slopes_LM_AP$LM_AP_Slope ~ slopes_LM_AR$LM_AVG_Slope)
#summary(lm(slopes_LM_AP$LM_AP_Slope ~ slopes_LM_AR$LM_AVG_Slope))


#plot(slopes_LM_SR_avg$LM_SR_Slope_Mean ~ slopes_LM_AR$LM_AVG_Slope)
#summary(lm(slopes_LM_SR_avg$LM_SR_Slope_Mean ~ slopes_LM_AR$LM_AVG_Slope))

ggplot()+
  geom_point(aes( x= slopes_LM_SR_avg$LM_SR_Slope_Mean , y = slopes_LM_AR$LM_AVG_Slope, col = slopes_LM_AP$Sample_Type, shape = slopes_LM_AP$Population ))


####1.4 Calculate difference Curve first ####

#We could use the df_avg_tp function for this

df_avg_diff<- df_avg_tp(data)%>%
  group_by('Location', 'Expt_No', 'Depth', 'Sample_Type', 'Population', 'Time_Range')%>%
  arrange('Location',
          'Expt_No',
          'Depth',
          factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
          factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
          factor(Time_Range, levels = c('T0_T3', 'T0_T6', 'T0_T9', 'T0_17', 'T0_T24')),
          as.numeric(Timepoint))

slopes_LM_Diff <- slope_lm_avg(df_avg_diff)

slopes_LM_Diff<- slopes_LM_Diff%>%
  arrange('Location',
          'Expt_No',
          'Depth',
          factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
          factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
          factor(Time_Range, levels = c('T0_T3', 'T0_T6', 'T0_T9', 'T0_17', 'T0_T24')))

plot(slopes_LM_Diff$LM_AVG_Slope ~ slopes_LM_AR$LM_AVG_Slope)

summary(lm(slopes_LM_Diff$LM_AVG_Slope ~ slopes_LM_AR$LM_AVG_Slope))


####1.5 Using LMER to calculate Diff Curve####
#I will use all teh points for this. So df_sr

lmer_model(data_sr)

slope_lm_ar_diff_lmer<- function(df_sr){ #takes SR dataframe as an input
  
  lm_vp<- list()
  slope_lm_sr_df<- data.frame()
  for (location in unique(df_sr$Location)){
    for (expt_no in unique(df_sr$Expt_No)){
      for (depth in unique(df_sr$Depth)){
        for (time in unique(df_sr$Time_Range)){
          for (viruses in unique(df_sr[df_sr$Microbe == "Viruses",]$Population)){
                
                
                df2<- df_sr[df_sr$Location == location,]
                df2<- df2[df2$Expt_No == expt_no,]
                df2<- df2[df2$Depth == depth,]
                df2<- df2[df2$Time_Range == time,]
                df2<- df2[df2$Population == viruses,]
                
                try(df3<- lmer_model(df2))#the output is a dataframe for VP, VPC and Diff together
                
                if (exists("df3")){
                for (prod in unique(df3$Sample_Type)){
                  
                  lm<- summary(lm(data = df3[df3$Sample_Type == prod,], Mean ~ as.numeric(Timepoint)))
                  print(lm)
                  slope<- c(location, expt_no, depth, time, viruses, prod, lm$coefficients[c(2,4)], lm$r.squared)
                  print(lm$coefficients[,3])
                  print(lm$coefficients[[2]])
                  lm_vp[[length(lm_vp)+1]] <- slope 
                  
                }
                } else{
                  for (prod in c("Diff", "VP", "VPC")){
                  slope<- c(location, expt_no, depth, time, viruses, prod, NA, NA, NA)
                  lm_vp[[length(lm_vp)+1]] <- slope 
                  }
                }
                
                
                 
          }     
        }
      }
    }
  }
  
  slope_lm_sr_df<- data.frame(t(sapply(lm_vp, c)))
  colnames(slope_lm_sr_df)<- c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'VP_Slope', 'VP_SE', 'VP_R_squared')
  slope_lm_sr_df[, c('VP_Slope', 'VP_SE', 'VP_R_squared')]<- lapply(slope_lm_sr_df[, c('VP_Slope', 'VP_SE', 'VP_squared')], as.numeric)
  return(slope_lm_sr_df)
  rm(df3)
}





slopes_LM_SR_Diff<- slope_lm_sr_diff(data_sr)

slopes_LM_SR_Diff<- slopes_LM_SR_Diff %>%
  arrange('Location',
          'Expt_No',
          'Depth',
          factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
          factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
          factor(Time_Range, levels = c('T0_T3', 'T0_T6', 'T0_T9', 'T0_17', 'T0_T24')))
  

plot(slopes_LM_SR_Diff$LM_Diff_Slope ~ slopes_LM_AR$LM_AVG_Slope)

summary(lm(slopes_LM_SR_Diff$LM_Diff_Slope~ slopes_LM_AR$LM_AVG_Slope))


ggplot()+
  geom_point(aes( x= slopes_LM_SR_Diff$LM_Diff_Slope, y = slopes_LM_AR$LM_AVG_Slope, col = slopes_LM_AP$Time_Range, shape = slopes_LM_AP$Population ))

#just the T0_T3 are off, which makes sense. But i need to code NAs in there instead of random values. 


a<- pivot_longer(slopes_LM_AP[,1:7], LM_AP_Slope, values_to = 'Slope', names_to = 'Slope_Type')

b<- pivot_longer(slopes_LM_SR_avg[,1:7], LM_SR_Slope_Mean, values_to = 'Slope', names_to = 'Slope_Type')

c<- pivot_longer(slopes_LM_AR[,1:7], LM_AVG_Slope, values_to = 'Slope', names_to = 'Slope_Type')

d<- pivot_longer(slopes_LM_Diff[,1:7], LM_AVG_Slope, values_to = 'Slope', names_to = 'Slope_Type')
d$Slope_Type <- 'LM_AVG_Diff'

e<- pivot_longer(slopes_LM_SR_Diff[,1:7], LM_Diff_Slope, values_to = 'Slope', names_to = 'Slope_Type')


lm_slopes<- full_join(a,b) %>%
  full_join(c)%>%
  full_join(d)%>%
  full_join(e)


ggplot( data = lm_slopes, aes(x = Slope_Type, y = Slope, col = Sample_Type, shape= Population))+
  geom_point(size = 2.5)+
  theme_bw()+
  xlab("Linear Regression Methods")+
  ylab("Slopes")+
  scale_color_lancet(alpha = 0.75)





lmer_model<- function(df, value = count){
  
  lmer_data<- data.frame()
  model_plots<- list()
  n<-0
  for (rep in c(1,2,3)){
    df$Replicate[df$Replicate == rep & df$Sample_Type == 'VPC'] <- rep+3
  }
  
  for (virus in unique(df$Population)){
    
    model<- lme4::lmer(data = df, count ~ Sample_Type*as.factor(Timepoint) + (1+ Sample_Type | Replicate))
    plot<- model_plot(model, df = df)
    emmeans<- emmeans::emmeans(model, ~ Sample_Type|as.factor(Timepoint))
    contrast<- pairs(emmeans) 
    df1<- data.frame(rep('c_Viruses', 6),
                     rep("Diff", 6),
                     summary(contrast)$Timepoint, 
                     -(summary(contrast)$estimate),
                     summary(contrast)$SE
    )
    colnames(df1)<- c("Population", "Sample_Type", "Timepoint", "Mean", "SE")
    st<- summary(emmeans)[1]
    tp<- summary(emmeans)[2]
    mean<- summary(emmeans)[3]
    se<- summary(emmeans)[4]
    pop<- rep('c_Viruses', 6)
    df2<- data.frame(pop,st,tp,mean,se)
    colnames(df2)<- c("Population", "Sample_Type", "Timepoint", "Mean", "SE")
    
    
    lmer_data<- rbind(lmer_data, df2, df1)
    # DF<- rbind(df2, df1) %>%
    #  arrange(Timepoint)
    #n<- n +1
    #lmer_data[[n]]<- DF
    #model_plots[[n]]<- plot
    
  }
  return(lmer_data)
  return(model_plots)
  
  
}

summary(lm(lmer_data[lmer_data$Sample_Type == "VP",]$Mean  ~ lmer_data[lmer_data$Sample_Type == "VP",]$Timepoint))

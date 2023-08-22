#####1.0 Installing Packages####
source("C:/Users/hisham.shaikh/OneDrive - UGent/Projects/FCM_R/ViralProduction_R/sourcesourcebaby.R")
library(tidyverse)
library(ggstatsplot)
library(ggsci)
#library(moments)


#We have the dataframes
#For LM we use no rep, per rep, and averaged. Therefore, find a way to 

simu_df_AR_Diff_function<- function(df = simu_df){
  simu_df_AR<- df %>% group_by(Expt_No, Sample_Type, Timepoint )%>%
    summarise(n = n(), mean = mean(Count), se = plotrix::std.error(Count))
  df5<- simu_df_AR
  
  VP<- df5[df5$Sample_Type == "VP",]
  VPC<- df5[df5$Sample_Type == "VPC",]
  Diff<- as.data.frame(cbind(VP[,c('Expt_No', 'Sample_Type', 'Timepoint')], VPC[,'mean'] -VP[,'mean'], VPC[,'se'] + VP[,'se']))
  Diff[,'Sample_Type']<- "Diff"
  
  simu_df_AR_Diff<- full_join(df5, Diff, by = NULL)
  rm(df)
  return(simu_df_AR_Diff)
}


#Calculating Difference Curve

simu_calc_diff_lm_AR<- function(df, variables){
  
  VP<- df[df$Sample_Type == "VP",]
  VPC<- df[df$Sample_Type == "VPC",]
  Diff<- as.data.frame(cbind(VP[,variables], VPC[,'VP'] -VP[,'VP'], VPC[,'VP_SE'] + VP[,'VP_SE']))
  Diff[,'Sample_Type']<- "Diff"
  colnames(Diff)[3]<- 'VP'
  colnames(Diff)[4]<- 'VP_SE'
  Diff$VP_R_Squared <- NA
  
  df<- full_join(df, Diff, by = NULL)
  return(df)
}


#lmer model fucntion for the simulation

simu_lmer_model<- function(df){
  
  lmer_data<- data.frame()
  #model_plots<- list()
  n<-0
  
  
  
  model<- lme4::lmer(data = df, Count ~ Sample_Type*as.factor(Timepoint) + (1  | Replicate))
  #plot<- model_plot(model, df = df)
  emmeans<- emmeans::emmeans(model, ~ Sample_Type|as.factor(Timepoint))
  contrast<- pairs(emmeans) 
  dataf1<- data.frame(rep("Diff", length(unique(df$Timepoint))),
                      summary(contrast)$Timepoint, 
                      -(summary(contrast)$estimate),
                      summary(contrast)$SE
  )
  colnames(dataf1)<- c("Sample_Type", "Timepoint", "Mean", "SE")
  st<- summary(emmeans)[1]
  tp<- summary(emmeans)[2]
  mean<- summary(emmeans)[3]
  se<- summary(emmeans)[4]
  dataf2<- data.frame(st,tp,mean,se)
  colnames(dataf2)<- c("Sample_Type", "Timepoint", "Mean", "SE")
  
  
  lmer_data<- rbind(lmer_data, dataf2, dataf1)
  # DF<- rbind(dataf2, dataf1) %>%
  #  arrange(Timepoint)
  #n<- n +1
  #lmer_data[[n]]<- DF
  #model_plots[[n]]<- plot
  
  
  return(lmer_data)
  return(model_plots)
  
  
}


#### Simulation functions####

#Functions for the 12 different methods for calculation of viral production data
{
  
  ####2.1 LM_AP####
  
  simu_lm_ap <- function(df = simu_df, variables = var){
    
    lm_vp<- list()
    slope<- data.frame()
    
    
    for (expt in length_df) {
      
      print(expt)
      
      for(type in unique(df$Sample_Type)){
        
        df1<- df %>% filter(Expt_No == expt,
                            Sample_Type == type)
        
        try(lm<- summary(lm(data = df1, Count ~ as.numeric(Timepoint))))
        slope1<- c(expt, type, lm$coefficients[c(2,4)], lm$r.squared)
        lm_vp[[length(lm_vp)+1]] <- slope1
        
      }
      
      print('done')
      rm(df1, expt, type, slope1, lm)
      
    }
    
    {
      
      slope<- data.frame(t(sapply(lm_vp, c)))
      colnames(slope)<- c('Expt_No', 'Sample_Type', 'VP', 'VP_SE', 'VP_R_Squared')
      slope$VP_Type<- 'LM_AP'
      slope<- slope%>% mutate(across(c('Expt_No', 'VP', 'VP_SE', 'VP_R_Squared'), as.numeric))
      
    }
    
    #calculate Lysogeny from here.
    
    LM_AP_output<- slope %>%
      simu_calc_diff_lm_AR(variables) %>%
      mutate(VP_Type = 'LM_AP')
    
    rm(slope, lm_vp)
    
    return(LM_AP_output)
    
  }
  
  
  ####2.2 LM_SR_AVG####
  
  simu_lm_sr_avg<- function(df = simu_df, variables = var){
    lm_vp<- list()
    slope<- data.frame()
    
    
    for (expt in length_df) {
      
      print(expt)
      
      for (rep in unique(df$Replicate)) {
        
        df1 <- df %>% filter(Expt_No == expt,
                             Replicate == rep)
        type<- unique(df1$Sample_Type)
        print(paste(expt, type, rep))
        
        try(lm<- summary(lm(data = df1, Count ~ as.numeric(Timepoint))))
        slope1<- c(expt, type, rep, lm$coefficients[c(2,4)], lm$r.squared)
        lm_vp[[length(lm_vp)+1]] <- slope1
        
      }
      
      print('done')
      rm(df1, expt, type, rep, slope1, lm)
      
    }
    
    {
      
      slope<- data.frame(t(sapply(lm_vp, c)))
      colnames(slope)<- c('Expt_No', 'Sample_Type', 'Replicate', 'VP', 'VP_SE', 'VP_R_Squared')
      slope$VP_Type<- 'LM_SR_AVG'
      slope<- slope%>% mutate(across(c('Expt_No', 'VP', 'VP_SE', 'VP_R_Squared'), as.numeric))
      
    }
    
    #calculate Lysogeny from here.
    
    LM_SR_AVG_output<- slope %>%
      group_by(Expt_No, Sample_Type) %>%
      summarise(n = n(), VP2 = mean(VP), VP_SE = plotrix::std.error(VP)) %>%
      rename(VP = VP2) %>%
      simu_calc_diff_lm_AR(variables) %>%
      mutate(VP_Type = 'LM_SR_AVG')
    
    rm(slope, lm_vp)
    
    return(LM_SR_AVG_output)
    
  }
  
  
  ####2.3 LM_AR ####
  
  simu_lm_ar <- function(df = simu_df, variables = var2){
    simu_df_AR<- df %>% group_by(Expt_No, Sample_Type, Timepoint )%>%
      summarise(n = n(), mean = mean(Count), se = plotrix::std.error(Count))
    
    lm_vp<- list()
    slope<- data.frame()
    
    for (expt in length_df){
      
      print(expt)
      
      for (type in unique(simu_df_AR$Sample_Type)){
        df1<- simu_df_AR %>% filter(Expt_No == expt,
                                    Sample_Type == type)
        
        try(lm<- summary(lm(data = df1, mean ~ as.numeric(Timepoint))))
        slope1<- c(expt, type, lm$coefficients[c(2,4)], lm$r.squared)
        lm_vp[[length(lm_vp)+1]] <- slope1
        
      }
      print('done')
      rm(df1, expt, type, slope1, lm)
    }
    
    {
      slope<- data.frame(t(sapply(lm_vp, c)))
      colnames(slope)<- c('Expt_No', 'Sample_Type', 'VP', 'VP_SE', 'VP_R_Squared')
      slope$VP_Type<- 'LM_AR'
      slope<- slope%>% mutate(across(c('Expt_No', 'VP', 'VP_SE', 'VP_R_Squared'), as.numeric))
      
    }
    
    LM_AR_output<- slope %>%
      group_by(Expt_No, Sample_Type) %>%
      simu_calc_diff_lm_AR(var) %>%
      mutate(VP_Type = 'LM_AR')
    
    
    rm( slope, lm_vp)
    
    return(LM_AR_output)
    
  }
  
  
  #####2.4 LM_AR_Diff ####
  
  simu_lm_ar_diff <- function(df = simu_df){
    
    simu_df_AR_Diff<- simu_df_AR_Diff_function(df)
    
    lm_vp<- list()
    slope<- data.frame()
    
    for (expt in length_df){
      
      print(expt)
      
      for (type in unique(simu_df_AR_Diff$Sample_Type)){
        df1<- simu_df_AR_Diff %>% filter(Expt_No == expt,
                                         Sample_Type == type)
        
        try(lm<- summary(lm(data = df1, mean ~ as.numeric(Timepoint))))
        slope1<- c(expt, type, lm$coefficients[c(2,4)], lm$r.squared)
        lm_vp[[length(lm_vp)+1]] <- slope1
        
      }
      print('done')
      rm(df1, expt, type, slope1, lm)
    }
    
    slope<- data.frame(t(sapply(lm_vp, c)))
    colnames(slope)<- c('Expt_No', 'Sample_Type', 'VP', 'VP_SE', 'VP_R_Squared')
    slope$VP_Type<- 'LM_AR_Diff'
    slope<- slope%>% mutate(across(c('Expt_No', 'VP', 'VP_SE', 'VP_R_Squared'), as.numeric))
    
    LM_AR_Diff_output<- slope %>%
      group_by(Expt_No, Sample_Type) 
    
    rm( slope, lm_vp)
    
    return(LM_AR_Diff_output)
    
  }
  
  
  ####2.5 LM_aR_Diff_LMER ####
  
  simu_lm_ar_diff_lmer <- function(df = simu_df){
    
    
    lm_vp<- list()
    slope<- data.frame()
    
    for (expt in length_df){
      
      print(expt)
      
      
      df1<- df %>% filter(Expt_No == expt)
      
      try(df2<- simu_lmer_model(df1))
      
      
      try({
        if (exists("df2")){
          for (type in unique(df2$Sample_Type)){
            
            lm<- summary(lm(data = df2[df2$Sample_Type == type,], Mean ~ as.numeric(Timepoint)))
            #print(lm)
            slope<- c(expt, type, lm$coefficients[c(2,4)], lm$r.squared)
            #print(lm$coefficients[,3])
            #print(lm$coefficients[[2]])
            lm_vp[[length(lm_vp)+1]] <- slope 
            
          }
        } else{
          for (type in c("Diff", "VP", "VPC")){
            slope<- c(expt, type, NA, NA, NA)
            lm_vp[[length(lm_vp)+1]] <- slope 
          }
        }
      })
      
      
      
      print('done')
      rm(df1, expt,type, lm)
    }
    {
      slope<- data.frame(t(sapply(lm_vp, c)))
      colnames(slope)<- c('Expt_No', 'Sample_Type', 'VP', 'VP_SE', 'VP_R_Squared')
      slope<- slope%>% mutate(across(c('Expt_No', 'VP', 'VP_SE', 'VP_R_Squared'), as.numeric))
    }
    LM_AR_Diff_LMER_output<- slope %>%
      group_by(Expt_No, Sample_Type)%>%
      mutate(VP_Type = 'LM_AR_Diff_LMER')
    
    rm( slope, lm_vp)
    
    return(LM_AR_Diff_LMER_output)
    
  }
  
  
  ####2.6 VPCL_SR_AVG
  
  simu_vpcl_sr_avg<-function(df = simu_df, variables = var){
    vipcal_vp<- list()
    vpcl<- data.frame()
    
    
    for (expt in length_df) {
      
      print(expt)
      
      for (rep in unique(df$Replicate)) {
        
        df1 <- df %>% filter(Expt_No == expt,
                             Replicate == rep)
        type<- unique(df1$Sample_Type)
        print(paste(expt, type, rep))
        
        
        try({
          p<- peaks(c(+10e+10, df1$Count, -10e+10))-1
          v<- valleys(c(+10e+10, df1$Count, -10e+10))-1
          
          if(identical(length(p), length(v))){
            print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
          }else{
            print("Number peaks and valleys are not the same. This will lead to erroneous viral production calculations")
          }
          
          if(length(p)==0){
            vp<- 0
            
          } else if (length(p)==1) {
            vp<- (df1$Count[p[1]] - df1$Count[v[1]])/(df1$Timepoint[p[1]] - df1$Timepoint[v[1]])
            
          } else if (length(p)==2) {
            vp<- ((df1$Count[p[1]] - df1$Count[v[1]])/(df1$Timepoint[p[1]] - df1$Timepoint[v[1]]) + 
                    (df1$Count[p[2]] - df1$Count[v[2]])/(df1$Timepoint[p[2]] - df1$Timepoint[v[2]]))/2
            
          } else if (length(p)==3) {
            vp<-  ((df1$Count[p[1]] - df1$Count[v[1]])/(df1$Timepoint[p[1]] - df1$Timepoint[v[1]]) + 
                     (df1$Count[p[2]] - df1$Count[v[2]])/(df1$Timepoint[p[2]] - df1$Timepoint[v[2]]) +
                     (df1$Count[p[3]] - df1$Count[v[3]])/(df1$Timepoint[p[3]] - df1$Timepoint[v[3]]))/3
          }
          
          
          vipcal<- c(expt, type, rep, vp)
          
          vipcal_vp[[length(vipcal_vp)+1]] <- vipcal
        })
        
        
      }
      
      print('done')
      rm(df1, expt, type, rep, vp)
      
    }
    
    {
      vpcl<- data.frame(t(sapply(vipcal_vp, c)))
      colnames(vpcl)<- c('Expt_No', 'Sample_Type', 'Replicate', 'VP')
      vpcl[, c('Expt_No', 'VP')]<- lapply(vpcl[, c('Expt_No', 'VP')], as.numeric)
      vpcl$VP_Type<- 'VPCL_SR_AVG'
      
    }
    
    #calculate Lysogeny from here.
    
    VPCL_SR_AVG_output<- vpcl %>%
      group_by(Expt_No, Sample_Type) %>%
      summarise(n = n(), VP2 = mean(VP), VP_SE = plotrix::std.error(VP)) %>%
      rename(VP = VP2) %>%
      simu_calc_diff_lm_AR(variables) %>%
      mutate(VP_Type = 'VPCL_SR_AVG') %>%
      select(-n)
    
    rm(vpcl, vipcal_vp)
    
    return(VPCL_SR_AVG_output)
    
  }
  
  
  ####2.7 VPCL_AR_No_SE####
  
  simu_vpcl_ar_no_se<-function(df = simu_df, variables = var){
    simu_df_AR<- df %>% group_by(Expt_No, Sample_Type, Timepoint )%>%
      summarise(n = n(), mean = mean(Count), se = plotrix::std.error(Count))
    
    
    vipcal_vp<- list()
    vpcl<- data.frame()
    
    
    for (expt in length_df) {
      
      print(expt)
      
      for (type in unique(simu_df_AR$Sample_Type)){
        df1 <- simu_df_AR %>% filter(Expt_No == expt,
                                     Sample_Type == type)
        print(paste(expt, type))
        
        
        try({
          p<- peaks(c(+10e+10, df1$mean, -10e+10))-1
          v<- valleys(c(+10e+10, df1$mean, -10e+10))-1
          
          if(identical(length(p), length(v))){
            print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
          }else{
            print("Number peaks and valleys are not the same. This will lead to erroneous viral production calculations")
          }
          
          if(length(p)==0){
            vp<- 0
            
          } else if (length(p)==1) {
            vp<- (df1$mean[p[1]] - df1$mean[v[1]])/(df1$Timepoint[p[1]] - df1$Timepoint[v[1]])
            
          } else if (length(p)==2) {
            vp<- ((df1$mean[p[1]] - df1$mean[v[1]])/(df1$Timepoint[p[1]] - df1$Timepoint[v[1]]) + 
                    (df1$mean[p[2]] - df1$mean[v[2]])/(df1$Timepoint[p[2]] - df1$Timepoint[v[2]]))/2
            
          } else if (length(p)==3) {
            vp<-  ((df1$mean[p[1]] - df1$mean[v[1]])/(df1$Timepoint[p[1]] - df1$Timepoint[v[1]]) + 
                     (df1$mean[p[2]] - df1$mean[v[2]])/(df1$Timepoint[p[2]] - df1$Timepoint[v[2]]) +
                     (df1$mean[p[3]] - df1$mean[v[3]])/(df1$Timepoint[p[3]] - df1$Timepoint[v[3]]))/3
          }
          
          
          vipcal<- c(expt, type, vp)
          
          vipcal_vp[[length(vipcal_vp)+1]] <- vipcal
        })
        
        
        
        print('done')
        
        
      }
      
      rm(df1, expt, type, vp)
      
    }
    
    {
      vpcl<- data.frame(t(sapply(vipcal_vp, c)))
      colnames(vpcl)<- c('Expt_No', 'Sample_Type', 'VP')
      vpcl[, c('Expt_No', 'VP')]<- lapply(vpcl[, c('Expt_No', 'VP')], as.numeric)
      vpcl$VP_Type<- 'VPCL_AR_No_SE'
      
    }
    
    #calculate Lysogeny from here.
    
    VPCL_AR_No_SE_output<- vpcl %>%
      group_by(Expt_No, Sample_Type) %>%
      summarise(n = n(), VP2 = mean(VP), VP_SE = plotrix::std.error(VP)) %>%
      rename(VP = VP2) %>%
      simu_calc_diff_lm_AR(variables) %>%
      select(-n)%>%
      mutate(VP_Type = 'VPCL_AR_No_SE')
    
    rm(vpcl, vipcal_vp)
    
    return(VPCL_AR_No_SE_output)
    
  }
  
  
  ####2.8 VPCL_AR_SE####
  
  simu_vpcl_ar_se<-function(df = simu_df, variables = var){
    simu_df_AR<- df %>% group_by(Expt_No, Sample_Type, Timepoint )%>%
      summarise(n = n(), mean = mean(Count), se = plotrix::std.error(Count))
    
    vipcal_vp<- list()
    vpcl<- data.frame()
    
    
    for (expt in length_df) {
      
      print(expt)
      
      for (type in unique(simu_df_AR$Sample_Type)){
        df1 <- simu_df_AR %>% filter(Expt_No == expt,
                                     Sample_Type == type)
        print(paste(expt, type))
        
        
        try({
          p<- peaks_se(c(+10e+10, df1$mean, -10e+10),
                       c(0, df1$se, 0))-1
          v<- valleys_se(c(+10e+10, df1$mean, -10e+10),
                         c(0, df1$se, 0))-1
          
          if(identical(length(p), length(v))){
            print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
          }else{
            print("Number peaks and valleys are not the same. This will lead to erroneous viral production calculations")
          }
          
          if(length(p)==0){
            vp<- 0
            
          } else if (length(p)==1) {
            vp<- (df1$mean[p[1]] - df1$mean[v[1]])/(df1$Timepoint[p[1]] - df1$Timepoint[v[1]])
            
          } else if (length(p)==2) {
            vp<- ((df1$mean[p[1]] - df1$mean[v[1]])/(df1$Timepoint[p[1]] - df1$Timepoint[v[1]]) + 
                    (df1$mean[p[2]] - df1$mean[v[2]])/(df1$Timepoint[p[2]] - df1$Timepoint[v[2]]))/2
            
          } else if (length(p)==3) {
            vp<-  ((df1$mean[p[1]] - df1$mean[v[1]])/(df1$Timepoint[p[1]] - df1$Timepoint[v[1]]) + 
                     (df1$mean[p[2]] - df1$mean[v[2]])/(df1$Timepoint[p[2]] - df1$Timepoint[v[2]]) +
                     (df1$mean[p[3]] - df1$mean[v[3]])/(df1$Timepoint[p[3]] - df1$Timepoint[v[3]]))/3
          }
          
          if(length(p)==0){
            se<- 0
            
          } else if (length(p)==1) {
            se<- (df1$se[p[1]] + df1$se[v[1]])/(df1$Timepoint[p[1]] - df1$Timepoint[v[1]])
            
          } else if (length(p)==2) {
            se<- ((df1$se[p[1]] + df1$se[v[1]])/(df1$Timepoint[p[1]] - df1$Timepoint[v[1]]) + 
                    (df1$se[p[2]] + df1$se[v[2]])/(df1$Timepoint[p[2]] - df1$Timepoint[v[2]]))/2
            
          } else if (length(p)==3) {
            se<-  ((df1$se[p[1]] + df1$se[v[1]])/(df1$Timepoint[p[1]] - df1$Timepoint[v[1]]) + 
                     (df1$se[p[2]] + df1$se[v[2]])/(df1$Timepoint[p[2]] - df1$Timepoint[v[2]]) +
                     (df1$se[p[3]] + df1$se[v[3]])/(df1$Timepoint[p[3]] - df1$Timepoint[v[3]]))/3
          }
          
          
          vipcal<- c(expt, type, vp, se)
          
          vipcal_vp[[length(vipcal_vp)+1]] <- vipcal
        })
        
        
        
        print('done')
        
        
      }
      
      rm(df1, expt, type, vp, se)
      
    }
    
    {
      vpcl<- data.frame(t(sapply(vipcal_vp, c)))
      colnames(vpcl)<- c('Expt_No', 'Sample_Type', 'VP', 'VP_SE')
      vpcl[, c('Expt_No', 'VP', 'VP_SE')]<- lapply(vpcl[, c('Expt_No', 'VP', 'VP_SE')], as.numeric)
      vpcl$VP_Type<- 'VPCL_AR_SE'
      
    }
    
    #calculate Lysogeny from here.
    
    VPCL_AR_SE_ouput<- vpcl %>%
      group_by(Expt_No, Sample_Type) %>%
      summarise(n = n(), VP2 = mean(VP), VP_SE = plotrix::std.error(VP)) %>%
      rename(VP = VP2) %>%
      simu_calc_diff_lm_AR(variables) %>%
      select(-n)%>%
      mutate(VP_Type = 'VPCL_AR_SE')
    
    rm(vpcl, vipcal_vp)
    
    return(VPCL_AR_SE_ouput)
    
  }
  
  
  ####2.9 VPCL_AR_Diff_No_SE####
  
  simu_vpcl_ar_diff_no_se<-function(df = simu_df, variables = var){
    
    simu_df_AR_Diff<- simu_df_AR_Diff_function(df)
    
    vipcal_vp<- list()
    vpcl<- data.frame()
    
    
    for (expt in length_df) {
      
      print(expt)
      
      for (type in unique(simu_df_AR_Diff$Sample_Type)){
        df1 <- simu_df_AR_Diff %>% filter(Expt_No == expt,
                                          Sample_Type == type)
        print(paste(expt, type))
        
        
        try({
          p<- peaks(c(+10e+10, df1$mean, -10e+10))-1
          v<- valleys(c(+10e+10, df1$mean, -10e+10))-1
          
          if(identical(length(p), length(v))){
            print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
          }else{
            print("Number peaks and valleys are not the same. This will lead to erroneous viral production calculations")
          }
          
          if(length(p)==0){
            vp<- 0
            
          } else if (length(p)==1) {
            vp<- (df1$mean[p[1]] - df1$mean[v[1]])/(df1$Timepoint[p[1]] - df1$Timepoint[v[1]])
            
          } else if (length(p)==2) {
            vp<- ((df1$mean[p[1]] - df1$mean[v[1]])/(df1$Timepoint[p[1]] - df1$Timepoint[v[1]]) + 
                    (df1$mean[p[2]] - df1$mean[v[2]])/(df1$Timepoint[p[2]] - df1$Timepoint[v[2]]))/2
            
          } else if (length(p)==3) {
            vp<-  ((df1$mean[p[1]] - df1$mean[v[1]])/(df1$Timepoint[p[1]] - df1$Timepoint[v[1]]) + 
                     (df1$mean[p[2]] - df1$mean[v[2]])/(df1$Timepoint[p[2]] - df1$Timepoint[v[2]]) +
                     (df1$mean[p[3]] - df1$mean[v[3]])/(df1$Timepoint[p[3]] - df1$Timepoint[v[3]]))/3
          }
          
          
          vipcal<- c(expt, type, vp)
          
          vipcal_vp[[length(vipcal_vp)+1]] <- vipcal
        })
        
        
        
        print('done')
        
        
      }
      
      rm(df1, expt, type, vp)
      
    }
    
    {
      vpcl<- data.frame(t(sapply(vipcal_vp, c)))
      colnames(vpcl)<- c('Expt_No', 'Sample_Type', 'VP')
      vpcl[, c('Expt_No', 'VP')]<- lapply(vpcl[, c('Expt_No', 'VP')], as.numeric)
      vpcl$VP_Type<- 'VPCL_AR_Diff_No_SE'
      vpcl$VP_SE<- NA
      vpcl$VP_R_Squared<- NA
      
    }
    
    #calculate Lysogeny from here.
    
    VPCL_AR_Diff_No_SE_output<- vpcl %>%
      group_by(Expt_No, Sample_Type) 
    
    rm(vpcl, vipcal_vp)
    
    return(VPCL_AR_Diff_No_SE_output)
    
  }
  
  
  ####2.10 VPCL_AR_Diff_SE####
  
  simu_vpcl_ar_diff_se<-function(df = simu_df, variables = var){
    
    simu_df_AR_Diff<- simu_df_AR_Diff_function(df)
    
    vipcal_vp<- list()
    vpcl<- data.frame()
    
    
    for (expt in length_df) {
      
      print(expt)
      
      for (type in unique(simu_df_AR_Diff$Sample_Type)){
        df1 <- simu_df_AR_Diff %>% filter(Expt_No == expt,
                                          Sample_Type == type)
        print(paste(expt, type))
        
        
        try({
          p<- peaks_se(c(+10e+10, df1$mean, -10e+10),
                       c(0, df1$se, 0))-1
          v<- valleys_se(c(+10e+10, df1$mean, -10e+10),
                         c(0, df1$se, 0))-1
          
          if(identical(length(p), length(v))){
            print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
          }else{
            print("Number peaks and valleys are not the same. This will lead to erroneous viral production calculations")
          }
          
          if(length(p)==0){
            vp<- 0
            
          } else if (length(p)==1) {
            vp<- (df1$mean[p[1]] - df1$mean[v[1]])/(df1$Timepoint[p[1]] - df1$Timepoint[v[1]])
            
          } else if (length(p)==2) {
            vp<- ((df1$mean[p[1]] - df1$mean[v[1]])/(df1$Timepoint[p[1]] - df1$Timepoint[v[1]]) + 
                    (df1$mean[p[2]] - df1$mean[v[2]])/(df1$Timepoint[p[2]] - df1$Timepoint[v[2]]))/2
            
          } else if (length(p)==3) {
            vp<-  ((df1$mean[p[1]] - df1$mean[v[1]])/(df1$Timepoint[p[1]] - df1$Timepoint[v[1]]) + 
                     (df1$mean[p[2]] - df1$mean[v[2]])/(df1$Timepoint[p[2]] - df1$Timepoint[v[2]]) +
                     (df1$mean[p[3]] - df1$mean[v[3]])/(df1$Timepoint[p[3]] - df1$Timepoint[v[3]]))/3
          }
          
          if(length(p)==0){
            se<- 0
            
          } else if (length(p)==1) {
            se<- (df1$se[p[1]] + df1$se[v[1]])/(df1$Timepoint[p[1]] - df1$Timepoint[v[1]])
            
          } else if (length(p)==2) {
            se<- ((df1$se[p[1]] + df1$se[v[1]])/(df1$Timepoint[p[1]] - df1$Timepoint[v[1]]) + 
                    (df1$se[p[2]] + df1$se[v[2]])/(df1$Timepoint[p[2]] - df1$Timepoint[v[2]]))/2
            
          } else if (length(p)==3) {
            se<-  ((df1$se[p[1]] + df1$se[v[1]])/(df1$Timepoint[p[1]] - df1$Timepoint[v[1]]) + 
                     (df1$se[p[2]] + df1$se[v[2]])/(df1$Timepoint[p[2]] - df1$Timepoint[v[2]]) +
                     (df1$se[p[3]] + df1$se[v[3]])/(df1$Timepoint[p[3]] - df1$Timepoint[v[3]]))/3
          }
          
          
          vipcal<- c(expt, type, vp, se)
          
          vipcal_vp[[length(vipcal_vp)+1]] <- vipcal
        })
        
        
        
        print('done')
        
        
      }
      
      rm(df1, expt, type, vp)
      
    }
    
    {
      vpcl<- data.frame(t(sapply(vipcal_vp, c)))
      colnames(vpcl)<- c('Expt_No', 'Sample_Type', 'VP', 'VP_SE')
      vpcl[, c('Expt_No', 'VP', 'VP_SE')]<- lapply(vpcl[, c('Expt_No', 'VP', 'VP_SE')], as.numeric)
      vpcl$VP_Type<- 'VPCL_AR_Diff_SE'
      vpcl$VP_R_Squared<- NA
      
    }
    
    #calculate Lysogeny from here.
    
    VPCL_AR_Diff_SE_output<- vpcl %>%
      group_by(Expt_No, Sample_Type) 
    
    rm(vpcl, vipcal_vp)
    
    return(VPCL_AR_Diff_SE_output)
    
  }
  
  
  ####2.11 VPCL_AR_Diff_LMER_No_SE####
  
  simu_vpcl_ar_diff_lmer_no_se <- function(df = simu_df){
    
    
    vipcal_vp<- list()
    vpcl<- data.frame()
    
    for (expt in length_df){
      
      print(expt)
      
      
      df1<- df %>% filter(Expt_No == expt)
      
      try(df2<- simu_lmer_model(df1))
      
      
      try({
        if (exists("df2")){
          
          for (type in unique(df2$Sample_Type)){
            
            
            df3 <- df2 %>% filter(Sample_Type == type)
            print(paste(expt, type))
            
            
            try({
              p<- peaks(c(+10e+10, df3$Mean, -10e+10))-1
              v<- valleys(c(+10e+10, df3$Mean, -10e+10))-1
              
              if(identical(length(p), length(v))){
                print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
              }else{
                print("Number peaks and valleys are not the same. This will lead to erroneous viral production calculations")
              }
              
              if(length(p)==0){
                vp<- 0
                
              } else if (length(p)==1) {
                vp<- (df3$Mean[p[1]] - df3$Mean[v[1]])/(df3$Timepoint[p[1]] - df3$Timepoint[v[1]])
                
              } else if (length(p)==2) {
                vp<- ((df3$Mean[p[1]] - df3$Mean[v[1]])/(df3$Timepoint[p[1]] - df3$Timepoint[v[1]]) + 
                        (df3$Mean[p[2]] - df3$Mean[v[2]])/(df3$Timepoint[p[2]] - df3$Timepoint[v[2]]))/2
                
              } else if (length(p)==3) {
                vp<-  ((df3$Mean[p[1]] - df3$Mean[v[1]])/(df3$Timepoint[p[1]] - df3$Timepoint[v[1]]) + 
                         (df3$Mean[p[2]] - df3$Mean[v[2]])/(df3$Timepoint[p[2]] - df3$Timepoint[v[2]]) +
                         (df3$Mean[p[3]] - df3$Mean[v[3]])/(df3$Timepoint[p[3]] - df3$Timepoint[v[3]]))/3
              }
              
              
              vipcal<- c(expt, type, vp)
              
              vipcal_vp[[length(vipcal_vp)+1]] <- vipcal
            })
            
            
          }
          
        }      else{
          for (type in c("Diff", "VP", "VPC")){
            vipcal<- c(expt, type, NA, NA, NA)
            vipcal_vp[[length(vipcal_vp)+1]] <- vipcal 
          }
        }
      })
      
      
      
      print('done')
      rm(df1, expt,type)
    }
    
    
    
    
    vpcl<- data.frame(t(sapply(vipcal_vp, c)))
    colnames(vpcl)<- c('Expt_No', 'Sample_Type', 'VP')
    vpcl<- vpcl%>% mutate(across(c('Expt_No', 'VP'), as.numeric))
    vpcl$VP_SE<- NA
    vpcl$VP_R_Squared<- NA
    vpcl$VP_Type<- 'VPCL_AR_Diff_LMER_No_SE'
    
    VPCL_AR_Diff_LMER_No_SE_output<- vpcl %>%
      group_by(Expt_No, Sample_Type) 
    
    rm( vpcl, vipcal_vp)
    
    return(VPCL_AR_Diff_LMER_No_SE_output)
    
  }
  
  
  ####2.12 VPCL_AR_Diff_LMER_No_SE####
  
  simu_vpcl_ar_diff_lmer_se <- function(df = simu_df){
    
    
    vipcal_vp<- list()
    vpcl<- data.frame()
    
    for (expt in length_df){
      
      print(expt)
      
      
      df1<- df %>% filter(Expt_No == expt)
      
      try(df2<- simu_lmer_model(df1))
      
      
      try({
        if (exists("df2")){
          
          for (type in unique(df2$Sample_Type)){
            
            
            df3 <- df2 %>% filter(Sample_Type == type)
            print(paste(expt, type))
            
            
            try({
              p<- peaks_se(c(+10e+10, df3$Mean, -10e+10),
                           c(0, df3$SE, 0))-1
              v<- valleys_se(c(+10e+10, df3$Mean, -10e+10),
                             c(0, df3$SE, 0))-1 
              if(identical(length(p), length(v))){
                print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
              }else{
                print("Number peaks and valleys are not the same. This will lead to erroneous viral production calculations")
              }
              
              if(length(p)==0){
                vp<- 0
                
              } else if (length(p)==1) {
                vp<- (df3$Mean[p[1]] - df3$Mean[v[1]])/(df3$Timepoint[p[1]] - df3$Timepoint[v[1]])
                
              } else if (length(p)==2) {
                vp<- ((df3$Mean[p[1]] - df3$Mean[v[1]])/(df3$Timepoint[p[1]] - df3$Timepoint[v[1]]) + 
                        (df3$Mean[p[2]] - df3$Mean[v[2]])/(df3$Timepoint[p[2]] - df3$Timepoint[v[2]]))/2
                
              } else if (length(p)==3) {
                vp<-  ((df3$Mean[p[1]] - df3$Mean[v[1]])/(df3$Timepoint[p[1]] - df3$Timepoint[v[1]]) + 
                         (df3$Mean[p[2]] - df3$Mean[v[2]])/(df3$Timepoint[p[2]] - df3$Timepoint[v[2]]) +
                         (df3$Mean[p[3]] - df3$Mean[v[3]])/(df3$Timepoint[p[3]] - df3$Timepoint[v[3]]))/3
              }
              
              if(length(p)==0){
                se<- 0
                
              } else if (length(p)==1) {
                se<- (df3$SE[p[1]] + df3$SE[v[1]])/(df3$Timepoint[p[1]] - df3$Timepoint[v[1]])
                
              } else if (length(p)==2) {
                se<- ((df3$SE[p[1]] + df3$SE[v[1]])/(df3$Timepoint[p[1]] - df3$Timepoint[v[1]]) + 
                        (df3$SE[p[2]] + df3$SE[v[2]])/(df3$Timepoint[p[2]] - df3$Timepoint[v[2]]))/2
                
              } else if (length(p)==3) {
                se<-  ((df3$se[p[1]] + df3$SE[v[1]])/(df3$Timepoint[p[1]] - df3$Timepoint[v[1]]) + 
                         (df3$SE[p[2]] + df3$SE[v[2]])/(df3$Timepoint[p[2]] - df3$Timepoint[v[2]]) +
                         (df3$SE[p[3]] + df3$SE[v[3]])/(df3$Timepoint[p[3]] - df3$Timepoint[v[3]]))/3
              }
              
              vipcal<- c(expt, type, vp, se)
              
              vipcal_vp[[length(vipcal_vp)+1]] <- vipcal
            })
            
            
          }
          
        }      else{
          for (type in c("Diff", "VP", "VPC")){
            vipcal<- c(expt, type, NA, NA, NA)
            vipcal_vp[[length(vipcal_vp)+1]] <- vipcal 
          }
        }
      })
      
      
      
      print('done')
      rm(df1, expt,type)
    }
    
    
    
    vpcl<- data.frame(t(sapply(vipcal_vp, c)))
    vpcl<- data.frame(t(sapply(vipcal_vp, c)))
    colnames(vpcl)<- c('Expt_No', 'Sample_Type', 'VP', 'VP_SE')
    vpcl<- vpcl%>% mutate(across(c('Expt_No', 'VP', 'VP_SE'), as.numeric))
    vpcl$VP_Type<- 'VPCL_AR_Diff_LMER_SE'
    vpcl$VP_R_Squared<- NA
    
    VPCL_AR_Diff_LMER_SE_output<- vpcl %>%
      group_by(Expt_No, Sample_Type) 
    
    rm(vpcl)
    
    return(VPCL_AR_Diff_LMER_SE_output)
    
  }
}

####Function list####

simu_calc_funct_list<- list(simu_lm_ap, simu_lm_sr_avg, simu_lm_ar,
                            simu_lm_ar_diff, simu_lm_ar_diff_lmer,
                            simu_vpcl_sr_avg, simu_vpcl_ar_no_se, 
                            simu_vpcl_ar_se, simu_vpcl_ar_diff_no_se, 
                            simu_vpcl_ar_diff_se, simu_vpcl_ar_diff_lmer_no_se,
                            simu_vpcl_ar_diff_lmer_se)


#####Simulation function####

simulation_function_parallel<- function(output_df = output_df_par, output_csv = T, output_dir = './simulations_par',
                               output_csv_name = 'simu_output_df.csv'){
  
  if (output_csv == T){
    #Exporting the dataframe to a csv
    if (dir.exists('./simulations')){
      print("The directory already exists")
    }else{
      print('The directory doesnt exist. Creating one.')
      dir.create('./simulations')
      if(dir.exists('./simulations')){
        'The `simulations` directory was created'
      }
    }
  }
  
  #create output_df 
  {
    output_df<- data.frame(matrix(ncol = 6, nrow =0))
    colnames(output_df)<- c('Expt_No', 'Sample_Type', 'VP', 'VP_SE', 'VP_R_Squared', 'VP_Type')
    output_df <- output_df %>%
      mutate_at(c('Sample_Type', 
                  'VP_Type'), as.character) %>%
      mutate_at(c('Expt_No', 'VP', 'VP_SE', 
                  'VP_R_Squared'), as.numeric)
    
  }
  
  
  #running the simulation
  dir.create('./simulations/logs')
  
  #parallel run. 
  cl<- parallel::makeCluster(6)
  doParallel::registerDoParallel(cl)
  
  
    
    foreach::foreach (mtd = 1:12,
                      .combine = rbind){
      
      print(paste("Started analysis using method ", mtd))
      #sink(paste0('./simulations/logs/simu_output_log_method', mtd,'.txt'))
      try(simu_calc_funct_list[[mtd]](simu_df))
      #sink()
      
   
  
  
  
  if(output_csv == T){
    write.csv(output_df, paste0(output_dir, '/', output_csv_name),
              row.names = F)
    
  }
  
  return(output_df)
  
    }
    parallel::stopCluster(cl)
}



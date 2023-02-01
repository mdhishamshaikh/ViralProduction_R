slope_lm_sr<- function(df_sr){ #takes SR dataframe as an input
  
  lm_vp<- list()
  slope_df <- data.frame()
  for (location in unique(df_sr$Location)){
    for (expt_no in unique(df_sr$Expt_No)){
      for (depth in unique(df_sr$Depth)){
        for (time in unique(df_sr$Time_Range)){
          for (viruses in unique(df_sr$count)[4:7]){
            for (prod in unique(df_sr$Sample_Type)){
              for (rep in unique(df_sr$Replicate)){
                
                df2<- df_sr[df_sr$Location == location,]
                df2<- df2[df2$Expt_No == expt_no,]
                df2<- df2[df2$Depth == depth,]
                df2<- df2[df2$Time_Range == time,]
                df2<- df2[df2$count == viruses,]
                df2<- df2[df2$Sample_Type == prod,]
                df2<- df2[df2$Replicate == rep,]
                
                lm<- lm(data = df2, value ~ Timepoint)
                print(summary(lm))
                slope<- c(location, expt_no, depth, time, viruses, prod, rep, lm$coefficients[[2]])
                lm_vp[[length(lm_vp)+1]] <- slope  
              }
            }
          }
        }
      }
    }
  }
  
  slope_lm_sr<- data.frame(t(sapply(lm_vp, c)))
  colnames(slope_lm_sr)<- c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'Replicate', 'LM_SR_slope')
  return(lm_vp)
}
# works to caluclate serparet replicate slopes fro VP and VPC.

#Averaging them now and calculating Diff

df_avg_slope<- function(df, keep_0.22 = F) {
  #only works if we remove 0.22
  DF<- df[df$Sample_Type != '0.22',]
  
  DF<- DF %>%
    #gather('c_Bacteria', 'c_HNA', 'c_LNA', 'c_Viruses', 'c_V1', 'c_V2', 'c_V3', key="count", value="value") %>%
    group_by(Location, Expt_No, Depth, Time_Range, Population, Sample_Type) %>%
    summarise(n =n(), mean=mean(as.numeric(LM_SR_Slope)), sd=sd(LM_SR_Slope)) #calculating means and sd
  
  if ('VPC' %in% DF$Sample_Type){
    colnames_mean<- c("VP", "VPC")
    colnames_sd<- c("VP", "VPC")
    
  } else {
    colnames_mean<- c("VP")
    colnames_sd<- c("VP")
    
  }
  
  DF_mean<- select(DF, -c('sd')) %>% #splitting the dataframe cause I haven't figure out how to spread teh table without adding NAs
    spread('Sample_Type', 'mean')
  if (length(colnames_mean)==2){
    colnames(DF_mean)[7:8]<- colnames_mean
  } else if (length(colnames_mean)==1){
    colnames(DF_mean)[7]<- colnames_mean
  } 
  
  DF_sd<- select(DF, -c('mean')) %>%
    spread('Sample_Type', 'sd')
  if (length(colnames_sd)==2){
    colnames(DF_sd)[7:8]<- colnames_sd
  } else if (length(colnames_sd)==1){
    colnames(DF_sd)[7]<- colnames_sd
  } 
  
  if ('VPC' %in% DF$Sample_Type){
    DF_mean$Diff <- with(DF_mean, VPC-VP) #calculating Diff mean
    DF_mean<- pivot_longer(DF_mean, cols = c("VP", "VPC", "Diff"), names_to= 'Sample_Type', values_to='slope_mean_value')
    DF_sd$Diff <- with(DF_sd, VPC+VP) #Calculating Diff sd, which is addition of the other sds
    DF_sd<- pivot_longer(DF_sd, cols = c("VP", "VPC", "Diff"), names_to='Sample_Type', values_to= 'slope_sd_value')
  }
  
  DF<- merge(DF_mean, DF_sd, by= c('Location', 'Expt_No', 'Depth',
                                   'Time_Range', 'Population', 'Sample_Type', 'n'))
  rm('DF_mean', 'DF_sd', 'colnames_mean', 'colnames_sd')
  
  
  return(DF)
}



####Making VP and Diff dataframe for VIPCAL average calculations #####

vpcl_df<- df_sr(NJ1)
vpcl_df<- vpcl_df[vpcl_df$Microbe == "Viruses",]  
vpcl_df$Sample_Type<- factor(vpcl_df$Sample_Type, levels = c("VPC", "VP"))
vpcl_df$Timepoint2<- as.factor(vpcl_df$Timepoint)

vpcl_df[vpcl_df$Sample_Type == "VPC",]$Replicate<- replace(vpcl_df[vpcl_df$Sample_Type == "VPC",]$Replicate, vpcl_df[vpcl_df$Sample_Type == "VPC",]$Replicate == c("1","2", "3"), c("4", "5", "6"))
#vpcl_df[order(vpcl_df$Replicate),]

vpcl_df

df<- 

lyso_model<- lmer(value~Sample_Type*Timepoint2 + (1 | Replicate), data = lmer_data)
summary(lyso_model)

warm.emm_lyso<- emmeans(lyso_model, ~ Sample_Type|Timepoint2)
warm.emm_lyso

tw.emm_lyso<- pairs(warm.emm_lyso)
tw.emm_lyso



summary(tw.emm_lyso)$SE
as.numeric(unique(lmer_data$Timepoint))


unique(vpcl_df$count)
dataf<-data.frame(matrix(ncol=2,nrow=0))

colnames(dataf)<- c('mean', 'sd')
         
         for (pop in unique(vpcl_df$count)){
  print(pop)
  DF<- vpcl_df[vpcl_df$count== pop,]
  model<-lmer(value~Sample_Type*Timepoint2 + (1|Replicate), data = DF)
  summary(model)
  
  warp.emm<- emmeans(model, ~Sample_Type|Timepoint2)
  warp.emm
  tw.emm<- pairs(warp.emm)
  tw.emm
  
  a<- list(summary(tw.emm)$estimate)
  b<- t(summary(tw.emm)$SE)
  
  dataf[nrow(dataf)+1,]<- c(a,b)
  }#need to make a dataframe before calculating timeframe and vipcal


vpcl_dff<- vpcl_df %>% 
  expand(unique(vpcl_df$Location), unique(vpcl_df$Expt_No), unique(vpcl_df$Depth), unique(vpcl_df$Sample_Type), unique(vpcl_df$count), unique(vpcl_df$Timepoint))
vpcl_dff<- as.datafr


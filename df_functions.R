

#Creating four dataframe functions.
#1. separate replicates
#2. averaged replicate
#3. SR with time ranges
#4. AVG with time ranges

#separate replicate dataframe
df_sr<- function(df, keep_0.22 = F){
  DF<- df%>%
    select(c('Location', 'Expt_No', 'Depth', 'Sample_Type', 'Timepoint', 'Replicate',
             'c_Bacteria', 'c_HNA', 'c_LNA', 'c_Viruses', 'c_V1', 'c_V2', 'c_V3'))%>%
    gather('c_Bacteria', 'c_HNA', 'c_LNA', 'c_Viruses', 'c_V1', 'c_V2', 'c_V3', key="count", value="value")%>%
    mutate(Microbe = if_else(count == 'c_Bacteria' | count == 'c_HNA' | count == 'c_LNA', "Bacteria", "Viruses"))%>%
    arrange('Location', 'Expt_No', 'Depth', 'Sample_Type', 'Timepoint', 'Replicate','count')
  if (keep_0.22 == F){
    DF<- DF[DF$Sample_Type != '0.22',]
  }
  return(DF)
}
NJ1_sr<- df_sr(NJ1)



df_avg<- function(df, keep_0.22 = F) {
  #only works if we remove 0.22
  DF<- df[df$Sample_Type != '0.22',]
  
  DF<- DF %>%
    gather('c_Bacteria', 'c_HNA', 'c_LNA', 'c_Viruses', 'c_V1', 'c_V2', 'c_V3', key="count", value="value") %>%
    group_by(Location, Expt_No, Depth, Sample_Type, Timepoint, count ) %>%
    summarise(n =n(), mean=mean(value), sd=sd(value)) #calculating means and sd
  
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
    DF_mean$Diff <- with(DF_mean, VPC-VP) #calcualting Diff mean
    DF_mean<- pivot_longer(DF_mean, cols = c("VP", "VPC", "Diff"), names_to= 'Sample_Type', values_to='mean_value')
    DF_sd$Diff <- with(DF_sd, VPC+VP) #Calculating Diff sd, which is addition of the other sds
    DF_sd<- pivot_longer(DF_sd, cols = c("VP", "VPC", "Diff"), names_to='Sample_Type', values_to= 'sd_value')
  }
  
  DF<- merge(DF_mean, DF_sd, by= c('Location', 'Expt_No', 'Depth',
                                      'Timepoint', 'count', 'n', 'Sample_Type')) 
  
  
  rm('DF_mean', 'DF_sd', 'colnames_mean', 'colnames_sd')
  
  
  return(DF)
}

NJ1_avg<- df_avg(NJ1)

#df average replciates timepoints

df_avg_tp<- function(df, keep_0.22 = F){
  DF<- df_avg(df)
  
  TP<- unique(DF$Timepoint) #temporary variable
  
  colnames<- c() #temporary variable
  
  for (col in 2: length(TP)){
    colnames[length(colnames)+1] <- paste("T", TP[1], "_T", TP[col], sep = "")
  }
  
  DF<- DF%>%
    mutate("T0_T3" = case_when(Timepoint == '0' ~ "T0:T3",
                               Timepoint == '3' ~ "T0:T3"))%>%
    
    mutate("T0_T6" = case_when(Timepoint == '0' ~ "T0:T6",
                               Timepoint == '3' ~ "T0:T6", 
                               Timepoint == '6' ~ "T0:T6"))%>%
    
    mutate("T0_T17" = case_when(Timepoint == '0' ~ "T0:T17",
                                Timepoint == '3' ~ "T0:T17", 
                                Timepoint == '6' ~ "T0:T17",
                                Timepoint == '17' ~ "T0:T17"))%>%
    
    mutate("T0_T20" = case_when(Timepoint == '0' ~ "T0:T20",
                                Timepoint == '3' ~ "T0:T20",
                                Timepoint == '6' ~ "T0:T20",
                                Timepoint == '17' ~ "T0:T20",
                                Timepoint == '20' ~ "T0:T20"))%>%
    
    mutate("T0_T24" = case_when(Timepoint == '0' ~ "T0:T24",
                                Timepoint == '3' ~ "T0:T24",
                                Timepoint == '6' ~ "T0:T24",
                                Timepoint == '17' ~ "T0:T24",
                                Timepoint == '20' ~ "T0:T24",
                                Timepoint == '24' ~ "T0:T24")) %>%
    pivot_longer(cols = colnames, names_to = "Time_Range", values_to = "Time_Time")%>%
    drop_na()
  
  return(DF)
}

NJ1_avg_tp<- df_avg_tp(NJ1)

df_sr_tp<- function(df, keep_0.22 = F){
  DF<- df_sr(df)
  
  TP<- unique(DF$Timepoint) #temporary variable
  
  colnames<- c() #temporary variable
  
  for (col in 2: length(TP)){
    colnames[length(colnames)+1] <- paste("T", TP[1], "_T", TP[col], sep = "")
  }
  
  DF<- DF%>%
    mutate("T0_T3" = case_when(Timepoint == '0' ~ "T0:T3",
                               Timepoint == '3' ~ "T0:T3"))%>%
    
    mutate("T0_T6" = case_when(Timepoint == '0' ~ "T0:T6",
                               Timepoint == '3' ~ "T0:T6", 
                               Timepoint == '6' ~ "T0:T6"))%>%
    
    mutate("T0_T17" = case_when(Timepoint == '0' ~ "T0:T17",
                                Timepoint == '3' ~ "T0:T17", 
                                Timepoint == '6' ~ "T0:T17",
                                Timepoint == '17' ~ "T0:T17"))%>%
    
    mutate("T0_T20" = case_when(Timepoint == '0' ~ "T0:T20",
                                Timepoint == '3' ~ "T0:T20",
                                Timepoint == '6' ~ "T0:T20",
                                Timepoint == '17' ~ "T0:T20",
                                Timepoint == '20' ~ "T0:T20"))%>%
    
    mutate("T0_T24" = case_when(Timepoint == '0' ~ "T0:T24",
                                Timepoint == '3' ~ "T0:T24",
                                Timepoint == '6' ~ "T0:T24",
                                Timepoint == '17' ~ "T0:T24",
                                Timepoint == '20' ~ "T0:T24",
                                Timepoint == '24' ~ "T0:T24")) %>%
    pivot_longer(cols = colnames, names_to = "Time_Range", values_to = "Time_Time")%>%
    drop_na()
  
  return(DF)
}
NJ1_sr_tp<- df_sr_tp(NJ1)

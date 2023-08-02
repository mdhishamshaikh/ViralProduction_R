source("vp_calc_functions.R")
library(tidyverse)
library(ggsci)
library(svglite)
library(readxl)

####1.0 Dataframe functions####

#Time ranges

tp<- function(DF){
  
  
  TP<- unique(as.numeric(DF$Timepoint))
  colnames<- c()
  for (col in 2: length(TP)){
    a<- paste("T", TP[1], "_T", TP[col], sep = "")
    colnames[length(colnames)+1]<- a
  }
  colvalues<- c()
  for (col in 2: length(TP)){
    a<- paste("T", TP[1], ":T", TP[col], sep = "")
    #a<- paste0("R", col-1)
    colvalues[length(colvalues)+1]<- a
  }
  
  ncol<- ncol(DF)
  DF[colnames]<- NA
  
  
  DF[,ncol+1]<- case_when(DF$Timepoint == TP[1] ~ colvalues[1],
                          DF$Timepoint == TP[2] ~ colvalues[1])
  
  
  DF[,ncol+2]<- case_when(DF$Timepoint == TP[1] ~ colvalues[2],
                          DF$Timepoint == TP[2] ~ colvalues[2],
                          DF$Timepoint == TP[3] ~ colvalues[2])
  
  DF[,ncol+3]<- case_when(DF$Timepoint == TP[1] ~ colvalues[3],
                          DF$Timepoint == TP[2] ~ colvalues[3],
                          DF$Timepoint == TP[3] ~ colvalues[3],
                          DF$Timepoint == TP[4] ~ colvalues[3])
  
  DF[,ncol+4]<- case_when(DF$Timepoint == TP[1] ~ colvalues[4],
                          DF$Timepoint == TP[2] ~ colvalues[4],
                          DF$Timepoint == TP[3] ~ colvalues[4],
                          DF$Timepoint == TP[4] ~ colvalues[4],
                          DF$Timepoint == TP[5] ~ colvalues[4])
  
  DF[,ncol+5]<- case_when(DF$Timepoint == TP[1] ~ colvalues[5],
                          DF$Timepoint == TP[2] ~ colvalues[5],
                          DF$Timepoint == TP[3] ~ colvalues[5],
                          DF$Timepoint == TP[4] ~ colvalues[5],
                          DF$Timepoint == TP[5] ~ colvalues[5],
                          DF$Timepoint == TP[6] ~ colvalues[5])

  
  DF<- DF %>%
    pivot_longer(cols = colnames, names_to = "Time_Range", values_to = "Time_Time")%>%
    drop_na()
  
  rm('colnames', 'colvalues', 'TP', 'a', 'ncol')
  return(DF)
}

#separate replicate dataframe

df_SR<- function(df, keep_0.22 = F){
  DF<- df%>%
    select(c('Location', 'Expt_No', 'Depth', 'Sample_Type', 'Timepoint', 'Replicate',
             'c_Bacteria', 'c_HNA', 'c_LNA', 'c_Viruses', 'c_V1', 'c_V2', 'c_V3'))%>%
    gather('c_Bacteria', 'c_HNA', 'c_LNA', 'c_Viruses', 'c_V1', 'c_V2', 'c_V3', key="Population", value="count")%>%
    mutate(Microbe = if_else(Population == 'c_Bacteria' | Population == 'c_HNA' | Population == 'c_LNA', "Bacteria", "Viruses"))%>%
    arrange('Location', 'Expt_No', 'Depth', 'Sample_Type', as.numeric(Timepoint), 'Replicate','Population')
  if (keep_0.22 == F){
    DF<- DF[DF$Sample_Type != '0.22',]
  }
  return(DF)
}

df_AVG<- function(df, keep_0.22 = F) {
  #only works if we remove 0.22
  DF<- df[df$Sample_Type != '0.22',]
  
  DF<- DF %>%
    gather('c_Bacteria', 'c_HNA', 'c_LNA', 'c_Viruses', 'c_V1', 'c_V2', 'c_V3', key="Population", value="count") %>%
    group_by(Location, Expt_No, Depth, Sample_Type, Timepoint, Population ) %>%
    summarise(n =n(), mean=mean(count), se=plotrix::std.error(count)) #calculating means and sd
 
  # I think this part is not necessary, since the spread function later on automatically gives VP and VPC as colnames for mean and se resp.
  if ('VPC' %in% DF$Sample_Type){
    colnames_mean<- c("VP", "VPC")
    colnames_se<- c("VP", "VPC")
    
  } else {
    colnames_mean<- c("VP")
    colnames_se<- c("VP")
    
  }
  
  DF_mean<- select(DF, -c('se')) %>% #splitting the dataframe cause I haven't figure out how to spread teh table without adding NAs
    spread('Sample_Type', 'mean')
  if (length(colnames_mean)==2){
    colnames(DF_mean)[7:8]<- colnames_mean
  } else if (length(colnames_mean)==1){
    colnames(DF_mean)[7]<- colnames_mean
  } 
  
  DF_se<- select(DF, -c('mean')) %>%
    spread('Sample_Type', 'se')
  if (length(colnames_se)==2){
    colnames(DF_se)[7:8]<- colnames_se
  } else if (length(colnames_se)==1){
    colnames(DF_se)[7]<- colnames_se
  } 
  
  if ('VPC' %in% DF$Sample_Type){
    DF_mean$Diff <- with(DF_mean, VPC-VP) #calculating Diff mean
    DF_mean<- pivot_longer(DF_mean, cols = c("VP", "VPC", "Diff"), names_to= 'Sample_Type', values_to='mean')
    DF_se$Diff <- with(DF_se, VPC+VP) #Calculating Diff se, which is addition of the other ses
    DF_se<- pivot_longer(DF_se, cols = c("VP", "VPC", "Diff"), names_to='Sample_Type', values_to= 'se')
  }
  
  DF<- merge(DF_mean, DF_se, by= c('Location', 'Expt_No', 'Depth',
                                   'Timepoint', 'Population', 'n', 'Sample_Type')) %>%
    mutate(Microbe = if_else(Population == 'c_Bacteria' | Population == 'c_HNA' | Population == 'c_LNA', "Bacteria", "Viruses"))%>%
    mutate(Subgroup = if_else(Population == 'c_Bacteria' | Population == 'c_Viruses', "Parent", "Subgroup"))
  rm('DF_mean', 'DF_se', 'colnames_mean', 'colnames_se')
  
  DF<- DF%>%
        arrange('Location',
            'Expt_No',
            'Depth',
            'Sample_Type',
            'Population',
            as.numeric(Timepoint))
   
  
  return(DF)
}


#df average replicates timepoints
# I don't think these extra two functions are necessary, you can just add one line of DF <- tp(DF) to add the timepoints i think?
# If I try this myself I get exactly the same dataframe as with this code

df_avg_tp<- function(df, keep_0.22 = F){
  DF3<- df_AVG(df)
  
  df_list<- list()
  
  for (location in unique(DF3$Location)){
    for (expt_no in unique(DF3$Expt_No)){
      for (depth in unique(DF3$Depth)){
        DF2<- DF3%>% dplyr::filter(Location == location & Expt_No == expt_no & Depth == depth )
        DF2<- tp(DF2)
        df_list[[length(df_list)+1]]<- DF2
        rm(DF2)
      }
    }
  }
  
  DF_tp<-  data.table::rbindlist(df_list)
  return(DF_tp)
}

df_sr_tp<- function(df, keep_0.22 = F){
  DF3<- df_SR(df)
  
  df_list<- list()
  
  # Are these forloops necessary because in each of these columns has only one unique value?
  for (location in unique(DF3$Location)){
    for (expt_no in unique(DF3$Expt_No)){
      for (depth in unique(DF3$Depth)){
        DF2<- DF3%>% dplyr::filter(Location == location & Expt_No == expt_no & Depth == depth )
        DF2<- tp(DF2)
        df_list[[length(df_list)+1]]<- DF2
        rm(DF2)
      }
    }
  }
  
  DF_tp<-  data.table::rbindlist(df_list)
  return(DF_tp)
}

####2.0 Peaks and Valleys####

peaks<- function(values){
  list<- c()
  
  for (len in 1:(length(values)-1)){
    d<- sign(values[len+1] - values[len])
    #print(d)
    list[[length(list)+1]]<- d
  }
  which(diff(as.numeric(list))<0) +2 -1
  
}

valleys<- function(values){
  list<- c()
  for (len in 1:(length(values)-1)){
    d<- sign(values[len+1] - values[len])
    #print(d)
    list[[length(list)+1]]<- d
  }
  which(diff(as.numeric(list))>0) +2 -1
}

peaks_se<- function(values, se){
  list<- c()
  
  for (len in 1:(length(values)-1)){
    d<- sign((values[len+1] - se[len+1]) - (values[len] + se[len]))
    #print(d)
    list[[length(list)+1]]<- d
    #print(list)
  }
  which(diff(as.numeric(list))<0) +2 -1
  
}

valleys_se<- function(values, se){
  list<- c()
  
  for (len in 1:(length(values)-1)){
    d<- sign((values[len+1] - se[len+1]) - (values[len] + se[len]))
    #print(d)
    list[[length(list)+1]]<- d
  }
  which(diff(as.numeric(list))>0) +2 -1
  
}

####3.0

#####Slope Functions####
slope_lm_sr<- function(df_sr){ #takes SR dataframe as an input
  
  lm_vp<- list()
  slope_lm_sr_df<- data.frame()
  for (location in unique(df_sr$Location)){
    for (expt_no in unique(df_sr$Expt_No)){
      for (depth in unique(df_sr$Depth)){
        for (viruses in unique(df_sr[df_sr$Microbe == "Viruses",]$Population)){
          for (prod in unique(df_sr$Sample_Type)){
            for (rep in unique(df_sr$Replicate)){
              
              
              df2<- df_sr %>% 
                filter(Location == location,
                       Expt_No == expt_no,
                       Depth == depth,
                       Population == viruses,
                       Sample_Type == prod,
                       Replicate == rep)
              
              
              for (time in unique(df2$Time_Range)){
                
                df3<- df2 %>% 
                  filter(Time_Range == time)
                
                try({  
                  lm<- summary(lm(data = df3, count ~ as.numeric(Timepoint)))
                  
                  #print(lm)
                  
                  slope<- c(location, expt_no, depth, time, viruses, prod, rep, lm$coefficients[c(2,4)], lm$r.squared)
                  
                  #print(lm$coefficients[,3])
                  #print(lm$coefficients[[2]])
                  lm_vp[[length(lm_vp)+1]] <- slope 
                  
                })
              }
            }
          }
        }
      }
    }
  }
  
  slope_lm_sr_df<- data.frame(t(sapply(lm_vp, c)))
  colnames(slope_lm_sr_df)<- c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'Replicate', 'VP_Slope', 'VP_SE', 'VP_R_Squared')
  slope_lm_sr_df[, c('VP_Slope', 'VP_SE', 'VP_R_Squared')]<- lapply(slope_lm_sr_df[, c('VP_Slope', 'VP_SE', 'VP_R_Squared')], as.numeric)
  return(slope_lm_sr_df)
}

slope_lm_avg<- function(df_avg){ #takes AR dataframe as an input
  
  lm_vp<- list()
  slope_lm_avg_df<- data.frame()
  for (location in unique(df_avg$Location)){
    for (expt_no in unique(df_avg$Expt_No)){
      for (depth in unique(df_avg$Depth)){
        for (viruses in unique(df_avg[df_avg$Microbe == "Viruses",]$Population)){
          for (prod in unique(df_avg$Sample_Type)){
            
            
            df2<- df_avg %>% 
              filter(Location == location,
                     Expt_No == expt_no,
                     Depth == depth,
                     Population == viruses,
                     Sample_Type == prod)
            
            
            for (time in unique(df2$Time_Range)){
              df3<- df2 %>% 
                filter(Time_Range == time)
              
              try({
                
                lm<- summary(lm(data = df3, mean ~ Timepoint))
                
                #print(lm)
                
                slope<- c(location, expt_no, depth, time, viruses, prod, lm$coefficients[c(2,4)], lm$r.squared)
                
                lm_vp[[length(lm_vp)+1]] <- slope  
                
              })
            }
          }
        }
      }
    }
  }
  
  
  slope_lm_avg_df<- data.frame(t(sapply(lm_vp, c)))
  colnames(slope_lm_avg_df)<- c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'VP_Slope', 'VP_SE', 'VP_R_Squared')
  slope_lm_avg_df[, c('VP_Slope', 'VP_SE', 'VP_R_Squared')]<- lapply(slope_lm_avg_df[, c('VP_Slope', 'VP_SE', 'VP_R_Squared')], as.numeric)
  return(slope_lm_avg_df)
}

slope_all_points<- function(df_sr){ #takes df_sr as the input
  lm_vp<- list()
  slope_lm_all_points_df<- data.frame()
  for (location in unique(df_sr$Location)){
    for (expt_no in unique(df_sr$Expt_No)){
      for (depth in unique(df_sr$Depth)){
        for (viruses in unique(df_sr[df_sr$Microbe == "Viruses",]$Population)){
          for (prod in unique(df_sr$Sample_Type)){
            
            df2<- df_sr %>% 
              filter(Location == location,
                     Expt_No == expt_no,
                     Depth == depth,
                     Population == viruses,
                     Sample_Type == prod)
            
            
            for (time in unique(df2$Time_Range)){
              
              df3<- df2 %>% 
                filter(Time_Range == time)
              
              try({
                
                lm<- summary(lm(data = df3, count ~ as.numeric(Timepoint)))
                
                #print(lm)
                
                slope<- c(location, expt_no, depth, time, viruses, prod, lm$coefficients[c(2,4)], lm$r.squared)
                lm_vp[[length(lm_vp)+1]] <- slope  
                
                
              })
            }
          }
        }
      }
    }
  }
  
  slope_lm_all_points_df<- data.frame(t(sapply(lm_vp, c)))
  colnames(slope_lm_all_points_df)<- c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'VP_Slope', 'VP_SE', 'VP_R_Squared')
  slope_lm_all_points_df[, c('VP_Slope', 'VP_SE', 'VP_R_Squared')]<- lapply(slope_lm_all_points_df[, c('VP_Slope', 'VP_SE', 'VP_R_Squared')], as.numeric)
  return(slope_lm_all_points_df)
}

slope_lm_ar_diff_lmer<- function(df_sr){ #takes SR dataframe as an input
  
  lm_vp<- list()
  slope_lm_ar_df<- data.frame()
  for (location in unique(df_sr$Location)){
    for (expt_no in unique(df_sr$Expt_No)){
      for (depth in unique(df_sr$Depth)){
        
        for (viruses in unique(df_sr[df_sr$Microbe == "Viruses",]$Population)){
          
          
          df2<- df_sr %>% 
            filter(Location == location,
                   Expt_No == expt_no,
                   Depth == depth,
                   Population == viruses
            )
          
          
          for (time in unique(df2$Time_Range)){
            
            df3<- df2 %>% 
              filter(Time_Range == time)
            
            try(df4<- lmer_model(df3))#the output is a dataframe for VP, VPC and Diff together
            try({
              if (exists("df4")){
                for (prod in unique(df4$Sample_Type)){
                  
                  lm<- summary(lm(data = df4[df4$Sample_Type == prod,], Mean ~ as.numeric(Timepoint)))
                  print(lm)
                  slope<- c(location, expt_no, depth, time, viruses, prod, lm$coefficients[c(2,4)], lm$r.squared)
                  #print(lm$coefficients[,3])
                  #print(lm$coefficients[[2]])
                  lm_vp[[length(lm_vp)+1]] <- slope 
                  
                }
              } else{
                for (prod in c("Diff", "VP", "VPC")){
                  slope<- c(location, expt_no, depth, time, viruses, prod, NA, NA, NA)
                  lm_vp[[length(lm_vp)+1]] <- slope 
                }
              }
            })
            
            
            
          }     
        }
      }
    }
  }
  
  slope_lm_ar_df<- data.frame(t(sapply(lm_vp, c)))
  colnames(slope_lm_ar_df)<- c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'VP_Slope', 'VP_SE', 'VP_R_Squared')
  slope_lm_ar_df[, c('VP_Slope', 'VP_SE', 'VP_R_Squared')]<- lapply(slope_lm_ar_df[, c('VP_Slope', 'VP_SE', 'VP_R_Squared')], as.numeric)
  return(slope_lm_ar_df)
  rm(df3)
}

####VIPCAL Functions####
vipcal_sr<- function(df_sr){ #takes SR dataframe as an input
  
  vipcal_vp<- list()
  vpcl_vp_df<- data.frame()
  for (location in unique(df_sr$Location)){
    for (expt_no in unique(df_sr$Expt_No)){
      for (depth in unique(df_sr$Depth)){
        for (viruses in unique(df_sr[df_sr$Microbe == "Viruses",]$Population)){
          for (prod in unique(df_sr$Sample_Type)){
            for (rep in unique(df_sr$Replicate)){
              
              
              df2<- df_sr %>% 
                filter(Location == location,
                       Expt_No == expt_no,
                       Depth == depth,
                       Population == viruses,
                       Sample_Type == prod,
                       Replicate == rep)
              
              
              for (time in unique(df2$Time_Range)){
                
                df3<- df2 %>% 
                  filter(Time_Range == time)
                
                try({
                  p<- peaks(c(+10e+10, df3$count, -10e+10))-1
                  v<- valleys(c(+10e+10, df3$count, -10e+10))-1
                  
                  if(identical(length(p), length(v))){
                    print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
                  }else{
                    print("Number peaks and valleys are not the same. This will lead to erroneous viral production calculations")
                  }
                  
                  if(length(p)==0){
                    vp<- 0
                    
                  } else if (length(p)==1) {
                    vp<- (df3$count[p[1]] - df3$count[v[1]])/(df3$Timepoint[p[1]] - df3$Timepoint[v[1]])
                    
                  } else if (length(p)==2) {
                    vp<- ((df3$count[p[1]] - df3$count[v[1]])/(df3$Timepoint[p[1]] - df3$Timepoint[v[1]]) + 
                            (df3$count[p[2]] - df3$count[v[2]])/(df3$Timepoint[p[2]] - df3$Timepoint[v[2]]))/2
                    
                  } else if (length(p)==3) {
                    vp<-  ((df3$count[p[1]] - df3$count[v[1]])/(df3$Timepoint[p[1]] - df3$Timepoint[v[1]]) + 
                             (df3$count[p[2]] - df3$count[v[2]])/(df3$Timepoint[p[2]] - df3$Timepoint[v[2]]) +
                             (df3$count[p[3]] - df3$count[v[3]])/(df3$Timepoint[p[3]] - df3$Timepoint[v[3]]))/3
                  }
                  
                  
                  vipcal<- c(location, expt_no, depth, time, viruses, prod, rep, vp)
                  
                  vipcal_vp[[length(vipcal_vp)+1]] <- vipcal
                })
              }
            }
          }
        }
      }
    }
  }
  
  vpcl_vp_df<- data.frame(t(sapply(vipcal_vp, c)))
  colnames(vpcl_vp_df)<- c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'Replicate', 'VP')
  vpcl_vp_df$VP<- as.numeric(vpcl_vp_df$VP)
  return(vpcl_vp_df)
}

vipcal_avg_se<- function(df_avg){ #takes avg_tp dataframe as an input
  
  vipcal_vp<- list()
  vpcl_vp_df<- data.frame()
  for (location in unique(df_avg$Location)){
    for (expt_no in unique(df_avg$Expt_No)){
      for (depth in unique(df_avg$Depth)){
        for (viruses in unique(df_avg[df_avg$Microbe == "Viruses",]$Population)){
          for (prod in unique(df_avg$Sample_Type)){
            
            
            
            df2<- df_avg %>% 
              filter(Location == location,
                     Expt_No == expt_no,
                     Depth == depth,
                     Population == viruses,
                     Sample_Type == prod)
            
            
            for (time in unique(df2$Time_Range)){
              
              df3<- df2 %>% 
                filter(Time_Range == time)
              
              try({
                p<- peaks_se(c(+10e+10, df3$mean, -10e+10),
                             c(0, df3$se, 0))-1
                v<- valleys_se(c(+10e+10, df3$mean, -10e+10),
                               c(0, df3$se, 0))-1
                
                if(identical(length(p), length(v))){
                  print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
                }else{
                  print("Number peaks and valleys are not the same. This will lead to erroneous viral production calculations")
                }
                
                if(length(p)==0){
                  vp<- 0
                  
                } else if (length(p)==1) {
                  vp<- (df3$mean[p[1]] - df3$mean[v[1]])/(df3$Timepoint[p[1]] - df3$Timepoint[v[1]])
                  
                } else if (length(p)==2) {
                  vp<- ((df3$mean[p[1]] - df3$mean[v[1]])/(df3$Timepoint[p[1]] - df3$Timepoint[v[1]]) + 
                          (df3$mean[p[2]] - df3$mean[v[2]])/(df3$Timepoint[p[2]] - df3$Timepoint[v[2]]))/2
                  
                } else if (length(p)==3) {
                  vp<-  ((df3$mean[p[1]] - df3$mean[v[1]])/(df3$Timepoint[p[1]] - df3$Timepoint[v[1]]) + 
                           (df3$mean[p[2]] - df3$mean[v[2]])/(df3$Timepoint[p[2]] - df3$Timepoint[v[2]]) +
                           (df3$mean[p[3]] - df3$mean[v[3]])/(df3$Timepoint[p[3]] - df3$Timepoint[v[3]]))/3
                }
                
                if(length(p)==0){
                  se<- 0
                  
                } else if (length(p)==1) {
                  se<- (df3$se[p[1]] + df3$se[v[1]])/(df3$Timepoint[p[1]] - df3$Timepoint[v[1]])
                  
                } else if (length(p)==2) {
                  se<- ((df3$se[p[1]] + df3$se[v[1]])/(df3$Timepoint[p[1]] - df3$Timepoint[v[1]]) + 
                          (df3$se[p[2]] + df3$se[v[2]])/(df3$Timepoint[p[2]] - df3$Timepoint[v[2]]))/2
                  
                } else if (length(p)==3) {
                  se<-  ((df3$se[p[1]] + df3$se[v[1]])/(df3$Timepoint[p[1]] - df3$Timepoint[v[1]]) + 
                           (df3$se[p[2]] + df3$se[v[2]])/(df3$Timepoint[p[2]] - df3$Timepoint[v[2]]) +
                           (df3$se[p[3]] + df3$se[v[3]])/(df3$Timepoint[p[3]] - df3$Timepoint[v[3]]))/3
                }
                
                
                vipcal<- c(location, expt_no, depth, time, viruses, prod, vp,se)
                
                vipcal_vp[[length(vipcal_vp)+1]] <- vipcal  
              })
            }
          }
        }
      }
    }
  }
  
  vpcl_vp_df<- data.frame(t(sapply(vipcal_vp, c)))
  colnames(vpcl_vp_df)<- c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population', 'Sample_Type',  'VP', 'VP_SE')
  vpcl_vp_df[, c('VP', 'VP_SE')]<- lapply(vpcl_vp_df[, c('VP', 'VP_SE')], as.numeric)
  
  
  return(vpcl_vp_df)
}

vipcal_avg<- function(df_avg){ #takes avg_tp dataframe as an input
  
  vipcal_vp<- list()
  vpcl_vp_df<- data.frame()
  for (location in unique(df_avg$Location)){
    for (expt_no in unique(df_avg$Expt_No)){
      for (depth in unique(df_avg$Depth)){
        for (viruses in unique(df_avg[df_avg$Microbe == "Viruses",]$Population)){
          for (prod in unique(df_avg$Sample_Type)){
            
            
            df2<- df_avg %>% 
              filter(Location == location,
                     Expt_No == expt_no,
                     Depth == depth,
                     Population == viruses,
                     Sample_Type == prod)
            
            
            for (time in unique(df2$Time_Range)){
              
              df3<- df2 %>% 
                filter(Time_Range == time)
              
              try({
                p<- peaks(c(+10e+10, df3$mean, -10e+10))-1
                v<- valleys(c(+10e+10, df3$mean, -10e+10))-1
                
                if(identical(length(p), length(v))){
                  print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
                }else{
                  print("Number peaks and valleys are not the same. This will lead to erroneous viral production calculations")
                }
                
                if(length(p)==0){
                  vp<- 0
                  
                } else if (length(p)==1) {
                  vp<- (df3$mean[p[1]] - df3$mean[v[1]])/(df3$Timepoint[p[1]] - df3$Timepoint[v[1]])
                  
                } else if (length(p)==2) {
                  vp<- ((df3$mean[p[1]] - df3$mean[v[1]])/(df3$Timepoint[p[1]] - df3$Timepoint[v[1]]) + 
                          (df3$mean[p[2]] - df3$mean[v[2]])/(df3$Timepoint[p[2]] - df3$Timepoint[v[2]]))/2
                  
                } else if (length(p)==3) {
                  vp<-  ((df3$mean[p[1]] - df3$mean[v[1]])/(df3$Timepoint[p[1]] - df3$Timepoint[v[1]]) + 
                           (df3$mean[p[2]] - df3$mean[v[2]])/(df3$Timepoint[p[2]] - df3$Timepoint[v[2]]) +
                           (df3$mean[p[3]] - df3$mean[v[3]])/(df3$Timepoint[p[3]] - df3$Timepoint[v[3]]))/3
                }
                
                
                vipcal<- c(location, expt_no, depth, time, viruses, prod, vp)
                
                vipcal_vp[[length(vipcal_vp)+1]] <- vipcal  
              })
            }
          }
        }
      }
    }
  }
  
  vpcl_vp_df<- data.frame(t(sapply(vipcal_vp, c)))
  colnames(vpcl_vp_df)<- c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population', 'Sample_Type',  'VP')
  vpcl_vp_df$VP<- as.numeric(vpcl_vp_df$VP)
  return(vpcl_vp_df)
}

vipcal_sr_diff_SE<- function(df_sr){ #takes SR dataframe as an input
  
  vipcal_vp<- list()
  vpcl_vp_df<- data.frame()
  for (location in unique(df_sr$Location)){
    for (expt_no in unique(df_sr$Expt_No)){
      for (depth in unique(df_sr$Depth)){
        for (viruses in unique(df_sr[df_sr$Microbe == "Viruses",]$Population)){
          
          
          df2<- df_sr %>% 
            filter(Location == location,
                   Expt_No == expt_no,
                   Depth == depth,
                   Population == viruses
            )
          
          
          for (time in unique(df2$Time_Range)){
            
            df3<- df2 %>% 
              filter(Time_Range == time)
            
            try(df4<- lmer_model(df3))#the output is a dataframe for VP, VPC and Diff together
            try({
              if (exists("df3")){
                for (prod in unique(df4$Sample_Type)){
                  df5<-df4 %>% filter(Sample_Type == prod)
                  
                  p<- peaks_se(c(+10e+10, df5$Mean, -10e+10),
                               c(0, df5$SE, 0))-1
                  v<- valleys_se(c(+10e+10, df5$Mean, -10e+10),
                                 c(0, df5$SE, 0))-1
                  
                  if(identical(length(p), length(v))){
                    print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
                  }else{
                    print("Number peaks and valleys are not the same. This will lead to erroneous viral production calculations")
                  }
                  print(paste(location, expt_no, depth, time, viruses, prod))
                  print(p)
                  print(v)
                  if(length(p)==0){
                    vp<- 0
                    
                  } else if (length(p)==1) {
                    vp<- (df5$Mean[p[1]] - df5$Mean[v[1]])/(df5$Timepoint[p[1]] - df5$Timepoint[v[1]])
                    
                  } else if (length(p)==2) {
                    vp<- ((df5$Mean[p[1]] - df5$Mean[v[1]])/(df5$Timepoint[p[1]] - df5$Timepoint[v[1]]) + 
                            (df5$Mean[p[2]] - df5$Mean[v[2]])/(df5$Timepoint[p[2]] - df5$Timepoint[v[2]]))/2
                    
                  } else if (length(p)==3) {
                    vp<-  ((df5$Mean[p[1]] - df5$Mean[v[1]])/(df5$Timepoint[p[1]] - df5$Timepoint[v[1]]) + 
                             (df5$Mean[p[2]] - df5$Mean[v[2]])/(df5$Timepoint[p[2]] - df5$Timepoint[v[2]]) +
                             (df5$Mean[p[3]] - df5$Mean[v[3]])/(df5$Timepoint[p[3]] - df5$Timepoint[v[3]]))/3
                  }
                  
                  
                  if(length(p)==0){
                    se<- 0
                    
                  } else if (length(p)==1) {
                    se<- (df5$SE[p[1]] + df5$SE[v[1]])/(df5$Timepoint[p[1]] - df5$Timepoint[v[1]])
                    
                  } else if (length(p)==2) {
                    se<- ((df5$SE[p[1]] + df5$SE[v[1]])/(df5$Timepoint[p[1]] - df5$Timepoint[v[1]]) + 
                            (df5$SE[p[2]] + df5$SE[v[2]])/(df5$Timepoint[p[2]] - df5$Timepoint[v[2]]))/2
                    
                  } else if (length(p)==3) {
                    se<-  ((df5$SE[p[1]] + df5$SE[v[1]])/(df5$Timepoint[p[1]] - df5$Timepoint[v[1]]) + 
                             (df5$SE[p[2]] + df5$SE[v[2]])/(df5$Timepoint[p[2]] - df5$Timepoint[v[2]]) +
                             (df5$SE[p[3]] + df5$SE[v[3]])/(df5$Timepoint[p[3]] - df5$Timepoint[v[3]]))/3
                  }
                  
                  vipcal<- c(location, expt_no, depth, time, viruses, prod, vp, se)
                  
                  vipcal_vp[[length(vipcal_vp)+1]] <- vipcal  
                  
                  rm(p,v)
                  
                  
                }
              } else{
                for (prod in unique(df3$Sample_Type)){
                  vipcal<- c(location, expt_no, depth, time, viruses, prod, NA, NA)
                  vipcal_vp[[length(vipcal_vp)+1]] <- vipcal 
                }
              }
              
            })
            
          }     
        }
      }
    }
  }
  
  vpcl_vp_df<- data.frame(t(sapply(vipcal_vp, c)))
  colnames(vpcl_vp_df)<- c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'VP', 'VP_SE')
  vpcl_vp_df$VP <- as.numeric(vpcl_vp_df$VP) 
  vpcl_vp_df$VP_SE <- as.numeric(vpcl_vp_df$VP_SE)
  return(vpcl_vp_df)
  rm(df3)
  rm(df4)
}

vipcal_sr_diff_no_SE<- function(df_sr){ #takes SR dataframe as an input
  
  vipcal_vp<- list()
  vpcl_vp_df<- data.frame()
  for (location in unique(df_sr$Location)){
    for (expt_no in unique(df_sr$Expt_No)){
      for (depth in unique(df_sr$Depth)){
        for (viruses in unique(df_sr[df_sr$Microbe == "Viruses",]$Population)){
          
          
          
          df2<- df_sr %>% 
            filter(Location == location,
                   Expt_No == expt_no,
                   Depth == depth,
                   Population == viruses)
          
          for (time in unique(df2$Time_Range)){
            
            df3<- df2 %>% 
              filter(Time_Range == time)
            
            try(df4<- lmer_model(df3))#the output is a dataframe for VP, VPC and Diff together
            try({
              if (exists("df4")){
                for (prod in unique(df4$Sample_Type)){
                  df5<-df4 %>% filter(Sample_Type == prod)
                  
                  p<- peaks(c(+10e+10, df5$Mean, -10e+10))-1
                  v<- valleys(c(+10e+10, df5$Mean, -10e+10))-1
                  
                  if(identical(length(p), length(v))){
                    print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
                  }else{
                    print("Number peaks and valleys are not the same. This will lead to erroneous viral production calculations")
                  }
                  
                  if(length(p)==0){
                    vp<- 0
                    
                  } else if (length(p)==1) {
                    vp<- (df5$Mean[p[1]] - df5$Mean[v[1]])/(df5$Timepoint[p[1]] - df5$Timepoint[v[1]])
                    
                  } else if (length(p)==2) {
                    vp<- ((df5$Mean[p[1]] - df5$Mean[v[1]])/(df5$Timepoint[p[1]] - df5$Timepoint[v[1]]) + 
                            (df5$Mean[p[2]] - df5$Mean[v[2]])/(df5$Timepoint[p[2]] - df5$Timepoint[v[2]]))/2
                    
                  } else if (length(p)==3) {
                    vp<-  ((df5$Mean[p[1]] - df5$Mean[v[1]])/(df5$Timepoint[p[1]] - df5$Timepoint[v[1]]) + 
                             (df5$Mean[p[2]] - df5$Mean[v[2]])/(df5$Timepoint[p[2]] - df5$Timepoint[v[2]]) +
                             (df5$Mean[p[3]] - df5$Mean[v[3]])/(df5$Timepoint[p[3]] - df5$Timepoint[v[3]]))/3
                  }
                  
                  
                  vipcal<- c(location, expt_no, depth, time, viruses, prod, vp)
                  
                  vipcal_vp[[length(vipcal_vp)+1]] <- vipcal  
                  
                  
                  
                  
                }
              } else{
                for (prod in unique(df3$Sample_Type)){
                  vipcal<- c(location, expt_no, depth, time, viruses, prod, NA)
                  vipcal_vp[[length(vipcal_vp)+1]] <- vipcal 
                }
              }
            })
            
            rm(p,v)
            
          }     
        }
      }
    }
  }
  
  vpcl_vp_df<- data.frame(t(sapply(vipcal_vp, c)))
  colnames(vpcl_vp_df)<- c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'VP')
  vpcl_vp_df$VP <- as.numeric(vpcl_vp_df$VP) 
  return(vpcl_vp_df)
  rm(df3, df4)
  
}

####To calculate lysogeny from all points data####
calc_diff_lm_AP<- function(df){
  
  colnames(df)[colnames(df) == 'VP_Slope'] <- 'VP'
  
  VP<- df[df$Sample_Type == "VP",]
  VPC<- df[df$Sample_Type == "VPC",]
  Diff<- as.data.frame(cbind(VP[,c('Location', 'Expt_No', 'Depth','Sample_Type', 'Population', 'Time_Range')], VPC[,'VP'] -VP[,'VP'], VPC[,'VP_SE'] + VP[,'VP_SE']))
  Diff[,4]<- "Diff"
  colnames(Diff)[7]<- 'VP'
  colnames(Diff)[8]<- 'VP_SE'
  Diff$VP_R_Squared <- NA
  
  df<- full_join(df, Diff, by = NULL)
  
  return(df)
}

calc_diff_lm_SR<- function(df){
  
  VP<- df[df$Sample_Type == "VP",]
  VPC<- df[df$Sample_Type == "VPC",]
  Diff<- as.data.frame(cbind(VP[,c('Location', 'Expt_No', 'Depth','Sample_Type', 'Population', 'Time_Range')],VPC[,'VP'] -VP[,'VP'], VPC[,'VP_SE'] + VP[,'VP_SE']))
  Diff[,4]<- "Diff"
  colnames(Diff)[7]<- 'LM_SR_Slope_Mean'
  colnames(Diff)[8]<- 'LM_SR_Slope_SE'
  
  df<- full_join(df, Diff, by = NULL)
  return(df)
}

calc_diff_lm_AR<- function(df){
  colnames(df)[colnames(df) == 'VP_Slope'] <- 'VP'
  VP<- df[df$Sample_Type == "VP",]
  VPC<- df[df$Sample_Type == "VPC",]
  Diff<- as.data.frame(cbind(VP[,c('Location', 'Expt_No', 'Depth','Sample_Type', 'Population', 'Time_Range')], VPC[,'VP'] -VP[,'VP'], VPC[,'VP_SE'] + VP[,'VP_SE']))
  Diff[,'Sample_Type']<- "Diff"
  colnames(Diff)[7]<- 'VP'
  colnames(Diff)[8]<- 'VP_SE'
  
  df<- full_join(df, Diff, by = NULL)
  return(df)
}

calc_diff_vpcl_AR<-function(df){
  
  VP<- df[df$Sample_Type == "VP",]
  VPC<- df[df$Sample_Type == "VPC",]
  Diff<- as.data.frame(cbind(VP[,c('Location', 'Expt_No', 'Depth','Sample_Type', 'Population', 'Time_Range')], VPC[,'VP'] -VP[,'VP']))
  Diff[,'Sample_Type']<- "Diff"
  colnames(Diff)[7]<- 'VP'
  
  df<- full_join(df, Diff, by = NULL)
  return(df)
}


####5.0 LMER Model####

lmer_model<- function(df, value = count){
  
  lmer_data<- data.frame()
  #model_plots<- list()
  n<-0
  for (rep in c(1,2,3)){
    df$Replicate[df$Replicate == rep & df$Sample_Type == 'VPC'] <- rep+3
  }
  
  for (virus in unique(df$Population)){
    
    model<- lme4::lmer(data = df, count ~ Sample_Type*as.factor(Timepoint) + (1  | Replicate))
    #plot<- model_plot(model, df = df)
    emmeans<- emmeans::emmeans(model, ~ Sample_Type|as.factor(Timepoint))
    contrast<- pairs(emmeans) 
    dataf1<- data.frame(rep(virus, length(unique(df$Timepoint))),
                        rep("Diff", length(unique(df$Timepoint))),
                        summary(contrast)$Timepoint, 
                        -(summary(contrast)$estimate),
                        summary(contrast)$SE
    )
    colnames(dataf1)<- c("Population", "Sample_Type", "Timepoint", "Mean", "SE")
    st<- summary(emmeans)[1]
    tp<- summary(emmeans)[2]
    mean<- summary(emmeans)[3]
    se<- summary(emmeans)[4]
    pop<- rep(virus, length(unique(df$Timepoint)))
    dataf2<- data.frame(pop,st,tp,mean,se)
    colnames(dataf2)<- c("Population", "Sample_Type", "Timepoint", "Mean", "SE")
    
    
    lmer_data<- rbind(lmer_data, dataf2, dataf1)
    # DF<- rbind(dataf2, dataf1) %>%
    #  arrange(Timepoint)
    #n<- n +1
    #lmer_data[[n]]<- DF
    #model_plots[[n]]<- plot
    
  }
  return(lmer_data)
  return(model_plots)
  
  
}

####Bacterial End point####

bacterial_endpoint_range<- function(data){ #Gives out an index of the timepoint where we should stop the assay
  DF<- df_AVG(data)%>%
    filter(Sample_Type == 'VP', Microbe == 'Bacteria')
  DF$Timepoint<- as.numeric(DF$Timepoint)
  gt<- c()
  bp<- c()
  timepoint<- unique(DF$Timepoint)
  len<-length(timepoint)
  for (bacteria in c('c_Bacteria', 'c_HNA', 'c_LNA')){
    for (x in 2:len){ 
      Ba<- DF[DF$Population== bacteria ,] %>%
        arrange(Timepoint)
      y<- (log10(2)*(Ba$Timepoint[x]-Ba$Timepoint[1]))/(log10(Ba$mean[x])-log10(Ba$mean[1]))
      
      gt[length(gt)+1]<- y
      #if (y > 24){
      # print("Low bacterial production") }
      #if (y < 0){
      # print("Low bacterial production") }
      # plot(plot, x=c(3,6,17,20,24))
      # abline(h=24)
      # abline(h=48, col=2)
      
      
    }}
  #Interesting! LNA bacterial growth isn't that significant
  
  Bacterial_GT <- gt[1:(len-1)] #Need to find a way to save these values.
  HNA_GT <- gt[len:((len*2)-2)]
  LNA_GT <- gt[((len*2)-1):((len*3)-3)]
  bp_df<- data.frame(Bacterial_GT, HNA_GT, LNA_GT) #in hours
  
  bp_endpoint<- intersect(which(Bacterial_GT>0), which(Bacterial_GT<24))[1] #tp at which it decreases below 24
  #use this as the index for highlighting on plot
  bp<- Bacterial_GT[bp_endpoint] #Bacterial production at the endpoint
  #intersect(which(HNA_GT>0), which(HNA_GT<24))[1]
  #intersect(which(LNA_GT>0), which(LNA_GT<24))[1]
  
  if (NA %in% bp_endpoint){
    bp_endpoint_range<- paste("T", timepoint[1], "_T", timepoint[length(timepoint)], sep = "")
  } else {
    bp_endpoint_range<- paste("T", timepoint[1], "_T", timepoint[bp_endpoint], sep = "")
  }
  
  return(bp_endpoint_range)
}

####1.0 VIPCAL SR ####
#Calculates VIPCAL for every replicate.

vipcal_sr<- function(df_sr){ #takes SR dataframe as an input
  
  vipcal_vp<- list()
  vpcl_vp_df<- data.frame()
  for (location in unique(df_sr$Location)){
    for (expt_no in unique(df_sr$Expt_No)){
      for (depth in unique(df_sr$Depth)){
        for (time in unique(df_sr$Time_Range)){
          for (viruses in unique(df_sr[df_sr$Microbe == "Viruses",]$Population)){
            for (prod in unique(df_sr$Sample_Type)){
              for (rep in unique(df_sr$Replicate)){
                
                vpcl_df2<- df_sr[df_sr$Location == location,]
                vpcl_df2<- vpcl_df2[vpcl_df2$Expt_No == expt_no,]
                vpcl_df2<- vpcl_df2[vpcl_df2$Depth == depth,]
                vpcl_df2<- vpcl_df2[vpcl_df2$Time_Range == time,]
                vpcl_df2<- vpcl_df2[vpcl_df2$Population == viruses,]
                vpcl_df2<- vpcl_df2[vpcl_df2$Sample_Type == prod,]
                vpcl_df2<- vpcl_df2[vpcl_df2$Replicate == rep,]
                
                
                p<- peaks(c(+10e+10, vpcl_df2$count, -10e+10))-1
                v<- valleys(c(+10e+10, vpcl_df2$count, -10e+10))-1
                
                if(identical(length(p), length(v))){
                  print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
                }else{
                  print("Number peaks and valleys are not the same. This will lead to erroneous viral production calculations")
                }
                
                if(length(p)==0){
                  vp<- NA
                  
                } else if (length(p)==1) {
                  vp<- (vpcl_df2$count[p[1]] - vpcl_df2$count[v[1]])/(vpcl_df2$Timepoint[p[1]] - vpcl_df2$Timepoint[v[1]])
                  
                } else if (length(p)==2) {
                  vp<- ((vpcl_df2$count[p[1]] - vpcl_df2$count[v[1]])/(vpcl_df2$Timepoint[p[1]] - vpcl_df2$Timepoint[v[1]]) + 
                          (vpcl_df2$count[p[2]] - vpcl_df2$count[v[2]])/(vpcl_df2$Timepoint[p[2]] - vpcl_df2$Timepoint[v[2]]))/2
                  
                } else if (length(p)==3) {
                  vp<-  ((vpcl_df2$count[p[1]] - vpcl_df2$count[v[1]])/(vpcl_df2$Timepoint[p[1]] - vpcl_df2$Timepoint[v[1]]) + 
                           (vpcl_df2$count[p[2]] - vpcl_df2$count[v[2]])/(vpcl_df2$Timepoint[p[2]] - vpcl_df2$Timepoint[v[2]]) +
                           (vpcl_df2$count[p[3]] - vpcl_df2$count[v[3]])/(vpcl_df2$Timepoint[p[3]] - vpcl_df2$Timepoint[v[3]]))/3
                }
                
                
                vipcal<- c(location, expt_no, depth, time, viruses, prod, rep, vp)
                
                vipcal_vp[[length(vipcal_vp)+1]] <- vipcal  
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
        for (time in unique(df_avg$Time_Range)){
          for (viruses in unique(df_avg[df_avg$Microbe == "Viruses",]$Population)){
            for (prod in unique(df_avg$Sample_Type)){
              
              
              vpcl_df2<- df_avg[df_avg$Location == location,]
              vpcl_df2<- vpcl_df2[vpcl_df2$Expt_No == expt_no,]
              vpcl_df2<- vpcl_df2[vpcl_df2$Depth == depth,]
              vpcl_df2<- vpcl_df2[vpcl_df2$Time_Range == time,]
              vpcl_df2<- vpcl_df2[vpcl_df2$Population == viruses,]
              vpcl_df2<- vpcl_df2[vpcl_df2$Sample_Type == prod,]
              
              
              p<- peaks_se(c(+10e+10, vpcl_df2$mean, -10e+10),
                           c(0, vpcl_df2$se, 0))-1
              v<- valleys_se(c(+10e+10, vpcl_df2$mean, -10e+10),
                             c(0, vpcl_df2$se, 0))-1
              
              if(identical(length(p), length(v))){
                print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
              }else{
                print("Number peaks and valleys are not the same. This will lead to erroneous viral production calculations")
              }
              
              if(length(p)==0){
                vp<- NA
                
              } else if (length(p)==1) {
                vp<- (vpcl_df2$mean[p[1]] - vpcl_df2$mean[v[1]])/(vpcl_df2$Timepoint[p[1]] - vpcl_df2$Timepoint[v[1]])
                
              } else if (length(p)==2) {
                vp<- ((vpcl_df2$mean[p[1]] - vpcl_df2$mean[v[1]])/(vpcl_df2$Timepoint[p[1]] - vpcl_df2$Timepoint[v[1]]) + 
                        (vpcl_df2$mean[p[2]] - vpcl_df2$mean[v[2]])/(vpcl_df2$Timepoint[p[2]] - vpcl_df2$Timepoint[v[2]]))/2
                
              } else if (length(p)==3) {
                vp<-  ((vpcl_df2$mean[p[1]] - vpcl_df2$mean[v[1]])/(vpcl_df2$Timepoint[p[1]] - vpcl_df2$Timepoint[v[1]]) + 
                         (vpcl_df2$mean[p[2]] - vpcl_df2$mean[v[2]])/(vpcl_df2$Timepoint[p[2]] - vpcl_df2$Timepoint[v[2]]) +
                         (vpcl_df2$mean[p[3]] - vpcl_df2$mean[v[3]])/(vpcl_df2$Timepoint[p[3]] - vpcl_df2$Timepoint[v[3]]))/3
              }
              
              if(length(p)==0){
                se<- NA
                
              } else if (length(p)==1) {
                se<- (vpcl_df2$se[p[1]] + vpcl_df2$se[v[1]])/(vpcl_df2$Timepoint[p[1]] - vpcl_df2$Timepoint[v[1]])
                
              } else if (length(p)==2) {
                se<- ((vpcl_df2$se[p[1]] + vpcl_df2$se[v[1]])/(vpcl_df2$Timepoint[p[1]] - vpcl_df2$Timepoint[v[1]]) + 
                        (vpcl_df2$se[p[2]] + vpcl_df2$se[v[2]])/(vpcl_df2$Timepoint[p[2]] - vpcl_df2$Timepoint[v[2]]))/2
                
              } else if (length(p)==3) {
                se<-  ((vpcl_df2$se[p[1]] + vpcl_df2$se[v[1]])/(vpcl_df2$Timepoint[p[1]] - vpcl_df2$Timepoint[v[1]]) + 
                         (vpcl_df2$se[p[2]] + vpcl_df2$se[v[2]])/(vpcl_df2$Timepoint[p[2]] - vpcl_df2$Timepoint[v[2]]) +
                         (vpcl_df2$se[p[3]] + vpcl_df2$se[v[3]])/(vpcl_df2$Timepoint[p[3]] - vpcl_df2$Timepoint[v[3]]))/3
              }
              
              
              vipcal<- c(location, expt_no, depth, time, viruses, prod, vp,se)
              
              vipcal_vp[[length(vipcal_vp)+1]] <- vipcal  
              
            }
          }
        }
      }
    }
  }
  
  vpcl_vp_df<- data.frame(t(sapply(vipcal_vp, c)))
  colnames(vpcl_vp_df)<- c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population', 'Sample_Type',  'VP_Mean', 'VP_SE')
  vpcl_vp_df[, c('VP_Mean', 'VP_SE')]<- lapply(vpcl_vp_df[, c('VP_Mean', 'VP_SE')], as.numeric)
  
  
  return(vpcl_vp_df)
}

vipcal_avg<- function(df_avg){ #takes avg_tp dataframe as an input
  
  vipcal_vp<- list()
  vpcl_vp_df<- data.frame()
  for (location in unique(df_avg$Location)){
    for (expt_no in unique(df_avg$Expt_No)){
      for (depth in unique(df_avg$Depth)){
        for (time in unique(df_avg$Time_Range)){
          for (viruses in unique(df_avg[df_avg$Microbe == "Viruses",]$Population)){
            for (prod in unique(df_avg$Sample_Type)){
              
              
              vpcl_df2<- df_avg[df_avg$Location == location,]
              vpcl_df2<- vpcl_df2[vpcl_df2$Expt_No == expt_no,]
              vpcl_df2<- vpcl_df2[vpcl_df2$Depth == depth,]
              vpcl_df2<- vpcl_df2[vpcl_df2$Time_Range == time,]
              vpcl_df2<- vpcl_df2[vpcl_df2$Population == viruses,]
              vpcl_df2<- vpcl_df2[vpcl_df2$Sample_Type == prod,]
              
              
              p<- peaks(c(+10e+10, vpcl_df2$mean, -10e+10))-1
              v<- valleys(c(+10e+10, vpcl_df2$mean, -10e+10))-1
              
              if(identical(length(p), length(v))){
                print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
              }else{
                print("Number peaks and valleys are not the same. This will lead to erroneous viral production calculations")
              }
              
              if(length(p)==0){
                vp<- NA
                
              } else if (length(p)==1) {
                vp<- (vpcl_df2$mean[p[1]] - vpcl_df2$mean[v[1]])/(vpcl_df2$Timepoint[p[1]] - vpcl_df2$Timepoint[v[1]])
                
              } else if (length(p)==2) {
                vp<- ((vpcl_df2$mean[p[1]] - vpcl_df2$mean[v[1]])/(vpcl_df2$Timepoint[p[1]] - vpcl_df2$Timepoint[v[1]]) + 
                        (vpcl_df2$mean[p[2]] - vpcl_df2$mean[v[2]])/(vpcl_df2$Timepoint[p[2]] - vpcl_df2$Timepoint[v[2]]))/2
                
              } else if (length(p)==3) {
                vp<-  ((vpcl_df2$mean[p[1]] - vpcl_df2$mean[v[1]])/(vpcl_df2$Timepoint[p[1]] - vpcl_df2$Timepoint[v[1]]) + 
                         (vpcl_df2$mean[p[2]] - vpcl_df2$mean[v[2]])/(vpcl_df2$Timepoint[p[2]] - vpcl_df2$Timepoint[v[2]]) +
                         (vpcl_df2$mean[p[3]] - vpcl_df2$mean[v[3]])/(vpcl_df2$Timepoint[p[3]] - vpcl_df2$Timepoint[v[3]]))/3
              }
              
              
              vipcal<- c(location, expt_no, depth, time, viruses, prod, vp)
              
              vipcal_vp[[length(vipcal_vp)+1]] <- vipcal  
              
            }
          }
        }
      }
    }
  }
  
  vpcl_vp_df<- data.frame(t(sapply(vipcal_vp, c)))
  colnames(vpcl_vp_df)<- c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population', 'Sample_Type',  'VP_Mean')
  vpcl_vp_df$VP_Mean<- as.numeric(vpcl_vp_df$VP_Mean)
  return(vpcl_vp_df)
}

vipcal_sr_diff_SE<- function(df_sr){ #takes SR dataframe as an input
  
  vipcal_vp<- list()
  vpcl_vp_df<- data.frame()
  for (location in unique(df_sr$Location)){
    for (expt_no in unique(df_sr$Expt_No)){
      for (depth in unique(df_sr$Depth)){
        for (time in unique(df_sr$Time_Range)){
          for (viruses in unique(df_sr[df_sr$Microbe == "Viruses",]$Population)){
            
            
            vpcl_df2<- df_sr[df_sr$Location == location,]
            vpcl_df2<- vpcl_df2[vpcl_df2$Expt_No == expt_no,]
            vpcl_df2<- vpcl_df2[vpcl_df2$Depth == depth,]
            vpcl_df2<- vpcl_df2[vpcl_df2$Time_Range == time,]
            vpcl_df2<- vpcl_df2[vpcl_df2$Population == viruses,]
            
            try(vpcl_df3<- lmer_model(vpcl_df2))#the output is a dataframe for VP, VPC and Diff together
            
            if (exists("vpcl_df3")){
              for (prod in unique(vpcl_df3$Sample_Type)){
                vpcl_df4<-vpcl_df3[vpcl_df3$Sample_Type == prod,]
                
                p<- peaks_sd(c(+10e+10, vpcl_df4$Mean, -10e+10),
                             c(0, vpcl_df4$SE, 0))-1
                v<- valleys_sd(c(+10e+10, vpcl_df4$Mean, -10e+10),
                               c(0, vpcl_df4$SE, 0))-1
                
                if(identical(length(p), length(v))){
                  print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
                }else{
                  print("Number peaks and valleys are not the same. This will lead to erroneous viral production calculations")
                }
                print(paste(location, expt_no, depth, time, viruses, prod))
                print(p)
                print(v)
                if(length(p)==0){
                  vp<- NA
                  
                } else if (length(p)==1) {
                  vp<- (vpcl_df4$Mean[p[1]] - vpcl_df4$Mean[v[1]])/(vpcl_df4$Timepoint[p[1]] - vpcl_df4$Timepoint[v[1]])
                  
                } else if (length(p)==2) {
                  vp<- ((vpcl_df4$Mean[p[1]] - vpcl_df4$Mean[v[1]])/(vpcl_df4$Timepoint[p[1]] - vpcl_df4$Timepoint[v[1]]) + 
                          (vpcl_df4$Mean[p[2]] - vpcl_df4$Mean[v[2]])/(vpcl_df4$Timepoint[p[2]] - vpcl_df4$Timepoint[v[2]]))/2
                  
                } else if (length(p)==3) {
                  vp<-  ((vpcl_df4$Mean[p[1]] - vpcl_df4$Mean[v[1]])/(vpcl_df4$Timepoint[p[1]] - vpcl_df4$Timepoint[v[1]]) + 
                           (vpcl_df4$Mean[p[2]] - vpcl_df4$Mean[v[2]])/(vpcl_df4$Timepoint[p[2]] - vpcl_df4$Timepoint[v[2]]) +
                           (vpcl_df4$Mean[p[3]] - vpcl_df4$Mean[v[3]])/(vpcl_df4$Timepoint[p[3]] - vpcl_df4$Timepoint[v[3]]))/3
                }
                
                if(length(p)==0){
                  se<- NA
                  
                } else if (length(p)==1) {
                  se<- (vpcl_df4$SE[p[1]] + vpcl_df4$SE[v[1]])/(vpcl_df4$Timepoint[p[1]] - vpcl_df4$Timepoint[v[1]])
                  
                } else if (length(p)==2) {
                  se<- ((vpcl_df4$SE[p[1]] + vpcl_df4$SE[v[1]])/(vpcl_df4$Timepoint[p[1]] - vpcl_df4$Timepoint[v[1]]) + 
                          (vpcl_df4$SE[p[2]] + vpcl_df4$SE[v[2]])/(vpcl_df4$Timepoint[p[2]] - vpcl_df4$Timepoint[v[2]]))/2
                  
                } else if (length(p)==3) {
                  se<-  ((vpcl_df4$SE[p[1]] + vpcl_df4$SE[v[1]])/(vpcl_df4$Timepoint[p[1]] - vpcl_df4$Timepoint[v[1]]) + 
                           (vpcl_df4$SE[p[2]] + vpcl_df4$SE[v[2]])/(vpcl_df4$Timepoint[p[2]] - vpcl_df4$Timepoint[v[2]]) +
                           (vpcl_df4$SE[p[3]] + vpcl_df4$SE[v[3]])/(vpcl_df4$Timepoint[p[3]] - vpcl_df4$Timepoint[v[3]]))/3
                }
                
                
                vipcal<- c(location, expt_no, depth, time, viruses, prod, vp,se)
                
                vipcal_vp[[length(vipcal_vp)+1]] <- vipcal   
                
                rm(p,v)
                
                
              }
            } else{
              for (prod in unique(vpcl_df2$Sample_Type)){
                vipcal<- c(location, expt_no, depth, time, viruses, prod, NA, NA)
                vipcal_vp[[length(vipcal_vp)+1]] <- vipcal 
              }
            }
            
            
            
          }     
        }
      }
    }
  }
  
  vpcl_vp_df<- data.frame(t(sapply(vipcal_vp, c)))
  colnames(vpcl_vp_df)<- c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population', 'Sample_Type',  'VP_Mean', 'VP_SE')
  vpcl_vp_df[, c('VP_Mean', 'VP_SE')]<- lapply(vpcl_vp_df[, c('VP_Mean', 'VP_SE')], as.numeric)
  
  
  return(vpcl_vp_df)
  rm(vpcl_df3)
  rm(vpcl_df4)
}

vipcal_sr_diff_no_SE<- function(df_sr){ #takes SR dataframe as an input
  
  vipcal_vp<- list()
  vpcl_vp_df<- data.frame()
  for (location in unique(df_sr$Location)){
    for (expt_no in unique(df_sr$Expt_No)){
      for (depth in unique(df_sr$Depth)){
        for (time in unique(df_sr$Time_Range)){
          for (viruses in unique(df_sr[df_sr$Microbe == "Viruses",]$Population)){
            
            
            vpcl_df2<- df_sr[df_sr$Location == location,]
            vpcl_df2<- vpcl_df2[vpcl_df2$Expt_No == expt_no,]
            vpcl_df2<- vpcl_df2[vpcl_df2$Depth == depth,]
            vpcl_df2<- vpcl_df2[vpcl_df2$Time_Range == time,]
            vpcl_df2<- vpcl_df2[vpcl_df2$Population == viruses,]
            
            try(vpcl_df3<- lmer_model(vpcl_df2))#the output is a dataframe for VP, VPC and Diff together
            
            if (exists("vpcl_df3")){
              for (prod in unique(vpcl_df3$Sample_Type)){
                vpcl_df4<-vpcl_df3[vpcl_df3$Sample_Type == prod,]
                
                p<- peaks(c(+10e+10, vpcl_df4$Mean, -10e+10))-1
                v<- valleys(c(+10e+10, vpcl_df4$Mean, -10e+10))-1
                
                if(identical(length(p), length(v))){
                  print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
                }else{
                  print("Number peaks and valleys are not the same. This will lead to erroneous viral production calculations")
                }
                
                if(length(p)==0){
                  vp<- NA
                  
                } else if (length(p)==1) {
                  vp<- (vpcl_df4$Mean[p[1]] - vpcl_df4$Mean[v[1]])/(vpcl_df4$Timepoint[p[1]] - vpcl_df4$Timepoint[v[1]])
                  
                } else if (length(p)==2) {
                  vp<- ((vpcl_df4$Mean[p[1]] - vpcl_df4$Mean[v[1]])/(vpcl_df4$Timepoint[p[1]] - vpcl_df4$Timepoint[v[1]]) + 
                          (vpcl_df4$Mean[p[2]] - vpcl_df4$Mean[v[2]])/(vpcl_df4$Timepoint[p[2]] - vpcl_df4$Timepoint[v[2]]))/2
                  
                } else if (length(p)==3) {
                  vp<-  ((vpcl_df4$Mean[p[1]] - vpcl_df4$Mean[v[1]])/(vpcl_df4$Timepoint[p[1]] - vpcl_df4$Timepoint[v[1]]) + 
                           (vpcl_df4$Mean[p[2]] - vpcl_df4$Mean[v[2]])/(vpcl_df4$Timepoint[p[2]] - vpcl_df4$Timepoint[v[2]]) +
                           (vpcl_df4$Mean[p[3]] - vpcl_df4$Mean[v[3]])/(vpcl_df4$Timepoint[p[3]] - vpcl_df4$Timepoint[v[3]]))/3
                }
                
                
                vipcal<- c(location, expt_no, depth, time, viruses, prod, vp)
                
                vipcal_vp[[length(vipcal_vp)+1]] <- vipcal  
                
                
                
                
              }
            } else{
              for (prod in unique(vpcl_df2$Sample_Type)){
                vipcal<- c(location, expt_no, depth, time, viruses, prod,  NA)
                vipcal_vp[[length(vipcal_vp)+1]] <- vipcal 
              }
            }
            
            rm(p,v)
            
          }     
        }
      }
    }
  }
  
  vpcl_vp_df<- data.frame(t(sapply(vipcal_vp, c)))
  colnames(vpcl_vp_df)<- c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'VP_Mean')
  vpcl_vp_df$VP_Mean <- as.numeric(vpcl_vp_df$VP_Mean) 
  return(vpcl_vp_df)
  rm(vpcl_df3, vpcl_df4)
  
}

#A script to install and load required packages, and to source global functions


####1. Installing and loading libraries####

#Installing BiocManager
{
  if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
}

#List of packages to install
packages_to_load<- c("tidyverse", 
                     "flowWorkspace",
                     "scales",
                     "readxl",
                     "emmeans",
                     "lme4")


#Checking if packages are already present. If absent, then installing packages from BiocManager 
for (pack in packages_to_load){
  if(!requireNamespace(pack))
    BiocManager::install(pack, force = T)
  }

#Loading libraries  
for (pack in packages_to_load){
   library(pack, character.only = T)
}





####2. Overview plot functions ####


####2.1 To create an overview table with time ranges. Needed for overview plotting####
overview_df_tp_avg<- function(df, keep_0.22 = F) {
  #only works if we remove 0.22
  DF<- df[df$Sample_Type != '0.22',]
  
  DF<- DF %>%
    gather('c_Bacteria', 'c_HNA', 'c_LNA', 'c_Viruses', 'c_V1', 'c_V2', 'c_V3', key="count", value="value") %>%
    #gathering all the count information in two columns
    group_by(Location, Expt_No, Depth, Timepoint, Sample_Type, count ) %>% #grouping the dataframe in this order
    summarise(n =n(), mean=mean(value), sd=sd(value)) #calculating means and sd
  
  if ('VPC' %in% DF$Sample_Type){ #Here, if we have VPC samples, it'll create necessary columns for it
    colnames_mean<- c("VP", "VPC")
    colnames_sd<- c("VP", "VPC")
    
  } else {
    colnames_mean<- c("VP")
    colnames_sd<- c("VP")
    
  }
  
  DF_mean<- select(DF, -c('sd')) %>% #splitting the dataframe cause I haven't figure out how to spread the table without adding NAs
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
  
  if ('VPC' %in% DF$Sample_Type){ #If VPC exists, then we need to calculate Diff curve. Here we use simple subtarction between mean VPC and mean VP, and add their SDs
    DF_mean$Diff <- with(DF_mean, VPC-VP) #calculating Diff mean
    DF_mean<- pivot_longer(DF_mean, cols = c("VP", "VPC", "Diff"), names_to= 'Sample_Type', values_to='mean_value')
    DF_sd$Diff <- with(DF_sd, VPC+VP) #Calculating Diff sd, which is addition of the other sds
    DF_sd<- pivot_longer(DF_sd, cols = c("VP", "VPC", "Diff"), names_to='Sample_Type', values_to= 'sd_value')
  }
  
  #Merging both mean and sd dataframes
  DF<- merge(DF_mean, DF_sd, by= c('Location', 'Expt_No', 'Depth', 
                                   'Timepoint', 'count', 'n', 'Sample_Type'))
  
  DF<-  mutate(DF, Subgroup = case_when(DF$count == 'c_Bacteria' ~ "Total",
                                        DF$count == 'c_Viruses' ~ "Total",
                                        DF$count == 'c_HNA' ~ "Bacteria",
                                        DF$count == 'c_LNA' ~ "Bacteria",
                                        DF$count == 'c_V1' ~ "Viruses",
                                        DF$count == 'c_V2' ~ "Viruses",
                                        DF$count == 'c_V3' ~ "Viruses")) %>%
    arrange(Timepoint)
  
  rm('DF_mean', 'DF_sd', 'colnames_mean', 'colnames_sd')
  
  
  
  
  #For the time ranges
  TP<- unique(DF$Timepoint) #identifying the timepoints present
  colnames<- c()
  for (col in 2: length(TP)){
    a<- paste("T", TP[1], "_T", TP[col], sep = "")
    colnames[length(colnames)+1]<- a
  }
  colvalues<- c()
  for (col in 2: length(TP)){
    a<- paste("T", TP[1], ":T", TP[col], sep = "")
    colvalues[length(colvalues)+1]<- a
  }
  
  ncol<- ncol(DF)
  DF[colnames]<- NA
  
  #This works for 6 timepoints
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


####2.2 To create an overview plot with time ranges####
overview_plots_tp_avg<- function(df, ...){
  n<- ggplot(df, aes(x= Timepoint, y= mean, color= count , shape=count))+
    geom_point(size= 1.5)+
    geom_smooth(size= 1.0, method = 'lm', se = F)+
    geom_line()+
    geom_hline(yintercept = 0, color= '#636363', size= 0.3, linetype = "dashed")+
    facet_grid(factor(Subgroup, levels = c("Total", "Bacteria", "Viruses")) +
                 factor(Sample_Type, levels = c("VP", "VPC", "Diff"))~
                 factor(Time_Time, levels = unique(df$Time_Time))  
    )+
    scale_color_manual(name= 'Populations',
                       labels=c("Total Bacteria", "HNA Bacteria", "LNA Bacteria",
                                "Total Viruses", "V1 Viruses", "V2 Viruses", "V3 Viruses"),
                       values= c(c_Bacteria = "#AD002A99", 
                                 c_HNA = "#00468B99",
                                 c_LNA = "#0099B499",
                                 c_Viruses = "#ED000099",
                                 c_V1 = "#1B191999",
                                 c_V2 = "#7E6148B2",
                                 c_V3 = "#9C9C9C"))+
    scale_shape_manual(name = 'Populations', 
                       values = c(16,16,16,17,17,17,17),
                       labels=c("Total Bacteria", "HNA Bacteria", "LNA Bacteria",
                                "Total Viruses", "V1 Viruses", "V2 Viruses", "V3 Viruses"))+
    scale_y_continuous(labels = unit_format(unit = NULL, scale = 1e-6),
                       breaks= seq(-4e+6,13e+6, 4e+6),
                       limits = c(-4e+6, 14e+6))+
    theme_bw()+
    geom_errorbar(aes(ymin=mean + se, ymax= mean - se), width = 0.5, size = 0.5)+
    scale_x_continuous(breaks = unique(df$Timepoint))+
      labs(title = paste(paste(unique(df$Location), "St", unique(df$Expt_No), "Depth", unique(df$Depth), sep = "_"), "- Viral Production Assay - Timepoint Sloughing"), subtitle = 'Overview - Bacterial and Viral counts for Lytic and Lysogenic inductions',
         x= 'Sampling Timepoints\n (in hours)', y='FCM Counts\n (in millions)')+
    theme(strip.text = element_text(face = "bold"),
          strip.background = element_rect( color = 'black', fill = 'white'),
          axis.title = element_text(face = 'bold'),
          title = element_text(face = 'bold'),
          legend.position = "bottom")+
    guides(color = guide_legend(nrow = 2, byrow = TRUE))

  return(n)
  
}


####2.3 For bacterial generation time end point####
bacterial_endpoint<- function (df){ 
  
  #creating the dataframe required for calculating bacterial generation time
  DF<- df[df$Sample_Type != '0.22',] #Lose 0.22 values
  
  DF<- gather(DF, 'c_Bacteria', 'c_HNA', 'c_LNA', 'c_Viruses', 'c_V1', 'c_V2', 'c_V3', key="count", value="value") %>%
    group_by(Location, Expt_No, Depth, Sample_Type, Timepoint, count ) %>%
    summarise(n =n(), mean=mean(value), sd=sd(value)) #calculating means and sd
  
  DF_mean<- DF[,1:8] %>% #splitting the datfram cause I haven't figure out how to spread teh table without adding NAs
    spread('Sample_Type', 'mean')
  colnames(DF_mean)[7:8]<- c("VP_mean", "VPC_mean")
  DF_mean$Diff_mean <- with(DF_mean, VPC_mean-VP_mean) #calculating Diff mean
  
  DF_sd<- DF[,c(1:7,9)] %>%
    spread('Sample_Type', 'sd')
  colnames(DF_sd)[7:8]<- c("VP_sd", "VPC_sd")
  DF_sd$Diff_sd <- with(DF_sd, VPC_sd+VP_sd) #Calculating Diff sd, whcih si addition of the other sds
  
  
  DF<- merge(DF_mean, DF_sd, by= c('Location', 'Expt_No', 'Depth',
                                      'Timepoint', 'count', 'n')) %>%
    mutate(Microbe = if_else(count == 'c_Bacteria' | count == 'c_HNA' | count == 'c_LNA', "Bacteria", "Viruses"))%>%
    mutate(Subgroup = if_else(count == 'c_Bacteria' | count == 'c_Viruses', "Parent", "Subgroup"))
  
  rm(DF_mean)
  rm(DF_sd)
  
#Calculating Bacterial Generation time
  gt<- c()
  bp<- c()
  
  for (bacteria in c('c_Bacteria', 'c_HNA', 'c_LNA')){
    for (x in 2:6){ #2:6 should be adjusted for more timepoints. Here it is for 6 tps.
      Ba<- DF[DF$count== bacteria,] %>%
        arrange(Timepoint)
      y<- (log10(2)*(Ba$Timepoint[x]-Ba$Timepoint[1]))/(log10(Ba$VP_mean[x])-log10(Ba$VP_mean[1]))
      
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
  
  Bacterial_GT <- gt[1:5] #Need to find a way to save these values.
  HNA_GT <- gt[6:10]
  LNA_GT <- gt[11:15]
  bp_df<- data.frame(Bacterial_GT, HNA_GT, LNA_GT) #in hours
  
  bp_endpoint<- intersect(which(Bacterial_GT>0), which(Bacterial_GT<24))[1] - 1 #use this as the index for highlighting on plot
  bp<- Bacterial_GT[bp_endpoint] #Bacterial production at the endpoint
  #intersect(which(HNA_GT>0), which(HNA_GT<24))[1]
  #intersect(which(LNA_GT>0), which(LNA_GT<24))[1]
  

  return(bp_endpoint)
  
}


####2.4 Highlighting bacterial end point on the overview plot####
bp_endpoint_highlight <- function(plot) {
  n<- ggplot_gtable(ggplot_build(plot))
  strip_both<- which(grepl('strip-', n$layout$name))
  fills<- c(rep("white", 14))
  fills[bp]<- "#E489A6" #Change the index to that of the time range after bacterial generation time
  k <- 1
  
  for (i in strip_both) {
    j<- which(grepl('rect', n$grobs[[i]]$grobs[[1]]$childrenOrder))
    n$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k<- k+1
  }
  #https://ojkftyk.blogspot.com/2019/01/r-ggplot2-change-colour-of-font-and.html
  n<- grid::grid.draw(n)
  return (n)
}




####2.5 Function that combines the functions from #2 ####
overview_vp_plots <- function(x){
  
  bp<- bacterial_endpoint(x)
  
  skeleton_plot<- overview_df_tp_avg(x)%>%
    overview_plots_tp_avg()
  
    skeleton_plot<- bp_endpoint_highlight(skeleton_plot)
    
  
  print(skeleton_plot)
}

 



####3.0 LMER plots function####
#writing a model plot function to make things easier.
model_plot<- function(model, df){ #lmer input for model, and a data frame for observed values
  a<- model
  b<- 
    par(mfrow= c(2,3))
  #1. Observed data Histogram
  hist(df$value,
       xlab = "Observed Values",
       main = "Observed Values"
  )
  
  #2. Predicted Data Histogram
  hist(predict(a),
       xlab = "Predicted Values",
       main = "Predicted Values"
  ) 
  
  #3. Predicted vs Observed Values Line plot
  plot(predict(a)~ df$value, 
       ylim = c(0, 1.5e+07), 
       xlim = c(0, 1.5e+07),
       xlab = "Observed values", 
       ylab = "Predicted Values", 
       main = "Predicted vs Observed Values")
  abline(0,1, lwd=2, col = "coral3")
  
  fit2<- lm(predict(model) ~ df$value)
  lgd <- c(
    paste("R^2 =", round(summary(fit2)$r.squared,3)),
    paste("Intercept =", format(coef(fit2)[1], scientific = T, digits = 4)),
    paste("Slope =", round(coef(fit2)[2],3))
  )
  legend("topleft", legend=lgd)
  abline(fit2, lwd=2)
  legend("bottomright", legend=c("predicted ~ observed", "1:1"), col=c(1,"coral3"), lty=1, lwd=2)
  rm(fit2, lgd)
  
  
  #4. Normalized Residuals vs Predicted Values
  plot(resid(a) ~ fitted(a), 
       xlab = "Predicted values", 
       ylab = "Normalized residuals", 
       main = "Normalized Residuals vs Predicted Values"
  )
  abline(h = 0, lty = 2, col = "coral3", lwd = 2)
  
  #5. Normalized Residuals vs Replicates
  boxplot(resid(a) ~ Replicate, 
          data = df, 
          xlab = "Replicate",
          ylab = "Normalized residuals",
          main = "Normalized Residuals vs Replicates"
  )
  abline(h = 0, lty = 2, col = "coral3", lwd = 2)
  
  #6. QQ-plot
  qqnorm(resid(a), pch = 16)
  qqline(resid(a), col = "coral3", lwd = 2)
  
}


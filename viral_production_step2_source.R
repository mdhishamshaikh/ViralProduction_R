### Source file for viral_production_step2.R ###

## Setup
# Installing BiocManager if not presented
{
  if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
}

# List of packages
packages_to_load<- c("tidyverse", "flowWorkspace", "scales",
                     "readxl", "emmeans", "lme4", "ggsci",
                     "svglite", "tidyr")

# If package not presented, install with BiocManager 
for(pack in packages_to_load){
  if(!requireNamespace(pack))
    BiocManager::install(pack, force = T, update = F)
}

# Load packages 
for (pack in packages_to_load){
  library(pack, character.only = T)
}

# A collection of functions that calculate viral production using different methods
source("vp_calc_functions.R")

## Functions needed for different calculation methods of viral production
# 1. Dataframe functions: Each method requires a specific type of dataframe as input
# An average of the replicates is taken or each replicate is evaluated separately 
# First of all, the different time points need to be defined
tp <- function(DF){
  
  TP <- unique(as.numeric(DF$Timepoint)) # Unique timepoints of assay
  
  # Determining colnames and colvalues, always start from T0 so iterate from 2nd element
  colnames<- c()
  
  for(col in 2: length(TP)){
    a <- paste("T", TP[1], "_T", TP[col], sep = "")
    colnames[length(colnames)+1] <- a
  }
  
  colvalues<- c()
  for(col in 2: length(TP)){
    a <- paste("T", TP[1], ":T", TP[col], sep = "")
    colvalues[length(colvalues)+1] <- a
  }
  
  # Adding columns to link colnames to colvalues
  ncol<- ncol(DF)
  DF[colnames]<- NA # Setting all to NA
  
  # 1st Case: Timepoint is 0 or 3
  DF[,ncol+1]<- case_when(DF$Timepoint == TP[1] ~ colvalues[1],
                          DF$Timepoint == TP[2] ~ colvalues[1])
  
  # 2nd Case: Timepoint is 0, 3 or 6
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
  
  DF <- DF %>%
    pivot_longer(cols = colnames, names_to = "Time_Range", values_to = "Time_Time")
  return(DF)
}
  
  





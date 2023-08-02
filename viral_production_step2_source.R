###---Source file for viral_production_step2.R---###
# General setup + collection of functions needed in methods of vp_calc_functions

### Setup
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

### Functions needed for different calculation methods of viral production (vp_calc_functions.R)

## 1. Dataframe functions: 

# Each method requires a specific type of dataframe as input: an average of the replicates is taken or each replicate is evaluated separately 

# First of all, the different time points need to be defined
tp <- function(DF){
  
  TP <- unique(as.numeric(DF$Timepoint)) # Unique timepoints of assay
  
  # Timepoint always starts from T = 0 => iteration starts from 2nd timepoint
  colnames<- c() # Of the form: "TX_TY"
  
  for(col in 2: length(TP)){
    a <- paste("T", TP[1], "_T", TP[col], sep = "") # Make current timepoint
    colnames[length(colnames)+1] <- a # Add it as column to vector colnames
  }
  
  colvalues<- c() # Of the form: "TX:TY"
  for(col in 2: length(TP)){
    a <- paste("T", TP[1], ":T", TP[col], sep = "") # Make current timepoint
    colvalues[length(colvalues)+1] <- a # Add it as column to vector colvalues
  }
  
  # Adding columns to link colnames to colvalues
  ncol <- ncol(DF)
  DF[colnames] <- NA # Initiate new columns
  
  # Simplified loop for adding the columns to dataframe
  for (i in 1:(length(TP)-1)) { 
    conditions <- DF$Timepoint %in% TP[1:(i+1)]
    DF[conditions, (ncol + i)] <- colvalues[i]
  }
  
  # We have added columns to each replicate with its timeframe
  # Change the format => increase number of rows and decrease the number of columns with pivot_longer()
  DF <- DF %>%
    pivot_longer(cols = colnames, names_to = "Time_Range", values_to = "Time_Time") %>%
    drop_na() # Drops all the NA's
  # When running only this function as DF <- tp(data), an empty dataframe will be returned because of the column "Comments", which consits of all NA's
  
  return(DF)
}

# Separate replicate dataframe
df_SR <- function(df, keep_0.22 = F){
  DF <- df %>%
    select(4:16) %>% # In the understanding that the input dataframe always has the same format
    # When looking up gather(), it suggested to switch to pivot_longer() function but when doing this, the arrange line doesn't work anymore
    # I left gather() in here but if the order of the rows is not that important => you could switch gather() to pivot_longer()
    #pivot_longer(cols = 7:13, names_to = 'Population', values_to = 'Count') %>%
    gather(7:13, key = 'Population', value = 'Count') %>% # Taking index of columns instead of name, again beginning format needs to be the same every time
    mutate(Microbe = if_else(Population %in% c('c_Bacteria', 'c_HNA', 'c_LNA'), 'Bacteria', 'Viruses')) %>% # Adding an extra column defining if replicate is from bacterial or viral origin
    arrange('Location', 'Expt_No', 'Depth', 'Sample_Type','Replicate','Population', as.numeric(Timepoint)) # Reorder the rows by the values of the selected columns
  
  # By default, only go further with VP and VPC replicates => filtering out 0.22 samples which are control samples
  if (keep_0.22 == F){
    DF <- DF[DF$Sample_Type != '0.22',]
  }
  
  # Adding timepoints
  DF <- tp(DF)
  
  return(DF)
}

# Average of replicates on same timepoint dataframe
# Since for the control samples (0.22) we only have timepoints 0 and 24, we need to filter these out
df_AVG <- function(df, keep_0.22 = F){
  # Removing control samples
  DF <- df[df$Sample_Type != '0.22',]
  
  # Calculating number, mean of replicates and standard error 
  DF <- DF %>%
    select(4:16) %>%
    gather(7:13, key = 'Population', value = 'Count') %>%
    group_by(Location, Expt_No, Depth, Sample_Type, Timepoint, Population) %>%
    summarise(n = n(), mean = mean(Count), se = plotrix::std.error(Count)) # plotrix::std.error calculates the standard error of the mean
    
  # Calculate the difference between VPC samples and VP samples since for some methods we will estimate the difference curve
  # We will consider the mean and se separately, so we will use two intermediate subdataframes
  DF_mean <- DF %>%
    select(-'se') %>%
    spread('Sample_Type', 'mean')
  
  DF_se <- DF %>%
    select(-'mean') %>%
    spread('Sample_Type', 'se')
  
  if ('VPC' %in% DF$Sample_Type){
    DF_mean$Diff <- with(DF_mean, VPC - VP)
    DF_mean <- pivot_longer(DF_mean, cols = c('VP', 'VPC', 'Diff'), names_to = 'Sample_Type', values_to = 'mean')
    DF_se$Diff <- with(DF_se, VPC + VP)
    DF_se <- pivot_longer(DF_se, cols = c('VP', 'VPC', 'Diff'), names_to = 'Sample_Type', values_to = 'se')
  }
  
  # Merging in one dataframe and adding columns Microbe and Subgroup
  DF <- merge(DF_mean, DF_se, by = c('Location', 'Expt_No', 'Depth',
                                     'Timepoint', 'Population', 'n', 'Sample_Type')) %>%
    mutate(Microbe = if_else(Population %in% c('c_Bacteria', 'c_HNA', 'c_LNA'), 'Bacteria', 'Viruses')) %>%
    mutate(Subgroup = if_else(Population %in% c('c_Bacteria', 'c_Viruses'), 'Parent', 'Subgroup')) %>%
    arrange('Location', 'Expt_No', 'Depth', 'Sample_Type', 'Population', as.numeric(Timepoint))
  
  # Adding timepoints
  DF <- tp(DF)
  
  return(DF)
}

## 2. Peaks and valleys
# Two main methods will be used for "Viral Production Analyses": Linear Regression vs VIPCAL
# Linear regression calculates the viral production using a slope
# On the other hand, VIPCAL uses average of increments. Therefore, peaks and valleys need to be determined based of the counts

# Determining peaks
peaks <- function(values){
  list <- c()

  for (i in 1:(length(values)-1)){
    d <- sign(values[i+1] - values[i]) 
    list[[length(list)+1]] <- d # Nested list with 3 possible values: -1,0,1
  }
  return(which(diff(as.numeric(list)) < 0) +1)
}
  
# Determining valleys
valleys <- function(values){
  list <- c()
  
  for (i in 1:(length(values)-1)){
    d <- sign(values[i+1] - values[i]) 
    list[[length(list)+1]] <- d # Nested list with 3 possible values: -1,0,1
  }
  return(which(diff(as.numeric(list)) > 0) +1)
}

# Determining standard error on both peaks and values
peaks_se <- function(values){
  list <- c()
  
  for (i in 1:(length(values)-1)){
    d <- sign((values[i+1] - se[i+1]) - (values[i] + se[i])) 
    list[[length(list)+1]] <- d # Nested list with 3 possible values: -1,0,1
  }
  return(which(diff(as.numeric(list)) < 0) +1)
}

valleys <- function(values){
  list <- c()
  
  for (i in 1:(length(values)-1)){
    d <- sign((values[i+1] - se[i+1]) - (values[i] + se[i])) 
    list[[length(list)+1]] <- d # Nested list with 3 possible values: -1,0,1
  }
  return(which(diff(as.numeric(list)) > 0) +1)
}

## 3. Slope Functions
# For Linear Regression






# For VIPCAL








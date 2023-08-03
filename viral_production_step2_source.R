###---Source file for viral_production_step2.R---###
# General setup + collection of functions needed in methods of vp_calc_functions

### Information
# Viral production can be divided in two phases: Lytic and Lysogenic
# VP replicates represent the lytic viral production, this can be calculated with linear regression (slope) or VIPCAL (average of increments)
# Analoge for VPC samples => they represent both lytic and lysogenic viral production
# To determine only the lysogenic viral production: in linear regression VPC slope - VP slope, in VIPCAL average of increments of the difference curve
# VIPCAL has his own issues => standard error has big influence => VIPCAL-SE takes this into account and only looks at  those increments that don't have overlapping range of SE
# If there is overlap we give a 0 => rather a zero than an uncertain number => VIPCAL-SE is more conservative

# VIPCAL = online tool for estimating lytically and lysogenicaaly produced viruses
# Lytic VP as the slope between two peaks that occur in viral abundance during incubation, following function holds: VP = [(Vmax1 - Vmin1) + (Vmax2 - Vmin2)] / (Tmax2 - Tmin1) => average over all
# Lysogenic VP computed from difference curve


### Setup
# Installing BiocManager if not presented
{
  if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
}

# List of packages
packages_to_load<- c("tidyverse", "flowWorkspace", "scales",
                     "readxl", "emmeans", "lme4", "ggsci",
                     "svglite", "tidyr", "data.table")

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
df_SR <- function(data, keep_0.22 = F){
  DF <- data %>%
    select(4:16) %>% # In the understanding that the input dataframe always has the same format
    # When looking up gather(), it suggested to switch to pivot_longer() function but when doing this, the arrange line doesn't work anymore
    # I left gather() in here but if the order of the rows is not that important => you could switch gather() to pivot_longer()
    #pivot_longer(cols = 7:13, names_to = 'Population', values_to = 'Count') %>%
    gather(7:13, key = 'Population', value = 'Count') %>% # Taking index of columns instead of name, again beginning format needs to be the same every time
    mutate(Microbe = if_else(Population %in% c('c_Bacteria', 'c_HNA', 'c_LNA'), 'Bacteria', 'Viruses')) %>% # Adding an extra column defining if replicate is from bacterial or viral origin
    arrange('Location', 'Station_Number', 'Depth', 'Sample_Type','Replicate','Population', as.numeric(Timepoint)) # Reorder the rows by the values of the selected columns
  
  # By default, only go further with VP and VPC replicates => filtering out 0.22 samples which are control samples
  if (keep_0.22 == F){
    DF <- DF[DF$Sample_Type != '0.22',]
  }
  
  # Adding timepoints
  for (location in unique(DF$Location)){
    for (station in unique(DF$Station_Number)){
      for (depth in unique(DF$Depth)){
        DF <- DF %>%
          filter(Location == location & Station_Number == station & Depth == depth) # Filter on location, station number and depth of sample
        
        DF <- tp(DF) # Add timepoints
      }
    }
  }
  
  # Change to data table for further analyses
  DF <- as.data.table(DF)
  
  return(DF)
}

# Average of replicates on same timepoint dataframe
# Since for the control samples (0.22) we only have timepoints 0 and 24, we need to filter these out
df_AVG <- function(data, keep_0.22 = F){
  # Removing control samples
  DF <- data[data$Sample_Type != '0.22',]
  
  # Calculating number, mean of replicates and standard error 
  DF <- DF %>%
    select(4:16) %>%
    gather(7:13, key = 'Population', value = 'Count') %>%
    group_by(Location, Station_Number, Depth, Sample_Type, Timepoint, Population) %>%
    summarise(n = n(), Mean = mean(Count), SE = plotrix::std.error(Count)) # plotrix::std.error calculates the standard error of the mean
    
  # Calculate the difference between VPC samples and VP samples since for some methods we will estimate the difference curve
  # We will consider the mean and se separately, so we will use two intermediate subdataframes
  DF_mean <- DF %>%
    select(-'SE') %>%
    spread('Sample_Type', 'Mean')
  
  DF_se <- DF %>%
    select(-'Mean') %>%
    spread('Sample_Type', 'SE')
  
  if ('VPC' %in% DF$Sample_Type){
    DF_mean$Diff <- with(DF_mean, VPC - VP)
    DF_mean <- pivot_longer(DF_mean, cols = c('VP', 'VPC', 'Diff'), names_to = 'Sample_Type', values_to = 'Mean')
    DF_se$Diff <- with(DF_se, VPC + VP)
    DF_se <- pivot_longer(DF_se, cols = c('VP', 'VPC', 'Diff'), names_to = 'Sample_Type', values_to = 'SE')
  }
  
  # Merging in one dataframe and adding columns Microbe and Subgroup
  DF <- merge(DF_mean, DF_se, by = c('Location', 'Station_Number', 'Depth',
                                     'Timepoint', 'Population', 'n', 'Sample_Type')) %>%
    mutate(Microbe = if_else(Population %in% c('c_Bacteria', 'c_HNA', 'c_LNA'), 'Bacteria', 'Viruses')) %>%
    mutate(Subgroup = if_else(Population %in% c('c_Bacteria', 'c_Viruses'), 'Parent', 'Subgroup')) %>%
    arrange('Location', 'Station_Number', 'Depth', 'Sample_Type', 'Population', as.numeric(Timepoint))
  
  # Adding timepoints
  for (location in unique(DF$Location)){
    for (station in unique(DF$Station_Number)){
      for (depth in unique(DF$Depth)){
        DF <- DF %>%
          filter(Location == location & Station_Number == station & Depth == depth)
        
        DF <- tp(DF)
      }
    }
  }
  # Change to data table for further analyses
  DF <- as.data.table(DF)
  
  return(DF)
}

## 2. Peaks and valleys
# Two main methods will be used for "Viral Production Analyses": Linear Regression vs VIPCAL
# Linear regression calculates the viral production using a slope
# On the other hand, VIPCAL uses average of increments. Therefore, peaks and valleys need to be determined based of the counts

# Determining peaks
peaks <- function(values){
  list <- c()

  for (i in 1:(length(values)-1)){ # -1 because we add +10e+10 and -10e+10 to values so that first and last count is not dismissed => if -1 is not presented, last element of list will be NA
    d <- sign(values[i+1] - values[i]) # Determine where the index changes by computing the next count with the current count
    list[length(list)+1] <- d # List with 3 possible values: -1,0,1
    # 0 = no change of count; -1 = counts go down: 1 = counts go up
  }
  # If value is -1, current index is a peak; If value is 1, next index is peak
  return(which(diff(list) < 0)) # Which() looks where the diff(list) is smaller then 0 => diff(list) determines the difference between consecutive values of the list => if difference is negative it means that next value of the list is lower then current => PEAK
}
  
# Determining valleys
valleys <- function(values){
  list <- c()
  
  for (i in 1:(length(values)-1)){
    d <- sign(values[i+1] - values[i]) 
    list[[length(list)+1]] <- d 
  }
  return(which(diff(as.numeric(list)) > 0))
}

# Taking the SE into account
peaks_se <- function(values){
  list <- c()
  
  for (i in 1:(length(values)-1)){
    d <- sign((values[i+1] - se[i+1]) - (values[i] + se[i])) 
    list[[length(list)+1]] <- d 
  }
  return(which(diff(as.numeric(list)) < 0))
}

valleys_se <- function(values){
  list <- c()
  
  for (i in 1:(length(values)-1)){
    d <- sign((values[i+1] - se[i+1]) - (values[i] + se[i])) 
    list[[length(list)+1]] <- d 
  }
  return(which(diff(as.numeric(list)) > 0))
}

## 3. LMER model
# Some of the methods use a difference curve to determine the lysogenic production => Lysogenic production = average of increments in difference curve
# This difference curve can be easily calculated by subtracting the VP slope of the VPC slope or
# LMER model can be used (Linear Mixed-Effects Model)
LMER_model <- function(DF){ # Dataframe as input 
  lmer_data <- data.frame() # Initialize result dataframe
  
  for (spec in unique(DF$Population)){
    # Fit linear mixed-effects model
    # lm_function of form: resp ~ expr
    # lmer_function of form: resp ~ FEexpr + (REexpr1 | factor1) + (REexpr2 | factor2) => contain special random-effects terms
    # The LMER model will take the different replicates into account as random-effects terms
    lmer_mod <- lmer(data = DF, Count ~ Sample_Type * as.factor(Timepoint) + (1 | Replicate)) # Does the count change in function of Sample_Type and Timepoint (interaction term) + random-effect: variability between different levels of Replicate => allowing variation between different Replicates
    
    # Compute estimate marginal means in LMER model => least-squares means
    esmmeans <- emmeans(lmer_mod, ~ Sample_Type | as.factor(Timepoint)) # VP and VPC samples
    
    res1 <- data.frame(rep(spec, length(unique(DF$Timepoint))),
                       summary(esmmeans)$Sample_Type,
                       summary(esmmeans)$Timepoint,
                       summary(esmmeans)$emmean,
                       summary(esmmeans)$SE)
    colnames(res1) <- c("Population", "Sample_Type", "Timepoint", "Mean", "SE")
    
    # Contrasts of least-squares means => VP - VPC = difference samples
    contrasts <- pairs(esmmeans)
    
    res2 <- data.frame(rep(spec, length(unique(DF$Timepoint))), # rep() = replicates the values
                       rep("Diff", length(unique(DF$Timepoint))),
                       summary(contrasts)$Timepoint,
                       -(summary(contrasts)$estimate),
                       summary(contrasts)$SE)
    colnames(res2) <- c("Population", "Sample_Type", "Timepoint", "Mean", "SE")
    
    # All together in one dataframe
    lmer_data <- rbind(lmer_data, res1, res2)
  }
  return(lmer_data)
}

## 4. Slope Functions
# Functions that will help calculate the slope between the different timepoints for linear regression or the increments for VIPCAL 

## 4.1 Linear Regression: 5 variants of methods => 4 different functions to determine slope
# F1: Slope of all points => no replicate treatment
slope_LM_allpoints <- function(DF_SR){ # Takes separate replicate dataframe as input, but won't filter on replicates
  lm_res <- list() # Initialize list for results
  
  for (location in unique(DF_SR$Location)){
    for (station in unique(DF_SR$Station_Number)){
      for (depth in unique(DF_SR$Depth)){
        for (virus in unique(DF_SR[DF_SR$Microbe == 'Viruses',]$Population)){ # Only want the virus samples, within virus samples look at different populations
          for (sample in unique(DF_SR$Sample_Type)){
            DF <- DF_SR %>%
              filter(Location == location & Station_Number == station & Depth == depth,
                     Population == virus, Sample_Type == sample) # Get all samples of one sample type, one population => Last iteration gives population c_V3, sample_type VPC
            
            # For each subselection of samples, another subselection based on Time_Range
            for (time in unique(DF$Time_Range)){
              DF2 <- DF %>%
                filter(Time_Range == time) # First iteration = c_Viruses, VP, Timepoint 0 and 3; 2nd iteration = c_Viruses, VP, Timepoint 0,3 and 6; ...
              
              # Fit linear model
              lm <- summary(lm(data = DF2, Count ~ as.numeric(Timepoint))) # Does the count of the current samples significantly changes over time?
              # Coefficient of Timepoint represents the average change in Count for each unit increase in Timepoint => if this coefficient is significant (pval < 0.05) => linear relation = Count significantly changes over Timepoint
              # Summary will provide other information such as standard errors, t- and p-values
              
              # Save the coefficient and its standard error together with the multiple R-squared value which is a measure of how well the linear regression model fits the data
              # A higher R-squared value indicates a better fit of the model to the data => the independent variable accounts for a larger proportion of the variability in the data
              # P-value = significance individual variable; R-squared = overall goodness-of-fit of linear model
              res <- c(location, station, depth, time, virus, sample, lm$coefficients[c(2,4)], lm$r.squared)
              lm_res[[length(lm_res)+1]] <- res
            }
          }
        }
      }
    }
  }
  # Changing nested list into dataframe
  S_LM_allpoints <- data.frame(t(sapply(lm_res, c))) # sapply will return a column vector of each nested list of lm_res, transposing this and setting to dataframe for result
  colnames(S_LM_allpoints) <- c('Location', 'Station_Number', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'VP_Slope', 'VP_SE', 'VP_R_Squared')
  
  # Results of linear model need to be changed to numeric instead of character
  S_LM_allpoints[, c('VP_Slope', 'VP_SE', 'VP_R_Squared')] <- lapply(S_LM_allpoints[, c('VP_Slope', 'VP_SE', 'VP_R_Squared')], as.numeric)
  
  return(S_LM_allpoints)
}

# F2: Slope with separate replicate treatment
slope_LM_rep <- function(DF_SR){ # Takes separate replicate dataframe as input, with filter on the replicates
  lm_res <- list()
  
  for (location in unique(DF_SR$Location)){
    for (station in unique(DF_SR$Station_Number)){
      for (depth in unique(DF_SR$Depth)){
        for (virus in unique(DF_SR[DF_SR$Microbe == 'Viruses',]$Population)){ 
          for (sample in unique(DF_SR$Sample_Type)){
            for (rep in unique(DF_SR$Replicate)){
              DF <- DF_SR %>%
                filter(Location == location & Station_Number == station & Depth == depth,
                       Population == virus, Sample_Type == sample, Replicate == rep)
              
              for (time in unique(DF$Time_Range)){
                DF2 <- DF %>%
                  filter(Time_Range == time)
                
                # Fit linear model
                lm <- summary(lm(data = DF2, Count ~ as.numeric(Timepoint)))
                res <- c(location, station, depth, time, virus, sample, rep, lm$coefficients[c(2,4)], lm$r.squared)
                lm_res[[length(lm_res)+1]] <- res
              }
            }
          }
        }
      }
    }
  }
  # Changing nested list into dataframe
  S_LM_rep <- data.frame(t(sapply(lm_res, c))) 
  colnames(S_LM_rep) <- c('Location', 'Station_Number', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'Replicate', 'VP_Slope', 'VP_SE', 'VP_R_Squared')
  
  # Results of linear model need to be changed to numeric instead of character
  S_LM_rep[, c('VP_Slope', 'VP_SE', 'VP_R_Squared')] <- lapply(S_LM_rep[, c('VP_Slope', 'VP_SE', 'VP_R_Squared')], as.numeric)
  
  return(S_LM_rep)
}

# F3: Slope with average replicate treatment
slope_LM_avg <- function(DF_AVG){ # Takes average replicate dataframe as input
  lm_res <- list()
  
  for (location in unique(DF_AVG$Location)){
    for (station in unique(DF_AVG$Station_Number)){
      for (depth in unique(DF_AVG$Depth)){
        for (virus in unique(DF_AVG[DF_AVG$Microbe == 'Viruses',]$Population)){ 
          for (sample in unique(DF_AVG$Sample_Type)){
            DF <- DF_AVG %>%
              filter(Location == location & Station_Number == station & Depth == depth,
                     Population == virus, Sample_Type == sample)
            
            for (time in unique(DF$Time_Range)){
              DF2 <- DF %>%
                filter(Time_Range == time)
              
              # Fit linear model
              lm <- summary(lm(data = DF2, Mean ~ as.numeric(Timepoint)))
              res <- c(location, station, depth, time, virus, sample, lm$coefficients[c(2,4)], lm$r.squared)
              lm_res[[length(lm_res)+1]] <- res
            }
          }
        }
      }
    }
  }
  # Changing nested list into dataframe
  S_LM_avg <- data.frame(t(sapply(lm_res, c))) 
  colnames(S_LM_avg) <- c('Location', 'Station_Number', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'VP_Slope', 'VP_SE', 'VP_R_Squared')
  
  # Results of linear model need to be changed to numeric instead of character
  S_LM_avg[, c('VP_Slope', 'VP_SE', 'VP_R_Squared')] <- lapply(S_LM_avg[, c('VP_Slope', 'VP_SE', 'VP_R_Squared')], as.numeric)
  
  return(S_LM_avg)
}  

# F4: Slope with LMER model for difference curve
slope_LM_avg_diff_lmer <- function(DF_SR){# Takes separate replicate dataframe as input, but we still work with averaged replicates => averaging is in LMER model
  lm_res <- list()
  
  for (location in unique(DF_SR$Location)){
    for (station in unique(DF_SR$Station_Number)){
      for (depth in unique(DF_SR$Depth)){
        for (virus in unique(DF_SR[DF_SR$Microbe == 'Viruses',]$Population)){
          DF <- DF_SR %>%
            filter(Location == location, Station_Number == station,
                   Depth == depth, Population == virus)
          
          for (time in unique(DF$Time_Range)){
            DF2 <- DF %>%
              filter(Time_Range == time)
            
            # Currently, still working with separate replicate dataframe => going to average by using LMER model and adding difference samples
            DF3 <- LMER_model(DF2)
            
            # Fit linear model
            for (sample in unique(DF3$Sample_Type)){
              lm <- summary(lm(data = DF3[DF3$Sample_Type == sample, ], Mean ~ as.numeric(Timepoint)))
              res <- c(location, station, depth, time, virus, sample, lm$coefficients[c(2,4)], lm$r.squared)
              lm_res[[length(lm_res)+1]] <- res
            }
          }
        }
      }
    }
  }
  # Changing nested list into dataframe
  S_LM_avg_diff_lmer <- data.frame(t(sapply(lm_res, c))) 
  colnames(S_LM_avg_diff_lmer) <- c('Location', 'Station_Number', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'VP_Slope', 'VP_SE', 'VP_R_Squared')
  
  # Results of linear model need to be changed to numeric instead of character
  S_LM_avg_diff_lmer[, c('VP_Slope', 'VP_SE', 'VP_R_Squared')] <- lapply(S_LM_avg_diff_lmer[, c('VP_Slope', 'VP_SE', 'VP_R_Squared')], as.numeric)
  
  return(S_LM_avg_diff_lmer)
}

## 4.2 VIPCAL: 7 variants of methods => 5 different functions to determine increments
# F1: Separate replicate treatment
vipcal_rep <- function(DF_SR){ # Takes separate replicate dataframe as input, with filter on the replicates
  vipcal_res <- list() # Initialize list for results
  
  for (location in unique(DF_SR$Location)){
    for (station in unique(DF_SR$Station_Number)){
      for (depth in unique(DF_SR$Depth)){
        for (virus in unique(DF_SR[DF_SR$Microbe == 'Viruses',]$Population)){ 
          for (sample in unique(DF_SR$Sample_Type)){
            for (rep in unique(DF_SR$Replicate)){
              DF <- DF_SR %>%
                filter(Location == location & Station_Number == station & Depth == depth,
                       Population == virus, Sample_Type == sample, Replicate == rep)
              
              for (time in unique(DF$Time_Range)) {
                DF2 <- DF %>%
                  filter(Time_Range == time)
                
                # VIPCAL: determing average of increments based on peaks and valleys
                
              }
            }
          }
        }
      }
    }
  }
}


# F2: Average replicate treatment

# F3: Average replicate treatment with SE

# F4: Average replicate treatment with LMER model for difference curve

# F5: Average replicate treatment with LMER model for difference curve with SE
























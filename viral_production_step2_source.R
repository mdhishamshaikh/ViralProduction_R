###---Source file for viral_production_step2.R---###
# General setup + collection of functions needed in methods of vp_calc_functions

### Information
# Viral production can be divided in two phases: Lytic and Lysogenic
# VP replicates represent the lytic viral production, this can be calculated with linear regression (slope) or VIPCAL (average of increments)
# Analoge for VPC samples => they represent both lytic and lysogenic viral production
# To determine only the lysogenic viral production: in linear regression VPC slope - VP slope, in VIPCAL average of increments of the difference curve
# VIPCAL has his own issues => standard error has big influence => VIPCAL-SE takes this into account and only looks at  those increments that don't have overlapping range of SE
# If there is overlap we give a 0 => rather a zero than an uncertain number => VIPCAL-SE is more conservative

# VIPCAL = online tool for estimating lyticaly and lysogenicaly produced viruses
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
    drop_na() # Drops all the rows which have a NA value
  # When running only this function as DF <- tp(data), an empty dataframe will be returned because of the column "Comments", which consists of NA
  
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
df_AVG <- function(data){
  # Removing control samples
  DF <- data[data$Sample_Type != '0.22',]
  
  # Calculating number, mean of replicates and standard error 
  DF <- DF %>%
    select(4:16) %>%
    gather(7:13, key = 'Population', value = 'Count') %>%
    group_by(Location, Station_Number, Depth, Sample_Type, Timepoint, Population) %>%
    summarise(n = n(), Mean = mean(Count), SE = plotrix::std.error(Count)) # plotrix::std.error calculates the standard error of the mean
    
  # Calculate the difference between VPC samples and VP samples since for some methods viral production is determined based on the difference curve
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
    # 0 = no change of count; -1 = counts go down; 1 = counts go up
  }
  return(which(diff(list) < 0)) # PEAK if negative difference
}
  
# Determining valleys
valleys <- function(values){
  list <- c()
  
  for (i in 1:(length(values)-1)){ 
    d <- sign(values[i+1] - values[i]) 
    list[length(list)+1] <- d 
    
  }
  return(which(diff(list) > 0)) # VALLEY if positive difference
}

# Taking the SE into account: only TRUE increments are kept
# If the SE between a peak and valley overlaps, can't guarantee that this difference is sufficient => keep those out
# Only peak of valley if there is a difference and standard errrors don't overlap
# By doing this also the created valley/peak by adding 10e10 and -10e10 is filtered out
# Negative evolution in viral production also filtered out => if you go from peak -> valley (downward) => only increments are kept
peaks_se <- function(values, se){
  list <- c()
  
  for (i in 1:(length(values)-1)){
    d <- sign((values[i+1] - se[i+1]) - (values[i] + se[i])) # Take SE into account by adding it to current count and substract it of next count
    list[length(list)+1] <- d # List with two possible values: -1 and 1
    # Only a positive sign if increments with no overlapping standard errors
  }
  return(which(diff(list) < 0)) # PEAK if negative difference
}

valleys_se <- function(values, se){
  list <- c()
  
  for (i in 1:(length(values)-1)){
    d <- sign((values[i+1] - se[i+1]) - (values[i] + se[i])) 
    list[length(list)+1] <- d 
  }
  return(which(diff(list) > 0)) # VALLEY if positive difference
}

## 3. LMER model
# To determine lysogenic production, a difference curve is used => Lysogenic production = average of increments in difference curve
# This difference curve can be easily calculated by subtracting the VP slope of the VPC slope or a
# LMER model can be used (Linear Mixed-Effects Model)
LMER_model <- function(DF){ # Dataframe as input 
  lmer_data <- data.frame() # Initialize result dataframe
  
  # The LMER model will take the different replicates into account as random-effects terms
  # Solely on replicates column, the model will have problems distinguishing VP and VPC samples: for both rep = 1, 2 or 3
  # In order to make sure the LMER model distinguishes between them, change replicate reference of VPC samples
  # By changing this, an effective separation between both VP and VPC samples are made. 
  # The LMER model can now tike into account the Sample_Type factor but also the interaction between Sample_Type and the new replicate reference
  for (rep in unique(DF$Replicate)){
    DF$Replicate[DF$Replicate == rep & DF$Sample_Type == 'VPC'] <- rep + 3
  }
  
  for (spec in unique(DF$Population)){
    # Fit linear mixed-effects model
    # lm_function of form: resp ~ expr
    # lmer_function of form: resp ~ FEexpr + (REexpr1 | factor1) + (REexpr2 | factor2) => contain special random-effects terms
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
  colnames(S_LM_allpoints) <- c('Location', 'Station_Number', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'VP', 'VP_SE', 'VP_R_Squared')
  
  # Results of linear model need to be changed to numeric instead of character
  S_LM_allpoints[, c('VP', 'VP_SE', 'VP_R_Squared')] <- lapply(S_LM_allpoints[, c('VP', 'VP_SE', 'VP_R_Squared')], as.numeric)
  
  return(S_LM_allpoints)
}

# F2: Slope with separate replicate treatment
slope_LM_sr <- function(DF_SR){ # Takes separate replicate dataframe as input, with filter on the replicates
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
  colnames(S_LM_rep) <- c('Location', 'Station_Number', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'Replicate', 'VP', 'VP_SE', 'VP_R_Squared')
  
  # Results of linear model need to be changed to numeric instead of character
  S_LM_rep[, c('VP', 'VP_SE', 'VP_R_Squared')] <- lapply(S_LM_rep[, c('VP', 'VP_SE', 'VP_R_Squared')], as.numeric)
  
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
  colnames(S_LM_avg) <- c('Location', 'Station_Number', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'VP', 'VP_SE', 'VP_R_Squared')
  
  # Results of linear model need to be changed to numeric instead of character
  S_LM_avg[, c('VP', 'VP_SE', 'VP_R_Squared')] <- lapply(S_LM_avg[, c('VP', 'VP_SE', 'VP_R_Squared')], as.numeric)
  
  return(S_LM_avg)
}  

# F4: Slope with LMER model for difference curve
slope_LM_avg_lmer <- function(DF_SR){# Takes separate replicate dataframe as input, but we still work with averaged replicates => averaging is in LMER model
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
  colnames(S_LM_avg_diff_lmer) <- c('Location', 'Station_Number', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'VP', 'VP_SE', 'VP_R_Squared')
  
  # Results of linear model need to be changed to numeric instead of character
  S_LM_avg_diff_lmer[, c('VP', 'VP_SE', 'VP_R_Squared')] <- lapply(S_LM_avg_diff_lmer[, c('VP', 'VP_SE', 'VP_R_Squared')], as.numeric)
  
  return(S_LM_avg_diff_lmer)
}

## 4.2 VIPCAL: 7 variants of methods => 5 different functions to determine increments
# F1: Separate replicate treatment
VIPCAL_sr <- function(DF_SR){ # Takes separate replicate dataframe as input, with filter on the replicates
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
                
                # VIPCAL
                # Determine peaks and valleys
                p <- peaks(c(+10e+10, DF2$Count, -10e+10)) # Adding to make sure that first and last element are not dismissed
                v <- valleys(c(+10e+10, DF2$Count, -10e+10))
                
                # Check if number of peaks and valleys correspond
                if (identical(length(p), length(v))){
                  print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
                } else {
                  print("Number of peaks and valleys identified differ, this will lead to erroneous viral production calculations")
                }
                
                # Calculate viral production: as the slope between the minimum (valley) and maximum (peak) viral abundance
                # Following function holds: {[(Viral_Count(peak1) - Viral_Count(valley1)) / (Time(peak1) - Time(valley1))] + [...] + [...]} / amount of peaks
                if (length(p) == 0){
                  vp <- 0 # No viral production if no peaks => constant viral abundance
                } else {
                  total_vp <- 0
                  for (i in 1:length(p)){ # Iterate through the peaks and calculate the viral production of peak and corresponding value
                    vp_i <- (DF2$Count[p[i]] - DF2$Count[v[i]]) / (DF2$Timepoint[p[i]] - DF2$Timepoint[v[i]])
                    total_vp <- total_vp + vp_i
                  }
                  vp <- total_vp / length(p) # Average of increments
                }
                res <- c(location, station, depth, time, virus, sample, rep, vp)
                vipcal_res[[length(vipcal_res)+1]] <- res
              }
            }
          }
        }
      }
    }
  }
  # Changing nested list into dataframe
  VPCL_sr <- data.frame(t(sapply(vipcal_res, c)))
  colnames(VPCL_sr)<- c('Location', 'Station_Number', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'Replicate', 'VP')
  VPCL_sr$VP <- as.numeric(VPCL_sr$VP) # Transform viral production to a numeric value instead a character
  
  return(VPCL_sr)
}

# F2: Average replicate treatment
VIPCAL_avg <- function(DF_AVG){ # Takes average replicate dataframe as input
  vipcal_res <- list()
  
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
              
              # VIPCAL
              # Determine peaks and valleys
              p <- peaks(c(+10e+10, DF2$Mean, -10e+10))
              v <- valleys(c(+10e+10, DF2$Mean, -10e+10))
              
              # Check if number of peaks and valleys correspond
              if (identical(length(p), length(v))){
                print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
              } else {
                print("Number of peaks and valleys identified differ, this will lead to erroneous viral production calculations")
              }
              
              # Calculate viral production: as the slope between the minimum (valley) and maximum (peak) viral abundance
              # Following function holds: {[(Viral_Count(peak1) - Viral_Count(valley1)) / (Time(peak1) - Time(valley1))] + [...] + [...]} / amount of peaks
              if (length(p) == 0){
                vp <- 0 
              } else {
                total_vp <- 0
                for (i in 1:length(p)){ 
                  vp_i <- (DF2$Mean[p[i]] - DF2$Mean[v[i]]) / (DF2$Timepoint[p[i]] - DF2$Timepoint[v[i]])
                  total_vp <- total_vp + vp_i
                }
                vp <- total_vp / length(p)
              }
              res <- c(location, station, depth, time, virus, sample, vp)
              vipcal_res[[length(vipcal_res)+1]] <- res
            }
          }
        }
      }
    }
  }
  # Changing nested list into dataframe
  VPCL_avg <- data.frame(t(sapply(vipcal_res, c)))
  colnames(VPCL_avg)<- c('Location', 'Station_Number', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'VP')
  VPCL_avg$VP <- as.numeric(VPCL_avg$VP)
  
  return(VPCL_avg)
}

# F3: Average replicate treatment with SE
VIPCAL_avg_se <- function(DF_AVG){ # Takes average replicate dataframe as input
  vipcal_res <- list()
  
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
              
              # VIPCAL
              # Determine peaks and valleys with SE taken into account
              p <- peaks_se(c(+10e+10, DF2$Mean, -10e+10),
                            c(0, DF2$SE, 0))
              v <- valleys_se(c(+10e+10, DF2$Mean, -10e+10),
                              c(0, DF2$SE, 0))
              
              # Check if number of peaks and valleys correspond
              if (identical(length(p), length(v))){
                print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
              } else {
                print("Number of peaks and valleys identified differ, this will lead to erroneous viral production calculations")
              }
              
              # Calculate viral production: as the slope between the minimum (valley) and maximum (peak) viral abundance
              # Following function holds: {[(Viral_Count(peak1) - Viral_Count(valley1)) / (Time(peak1) - Time(valley1))] + [...] + [...]} / amount of peaks
              # Also calculating standard error on the viral production calculation
              if (length(p) == 0){
                vp <- 0 
                se <- 0
              } else {
                total_vp <- 0
                total_se <- 0
                for (i in 1:length(p)){ 
                  vp_i <- (DF2$Mean[p[i]] - DF2$Mean[v[i]]) / (DF2$Timepoint[p[i]] - DF2$Timepoint[v[i]])
                  total_vp <- total_vp + vp_i
                  
                  se_i <- (DF2$SE[p[i]] + DF2$SE[v[i]]) / (DF2$Timepoint[p[i]] - DF2$Timepoint[v[i]])
                  total_se <- total_se + se_i
                }
                vp <- total_vp / length(p)
                se <- total_se / length(p)
              }
              res <- c(location, station, depth, time, virus, sample, vp, se)
              vipcal_res[[length(vipcal_res)+1]] <- res
            }
          }
        }
      }
    }
  }
  # Changing nested list into dataframe
  VPCL_avg_se <- data.frame(t(sapply(vipcal_res, c)))
  colnames(VPCL_avg_se)<- c('Location', 'Station_Number', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'VP', 'VP_SE')
  VPCL_avg_se[, c('VP', 'VP_SE')] <- lapply(VPCL_avg_se[, c('VP', 'VP_SE')], as.numeric)
  
  return(VPCL_avg_se)
}

# F4: Average replicate treatment with LMER model for difference curve
VIPCAL_avg_lmer <- function(DF_SR){ # Takes separate replicate dataframe as input, but we still work with averaged replicates => averaging is in LMER model
  vipcal_res <- list()
  
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
            
            # VIPCAL
            for (sample in unique(DF3$Sample_Type)){
              DF4 <- DF3 %>%
                filter(Sample_Type == sample)
              
              # Determine peaks and valleys
              p <- peaks(c(+10e+10, DF4$Mean, -10e+10))
              v <- valleys(c(+10e+10, DF4$Mean, -10e+10))
              
              # Check if number of peaks and valleys correspond
              if (identical(length(p), length(v))){
                print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
              } else {
                print("Number of peaks and valleys identified differ, this will lead to erroneous viral production calculations")
              }
              
              # Calculate viral production: as the slope between the minimum (valley) and maximum (peak) viral abundance
              # Following function holds: {[(Viral_Count(peak1) - Viral_Count(valley1)) / (Time(peak1) - Time(valley1))] + [...] + [...]} / amount of peaks
              if (length(p) == 0){
                vp <- 0 
              } else {
                total_vp <- 0
                for (i in 1:length(p)){ 
                  vp_i <- (DF4$Mean[p[i]] - DF4$Mean[v[i]]) / (DF4$Timepoint[p[i]] - DF4$Timepoint[v[i]])
                  total_vp <- total_vp + vp_i
                }
                vp <- total_vp / length(p)
              }
              res <- c(location, station, depth, time, virus, sample, vp)
              vipcal_res[[length(vipcal_res)+1]] <- res
            }
          }
        }
      }
    }
  }
  # Changing nested list into dataframe
  VPCL_avg_lmer <- data.frame(t(sapply(vipcal_res, c)))
  colnames(VPCL_avg_lmer)<- c('Location', 'Station_Number', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'VP')
  VPCL_avg_lmer$VP <- as.numeric(VPCL_avg_lmer$VP)
  
  return(VPCL_avg_lmer)
}

# F5: Average replicate treatment with LMER model for difference curve with SE
VIPCAL_avg_lmer_se <- function(DF_SR){ # Takes separate replicate dataframe as input, but we still work with averaged replicates => averaging is in LMER model
  vipcal_res <- list()
  
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
            
            # VIPCAL
            for (sample in unique(DF3$Sample_Type)){
              DF4 <- DF3 %>%
                filter(Sample_Type == sample)
              
              # Determine peaks and valleys
              p <- peaks_se(c(+10e+10, DF4$Mean, -10e+10),
                            c(0, DF4$SE, 0))
              v <- valleys_se(c(+10e+10, DF4$Mean, -10e+10),
                              c(0, DF4$SE, 0))
              
              # Check if number of peaks and valleys correspond
              if (identical(length(p), length(v))){
                print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
              } else {
                print("Number of peaks and valleys identified differ, this will lead to erroneous viral production calculations")
              }
              
              # Calculate viral production: as the slope between the minimum (valley) and maximum (peak) viral abundance
              # Following function holds: {[(Viral_Count(peak1) - Viral_Count(valley1)) / (Time(peak1) - Time(valley1))] + [...] + [...]} / amount of peaks
              # Also calculating standard error on the viral production calculation
              if (length(p) == 0){
                vp <- 0 
                se <- 0
              } else {
                total_vp <- 0
                total_se <- 0
                for (i in 1:length(p)){ 
                  vp_i <- (DF4$Mean[p[i]] - DF4$Mean[v[i]]) / (DF4$Timepoint[p[i]] - DF4$Timepoint[v[i]])
                  total_vp <- total_vp + vp_i
                  
                  se_i <- (DF4$SE[p[i]] + DF4$SE[v[i]]) / (DF4$Timepoint[p[i]] - DF4$Timepoint[v[i]])
                  total_se <- total_se + se_i
                }
                vp <- total_vp / length(p)
                se <- total_se / length(p)
              }
              res <- c(location, station, depth, time, virus, sample, vp, se)
              vipcal_res[[length(vipcal_res)+1]] <- res
            }
          }
        }
      }
    }
  }
  # Changing nested list into dataframe
  VPCL_avg_lmer_se <- data.frame(t(sapply(vipcal_res, c)))
  colnames(VPCL_avg_lmer_se)<- c('Location', 'Station_Number', 'Depth', 'Time_Range', 'Population', 'Sample_Type', 'VP', 'VP_SE')
  VPCL_avg_lmer_se[, c('VP', 'VP_SE')] <- lapply(VPCL_avg_lmer_se[, c('VP', 'VP_SE')], as.numeric)
  
  return(VPCL_avg_lmer_se)
}

## 5. Lysogenic production
# Above functions provide to calculate the viral production for each of the samples
# VP samples represent lytic viral production there where VPC samples represent lytic and lysogenic viral production
# Next function will make it able to calculate the lysogenic viral production for the methods that don't use a difference curve

# I think it is possible to make one function instead of 4 and just add some variables which distinguish the different cases
calc_DIFF <- function(DF, VIPCAL = F, SE = F){ # Input dataframe consist of output of slope functions for Linear regression or VIPCAL functions of VIPCAL
  # Determine difference samples by subtraction
  # For linear model always the same, for VIPCAL variants with and without SE. R_Squared values are only presented with linear regression, not with VIPCAL
  if (VIPCAL == T){
    if (SE == F){
      diff <- DF %>%
        group_by(Location, Station_Number, Depth, Time_Range, Population) %>%
        summarize(
          VP = sum(VP[Sample_Type == "VPC"]) - sum(VP[Sample_Type == "VP"]), # Calculate viral production for difference samples
          Sample_Type = "Diff")
    }else if (SE == T){
      diff <- DF %>%
        group_by(Location, Station_Number, Depth, Time_Range, Population) %>%
        summarize(
          VP = sum(VP[Sample_Type == "VPC"]) - sum(VP[Sample_Type == "VP"]), # Calculate viral production for difference samples
          VP_SE = sum(VP_SE[Sample_Type == "VPC"]) + sum(VP_SE[Sample_Type == "VP"]), # Calculate the standard error on the viral production
          Sample_Type = "Diff")
    }
  }else{ 
    diff <- DF %>%
      group_by(Location, Station_Number, Depth, Time_Range, Population) %>%
      summarize(
        VP = sum(VP[Sample_Type == "VPC"]) - sum(VP[Sample_Type == "VP"]), # Calculate viral production for difference samples
        VP_SE = sum(VP_SE[Sample_Type == "VPC"]) + sum(VP_SE[Sample_Type == "VP"]), # Calculate the standard error on the viral production
        VP_R_Squared = NA,
        Sample_Type = "Diff")
  }
  # Join diff samples to VP and VPC samples
  joined_df <- full_join(DF, diff, by = NULL)
  
  return(joined_df)
}

## 6. Different methods for calculating viral production
# The viral production is calculated by the means of two main methods: Linear Regression vs VIPCAL
# Of these methods different variants are applied considering the replicate treatment, standard error taken into account or difference curve estimation
# In total, there are 5 variants of linear regression and 7 of VIPCAL

# 6.1 LM-1: Linear regression with no replicate treatment (LM_AP)
LM_1 <- function(data){
  # Create correct type of dataframe
  DF <- df_SR(data)
  
  # Calculate viral production
  DF_VP <- slope_LM_allpoints(DF) %>% 
    calc_DIFF()
  
  # Arrange dataframe
  DF_LM1 <- DF_VP %>%
    group_by(Location, Station_Number, Depth, Time_Range, Population, Sample_Type) %>%
    arrange('Location', 'Station_Number', 'Depth',
            factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
            factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
            factor(Time_Range, levels = c("T0_T3", "T0_T6", "T0_T17", "T0_T20", "T0_T24"))) %>% # Added cause order for diff samples was gone (T17, T20, T24, T3, T6)
    mutate(VP_Type = 'LM_AP')
  
  return(DF_LM1)
}

# 6.2 LM-2: Linear regression with separate replicate treatment (LM_SR)
# In the end, we are interested in looking at lytic and lysogenic viral production. To accomplish this, difference samples are needed (lysogenic production)
# This method uses a separate replicate treatment, for each sample (VP and VPC) it distinguishes between the separate replicates
# Here, an additional parameter 'avg' is added. If avg = F: the output will consist only the VP and VPC samples with separate replicate treatment; if avg = T: after determining slopes with separate replicate treatment, data is averaged and difference samples are calculated
LM_2 <- function(data, avg = F){
  # Create correct type of dataframe
  DF <- df_SR(data)
  
  # Calculate viral production
  DF_VP <- slope_LM_sr(DF)
  
  # Arrange dataframe
  if (avg == F){
    DF_LM2 <- DF_VP %>%
      group_by(Location, Station_Number, Depth, Time_Range, Population, Sample_Type) %>%
      arrange('Location', 'Station_Number', 'Depth',
              factor(Sample_Type, levels = c('VP', 'VPC')),
              factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3'))) %>%
      mutate(VP_Type = 'LM_SR')
  } else if (avg == T){
    DF_LM2 <- DF_VP %>%
      group_by(Location, Station_Number, Depth, Time_Range, Population, Sample_Type) %>%
      summarise(VP_mean = mean(VP), VP_SE = plotrix::std.error(VP), VP_R_Squared = mean(VP_R_Squared)) %>%
      rename(VP = VP_mean) %>% # Necessary to name VP_mean in previous step so that SE is calculated correctly
      calc_DIFF() %>%
      arrange('Location', 'Station_Number', 'Depth',
              factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
              factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
              factor(Time_Range, levels = c("T0_T3", "T0_T6", "T0_T17", "T0_T20", "T0_T24"))) %>% 
      mutate(VP_Type = 'LM_SR_AVG')
  }
  
  return(DF_LM2)
}

# 6.3 LM-3: Linear regression with averaged replicate treatment (LM_AR)
LM_3 <- function(data){
  # Create correct type of dataframe
  DF <- df_AVG(data) %>%
    subset(Sample_Type != 'Diff') # LM-3 has no difference curve estimation
  
  # Calculate viral production
  DF_VP <- slope_LM_avg(DF) %>%
    calc_DIFF()
  
  # Arrange dataframe
  DF_LM3 <- DF_VP %>%
    group_by(Location, Station_Number, Depth, Time_Range, Population, Sample_Type) %>%
    arrange('Location', 'Station_Number', 'Depth',
            factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
            factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
            factor(Time_Range, levels = c("T0_T3", "T0_T6", "T0_T17", "T0_T20", "T0_T24"))) %>% 
    mutate(VP_Type = 'LM_AR')
  
  return(DF_LM3)
}

# 6.4 LM-4: Linear regression with averaged replicate treatment and difference curve estimation by subtraction (LM_AR_DIFF)
LM_4 <- function(data){
  # Create correct type of dataframe
  DF <- df_AVG(data) 
  
  # Calculate viral production
  DF_VP <- slope_LM_avg(DF)
  
  # Arrange dataframe
  DF_LM4 <- DF_VP %>%
    group_by(Location, Station_Number, Depth, Time_Range, Population, Sample_Type) %>%
    arrange('Location', 'Station_Number', 'Depth',
            factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
            factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
            factor(Time_Range, levels = c("T0_T3", "T0_T6", "T0_T17", "T0_T20", "T0_T24"))) %>% 
    mutate(VP_Type = 'LM_AR_DIFF')
  
  return(DF_LM4)
}

# 6.5 LM-5: Linear regression with averaged replicate treatment and difference curve estimation by LMER (LM_AR_DIFF_LMER)
LM_5 <- function(data){
  # Create correct type of dataframe
  DF <- df_SR(data) 
  
  # Calculate viral production
  DF_VP <- slope_LM_avg_lmer(DF)
  
  # Arrange dataframe
  DF_LM5 <- DF_VP %>%
    group_by(Location, Station_Number, Depth, Time_Range, Population, Sample_Type) %>%
    arrange('Location', 'Station_Number', 'Depth',
            factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
            factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
            factor(Time_Range, levels = c("T0_T3", "T0_T6", "T0_T17", "T0_T20", "T0_T24"))) %>% 
    mutate(VP_Type = 'LM_AR_DIFF_LMER')
  
  return(DF_LM5)
}

# 6.6 VPCL-1: VIPCAL with separate replicate treatment (VPCL_SR)
# Here the same as with LM-2, an additional parameter to determine if you want to see the separate replicate results or the averaged ones
VPCL_1 <- function(data, avg = F){
  # Create correct type of dataframe
  DF <- df_SR(data)
  
  # Calculate viral production
  DF_VP <- VIPCAL_sr(DF)
  
  # Arrange dataframe
  if (avg == F){
    DF_VPCL1 <- DF_VP %>%
      group_by(Location, Station_Number, Depth, Time_Range, Population, Sample_Type) %>%
      arrange('Location', 'Station_Number', 'Depth',
              factor(Sample_Type, levels = c('VP', 'VPC')),
              factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3'))) %>%
      mutate(VP_Type = 'VPCL_SR')
  }else if (avg == T){
    DF_VPCL1 <- DF_VP %>%
      group_by(Location, Station_Number, Depth, Time_Range, Population, Sample_Type) %>%
      summarise(VP_mean = mean(VP), VP_SE = plotrix::std.error(VP)) %>%
      rename(VP = VP_mean) %>% # Necessary to name VP_mean in previous step so that SE is calculated correctly
      calc_DIFF(VIPCAL = T) %>%
      arrange('Location', 'Station_Number', 'Depth',
              factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
              factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
              factor(Time_Range, levels = c("T0_T3", "T0_T6", "T0_T17", "T0_T20", "T0_T24"))) %>% 
      mutate(VP_Type = 'VPCL_SR_AVG')
  }
  
  return(DF_VPCL1)
}

# 6.7 VPCL-2: VIPCAL with averaged replicate treatment (VPCL_AR)
VPCL_2 <- function(data){
  # Create correct type of dataframe
  DF <- df_AVG(data) %>%
    subset(Sample_Type != 'Diff') # VPCL-2 has no difference curve estimation
  
  # Calculate viral production
  DF_VP <- VIPCAL_avg(DF) %>%
    calc_DIFF(VIPCAL = T)
  
  # Arrange dataframe
  VPCL2 <- DF_VP %>%
    group_by(Location, Station_Number, Depth, Time_Range, Population, Sample_Type) %>%
    arrange('Location', 'Station_Number', 'Depth',
            factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
            factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
            factor(Time_Range, levels = c("T0_T3", "T0_T6", "T0_T17", "T0_T20", "T0_T24"))) %>% 
    mutate(VP_Type = 'VPCL_AR')
  
  return(VPCL2)
}

# 6.8 VPCL-3: VIPCAL with averaged replicate treatment with SE (VPCL_AR_SE)
VPCL_3 <- function(data){
  # Create correct type of dataframe
  DF <- df_AVG(data) %>%
    subset(Sample_Type != 'Diff') # VPCL-3 has no difference curve estimation
  
  # Calculate viral production
  DF_VP <- VIPCAL_avg_se(DF) %>%
    calc_DIFF(VIPCAL = T, SE = T)
  
  # Arrange dataframe
  VPCL3 <- DF_VP %>%
    group_by(Location, Station_Number, Depth, Time_Range, Population, Sample_Type) %>%
    arrange('Location', 'Station_Number', 'Depth',
            factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
            factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
            factor(Time_Range, levels = c("T0_T3", "T0_T6", "T0_T17", "T0_T20", "T0_T24"))) %>% 
    mutate(VP_Type = 'VPCL_AR_SE')
  
  return(VPCL3)
}

# 6.9 VPCL-4: VIPCAL with averaged replicate treatment with difference curve estimation by subtraction (VPCL_AR_DIFF)
VPCL_4 <- function(data){
  # Create correct type of dataframe
  DF <- df_AVG(data)
  
  # Calculate viral production
  DF_VP <- VIPCAL_avg(DF)
  
  # Arrange dataframe
  VPCL4 <- DF_VP %>%
    group_by(Location, Station_Number, Depth, Time_Range, Population, Sample_Type) %>%
    arrange('Location', 'Station_Number', 'Depth',
            factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
            factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
            factor(Time_Range, levels = c("T0_T3", "T0_T6", "T0_T17", "T0_T20", "T0_T24"))) %>% 
    mutate(VP_Type = 'VPCL_AR_DIFF')
  
  return(VPCL4)
}

# 6.10 VPCL-5: VIPCAL with averaged replicate treatment with difference curve estimation by subtraction and with SE (VPCL_AR_DIFF_SE)
VPCL_5 <- function(data){
  # Create correct type of dataframe
  DF <- df_AVG(data)
  
  # Calculate viral production
  DF_VP <- VIPCAL_avg_se(DF)
  
  # Arrange dataframe
  VPCL5 <- DF_VP %>%
    group_by(Location, Station_Number, Depth, Time_Range, Population, Sample_Type) %>%
    arrange('Location', 'Station_Number', 'Depth',
            factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
            factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
            factor(Time_Range, levels = c("T0_T3", "T0_T6", "T0_T17", "T0_T20", "T0_T24"))) %>% 
    mutate(VP_Type = 'VPCL_AR_DIFF_SE')
  
  return(VPCL5)
}

# 6.11 VPCL-6: VIPCAL with averaged replicate treatment with difference curve estimation by LMER (VPCL_AR_DIFF_LMER)
VPCL_6 <- function(data){
  # Create correct type of dataframe
  DF <- df_SR(data)
  
  # Calculate viral production
  DF_VP <- VIPCAL_avg_lmer(DF)
  
  # Arrange dataframe
  VPCL6 <- DF_VP %>%
    group_by(Location, Station_Number, Depth, Time_Range, Population, Sample_Type) %>%
    arrange('Location', 'Station_Number', 'Depth',
            factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
            factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
            factor(Time_Range, levels = c("T0_T3", "T0_T6", "T0_T17", "T0_T20", "T0_T24"))) %>% 
    mutate(VP_Type = 'VPCL_AR_DIFF_LMER')
  
  return(VPCL6)
}

# 6.12 VPCL-7: VIPCAL with averaged replicate treatment with difference curve estimation by LMER and with SE (VPCL_AR_DIFF_LMER)
VPCL_7 <- function(data){
  # Create correct type of dataframe
  DF <- df_SR(data)
  
  # Calculate viral production
  DF_VP <- VIPCAL_avg_lmer_se(DF)
  
  # Arrange dataframe
  VPCL7 <- DF_VP %>%
    group_by(Location, Station_Number, Depth, Time_Range, Population, Sample_Type) %>%
    arrange('Location', 'Station_Number', 'Depth',
            factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
            factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
            factor(Time_Range, levels = c("T0_T3", "T0_T6", "T0_T17", "T0_T20", "T0_T24"))) %>% 
    mutate(VP_Type = 'VPCL_AR_DIFF_LMER_SE')
  
  return(VPCL7)
}

# List of all different methods
calculate_VP_list <- list(LM_1, LM_2, LM_3, LM_4, LM_5,
                          VPCL_1, VPCL_2, VPCL_3, VPCL_4,
                          VPCL_5, VPCL_6, VPCL_7)
names(calculate_VP_list) <- c('LM_1', 'LM_2', 'LM_3', 'LM_4', 'LM_5',
                              'VPCL_1', 'VPCL_2', 'VPCL_3', 'VPCL_4',
                              'VPCL_5', 'VPCL_6', 'VPCL_7')

# 7. Bacterial Endpoint
# Estimating moment where we should stop the assay
# Bacterial growth can induce an increase in collision rates
bacterial_endpoint <- function(data){ # Returns timepoint, where we should stop the assay
  
}


### Viral Production Assay: Step 2 ###

## 1. Importing used functions, these can be found in viral_production_step2_source.R
source("viral_production_step2_source.R")

## 2. Importing data from Step 1
data <- read.csv('NJ1.csv') %>%
  rename(Station_Number = Expt_No) # Changing column name to something more appropriate => you can also use index of column (5)

data_all <- read.csv('NJ2020.csv') %>%
  rename(Station_Number = Expt_No)

df_abundance <- read.csv('NJ2020_abundance.csv') %>% # Consist of the abundances of all populations in original seawater sample for each experiment
  rename(Station_Number = Expt_No)

## 3. Calculating viral production
# Main function for viral production calculation
# Different variables are presented to adjust for the desired output
calc_VP <- function(data, output_dir = '', method = c(1:12), write_csv = T, SR_calc = T, bp_endpoint = T){
  
  ## 1. Check if output directory is valid (variable: output_dir)
  if (output_dir == ''){
    print('No output directory is given!')
    stop('Please define output directory before proceeding.')
    
  }else if (file.exists(output_dir)){
    print(paste0('The ', output_dir, ' folder already exists!'))
    stop('Please define another output directory before proceeding.')
    
  }else {
    dir.create(output_dir)
  }
  
  # Storing all errors and warnings in global environment
  .GlobalEnv$calc_vp_error_list <- list()
  .GlobalEnv$calc_vp_warn_list <- list()
  
  ## 2. Calculate viral production
  # Create output dataframe
  res_path <- paste0('./', output_dir, '/') # output.dir needs to be in current working directory!
  output_df <- data.frame(Location = character(), Station_Number = character(), Depth = character(),
                          Time_Range = character(), Population = character(), Sample_Type = character(),
                          VP = numeric(), VP_SE = numeric(), VP_R_Squared = numeric(), VP_Type = character())
  
  # Run methods
  # Method is a vector consisting of numbers 1 to 12, this represents all the possible methods => see names(calculate_VP_list)
  # If only some of the functions are wanted, adjust variable method 
  for (mtd in method){
    # tryCatch() provides structured way to handle errors and warnings
    tryCatch(
      expr = { 
        print(paste0('Processing using method: ', names(calculate_VP_list)[mtd]))
        
        res <- calculate_VP_list[[mtd]](data)
        output_df <- output_df %>%
          full_join(res)
      },
      
      error = function(e){ # e represents a possible error
        err <- paste(Sys.time(), paste0('Error in analysis using method ', names(calculate_VP_list)[mtd], ':'), e)
        print(err)
        calc_vp_error_list[[length(calc_vp_error_list) + 1]] <<- err
      },
      
      warning = function(w){ # w represents a possible warning
        warn <- paste(Sys.time(), paste0('Warning in analysis using method ', names(calculate_VP_list)[mtd], ':'), w)
        print(warn)
        calc_vp_warn_list[[length(calc_vp_warn_list) + 1]] <<- warn
      },
      
      finally = {
        print(paste0('Analysis done for method ', mtd, '. Please check calc_vp_error_list and calc_vp_warn_list for any error or warnings.'))
      }
    )
  }
  
  # Arrange output dataframe
  output_df <- output_df %>%
    arrange('Location', 'Station_Number', 'Depth',
            factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
            factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')))
  
  # Write results in csv. If no csv is wanted, set write_csv to F
  if (write_csv == T){
    write.csv(output_df, file.path(res_path, 'vp_calc_ALL.csv'), row.names = F)
  }
  
  # Save only the values of the whole assay (T0_T24)
  output_24_df <- output_df %>%
    filter(Time_Range == 'T0_T24')
  
  write.csv(output_24_df, file.path(res_path, 'vp_calc_24.csv'), row.names = F)
  
  ## 3. Separate replicate treatment results
  # LM-2 and VPCL-1 have separate replicate treatment. In the end, the samples are averaged and difference is measured since we are interested in lytic and lysogenic production
  # If results of separate replicate treatment are not wished, set SR_calc to F
  if (SR_calc == T){
    # Create output dataframe
    output_df_SR <- data.frame(Location = character(), Station_Number = character(), Depth = character(),
                               Time_Range = character(), Population = character(), Sample_Type = character(),
                               Replicate = character(), VP = numeric(), VP_SE = numeric(), 
                               VP_R_Squared = numeric(), VP_Type = character())
    # Select correct methods
    method_names <- c('LM_2', 'VPCL_1')
    method_SR <- which(names(calculate_VP_list) %in% method_names)
    
    # Run methods
    for (mtd in method_SR){
      tryCatch(
        expr = { 
          print(paste0('Processing using method: ', names(calculate_VP_list)[mtd], ', only separate replicate results.'))
          
          res <- calculate_VP_list[[mtd]](data, avg = F)
          output_df_SR <- output_df_SR %>%
            full_join(res)
        },
        
        error = function(e){ # e represents a possible error
          err <- paste(Sys.time(), paste0('Error in analysis using method ', mtd, ':'), e)
          print(err)
          calc_vp_error_list[[length(calc_vp_error_list) + 1]] <<- err
        },
        
        warning = function(w){ # w represents a possible warning
          warn <- paste(Sys.time(), paste0('Warning in analysis using method ', mtd, ':'), w)
          print(warn)
          calc_vp_warn_list[[length(calc_vp_warn_list) + 1]] <<- warn
        },
        
        finally = {
          print(paste0('Analysis done for method ', mtd, '. Please check calc_vp_error_list and calc_vp_warn_list for any error or warnings.'))
        }
      )
    }
    
    # Arrange output dataframe
    output_df_SR <- output_df_SR %>%
      arrange('Location', 'Station_Number', 'Depth',
              factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
              factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')))
    
    # Write results in csv. If no csv is wanted, set write_csv to F
    if (write_csv == T){
      write.csv(output_df_SR, file.path(res_path, 'vp_calc_SR.csv'), row.names = F)
    } 
  }
  
  ## 4. Bacterial endpoint
  # To take the bacterial endpoint into account and stop the assay earlier, avoiding biased results for VP samples. If not wanted, set bp_endpoint to F
  if (bp_endpoint == T){
    # Create output dataframe
    output_df_bp <- data.frame(Location = character(), Station_Number = character(), Depth = character(),
                               Time_Range = character(), Population = character(), Sample_Type = character(),
                               VP = numeric(), VP_SE = numeric(), VP_R_Squared = numeric(), VP_Type = character())
    
    # Determine bacterial endpoint
    bp_df <- data %>%
      unite(c('Location', 'Station_Number', 'Depth'), col = 'tag', remove = F)
    endpoint_list <- list()
    
    for (combi_tag in unique(bp_df$tag)){
      bp_df1 <- bp_df %>%
        filter(tag == combi_tag)
      
      bp_range <- bacterial_endpoint(bp_df1)
      endpoint_list[[length(endpoint_list)+1]] <- c(combi_tag, bp_range)
    }
    
    # Select results from output_df based on the bp_range and add to output dataframe
    for (i in 1:length(unique(bp_df$tag))){
      res <- output_df %>%
        unite('tag', c('Location', 'Station_Number', 'Depth'), remove = F) %>%
        filter(tag == endpoint_list[[i]][1] & Time_Range == endpoint_list[[i]][2]) %>%
        select(-tag)
      
      output_df_bp <- output_df_bp %>%
        full_join(res)
    }
    
    # Arrange output dataframe
    output_df_bp <- output_df_bp %>%
      arrange('Location', 'Station_Number', 'Depth',
              factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
              factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')))
    
    # Write results in csv. If no csv is wanted, set write_csv to F
    if (write_csv == T){
      write.csv(output_df_bp, file.path(res_path, 'vp_calc_BP.csv'), row.names = F)
    } 
  }
  
  return(output_df)
}

# Running calc_VP function
print(names(calculate_VP_list)) # Order of different methods possible to calculate viral production
#vp_calc_NJ1 <- calc_VP(data, output_dir = 'vp_calc_NJ1')
vp_calc_NJ2020 <- calc_VP(data_all, output_dir = 'vp_calc_NJ2020')

## 4. Analyzing results
# Three different input dataframes: 
# 1. NJ2020: output of Step 1 => viral and bacterial abundances from flow cytometer
# 2. vp_calc_NJ2020: output of Step 2 => viral production calculations based on different methods => VP = [VLP per mL per h] (VLP = virus-like particles)
# 3. NJ2020_abundance: consists of the viral and bacterial abundances of the original sample (seawater, WW_sample)
analyze_vpres <- function(vpres, data, abundance){
  
  # 1. Correct viral production: VP_values are based of T0 which consists of 40% bacteria, ideally this should be 100% like in the original sample
  # Current values represent thus 40%, factorize so that they represent 100%. Also, with propagation of error, take the SE into account
  vpres_corrected <- vpres %>%
    unite(c('Location', 'Station_Number', 'Depth'), col = 'tag', remove = F) %>%
    mutate(c_VP = VP / 0.4363207,
           c_VP_SE = abs(c_VP) * (VP_SE / abs(VP))) %>%
    mutate_at(c('Station_Number'), as.integer) %>%
    filter(Sample_Type != 'VPC') %>% # Lytic fase = VP_samples; Lysogenic fase = Diff_Samples
    arrange('tag', 'Sample_Type', 'Population', 'Time_Range', 'VP_Type')
  
  # 2. Lytic and lysogenic viral production: c_VP corresponds to the amount of virus-like particles (VLP) per mL per hour in the given Time_Range
  # For each Time_Range, the percentage of infected cells will be calculated. VP_samples represent the lytic fase (% lytically infected cells), Diff_samples represent the lysogenic fase (% lysogenic cells)
  # Therefore, the burst size is needed: Burst size = the number of viral particles released from an infected host cell during the lytic cycle of viral replication => how many new viral particles will be released when cell bursts
  # A higher burst size will result in a lower percentage infected cells since more phages burst
  
  # Extract average bacterial abundance at T0 from original data
  B0_df <- df_AVG(data) %>% # Average data over replicates and add tag => necessary to get the bacterial abundance at T0 for each sample
    unite(c('Location', 'Station_Number', 'Depth'), col = 'tag', remove = F) %>%
    filter(Population == 'c_Bacteria', Timepoint == 0)
  
  # Add bacterial abundance at timepoint 0 as column B_T0 to dataframe
  perc_df <- vpres_corrected %>%
    left_join(select(B0_df, tag, Sample_Type, Mean), by = c('tag', 'Sample_Type')) %>%
    rename(B_T0 = Mean) %>%
    distinct() # Remove duplicate rows (created by left join)
  
  # Calculate percentage lytic and lysogenic cells for each burst size
  BS <- c(10,25,40) # Quick search online only returned burst sizes of laboratory environment => 25 as mid burst size and 15 below and above for range of three different burst sizes
  
  for (bs in BS){
    col_name <- paste0('%Cells_BS_', bs)
    perc_df[[col_name]] <- perc_df$VP * (100 / (perc_df$B_T0 * bs))
  }
  
  # 3. Lysis & Lysogenic rate of bacteria:
  # Lysis rate = rate at which bacterial cells rupture; Lysogenic rate = rate at which bacterial cells become lysogenized
  # First, calculate the lytic and lysogenic viral production in the original sample: VP * (B_OS / B_T0) with B_OS the bacterial abundance in the original sample
  
  rate_df <- perc_df %>%
    left_join(select(abundance, Station_Number, Total_Bacteria), by = c('Station_Number')) %>%
    rename(B_OS = Total_Bacteria) %>%
    mutate(c_VP_OS = c_VP * (B_OS / B_T0))
  
  # Calculate the rate for each burst size
  for (bs in BS){
    col_name <- paste0('Rate_BS_', bs)
    rate_df[[col_name]] <- rate_df$c_VP_OS / bs
  }
  
  # 4. Percentage of bacterial production lysed and bacterial loss per day: 
  # Bacterial production lysed = quantity of bacterial biomass or production that undergoes lysis due to viral infection
  # Bacterial loss per day = rate at which bacteria are removed due to viral lysis
  BSP_OS = 0.0027e6 # From DEMO VIPCAL; Bacterial secondary production in original sample
  
  bacterial_df <- rate_df %>%
    mutate('%BP_Lysed_BS_10' = ifelse(Sample_Type == 'Diff', NA, Rate_BS_10 / BSP_OS),
           '%BP_Lysed_BS_25' = ifelse(Sample_Type == 'Diff', NA,Rate_BS_25 / BSP_OS),
           '%BP_Lysed_BS_40' = ifelse(Sample_Type == 'Diff', NA,Rate_BS_40 / BSP_OS),
           '%B_Loss_BS_10' = ifelse(Sample_Type == 'Diff', NA,((Rate_BS_10 * 100) / B_OS) * 24),
           '%B_Loss_BS_25' = ifelse(Sample_Type == 'Diff', NA,((Rate_BS_25 * 100) / B_OS) * 24),
           '%B_Loss_BS_40' = ifelse(Sample_Type == 'Diff', NA,((Rate_BS_40 * 100) / B_OS) * 24))
  
  # 5. Viral turnover time: Time it takes for new viral particles to be produced
  # Just like the bacterial abundance in the original sample, need of adding the viral abundance of the original sample
  viral_abundance <- df_abundance %>%
    rename(c_Bacteria = Total_Bacteria, c_Viruses = Total_Viruses, c_V1 = V1, c_V2 = V2, c_V3 = V3) %>%
    pivot_longer(cols = c(7:13), names_to = 'Population', values_to = 'Abundance')
  
  viral_df <- bacterial_df %>%
    left_join(select(viral_abundance, Station_Number, Population, Abundance), by = c('Station_Number', 'Population')) %>%
    rename(V_OS = Abundance) %>%
    mutate(V_TT = c_VP_OS / V_OS)
  
  # 6. Nutrient release in bacteria and viruses for C, N and P: 
  # In literature, search for average content of these nutrients in bacteria and viruses
  
  # Data for bacteria from fjorden in Norway just behind North Sea (October 1993) => native aquatic sample => native bacteria
  # 1 fg = 10-15 g
  C_B <- 19e-15
  N_B <- 5e-15
  P_B <- 0.8e-15
  
  # Based on sequence information of bacteriophages => elemental composition of Synechococcus ssp (only marine virus in article)
  # Number of atoms are given in article => transform to amount of gram
  C_V <- (1664612 / 6.022e23) * 12.01 # ([atoms] / [atoms/mol]) * [g/mol]
  N_V <- (563058 / 6.022e23) * 14.01
  P_V <- (92428 / 6.022e23) * 30.97
  
  nutrient_df <- viral_df %>%
    mutate(B_DOC_BS_10 = Rate_BS_10 * C_B,
           B_DOC_BS_25 = Rate_BS_25 * C_B,
           B_DOC_BS_40 = Rate_BS_40 * C_B,
           B_DON_BS_10 = Rate_BS_10 * N_B,
           B_DON_BS_25 = Rate_BS_25 * N_B,
           B_DON_BS_40 = Rate_BS_40 * N_B,
           B_DOP_BS_10 = Rate_BS_10 * P_B,
           B_DOP_BS_25 = Rate_BS_25 * P_B,
           B_DOP_BS_40 = Rate_BS_40 * P_B,
           V_DOC = V_TT * C_V,
           V_DON = V_TT * N_V,
           V_DOP = V_TT * P_V)

  return(nutrient_df)
}

analyze_VP_NJ2020 <- analyze_vpres(vp_calc_NJ2020, data_all, df_abundance)





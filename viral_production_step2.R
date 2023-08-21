### Viral Production Assay: Step 2 ###

## 1. Importing used functions, these can be found in viral_production_step2_source.R
source("viral_production_step2_source.R")

## 2. Importing data from Step 1
data <- read.csv('NJ2020.csv') %>%
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
    .GlobalEnv$output_dir <- output_dir
    dir.create(output_dir)
  }
  
  # Storing all errors and warnings in global environment
  .GlobalEnv$calc_vp_error_list <- list()
  .GlobalEnv$calc_vp_warn_list <- list()
  
  ## 2. Calculate viral production
  # Create output dataframe
  res_path <- paste0('./', output_dir, '/') # output.dir needs to be in current working directory!
  output_df <- data.frame(Location = character(), Station_Number = integer(), Depth = character(),
                          Time_Range = character(), Population = character(), Sample_Type = character(),
                          VP = numeric(), abs_VP = numeric(), VP_SE = numeric(), VP_R_Squared = numeric(), VP_Type = character())
  
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
  
  if (write_csv == T){
    write.csv(output_24_df, file.path(res_path, 'vp_calc_24.csv'), row.names = F)
  }
  
  ## 3. Separate replicate treatment results
  # LM-2 and VPCL-1 have separate replicate treatment. In the end, the samples are averaged and difference is measured since we are interested in lytic and lysogenic production
  # If results of separate replicate treatment are not wished, set SR_calc to F
  if (SR_calc == T){
    # Create output dataframe
    output_df_SR <- data.frame(Location = character(), Station_Number = integer(), Depth = character(),
                               Time_Range = character(), Population = character(), Sample_Type = character(),
                               Replicate = character(), VP = numeric(), abs_VP = numeric(), VP_SE = numeric(), 
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
    output_df_bp <- data.frame(Location = character(), Station_Number = integer(), Depth = character(),
                               Time_Range = character(), Population = character(), Sample_Type = character(),
                               VP = numeric(), abs_VP = numeric(), VP_SE = numeric(), VP_R_Squared = numeric(), VP_Type = character())
    
    # Determine bacterial endpoint
    bp_df <- data %>%
      unite(all_of(c('Location', 'Station_Number', 'Depth')), col = 'tag', remove = F)
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
        unite('tag', all_of(c('Location', 'Station_Number', 'Depth')), remove = F) %>%
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
vp_calc_NJ2020 <- calc_VP(data, output_dir = 'vp_calc_NJ2020')

## 4. Analyzing results
# Three different input dataframes: 
# 1. data: output of Step 1 => viral and bacterial abundances from flow cytometer
# 2. vpres: output of Step 2 => viral production calculations based on different methods => VP = [VLP per mL per h] (VLP = virus-like particles)
# 3. abundance: consists of the viral and bacterial abundances of the original sample (seawater, WW_sample)
analyze_vpres <- function(vpres, data, abundance, BS = c(), BSP = NULL, nutrient_content_B = list(), 
                          nutrient_content_V = list(), write_output = T){
  
  ## Setup
  # Check if output folder of step 2 is there (Step 2 needs to be runned)
  if (!file.exists(output_dir)){
    print(paste0('The ', output_dir, ' folder does not exists!'))
    stop('Please run step 2 before proceeding.')
  }
  
  res_path <- paste0('./', output_dir, '/') # output.dir needs to be in current working directory!
  
  # Bacterial abundance at T0
  B0_df <- df_AVG(data) %>%
    filter(Timepoint == 0, Population == 'c_Bacteria', Sample_Type == 'VP') %>%
    select(all_of(c('Station_Number', 'Timepoint', 'Population', 'Sample_Type', 'Mean'))) %>%
    distinct() # Remove duplicate rows
  
  # Bacterial and viral abundance in the original sample
  OS_df <- abundance %>%
    select(all_of(c('Station_Number', 'Total_Bacteria', 'Total_Viruses')))
  
  # Join both to vpres dataframe
  vpres_df <- vpres %>%
    left_join(select(B0_df, Station_Number, Mean), by = c('Station_Number')) %>%
    rename(B_0 = Mean) %>%
    left_join(select(OS_df, Station_Number, Total_Bacteria, Total_Viruses), by = c('Station_Number')) %>%
    rename(B_OS = Total_Bacteria, V_OS = Total_Viruses)
  
  # Give default values if burst size and bacterial secondary production are not given
  if (is_empty(BS)){
    BS <- c(10,25,40) 
  }
  if (is.null(BSP)){
    BSP <- 0.0027e6 # Is in the same unit as the viral production values => VLP cells / mL h
  }
  if (is_empty(nutrient_content_B)){
    # Article: Fagerbakke, KM & Heldal, Mikal & Norland, S. (1996). Content of Carbon, Nitrogen, Oxygen, Sulfur and Phosphorus in Native Aquatic and Cultured Bacteria. Aquatic Microbial Ecology - AQUAT MICROB ECOL. 10. 15-27. 10.3354/ame010015
    nutrient_content_B <- list(C = 19e-15, N = 5e-15, V = 0.8e-15) # unit = g (C/N/P) / cell 
  }
  if (is_empty(nutrient_content_V)){
    # Article: Jover, L., Effler, T., Buchan, A. et al. The elemental composition of virus particles: implications for marine biogeochemical cycles. Nat Rev Microbiol 12, 519–528 (2014). https://doi.org/10.1038/nrmicro3289
    C_V <- (1664612 / 6.022e23) * 12.01 # ([atoms] / [atoms/mol]) * [g/mol]
    N_V <- (563058 / 6.022e23) * 14.01
    P_V <- (92428 / 6.022e23) * 30.97
    nutrient_content_V <- list(C = C_V, N = N_V, P = P_V) # unit = g (C/N/P) / cell 
  }
  
  ## 1. Correct values:
  # Viral production values are calculated based on T0_concentrate (after filtration), the amount of bacteria does not match the amount of the original sample (OS)
  # Correction of the viral production rate (VP) and absolute viral production (abs_VP) by multiplying by B_OS/B_0
  # With propagation of error, the standard error is also taken into account
  vpres_corrected <- vpres_df %>%
    mutate(c_VP = VP * (B_OS / B_0),
           c_abs_VP = abs_VP * (B_OS / B_0),
           c_VP_SE = abs(c_VP) * (VP_SE / abs(VP)))
  
  ## 2. Lytic and lysogenic viral production:
  # abs_VP corresponds to the absolute viral production at that point of the experiment, VP takes the time range in account and refers to the mean viral production or viral production rate in that time range
  # Given the burst size, the percentage of cells can be determined (VP_samples = percentage of lytically infected cells; Diff_samples = percentage of lysogenic cells; VPC_samples = both)
  # A higher burst size, meaning that the bacteriophages will burst, will result in a lower percentage cells
  
  ## 3. Lysis & Lysogenic rate of bacteria:
  # The rate at which bacterial cells rupture (lytic) or become lysogenized (lysogenic)
  # Need to use the viral production rate of the original sample (c_VP)
  
  ## 4. Percentage of bacterial production lysed and bacterial loss per day:
  # The quantity of bacterial biomass that undergoes lysis, depends on bacterial secondary production, and the rate at which bacteria are removed due to viral lysis
  
  for (bs in BS){
    ## 2.
    col_name <- paste0('%Cells_BS_', bs)
    vpres_corrected[[col_name]] <- vpres_corrected$abs_VP * (100 / (vpres_corrected$B_0 * bs))
    
    ## 3.
    col_name1 <- paste0('Rate_BS_', bs)
    vpres_corrected[[col_name1]] <- vpres_corrected$c_VP / bs
    
    ## 4.
    col_name2 <- paste0('%BP_Lysed_BS_', bs)
    vpres_corrected[[col_name2]] <- vpres_corrected[[col_name1]] / BSP
    
    col_name3 <- paste0('%B_Loss_BS_', bs)
    vpres_corrected[[col_name3]] <- ((vpres_corrected[[col_name1]] * 100) / vpres_corrected$B_OS) * 24
  }
  
  ## 5. Viral turnover time: 
  # Given x is the population of viruses right now, the viral turnover time describes the time it takes to replace x by new viruses => due to lysis new viral particles come free
  vpres_corrected <- vpres_corrected %>%
    mutate(V_TT = c_VP / V_OS)
  
  ## 6. Nutrient release in bacteria and viruses for C, N and P:
  for (bs in BS){
    for (nutrient in 1:length(nutrient_content_B)){
      current_rate_column <- paste0('Rate_BS_', bs)
      col_name <- paste0('DO', names(nutrient_content_B[nutrient]), '_B_BS_', bs)
      vpres_corrected[[col_name]] <- vpres_corrected[[current_rate_column]] * nutrient_content_B[[nutrient]]
    }
  }
  
  for (nutrient in 1:length(nutrient_content_V)){
    col_name <- paste0('DO', names(nutrient_content_V[nutrient]), '_V')
    vpres_corrected[[col_name]] <- vpres_corrected[[current_rate_column]] * nutrient_content_V[[nutrient]]
  }
  
  # Write results in csv. If no csv is wanted, set write_csv to F
  if (write_output == T){
    
    # Create data dictionary
    data_dictionary <- data.frame(
      Variable = colnames(vpres_corrected),
      Unit = c('/', '/', 'm', 'h', '/', '/', '#VLP (virus-like particles)/mLh', '#VLP/mL', '/', '/', '/',
               '#VLP/mL', '#VLP/mL', '#VLP/mL', '#VLP/mLh', '#VLP/mL', '/', '%', '#VLP/mLh', '%', '%',
               '%', '#VLP/mLh', '%', '%', '%', '#VLP/mLh', '%', '%', '1/h', 'g C/mLh', 'g N/mLh', 'g P/mLh',
               'g C/mLh', 'g N/mLh', 'g P/mLh', 'g C/mLh', 'g N/mLh', 'g P/mLh', 'g C/mLh', 'g N/mLh', 'g P/mLh'),
      Description = c(
        'Location of the experiment',
        'Number of station where experiment is conducted',
        'Depth at which experiment is performed',
        'Timepoints of sampling data, expressed in a time range: starting at T = 0 until time of measuring X => T0_TX',
        'Population Types: c_Viruses covers the entire virus population, while c_V1, c_V2, and c_V3 represent subpopulations',
        'Sample Types: VP rfor lytic viral production, Diff for lysogenic viral production, VPC for both',
        'Viral production rate: the mean viral production rate during the current Time_Range',
        'Absolute viral production at current timepoint',
        'The standard error on the viral production rate',
        'R-squared value: goodness of fit of the linear regression model',
        'Calculation method of viral production',
        'Bacterial abundance at the beginning of the experiment (T0)',
        'Bacterial abundance in the original sample',
        'Viral abundance in the original sample',
        'Corrected viral production rate: viral production rate in the original sample',
        'Corrected absolute viral production: absolute viral production in the original sample',
        'Corrected standard error on the viral production rate',
        'Percentage of cells for given burst size: % lytically infected cells for VP samples, % lysogenic cells for Diff samples',
        'Rate of bacteria for given burst size: lysis rate of bacteria for VP samples, lysogenic rate of bacteria for Diff samples',
        'Percentage of bacterial production lysed: the quantity of bacterial biomass that undergoes lysis',
        'Percentage of bacterial loss per day: the rate at which bacteria are removed due to viral lysis',
        'Percentage of cells for given burst size: % lytically infected cells for VP samples, % lysogenic cells for Diff samples',
        'Rate of bacteria for given burst size: lysis rate of bacteria for VP samples, lysogenic rate of bacteria for Diff samples',
        'Percentage of bacterial production lysed: the quantity of bacterial biomass that undergoes lysis',
        'Percentage of bacterial loss per day: the rate at which bacteria are removed due to viral lysis',
        'Percentage of cells for given burst size: % lytically infected cells for VP samples, % lysogenic cells for Diff samples',
        'Rate of bacteria for given burst size: lysis rate of bacteria for VP samples, lysogenic rate of bacteria for Diff samples',
        'Percentage of bacterial production lysed: the quantity of bacterial biomass that undergoes lysis',
        'Percentage of bacterial loss per day: the rate at which bacteria are removed due to viral lysis',
        'Viral turnover time: time to replacte the current virus population by new viruses',
        'Dissolved organic carbon release of bacteria for given burst size',
        'Dissolved organic nitrogen release of bacteria for given burst size',
        'Dissolved organic phosphorous release of bacteria for given burst size',
        'Dissolved organic carbon release of bacteria for given burst size',
        'Dissolved organic nitrogen release of bacteria for given burst size',
        'Dissolved organic phosphorous release of bacteria for given burst size',
        'Dissolved organic carbon release of bacteria for given burst size',
        'Dissolved organic nitrogen release of bacteria for given burst size',
        'Dissolved organic phosphorous release of bacteria for given burst size',
        'Dissolved organic carbon release of viruses',
        'Dissolved organic nitrogen release of viruses',
        'Dissolved organic phosphorous release of viruses'
      )
    )
    
    # Create Workbook
    wb <- createWorkbook()
    
    # Add data and dictionary as separate sheets
    addWorksheet(wb, 'Data')
    writeData(wb, sheet = 'Data', vpres_corrected)
    
    addWorksheet(wb, 'Data dictionary')
    writeData(wb, sheet = 'Data dictionary', data_dictionary)
    
    # Set column width and bold column names
    setColWidths(wb, sheet = "Data dictionary", cols = 1:ncol(data_dictionary), widths = "auto")
    
    bold_style <- createStyle(textDecoration = c("bold"))
    addStyle(wb, "Data", style = bold_style, rows = 1, cols = 1:ncol(vpres_corrected))
    addStyle(wb, "Data dictionary", style = bold_style, rows = 1, cols = 1:ncol(data_dictionary))
    
    # Save workbook as excel file
    saveWorkbook(wb, file = paste0(res_path, 'vp_calc_ANALYZED.xlsx'), overwrite = T)
  }

  return(vpres_corrected)
}

analyze_VP_NJ2020 <- analyze_vpres(vp_calc_NJ2020, data, df_abundance)





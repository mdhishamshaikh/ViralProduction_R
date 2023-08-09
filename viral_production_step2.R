### Viral Production Assay: Step 2 ###

## 1. Importing used functions, these can be found in viral_production_step2_source.R
source("viral_production_step2_source.R")

## 2. Importing data from Step 1
data <- read.csv('NJ1.csv')
names(data)[names(data) == 'Expt_No'] <- 'Station_Number' # Changing column name to something more appropriate => you can also use index of column (5)

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
  
  # Write results in csv. If no csv is wanted, set write_csv to F
  if (write_csv == T){
    write.csv(output_df, file.path(res_path, 'vp_calc_ALL.csv'), row.names = F)
  }
  
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
    
    # Write results in csv. If no csv is wanted, set write_csv to F
    if (write_csv == T){
      write.csv(output_df_bp, file.path(res_path, 'vp_calc_BP.csv'), row.names = F)
    } 
  }
  
  return(output_df)
}

# Running calc_VP function
print(names(calculate_VP_list)) # Order of different methods possible to calculate viral production
vp_calc_NJ1 <- calc_VP(data, output_dir = 'vp_calc_NJ1')

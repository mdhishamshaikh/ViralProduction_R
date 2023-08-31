calculate_viral_production <- function(data,
                                       methods = c(1:12),
                                       SR_calc = TRUE,
                                       BP_endpoint = TRUE,
                                       write_csv = TRUE,
                                       output_dir = ''){
  ## 1. Checks
  # Check for valid output directory if csv files needs to be written
  if (write_csv == T){
    if (output_dir == ''){
      stop('No output directory is given, please define output directory before proceeding!')
      
    } else if (file.exists(output_dir)){
      stop('The output folder already exists, please define another output directory before proceeding!')
      
    } else{
      .GlobalEnv$output_dir <- output_dir
      dir.create(output_dir)
    }
  }
  
  ## 2. Setup
  # Store all errors and warnings in list
  .GlobalEnv$calc_vp_error_list <- list()
  .GlobalEnv$calc_vp_warn_list <- list()
  
  # Import all the possible methods and check the populations to analyze
  vp_list_of_methods()
  vp_check_populations(data)
  
  ## 3. Calculate viral production
  vp_results_path <- paste0(output_dir, '/')
  vp_results_output_df <- data.frame(tag = character(), Location = character(), Station_Number = integer(), 
                                     Depth = integer(), Time_Range = character(), Population = character(), 
                                     Sample_Type = character(), VP = numeric(), abs_VP = numeric(), VP_SE = numeric(), 
                                     VP_R_Squared = numeric(), VP_Type = character())
  
  for (method in methods){
    tryCatch(
      expr = {
        print(paste0('Processing using method: ', names(.GlobalEnv$list_of_methods)[method]))
        vp_results <- .GlobalEnv$list_of_methods[[method]](data)
        
        vp_results_output_df <- dplyr::full_join(vp_results_output_df, vp_results)
      },
      
      error = function(e){
        err <- paste(Sys.time(), paste0('Error in analysis using method ', names(.GlobalEnv$list_of_methods)[method], ':'), e)
        .GlobalEnv$calc_vp_error_list[[length(.GlobalEnv$calc_vp_error_list) + 1]] <<- err
        print(err)
      },
      
      warning = function(w){
        warn <- paste(Sys.time(), paste0('Warning in analysis using method ', names(.GlobalEnv$list_of_methods)[method], ':'), w)
        .GlobalEnv$calc_vp_warn_list[[length(.GlobalEnv$calc_vp_warn_list) + 1]] <<- warn
        print(warn)
      },
      
      finally = {
        print(paste0('Analysis done for method ', names(.GlobalEnv$list_of_methods)[method], '. Please check calc_vp_error_list and calc_vp_warn_list for any error or warnings.'))
      }
    )
  }
  
  vp_results_output_df <- vp_results_output_df %>%
    dplyr::arrange('tag',
                   factor(.data$Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                   factor(.data$Population, levels = .GlobalEnv$populations_to_analyze[which(startsWith(.GlobalEnv$populations_to_analyze, 'c_V'))])) %>%
    dplyr::select(-.data$tag)
  
  vp_results_output_T24_df <- vp_results_output_df %>%
    dplyr::filter(.data$Time_Range == 'T0_T24')
  
  if (write_csv == T){
    utils::write.csv(vp_results_output_df, file.path(vp_results_path, 'vp_results_ALL.csv'), row.names = F)
    utils::write.csv(vp_results_output_T24_df, file.path(vp_results_path, 'vp_results_24.csv'), row.names = F)
  }
  
  ## 4. Results of separate replicate treatment
  
  ## 5. Results with bacterial endpoint taken into account
  
  
  
  return(vp_results_output_df)
}
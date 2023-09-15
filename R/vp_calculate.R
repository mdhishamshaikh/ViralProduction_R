#' Wrapper function to calculate viral production
#' 
#' @description
#' Wrapper function that performs viral production calculation. Given the output of the flow cytometry step, viral production
#' is calculated for all methods of linear regression and VIPCAL. Based on the parameters, different data frames will be available
#' in the global environment and written as csv files.
#' 
#' 1. `vp_results_ALL.csv`: Contains the viral production results for all samples.
#' 2. `vp_results_T24.csv`: Contains the viral production results for all samples at the end of the assay (T0_T24).
#' 3. `vp_results_SR.csv`: Contains the viral production results for all samples of the separate replicate treatment, VP and VPC samples only.
#' 4. `vp_results_BP.csv`: Contains the viral production results for all samples with the bacterial endpoint into account. Results for that time range will be stored in data frame.
#' 
#' More details about the calculation methods:
#' 
#' - Linear regression variants: [viralprod::vp_methods_linear]
#' - VIPCAL variants: [viralprod::vp_methods_VIPCAL]
#'
#' @param x Data frame with the output of the flow cytometry, has to have the `viralprod` class.
#' @param ... Arguments passed on to the next function.
#' @param methods Integer vector with the indexes of `list_of_methods`. Indexes determine which methods will be performed within the viral production calculation. 
#' @param SR_calc If \code{TRUE}, separate replicate treatment results will be stored and saved as an separate data frame. 
#' Set to \code{FALSE}, if separate replicate treatment results are not wanted. (Default = \code{TRUE})
#' @param BP_endpoint If \code{TRUE}, the bacterial endpoint will be taken into account and only those results will be saved in a new
#' data frame. If not wanted, set to \code{FALSE}. (Default = \code{TRUE})
#' @param write_output If \code{TRUE}, the output data frames will be saved as csv files in a folder specified by `output_dir`.
#' If no csv files are wanted, set to \code{FALSE}. (Default = \code{TRUE})
#' @param output_dir String that refers to the location of folder to save the data frames as csv files.
#'
#' @return Depending on setting of parameters, different data frames with the viral production calculation will be available in the global environment.
#' @export
#' 
#' @name vp_calculate
#' @rdname vp_calculate
#'
#' @examples \dontrun{
#' # Setup
#' data_NJ2020_all <- read.csv(system.file('extdata', 
#' 'NJ2020_Station_2_and_6_all_populations.csv', package = "viralprod"))
#' 
#' # Perform
#' # Add S3 class
#' x <- vp_class_count_data(data_NJ2020_all)
#' 
#' # Default method
#' vp_calculate(x = data_NJ2020_all, write_output = F)
#' 
#' # S3 class, viralprod, method
#' # Write output files
#' vp_calculate(x, 
#' output_dir = paste0(system.file(“extdata”, package = “viralprod”), 
#' “/NJ2020_vp_results”))
#' 
#' # No output files
#' vp_calculate(x, write_output = F)
#' 
#' # No bacterial endpoint and separate replicate treatment results
#' vp_calculate(x, write_output = F, SR_calc = F, BP_endpoint = F)
#' 
#' # Sub selection of the methods
#' vp_calculate(x, write_output = F, methods = c(2,3,8,12))
#' }
vp_calculate <- function(x, ...){
  UseMethod("vp_calculate")
}


#' @export
#' @rdname vp_calculate
vp_calculate.default <- function(x, ...){
  stop('The input data frame has not the correct S3 class! You can check your data frame with `vp_class_count_data()`.')
}


#' @export
#' @rdname vp_calculate
vp_calculate.viralprod <- function(x, ...,
                                   methods = c(1:12),
                                   SR_calc = TRUE,
                                   BP_endpoint = TRUE,
                                   write_output = TRUE,
                                   output_dir = ''){
  ## 1. Checks
  # Check for valid output directory if csv files needs to be written
  if (write_output == T){
    if (output_dir == ''){
      stop('No output directory is given, please define output directory before proceeding!')
      
    } else if (file.exists(output_dir)){
      # Give the user an option if he wants to continue and possibly overwrite results
      user_choice <- utils::menu(
        c("Continue and possibly overwrite existing results",
          "Stop and define another output directory"),
        title = "Warning: The output folder already exists!",
        graphics = FALSE)
      
      if (user_choice == 1){
        print('Continuing with the existing output directory, existing results will possibly be overwritten!')
        vp_results_path <- paste0(output_dir, '/')
        
      } else {
        stop('Process stopped, please define another output directory!')
      }
    } else{
      dir.create(output_dir)
      vp_results_path <- paste0(output_dir, '/')
    }
  }
  
  # Check if data frame isn't empty
  if (rlang::is_empty(x)){
    stop('Not able to proceed since input data frames is empty!')
  }
  
  ## 2. Setup
  # Store all errors and warnings in list
  calc_vp_error_list <- list()
  calc_vp_warn_list <- list()
  
  # Import all the possible methods and check the populations to analyze
  vp_list_of_methods()
  vp_check_populations(x)
  
  ## 3. Calculate viral production
  vp_results_output_df <- data.frame(tag = character(), Location = character(), Station_Number = integer(), 
                                     Depth = integer(), Time_Range = character(), Population = character(), 
                                     Sample_Type = character(), VP = numeric(), abs_VP = numeric(), 
                                     VP_SE = numeric(), VP_R_Squared = numeric(), VP_Method = character())
  
  for (method in methods){
    tryCatch(
      expr = {
        print(paste0('Processing using method: ', names(.GlobalEnv$list_of_methods)[method]))
        vp_results <- .GlobalEnv$list_of_methods[[method]](x)
        
        vp_results_output_df <- dplyr::full_join(vp_results_output_df, vp_results)
      },
      
      error = function(e){
        err <- paste(Sys.time(), paste0('Error in analysis using method ', names(.GlobalEnv$list_of_methods)[method], ':'), e)
        calc_vp_error_list[[length(calc_vp_error_list) + 1]] <<- err
        print(err)
      },
      
      warning = function(w){
        warn <- paste(Sys.time(), paste0('Warning in analysis using method ', names(.GlobalEnv$list_of_methods)[method], ':'), w)
        calc_vp_warn_list[[length(calc_vp_warn_list) + 1]] <<- warn
        print(warn)
      },
      
      finally = {
        print(paste0('Analysis done for method ', names(.GlobalEnv$list_of_methods)[method], '. Please check calc_vp_error_list and calc_vp_warn_list for any error or warnings.'))
      }
    )
  }
  .GlobalEnv$vp_results_output_df <- vp_results_output_df %>%
    dplyr::arrange('tag',
                   factor(.data$Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                   factor(.data$Population, levels = .GlobalEnv$populations_to_analyze[grep('c_V', .GlobalEnv$populations_to_analyze)])) %>%
    dplyr::select(-'tag')
  
  .GlobalEnv$vp_results_output_T24_df <- .GlobalEnv$vp_results_output_df %>%
    dplyr::filter(.data$Time_Range == 'T0_T24')
  
  if (write_output == T){
    utils::write.csv(.GlobalEnv$vp_results_output_df, file.path(vp_results_path, 'vp_results_ALL.csv'), row.names = F)
    utils::write.csv(.GlobalEnv$vp_results_output_T24_df, file.path(vp_results_path, 'vp_results_24.csv'), row.names = F)
  }
  
  ## 4. Results of separate replicate treatment
  if (SR_calc == T){
    vp_results_output_SR_df <- data.frame(tag = character(), Location = character(), Station_Number = integer(), 
                                          Depth = integer(), Time_Range = character(), Population = character(), 
                                          Sample_Type = character(), VP = numeric(), abs_VP = numeric(), VP_SE = numeric(), 
                                          VP_R_Squared = numeric(), VP_Method = character())
    
    methods_SR <- grep("separate_replicates", .GlobalEnv$list_of_methods)
    
    for (method in methods_SR){
      tryCatch(
        expr = {
          print(paste0('Processing using method: ', names(.GlobalEnv$list_of_methods)[method], ', only separate replicate results'))
          vp_results <- .GlobalEnv$list_of_methods[[method]](x, AVG = F)
          
          vp_results_output_SR_df <- dplyr::full_join(vp_results_output_SR_df, vp_results)
        },
        
        error = function(e){
          err <- paste(Sys.time(), paste0('Error in analysis using method ', names(.GlobalEnv$list_of_methods)[method], ':'), e)
          calc_vp_error_list[[length(calc_vp_error_list) + 1]] <<- err
          print(err)
        },
        
        warning = function(w){
          warn <- paste(Sys.time(), paste0('Warning in analysis using method ', names(.GlobalEnv$list_of_methods)[method], ':'), w)
          calc_vp_warn_list[[length(calc_vp_warn_list) + 1]] <<- warn
          print(warn)
        },
        
        finally = {
          print(paste0('Analysis done for method ', names(.GlobalEnv$list_of_methods)[method], '. Please check calc_vp_error_list and calc_vp_warn_list for any error or warnings.'))
        }
      )
    }
    .GlobalEnv$calc_vp_error_list <- calc_vp_error_list
    .GlobalEnv$calc_vp_warn_list <- calc_vp_warn_list
    
    .GlobalEnv$vp_results_output_SR_df <- vp_results_output_SR_df %>%
      dplyr::arrange('tag',
                     factor(.data$Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                     factor(.data$Population, levels = .GlobalEnv$populations_to_analyze[grep('c_V', .GlobalEnv$populations_to_analyze)])) %>%
      dplyr::select(-'tag')
    
    if (write_output == T){
      utils::write.csv(.GlobalEnv$vp_results_output_SR_df, file.path(vp_results_path, 'vp_results_SR.csv'), row.names = F)
    }
  }
  
  ## 5. Results with bacterial endpoint taken into account
  if (BP_endpoint == T){
    vp_results_output_BP_df <- data.frame(tag = character(), Location = character(), Station_Number = integer(), 
                                          Depth = integer(), Time_Range = character(), Population = character(), 
                                          Sample_Type = character(), VP = numeric(), abs_VP = numeric(), VP_SE = numeric(), 
                                          VP_R_Squared = numeric(), VP_Method = character())
    
    BP_result_list <- list()
    data_with_tag <- x %>%
      tidyr::unite(dplyr::all_of(c('Location', 'Station_Number', 'Depth')), col = 'tag', remove = F)
    
    for (combi_tag in unique(data_with_tag$tag)){
      BP_DF <- data_with_tag %>%
        dplyr::filter(.data$tag == combi_tag)
      
      BP_res <- vp_bacterial_endpoint(BP_DF)
      BP_result_list[[length(BP_result_list) + 1]] <- c(combi_tag, BP_res)
    }
    
    for (index in 1:length(unique(data_with_tag$tag))){
      vp_results_output_BP_index <- .GlobalEnv$vp_results_output_df %>%
        tidyr::unite('tag', dplyr::all_of(c('Location', 'Station_Number', 'Depth')), remove = F) %>%
        dplyr::filter(.data$tag == BP_result_list[[index]][1] & .data$Time_Range == BP_result_list[[index]][2])
      
      vp_results_output_BP_df <- dplyr::full_join(vp_results_output_BP_df, vp_results_output_BP_index)
    }
    
    .GlobalEnv$vp_results_output_BP_df <- vp_results_output_BP_df %>%
      dplyr::arrange('tag',
                     factor(.data$Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                     factor(.data$Population, levels = .GlobalEnv$populations_to_analyze[grep('c_V', .GlobalEnv$populations_to_analyze)])) %>%
      dplyr::select(-'tag')
    
    if (write_output == T){
      utils::write.csv(.GlobalEnv$vp_results_output_BP_df, file.path(vp_results_path, 'vp_results_BP.csv'), row.names = F)
    }
  }
}

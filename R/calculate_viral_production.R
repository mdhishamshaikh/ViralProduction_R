#' Wrapper function that calculates viral production
#' 
#' @description
#' Wrapper function that does everything at once. Given the output of the flow cytometry step, viral production
#' is calculated for all the methods of linear regression and VIPCAL. The output data frames are written as csv files,
#' four different csv files are saved:
#' 
#' 1. `vp_results_ALL.csv`: contains the results of all the samples for each time range.
#' 2. `vp_results_T24.csv`: contains the results of all the samples for the end of the assay (T0_T24).
#' 3. `vp_results_SR.csv`: contains the results of the separate replicate treatment of all the samples for each time range.
#' 4. `vp_results_BP.csv`: contains the results of all the samples but takes the bacterial endpoint into account, bacterial endpoint decides time range of the results.
#' 
#' More details about each of the methods used:
#' 
#' - Linear regression variants: [viralprod::vp_methods_linear]
#' - VIPCAL variants: [viralprod::vp_methods_VIPCAL]
#'
#' @param data Data frame with the output of the flow cytometry.
#' @param methods Integer vector with the indexes of `list_of_methods`. Indexes determine which methods will be used for viral production calculation. 
#' @param SR_calc If \code{TRUE}, separate replicate treatment results will be performed and saved as an separate data frame. 
#' Set to \code{FALSE}, if separate replicate treatment results are not wanted. (Default = \code{TRUE})
#' @param BP_endpoint If \code{TRUE}, the bacterial endpoint will be taken into account and only those results will be saved in a new
#' data frame. If not wanted, set to \code{FALSE}. (Default = \code{TRUE})
#' @param write_csv If \code{TRUE}, the output data frames will be saved as csv files in a folder specified by `output_dir`.
#' If no csv files are wanted, set to \code{FALSE}. (Default = \code{TRUE})
#' @param output_dir String that refers to the location of folder to save the data frames as csv files.
#'
#' @return Depending on setting of parameters, different data frames with the viral production calculation will be available in the global environment.
#' @export
#' 
#' @name calculate_viral_production
#' @rdname calculate_viral_production
#'
#' @examples \dontrun{
#' data_NJ2020_all <- read.csv(system.file('extdata', 
#' 'NJ2020_Station_2_and_6_all_populations.csv', package = "viralprod"))
#' 
#' calculate_viral_production(data_NJ2020_all, 
#' output_dir = paste0(system.file(“extdata”, package = “viralprod”), 
#' “/vp_results_NJ2020”))
#' 
#' calculate_viral_production(data_NJ2020_all, write_csv = F)
#' 
#' calculate_viral_production(data_NJ2020_all, write_csv = F, SR_calc = F, BP_endpoint = F)
#' 
#' calculate_viral_production(data_NJ2020_all, write_csv = F, methods = c(2,3,8,12))
#' }
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
  .GlobalEnv$vp_results_output_df <- data.frame(tag = character(), Location = character(), Station_Number = integer(), 
                                                Depth = integer(), Time_Range = character(), Population = character(), 
                                                Sample_Type = character(), VP = numeric(), abs_VP = numeric(), VP_SE = numeric(), 
                                                VP_R_Squared = numeric(), VP_Method = character())
  
  for (method in methods){
    tryCatch(
      expr = {
        print(paste0('Processing using method: ', names(.GlobalEnv$list_of_methods)[method]))
        vp_results <- .GlobalEnv$list_of_methods[[method]](data)
        
        .GlobalEnv$vp_results_output_df <- dplyr::full_join(.GlobalEnv$vp_results_output_df, vp_results)
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
  
  .GlobalEnv$vp_results_output_df <- .GlobalEnv$vp_results_output_df %>%
    dplyr::arrange('tag',
                   factor(.data$Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                   factor(.data$Population, levels = .GlobalEnv$populations_to_analyze[grep('c_V', .GlobalEnv$populations_to_analyze)])) %>%
    dplyr::select(-'tag')
  
  .GlobalEnv$vp_results_output_T24_df <- .GlobalEnv$vp_results_output_df %>%
    dplyr::filter(.data$Time_Range == 'T0_T24')
  
  if (write_csv == T){
    utils::write.csv(.GlobalEnv$vp_results_output_df, file.path(vp_results_path, 'vp_results_ALL.csv'), row.names = F)
    utils::write.csv(.GlobalEnv$vp_results_output_T24_df, file.path(vp_results_path, 'vp_results_24.csv'), row.names = F)
  }
  
  ## 4. Results of separate replicate treatment
  if (SR_calc == T){
    .GlobalEnv$vp_results_output_SR_df <- data.frame(tag = character(), Location = character(), Station_Number = integer(), 
                                                     Depth = integer(), Time_Range = character(), Population = character(), 
                                                     Sample_Type = character(), VP = numeric(), abs_VP = numeric(), VP_SE = numeric(), 
                                                     VP_R_Squared = numeric(), VP_Method = character())
    
    methods_SR <- grep("separate_replicates", .GlobalEnv$list_of_methods)
    
    for (method in methods_SR){
      tryCatch(
        expr = {
          print(paste0('Processing using method: ', names(.GlobalEnv$list_of_methods)[method], ', only separate replicate results'))
          vp_results <- .GlobalEnv$list_of_methods[[method]](data, AVG = F)
          
          .GlobalEnv$vp_results_output_SR_df <- dplyr::full_join(.GlobalEnv$vp_results_output_SR_df, vp_results)
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
    
    .GlobalEnv$vp_results_output_SR_df <- .GlobalEnv$vp_results_output_SR_df %>%
      dplyr::arrange('tag',
                     factor(.data$Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                     factor(.data$Population, levels = .GlobalEnv$populations_to_analyze[grep('c_V', .GlobalEnv$populations_to_analyze)])) %>%
      dplyr::select(-'tag')
    
    if (write_csv == T){
      utils::write.csv(.GlobalEnv$vp_results_output_SR_df, file.path(vp_results_path, 'vp_results_SR.csv'), row.names = F)
    }
  }
  
  ## 5. Results with bacterial endpoint taken into account
  if (BP_endpoint == T){
    .GlobalEnv$vp_results_output_BP_df <- data.frame(tag = character(), Location = character(), Station_Number = integer(), 
                                                     Depth = integer(), Time_Range = character(), Population = character(), 
                                                     Sample_Type = character(), VP = numeric(), abs_VP = numeric(), VP_SE = numeric(), 
                                                     VP_R_Squared = numeric(), VP_Method = character())
    
    BP_result_list <- list()
    data_with_tag <- data %>%
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
      
      .GlobalEnv$vp_results_output_BP_df <- dplyr::full_join(.GlobalEnv$vp_results_output_BP_df, vp_results_output_BP_index)
    }
    
    .GlobalEnv$vp_results_output_BP_df <- .GlobalEnv$vp_results_output_BP_df %>%
      dplyr::arrange('tag',
                     factor(.data$Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                     factor(.data$Population, levels = .GlobalEnv$populations_to_analyze[grep('c_V', .GlobalEnv$populations_to_analyze)])) %>%
      dplyr::select(-'tag')
    
    if (write_csv == T){
      utils::write.csv(.GlobalEnv$vp_results_output_BP_df, file.path(vp_results_path, 'vp_results_BP.csv'), row.names = F)
    }
  }
}

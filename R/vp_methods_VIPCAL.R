#' VIPCAL variants for viral production calculation
#' 
#' @description
#' Two main methods are used to determine the viral production rate in the viral reduction assay: `Linear Regression` vs
#' `VIPCAL`. Different variants of both methods are performed: variance in replicate treatment, use of standard error
#' and type of estimation of the difference curve is presented. Next, the different variants of `VIPCAL`
#' which examines the average of increments as a benchmark for viral production. 
#' The following general step-by-step plan is carried out by each of the variants:
#'    
#' 1. Constructing correct data frame depending on replicate treatment.
#' 2. Calculate viral production depending on replicate treatment, standard error use and difference curve estimation.
#' If no difference curve estimation, calculate lysogenic viral production.
#' 3. Arrange output data frame.
#' 
#' `vp_VIPCAL_separate_replicates`: separate replicate treatment, no SE taken into account, no difference curve estimation.
#' 
#' `vp_VIPCAL_average_replicates`: average replicate treatment, no SE taken into account, no difference curve estimation.
#' 
#' `vp_VIPCAL_average_replicates_SE`: average replicate treatment, SE taken into account, no difference curve estimation.
#' 
#' `vp_VIPCAL_average_replicates_diff`: average replicate treatment, no SE taken into account, difference curve estimation by subtraction.
#' 
#' `vp_VIPCAL_average_replicates_diff_SE`: average replicate treatment, SE taken into account, difference curve estimation by subtraction.
#' 
#' `vp_VIPCAL_average_replicates_diff_LMER`: average replicate treatment, no SE taken into account, difference curve estimation by LMER model.
#' 
#' `vp_VIPCAL_average_replicates_diff_LMER_SE`: average replicate treatment, SE taken into account, difference curve estimation by LMER model.
#' 
#' More details about the used functions:
#' 
#' - Constructing data frames: [viralprod::vp_dataframes], [viralprod::vp_check_populations] 
#' - Calculating viral production with VIPCAL: [viralprod::determine_vp_VIPCAL]
#' - LMER model: [viralprod::vp_LMER_model]
#' - Calculating lysogenic production: [viralprod::vp_calculate_difference_samples]
#' 
#' @param data Data frame with the output of the flow cytometry.
#' @param AVG One of the methods uses an separate replicate treatment. In the end, interested in the lytic and lysogenic
#' viral production. For the last one, difference samples are needed. If \code{TRUE}, average over replicates after determining viral production with
#' separate replicate treatment and lysogenic viral production is calculated. If \code{FALSE}, no averaging over replicates,
#' output will consists of only VP and VPC samples with separate replicate treatment. (Default = \code{TRUE})
#'
#' @return Data frame with the viral production rate and the absolute viral production for each population at the different time points of the assay.
#' 
#' @name vp_methods_VIPCAL
#' @rdname vp_methods_VIPCAL
#'
#' @examples \dontrun{
#' data_NJ2020_all <- read.csv(system.file('extdata', 
#' 'NJ2020_Station_2_and_6_all_populations.csv', package = "viralprod"))
#' 
#' vp_check_populations(data_NJ2020_all)
#' 
#' vp_VIPCAL_separate_replicates(data_NJ2020_all)
#' vp_VIPCAL_separate_replicates(data_NJ2020_all, AVG = F)
#' 
#' vp_VIPCAL_average_replicates(data_NJ2020_all)
#' 
#' vp_VIPCAL_average_replicates_SE(data_NJ2020_all)
#' 
#' vp_VIPCAL_average_replicates_diff(data_NJ2020_all)
#' 
#' vp_VIPCAL_average_replicates_diff_SE(data_NJ2020_all)
#' 
#' vp_VIPCAL_average_replicates_diff_LMER(data_NJ2020_all)
#' 
#' vp_VIPCAL_average_replicates_diff_LMER_SE(data_NJ2020_all)
#' }
vp_VIPCAL_separate_replicates <- function(data,
                                          AVG = T){
  separate_replicate_dataframe_with_timepoints <- vp_separate_replicate_dataframe(data)
  
  determine_viral_production_dataframe <- determine_vp_VIPCAL_separate_replicates(separate_replicate_dataframe_with_timepoints)
  
  if (AVG == F){
    viral_production_VIPCAL <- determine_viral_production_dataframe %>%
      dplyr::group_by(.data$tag, .data$Time_Range, .data$Population, .data$Sample_Type) %>%
      dplyr::arrange('tag',
                     factor(.data$Sample_Type, levels = c('VP', 'VPC')),
                     factor(.data$Population, levels = .GlobalEnv$populations_to_analyze[grep('c_V', .GlobalEnv$populations_to_analyze)])) %>%
      dplyr::mutate(VP_Method = 'VPCL_SR') %>%
      dplyr::select('tag', 'Location', 'Station_Number', 'Depth', 'Replicate', dplyr::everything())
  } else {
    viral_production_VIPCAL <- determine_viral_production_dataframe %>%
      dplyr::group_by(.data$tag, .data$Location, .data$Station_Number, .data$Depth, .data$Time_Range, .data$Population, .data$Sample_Type) %>%
      dplyr::summarise(
        VP_mean = mean(.data$VP),
        abs_VP_mean = mean(.data$abs_VP),
        VP_SE = plotrix::std.error(.data$VP)) %>%
      dplyr::rename(VP = "VP_mean", abs_VP = "abs_VP_mean") %>%
      vp_calculate_difference_samples(VIPCAL = T) %>%
      dplyr::arrange('tag',
                     factor(.data$Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                     factor(.data$Population, levels = .GlobalEnv$populations_to_analyze[grep('c_V', .GlobalEnv$populations_to_analyze)])) %>%
      dplyr::mutate(VP_Method = 'VPCL_SR_AVG') %>%
      dplyr::select('tag', 'Location', 'Station_Number', 'Depth', dplyr::everything())
  }
  
  return(viral_production_VIPCAL)
}


#' @rdname vp_methods_VIPCAL
vp_VIPCAL_average_replicates <- function(data){
  average_replicate_dataframe_with_timepoints <- vp_average_replicate_dataframe(data)
  
  average_replicate_dataframe_no_diff_curve <- average_replicate_dataframe_with_timepoints %>%
    subset(average_replicate_dataframe_with_timepoints$Sample_Type != 'Diff')
  
  determine_viral_production_dataframe <- determine_vp_VIPCAL_average_replicates(average_replicate_dataframe_no_diff_curve) %>%
    vp_calculate_difference_samples(VIPCAL = T)
  
  viral_production_VIPCAL <- determine_viral_production_dataframe %>%
    dplyr::group_by(.data$tag, .data$Time_Range, .data$Population, .data$Sample_Type) %>%
    dplyr::arrange('tag',
                   factor(.data$Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                   factor(.data$Population, levels = .GlobalEnv$populations_to_analyze[grep('c_V', .GlobalEnv$populations_to_analyze)])) %>%
    dplyr::mutate(VP_Method = 'VPCL_AR') %>%
    dplyr::select('tag', 'Location', 'Station_Number', 'Depth', dplyr::everything())
  
  return(viral_production_VIPCAL)
}


#' @rdname vp_methods_VIPCAL
vp_VIPCAL_average_replicates_SE <- function(data){
  average_replicate_dataframe_with_timepoints <- vp_average_replicate_dataframe(data)
  
  average_replicate_dataframe_no_diff_curve <- average_replicate_dataframe_with_timepoints %>%
    subset(average_replicate_dataframe_with_timepoints$Sample_Type != 'Diff')
  
  determine_viral_production_dataframe <- determine_vp_VIPCAL_average_replicates_SE(average_replicate_dataframe_no_diff_curve) %>%
    vp_calculate_difference_samples(VIPCAL = T, SE = T)
  
  viral_production_VIPCAL <- determine_viral_production_dataframe %>%
    dplyr::group_by(.data$tag, .data$Time_Range, .data$Population, .data$Sample_Type) %>%
    dplyr::arrange('tag',
                   factor(.data$Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                   factor(.data$Population, levels = .GlobalEnv$populations_to_analyze[grep('c_V', .GlobalEnv$populations_to_analyze)])) %>%
    dplyr::mutate(VP_Method = 'VPCL_AR_SE') %>%
    dplyr::select('tag', 'Location', 'Station_Number', 'Depth', dplyr::everything())
  
  return(viral_production_VIPCAL)
}


#' @rdname vp_methods_VIPCAL
vp_VIPCAL_average_replicates_diff <- function(data){
  average_replicate_dataframe_with_timepoints <- vp_average_replicate_dataframe(data)
  
  determine_viral_production_dataframe <- determine_vp_VIPCAL_average_replicates(average_replicate_dataframe_with_timepoints)
  
  viral_production_VIPCAL <- determine_viral_production_dataframe %>%
    dplyr::group_by(.data$tag, .data$Time_Range, .data$Population, .data$Sample_Type) %>%
    dplyr::arrange('tag',
                   factor(.data$Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                   factor(.data$Population, levels = .GlobalEnv$populations_to_analyze[grep('c_V', .GlobalEnv$populations_to_analyze)])) %>%
    dplyr::mutate(VP_Method = 'VPCL_AR_DIFF') %>%
    dplyr::select('tag', 'Location', 'Station_Number', 'Depth', dplyr::everything())
  
  return(viral_production_VIPCAL)
}


#' @rdname vp_methods_VIPCAL
vp_VIPCAL_average_replicates_diff_SE <- function(data){
  average_replicate_dataframe_with_timepoints <- vp_average_replicate_dataframe(data)
  
  determine_viral_production_dataframe <- determine_vp_VIPCAL_average_replicates_SE(average_replicate_dataframe_with_timepoints)
  
  viral_production_VIPCAL <- determine_viral_production_dataframe %>%
    dplyr::group_by(.data$tag, .data$Time_Range, .data$Population, .data$Sample_Type) %>%
    dplyr::arrange('tag',
                   factor(.data$Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                   factor(.data$Population, levels = .GlobalEnv$populations_to_analyze[grep('c_V', .GlobalEnv$populations_to_analyze)])) %>%
    dplyr::mutate(VP_Method = 'VPCL_AR_DIFF_SE') %>%
    dplyr::select('tag', 'Location', 'Station_Number', 'Depth', dplyr::everything())
  
  return(viral_production_VIPCAL)
}


#' @rdname vp_methods_VIPCAL
vp_VIPCAL_average_replicates_diff_LMER <- function(data){
  separate_replicate_dataframe_with_timepoints <- vp_separate_replicate_dataframe(data)
  
  determine_viral_production_dataframe <- determine_vp_VIPCAL_LMER_model(separate_replicate_dataframe_with_timepoints)
  
  viral_production_VIPCAL <- determine_viral_production_dataframe %>%
    dplyr::group_by(.data$tag, .data$Time_Range, .data$Population, .data$Sample_Type) %>%
    dplyr::arrange('tag',
                   factor(.data$Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                   factor(.data$Population, levels = .GlobalEnv$populations_to_analyze[grep('c_V', .GlobalEnv$populations_to_analyze)])) %>%
    dplyr::mutate(VP_Method = 'VPCL_AR_DIFF_LMER') %>%
    dplyr::select('tag', 'Location', 'Station_Number', 'Depth', dplyr::everything())
  
  return(viral_production_VIPCAL)
}


#' @rdname vp_methods_VIPCAL
vp_VIPCAL_average_replicates_diff_LMER_SE <- function(data){
  separate_replicate_dataframe_with_timepoints <- vp_separate_replicate_dataframe(data)
  
  determine_viral_production_dataframe <- determine_vp_VIPCAL_LMER_model_SE(separate_replicate_dataframe_with_timepoints)
  
  viral_production_VIPCAL <- determine_viral_production_dataframe %>%
    dplyr::group_by(.data$tag, .data$Time_Range, .data$Population, .data$Sample_Type) %>%
    dplyr::arrange('tag',
                   factor(.data$Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                   factor(.data$Population, levels = .GlobalEnv$populations_to_analyze[grep('c_V', .GlobalEnv$populations_to_analyze)])) %>%
    dplyr::mutate(VP_Method = 'VPCL_AR_DIFF_LMER_SE') %>%
    dplyr::select('tag', 'Location', 'Station_Number', 'Depth', dplyr::everything())
  
  return(viral_production_VIPCAL)
}

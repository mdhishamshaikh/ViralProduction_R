#' Viral production calculation by linear regression
#' 
#' @description
#' Two main methods are used to determine the viral production rate over time in our assay: `Linear Regression` vs
#' `VIPCAL`. Different variants of both methods are performed, variance in replicate treatment, use of standard error
#' or type of estimation of difference curve is presented. Next, the different variants of `Linear Regression`,
#' which uses the slope between the viral counts to determine the viral production rate. The following general 
#' step-by-step plan is carried out by each of the variants:
#'    
#' 1. Constructing correct data frame depending on replicate treatment.
#' 2. Calculate viral production depending on replicate treatment, standard error use and difference curve estimation.
#' If no difference curve estimation, calculate lysogenic viral production.
#' 3. Arrange output data frame.
#'    
#' `vp_linear_allpoints` considers all points: no replicate treatment, no SE taken into account, no difference curve estimation.
#' 
#' `vp_linear_separate_replicates`: separate replicate treatment, no SE taken into account, no difference curve estimation. 
#' 
#' `vp_linear_average_replicates`: average replicate treatment, no SE taken into account, no difference curve estimation. 
#' 
#' `vp_linear_average_replicates_diff`: average replicate treatment, no SE taken into account, difference curve estimation by subtraction.
#' 
#' `vp_linear_average_replicates_diff_LMER`: average replicate treatment, no SE taken into account, difference curve estimation by LMER model.
#' 
#' More details about the used functions:
#' 
#' - Constructing data frames: [viralprod::vp_dataframes]
#' - Calculating viral production with linear regression: [viralprod::determine_vp_LM]
#' - LMER model: [viralprod::vp_LMER_model]
#' - Calculating lyosgenic production: [viralprod::vp_calculate_difference_samples]
#' 
#' @param data Data frame with the output of the flow cytometer.
#' @param AVG Interested in the lytic and lysogenic viral production. To study the lysogenic viral production, 
#' in need of difference samples. If \code{TRUE}, average over replicates after determining viral production with
#' separate replicate treatment and lysogenic viral production is calculated. If \code{FALSE}, no averaging over replicates,
#' output will consists of only VP and VPC samples with separate replicate treatment. (Default = \code{TRUE})
#'
#' @return Data frame with the viral production rate and the absolute viral production for each population at given time range of the assay.
#' 
#' @name vp_methods_linear
#' @rdname vp_methods_LM
#'
#' @examples \dontrun{
#' data_NJ2020 <- read.csv(system.file('extdata', 'NJ2020_subset.csv', package = "viralprod"))
#' 
#' vp_linear_allpoints(data_NJ2020)
#' 
#' vp_linear_separate_replicates(data_NJ2020)
#' vp_linear_separate_replicates(data_NJ2020, AVG = F)
#' 
#' vp_linear_average_replicates(data_NJ2020)
#' 
#' vp_linear_average_replicates_diff(data_NJ2020)
#' 
#' vp_linear_average_replicates_diff_LMER(data_NJ2020)
#' }
vp_linear_allpoints <- function(data){
  separate_replicate_dataframe_with_timepoints <- vp_separate_replicate_dataframe(data)
  
  determine_viral_production_dataframe <- determine_vp_linear_allpoints(separate_replicate_dataframe_with_timepoints) %>%
    vp_calculate_difference_samples()
  
  viral_production_LM <- determine_viral_production_dataframe %>%
    dplyr::group_by(.data$tag, .data$Time_Range, .data$Population, .data$Sample_Type) %>%
    dplyr::arrange('tag',
                   factor(.data$Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                   factor(.data$Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3'))) %>%
    dplyr::mutate(VP_method = 'LM_ALLPOINTS') %>%
    dplyr::select(.data$tag, .data$Location, .data$Station_Number, .data$Depth, dplyr::everything())
  
  return(viral_production_LM)
}


#' @rdname vp_methods_LM
vp_linear_separate_replicates <- function(data, AVG = TRUE){
  separate_replicate_dataframe_with_timepoints <- vp_separate_replicate_dataframe(data)
  
  determine_viral_production_dataframe <- determine_vp_linear_separate_replicates(separate_replicate_dataframe_with_timepoints)
  
  if (AVG == F){
    viral_production_LM <- determine_viral_production_dataframe %>%
      dplyr::group_by(.data$tag, .data$Time_Range, .data$Population, .data$Sample_Type) %>%
      dplyr::arrange('tag',
                     factor(.data$Sample_Type, levels = c('VP', 'VPC')),
                     factor(.data$Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3'))) %>%
      dplyr::mutate(VP_method = 'LM_SR') %>%
      dplyr::select(.data$tag, .data$Location, .data$Station_Number, .data$Depth, dplyr::everything())
  } else {
    viral_production_LM <- determine_viral_production_dataframe %>%
      dplyr::group_by(.data$tag, .data$Location, .data$Station_Number, .data$Depth, .data$Time_Range, .data$Population, .data$Sample_Type) %>%
      dplyr::summarise(
        VP_mean = mean(.data$VP),
        abs_VP_mean = mean(.data$abs_VP),
        VP_SE = plotrix::std.error(.data$VP),
        VP_R_Squared = mean(.data$VP_R_Squared)) %>%
      dplyr::rename(VP = .data$VP_mean, abs_VP = .data$abs_VP_mean) %>%
      vp_calculate_difference_samples() %>%
      dplyr::arrange('tag',
                     factor(.data$Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                     factor(.data$Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3'))) %>%
      dplyr::mutate(VP_method = 'LM_SR_AVG') %>%
      dplyr::select(.data$tag, .data$Location, .data$Station_Number, .data$Depth, dplyr::everything())
  }
  
  return(viral_production_LM)
}


#' @rdname vp_methods_LM
vp_linear_average_replicates <- function(data){
  average_replicate_dataframe_with_timepoints <- vp_average_replicate_dataframe(data)
  
  average_replicate_dataframe_no_diff_curve <- average_replicate_dataframe_with_timepoints %>%
    subset(average_replicate_dataframe_with_timepoints$Sample_Type != 'Diff')
  
  determine_viral_production_dataframe <- determine_vp_linear_average_replicates(average_replicate_dataframe_no_diff_curve) %>%
    vp_calculate_difference_samples()
  
  viral_production_LM <- determine_viral_production_dataframe %>%
    dplyr::group_by(.data$tag, .data$Time_Range, .data$Population, .data$Sample_Type) %>%
    dplyr::arrange('tag',
                   factor(.data$Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                   factor(.data$Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3'))) %>%
    dplyr::mutate(VP_method = 'LM_AR') %>%
    dplyr::select(.data$tag, .data$Location, .data$Station_Number, .data$Depth, dplyr::everything())
  
  return(viral_production_LM)
}
  

#' @rdname vp_methods_LM
vp_linear_average_replicates_diff <- function(data){
  average_replicate_dataframe_with_timepoints <- vp_average_replicate_dataframe(data)
  
  determine_viral_production_dataframe <- determine_vp_linear_average_replicates(average_replicate_dataframe_with_timepoints)
  
  viral_production_LM <- determine_viral_production_dataframe %>%
    dplyr::group_by(.data$tag, .data$Time_Range, .data$Population, .data$Sample_Type) %>%
    dplyr::arrange('tag',
                   factor(.data$Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                   factor(.data$Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3'))) %>%
    dplyr::mutate(VP_method = 'LM_AR_DIFF') %>%
    dplyr::select(.data$tag, .data$Location, .data$Station_Number, .data$Depth, dplyr::everything())
  
  return(viral_production_LM)
}


#' @rdname vp_methods_LM
vp_linear_average_replicates_diff_LMER <- function(data){
  separate_replicate_dataframe_with_timepoints <- vp_separate_replicate_dataframe(data)
  
  determine_viral_production_dataframe <- determine_vp_linear_LMER_model(separate_replicate_dataframe_with_timepoints)
  
  viral_production_LM <- determine_viral_production_dataframe %>%
    dplyr::group_by(.data$tag, .data$Time_Range, .data$Population, .data$Sample_Type) %>%
    dplyr::arrange('tag',
                   factor(.data$Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                   factor(.data$Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3'))) %>%
    dplyr::mutate(VP_method = 'LM_AR_DIFF_LMER') %>%
    dplyr::select(.data$tag, .data$Location, .data$Station_Number, .data$Depth, dplyr::everything())
  
  return(viral_production_LM)
}
  
#' Viral production calculation by VIPCAL
#' 
#' @description
#' Two main methods are used to determine the viral production rate over time in our assay: `Linear Regression` vs
#' `VIPCAL`. Different variants of both methods are performed, variance in replicate treatment, use of standard error
#' or type of estimation of difference curve is presented. Next, the different variants of `VIPCAL`,
#' which uses the average of increments between the viral counts to determine the viral production rate. 
#' The following general step-by-step plan is carried out by each of the variants:
#'    
#' 1. Constructing correct data frame depending on replicate treatment.
#' 2. Calculate viral production depending on replicate treatment, standard error use and difference curve estimation.
#' If no difference curve estimation, calculate lysogenic viral production.
#' 3. Arrange output dataframe.
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
#' - Constructing data frames: [viralprod::vp_dataframes]
#' - Calculating viral production with VIPCAL: [viralprod::determine_vp_VIPCAL]
#' - LMER model: [viralprod::vp_LMER_model]
#' - Calculating lyosgenic production: [viralprod::vp_calculate_difference_samples]
#' 
#' @param data Dataframe with the output of the flow cytometer (Step 1).
#' @param AVG Interested in the lytic and lysogenic viral production. To study the lysogenic viral production, 
#' in need of difference samples. If \code{TRUE}, average over replicates after determining viral production with
#' separate replicate treatment and lysogenic viral production is calculated. If \code{FALSE}, no averaging over replicates,
#' output will consits of only VP and VPC samples with separate replicate treatment. (Default = \code{TRUE})
#'
#' @return Dataframe with the viral production rate and the absolute viral production for each population at given time range of the assay.
#' @export
#' 
#' @name vp_methods_VIPCAL
#' @rdname vp_methods_VIPCAL
#'
#' @examples \dontrun{
#' data_NJ2020 <- read.csv(system.file('extdata', 'NJ2020_subset.csv', package = "viralprod"))
#' 
#' vp_VIPCAL_separate_replicates(data_NJ2020)
#' vp_VIPCAL_separate_replicates(data_NJ2020, AVG = F)
#' 
#' vp_VIPCAL_average_replicates(data_NJ2020)
#' 
#' vp_VIPCAL_average_replicates_SE(data_NJ2020)
#' 
#' vp_VIPCAL_average_replicates_diff(data_NJ2020)
#' 
#' vp_VIPCAL_average_replicates_diff_SE(data_NJ2020)
#' 
#' vp_VIPCAL_average_replicates_diff_LMER(data_NJ2020)
#' 
#' vp_VIPCAL_average_replicates_diff_LMER_SE(data_NJ2020)
#' }
vp_VIPCAL_separate_replicates <- function(data, AVG = T){
  separate_replicate_dataframe_with_timepoints <- vp_separate_replicate_dataframe(data)
  
  determine_viral_production_dataframe <- determine_vp_VIPCAL_separate_replicates(separate_replicate_dataframe_with_timepoints)
  
  if (AVG == F){
    viral_production_VIPCAL <- determine_viral_production_dataframe %>%
      dplyr::group_by(.data$tag, .data$Time_Range, .data$Population, .data$Sample_Type) %>%
      dplyr::arrange('tag',
                     factor(.data$Sample_Type, levels = c('VP', 'VPC')),
                     factor(.data$Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3'))) %>%
      dplyr::mutate(VP_method = 'VPCL_SR') %>%
      dplyr::select(.data$tag, .data$Location, .data$Station_Number, .data$Depth, dplyr::everything())
  } else {
    viral_production_VIPCAL <- determine_viral_production_dataframe %>%
      dplyr::group_by(.data$tag, .data$Location, .data$Station_Number, .data$Depth, .data$Time_Range, .data$Population, .data$Sample_Type) %>%
      dplyr::summarise(
        VP_mean = mean(.data$VP),
        abs_VP_mean = mean(.data$abs_VP),
        VP_SE = plotrix::std.error(.data$VP)) %>%
      dplyr::rename(VP = .data$VP_mean, abs_VP = .data$abs_VP_mean) %>%
      vp_calculate_difference_samples(VIPCAL = T) %>%
      dplyr::arrange('tag',
                     factor(.data$Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                     factor(.data$Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3'))) %>%
      dplyr::mutate(VP_method = 'VPCL_SR_AVG') %>%
      dplyr::select(.data$tag, .data$Location, .data$Station_Number, .data$Depth, dplyr::everything())
  }
  
  return(viral_production_VIPCAL)
}


#' @export
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
                   factor(.data$Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3'))) %>%
    dplyr::mutate(VP_method = 'VPCL_AR') %>%
    dplyr::select(.data$tag, .data$Location, .data$Station_Number, .data$Depth, dplyr::everything())
  
  return(viral_production_VIPCAL)
}


#' @export
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
                   factor(.data$Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3'))) %>%
    dplyr::mutate(VP_method = 'VPCL_AR_SE') %>%
    dplyr::select(.data$tag, .data$Location, .data$Station_Number, .data$Depth, dplyr::everything())
  
  return(viral_production_VIPCAL)
}


#' @export
#' @rdname vp_methods_VIPCAL
vp_VIPCAL_average_replicates_diff <- function(data){
  average_replicate_dataframe_with_timepoints <- vp_average_replicate_dataframe(data)
  
  determine_viral_production_dataframe <- determine_vp_VIPCAL_average_replicates(average_replicate_dataframe_with_timepoints)
  
  viral_production_VIPCAL <- determine_viral_production_dataframe %>%
    dplyr::group_by(.data$tag, .data$Time_Range, .data$Population, .data$Sample_Type) %>%
    dplyr::arrange('tag',
                   factor(.data$Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                   factor(.data$Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3'))) %>%
    dplyr::mutate(VP_method = 'VPCL_AR_DIFF') %>%
    dplyr::select(.data$tag, .data$Location, .data$Station_Number, .data$Depth, dplyr::everything())
  
  return(viral_production_VIPCAL)
}


#' @export
#' @rdname vp_methods_VIPCAL
vp_VIPCAL_average_replicates_diff_SE <- function(data){
  average_replicate_dataframe_with_timepoints <- vp_average_replicate_dataframe(data)
  
  determine_viral_production_dataframe <- determine_vp_VIPCAL_average_replicates_SE(average_replicate_dataframe_with_timepoints)
  
  viral_production_VIPCAL <- determine_viral_production_dataframe %>%
    dplyr::group_by(.data$tag, .data$Time_Range, .data$Population, .data$Sample_Type) %>%
    dplyr::arrange('tag',
                   factor(.data$Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                   factor(.data$Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3'))) %>%
    dplyr::mutate(VP_method = 'VPCL_AR_DIFF_SE') %>%
    dplyr::select(.data$tag, .data$Location, .data$Station_Number, .data$Depth, dplyr::everything())
  
  return(viral_production_VIPCAL)
}

#' @export
#' @rdname vp_methods_VIPCAL
vp_VIPCAL_average_replicates_diff_LMER <- function(data){
  separate_replicate_dataframe_with_timepoints <- vp_separate_replicate_dataframe(data)
  
  determine_viral_production_dataframe <- determine_vp_VIPCAL_LMER_model(separate_replicate_dataframe_with_timepoints)
  
  viral_production_VIPCAL <- determine_viral_production_dataframe %>%
    dplyr::group_by(.data$tag, .data$Time_Range, .data$Population, .data$Sample_Type) %>%
    dplyr::arrange('tag',
                   factor(.data$Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                   factor(.data$Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3'))) %>%
    dplyr::mutate(VP_method = 'VPCL_AR_DIFF_LMER') %>%
    dplyr::select(.data$tag, .data$Location, .data$Station_Number, .data$Depth, dplyr::everything())
  
  return(viral_production_VIPCAL)
}


#' @export
#' @rdname vp_methods_VIPCAL
vp_VIPCAL_average_replicates_diff_LMER_SE <- function(data){
  separate_replicate_dataframe_with_timepoints <- vp_separate_replicate_dataframe(data)
  
  determine_viral_production_dataframe <- determine_vp_VIPCAL_LMER_model_SE(separate_replicate_dataframe_with_timepoints)
  
  viral_production_VIPCAL <- determine_viral_production_dataframe %>%
    dplyr::group_by(.data$tag, .data$Time_Range, .data$Population, .data$Sample_Type) %>%
    dplyr::arrange('tag',
                   factor(.data$Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                   factor(.data$Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3'))) %>%
    dplyr::mutate(VP_method = 'VPCL_AR_DIFF_LMER_SE') %>%
    dplyr::select(.data$tag, .data$Location, .data$Station_Number, .data$Depth, dplyr::everything())
  
  return(viral_production_VIPCAL)
}

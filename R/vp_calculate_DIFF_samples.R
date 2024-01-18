#' Calculate lysogenic viral production
#' 
#' @description
#' The viral reduction assay has two types of samples: VP and VPC. In VP samples, count of bacteriophages in the lytic phase
#' can be measured (`lytic viral production`), while in VPC samples, count of bacteriophages in both the lytic and
#' lysogenic phase can be measured (`lytic + lysogenic viral production`) since treatment with antibiotic `Mitomycin-C`.
#' 
#' Bacteriophages in the lysogenic phase integrate with the genome of the bacteria and can't be measured normally. Mitomycin-C
#' inhibits DNA synthesis in the bacteria, therefore the bacteriophage needs to go into the lytic phase and measurement is 
#' possible. The details around the different calculation methods of viral production, being it either with linear regression
#' or VIPCAL, are available on: [viralprod::determine_vp_linear_regression] and [viralprod::determine_vp_VIPCAL]. 
#' 
#' Some of the variants estimate the difference curve by subtraction or LMER model to calculate the lysogenic viral production. 
#' If there is no estimation of the difference curve, lysogenic viral production is calculated afterwards as the difference 
#' between the viral production values of VPC samples and VP samples.
#'
#' @param DF Data frame containing the viral production for VP and VPC samples calculated with linear regression or VIPCAL.
#' @param VIPCAL If \code{FALSE}, viral production is calculated with linear regression. If viral production is calculated with VIPCAL, set to \code{TRUE}. (Default = \code{FALSE})
#' @param SE If \code{FALSE}, viral production is calculated with VIPCAL without taking the standard error into account. If VIPCAL-SE is used, set to \code{TRUE}. (Default = \code{FALSE})
#'
#' @return Data frame with the viral production rate and the absolute viral production for each population at the different time points of the assay. 
#' 
#' @name vp_calculate_difference_samples
#' @rdname vp_calc_DIFF
#'
#' @examples \dontrun{
#' data_NJ2020_all <- read.csv(system.file('extdata', 
#' 'NJ2020_Station_2_and_6_all_populations.csv', package = "viralprod"))
#' 
#' vp_check_populations(data_NJ2020_all)
#' 
#' # For average replicate treatment, no difference curve estimation by subtraction
#' # Take subset of data frame
#' DF_SR <- vp_separate_replicate_dataframe(data_NJ2020_all)
#' DF_AVG <- vp_average_replicate_dataframe(data_NJ2020_all) %>% subset(Sample_Type != 'Diff') 
#' 
#' vp_linear_allpoints <- determine_vp_linear_allpoints(DF_SR)
#' vp_VIPCAL_AR <- determine_vp_VIPCAL_average_replicates(DF_AVG)
#' vp_VIPCAL_AR_SE <- determine_vp_VIPCAL_average_replicates_SE(DF_AVG)
#' 
#' vp_calculate_difference_samples(vp_linear_allpoints)
#' vp_calculate_difference_samples(vp_VIPCAL_AR, VIPCAL = T)
#' vp_calculate_difference_samples(vp_VIPCAL_AR_SE, VIPCAL = T, SE = T)
#' }
vp_calculate_difference_samples <- function(DF, 
                                            VIPCAL = FALSE, 
                                            SE = FALSE){
  if (VIPCAL == T){
    if (SE == T){
      difference_samples <- DF %>%
        dplyr::group_by(.data$tag, .data$Location, .data$Station_Number, .data$Depth, .data$Time_Range, .data$Population) %>%
        dplyr::summarise(
          VP = sum(.data$VP[.data$Sample_Type == 'VPC']) - sum(.data$VP[.data$Sample_Type == "VP"]),
          abs_VP = sum(.data$abs_VP[.data$Sample_Type == 'VPC']) - sum(.data$abs_VP[.data$Sample_Type == "VP"]),
          VP_SE = sqrt((sum(.data$VP_SE[.data$Sample_Type == "VPC"]))^2 + (sum(.data$VP_SE[.data$Sample_Type == "VP"]))^2), 
          Sample_Type = "Diff")
    } else {
      difference_samples <- DF %>%
        dplyr::group_by(.data$tag, .data$Location, .data$Station_Number, .data$Depth, .data$Time_Range, .data$Population) %>%
        dplyr::summarise(
          VP = sum(.data$VP[.data$Sample_Type == 'VPC']) - sum(.data$VP[.data$Sample_Type == "VP"]),
          abs_VP = sum(.data$abs_VP[.data$Sample_Type == 'VPC']) - sum(.data$abs_VP[.data$Sample_Type == "VP"]), 
          Sample_Type = "Diff")
    }
  } else {
    difference_samples <- DF %>%
      dplyr::group_by(.data$tag, .data$Location, .data$Station_Number, .data$Depth, .data$Time_Range, .data$Population) %>%
      dplyr::summarise(
        VP = sum(.data$VP[.data$Sample_Type == 'VPC']) - sum(.data$VP[.data$Sample_Type == "VP"]),
        abs_VP = sum(.data$abs_VP[.data$Sample_Type == 'VPC']) - sum(.data$abs_VP[.data$Sample_Type == "VP"]),
        VP_SE = sqrt((sum(.data$VP_SE[.data$Sample_Type == "VPC"]))^2 + (sum(.data$VP_SE[.data$Sample_Type == "VP"]))^2), 
        VP_R_Squared = NA,
        Sample_Type = "Diff")
  }
  full_dataframe <- dplyr::full_join(DF, difference_samples, by = NULL)
  
  return(full_dataframe)
}

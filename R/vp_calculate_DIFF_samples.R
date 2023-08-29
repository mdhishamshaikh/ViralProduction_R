#' Calculate lysogenic viral production
#' 
#' Viral production can be divided in two phases: lytic and lysogenic viral production. The VP samples
#' cover the lytic viral production, VPC samples include both phases since treatment with the antibiotic `mitomycin C`. 
#' If the bacteriophage is in the lysogenic phase, it integrates with the genome of the bacteria and 
#' can't be measured. Since mitomycin C inhibits DNA synthesis in the bacteria, bacteriophages go into the lytic 
#' phase and measurement is possible. Viral production is calculated by different variants of `Linear Regression`
#' or `VIPCAL`. Some of these variants have the calculation of the lysogenic viral production included by estimating
#' a difference curve by subtraction or LMER model. For the other variants, calculate lysogenic viral production 
#' afterwards.
#'
#' @param DF Dataframe containing the viral production for VP and VPC samples calculated with linear regression or VIPCAL.
#' @param VIPCAL If \code{FALSE}, viral production is calculated with linear regression. If viral production is calculated with VIPCAL, set to \code{TRUE}. (Default = \code{FALSE})
#' @param SE If \code{FALSE}, viral production is calculated with VIPCAL without taking the standard error into account. If VIPCAL takes the SE into account, set to \code{TRUE}. (Default = \code{FALSE})
#'
#' @return Dataframe with the viral production added for lysogenic phase. In the column Sample_Type a new sample, `Diff` is introduced that represents the lysogenic viral production. 
#' 
#' @name vp_calculate_difference_samples
#' @rdname vp_calc_DIFF
#'
#' @examples \dontrun{
#' data_NJ2020 <- read.csv(system.file('extdata', 'NJ2020_subset.csv', package = "viralprod"))
#' DF_SR <- vp_separate_replicate_dataframe(data_NJ2020)
#' DF_AVG <- vp_average_replicate_dataframe(data_NJ2020) %>% subset(Sample_Type != 'Diff')
#' 
#' vp_linear_allpoints <- determine_vp_linear_allpoints(DF_SR)
#' vp_VIPCAL_AR <- determine_vp_VIPCAL_average_replicates(DF_AVG)
#' vp_VIPCAL_AR_SE <- determine_vp_VIPCAL_average_replicates_SE(DF_AVG)
#' 
#' vp_calculate_difference_samples(vp_linear_allpoints)
#' vp_calculate_difference_samples(vp_VIPCAL_AR, VIPCAL = T)
#' vp_calculate_difference_samples(vp_VIPCAL_AR_SE, VIPCAL = T, SE = T)
#' }
vp_calculate_difference_samples <- function(DF, VIPCAL = FALSE, SE = FALSE){
  if (VIPCAL == T){
    if (SE == T){
      difference_samples <- DF %>%
        dplyr::group_by(.data$tag, .data$Location, .data$Station_Number, .data$Depth, .data$Time_Range, .data$Population) %>%
        dplyr::summarise(
          VP = sum(.data$VP[.data$Sample_Type == 'VPC']) - sum(.data$VP[.data$Sample_Type == "VP"]),
          abs_VP = sum(.data$abs_VP[.data$Sample_Type == 'VPC']) - sum(.data$abs_VP[.data$Sample_Type == "VP"]),
          VP_SE = sum(.data$VP_SE[.data$Sample_Type == "VPC"]) + sum(.data$VP_SE[.data$Sample_Type == "VP"]), 
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
        VP_SE = sum(.data$VP_SE[.data$Sample_Type == "VPC"]) + sum(.data$VP_SE[.data$Sample_Type == "VP"]), 
        VP_R_Squared = NA,
        Sample_Type = "Diff")
  }
  full_dataframe <- dplyr::full_join(DF, difference_samples, by = NULL)
  
  return(full_dataframe)
}
#' End-to-End viral production analysis
#' 
#' @description
#' Wrapper function that includes everything from the `viralprod` package. Calculation, analyzing and
#' visualization of viral production based on the input data from the flow cytometry and the abundances
#' of bacterial and virus populations in the original sample. All outputs will be saved together in the 
#' folder specified by the argument `output_dir`.
#' 
#' More details about the previous performed steps can be found on:
#' 
#' 1. Calculate viral production step: [viralprod::vp_calculate]
#' 2. Analyze viral production step: [viralprod::vp_analyze]
#' 3. Visualize viral production data: [viralprod::vp_visualize] 
#' 
#' @param data Data frame with the output of the flow cytometry.
#' @param original_abundances Data frame with the abundances of bacterial and virus population in the original sample.
#' @param methods Integer vector with the indexes of `list_of_methods`. Indexes determine which methods will be performed within the viral production calculation. 
#' @param SR_calc If \code{TRUE}, separate replicate treatment results will be stored and saved as an separate data frame. 
#' Set to \code{FALSE}, if separate replicate treatment results are not wanted. (Default = \code{TRUE})
#' @param BP_endpoint If \code{TRUE}, the bacterial endpoint will be taken into account and only those results will be saved in a new
#' data frame. If not wanted, set to \code{FALSE}. (Default = \code{TRUE})
#' @param burst_sizes Vector with three different burst sizes. The burst size refers to the number of new viral particles released from an infected bacterial cell. 
#' @param bacterial_secondary_production Value for the bacterial secondary production, how much new bacterial biomass is produced as a result of bacterial growth and reproduction. 
#' @param nutrient_content_bacteria List with the amount of organic carbon, nitrogen and phosphor released by a bacteria, preferred a aquatic, North Sea bacteria.
#' @param nutrient_content_viruses List with the amount of organic carbon, nitrogen and phosphor released by a marine virus (bacteriophage).
#' @param write_output If \code{TRUE}, the output data frames will be saved as csv files in the same folder of the viral production calculation.
#' If no csv files are wanted, set to \code{FALSE}. (Default = \code{TRUE})
#' @param output_dir String that refers to the location of folder to save the data frames as csv files. 
#'
#' @return All outputs will be saved in the folder specified by `output_dir`.
#' @export
#' 
#' @name vp_end_to_end
#' @rdname vp_end_to_end
#'
#' @examples \dontrun{
#' # Setup
#' data_NJ2020_all <- read.csv(system.file('extdata', 
#' 'NJ2020_Station_2_and_6_all_populations.csv', package = "viralprod"))
#' 
#' original_abundances_NJ2020 <- read.csv(system.file('extdata',
#' 'NJ2020_original_abundances.csv', package = "viralprod"))
#' 
#' # Perform
#' vp_end_to_end <- function(
#' data_NJ2020_all,
#' original_abundances_NJ2020,
#' output_dir = paste0(system.file(“extdata”, package = “viralprod”), 
#' “/vp_results_NJ2020_all”))
#' }
vp_end_to_end <- function(data = data.frame(),
                          original_abundances = data.frame(),
                          methods = c(1:12),
                          SR_calc = TRUE,
                          BP_endpoint = TRUE,
                          burst_sizes = c(),
                          bacterial_secondary_production = NULL,
                          nutrient_content_bacteria = list(),
                          nutrient_content_viruses = list(),
                          write_output = TRUE,
                          output_dir = ''){
  ## 1. Check the input data frames
  data_viralprod <- vp_class_count_data(data)
  original_abundances_viralprod <- vp_class_ori_abu(original_abundances)
  
  ## 2. Calculate viral production
  vp_calculate(x = data_viralprod,
                             methods = methods,
                             SR_calc = SR_calc,
                             BP_endpoint = BP_endpoint,
                             write_output = write_output,
                             output_dir = output_dir)
  
  ## 2. Analyze viral production
  vp_analyze(x = data_viralprod,
                           vp_results = .GlobalEnv$vp_results_output_df,
                           original_abundances = original_abundances_viralprod,
                           burst_sizes = burst_sizes,
                           bacterial_secondary_production = bacterial_secondary_production,
                           nutrient_content_bacteria = nutrient_content_bacteria,
                           nutrient_content_viruses = nutrient_content_viruses,
                           write_output = write_output,
                           output_dir = output_dir)
  
  ## 3. Visualize viral production
  vp_visualize(x = data_viralprod,
                             vp_results = .GlobalEnv$vp_results_output_df,
                             original_abundances = original_abundances_viralprod,
                             burst_sizes = burst_sizes,
                             bacterial_secondary_production = bacterial_secondary_production,
                             nutrient_content_bacteria = nutrient_content_bacteria,
                             nutrient_content_viruses = nutrient_content_viruses,
                             write_output = write_output,
                             output_dir = output_dir)
}

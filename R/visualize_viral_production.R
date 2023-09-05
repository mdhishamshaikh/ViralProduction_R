#' Wrapper function to visualize viral production data
#' 
#' @description
#' A major step in data analyses is `Data Visualization`. Some suggestions to visualize the flow cytometry, viral
#' production or analyzed viral production data are given. The different visualizations will be performed, stored  
#' and saved as PDFs if wanted.
#' 
#' More details about the previous performed steps can be found on:
#' 
#' 1. Calculate viral production step: [viralprod::calculate_viral_production]
#' 2. Analyze viral production step: [viralprod::analyze_viral_production]
#' 3. Visualizations of viral production data: [viralprod::vp_visuals] 
#' 
#' @param data Data frame with the output of the flow cytometry.
#' @param original_abundances Data frame with the abundances of bacterial and virus population in the original sample.
#' @param write_csv If \code{TRUE}, the output data frames will be saved as csv files in the same folder of the viral production calculation.
#' If no csv files are wanted, set to \code{FALSE}. (Default = \code{TRUE})
#' @param output_dir String that refers to the location of folder to save the data frames as csv files. 
#' 
#' @return Plot objects will be stored in variable `plot_list` in the global environment. A subfolder, `Figures`, will be created in the given output folder for storing the figures as PDFs.
#' @export
#' 
#' @name visualize_viral_production
#' @rdname visualize_viral_production
#'
#' @examples \dontrun{
#' # Setup
#' data_NJ2020_all <- read.csv(system.file('extdata', 
#' 'NJ2020_Station_2_and_6_all_populations.csv', package = "viralprod"))
#' 
#' original_abundances_NJ2020 <- read.csv(system.file('extdata',
#' 'NJ2020_original_abundances.csv', package = "viralprod"))
#' 
#' #' # Perform
#' # Write output files
#' visualize_viral_production(data_NJ2020_all, original_abundances_NJ2020, 
#' output_dir = paste0(system.file(“extdata”, package = “viralprod”), 
#' “/NJ2020_vp_results”))
#' 
#' # No output files
#' visualize_viral_production(data_NJ2020_all, original_abundances_NJ2020, write_csv = F)
#' }
visualize_viral_production <- function(data = data.frame(),
                                       original_abundances = data.frame(),
                                       write_csv = TRUE,
                                       output_dir = ''){
  ## 1. Checks
  # Check for valid output directory if PDFs needs to be written
  if (write_csv == T){
    if(output_dir == ''){
      stop('No output directory is given, please define output directory before proceeding!')
      
    } else if(!file.exists(output_dir)){
      # Give the user an option if he wants to continue and store figures in separate folder
      user_choice <- utils::menu(
        c("Continue and store figures in folder by your choice",
          "Stop and define same output directory of calculation results"),
        title = message("Warning: The output folder does not exists!"),
        graphics = FALSE)
      
      if (user_choice == 1){
        message('Continuing with the given output directory, analyzing results will be stored in separate folder!')
        dir.create(output_dir)
        visualize_vp_results_path <- paste0(output_dir, '/Figures/')
        
      } else {
        stop('Process stopped, please define another output directory!')
      }
    } else {
      visualize_vp_results_path <- paste0(output_dir, '/Figures/')
    }
  }
  
  # Check if both data frames aren't empty
  if (any(sapply(list(data, original_abundances), function(df) rlang::is_empty(df)))){
    stop('Not able to proceed analyzing since one of the input data frames is empty!')
  }
  
  ## 2. Visuals
  .GlobalEnv$plot_list <- list()
  
  plot_overview_counts_over_time(data)
  plot_collision_rates(data, original_abundances)
  
  calculate_viral_production(data, write_csv = F)
  plot_comparison_methods(.GlobalEnv$vp_results_output_df)
  
  analyze_viral_production(.GlobalEnv$vp_results_output_BP_df, data, original_abundances, write_csv = F)
  plot_percentage_cells(.GlobalEnv$analyzed_vp_results_df)
  
  analyze_viral_production(.GlobalEnv$vp_results_output_T24_df, data, original_abundances, write_csv = F)
  plot_nutrient_release(.GlobalEnv$analyzed_vp_results_df)
  
  ## 3. Save figures as PDFs
  
}

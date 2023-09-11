#' Wrapper function to visualize viral production data
#' 
#' @description
#' A major step in data analyses is `data visualization`. Some suggestions to visualize the flow cytometry, viral
#' production or analyzed viral production data are given. The different visualizations will be performed, stored  
#' and saved as PDFs if wanted.
#' 
#' More details about the previous performed steps can be found on:
#' 
#' 1. Calculate viral production step: [viralprod::calculate_viral_production]
#' 2. Analyze viral production step: [viralprod::analyze_viral_production]
#' 3. Visualizations of viral production data: [viralprod::vp_visuals] 
#' 
#' @param x Data frame with the output of the flow cytometry, has to have the `viralprod` class.
#' @param ... Arguments passed on to the next function.
#' @param vp_results Data frame with the viral production calculation results, available in global environment. See [viralprod::calculate_viral_production] for more details.
#' @param original_abundances Data frame with the abundances of bacterial and virus population in the original sample, has to have the `viralprod_analyze` class.
#' @param burst_sizes Vector with three different burst sizes. The burst size refers to the number of new viral particles released from an infected bacterial cell. 
#' @param bacterial_secondary_production Value for the bacterial secondary production, how much new bacterial biomass is produced as a result of bacterial growth and reproduction. 
#' @param nutrient_content_bacteria List with the amount of organic carbon, nitrogen and phosphor released by a bacteria, preferred a aquatic, North Sea bacteria.
#' @param nutrient_content_viruses List with the amount of organic carbon, nitrogen and phosphor released by a marine virus (bacteriophage).
#' @param write_output If \code{TRUE}, the output data frames will be saved as csv files in the same folder of the viral production calculation.
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
#' x <- new_viralprod_class(data_NJ2020_all)
#' y <- new_viralprod_class_2(original_abundances_NJ2020)
#' calculate_viral_production(x, write_output = F)
#' 
#' # Perform
#' # Default method
#' visualize_viral_production(x = data_NJ2020_all, vp_results = vp_results_output_df, 
#' original_abundances = original_abundances_NJ2020, write_output = F)
#' 
#' # S3 class, viralprod, method
#' # Error expected when original abundances don't have correct class
#' visualize_viral_production(x, vp_results = vp_results_output_df, 
#' original_abundances = original_abundances_NJ2020, write_output = F)
#' 
#' # Write output files
#' visualize_viral_production(x, vp_results = vp_results_output_df, 
#' original_abundances = y, 
#' output_dir = paste0(system.file(“extdata”, package = “viralprod”), 
#' “/NJ2020_vp_results”))
#' 
#' # No output files
#' visualize_viral_production(x, vp_results = vp_results_output_df, 
#' original_abundances = y, write_output = F)
#' 
#' # Set own parameter values
#' visualize_viral_production(x, vp_results = vp_results_output_df, 
#' original_abundances = y,
#' burst_sizes = c(15,30,50), bacterial_secondary_production = 1000, 
#' nutrient_content_bacteria = list(C = 20, N = 15, P = 5),
#' nutrient_content_virus = list(C = 5, N = 3, P = 1)) 
#' }
visualize_viral_production <- function(x, ...){
  UseMethod("visualize_viral_production")
}


#' @export
#' @rdname visualize_viral_production
visualize_viral_production.default <- function(x, ...){
  message('The input data frame has not the right format, please adjust before moving on! You can check your data frame with `new_viralprod_class()`.')
}


#' @export
#' @rdname visualize_viral_production
visualize_viral_production.viralprod <- function(x, ...,
                                                 vp_results = data.frame(),
                                                 original_abundances = data.frame(),
                                                 burst_sizes = c(),
                                                 bacterial_secondary_production = NULL,
                                                 nutrient_content_bacteria = list(),
                                                 nutrient_content_viruses = list(),
                                                 write_output = TRUE,
                                                 output_dir = ''){
  ## 1. Checks
  # Check if original_abundances have correct class
  checkmate::assert_class(original_abundances, 'viralprod_analyze')
  
  # Check for valid output directory if PDFs needs to be written
  if (write_output == T){
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
        visualize_vp_results_path <- paste0(output_dir, '/')
        
      } else {
        stop('Process stopped, please define another output directory!')
      }
    } else {
      visualize_vp_results_path <- paste0(output_dir, '/Figures/')
      if (!file.exists(visualize_vp_results_path)){
        dir.create(visualize_vp_results_path)
      }
    }
  }
  
  # Check if all three data frames aren't empty
  if (any(sapply(list(x, vp_results, original_abundances), function(df) rlang::is_empty(df)))){
    stop('Not able to proceed analyzing since one of the input data frames is empty!')
  }
  
  ## 2. Visuals
  .GlobalEnv$plot_list <- list()
  
  plot_overview_counts_over_time(x)
  plot_collision_rates(x, original_abundances)
  plot_comparison_methods(vp_results)
  plot_VIPCAL_vs_VIPCAL_SE(vp_results)
  
  analyze_viral_production(x,
                           vp_results = .GlobalEnv$vp_results_output_BP_df,
                           original_abundances = original_abundances,
                           burst_sizes = burst_sizes,
                           bacterial_secondary_production = bacterial_secondary_production,
                           nutrient_content_bacteria = nutrient_content_bacteria,
                           nutrient_content_viruses = nutrient_content_viruses,
                           write_output = F)
  plot_percentage_cells(.GlobalEnv$analyzed_vp_results_df)
  
  analyze_viral_production(x,
                           vp_results = .GlobalEnv$vp_results_output_T24_df,
                           original_abundances = original_abundances,
                           burst_sizes = burst_sizes,
                           bacterial_secondary_production = bacterial_secondary_production,
                           nutrient_content_bacteria = nutrient_content_bacteria,
                           nutrient_content_viruses = nutrient_content_viruses,
                           write_output = F)
  plot_nutrient_release(.GlobalEnv$analyzed_vp_results_df)
  
  ## 3. Save figures as PDFs
  if(write_output == T){
    for (plot in 1:length(.GlobalEnv$plot_list)){
      plot_name <- names(.GlobalEnv$plot_list)[plot]
      file_path <- paste0(visualize_vp_results_path, paste0(plot_name, '.pdf'))
      ggplot2::ggsave(filename = file_path, 
                      plot = .GlobalEnv$plot_list[[plot]]$plot_object, 
                      width = .GlobalEnv$plot_list[[plot]]$width, 
                      height = .GlobalEnv$plot_list[[plot]]$height)
    }
  }
}

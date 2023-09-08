#' List of methods used for viral production calculation
#' 
#' @description
#' A list of all the different variants of linear regression and VIPCAL available to calculate the viral
#' production rate.
#' 
#' More details on the different methods:
#' 
#' - Variants of linear regression: [viralprod::vp_methods_linear]
#' - Variants of VIPCAL: [viralprod::vp_methods_VIPCAL]
#' 
#' @return A list, with the different variants of linear regression and VIPCAL to calculate viral production, will be available in the global environment.
#' 
#' @name vp_list_of_methods
#' @rdname vp_list_of_methods
#' 
#' @export
#' 
#' @examples \dontrun{
#' viralprod::vp_list_of_methods()
#' 
#' # If you want to run for example the first method
#' data_NJ2020_all <- read.csv(system.file('extdata', 
#' 'NJ2020_Station_2_and_6_all_populations.csv', package = "viralprod"))
#' 
#' vp_check_populations(data_NJ2020_all)
#' list_of_methods[[1]](data_NJ2020_all)
#' }
vp_list_of_methods <- function(){
  methods_to_calculate_viral_production <- list(
    "vp_linear_allpoints",
    "vp_linear_separate_replicates",
    "vp_linear_average_replicates",
    "vp_linear_average_replicates_diff",
    "vp_linear_average_replicates_diff_LMER",
    "vp_VIPCAL_separate_replicates",
    "vp_VIPCAL_average_replicates",
    "vp_VIPCAL_average_replicates_SE",
    "vp_VIPCAL_average_replicates_diff",
    "vp_VIPCAL_average_replicates_diff_SE",
    "vp_VIPCAL_average_replicates_diff_LMER",
    "vp_VIPCAL_average_replicates_diff_LMER_SE"
  )
  methods_list <- list()
  
  for (method in methods_to_calculate_viral_production){
    method_function <- get(method)
    methods_list[[length(methods_list) + 1]] <- method_function
  }
  names(methods_list) <- methods_to_calculate_viral_production
  
  .GlobalEnv$list_of_methods <- methods_list
}
